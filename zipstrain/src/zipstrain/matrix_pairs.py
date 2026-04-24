from __future__ import annotations

import gc
import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import duckdb
import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

from zipstrain import __version__
from zipstrain import compare as cp

PROFILE_REQUIRED_COLUMNS = {"chrom", "pos", "gene", "genome", "A", "T", "C", "G"}
COUNT_DTYPES = {
    "uint16": np.uint16,
    "uint32": np.uint32,
}
MATRIX_BUILD_MIN_COV = 5
FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS = "allele_presence_after_cov_filter"
MATRIX_PAIR_BACKENDS = ("numpy", "torch", "torch-cpu", "torch-cuda", "torch-mps")
BuildProgressCallback = Callable[[dict[str, object]], None]
CompareProgressCallback = Callable[[dict[str, object]], None]
MATRIX_COMPARISON_SUPPORTED_CALCULATIONS = ("ani", "ibs")
MATRIX_BUILD_FIXED_HEADROOM_BYTES = 128 * 1024 * 1024
MATRIX_BUILD_TEMP_BYTES_PER_POSITION = 32
MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES = 256 * 1024 * 1024
MATRIX_BUILD_MIN_COMMIT_BATCH_BYTES = 64 * 1024 * 1024
MATRIX_BUILD_DUCKDB_MEMORY_FRACTION = 0.25


def parse_matrix_calculations(calculate: Optional[str] = None) -> tuple[str, ...]:
    calculations = cp.parse_genome_calculations(calculate)
    supported = tuple(metric for metric in MATRIX_COMPARISON_SUPPORTED_CALCULATIONS if metric in calculations)
    if supported:
        return supported
    raise ValueError("Matrix compare currently supports only ani and ibs calculations.")


def matrix_metric_output_columns(calculate: Optional[str] = None) -> list[str]:
    columns = ["sample_1", "sample_2", "genome"]
    calculations = parse_matrix_calculations(calculate)
    if "ani" in calculations:
        columns.extend(["total_positions", "share_allele_pos", "genome_pop_ani"])
    if "ibs" in calculations:
        columns.append("max_consecutive_length")
    return columns


def matrix_pair_output_schema(calculate: Optional[str] = None) -> pa.Schema:
    fields = [
        pa.field("sample_1", pa.string()),
        pa.field("sample_2", pa.string()),
        pa.field("genome", pa.string()),
    ]
    calculations = parse_matrix_calculations(calculate)
    if "ani" in calculations:
        fields.extend(
            [
                pa.field("total_positions", pa.int64()),
                pa.field("share_allele_pos", pa.int64()),
                pa.field("genome_pop_ani", pa.float64()),
            ]
        )
    if "ibs" in calculations:
        fields.append(pa.field("max_consecutive_length", pa.int64()))
    return pa.schema(fields)


@dataclass(frozen=True)
class ScaffoldSpec:
    scaffold_idx: int
    genome: str
    chrom: str
    index_base: int
    vector_length: int
    min_pos: int
    max_pos: int


@dataclass(frozen=True)
class MatrixDbSummary:
    output_file: Path
    profile_files: int
    sample_count: int
    scaffold_count: int
    stored_rows: int
    count_dtype: str
    genome_scope: str
    memory_limit_gb: float


@dataclass(frozen=True)
class MatrixCompareSummary:
    output_file: Path
    requested_pairs: int
    written_rows: int
    scaffold_count: int
    genome_count: int
    anchor_groups: int
    target_chunks: int
    min_cov: int
    genome_scope: str
    backend: str
    device: str
    memory_limit_gb: float
    position_tile_size: int


@dataclass(frozen=True)
class MatrixBuildMemoryPlan:
    duckdb_memory_limit_bytes: int
    commit_batch_bytes: int
    estimated_python_peak_bytes: int
    limiting_scaffold: str


class TorchBackendMissingError(ImportError):
    pass


class MatrixPairComputeBackend:
    def __init__(self, backend: str):
        if backend not in MATRIX_PAIR_BACKENDS:
            raise ValueError(
                f"Unsupported backend '{backend}'. Choose one of {', '.join(MATRIX_PAIR_BACKENDS)}."
            )
        self.requested = backend
        self.kind = "numpy"
        self.device = "cpu"
        self.torch = None
        if backend == "numpy":
            return
        self.kind = "torch"
        self.torch = self._import_torch()
        self.device = self._resolve_torch_device(backend)

    def _import_torch(self):
        try:
            return importlib.import_module("torch")
        except ImportError as exc:
            raise TorchBackendMissingError(
                "Torch backend requested but 'torch' is not installed. "
                'Install with: pip install "zipstrain[matrix]". '
                "For CUDA, install the matching Torch wheel from pytorch.org."
            ) from exc

    def _torch_mps_available(self) -> bool:
        return bool(
            hasattr(self.torch, "backends")
            and hasattr(self.torch.backends, "mps")
            and self.torch.backends.mps.is_available()
        )

    def _resolve_torch_device(self, backend: str) -> str:
        if backend == "torch-cpu":
            return "cpu"
        if backend == "torch-cuda":
            if not self.torch.cuda.is_available():
                raise RuntimeError("Torch CUDA backend requested but CUDA is not available.")
            return "cuda"
        if backend == "torch-mps":
            if not self._torch_mps_available():
                raise RuntimeError("Torch MPS backend requested but MPS is not available.")
            return "mps"
        if self.torch.cuda.is_available():
            return "cuda"
        if self._torch_mps_available():
            return "mps"
        return "cpu"

def _looks_like_profile_parquet(path: Path) -> bool:
    try:
        names = set(pq.read_schema(path).names)
    except Exception:
        return False
    return PROFILE_REQUIRED_COLUMNS.issubset(names)


def discover_profile_parquets(profile_dir: Path) -> list[Path]:
    if not profile_dir.exists():
        raise FileNotFoundError(f"Profile directory does not exist: {profile_dir}")
    if not profile_dir.is_dir():
        raise NotADirectoryError(f"Profile path is not a directory: {profile_dir}")
    candidates = sorted(path for path in profile_dir.glob("*.parquet") if path.is_file())
    profiles = [path for path in candidates if _looks_like_profile_parquet(path)]
    if not profiles:
        raise ValueError(
            f"No classic ZipStrain profile parquets were found in {profile_dir}. "
            "Expected parquet files with columns chrom,pos,gene,genome,A,T,C,G."
        )
    sample_names = [path.stem for path in profiles]
    if len(sample_names) != len(set(sample_names)):
        raise ValueError("Profile file stems must be unique to derive unique sample names.")
    return profiles


def _collect_scaffold_specs(profile_paths: list[Path], genome: Optional[str] = None) -> list[ScaffoldSpec]:
    scope = pl.lit(True)
    if genome is not None:
        scope = scope & (pl.col("genome") == genome)
    scaffolds: dict[tuple[str, str], dict[str, int]] = {}
    for profile_path in profile_paths:
        frame = (
            pl.scan_parquet(profile_path)
            .filter(scope)
            .select("genome", "chrom", "pos")
            .group_by(["genome", "chrom"])
            .agg(
                pl.col("pos").min().cast(pl.Int64).alias("min_pos"),
                pl.col("pos").max().cast(pl.Int64).alias("max_pos"),
            )
            .collect(engine="streaming")
        )
        for row in frame.iter_rows(named=True):
            key = (str(row["genome"]), str(row["chrom"]))
            min_pos = int(row["min_pos"])
            max_pos = int(row["max_pos"])
            current = scaffolds.get(key)
            if current is None:
                scaffolds[key] = {"min_pos": min_pos, "max_pos": max_pos}
            else:
                current["min_pos"] = min(current["min_pos"], min_pos)
                current["max_pos"] = max(current["max_pos"], max_pos)
    if not scaffolds:
        detail = genome if genome is not None else "all"
        raise ValueError(f"No profile rows found for genome scope: {detail}")
    specs: list[ScaffoldSpec] = []
    for scaffold_idx, ((genome_name, chrom), bounds) in enumerate(sorted(scaffolds.items())):
        index_base = int(bounds["min_pos"])
        max_pos = int(bounds["max_pos"])
        specs.append(
            ScaffoldSpec(
                scaffold_idx=scaffold_idx,
                genome=genome_name,
                chrom=chrom,
                index_base=index_base,
                vector_length=max_pos - index_base + 1,
                min_pos=index_base,
                max_pos=max_pos,
            )
        )
    return specs


def _collect_scaffold_specs_from_bed(
    profile_paths: list[Path],
    bed_file: Path,
    genome: Optional[str] = None,
) -> list[ScaffoldSpec]:
    bed_path = Path(bed_file)
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file does not exist: {bed_path}")
    if not bed_path.is_file():
        raise FileNotFoundError(f"BED path is not a file: {bed_path}")

    bed_spans = (
        pl.scan_csv(bed_path, separator="\t", has_header=False)
        .select(
            pl.col("column_1").cast(pl.Utf8).alias("chrom"),
            pl.col("column_2").cast(pl.Int64).alias("start"),
            pl.col("column_3").cast(pl.Int64).alias("end"),
        )
        .group_by("chrom")
        .agg(
            pl.col("start").min().alias("min_start"),
            pl.col("end").max().alias("max_end"),
        )
        .collect(engine="streaming")
    )

    if bed_spans.height == 0:
        raise ValueError(f"No scaffold intervals found in BED file: {bed_path}")

    scope = pl.lit(True)
    if genome is not None:
        scope = scope & (pl.col("genome") == genome)

    profile_pairs = (
        pl.scan_parquet([str(path) for path in profile_paths])
        .filter(scope)
        .select(
            pl.col("genome").cast(pl.Utf8).alias("genome"),
            pl.col("chrom").cast(pl.Utf8).alias("chrom"),
        )
        .unique()
        .collect(engine="streaming")
    )

    if profile_pairs.height == 0:
        detail = genome if genome is not None else "all"
        raise ValueError(f"No profile rows found for genome scope: {detail}")

    specs_df = profile_pairs.join(bed_spans, on="chrom", how="left")
    missing_scaffolds = (
        specs_df
        .filter(pl.col("max_end").is_null() | pl.col("min_start").is_null())
        .select("chrom")
        .unique()
        .get_column("chrom")
        .to_list()
    )
    if missing_scaffolds:
        raise ValueError(
            "BED file is missing scaffold intervals for the following profile scaffolds: "
            + ", ".join(sorted(str(name) for name in missing_scaffolds))
        )

    specs: list[ScaffoldSpec] = []
    for scaffold_idx, row in enumerate(
        specs_df.sort(["genome", "chrom"]).iter_rows(named=True)
    ):
        min_start = int(row["min_start"])
        max_end = int(row["max_end"])
        # BED intervals are 0-based, half-open; classic profile positions are 1-based.
        index_base = min_start + 1
        vector_length = max_end - min_start
        if vector_length <= 0:
            raise ValueError(
                f"Invalid BED span for scaffold {row['chrom']}: start={min_start}, end={max_end}"
            )
        specs.append(
            ScaffoldSpec(
                scaffold_idx=scaffold_idx,
                genome=str(row["genome"]),
                chrom=str(row["chrom"]),
                index_base=index_base,
                vector_length=vector_length,
                min_pos=index_base,
                max_pos=max_end,
            )
        )
    return specs


def _pack_matrix(matrix: np.ndarray) -> bytes:
    return np.ascontiguousarray(matrix).tobytes()


def _unpack_matrix(blob: bytes, dtype_name: str, shape: tuple[int, int]) -> np.ndarray:
    dtype = COUNT_DTYPES[dtype_name]
    out = np.frombuffer(blob, dtype=dtype)
    expected = shape[0] * shape[1]
    if out.size != expected:
        raise ValueError(f"Matrix blob has {out.size} values, expected {expected}")
    return out.reshape(shape)


def _memory_limit_bytes(memory_limit_gb: float) -> int:
    if memory_limit_gb <= 0:
        raise ValueError("memory_limit_gb must be > 0")
    return int(memory_limit_gb * (1024 ** 3))


def _commit_batch_bytes(commit_batch_gb: float) -> int:
    if commit_batch_gb <= 0:
        raise ValueError("commit_batch_gb must be > 0")
    return int(commit_batch_gb * (1024 ** 3))


def _estimate_sample_scaffold_bytes(vector_length: int, dtype_name: str) -> int:
    return vector_length * 4 * np.dtype(COUNT_DTYPES[dtype_name]).itemsize


def _estimate_builder_python_peak_bytes(vector_length: int, dtype_name: str) -> int:
    matrix_bytes = _estimate_sample_scaffold_bytes(vector_length, dtype_name)
    temp_bytes = vector_length * MATRIX_BUILD_TEMP_BYTES_PER_POSITION
    return matrix_bytes * 2 + temp_bytes + MATRIX_BUILD_FIXED_HEADROOM_BYTES


def _format_memory_bytes(num_bytes: int) -> str:
    return f"{num_bytes / (1024 ** 3):.2f} GB"


def _duckdb_memory_limit_setting(memory_limit_bytes: int) -> str:
    memory_limit_mib = max(1, int(memory_limit_bytes // (1024 ** 2)))
    return f"{memory_limit_mib}MiB"


def _plan_matrix_build_memory(
    scaffolds: list[ScaffoldSpec],
    count_dtype: str,
    memory_limit_bytes: int,
    commit_batch_gb: Optional[float] = None,
) -> MatrixBuildMemoryPlan:
    if not scaffolds:
        raise ValueError("No scaffolds were provided for matrix DB build.")

    limiting_spec = max(
        scaffolds,
        key=lambda spec: _estimate_builder_python_peak_bytes(spec.vector_length, count_dtype),
    )
    estimated_python_peak_bytes = _estimate_builder_python_peak_bytes(
        limiting_spec.vector_length,
        count_dtype,
    )

    if estimated_python_peak_bytes + MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES > memory_limit_bytes:
        raise MemoryError(
            f"Scaffold {limiting_spec.chrom} is estimated to require about "
            f"{_format_memory_bytes(estimated_python_peak_bytes)} of Python-side working memory, "
            f"which leaves less than {_format_memory_bytes(MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES)} "
            f"for DuckDB under the total limit of {_format_memory_bytes(memory_limit_bytes)}. "
            "Increase --memory-limit-gb or narrow the genome scope."
        )

    available_for_duckdb = max(
        MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES,
        memory_limit_bytes - estimated_python_peak_bytes,
    )
    duckdb_memory_limit_bytes = max(
        MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES,
        min(
            available_for_duckdb,
            int(memory_limit_bytes * MATRIX_BUILD_DUCKDB_MEMORY_FRACTION),
        ),
    )

    if commit_batch_gb is None:
        commit_batch_bytes = max(
            MATRIX_BUILD_MIN_COMMIT_BATCH_BYTES,
            min(
                duckdb_memory_limit_bytes // 2,
                available_for_duckdb // 4,
            ),
        )
    else:
        commit_batch_bytes = _commit_batch_bytes(commit_batch_gb)
        commit_batch_bytes = min(
            commit_batch_bytes,
            max(MATRIX_BUILD_MIN_COMMIT_BATCH_BYTES, duckdb_memory_limit_bytes // 2),
        )

    return MatrixBuildMemoryPlan(
        duckdb_memory_limit_bytes=duckdb_memory_limit_bytes,
        commit_batch_bytes=commit_batch_bytes,
        estimated_python_peak_bytes=estimated_python_peak_bytes,
        limiting_scaffold=limiting_spec.chrom,
    )


def _init_matrix_db_schema(conn: duckdb.DuckDBPyConnection) -> None:
    conn.execute(
        """
        CREATE TABLE matrix_db_metadata (
          key VARCHAR PRIMARY KEY,
          value VARCHAR NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_db_samples (
          sample_idx INTEGER PRIMARY KEY,
          sample_name VARCHAR NOT NULL UNIQUE,
          profile_path VARCHAR NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_db_scaffolds (
          scaffold_idx INTEGER PRIMARY KEY,
          genome VARCHAR NOT NULL,
          chrom VARCHAR NOT NULL,
          index_base BIGINT NOT NULL,
          vector_length BIGINT NOT NULL,
          min_pos BIGINT NOT NULL,
          max_pos BIGINT NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_db_sample_scaffold_matrices (
          sample_idx INTEGER NOT NULL,
          scaffold_idx INTEGER NOT NULL,
          count_dtype VARCHAR NOT NULL,
          matrix_rows BIGINT NOT NULL,
          matrix_cols INTEGER NOT NULL,
          matrix_blob BLOB NOT NULL,
          PRIMARY KEY (sample_idx, scaffold_idx)
        )
        """
    )


def _open_matrix_build_connection(
    output_file: Path,
    duckdb_memory_limit_bytes: int,
) -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect(str(output_file))
    conn.execute("SET preserve_insertion_order=false")
    conn.execute(f"SET memory_limit='{_duckdb_memory_limit_setting(duckdb_memory_limit_bytes)}'")
    return conn


def _restart_matrix_build_transaction(
    conn: duckdb.DuckDBPyConnection,
    output_file: Path,
    duckdb_memory_limit_bytes: int,
) -> duckdb.DuckDBPyConnection:
    conn.execute("COMMIT")
    conn.close()
    gc.collect()
    conn = _open_matrix_build_connection(output_file, duckdb_memory_limit_bytes)
    conn.execute("BEGIN")
    return conn


def _load_profile_scaffold_matrix(
    profile_path: Path,
    spec: ScaffoldSpec,
    count_dtype: str,
    min_cov: int,
) -> Optional[np.ndarray]:
    np_dtype = COUNT_DTYPES[count_dtype]
    frame = (
        pl.scan_parquet(profile_path)
        .filter((pl.col("genome") == spec.genome) & (pl.col("chrom") == spec.chrom))
        .select("pos", "A", "T", "C", "G")
        .collect(engine="streaming")
    )
    if frame.height == 0:
        return None
    coverage = (
        frame["A"].cast(pl.Int64)
        + frame["T"].cast(pl.Int64)
        + frame["C"].cast(pl.Int64)
        + frame["G"].cast(pl.Int64)
    )
    keep_mask = coverage >= min_cov
    if not keep_mask.any():
        return None
    frame = frame.filter(keep_mask)
    matrix = np.zeros((spec.vector_length, 4), dtype=np_dtype)
    pos = frame["pos"].to_numpy().astype(np.int64) - spec.index_base
    matrix[pos, 0] = (frame["A"].to_numpy() > 0).astype(np_dtype, copy=False)
    matrix[pos, 1] = (frame["T"].to_numpy() > 0).astype(np_dtype, copy=False)
    matrix[pos, 2] = (frame["C"].to_numpy() > 0).astype(np_dtype, copy=False)
    matrix[pos, 3] = (frame["G"].to_numpy() > 0).astype(np_dtype, copy=False)
    return matrix


def build_matrix_db(
    profile_dir: Path,
    output_file: Path,
    genome: str = "all",
    count_dtype: str = "uint16",
    memory_limit_gb: float = 16.0,
    commit_batch_gb: Optional[float] = None,
    bed_file: Optional[Path] = None,
    progress_callback: Optional[BuildProgressCallback] = None,
) -> MatrixDbSummary:
    if count_dtype not in COUNT_DTYPES:
        raise ValueError(f"Unsupported count dtype '{count_dtype}'. Choose one of {', '.join(COUNT_DTYPES)}.")
    profile_paths = discover_profile_parquets(profile_dir)
    genome_scope = None if genome == "all" else genome
    if bed_file is None:
        scaffolds = _collect_scaffold_specs(profile_paths=profile_paths, genome=genome_scope)
    else:
        scaffolds = _collect_scaffold_specs_from_bed(
            profile_paths=profile_paths,
            bed_file=Path(bed_file),
            genome=genome_scope,
        )
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    memory_plan = _plan_matrix_build_memory(
        scaffolds=scaffolds,
        count_dtype=count_dtype,
        memory_limit_bytes=memory_limit_bytes,
        commit_batch_gb=commit_batch_gb,
    )
    commit_batch_bytes = memory_plan.commit_batch_bytes

    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    total_work = len(profile_paths) * len(scaffolds)
    completed_work = 0

    if progress_callback is not None:
        progress_callback(
            {
                "phase": "start",
                "completed": completed_work,
                "total": total_work,
                "sample_name": "",
                "genome": genome_scope or "all",
                "scaffold": "",
                "stored_rows": 0,
            }
        )

    conn = _open_matrix_build_connection(output_file, memory_plan.duckdb_memory_limit_bytes)
    build_succeeded = False
    try:
        _init_matrix_db_schema(conn)
        conn.execute("BEGIN")
        metadata_rows = [
            ("zipstrain_version", __version__),
            ("profiles_dir", str(profile_dir.resolve())),
            ("profile_format", "classic_zipstrain_profile_parquet"),
            ("genome_scope", genome_scope or "all"),
            ("count_dtype", count_dtype),
            ("layout", "per_sample_per_scaffold_dense_matrix"),
            ("matrix_value_semantics", FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS),
            ("coverage_filter_min_cov", str(MATRIX_BUILD_MIN_COV)),
            ("memory_limit_gb", str(memory_limit_gb)),
            ("duckdb_memory_limit_gb", f"{memory_plan.duckdb_memory_limit_bytes / (1024 ** 3):.6f}"),
            ("effective_commit_batch_gb", f"{commit_batch_bytes / (1024 ** 3):.6f}"),
            ("estimated_python_peak_gb", f"{memory_plan.estimated_python_peak_bytes / (1024 ** 3):.6f}"),
        ]
        sample_rows = [(idx, path.stem, str(path.resolve())) for idx, path in enumerate(profile_paths)]
        scaffold_rows = [
            (
                spec.scaffold_idx,
                spec.genome,
                spec.chrom,
                spec.index_base,
                spec.vector_length,
                spec.min_pos,
                spec.max_pos,
            )
            for spec in scaffolds
        ]
        conn.executemany("INSERT INTO matrix_db_metadata VALUES (?, ?)", metadata_rows)
        conn.executemany("INSERT INTO matrix_db_samples VALUES (?, ?, ?)", sample_rows)
        conn.executemany(
            "INSERT INTO matrix_db_scaffolds VALUES (?, ?, ?, ?, ?, ?, ?)",
            scaffold_rows,
        )
        conn = _restart_matrix_build_transaction(conn, output_file, memory_plan.duckdb_memory_limit_bytes)

        stored_rows = 0
        batch_bytes = 0
        for sample_idx, profile_path in enumerate(profile_paths):
            sample_name = profile_path.stem
            for spec in scaffolds:
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "processing",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": spec.genome,
                            "scaffold": spec.chrom,
                            "stored_rows": stored_rows,
                        }
                    )
                estimated_python_peak = _estimate_builder_python_peak_bytes(spec.vector_length, count_dtype)
                if estimated_python_peak + memory_plan.duckdb_memory_limit_bytes > memory_limit_bytes:
                    raise MemoryError(
                        f"Scaffold {spec.chrom} for sample {profile_path.name} is estimated to use about "
                        f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory plus "
                        f"{_format_memory_bytes(memory_plan.duckdb_memory_limit_bytes)} for DuckDB, "
                        f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                    )
                matrix = _load_profile_scaffold_matrix(
                    profile_path=profile_path,
                    spec=spec,
                    count_dtype=count_dtype,
                    min_cov=MATRIX_BUILD_MIN_COV,
                )
                if matrix is None:
                    completed_work += 1
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "advance",
                                "completed": completed_work,
                                "total": total_work,
                                "sample_name": sample_name,
                                "genome": spec.genome,
                                "scaffold": spec.chrom,
                                "stored_rows": stored_rows,
                            }
                        )
                    continue
                matrix_blob = _pack_matrix(matrix)
                conn.execute(
                    "INSERT INTO matrix_db_sample_scaffold_matrices VALUES (?, ?, ?, ?, ?, ?)",
                    [
                        sample_idx,
                        spec.scaffold_idx,
                        count_dtype,
                        matrix.shape[0],
                        matrix.shape[1],
                        matrix_blob,
                    ],
                )
                batch_bytes += len(matrix_blob)
                del matrix_blob
                del matrix
                stored_rows += 1
                completed_work += 1
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "advance",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": spec.genome,
                            "scaffold": spec.chrom,
                            "stored_rows": stored_rows,
                        }
                    )
                if batch_bytes >= commit_batch_bytes:
                    conn = _restart_matrix_build_transaction(
                        conn,
                        output_file,
                        memory_plan.duckdb_memory_limit_bytes,
                    )
                    batch_bytes = 0
            if batch_bytes > 0:
                conn = _restart_matrix_build_transaction(
                    conn,
                    output_file,
                    memory_plan.duckdb_memory_limit_bytes,
                )
                batch_bytes = 0
        conn.execute("COMMIT")
        build_succeeded = True
    except Exception:
        try:
            conn.execute("ROLLBACK")
        except Exception:
            pass
        raise
    finally:
        conn.close()
        if not build_succeeded and output_file.exists():
            output_file.unlink(missing_ok=True)

    if progress_callback is not None:
        progress_callback(
            {
                "phase": "done",
                "completed": completed_work,
                "total": total_work,
                "sample_name": "",
                "genome": genome_scope or "all",
                "scaffold": "",
                "stored_rows": stored_rows,
            }
        )

    return MatrixDbSummary(
        output_file=output_file,
        profile_files=len(profile_paths),
        sample_count=len(profile_paths),
        scaffold_count=len(scaffolds),
        stored_rows=stored_rows,
        count_dtype=count_dtype,
        genome_scope=genome_scope or "all",
        memory_limit_gb=memory_limit_gb,
    )


def _load_matrix_db_metadata(conn: duckdb.DuckDBPyConnection) -> dict[str, str]:
    return {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM matrix_db_metadata").fetchall()}


def _load_matrix_db_samples(conn: duckdb.DuckDBPyConnection) -> list[tuple[int, str]]:
    return [
        (int(sample_idx), str(sample_name))
        for sample_idx, sample_name in conn.execute(
            "SELECT sample_idx, sample_name FROM matrix_db_samples ORDER BY sample_idx"
        ).fetchall()
    ]


def _load_matrix_db_scaffolds(conn: duckdb.DuckDBPyConnection, genome: Optional[str] = None) -> list[ScaffoldSpec]:
    if genome is None:
        rows = conn.execute(
            """
            SELECT scaffold_idx, genome, chrom, index_base, vector_length, min_pos, max_pos
            FROM matrix_db_scaffolds
            ORDER BY scaffold_idx
            """
        ).fetchall()
    else:
        rows = conn.execute(
            """
            SELECT scaffold_idx, genome, chrom, index_base, vector_length, min_pos, max_pos
            FROM matrix_db_scaffolds
            WHERE genome = ?
            ORDER BY scaffold_idx
            """,
            [genome],
        ).fetchall()
    return [
        ScaffoldSpec(
            scaffold_idx=int(scaffold_idx),
            genome=str(genome_name),
            chrom=str(chrom),
            index_base=int(index_base),
            vector_length=int(vector_length),
            min_pos=int(min_pos),
            max_pos=int(max_pos),
        )
        for scaffold_idx, genome_name, chrom, index_base, vector_length, min_pos, max_pos in rows
    ]


def _load_sample_scaffold_matrices(
    conn: duckdb.DuckDBPyConnection,
    scaffold_idx: int,
    sample_ids: list[int],
    vector_length: int,
) -> tuple[str, dict[int, np.ndarray]]:
    if not sample_ids:
        raise ValueError("sample_ids must not be empty")
    rows = conn.execute(
        f"""
        SELECT sample_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob
        FROM matrix_db_sample_scaffold_matrices
        WHERE scaffold_idx = ? AND sample_idx IN ({','.join(['?'] * len(sample_ids))})
        ORDER BY sample_idx
        """,
        [scaffold_idx, *sample_ids],
    ).fetchall()
    dtype_name = ""
    matrices: dict[int, np.ndarray] = {}
    for sample_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in rows:
        dtype_name = str(count_dtype)
        matrices[int(sample_idx)] = _unpack_matrix(
            bytes(matrix_blob),
            dtype_name,
            (int(matrix_rows), int(matrix_cols)),
        )
    if dtype_name and any(matrix.shape[0] != vector_length for matrix in matrices.values()):
        raise ValueError(f"Stored matrix rows do not match scaffold vector length for scaffold_idx={scaffold_idx}")
    return dtype_name, matrices


def _load_scaffold_count_dtype(
    conn: duckdb.DuckDBPyConnection,
    scaffold_idx: int,
) -> str:
    row = conn.execute(
        """
        SELECT count_dtype
        FROM matrix_db_sample_scaffold_matrices
        WHERE scaffold_idx = ?
        LIMIT 1
        """,
        [scaffold_idx],
    ).fetchone()
    return "" if row is None else str(row[0])


def _plan_chunk_sizes(
    vector_length: int,
    remaining_targets: int,
    dtype_name: str,
    memory_limit_bytes: int,
    backend_kind: str,
    position_tile_size: Optional[int] = None,
) -> tuple[int, int]:
    dtype_bytes = np.dtype(COUNT_DTYPES[dtype_name]).itemsize
    reserve = int(memory_limit_bytes * 0.15)
    budget = max(memory_limit_bytes - reserve, 64 * 1024 * 1024)

    per_position_anchor = 4 * dtype_bytes + 8
    per_position_target = 4 * dtype_bytes + 8
    if backend_kind == "torch":
        per_position_anchor += 4 * dtype_bytes
        per_position_target += 4 * dtype_bytes

    if position_tile_size is not None:
        tile_size = min(vector_length, position_tile_size)
        max_targets = max(
            1,
            int((budget / max(tile_size, 1) - per_position_anchor) // max(per_position_target, 1)),
        )
        return min(remaining_targets, max_targets), tile_size

    max_targets_full = max(
        1,
        int((budget / max(vector_length, 1) - per_position_anchor) // max(per_position_target, 1)),
    )
    if max_targets_full >= 1:
        return min(remaining_targets, max_targets_full), vector_length

    tile_size = max(1, int(budget // max(per_position_anchor + per_position_target, 1)))
    return 1, min(vector_length, tile_size)


def _compare_tile_presence_numpy(
    anchor_matrix: np.ndarray,
    target_matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    anchor_presence = anchor_matrix.astype(np.int8, copy=False)
    target_presence = target_matrices.astype(np.int8, copy=False)

    anchor_cov = anchor_presence.max(axis=1).astype(np.int32, copy=False)
    target_cov = target_presence.max(axis=1).astype(np.int32, copy=False)
    totals = (anchor_cov @ target_cov).astype(np.int64, copy=False)

    shared_signal = np.matmul(anchor_presence[:, np.newaxis, :], target_presence).squeeze(1) > 0
    shared = shared_signal.sum(axis=0, dtype=np.int64)
    return totals, shared


def _compare_tile_presence_torch_tensors(
    torch_module,
    anchor_t,
    targets_t,
):
    anchor_cov = anchor_t.amax(dim=1)
    target_cov = targets_t.amax(dim=1)
    totals = torch_module.matmul(anchor_cov.unsqueeze(0), target_cov).squeeze(0).to(torch_module.int64)

    shared_scores = torch_module.matmul(anchor_t.unsqueeze(1), targets_t).squeeze(1)
    shared = (shared_scores > 0).sum(dim=0, dtype=torch_module.int64)
    return totals, shared


def _compare_tile_presence_numpy_with_mask(
    anchor_matrix: np.ndarray,
    target_matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    anchor_presence = anchor_matrix.astype(np.int8, copy=False)
    target_presence = target_matrices.astype(np.int8, copy=False)

    anchor_cov = anchor_presence.max(axis=1).astype(np.int32, copy=False)
    target_cov = target_presence.max(axis=1).astype(np.int32, copy=False)
    totals = (anchor_cov @ target_cov).astype(np.int64, copy=False)

    shared_mask = np.matmul(anchor_presence[:, np.newaxis, :], target_presence).squeeze(1) > 0
    shared = shared_mask.sum(axis=0, dtype=np.int64)
    return totals, shared, shared_mask


def _shared_mask_presence_torch_tensors(
    torch_module,
    anchor_t,
    targets_t,
):
    shared_scores = torch_module.matmul(anchor_t.unsqueeze(1), targets_t).squeeze(1)
    return shared_scores > 0


def _compare_tile_presence_torch_tensors_with_mask(
    torch_module,
    anchor_t,
    targets_t,
):
    anchor_cov = anchor_t.amax(dim=1)
    target_cov = targets_t.amax(dim=1)
    totals = torch_module.matmul(anchor_cov.unsqueeze(0), target_cov).squeeze(0).to(torch_module.int64)
    shared_mask = _shared_mask_presence_torch_tensors(
        torch_module=torch_module,
        anchor_t=anchor_t,
        targets_t=targets_t,
    )
    shared = shared_mask.sum(dim=0, dtype=torch_module.int64)
    return totals, shared, shared_mask


def _update_ibs_numpy(
    shared_mask: np.ndarray,
    current_runs: np.ndarray,
    max_runs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if shared_mask.size == 0:
        return current_runs, max_runs
    next_runs = np.zeros_like(current_runs)
    updated_max = max_runs.copy()
    for col_idx in range(shared_mask.shape[1]):
        col_mask = shared_mask[:, col_idx]
        if not bool(col_mask.any()):
            continue
        padded = np.concatenate(([False], col_mask.astype(bool, copy=False), [False]))
        diff = np.diff(padded.astype(np.int8, copy=False))
        starts = np.flatnonzero(diff == 1)
        ends = np.flatnonzero(diff == -1)
        lengths = (ends - starts).astype(np.int64, copy=False)
        if col_mask[0]:
            lengths[0] += current_runs[col_idx]
        if lengths.size:
            updated_max[col_idx] = max(updated_max[col_idx], int(lengths.max()))
        if col_mask[-1]:
            next_runs[col_idx] = lengths[-1]
    return next_runs, updated_max


def _shared_mask_to_numpy(shared_mask) -> np.ndarray:
    if isinstance(shared_mask, np.ndarray):
        return shared_mask
    if hasattr(shared_mask, "detach"):
        return shared_mask.detach().cpu().numpy()
    if hasattr(shared_mask, "cpu"):
        return shared_mask.cpu().numpy()
    return np.asarray(shared_mask)


def _make_arrow_table(
    sample_1: str,
    sample_2: list[str],
    genome: str,
    calculations: tuple[str, ...],
    total_positions: Optional[np.ndarray] = None,
    share_allele_pos: Optional[np.ndarray] = None,
    max_consecutive_length: Optional[np.ndarray] = None,
) -> pa.Table:
    arrays = [
        pa.array([sample_1] * len(sample_2), type=pa.string()),
        pa.array(sample_2, type=pa.string()),
        pa.array([genome] * len(sample_2), type=pa.string()),
    ]
    if "ani" in calculations:
        if total_positions is None or share_allele_pos is None:
            raise ValueError("ANI output requested but total_positions/share_allele_pos were not provided.")
        arrays.extend(
            [
                pa.array(total_positions.tolist(), type=pa.int64()),
                pa.array(share_allele_pos.tolist(), type=pa.int64()),
                pa.array((share_allele_pos / total_positions * 100.0).tolist(), type=pa.float64()),
            ]
        )
    if "ibs" in calculations:
        if max_consecutive_length is None:
            raise ValueError("IBS output requested but max_consecutive_length was not provided.")
        arrays.append(pa.array(max_consecutive_length.tolist(), type=pa.int64()))
    return pa.Table.from_arrays(
        arrays,
        schema=matrix_pair_output_schema(calculations),
    )


def _prepare_torch_matrix(
    compute_backend: MatrixPairComputeBackend,
    matrix: np.ndarray,
    matrix_value_semantics: str,
):
    tensor = compute_backend.torch.from_numpy(np.ascontiguousarray(matrix))
    if compute_backend.device == "cuda" and hasattr(tensor, "pin_memory"):
        tensor = tensor.pin_memory()
        kwargs: dict[str, object] = {"device": compute_backend.device, "non_blocking": True}
    else:
        kwargs = {"device": compute_backend.device}
    if matrix_value_semantics == FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
        kwargs["dtype"] = compute_backend.torch.float32
    return tensor.to(**kwargs)


def _compare_anchor_against_target_chunk_torch(
    compute_backend: MatrixPairComputeBackend,
    anchor_torch,
    target_torch,
    vector_length: int,
    tile_size: int,
    matrix_value_semantics: str,
    min_cov: int,
    need_ibs: bool = False,
) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    target_count = int(target_torch.shape[2])
    chunk_totals_torch = compute_backend.torch.zeros(
        target_count,
        dtype=compute_backend.torch.int64,
        device=compute_backend.device,
    )
    chunk_shared_torch = compute_backend.torch.zeros(
        target_count,
        dtype=compute_backend.torch.int64,
        device=compute_backend.device,
    )
    current_runs = None
    max_runs = None
    if need_ibs:
        current_runs = np.zeros(target_count, dtype=np.int64)
        max_runs = np.zeros(target_count, dtype=np.int64)
    for tile_start in range(0, vector_length, tile_size):
        tile_stop = min(vector_length, tile_start + tile_size)
        if need_ibs:
            total_inc, shared_inc, shared_mask = _compare_tile_presence_torch_tensors_with_mask(
                torch_module=compute_backend.torch,
                anchor_t=anchor_torch[tile_start:tile_stop, :],
                targets_t=target_torch[tile_start:tile_stop, :, :],
            )
        else:
            total_inc, shared_inc = _compare_tile_presence_torch_tensors(
                torch_module=compute_backend.torch,
                anchor_t=anchor_torch[tile_start:tile_stop, :],
                targets_t=target_torch[tile_start:tile_stop, :, :],
            )
        chunk_totals_torch += total_inc
        chunk_shared_torch += shared_inc
        if need_ibs:
            current_runs, max_runs = _update_ibs_numpy(
                shared_mask=_shared_mask_to_numpy(shared_mask),
                current_runs=current_runs,
                max_runs=max_runs,
            )
    return chunk_totals_torch.cpu().numpy(), chunk_shared_torch.cpu().numpy(), max_runs


def _matrix_compare_reuse_target_chunks_torch(
    conn: duckdb.DuckDBPyConnection,
    output_file: Path,
    min_cov: int,
    genome_scope: Optional[str],
    metadata: dict[str, str],
    samples: list[tuple[int, str]],
    scaffolds: list[ScaffoldSpec],
    requested_pairs: int,
    memory_limit_bytes: int,
    memory_limit_gb: float,
    position_tile_size: Optional[int],
    compute_backend: MatrixPairComputeBackend,
    matrix_value_semantics: str,
    calculations: tuple[str, ...],
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    total_work = requested_pairs * len(scaffolds)
    completed_work = 0
    written_rows = 0
    target_chunks = 0
    anchor_groups = max(len(samples) - 1, 0)
    dtype_name = str(metadata.get("count_dtype", "uint16"))
    default_tile_size = 0
    max_vector_length = max(spec.vector_length for spec in scaffolds)
    global_block_size, _ = _plan_chunk_sizes(
        vector_length=max_vector_length,
        remaining_targets=len(samples),
        dtype_name=dtype_name,
        memory_limit_bytes=memory_limit_bytes,
        backend_kind=compute_backend.kind,
        position_tile_size=position_tile_size,
    )
    global_block_size = max(1, global_block_size)
    samples_by_id = {sample_idx: sample_name for sample_idx, sample_name in samples}
    scaffolds_by_genome: dict[str, list[ScaffoldSpec]] = {}
    for spec in scaffolds:
        scaffolds_by_genome.setdefault(spec.genome, []).append(spec)

    if progress_callback is not None:
        progress_callback(
            {
                "phase": "start",
                "completed": completed_work,
                "total": total_work,
                "anchor_name": "",
                "genome": genome_scope or metadata.get("genome_scope", "all"),
                "scaffold": "",
                "targets_completed": 0,
                "targets_total": 0,
                "target_chunks": target_chunks,
            }
        )

    output_schema = matrix_pair_output_schema(calculations)
    with pq.ParquetWriter(output_file, output_schema, compression="zstd") as writer:
        for genome_name in sorted(scaffolds_by_genome):
            genome_scaffolds = scaffolds_by_genome[genome_name]
            for block_start in range(0, len(samples), global_block_size):
                block_rows = samples[block_start:block_start + global_block_size]
                if not block_rows:
                    continue
                external_rows = samples[:block_start]
                pair_count_for_block = len(external_rows) * len(block_rows) + (len(block_rows) * (len(block_rows) - 1)) // 2
                if pair_count_for_block == 0:
                    continue

                block_ids = np.array([sample_idx for sample_idx, _sample_name in block_rows], dtype=np.int64)
                block_names = [sample_name for _sample_idx, sample_name in block_rows]
                totals_external = np.zeros((len(external_rows), len(block_rows)), dtype=np.int64)
                shared_external = np.zeros((len(external_rows), len(block_rows)), dtype=np.int64)
                totals_internal = np.zeros((len(block_rows), len(block_rows)), dtype=np.int64)
                shared_internal = np.zeros((len(block_rows), len(block_rows)), dtype=np.int64)
                ibs_external = np.zeros((len(external_rows), len(block_rows)), dtype=np.int64) if "ibs" in calculations else None
                ibs_internal = np.zeros((len(block_rows), len(block_rows)), dtype=np.int64) if "ibs" in calculations else None

                for spec in genome_scaffolds:
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "processing",
                                "completed": completed_work,
                                "total": total_work,
                                "anchor_name": block_names[0],
                                "genome": spec.genome,
                                "scaffold": spec.chrom,
                                "targets_completed": 0,
                                "targets_total": len(block_ids),
                                "target_chunks": target_chunks,
                            }
                        )
                    scaffold_dtype = _load_scaffold_count_dtype(conn, spec.scaffold_idx)
                    if not scaffold_dtype:
                        completed_work += pair_count_for_block
                        target_chunks += 1
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": block_names[0],
                                    "genome": spec.genome,
                                    "scaffold": spec.chrom,
                                    "targets_completed": len(block_ids),
                                    "targets_total": len(block_ids),
                                    "target_chunks": target_chunks,
                                }
                            )
                        continue

                    zero_matrix = np.zeros((spec.vector_length, 4), dtype=COUNT_DTYPES[dtype_name])
                    tile_targets, tile_size = _plan_chunk_sizes(
                        vector_length=spec.vector_length,
                        remaining_targets=len(block_ids),
                        dtype_name=dtype_name,
                        memory_limit_bytes=memory_limit_bytes,
                        backend_kind=compute_backend.kind,
                        position_tile_size=position_tile_size,
                    )
                    default_tile_size = tile_size
                    if tile_targets < len(block_ids):
                        raise RuntimeError(
                            "Internal target block size exceeded the planned torch chunk capacity. "
                            "This indicates the global block planner is inconsistent."
                        )
                    _dtype, loaded_targets = _load_sample_scaffold_matrices(
                        conn,
                        scaffold_idx=spec.scaffold_idx,
                        sample_ids=block_ids.tolist(),
                        vector_length=spec.vector_length,
                    )
                    target_matrices = np.stack(
                        [loaded_targets.get(int(sample_idx), zero_matrix) for sample_idx in block_ids],
                        axis=2,
                    )
                    target_torch = _prepare_torch_matrix(
                        compute_backend=compute_backend,
                        matrix=target_matrices,
                        matrix_value_semantics=matrix_value_semantics,
                    )
                    processed_pairs_for_block = 0

                    for external_pos, (anchor_idx, _anchor_name) in enumerate(external_rows):
                        _anchor_dtype, loaded_anchor = _load_sample_scaffold_matrices(
                            conn,
                            scaffold_idx=spec.scaffold_idx,
                            sample_ids=[anchor_idx],
                            vector_length=spec.vector_length,
                        )
                        anchor_matrix = loaded_anchor.get(anchor_idx, zero_matrix)
                        anchor_torch = _prepare_torch_matrix(
                            compute_backend=compute_backend,
                            matrix=anchor_matrix,
                            matrix_value_semantics=matrix_value_semantics,
                        )
                        total_inc, shared_inc, ibs_inc = _compare_anchor_against_target_chunk_torch(
                            compute_backend=compute_backend,
                            anchor_torch=anchor_torch,
                            target_torch=target_torch,
                            vector_length=spec.vector_length,
                            tile_size=tile_size,
                            matrix_value_semantics=matrix_value_semantics,
                            min_cov=min_cov,
                            need_ibs="ibs" in calculations,
                        )
                        totals_external[external_pos, :] += total_inc
                        shared_external[external_pos, :] += shared_inc
                        if ibs_external is not None and ibs_inc is not None:
                            ibs_external[external_pos, :] = np.maximum(ibs_external[external_pos, :], ibs_inc)
                        processed_pairs_for_block += len(block_ids)
                        completed_work += len(block_ids)
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": samples_by_id[anchor_idx],
                                    "genome": spec.genome,
                                    "scaffold": spec.chrom,
                                    "targets_completed": external_pos + 1,
                                    "targets_total": len(external_rows) + len(block_ids) - 1,
                                    "target_chunks": target_chunks,
                                }
                            )

                    for local_anchor_pos in range(len(block_ids) - 1):
                        total_inc, shared_inc, ibs_inc = _compare_anchor_against_target_chunk_torch(
                            compute_backend=compute_backend,
                            anchor_torch=target_torch[:, :, local_anchor_pos],
                            target_torch=target_torch[:, :, local_anchor_pos + 1:],
                            vector_length=spec.vector_length,
                            tile_size=tile_size,
                            matrix_value_semantics=matrix_value_semantics,
                            min_cov=min_cov,
                            need_ibs="ibs" in calculations,
                        )
                        totals_internal[local_anchor_pos, local_anchor_pos + 1:] += total_inc
                        shared_internal[local_anchor_pos, local_anchor_pos + 1:] += shared_inc
                        if ibs_internal is not None and ibs_inc is not None:
                            ibs_internal[local_anchor_pos, local_anchor_pos + 1:] = np.maximum(
                                ibs_internal[local_anchor_pos, local_anchor_pos + 1:],
                                ibs_inc,
                            )
                        remaining_targets = len(block_ids) - local_anchor_pos - 1
                        processed_pairs_for_block += remaining_targets
                        completed_work += remaining_targets
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": block_names[local_anchor_pos],
                                    "genome": spec.genome,
                                    "scaffold": spec.chrom,
                                    "targets_completed": local_anchor_pos + 1,
                                    "targets_total": len(block_ids) - 1,
                                    "target_chunks": target_chunks,
                                }
                            )

                    if processed_pairs_for_block != pair_count_for_block:
                        raise RuntimeError(
                            "Torch block compare progress accounting drifted from the expected pair count."
                        )
                    target_chunks += 1
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "advance",
                                "completed": completed_work,
                                "total": total_work,
                                "anchor_name": block_names[0],
                                "genome": spec.genome,
                                "scaffold": spec.chrom,
                                "targets_completed": len(block_ids),
                                "targets_total": len(block_ids),
                                "target_chunks": target_chunks,
                            }
                        )

                for external_pos, (anchor_idx, _anchor_name) in enumerate(external_rows):
                    totals = totals_external[external_pos]
                    shared = shared_external[external_pos]
                    mask = totals > 0
                    if not mask.any():
                        continue
                    table = _make_arrow_table(
                        sample_1=samples_by_id[anchor_idx],
                        sample_2=[block_names[idx] for idx, keep in enumerate(mask) if keep],
                        genome=genome_name,
                        calculations=calculations,
                        total_positions=totals[mask] if "ani" in calculations else None,
                        share_allele_pos=shared[mask] if "ani" in calculations else None,
                        max_consecutive_length=ibs_external[external_pos, :][mask] if ibs_external is not None else None,
                    )
                    writer.write_table(table)
                    written_rows += table.num_rows

                for local_anchor_pos, anchor_name in enumerate(block_names[:-1]):
                    totals = totals_internal[local_anchor_pos, local_anchor_pos + 1:]
                    shared = shared_internal[local_anchor_pos, local_anchor_pos + 1:]
                    mask = totals > 0
                    if not mask.any():
                        continue
                    later_names = block_names[local_anchor_pos + 1:]
                    table = _make_arrow_table(
                        sample_1=anchor_name,
                        sample_2=[later_names[idx] for idx, keep in enumerate(mask) if keep],
                        genome=genome_name,
                        calculations=calculations,
                        total_positions=totals[mask] if "ani" in calculations else None,
                        share_allele_pos=shared[mask] if "ani" in calculations else None,
                        max_consecutive_length=ibs_internal[local_anchor_pos, local_anchor_pos + 1:][mask] if ibs_internal is not None else None,
                    )
                    writer.write_table(table)
                    written_rows += table.num_rows

    if progress_callback is not None:
        progress_callback(
            {
                "phase": "done",
                "completed": completed_work,
                "total": total_work,
                "anchor_name": "",
                "genome": genome_scope or metadata.get("genome_scope", "all"),
                "scaffold": "",
                "targets_completed": 0,
                "targets_total": 0,
                "target_chunks": target_chunks,
            }
        )

    return MatrixCompareSummary(
        output_file=output_file,
        requested_pairs=requested_pairs,
        written_rows=written_rows,
        scaffold_count=len(scaffolds),
        genome_count=len(scaffolds_by_genome),
        anchor_groups=anchor_groups,
        target_chunks=target_chunks,
        min_cov=min_cov,
        genome_scope=genome_scope or metadata.get("genome_scope", "all"),
        backend=compute_backend.kind,
        device=compute_backend.device,
        memory_limit_gb=memory_limit_gb,
        position_tile_size=position_tile_size or default_tile_size,
    )


def matrix_compare(
    matrix_db_file: Path,
    output_file: Path,
    min_cov: int = 5,
    genome: str = "all",
    memory_limit_gb: float = 16.0,
    position_tile_size: Optional[int] = None,
    backend: str = "numpy",
    calculate: Optional[str] = "all",
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    if min_cov < 1:
        raise ValueError("min_cov must be >= 1")
    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    compute_backend = MatrixPairComputeBackend(backend)
    calculations = parse_matrix_calculations(calculate)
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    conn = duckdb.connect(str(matrix_db_file), read_only=True)
    try:
        metadata = _load_matrix_db_metadata(conn)
        matrix_value_semantics = metadata.get("matrix_value_semantics")
        if matrix_value_semantics != FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
            raise ValueError(
                "This matrix DB uses an unsupported legacy storage layout. "
                "Rebuild it with 'zipstrain utilities build-matrix-db'."
            )
        filtered_min_cov = metadata.get("coverage_filter_min_cov")
        required_min_cov = int(filtered_min_cov or MATRIX_BUILD_MIN_COV)
        if min_cov != required_min_cov:
            raise ValueError(
                f"Matrix DB was built with fixed coverage filter min_cov={required_min_cov}. "
                f"Requested min_cov={min_cov}. Rebuild the matrix DB or use --min-cov {required_min_cov}."
            )
        samples = _load_matrix_db_samples(conn)
        requested_pairs = len(samples) * (len(samples) - 1) // 2
        if requested_pairs == 0:
            if progress_callback is not None:
                progress_callback(
                    {
                        "phase": "start",
                        "completed": 0,
                        "total": 0,
                        "anchor_name": "",
                        "genome": genome,
                        "scaffold": "",
                        "targets_completed": 0,
                        "targets_total": 0,
                        "target_chunks": 0,
                    }
                )
            output_schema = matrix_pair_output_schema(calculations)
            with pq.ParquetWriter(output_file, output_schema, compression="zstd") as writer:
                writer.write_table(pa.Table.from_arrays([pa.array([], type=f.type) for f in output_schema], schema=output_schema))
            if progress_callback is not None:
                progress_callback(
                    {
                        "phase": "done",
                        "completed": 0,
                        "total": 0,
                        "anchor_name": "",
                        "genome": genome,
                        "scaffold": "",
                        "targets_completed": 0,
                        "targets_total": 0,
                        "target_chunks": 0,
                    }
                )
            return MatrixCompareSummary(
                output_file=output_file,
                requested_pairs=0,
                written_rows=0,
                scaffold_count=0,
                genome_count=0,
                anchor_groups=0,
                target_chunks=0,
                min_cov=min_cov,
                genome_scope=genome,
                backend=compute_backend.kind,
                device=compute_backend.device,
                memory_limit_gb=memory_limit_gb,
                position_tile_size=position_tile_size or 0,
            )

        genome_scope = None if genome == "all" else genome
        scaffolds = _load_matrix_db_scaffolds(conn, genome=genome_scope)
        if not scaffolds:
            raise ValueError(f"No scaffolds found for genome scope: {genome}")
        if compute_backend.kind == "torch":
            return _matrix_compare_reuse_target_chunks_torch(
                conn=conn,
                output_file=output_file,
                min_cov=min_cov,
                genome_scope=genome_scope,
                metadata=metadata,
                samples=samples,
                scaffolds=scaffolds,
                requested_pairs=requested_pairs,
                memory_limit_bytes=memory_limit_bytes,
                memory_limit_gb=memory_limit_gb,
                position_tile_size=position_tile_size,
                compute_backend=compute_backend,
                matrix_value_semantics=matrix_value_semantics,
                calculations=calculations,
                progress_callback=progress_callback,
            )

        total_work = requested_pairs * len(scaffolds)
        completed_work = 0
        written_rows = 0
        target_chunks = 0
        default_tile_size = 0
        anchor_groups = 0

        if progress_callback is not None:
            progress_callback(
                {
                    "phase": "start",
                    "completed": completed_work,
                    "total": total_work,
                    "anchor_name": "",
                    "genome": genome_scope or "all",
                    "scaffold": "",
                    "targets_completed": 0,
                    "targets_total": 0,
                    "target_chunks": target_chunks,
                }
            )

        output_schema = matrix_pair_output_schema(calculations)
        with pq.ParquetWriter(output_file, output_schema, compression="zstd") as writer:
            for anchor_offset, (sample_1_idx, sample_1_name) in enumerate(samples[:-1]):
                target_sample_rows = samples[anchor_offset + 1 :]
                if not target_sample_rows:
                    continue
                anchor_groups += 1
                target_ids_all = np.array([sample_idx for sample_idx, _sample_name in target_sample_rows], dtype=np.int64)
                target_names_all = [sample_name for _sample_idx, sample_name in target_sample_rows]
                totals_by_genome: dict[str, np.ndarray] = {
                    spec.genome: np.zeros(len(target_ids_all), dtype=np.int64) for spec in scaffolds
                }
                shared_by_genome: dict[str, np.ndarray] = {
                    spec.genome: np.zeros(len(target_ids_all), dtype=np.int64) for spec in scaffolds
                }
                ibs_by_genome: dict[str, np.ndarray] = (
                    {spec.genome: np.zeros(len(target_ids_all), dtype=np.int64) for spec in scaffolds}
                    if "ibs" in calculations
                    else {}
                )

                for spec in scaffolds:
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "processing",
                                "completed": completed_work,
                                "total": total_work,
                                "anchor_name": sample_1_name,
                                "genome": spec.genome,
                                "scaffold": spec.chrom,
                                "targets_completed": 0,
                                "targets_total": len(target_ids_all),
                                "target_chunks": target_chunks,
                            }
                        )
                    dtype_name = _load_scaffold_count_dtype(conn, spec.scaffold_idx)
                    if not dtype_name:
                        completed_work += len(target_ids_all)
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": sample_1_name,
                                    "genome": spec.genome,
                                    "scaffold": spec.chrom,
                                    "targets_completed": len(target_ids_all),
                                    "targets_total": len(target_ids_all),
                                    "target_chunks": target_chunks,
                                }
                            )
                        continue
                    zero_matrix = np.zeros((spec.vector_length, 4), dtype=COUNT_DTYPES[dtype_name])
                    _anchor_dtype, anchor_loaded = _load_sample_scaffold_matrices(
                        conn,
                        scaffold_idx=spec.scaffold_idx,
                        sample_ids=[sample_1_idx],
                        vector_length=spec.vector_length,
                    )
                    anchor_matrix = anchor_loaded.get(sample_1_idx, zero_matrix)
                    anchor_torch = None
                    if compute_backend.kind == "torch":
                        anchor_array = np.ascontiguousarray(anchor_matrix)
                        anchor_torch = compute_backend.torch.from_numpy(anchor_array).to(
                            device=compute_backend.device,
                            dtype=compute_backend.torch.float32,
                        )

                    target_offset = 0
                    while target_offset < len(target_ids_all):
                        max_targets, tile_size = _plan_chunk_sizes(
                            vector_length=spec.vector_length,
                            remaining_targets=len(target_ids_all) - target_offset,
                            dtype_name=dtype_name,
                            memory_limit_bytes=memory_limit_bytes,
                            backend_kind=compute_backend.kind,
                            position_tile_size=position_tile_size,
                        )
                        default_tile_size = tile_size
                        chunk_ids = target_ids_all[target_offset: target_offset + max_targets]
                        _chunk_dtype, loaded_targets = _load_sample_scaffold_matrices(
                            conn,
                            scaffold_idx=spec.scaffold_idx,
                            sample_ids=chunk_ids.tolist(),
                            vector_length=spec.vector_length,
                        )
                        target_matrices = np.stack(
                            [loaded_targets.get(int(sample_idx), zero_matrix) for sample_idx in chunk_ids],
                            axis=2,
                        )
                        totals = totals_by_genome[spec.genome][target_offset: target_offset + max_targets]
                        shared = shared_by_genome[spec.genome][target_offset: target_offset + max_targets]
                        ibs_max = ibs_by_genome[spec.genome][target_offset: target_offset + max_targets] if "ibs" in calculations else None
                        current_runs = np.zeros(len(chunk_ids), dtype=np.int64) if "ibs" in calculations else None
                        for tile_start in range(0, spec.vector_length, tile_size):
                            tile_stop = min(spec.vector_length, tile_start + tile_size)
                            anchor_tile = anchor_matrix[tile_start:tile_stop, :]
                            target_tile = target_matrices[tile_start:tile_stop, :, :]
                            if "ibs" in calculations:
                                total_inc, shared_inc, shared_mask = _compare_tile_presence_numpy_with_mask(
                                    anchor_matrix=anchor_tile,
                                    target_matrices=target_tile,
                                )
                            else:
                                total_inc, shared_inc = _compare_tile_presence_numpy(
                                    anchor_matrix=anchor_tile,
                                    target_matrices=target_tile,
                                )
                            totals += total_inc
                            shared += shared_inc
                            if "ibs" in calculations and current_runs is not None and ibs_max is not None:
                                current_runs, updated_max = _update_ibs_numpy(
                                    shared_mask=shared_mask,
                                    current_runs=current_runs,
                                    max_runs=ibs_max,
                                )
                                ibs_max[:] = updated_max
                        target_offset += max_targets
                        completed_work += max_targets
                        target_chunks += 1
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": sample_1_name,
                                    "genome": spec.genome,
                                    "scaffold": spec.chrom,
                                    "targets_completed": target_offset,
                                    "targets_total": len(target_ids_all),
                                    "target_chunks": target_chunks,
                                }
                            )

                for genome_name in sorted(totals_by_genome):
                    totals = totals_by_genome[genome_name]
                    shared = shared_by_genome[genome_name]
                    mask = totals > 0
                    if not mask.any():
                        continue
                    table = _make_arrow_table(
                        sample_1=sample_1_name,
                        sample_2=[name for idx, name in enumerate(target_names_all) if mask[idx]],
                        genome=genome_name,
                        calculations=calculations,
                        total_positions=totals[mask] if "ani" in calculations else None,
                        share_allele_pos=shared[mask] if "ani" in calculations else None,
                        max_consecutive_length=ibs_by_genome[genome_name][mask] if "ibs" in calculations else None,
                    )
                    writer.write_table(table)
                    written_rows += table.num_rows

        if progress_callback is not None:
            progress_callback(
                {
                    "phase": "done",
                    "completed": completed_work,
                    "total": total_work,
                    "anchor_name": "",
                    "genome": genome_scope or metadata.get("genome_scope", "all"),
                    "scaffold": "",
                    "targets_completed": 0,
                    "targets_total": 0,
                    "target_chunks": target_chunks,
                }
            )

        return MatrixCompareSummary(
            output_file=output_file,
            requested_pairs=requested_pairs,
            written_rows=written_rows,
            scaffold_count=len(scaffolds),
            genome_count=len({spec.genome for spec in scaffolds}),
            anchor_groups=anchor_groups,
            target_chunks=target_chunks,
            min_cov=min_cov,
            genome_scope=genome_scope or metadata.get("genome_scope", "all"),
            backend=compute_backend.kind,
            device=compute_backend.device,
            memory_limit_gb=memory_limit_gb,
            position_tile_size=position_tile_size or default_tile_size,
        )
    finally:
        conn.close()
