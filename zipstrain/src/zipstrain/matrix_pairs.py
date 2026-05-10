from __future__ import annotations

import asyncio
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import gc
import importlib
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import duckdb
import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

from zipstrain import compare as cp

PROFILE_REQUIRED_COLUMNS = {"chrom", "pos", "gene", "genome", "A", "T", "C", "G"}
COUNT_DTYPES = {
    "uint16": np.uint16,
    "uint32": np.uint32,
}
MATRIX_BUILD_MIN_COV = 5
FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS = "allele_presence_after_cov_filter"
CURRENT_MATRIX_DB_LAYOUT = "per_sample_per_genome_dense_matrix"
MATRIX_PAIR_BACKENDS = ("numpy", "torch", "torch-cpu", "torch-cuda", "torch-mps")
MATRIX_IO_EXECUTOR_KINDS = ("thread", "process")
BuildProgressCallback = Callable[[dict[str, object]], None]
CompareProgressCallback = Callable[[dict[str, object]], None]
MATRIX_COMPARISON_SUPPORTED_CALCULATIONS = ("ani", "ibs")
MATRIX_BUILD_FIXED_HEADROOM_BYTES = 128 * 1024 * 1024
MATRIX_BUILD_TEMP_BYTES_PER_POSITION = 32
MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES = 256 * 1024 * 1024
MATRIX_BUILD_MIN_COMMIT_BATCH_BYTES = 64 * 1024 * 1024
MATRIX_BUILD_DUCKDB_MEMORY_FRACTION = 0.25
MATRIX_COMPARE_CHECKPOINT_BATCH_UNITS = 4
MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS = 16
MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE_DEFAULT = 1
MPS_MAX_GRAPH_TENSOR_ELEMENTS = (2**31) - 1
MPS_GRAPH_TENSOR_HEADROOM_FRACTION = 0.9


def _emit_matrix_compare_writer_log(
    *,
    start_time: float,
    completed: int,
    total: int,
    batch_pairs: int,
    batch_rows: int,
    anchor_name: str,
    genome: str,
    targets_completed: int,
    targets_total: int,
    target_chunks: int,
) -> None:
    percent = (completed / max(total, 1)) * 100 if total > 0 else 0.0
    fields = {
        "elapsed": f"{max(time.monotonic() - start_time, 0.0):.1f}s",
        "completed": completed,
        "total": total,
        "percent": f"{percent:.1f}",
        "batch_pairs": batch_pairs,
        "batch_rows": batch_rows,
        "anchor_name": anchor_name,
        "genome": genome,
        "targets_completed": targets_completed,
        "targets_total": targets_total,
        "target_chunks": target_chunks,
    }
    payload = " ".join(f"{key}={value}" for key, value in fields.items())
    sys.stderr.write(f"MATRIX-COMPARE PROGRESS {payload}\n")
    sys.stderr.flush()

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


def matrix_compare_result_db_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("sample_idx_1", pa.int64()),
            pa.field("sample_idx_2", pa.int64()),
            pa.field("sample_1", pa.string()),
            pa.field("sample_2", pa.string()),
            pa.field("genome_idx", pa.int64()),
            pa.field("genome", pa.string()),
            pa.field("total_positions", pa.int64()),
            pa.field("share_allele_pos", pa.int64()),
            pa.field("genome_pop_ani", pa.float64()),
            pa.field("max_consecutive_length", pa.int64()),
        ]
    )


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
class GenomeSpec:
    genome_idx: int
    genome: str
    matrix_length: int
    true_length: int
    scaffold_count: int


@dataclass(frozen=True)
class GenomeScaffoldOffset:
    genome_idx: int
    scaffold_ordinal: int
    genome: str
    chrom: str
    axis_start: int
    axis_end: int
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
class MatrixDbAppendSummary:
    output_file: Path
    appended_profile_files: int
    appended_sample_count: int
    total_sample_count: int
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


def _build_genome_specs(
    scaffolds: list[ScaffoldSpec],
) -> tuple[list[GenomeSpec], list[GenomeScaffoldOffset]]:
    genomes: list[GenomeSpec] = []
    offsets: list[GenomeScaffoldOffset] = []
    grouped: dict[str, list[ScaffoldSpec]] = {}
    for spec in scaffolds:
        grouped.setdefault(spec.genome, []).append(spec)

    for genome_idx, genome_name in enumerate(sorted(grouped)):
        genome_scaffolds = sorted(grouped[genome_name], key=lambda spec: (spec.chrom, spec.index_base))
        axis_cursor = 0
        true_length = 0
        for ordinal, spec in enumerate(genome_scaffolds):
            axis_start = axis_cursor
            axis_end = axis_start + spec.vector_length - 1
            offsets.append(
                GenomeScaffoldOffset(
                    genome_idx=genome_idx,
                    scaffold_ordinal=ordinal,
                    genome=genome_name,
                    chrom=spec.chrom,
                    axis_start=axis_start,
                    axis_end=axis_end,
                    index_base=spec.index_base,
                    vector_length=spec.vector_length,
                    min_pos=spec.min_pos,
                    max_pos=spec.max_pos,
                )
            )
            axis_cursor = axis_end + 1
            true_length += spec.vector_length
            if ordinal < len(genome_scaffolds) - 1:
                axis_cursor += 1
        genomes.append(
            GenomeSpec(
                genome_idx=genome_idx,
                genome=genome_name,
                matrix_length=axis_cursor,
                true_length=true_length,
                scaffold_count=len(genome_scaffolds),
            )
        )
    return genomes, offsets


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


def _estimate_builder_python_peak_bytes(matrix_length: int, dtype_name: str) -> int:
    matrix_bytes = _estimate_sample_scaffold_bytes(matrix_length, dtype_name)
    temp_bytes = matrix_length * MATRIX_BUILD_TEMP_BYTES_PER_POSITION
    return matrix_bytes * 2 + temp_bytes + MATRIX_BUILD_FIXED_HEADROOM_BYTES


def _format_memory_bytes(num_bytes: int) -> str:
    return f"{num_bytes / (1024 ** 3):.2f} GB"


def _duckdb_memory_limit_setting(memory_limit_bytes: int) -> str:
    memory_limit_mib = max(1, int(memory_limit_bytes // (1024 ** 2)))
    return f"{memory_limit_mib}MiB"


def _plan_matrix_build_memory(
    genomes: list[GenomeSpec],
    count_dtype: str,
    memory_limit_bytes: int,
    commit_batch_gb: Optional[float] = None,
) -> MatrixBuildMemoryPlan:
    if not genomes:
        raise ValueError("No genomes were provided for matrix DB build.")

    limiting_spec = max(
        genomes,
        key=lambda spec: _estimate_builder_python_peak_bytes(spec.matrix_length, count_dtype),
    )
    estimated_python_peak_bytes = _estimate_builder_python_peak_bytes(
        limiting_spec.matrix_length,
        count_dtype,
    )

    if estimated_python_peak_bytes + MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES > memory_limit_bytes:
        raise MemoryError(
            f"Genome {limiting_spec.genome} is estimated to require about "
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
        limiting_scaffold=limiting_spec.genome,
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
        CREATE TABLE matrix_db_genomes (
          genome_idx INTEGER PRIMARY KEY,
          genome VARCHAR NOT NULL UNIQUE,
          matrix_length BIGINT NOT NULL,
          true_length BIGINT NOT NULL,
          scaffold_count INTEGER NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_db_genome_scaffolds (
          genome_idx INTEGER NOT NULL,
          scaffold_ordinal INTEGER NOT NULL,
          genome VARCHAR NOT NULL,
          chrom VARCHAR NOT NULL,
          axis_start BIGINT NOT NULL,
          axis_end BIGINT NOT NULL,
          index_base BIGINT NOT NULL,
          vector_length BIGINT NOT NULL,
          min_pos BIGINT NOT NULL,
          max_pos BIGINT NOT NULL,
          PRIMARY KEY (genome_idx, scaffold_ordinal)
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_db_sample_genome_matrices (
          sample_idx INTEGER NOT NULL,
          genome_idx INTEGER NOT NULL,
          count_dtype VARCHAR NOT NULL,
          matrix_rows BIGINT NOT NULL,
          matrix_cols INTEGER NOT NULL,
          matrix_blob BLOB NOT NULL,
          PRIMARY KEY (sample_idx, genome_idx)
        )
        """
    )


def _init_matrix_compare_db_schema(conn: duckdb.DuckDBPyConnection) -> None:
    conn.execute(
        """
        CREATE TABLE matrix_compare_metadata (
          key VARCHAR PRIMARY KEY,
          value VARCHAR NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_compare_samples (
          sample_idx INTEGER PRIMARY KEY,
          sample_name VARCHAR NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_compare_genomes (
          genome_idx INTEGER PRIMARY KEY,
          genome VARCHAR NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_compare_completed_pair_genomes (
          sample_idx_1 INTEGER NOT NULL,
          sample_idx_2 INTEGER NOT NULL,
          genome_idx INTEGER NOT NULL,
          PRIMARY KEY (sample_idx_1, sample_idx_2, genome_idx)
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE matrix_compare_results (
          sample_idx_1 INTEGER NOT NULL,
          sample_idx_2 INTEGER NOT NULL,
          sample_1 VARCHAR NOT NULL,
          sample_2 VARCHAR NOT NULL,
          genome_idx INTEGER NOT NULL,
          genome VARCHAR NOT NULL,
          total_positions BIGINT,
          share_allele_pos BIGINT,
          genome_pop_ani DOUBLE,
          max_consecutive_length BIGINT
        )
        """
    )


def _compare_db_metadata_rows(
    *,
    matrix_metadata: dict[str, str],
    genome_scope: Optional[str],
    calculations: tuple[str, ...],
) -> list[tuple[str, str]]:
    return [
        ("matrix_db_layout", matrix_metadata.get("layout", "")),
        ("matrix_value_semantics", matrix_metadata.get("matrix_value_semantics", "")),
        ("matrix_count_dtype", matrix_metadata.get("count_dtype", "")),
        ("genome_scope", genome_scope or "all"),
        ("calculate", "+".join(calculations)),
    ]


def _validate_matrix_compare_db_metadata(
    *,
    compare_metadata: dict[str, str],
    matrix_metadata: dict[str, str],
    genome_scope: Optional[str],
    calculations: tuple[str, ...],
) -> None:
    expected = dict(
        _compare_db_metadata_rows(
            matrix_metadata=matrix_metadata,
            genome_scope=genome_scope,
            calculations=calculations,
        )
    )
    mismatches: list[str] = []
    for key, expected_value in expected.items():
        actual = compare_metadata.get(key)
        if actual != expected_value:
            mismatches.append(f"{key} expected={expected_value!r} actual={actual!r}")
    if mismatches:
        raise ValueError(
            "Existing matrix compare DB is incompatible with this run: "
            + "; ".join(mismatches)
        )


def _matrix_compare_table_exists(
    compare_conn: duckdb.DuckDBPyConnection,
    table_name: str,
) -> bool:
    return bool(
        compare_conn.execute(
            """
            SELECT count(*)
            FROM information_schema.tables
            WHERE table_schema = current_schema()
              AND table_name = ?
            """,
            [table_name],
        ).fetchone()[0]
    )


def _ensure_matrix_compare_completed_table(
    compare_conn: duckdb.DuckDBPyConnection,
) -> None:
    if _matrix_compare_table_exists(compare_conn, "matrix_compare_completed_pair_genomes"):
        return
    compare_conn.execute(
        """
        CREATE TABLE matrix_compare_completed_pair_genomes (
          sample_idx_1 INTEGER NOT NULL,
          sample_idx_2 INTEGER NOT NULL,
          genome_idx INTEGER NOT NULL,
          PRIMARY KEY (sample_idx_1, sample_idx_2, genome_idx)
        )
        """
    )
    compare_conn.execute(
        """
        INSERT OR IGNORE INTO matrix_compare_completed_pair_genomes
        SELECT DISTINCT sample_idx_1, sample_idx_2, genome_idx
        FROM matrix_compare_results
        """
    )


def _prepare_matrix_compare_db(
    *,
    output_file: Path,
    matrix_metadata: dict[str, str],
    samples: list[tuple[int, str]],
    genomes: list[GenomeSpec],
    genome_scope: Optional[str],
    calculations: tuple[str, ...],
) -> duckdb.DuckDBPyConnection:
    existed_before = output_file.exists() and output_file.stat().st_size > 0
    compare_conn = duckdb.connect(str(output_file))
    compare_conn.execute("SET preserve_insertion_order=false")
    if not existed_before:
        _init_matrix_compare_db_schema(compare_conn)
        compare_conn.executemany(
            "INSERT INTO matrix_compare_metadata VALUES (?, ?)",
            _compare_db_metadata_rows(
                matrix_metadata=matrix_metadata,
                genome_scope=genome_scope,
                calculations=calculations,
            ),
        )
    else:
        if not _matrix_compare_table_exists(compare_conn, "matrix_compare_metadata"):
            raise ValueError(
                "Existing matrix compare DB is missing metadata and cannot be resumed safely. "
                "Remove it and rerun matrix-compare."
            )
        compare_metadata = {
            str(k): str(v)
            for k, v in compare_conn.execute("SELECT key, value FROM matrix_compare_metadata").fetchall()
        }
        _validate_matrix_compare_db_metadata(
            compare_metadata=compare_metadata,
            matrix_metadata=matrix_metadata,
            genome_scope=genome_scope,
            calculations=calculations,
        )
        _ensure_matrix_compare_completed_table(compare_conn)

    compare_conn.executemany(
        "INSERT OR REPLACE INTO matrix_compare_samples VALUES (?, ?)",
        [(sample_idx, sample_name) for sample_idx, sample_name in samples],
    )
    compare_conn.executemany(
        "INSERT OR REPLACE INTO matrix_compare_genomes VALUES (?, ?)",
        [(spec.genome_idx, spec.genome) for spec in genomes],
    )
    compare_conn.commit()
    return compare_conn


def _load_matrix_compare_resume_state(
    compare_conn: duckdb.DuckDBPyConnection,
    sample_count: int,
    genome_ids: list[int],
) -> tuple[dict[int, set[tuple[int, int]]], int, int]:
    total_pairs = sample_count * (sample_count - 1) // 2
    if total_pairs == 0 or not genome_ids:
        return {genome_idx: set() for genome_idx in genome_ids}, 0, 0

    completed_rows = compare_conn.execute(
        """
        SELECT sample_idx_1, sample_idx_2, genome_idx
        FROM matrix_compare_completed_pair_genomes
        WHERE sample_idx_1 < ? AND sample_idx_2 < ?
        """,
        [sample_count, sample_count],
    ).fetchall()
    completed_by_genome: dict[int, set[tuple[int, int]]] = {
        genome_idx: set() for genome_idx in genome_ids
    }
    completed_pair_counts: dict[tuple[int, int], int] = {}
    valid_genome_ids = set(genome_ids)
    completed_pair_genomes = 0
    for sample_idx_1, sample_idx_2, genome_idx in completed_rows:
        genome_idx = int(genome_idx)
        if genome_idx not in valid_genome_ids:
            continue
        pair = (int(sample_idx_1), int(sample_idx_2))
        genome_pairs = completed_by_genome[genome_idx]
        if pair in genome_pairs:
            continue
        genome_pairs.add(pair)
        completed_pair_counts[pair] = completed_pair_counts.get(pair, 0) + 1
        completed_pair_genomes += 1

    fully_completed_pairs = sum(
        1 for count in completed_pair_counts.values() if count == len(genome_ids)
    )
    remaining_pairs = max(0, total_pairs - fully_completed_pairs)
    remaining_pair_genomes = max(0, total_pairs * len(genome_ids) - completed_pair_genomes)
    return completed_by_genome, remaining_pairs, remaining_pair_genomes


def _insert_matrix_compare_result_table(
    compare_conn: duckdb.DuckDBPyConnection,
    table: pa.Table,
) -> None:
    if table.num_rows == 0:
        return
    compare_conn.register("_matrix_compare_result_batch", table)
    try:
        compare_conn.execute(
            """
            INSERT INTO matrix_compare_results
            SELECT
              sample_idx_1,
              sample_idx_2,
              sample_1,
              sample_2,
              genome_idx,
              genome,
              total_positions,
              share_allele_pos,
              genome_pop_ani,
              max_consecutive_length
            FROM _matrix_compare_result_batch
            """
        )
    finally:
        compare_conn.unregister("_matrix_compare_result_batch")


def _mark_completed_pair_genomes(
    compare_conn: duckdb.DuckDBPyConnection,
    completed_pair_genomes: list[tuple[int, int, int]],
) -> None:
    if not completed_pair_genomes:
        return
    table = pa.table(
        {
            "sample_idx_1": pa.array(
                np.fromiter((row[0] for row in completed_pair_genomes), dtype=np.int64),
                type=pa.int64(),
            ),
            "sample_idx_2": pa.array(
                np.fromiter((row[1] for row in completed_pair_genomes), dtype=np.int64),
                type=pa.int64(),
            ),
            "genome_idx": pa.array(
                np.fromiter((row[2] for row in completed_pair_genomes), dtype=np.int64),
                type=pa.int64(),
            ),
        }
    )
    compare_conn.register("_matrix_compare_completed_pair_genome_batch", table)
    try:
        compare_conn.execute(
            """
            INSERT OR REPLACE INTO matrix_compare_completed_pair_genomes
            SELECT DISTINCT sample_idx_1, sample_idx_2, genome_idx
            FROM _matrix_compare_completed_pair_genome_batch
            """
        )
    finally:
        compare_conn.unregister("_matrix_compare_completed_pair_genome_batch")


def _completed_pair_rows_from_payload(
    payload: dict[str, object],
) -> list[tuple[int, int, int]]:
    sample_1_idx = int(payload["sample_1_idx"])
    genome_idx = int(payload["genome_idx"])
    sample_2_idx = np.asarray(payload["sample_2_idx"], dtype=np.int64)
    return [
        (sample_1_idx, int(sample_2_idx_value), genome_idx)
        for sample_2_idx_value in sample_2_idx.tolist()
    ]


def _write_matrix_compare_batch(
    output_file: Path,
    batch_tables: list[pa.Table],
    batch_completed_rows: list[tuple[int, int, int]],
) -> int:
    batch_rows = 0
    write_conn = duckdb.connect(str(output_file))
    write_conn.execute("SET preserve_insertion_order=false")
    try:
        write_conn.execute("BEGIN")
        try:
            for table in batch_tables:
                _insert_matrix_compare_result_table(write_conn, table)
                batch_rows += table.num_rows
            _mark_completed_pair_genomes(write_conn, batch_completed_rows)
            write_conn.execute("COMMIT")
        except Exception:
            write_conn.execute("ROLLBACK")
            raise
    finally:
        write_conn.close()
    return batch_rows


def _write_matrix_compare_payload_batch(
    output_file: Path,
    batch_payloads: list[dict[str, object]],
    log_context: Optional[dict[str, object]] = None,
) -> int:
    batch_rows = 0
    completed_rows: list[tuple[int, int, int]] = []
    write_conn = duckdb.connect(str(output_file))
    write_conn.execute("SET preserve_insertion_order=false")
    try:
        write_conn.execute("BEGIN")
        try:
            for payload in batch_payloads:
                completed_rows.extend(_completed_pair_rows_from_payload(payload))
                table = _make_arrow_table_from_compare_payload(payload)
                if table is None:
                    continue
                _insert_matrix_compare_result_table(write_conn, table)
                batch_rows += table.num_rows
            _mark_completed_pair_genomes(write_conn, completed_rows)
            write_conn.execute("COMMIT")
        except Exception:
            write_conn.execute("ROLLBACK")
            raise
    finally:
        write_conn.close()
    if log_context is not None:
        _emit_matrix_compare_writer_log(
            start_time=float(log_context["start_time"]),
            completed=int(log_context["completed"]),
            total=int(log_context["total"]),
            batch_pairs=int(log_context["batch_pairs"]),
            batch_rows=batch_rows,
            anchor_name=str(log_context["anchor_name"]),
            genome=str(log_context["genome"]),
            targets_completed=int(log_context["targets_completed"]),
            targets_total=int(log_context["targets_total"]),
            target_chunks=int(log_context["target_chunks"]),
        )
    return batch_rows


def _validate_matrix_db_appendable(metadata: dict[str, str]) -> tuple[str, str]:
    layout = metadata.get("layout", "")
    if layout != CURRENT_MATRIX_DB_LAYOUT:
        raise ValueError(
            "Matrix DB layout is not append-compatible with this builder. "
            "Rebuild the matrix DB with the current builder."
        )
    semantics = metadata.get("matrix_value_semantics", "")
    if semantics != FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
        raise ValueError(
            "Matrix DB value semantics are incompatible with append. "
            "Rebuild the matrix DB with the current builder."
        )
    min_cov = metadata.get("coverage_filter_min_cov")
    if min_cov != str(MATRIX_BUILD_MIN_COV):
        raise ValueError(
            "Matrix DB coverage filter is incompatible with append. "
            "Rebuild the matrix DB with the current builder."
        )
    count_dtype = metadata.get("count_dtype", "")
    if count_dtype not in COUNT_DTYPES:
        raise ValueError(
            f"Matrix DB count dtype '{count_dtype}' is not supported for append."
        )
    genome_scope = metadata.get("genome_scope", "all")
    return count_dtype, genome_scope


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


def _load_profile_genome_matrix(
    profile_path: Path,
    genome_spec: GenomeSpec,
    genome_offsets: list[GenomeScaffoldOffset],
    count_dtype: str,
    min_cov: int,
) -> np.ndarray:
    np_dtype = COUNT_DTYPES[count_dtype]
    matrix = np.zeros((genome_spec.matrix_length, 4), dtype=np_dtype)
    frame = (
        pl.scan_parquet(profile_path)
        .filter(pl.col("genome") == genome_spec.genome)
        .select("chrom", "pos", "A", "T", "C", "G")
        .collect(engine="streaming")
    )
    if frame.height == 0:
        return matrix
    coverage = (
        frame["A"].cast(pl.Int64)
        + frame["T"].cast(pl.Int64)
        + frame["C"].cast(pl.Int64)
        + frame["G"].cast(pl.Int64)
    )
    frame = frame.filter(coverage >= min_cov)
    if frame.height == 0:
        return matrix
    for offset in genome_offsets:
        scaffold_frame = frame.filter(pl.col("chrom") == offset.chrom)
        if scaffold_frame.height == 0:
            continue
        axis_pos = (
            scaffold_frame["pos"].to_numpy().astype(np.int64, copy=False)
            - offset.index_base
            + offset.axis_start
        )
        matrix[axis_pos, 0] = (scaffold_frame["A"].to_numpy() > 0).astype(np_dtype, copy=False)
        matrix[axis_pos, 1] = (scaffold_frame["T"].to_numpy() > 0).astype(np_dtype, copy=False)
        matrix[axis_pos, 2] = (scaffold_frame["C"].to_numpy() > 0).astype(np_dtype, copy=False)
        matrix[axis_pos, 3] = (scaffold_frame["G"].to_numpy() > 0).astype(np_dtype, copy=False)
    return matrix


def _validate_profile_against_matrix_contract(
    profile_path: Path,
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
) -> None:
    genome_names = {spec.genome for spec in genomes}
    offsets_by_key = {(offset.genome, offset.chrom): offset for offset in genome_scaffolds}
    grouped = (
        pl.scan_parquet(profile_path)
        .select("genome", "chrom", "pos")
        .group_by(["genome", "chrom"])
        .agg(
            pl.col("pos").min().cast(pl.Int64).alias("min_pos"),
            pl.col("pos").max().cast(pl.Int64).alias("max_pos"),
        )
        .collect(engine="streaming")
    )
    if grouped.height == 0:
        return

    unknown_genomes: set[str] = set()
    unknown_scaffolds: set[str] = set()
    out_of_range: list[str] = []
    for row in grouped.iter_rows(named=True):
        genome_name = str(row["genome"])
        chrom = str(row["chrom"])
        min_pos = int(row["min_pos"])
        max_pos = int(row["max_pos"])
        if genome_name not in genome_names:
            unknown_genomes.add(genome_name)
            continue
        offset = offsets_by_key.get((genome_name, chrom))
        if offset is None:
            unknown_scaffolds.add(f"{genome_name}:{chrom}")
            continue
        if min_pos < offset.min_pos or max_pos > offset.max_pos:
            out_of_range.append(
                f"{genome_name}:{chrom} profile_pos_range={min_pos}-{max_pos} "
                f"matrix_pos_range={offset.min_pos}-{offset.max_pos}"
            )

    errors: list[str] = []
    if unknown_genomes:
        errors.append("unknown genomes " + ", ".join(sorted(unknown_genomes)))
    if unknown_scaffolds:
        errors.append("unknown scaffolds " + ", ".join(sorted(unknown_scaffolds)))
    if out_of_range:
        errors.append("out-of-range positions " + "; ".join(sorted(out_of_range)))
    if errors:
        raise ValueError(
            f"Profile {profile_path.name} does not match the existing matrix DB contract: "
            + " | ".join(errors)
        )


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
    genomes, genome_scaffolds = _build_genome_specs(scaffolds)
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    memory_plan = _plan_matrix_build_memory(
        genomes=genomes,
        count_dtype=count_dtype,
        memory_limit_bytes=memory_limit_bytes,
        commit_batch_gb=commit_batch_gb,
    )

    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    total_work = len(profile_paths) * len(genomes)
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
            ("profiles_dir", str(profile_dir.resolve())),
            ("profile_format", "classic_zipstrain_profile_parquet"),
            ("genome_scope", genome_scope or "all"),
            ("count_dtype", count_dtype),
            ("layout", CURRENT_MATRIX_DB_LAYOUT),
            ("matrix_value_semantics", FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS),
            ("coverage_filter_min_cov", str(MATRIX_BUILD_MIN_COV)),
            ("memory_limit_gb", str(memory_limit_gb)),
            ("duckdb_memory_limit_gb", f"{memory_plan.duckdb_memory_limit_bytes / (1024 ** 3):.6f}"),
            ("effective_commit_batch_gb", f"{memory_plan.commit_batch_bytes / (1024 ** 3):.6f}"),
            ("estimated_python_peak_gb", f"{memory_plan.estimated_python_peak_bytes / (1024 ** 3):.6f}"),
            ("separator_rows_between_scaffolds", "1"),
        ]
        sample_rows = [(idx, path.stem, str(path.resolve())) for idx, path in enumerate(profile_paths)]
        genome_rows = [
            (
                spec.genome_idx,
                spec.genome,
                spec.matrix_length,
                spec.true_length,
                spec.scaffold_count,
            )
            for spec in genomes
        ]
        scaffold_rows = [
            (
                spec.genome_idx,
                spec.scaffold_ordinal,
                spec.genome,
                spec.chrom,
                spec.axis_start,
                spec.axis_end,
                spec.index_base,
                spec.vector_length,
                spec.min_pos,
                spec.max_pos,
            )
            for spec in genome_scaffolds
        ]
        conn.executemany("INSERT INTO matrix_db_metadata VALUES (?, ?)", metadata_rows)
        conn.executemany("INSERT INTO matrix_db_samples VALUES (?, ?, ?)", sample_rows)
        conn.executemany(
            "INSERT INTO matrix_db_genomes VALUES (?, ?, ?, ?, ?)",
            genome_rows,
        )
        conn.executemany(
            "INSERT INTO matrix_db_genome_scaffolds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            scaffold_rows,
        )
        conn = _restart_matrix_build_transaction(conn, output_file, memory_plan.duckdb_memory_limit_bytes)

        stored_rows = 0
        scaffolds_by_genome_idx: dict[int, list[GenomeScaffoldOffset]] = {}
        for offset in genome_scaffolds:
            scaffolds_by_genome_idx.setdefault(offset.genome_idx, []).append(offset)
        for sample_idx, profile_path in enumerate(profile_paths):
            sample_name = profile_path.stem
            for genome_spec in genomes:
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "processing",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": genome_spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )
                estimated_python_peak = _estimate_builder_python_peak_bytes(genome_spec.matrix_length, count_dtype)
                if estimated_python_peak + memory_plan.duckdb_memory_limit_bytes > memory_limit_bytes:
                    raise MemoryError(
                        f"Genome {genome_spec.genome} for sample {profile_path.name} is estimated to use about "
                        f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory plus "
                        f"{_format_memory_bytes(memory_plan.duckdb_memory_limit_bytes)} for DuckDB, "
                        f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                    )
                matrix = _load_profile_genome_matrix(
                    profile_path=profile_path,
                    genome_spec=genome_spec,
                    genome_offsets=scaffolds_by_genome_idx[genome_spec.genome_idx],
                    count_dtype=count_dtype,
                    min_cov=MATRIX_BUILD_MIN_COV,
                )
                conn.execute(
                    """
                    INSERT INTO matrix_db_sample_genome_matrices
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    [
                        sample_idx,
                        genome_spec.genome_idx,
                        count_dtype,
                        int(matrix.shape[0]),
                        int(matrix.shape[1]),
                        _pack_matrix(matrix),
                    ],
                )
                stored_rows += 1
                completed_work += 1
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "advance",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": genome_spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )
                conn = _restart_matrix_build_transaction(
                    conn,
                    output_file,
                    memory_plan.duckdb_memory_limit_bytes,
                )
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


def append_matrix_db(
    profile_dir: Path,
    matrix_db_file: Path,
    memory_limit_gb: float = 16.0,
    progress_callback: Optional[BuildProgressCallback] = None,
) -> MatrixDbAppendSummary:
    profile_paths = discover_profile_parquets(profile_dir)
    matrix_db_file = Path(matrix_db_file).resolve()
    if not matrix_db_file.exists():
        raise FileNotFoundError(f"Matrix DB file does not exist: {matrix_db_file}")

    read_conn = duckdb.connect(str(matrix_db_file), read_only=True)
    try:
        metadata = _load_matrix_db_metadata(read_conn)
        count_dtype, genome_scope = _validate_matrix_db_appendable(metadata)
        genomes = _load_matrix_db_genomes(read_conn)
        genome_scaffolds = _load_matrix_db_genome_scaffolds(read_conn)
        existing_samples = _load_matrix_db_samples(read_conn)
    finally:
        read_conn.close()

    existing_sample_names = {sample_name for _sample_idx, sample_name in existing_samples}
    duplicate_sample_names = sorted(path.stem for path in profile_paths if path.stem in existing_sample_names)
    if duplicate_sample_names:
        raise ValueError(
            "Cannot append profiles whose sample names already exist in the matrix DB: "
            + ", ".join(duplicate_sample_names)
        )

    for profile_path in profile_paths:
        _validate_profile_against_matrix_contract(
            profile_path=profile_path,
            genomes=genomes,
            genome_scaffolds=genome_scaffolds,
        )

    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    memory_plan = _plan_matrix_build_memory(
        genomes=genomes,
        count_dtype=count_dtype,
        memory_limit_bytes=memory_limit_bytes,
    )

    next_sample_idx = max((sample_idx for sample_idx, _sample_name in existing_samples), default=-1) + 1
    appended_samples = [
        (next_sample_idx + idx, profile_path.stem, str(profile_path.resolve()), profile_path)
        for idx, profile_path in enumerate(profile_paths)
    ]
    total_work = len(appended_samples) * len(genomes)
    completed_work = 0
    stored_rows = 0

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

    scaffolds_by_genome_idx: dict[int, list[GenomeScaffoldOffset]] = {}
    for offset in genome_scaffolds:
        scaffolds_by_genome_idx.setdefault(offset.genome_idx, []).append(offset)

    conn = _open_matrix_build_connection(matrix_db_file, memory_plan.duckdb_memory_limit_bytes)
    append_succeeded = False
    appended_sample_indices = [sample_idx for sample_idx, _sample_name, _profile_path, _profile_file in appended_samples]
    try:
        conn.execute("BEGIN")
        conn.executemany(
            "INSERT INTO matrix_db_samples VALUES (?, ?, ?)",
            [(sample_idx, sample_name, profile_path_str) for sample_idx, sample_name, profile_path_str, _profile_path in appended_samples],
        )
        conn = _restart_matrix_build_transaction(conn, matrix_db_file, memory_plan.duckdb_memory_limit_bytes)

        for sample_idx, sample_name, _profile_path_str, profile_path in appended_samples:
            for genome_spec in genomes:
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "processing",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": genome_spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )
                estimated_python_peak = _estimate_builder_python_peak_bytes(genome_spec.matrix_length, count_dtype)
                if estimated_python_peak + memory_plan.duckdb_memory_limit_bytes > memory_limit_bytes:
                    raise MemoryError(
                        f"Genome {genome_spec.genome} for sample {profile_path.name} is estimated to use about "
                        f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory plus "
                        f"{_format_memory_bytes(memory_plan.duckdb_memory_limit_bytes)} for DuckDB, "
                        f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                    )
                matrix = _load_profile_genome_matrix(
                    profile_path=profile_path,
                    genome_spec=genome_spec,
                    genome_offsets=scaffolds_by_genome_idx[genome_spec.genome_idx],
                    count_dtype=count_dtype,
                    min_cov=MATRIX_BUILD_MIN_COV,
                )
                conn.execute(
                    """
                    INSERT INTO matrix_db_sample_genome_matrices
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    [
                        sample_idx,
                        genome_spec.genome_idx,
                        count_dtype,
                        int(matrix.shape[0]),
                        int(matrix.shape[1]),
                        _pack_matrix(matrix),
                    ],
                )
                stored_rows += 1
                completed_work += 1
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "advance",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": genome_spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )
                conn = _restart_matrix_build_transaction(
                    conn,
                    matrix_db_file,
                    memory_plan.duckdb_memory_limit_bytes,
                )
        conn.execute("COMMIT")
        append_succeeded = True
    except Exception:
        try:
            conn.execute("ROLLBACK")
        except Exception:
            pass
        raise
    finally:
        conn.close()
        if not append_succeeded and appended_sample_indices:
            cleanup_conn = duckdb.connect(str(matrix_db_file))
            compare_conn.commit()
            compare_conn.execute("BEGIN")
            try:
                cleanup_conn.execute(
                    f"DELETE FROM matrix_db_sample_genome_matrices WHERE sample_idx IN ({','.join(['?'] * len(appended_sample_indices))})",
                    appended_sample_indices,
                )
                cleanup_conn.execute(
                    f"DELETE FROM matrix_db_samples WHERE sample_idx IN ({','.join(['?'] * len(appended_sample_indices))})",
                    appended_sample_indices,
                )
            finally:
                cleanup_conn.close()

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

    return MatrixDbAppendSummary(
        output_file=matrix_db_file,
        appended_profile_files=len(profile_paths),
        appended_sample_count=len(profile_paths),
        total_sample_count=len(existing_samples) + len(profile_paths),
        scaffold_count=len(genome_scaffolds),
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


def _load_matrix_db_genomes(conn: duckdb.DuckDBPyConnection, genome: Optional[str] = None) -> list[GenomeSpec]:
    if genome is None:
        rows = conn.execute(
            """
            SELECT genome_idx, genome, matrix_length, true_length, scaffold_count
            FROM matrix_db_genomes
            ORDER BY genome_idx
            """
        ).fetchall()
    else:
        rows = conn.execute(
            """
            SELECT genome_idx, genome, matrix_length, true_length, scaffold_count
            FROM matrix_db_genomes
            WHERE genome = ?
            ORDER BY genome_idx
            """,
            [genome],
        ).fetchall()
    return [
        GenomeSpec(
            genome_idx=int(genome_idx),
            genome=str(genome_name),
            matrix_length=int(matrix_length),
            true_length=int(true_length),
            scaffold_count=int(scaffold_count),
        )
        for genome_idx, genome_name, matrix_length, true_length, scaffold_count in rows
    ]


def _load_matrix_db_genome_scaffolds(
    conn: duckdb.DuckDBPyConnection,
    genome_idx: Optional[int] = None,
    genome: Optional[str] = None,
) -> list[GenomeScaffoldOffset]:
    if genome_idx is None and genome is None:
        rows = conn.execute(
            """
            SELECT genome_idx, scaffold_ordinal, genome, chrom, axis_start, axis_end, index_base, vector_length, min_pos, max_pos
            FROM matrix_db_genome_scaffolds
            ORDER BY genome_idx, scaffold_ordinal
            """
        ).fetchall()
    elif genome_idx is not None:
        rows = conn.execute(
            """
            SELECT genome_idx, scaffold_ordinal, genome, chrom, axis_start, axis_end, index_base, vector_length, min_pos, max_pos
            FROM matrix_db_genome_scaffolds
            WHERE genome_idx = ?
            ORDER BY scaffold_ordinal
            """,
            [genome_idx],
        ).fetchall()
    else:
        rows = conn.execute(
            """
            SELECT genome_idx, scaffold_ordinal, genome, chrom, axis_start, axis_end, index_base, vector_length, min_pos, max_pos
            FROM matrix_db_genome_scaffolds
            WHERE genome = ?
            ORDER BY scaffold_ordinal
            """,
            [genome],
        ).fetchall()
    return [
        GenomeScaffoldOffset(
            genome_idx=int(genome_idx_val),
            scaffold_ordinal=int(scaffold_ordinal),
            genome=str(genome_name),
            chrom=str(chrom),
            axis_start=int(axis_start),
            axis_end=int(axis_end),
            index_base=int(index_base),
            vector_length=int(vector_length),
            min_pos=int(min_pos),
            max_pos=int(max_pos),
        )
        for genome_idx_val, scaffold_ordinal, genome_name, chrom, axis_start, axis_end, index_base, vector_length, min_pos, max_pos in rows
    ]


def _load_sample_genome_matrices(
    conn: duckdb.DuckDBPyConnection,
    genome_idx: int,
    sample_ids: list[int],
    matrix_length: int,
    dtype_name: str,
) -> dict[int, np.ndarray]:
    if not sample_ids:
        raise ValueError("sample_ids must not be empty")
    np_dtype = COUNT_DTYPES[dtype_name]
    rows = conn.execute(
        f"""
        SELECT sample_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob
        FROM matrix_db_sample_genome_matrices
        WHERE genome_idx = ? AND sample_idx IN ({','.join(['?'] * len(sample_ids))})
        ORDER BY sample_idx
        """,
        [genome_idx, *sample_ids],
    ).fetchall()
    if not rows:
        return {}
    matrices: dict[int, np.ndarray] = {}
    for sample_idx, stored_dtype, matrix_rows, matrix_cols, matrix_blob in rows:
        rows_int = int(matrix_rows)
        cols_int = int(matrix_cols)
        if rows_int != matrix_length or cols_int != 4:
            raise ValueError(
                f"Stored genome matrix shape mismatch for genome_idx={genome_idx}, sample_idx={sample_idx}: "
                f"expected ({matrix_length}, 4), found ({rows_int}, {cols_int})"
            )
        if str(stored_dtype) != dtype_name:
            raise ValueError(
                f"Stored genome matrix dtype mismatch for genome_idx={genome_idx}, sample_idx={sample_idx}: "
                f"expected {dtype_name}, found {stored_dtype}"
            )
        matrices[int(sample_idx)] = _unpack_matrix(
            bytes(matrix_blob),
            dtype_name,
            (rows_int, cols_int),
        ).astype(np_dtype, copy=False)
    return matrices


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


def _effective_torch_tile_size(
    compute_backend: MatrixPairComputeBackend,
    requested_tile_size: int,
    target_count: int,
) -> int:
    tile_size = max(1, int(requested_tile_size))
    if compute_backend.device != "mps":
        return tile_size
    max_elements = max(
        1,
        int(MPS_MAX_GRAPH_TENSOR_ELEMENTS * MPS_GRAPH_TENSOR_HEADROOM_FRACTION),
    )
    per_position_target_elements = 4 * max(int(target_count), 1)
    max_tile_for_targets = max(1, max_elements // per_position_target_elements)
    max_tile_for_anchor = max(1, max_elements // 4)
    return min(tile_size, max_tile_for_targets, max_tile_for_anchor)


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
    sample_1_idx: int,
    sample_1: str,
    sample_2_idx: list[int],
    sample_2: list[str],
    genome_idx: int,
    genome: str,
    calculations: tuple[str, ...],
    total_positions: Optional[np.ndarray] = None,
    share_allele_pos: Optional[np.ndarray] = None,
    max_consecutive_length: Optional[np.ndarray] = None,
) -> pa.Table:
    row_count = len(sample_2)
    if len(sample_2_idx) != row_count:
        raise ValueError("sample_2_idx and sample_2 must have the same length.")

    total_values: list[Optional[int]] = [None] * row_count
    shared_values: list[Optional[int]] = [None] * row_count
    ani_values: list[Optional[float]] = [None] * row_count
    if "ani" in calculations:
        if total_positions is None or share_allele_pos is None:
            raise ValueError("ANI output requested but total_positions/share_allele_pos were not provided.")
        total_values = total_positions.astype(np.int64, copy=False).tolist()
        shared_values = share_allele_pos.astype(np.int64, copy=False).tolist()
        ani_values = (share_allele_pos / total_positions * 100.0).tolist()

    ibs_values: list[Optional[int]] = [None] * row_count
    if "ibs" in calculations:
        if max_consecutive_length is None:
            raise ValueError("IBS output requested but max_consecutive_length was not provided.")
        ibs_values = max_consecutive_length.astype(np.int64, copy=False).tolist()

    return pa.Table.from_arrays(
        [
            pa.array([sample_1_idx] * row_count, type=pa.int64()),
            pa.array(sample_2_idx, type=pa.int64()),
            pa.array([sample_1] * row_count, type=pa.string()),
            pa.array(sample_2, type=pa.string()),
            pa.array([genome_idx] * row_count, type=pa.int64()),
            pa.array([genome] * row_count, type=pa.string()),
            pa.array(total_values, type=pa.int64()),
            pa.array(shared_values, type=pa.int64()),
            pa.array(ani_values, type=pa.float64()),
            pa.array(ibs_values, type=pa.int64()),
        ],
        schema=matrix_compare_result_db_schema(),
    )


def _make_arrow_table_from_compare_payload(
    payload: dict[str, object],
) -> Optional[pa.Table]:
    calculations = tuple(payload["calculations"])
    total_positions_all = np.asarray(payload["total_positions"], dtype=np.int64)
    mask = total_positions_all > 0
    if not bool(mask.any()):
        return None

    sample_2_idx_all = np.asarray(payload["sample_2_idx"], dtype=np.int64)
    sample_2_all = list(payload["sample_2"])
    share_allele_pos_all = np.asarray(payload["share_allele_pos"], dtype=np.int64)
    max_consecutive_all = payload.get("max_consecutive_length")
    if max_consecutive_all is not None:
        max_consecutive_all = np.asarray(max_consecutive_all, dtype=np.int64)

    return _make_arrow_table(
        sample_1_idx=int(payload["sample_1_idx"]),
        sample_1=str(payload["sample_1"]),
        sample_2_idx=sample_2_idx_all[mask].astype(np.int64, copy=False).tolist(),
        sample_2=[sample_2_all[idx] for idx, keep in enumerate(mask) if keep],
        genome_idx=int(payload["genome_idx"]),
        genome=str(payload["genome"]),
        calculations=calculations,
        total_positions=total_positions_all[mask] if "ani" in calculations else None,
        share_allele_pos=share_allele_pos_all[mask] if "ani" in calculations else None,
        max_consecutive_length=max_consecutive_all[mask] if max_consecutive_all is not None else None,
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


def _load_anchor_queue_batch_for_torch(
    matrix_db_file: Path,
    genome_idx: int,
    batch_rows: list[tuple[int, str]],
    matrix_length: int,
    dtype_name: str,
    zero_matrix: np.ndarray,
) -> list[tuple[int, str, np.ndarray]]:
    if not batch_rows:
        return []
    sample_ids = [sample_idx for sample_idx, _sample_name in batch_rows]
    read_conn = duckdb.connect(str(matrix_db_file), read_only=True)
    try:
        loaded_anchors = _load_sample_genome_matrices(
            read_conn,
            genome_idx=genome_idx,
            sample_ids=sample_ids,
            matrix_length=matrix_length,
            dtype_name=dtype_name,
        )
    finally:
        read_conn.close()
    return [
        (
            sample_idx,
            sample_name,
            loaded_anchors.get(sample_idx, zero_matrix),
        )
        for sample_idx, sample_name in batch_rows
    ]


def _load_target_queue_block_for_torch(
    matrix_db_file: Path,
    genome_idx: int,
    block_rows: list[tuple[int, str]],
    matrix_length: int,
    dtype_name: str,
) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray]:
    block_ids = np.array([sample_idx for sample_idx, _sample_name in block_rows], dtype=np.int64)
    block_names = [sample_name for _sample_idx, sample_name in block_rows]
    zero_matrix = np.zeros((matrix_length, 4), dtype=COUNT_DTYPES[dtype_name])
    read_conn = duckdb.connect(str(matrix_db_file), read_only=True)
    try:
        loaded_targets = _load_sample_genome_matrices(
            read_conn,
            genome_idx=genome_idx,
            sample_ids=block_ids.tolist(),
            matrix_length=matrix_length,
            dtype_name=dtype_name,
        )
    finally:
        read_conn.close()
    target_matrices = np.stack(
        [loaded_targets.get(int(sample_idx), zero_matrix) for sample_idx in block_ids],
        axis=2,
    )
    return block_ids, block_names, zero_matrix, target_matrices


def _load_target_prefetch_unit_for_torch(
    unit_index: int,
    matrix_db_file: Path,
    genome_idx: int,
    block_rows: list[tuple[int, str]],
    matrix_length: int,
    dtype_name: str,
) -> tuple[int, tuple[np.ndarray, list[str], np.ndarray, np.ndarray]]:
    return (
        unit_index,
        _load_target_queue_block_for_torch(
            matrix_db_file=matrix_db_file,
            genome_idx=genome_idx,
            block_rows=block_rows,
            matrix_length=matrix_length,
            dtype_name=dtype_name,
        ),
    )


def _build_matrix_io_executor(kind: str):
    if kind == "thread":
        return ThreadPoolExecutor(max_workers=1)
    if kind == "process":
        return ProcessPoolExecutor(max_workers=1)
    raise ValueError(f"Unsupported matrix I/O executor kind '{kind}'.")


def _compare_anchor_against_target_chunk_torch(
    compute_backend: MatrixPairComputeBackend,
    anchor_torch,
    target_torch,
    vector_length: int,
    tile_size: int,
    matrix_value_semantics: str,
    need_ibs: bool = False,
) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    target_count = int(target_torch.shape[2])
    tile_size = _effective_torch_tile_size(
        compute_backend=compute_backend,
        requested_tile_size=tile_size,
        target_count=target_count,
    )
    chunk_totals_torch, chunk_shared_torch, max_runs = _compare_anchor_against_target_chunk_torch_device(
        compute_backend=compute_backend,
        anchor_torch=anchor_torch,
        target_torch=target_torch,
        vector_length=vector_length,
        tile_size=tile_size,
        matrix_value_semantics=matrix_value_semantics,
        need_ibs=need_ibs,
    )
    combined_np = _download_torch_result_tensor_batch(
        compute_backend=compute_backend,
        totals_tensors=[chunk_totals_torch],
        shared_tensors=[chunk_shared_torch],
    )
    return combined_np[0, 0, :target_count].copy(), combined_np[0, 1, :target_count].copy(), max_runs


def _compare_anchor_against_target_chunk_torch_device(
    compute_backend: MatrixPairComputeBackend,
    anchor_torch,
    target_torch,
    vector_length: int,
    tile_size: int,
    matrix_value_semantics: str,
    need_ibs: bool = False,
):
    target_count = int(target_torch.shape[2])
    tile_size = _effective_torch_tile_size(
        compute_backend=compute_backend,
        requested_tile_size=tile_size,
        target_count=target_count,
    )
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
    return chunk_totals_torch, chunk_shared_torch, max_runs


def _download_torch_result_tensor_batch(
    compute_backend: MatrixPairComputeBackend,
    totals_tensors: list,
    shared_tensors: list,
) -> np.ndarray:
    if not totals_tensors:
        return np.zeros((0, 2, 0), dtype=np.int64)
    if len(totals_tensors) != len(shared_tensors):
        raise ValueError("totals_tensors and shared_tensors must have the same length")
    torch_module = compute_backend.torch
    lengths = [int(tensor.shape[0]) for tensor in totals_tensors]
    max_length = max(lengths)
    padded_totals = []
    padded_shared = []
    for totals_tensor, shared_tensor, length in zip(totals_tensors, shared_tensors, lengths):
        if int(shared_tensor.shape[0]) != length:
            raise ValueError("totals and shared tensors must have matching lengths")
        pad_width = max_length - length
        if pad_width > 0:
            padded_totals.append(torch_module.nn.functional.pad(totals_tensor, (0, pad_width)))
            padded_shared.append(torch_module.nn.functional.pad(shared_tensor, (0, pad_width)))
        else:
            padded_totals.append(totals_tensor)
            padded_shared.append(shared_tensor)
    stacked_totals = torch_module.stack(padded_totals, dim=0)
    stacked_shared = torch_module.stack(padded_shared, dim=0)
    combined = torch_module.stack((stacked_totals, stacked_shared), dim=1)
    if compute_backend.device == "cuda":
        combined_cpu = torch_module.empty(
            combined.shape,
            dtype=combined.dtype,
            device="cpu",
            pin_memory=True,
        )
        combined_cpu.copy_(combined, non_blocking=True)
        torch_module.cuda.current_stream().synchronize()
        combined_np = combined_cpu.numpy()
    else:
        combined_np = combined.cpu().numpy()
    return combined_np


async def _matrix_compare_reuse_target_chunks_torch_async(
    matrix_db_file: Path,
    conn: duckdb.DuckDBPyConnection,
    compare_conn: Optional[duckdb.DuckDBPyConnection],
    output_file: Path,
    min_cov: int,
    genome_scope: Optional[str],
    metadata: dict[str, str],
    samples: list[tuple[int, str]],
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
    completed_pairs_by_genome: dict[int, set[tuple[int, int]]],
    requested_pairs: int,
    remaining_work: int,
    memory_limit_bytes: int,
    memory_limit_gb: float,
    anchor_queue_size: int,
    target_queue_size: int,
    result_transfer_batch_size: int,
    loader_executor_kind: str,
    writer_executor_kind: str,
    position_tile_size: Optional[int],
    compute_backend: MatrixPairComputeBackend,
    matrix_value_semantics: str,
    calculations: tuple[str, ...],
    emit_writer_logs: bool,
    run_start_time: float,
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    total_work = remaining_work
    completed_work = 0
    written_rows = 0
    target_chunks = 0
    anchor_groups = max(len(samples) - 1, 0)
    dtype_name = str(metadata.get("count_dtype", "uint16"))
    default_tile_size = 0
    max_vector_length = max(spec.matrix_length for spec in genomes)
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
    loop = asyncio.get_running_loop()
    target_loader_executor = _build_matrix_io_executor(loader_executor_kind)
    anchor_loader_executor = _build_matrix_io_executor(loader_executor_kind)
    writer_executor = _build_matrix_io_executor(writer_executor_kind)
    pending_write_futures: deque[tuple[asyncio.Future[int], list[dict[str, object]]]] = deque()

    def load_target_block_sync(
        spec: GenomeSpec,
        block_rows: list[tuple[int, str]],
    ) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray]:
        block_ids = np.array([sample_idx for sample_idx, _sample_name in block_rows], dtype=np.int64)
        block_names = [sample_name for _sample_idx, sample_name in block_rows]
        zero_matrix = np.zeros((spec.matrix_length, 4), dtype=COUNT_DTYPES[dtype_name])
        loaded_targets = _load_sample_genome_matrices(
            conn,
            genome_idx=spec.genome_idx,
            sample_ids=block_ids.tolist(),
            matrix_length=spec.matrix_length,
            dtype_name=dtype_name,
        )
        target_matrices = np.stack(
            [loaded_targets.get(int(sample_idx), zero_matrix) for sample_idx in block_ids],
            axis=2,
        )
        return block_ids, block_names, zero_matrix, target_matrices

    async def drain_writer_results(force_one: bool = False, force_all: bool = False) -> None:
        nonlocal written_rows, completed_work
        while pending_write_futures:
            batch_future, batch_progress = pending_write_futures[0]
            if force_all or force_one:
                batch_rows = await batch_future
            elif not batch_future.done():
                break
            else:
                batch_rows = batch_future.result()
            pending_write_futures.popleft()
            written_rows += batch_rows
            for event in batch_progress:
                completed_work += int(event["delta"])
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "advance",
                            "completed": completed_work,
                            "total": total_work,
                            "anchor_name": str(event["anchor_name"]),
                            "genome": str(event["genome"]),
                            "scaffold": "",
                            "targets_completed": int(event["targets_completed"]),
                            "targets_total": int(event["targets_total"]),
                            "target_chunks": int(event["target_chunks"]),
                        }
                    )
            if force_one:
                force_one = False

    ordered_genomes = sorted(genomes, key=lambda item: item.genome)
    work_units = [
        (block_start, samples[block_start:block_start + global_block_size], spec)
        for block_start in range(0, len(samples), global_block_size)
        for spec in ordered_genomes
        if samples[block_start:block_start + global_block_size]
    ]
    target_queue_capacity = max(0, target_queue_size - 1)
    target_queue: deque[tuple[int, tuple[np.ndarray, list[str], np.ndarray, np.ndarray]]] = deque()
    next_target_prefetch_idx = 1
    pending_target_future: Optional[asyncio.Future[tuple[int, tuple[np.ndarray, list[str], np.ndarray, np.ndarray]]]] = None

    def submit_target_prefetch(current_unit_index: int) -> None:
        nonlocal next_target_prefetch_idx, pending_target_future
        if target_queue_capacity <= 0 or pending_target_future is not None:
            return
        if next_target_prefetch_idx <= current_unit_index:
            next_target_prefetch_idx = current_unit_index + 1
        if next_target_prefetch_idx >= len(work_units):
            return
        if len(target_queue) >= target_queue_capacity:
            return
        unit_index = next_target_prefetch_idx
        next_target_prefetch_idx += 1
        _prefetch_block_start, prefetch_block_rows, prefetch_spec = work_units[unit_index]
        pending_target_future = loop.run_in_executor(
            target_loader_executor,
            _load_target_prefetch_unit_for_torch,
            unit_index,
            matrix_db_file,
            prefetch_spec.genome_idx,
            prefetch_block_rows,
            prefetch_spec.matrix_length,
            dtype_name,
        )

    async def drain_target_prefetch(force: bool = False) -> None:
        nonlocal pending_target_future
        if pending_target_future is None:
            return
        if force:
            prefetched = await pending_target_future
        else:
            if not pending_target_future.done():
                return
            prefetched = pending_target_future.result()
        target_queue.append(prefetched)
        pending_target_future = None

    try:
        for unit_index, (block_start, block_rows, spec) in enumerate(work_units):
            await drain_target_prefetch(force=False)
            if target_queue and target_queue[0][0] == unit_index:
                block_ids, block_names, zero_matrix, target_matrices = target_queue.popleft()[1]
            else:
                block_ids, block_names, zero_matrix, target_matrices = load_target_block_sync(spec, block_rows)
            submit_target_prefetch(unit_index)
            completed_pairs_in_block = completed_pairs_by_genome.get(spec.genome_idx, set())
            external_rows = []
            external_missing_positions: dict[int, list[int]] = {}
            for anchor_idx, anchor_name in samples[:block_start]:
                missing_positions = [
                    pos
                    for pos, target_idx in enumerate(block_ids.tolist())
                    if (anchor_idx, int(target_idx)) not in completed_pairs_in_block
                ]
                if missing_positions:
                    external_rows.append((anchor_idx, anchor_name))
                    external_missing_positions[anchor_idx] = missing_positions

            internal_missing_positions: dict[int, list[int]] = {}
            for local_anchor_pos in range(len(block_ids) - 1):
                anchor_idx = int(block_ids[local_anchor_pos])
                missing_positions = [
                    later_pos - (local_anchor_pos + 1)
                    for later_pos in range(local_anchor_pos + 1, len(block_ids))
                    if (anchor_idx, int(block_ids[later_pos])) not in completed_pairs_in_block
                ]
                if missing_positions:
                    internal_missing_positions[local_anchor_pos] = missing_positions

            pair_count_for_block = sum(len(v) for v in external_missing_positions.values()) + sum(
                len(v) for v in internal_missing_positions.values()
            )
            if pair_count_for_block == 0:
                continue

            genome_name = spec.genome
            tile_targets, tile_size = _plan_chunk_sizes(
                vector_length=spec.matrix_length,
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
            target_torch = _prepare_torch_matrix(
                compute_backend=compute_backend,
                matrix=target_matrices,
                matrix_value_semantics=matrix_value_semantics,
            )
            processed_pairs_for_block = 0
            target_chunks += 1
            pending_device_results: list[dict[str, object]] = []
            pending_payloads: list[dict[str, object]] = []
            pending_progress: list[dict[str, object]] = []

            def submit_write_batch(
                batch_payloads: list[dict[str, object]],
                batch_progress: list[dict[str, object]],
            ) -> None:
                batch_pairs = sum(int(event["delta"]) for event in batch_progress)
                last_event = batch_progress[-1] if batch_progress else {
                    "anchor_name": "",
                    "genome": genome_name,
                    "targets_completed": 0,
                    "targets_total": 0,
                    "target_chunks": target_chunks,
                }
                log_context = None
                if emit_writer_logs and batch_progress:
                    log_context = {
                        "start_time": run_start_time,
                        "completed": completed_work + batch_pairs,
                        "total": total_work,
                        "batch_pairs": batch_pairs,
                        "anchor_name": str(last_event["anchor_name"]),
                        "genome": str(last_event["genome"]),
                        "targets_completed": int(last_event["targets_completed"]),
                        "targets_total": int(last_event["targets_total"]),
                        "target_chunks": int(last_event["target_chunks"]),
                    }
                pending_write_futures.append(
                    (
                        loop.run_in_executor(
                            writer_executor,
                            _write_matrix_compare_payload_batch,
                            output_file,
                            batch_payloads,
                            log_context,
                        ),
                        batch_progress,
                    )
                )

            async def flush_pending_device_results() -> None:
                nonlocal pending_device_results, pending_payloads
                if not pending_device_results:
                    return
                batch_device_results = pending_device_results
                pending_device_results = []
                combined_np = _download_torch_result_tensor_batch(
                    compute_backend=compute_backend,
                    totals_tensors=[item["total_positions"] for item in batch_device_results],
                    shared_tensors=[item["share_allele_pos"] for item in batch_device_results],
                )
                for result_idx, item in enumerate(batch_device_results):
                    valid_count = int(item["valid_count"])
                    pending_payloads.append(
                        {
                            "sample_1_idx": int(item["sample_1_idx"]),
                            "sample_1": str(item["sample_1"]),
                            "sample_2_idx": np.asarray(item["sample_2_idx"], dtype=np.int64),
                            "sample_2": list(item["sample_2"]),
                            "genome_idx": int(item["genome_idx"]),
                            "genome": str(item["genome"]),
                            "calculations": tuple(item["calculations"]),
                            "total_positions": combined_np[result_idx, 0, :valid_count].astype(np.int64, copy=True),
                            "share_allele_pos": combined_np[result_idx, 1, :valid_count].astype(np.int64, copy=True),
                            "max_consecutive_length": item["max_consecutive_length"],
                        }
                    )

            async def flush_pending_block_units() -> None:
                nonlocal pending_payloads, pending_progress, completed_work
                await flush_pending_device_results()
                if not pending_payloads and not pending_progress:
                    return
                batch_payloads = pending_payloads
                batch_progress = pending_progress
                pending_payloads = []
                pending_progress = []
                if batch_payloads:
                    if pending_write_futures:
                        await drain_writer_results(force_one=True)
                    submit_write_batch(batch_payloads, batch_progress)
                    return
                for event in batch_progress:
                    completed_work += int(event["delta"])
                if batch_progress and progress_callback is not None:
                    last_event = batch_progress[-1]
                    progress_callback(
                        {
                            "phase": "advance",
                            "completed": completed_work,
                            "total": total_work,
                            "anchor_name": str(last_event["anchor_name"]),
                            "genome": str(last_event["genome"]),
                            "scaffold": "",
                            "targets_completed": int(last_event["targets_completed"]),
                            "targets_total": int(last_event["targets_total"]),
                            "target_chunks": int(last_event["target_chunks"]),
                        }
                    )

            anchor_queue: deque[tuple[int, str, np.ndarray]] = deque()
            next_anchor_offset = 0
            pending_anchor_future: Optional[asyncio.Future[list[tuple[int, str, np.ndarray]]]] = None

            def submit_anchor_prefetch() -> None:
                nonlocal next_anchor_offset, pending_anchor_future
                if pending_anchor_future is not None:
                    return
                if next_anchor_offset >= len(external_rows):
                    return
                queue_deficit = anchor_queue_size - len(anchor_queue)
                if queue_deficit <= 0:
                    return
                batch_rows = external_rows[next_anchor_offset:next_anchor_offset + queue_deficit]
                if not batch_rows:
                    return
                next_anchor_offset += len(batch_rows)
                pending_anchor_future = loop.run_in_executor(
                    anchor_loader_executor,
                    _load_anchor_queue_batch_for_torch,
                    matrix_db_file,
                    spec.genome_idx,
                    batch_rows,
                    spec.matrix_length,
                    dtype_name,
                    zero_matrix,
                )

            async def drain_anchor_prefetch(force: bool = False) -> None:
                nonlocal pending_anchor_future
                if pending_anchor_future is None:
                    return
                if force:
                    prefetched = await pending_anchor_future
                else:
                    if not pending_anchor_future.done():
                        return
                    prefetched = pending_anchor_future.result()
                anchor_queue.extend(prefetched)
                pending_anchor_future = None

            submit_anchor_prefetch()
            await drain_anchor_prefetch(force=True)
            for external_pos in range(len(external_rows)):
                await drain_writer_results()
                if not anchor_queue:
                    await drain_anchor_prefetch(force=True)
                if not anchor_queue:
                    raise RuntimeError("Anchor queue unexpectedly empty during torch matrix compare.")
                anchor_idx, _anchor_name, anchor_matrix = anchor_queue.popleft()
                submit_anchor_prefetch()
                missing_positions = external_missing_positions[anchor_idx]
                anchor_torch = _prepare_torch_matrix(
                    compute_backend=compute_backend,
                    matrix=anchor_matrix,
                    matrix_value_semantics=matrix_value_semantics,
                )
                total_inc_torch, shared_inc_torch, ibs_inc = _compare_anchor_against_target_chunk_torch_device(
                    compute_backend=compute_backend,
                    anchor_torch=anchor_torch,
                    target_torch=target_torch,
                    vector_length=spec.matrix_length,
                    tile_size=tile_size,
                    matrix_value_semantics=matrix_value_semantics,
                    need_ibs="ibs" in calculations,
                )
                await drain_anchor_prefetch(force=False)
                pending_device_results.append(
                    {
                        "sample_1_idx": anchor_idx,
                        "sample_1": samples_by_id[anchor_idx],
                        "sample_2_idx": block_ids[missing_positions].astype(np.int64, copy=True),
                        "sample_2": [block_names[pos] for pos in missing_positions],
                        "genome_idx": spec.genome_idx,
                        "genome": genome_name,
                        "calculations": calculations,
                        "total_positions": total_inc_torch[missing_positions],
                        "share_allele_pos": shared_inc_torch[missing_positions],
                        "max_consecutive_length": None
                        if ibs_inc is None
                        else ibs_inc[missing_positions].astype(np.int64, copy=True),
                        "valid_count": len(missing_positions),
                    }
                )
                processed_pairs_for_block += len(missing_positions)
                pending_progress.append(
                    {
                        "delta": len(missing_positions),
                        "anchor_name": samples_by_id[anchor_idx],
                        "genome": spec.genome,
                        "targets_completed": external_pos + 1,
                        "targets_total": len(external_rows) + len(block_ids) - 1,
                        "target_chunks": target_chunks,
                    }
                )
                if len(pending_device_results) >= result_transfer_batch_size:
                    await flush_pending_device_results()
                if len(pending_progress) >= MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS:
                    await flush_pending_block_units()

            for local_anchor_pos in range(len(block_ids) - 1):
                await drain_writer_results()
                missing_positions = internal_missing_positions.get(local_anchor_pos, [])
                if not missing_positions:
                    continue
                total_inc_torch, shared_inc_torch, ibs_inc = _compare_anchor_against_target_chunk_torch_device(
                    compute_backend=compute_backend,
                    anchor_torch=target_torch[:, :, local_anchor_pos],
                    target_torch=target_torch[:, :, local_anchor_pos + 1:],
                    vector_length=spec.matrix_length,
                    tile_size=tile_size,
                    matrix_value_semantics=matrix_value_semantics,
                    need_ibs="ibs" in calculations,
                )
                later_ids = block_ids[local_anchor_pos + 1:]
                later_names = block_names[local_anchor_pos + 1:]
                pending_device_results.append(
                    {
                        "sample_1_idx": int(block_ids[local_anchor_pos]),
                        "sample_1": block_names[local_anchor_pos],
                        "sample_2_idx": later_ids[missing_positions].astype(np.int64, copy=True),
                        "sample_2": [later_names[pos] for pos in missing_positions],
                        "genome_idx": spec.genome_idx,
                        "genome": genome_name,
                        "calculations": calculations,
                        "total_positions": total_inc_torch[missing_positions],
                        "share_allele_pos": shared_inc_torch[missing_positions],
                        "max_consecutive_length": None
                        if ibs_inc is None
                        else ibs_inc[missing_positions].astype(np.int64, copy=True),
                        "valid_count": len(missing_positions),
                    }
                )
                processed_pairs_for_block += len(missing_positions)
                pending_progress.append(
                    {
                        "delta": len(missing_positions),
                        "anchor_name": block_names[local_anchor_pos],
                        "genome": spec.genome,
                        "targets_completed": local_anchor_pos + 1,
                        "targets_total": len(block_ids) - 1,
                        "target_chunks": target_chunks,
                    }
                )
                if len(pending_device_results) >= result_transfer_batch_size:
                    await flush_pending_device_results()
                if len(pending_progress) >= MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS:
                    await flush_pending_block_units()

            await flush_pending_block_units()
            await drain_writer_results()

            if processed_pairs_for_block != pair_count_for_block:
                raise RuntimeError(
                    "Torch block compare progress accounting drifted from the expected remaining pair count."
                )
            await drain_target_prefetch(force=False)
            submit_target_prefetch(unit_index)

        await drain_writer_results(force_all=True)
    finally:
        target_loader_executor.shutdown(wait=True)
        anchor_loader_executor.shutdown(wait=True)
        writer_executor.shutdown(wait=True)

    return MatrixCompareSummary(
        output_file=output_file,
        requested_pairs=requested_pairs,
        written_rows=written_rows,
        scaffold_count=len(genome_scaffolds),
        genome_count=len(genomes),
        anchor_groups=anchor_groups,
        target_chunks=target_chunks,
        min_cov=min_cov,
        genome_scope=genome_scope or metadata.get("genome_scope", "all"),
        backend=compute_backend.kind,
        device=compute_backend.device,
        memory_limit_gb=memory_limit_gb,
        position_tile_size=position_tile_size or default_tile_size,
    )


def _matrix_compare_reuse_target_chunks_torch(
    matrix_db_file: Path,
    conn: duckdb.DuckDBPyConnection,
    compare_conn: Optional[duckdb.DuckDBPyConnection],
    output_file: Path,
    min_cov: int,
    genome_scope: Optional[str],
    metadata: dict[str, str],
    samples: list[tuple[int, str]],
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
    completed_pairs_by_genome: dict[int, set[tuple[int, int]]],
    requested_pairs: int,
    remaining_work: int,
    memory_limit_bytes: int,
    memory_limit_gb: float,
    anchor_queue_size: int,
    target_queue_size: int,
    result_transfer_batch_size: int,
    loader_executor_kind: str,
    writer_executor_kind: str,
    position_tile_size: Optional[int],
    compute_backend: MatrixPairComputeBackend,
    matrix_value_semantics: str,
    calculations: tuple[str, ...],
    emit_writer_logs: bool,
    run_start_time: float,
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    return asyncio.run(
        _matrix_compare_reuse_target_chunks_torch_async(
            matrix_db_file=matrix_db_file,
            conn=conn,
            compare_conn=compare_conn,
            output_file=output_file,
            min_cov=min_cov,
            genome_scope=genome_scope,
            metadata=metadata,
            samples=samples,
            genomes=genomes,
            genome_scaffolds=genome_scaffolds,
            completed_pairs_by_genome=completed_pairs_by_genome,
            requested_pairs=requested_pairs,
            remaining_work=remaining_work,
            memory_limit_bytes=memory_limit_bytes,
            memory_limit_gb=memory_limit_gb,
            anchor_queue_size=anchor_queue_size,
            target_queue_size=target_queue_size,
            result_transfer_batch_size=result_transfer_batch_size,
            loader_executor_kind=loader_executor_kind,
            writer_executor_kind=writer_executor_kind,
            position_tile_size=position_tile_size,
            compute_backend=compute_backend,
            matrix_value_semantics=matrix_value_semantics,
            calculations=calculations,
            emit_writer_logs=emit_writer_logs,
            run_start_time=run_start_time,
            progress_callback=progress_callback,
        )
    )


def matrix_compare(
    matrix_db_file: Path,
    output_file: Path,
    min_cov: int = 5,
    genome: str = "all",
    memory_limit_gb: float = 16.0,
    anchor_queue_size: int = 1,
    target_queue_size: int = 1,
    result_transfer_batch_size: int = MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE_DEFAULT,
    loader_executor_kind: str = "thread",
    writer_executor_kind: str = "thread",
    position_tile_size: Optional[int] = None,
    backend: str = "numpy",
    calculate: Optional[str] = "all",
    emit_writer_logs: bool = False,
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    if min_cov < 1:
        raise ValueError("min_cov must be >= 1")
    if anchor_queue_size < 1:
        raise ValueError("anchor_queue_size must be >= 1")
    if target_queue_size < 1:
        raise ValueError("target_queue_size must be >= 1")
    if result_transfer_batch_size < 1:
        raise ValueError("result_transfer_batch_size must be >= 1")
    if loader_executor_kind not in MATRIX_IO_EXECUTOR_KINDS:
        raise ValueError(
            f"loader_executor_kind must be one of {', '.join(MATRIX_IO_EXECUTOR_KINDS)}"
        )
    if writer_executor_kind not in MATRIX_IO_EXECUTOR_KINDS:
        raise ValueError(
            f"writer_executor_kind must be one of {', '.join(MATRIX_IO_EXECUTOR_KINDS)}"
        )
    matrix_db_file = matrix_db_file.resolve()
    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    compute_backend = MatrixPairComputeBackend(backend)
    calculations = parse_matrix_calculations(calculate)
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    run_start_time = time.monotonic()
    conn = duckdb.connect(str(matrix_db_file), read_only=True)
    compare_conn: Optional[duckdb.DuckDBPyConnection] = None
    try:
        metadata = _load_matrix_db_metadata(conn)
        matrix_value_semantics = metadata.get("matrix_value_semantics")
        if matrix_value_semantics != FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
            raise ValueError(
                "This matrix DB uses an unsupported legacy storage layout. "
                "Rebuild it with 'zipstrain utilities build-matrix-db'."
            )
        samples = _load_matrix_db_samples(conn)
        genome_scope = None if genome == "all" else genome
        genomes = _load_matrix_db_genomes(conn, genome=genome_scope)
        if not genomes:
            raise ValueError(f"No scaffolds found for genome scope: {genome}")
        genome_scaffolds = _load_matrix_db_genome_scaffolds(conn, genome=genome_scope)
        compare_conn = _prepare_matrix_compare_db(
            output_file=output_file,
            matrix_metadata=metadata,
            samples=samples,
            genomes=genomes,
            genome_scope=genome_scope,
            calculations=calculations,
        )
        completed_pairs_by_genome, requested_pairs, total_work = _load_matrix_compare_resume_state(
            compare_conn,
            sample_count=len(samples),
            genome_ids=[spec.genome_idx for spec in genomes],
        )
        if total_work == 0:
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

        if compute_backend.kind == "torch":
            compare_conn.close()
            compare_conn = None
            return _matrix_compare_reuse_target_chunks_torch(
                matrix_db_file=matrix_db_file,
                conn=conn,
                compare_conn=compare_conn,
                output_file=output_file,
                min_cov=min_cov,
                genome_scope=genome_scope,
                metadata=metadata,
                samples=samples,
                genomes=genomes,
                genome_scaffolds=genome_scaffolds,
                completed_pairs_by_genome=completed_pairs_by_genome,
                requested_pairs=requested_pairs,
                remaining_work=total_work,
                memory_limit_bytes=memory_limit_bytes,
                memory_limit_gb=memory_limit_gb,
                anchor_queue_size=anchor_queue_size,
                target_queue_size=target_queue_size,
                result_transfer_batch_size=result_transfer_batch_size,
                loader_executor_kind=loader_executor_kind,
                writer_executor_kind=writer_executor_kind,
                position_tile_size=position_tile_size,
                compute_backend=compute_backend,
                matrix_value_semantics=matrix_value_semantics,
                calculations=calculations,
                emit_writer_logs=emit_writer_logs,
                run_start_time=run_start_time,
                progress_callback=progress_callback,
            )

        completed_work = 0
        written_rows = 0
        target_chunks = 0
        default_tile_size = 0
        anchor_groups = 0

        dtype_name = str(metadata.get("count_dtype", "uint16"))
        for anchor_offset, (sample_1_idx, sample_1_name) in enumerate(samples[:-1]):
            target_sample_rows = samples[anchor_offset + 1 :]
            if not target_sample_rows:
                continue
            anchor_groups += 1
            for spec in genomes:
                completed_pairs_for_genome = completed_pairs_by_genome.get(spec.genome_idx, set())
                remaining_target_rows = [
                    row
                    for row in target_sample_rows
                    if (sample_1_idx, row[0]) not in completed_pairs_for_genome
                ]
                if not remaining_target_rows:
                    continue
                target_ids_all = np.array([sample_idx for sample_idx, _sample_name in remaining_target_rows], dtype=np.int64)
                target_names_all = [sample_name for _sample_idx, sample_name in remaining_target_rows]
                zero_matrix = np.zeros((spec.matrix_length, 4), dtype=COUNT_DTYPES[dtype_name])
                anchor_loaded = _load_sample_genome_matrices(
                    conn,
                    genome_idx=spec.genome_idx,
                    sample_ids=[sample_1_idx],
                    matrix_length=spec.matrix_length,
                    dtype_name=dtype_name,
                )
                anchor_matrix = anchor_loaded.get(sample_1_idx, zero_matrix)
                pending_tables: list[pa.Table] = []
                pending_completed_rows: list[tuple[int, int, int]] = []
                pending_progress: list[dict[str, int]] = []

                def flush_pending_numpy_chunks() -> None:
                    nonlocal written_rows, completed_work
                    if not pending_tables and not pending_completed_rows and not pending_progress:
                        return
                    batch_rows = 0
                    batch_pairs = sum(progress_info["delta"] for progress_info in pending_progress)
                    if pending_tables:
                        compare_conn.commit()
                        compare_conn.execute("BEGIN")
                        try:
                            for table in pending_tables:
                                _insert_matrix_compare_result_table(compare_conn, table)
                                batch_rows += table.num_rows
                            _mark_completed_pair_genomes(compare_conn, pending_completed_rows)
                            compare_conn.execute("COMMIT")
                        except Exception:
                            compare_conn.execute("ROLLBACK")
                            raise
                        written_rows += batch_rows
                    elif pending_completed_rows:
                        compare_conn.commit()
                        compare_conn.execute("BEGIN")
                        try:
                            _mark_completed_pair_genomes(compare_conn, pending_completed_rows)
                            compare_conn.execute("COMMIT")
                        except Exception:
                            compare_conn.execute("ROLLBACK")
                            raise
                    for progress_info in pending_progress:
                        completed_work += progress_info["delta"]
                    if pending_progress:
                        last_progress = pending_progress[-1]
                        if emit_writer_logs:
                            _emit_matrix_compare_writer_log(
                                start_time=run_start_time,
                                completed=completed_work,
                                total=total_work,
                                batch_pairs=batch_pairs,
                                batch_rows=batch_rows,
                                anchor_name=sample_1_name,
                                genome=spec.genome,
                                targets_completed=last_progress["targets_completed"],
                                targets_total=len(target_ids_all),
                                target_chunks=target_chunks,
                            )
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": sample_1_name,
                                    "genome": spec.genome,
                                    "scaffold": "",
                                    "targets_completed": last_progress["targets_completed"],
                                    "targets_total": len(target_ids_all),
                                    "target_chunks": target_chunks,
                                }
                            )
                    pending_tables.clear()
                    pending_completed_rows.clear()
                    pending_progress.clear()

                target_offset = 0
                while target_offset < len(target_ids_all):
                    max_targets, tile_size = _plan_chunk_sizes(
                        vector_length=spec.matrix_length,
                        remaining_targets=len(target_ids_all) - target_offset,
                        dtype_name=dtype_name,
                        memory_limit_bytes=memory_limit_bytes,
                        backend_kind=compute_backend.kind,
                        position_tile_size=position_tile_size,
                    )
                    default_tile_size = tile_size
                    chunk_ids = target_ids_all[target_offset: target_offset + max_targets]
                    chunk_names = target_names_all[target_offset: target_offset + max_targets]
                    loaded_targets = _load_sample_genome_matrices(
                        conn,
                        genome_idx=spec.genome_idx,
                        sample_ids=chunk_ids.tolist(),
                        matrix_length=spec.matrix_length,
                        dtype_name=dtype_name,
                    )
                    target_matrices = np.stack(
                        [loaded_targets.get(int(sample_idx), zero_matrix) for sample_idx in chunk_ids],
                        axis=2,
                    )
                    totals_chunk = np.zeros(len(chunk_ids), dtype=np.int64)
                    shared_chunk = np.zeros(len(chunk_ids), dtype=np.int64)
                    ibs_max = np.zeros(len(chunk_ids), dtype=np.int64) if "ibs" in calculations else None
                    current_runs = np.zeros(len(chunk_ids), dtype=np.int64) if "ibs" in calculations else None
                    for tile_start in range(0, spec.matrix_length, tile_size):
                        tile_stop = min(spec.matrix_length, tile_start + tile_size)
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
                        totals_chunk += total_inc
                        shared_chunk += shared_inc
                        if "ibs" in calculations and current_runs is not None and ibs_max is not None:
                            current_runs, updated_max = _update_ibs_numpy(
                                shared_mask=shared_mask,
                                current_runs=current_runs,
                                max_runs=ibs_max,
                            )
                            ibs_max[:] = updated_max
                    mask = totals_chunk > 0
                    if mask.any():
                        pending_tables.append(
                            _make_arrow_table(
                                sample_1_idx=sample_1_idx,
                                sample_1=sample_1_name,
                                sample_2_idx=[int(chunk_ids[idx]) for idx, keep in enumerate(mask) if keep],
                                sample_2=[name for idx, name in enumerate(chunk_names) if mask[idx]],
                                genome_idx=spec.genome_idx,
                                genome=spec.genome,
                                calculations=calculations,
                                total_positions=totals_chunk[mask] if "ani" in calculations else None,
                                share_allele_pos=shared_chunk[mask] if "ani" in calculations else None,
                                max_consecutive_length=ibs_max[mask] if "ibs" in calculations else None,
                            )
                        )
                    pending_completed_rows.extend(
                        (sample_1_idx, int(target_idx), spec.genome_idx)
                        for target_idx in chunk_ids.tolist()
                    )
                    target_offset += max_targets
                    target_chunks += 1
                    pending_progress.append(
                        {
                            "delta": len(chunk_ids),
                            "targets_completed": target_offset,
                        }
                    )
                    if len(pending_progress) >= MATRIX_COMPARE_CHECKPOINT_BATCH_UNITS:
                        flush_pending_numpy_chunks()

                flush_pending_numpy_chunks()

        return MatrixCompareSummary(
            output_file=output_file,
            requested_pairs=requested_pairs,
            written_rows=written_rows,
            scaffold_count=len(genome_scaffolds),
            genome_count=len(genomes),
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
        if compare_conn is not None:
            compare_conn.close()


def export_matrix_compare_parquet(
    matrix_compare_db_file: Path,
    output_file: Path,
) -> Path:
    matrix_compare_db_file = Path(matrix_compare_db_file).resolve()
    output_file = Path(output_file).resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    conn = duckdb.connect(str(matrix_compare_db_file), read_only=True)
    try:
        compare_metadata = {
            str(k): str(v)
            for k, v in conn.execute("SELECT key, value FROM matrix_compare_metadata").fetchall()
        }
        calculate = compare_metadata.get("calculate", "ani")
        columns = matrix_metric_output_columns(calculate)
        query = f"""
            SELECT {', '.join(columns)}
            FROM matrix_compare_results
            ORDER BY sample_idx_1, sample_idx_2, genome_idx
        """
        conn.execute(
            f"COPY ({query}) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
            [str(output_file)],
        )
    finally:
        conn.close()
    return output_file
