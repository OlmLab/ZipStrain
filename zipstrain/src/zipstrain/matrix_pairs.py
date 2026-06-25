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
CURRENT_MATRIX_HDF5_LAYOUT = "per_genome_sample_major_dense_matrix_hdf5"
CURRENT_MATRIX_HDF5_SPARSE_LAYOUT = "per_genome_sample_major_sparse_indices_matrix_hdf5"
MATRIX_HDF5_FILE_VERSION = "1"
MATRIX_PAIR_BACKENDS = ("numpy", "torch", "torch-cpu", "torch-cuda", "torch-mps")
MATRIX_IO_EXECUTOR_KINDS = ("thread", "process")
MATRIX_HDF5_SUFFIXES = (".h5", ".hdf5", ".hd5")
HDF5_FILE_SIGNATURE = b"\x89HDF\r\n\x1a\n"
MATRIX_HDF5_CONTRACT_GENOMES_GROUP = "contract_genomes"
MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP = "contract_genome_scaffolds"
BuildProgressCallback = Callable[[dict[str, object]], None]
CompareProgressCallback = Callable[[dict[str, object]], None]
MATRIX_COMPARISON_SUPPORTED_CALCULATIONS = ("ani", "ibs", "gene")
MATRIX_COMPARISON_CALCULATION_ALIASES = {
    "ani": "ani",
    "popani": "ani",
    "ibs": "ibs",
    "max_block": "ibs",
    "max_consecutive_length": "ibs",
    "gene": "gene",
    "genes": "gene",
    "gene_ani": "gene",
}
MATRIX_BUILD_FIXED_HEADROOM_BYTES = 128 * 1024 * 1024
MATRIX_BUILD_TEMP_BYTES_PER_POSITION = 32
MATRIX_BUILD_MIN_DUCKDB_MEMORY_BYTES = 256 * 1024 * 1024
MATRIX_BUILD_MIN_COMMIT_BATCH_BYTES = 64 * 1024 * 1024
MATRIX_BUILD_DUCKDB_MEMORY_FRACTION = 0.25
MATRIX_COMPARE_CHECKPOINT_BATCH_UNITS = 4
MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS = 16
MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE_DEFAULT = 1
MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT = 128.0


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

def parse_matrix_calculations(
    calculate: Optional[str] = None,
    *,
    include_gene_from_all: bool = False,
) -> tuple[str, ...]:
    """Normalize matrix compare metric selections into a stable tuple.

    Parameters
    ----------
    calculate:
        User-facing metric selection string such as ``ani``, ``ani+ibs``,
        ``+gene``, or ``all``. Commas and plus signs are both accepted as
        separators.
    include_gene_from_all:
        When ``True``, interpret ``all`` as ``ani+ibs+gene``. This is used for
        matrix stores that already contain gene annotations.

    Returns
    -------
    tuple[str, ...]
        A tuple ordered according to ZipStrain's supported matrix metrics.
    """
    raw_tokens: list[str] = []
    if calculate is None:
        raw_tokens = ["all"]
    elif isinstance(calculate, str):
        for plus_part in calculate.split("+"):
            for comma_part in plus_part.split(","):
                token = comma_part.strip().lower()
                if token:
                    raw_tokens.append(token)
    else:
        raw_tokens = [str(calculate).strip().lower()]

    if not raw_tokens:
        raw_tokens = ["all"]

    normalized: set[str] = set()
    if "all" in raw_tokens:
        normalized.update(("ani", "ibs"))
        if include_gene_from_all:
            normalized.add("gene")

    for token in raw_tokens:
        if token == "all":
            continue
        mapped = MATRIX_COMPARISON_CALCULATION_ALIASES.get(token)
        if mapped is None:
            supported = "all," + ",".join(MATRIX_COMPARISON_SUPPORTED_CALCULATIONS)
            raise ValueError(f"Unsupported matrix calculation '{token}'. Supported values: {supported}")
        normalized.add(mapped)

    if "gene" in normalized:
        normalized.add("ani")

    supported = tuple(metric for metric in MATRIX_COMPARISON_SUPPORTED_CALCULATIONS if metric in normalized)
    if supported:
        return supported
    raise ValueError("Matrix compare requires at least one supported calculation.")


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


def matrix_compare_gene_result_db_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("sample_idx_1", pa.int64()),
            pa.field("sample_idx_2", pa.int64()),
            pa.field("sample_1", pa.string()),
            pa.field("sample_2", pa.string()),
            pa.field("genome_idx", pa.int64()),
            pa.field("genome", pa.string()),
            pa.field("gene", pa.string()),
            pa.field("gene_pop_ani", pa.float64()),
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
class GeneRangeSpec:
    gene_idx: int
    gene: str
    genome_idx: int
    genome: str
    chrom: str
    axis_start: int
    axis_end: int


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
class MatrixHdf5Summary:
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
    ignored_genome_count: int
    count_dtype: str
    genome_scope: str
    memory_limit_gb: float


@dataclass(frozen=True)
class MatrixHdf5AppendSummary:
    output_file: Path
    appended_profile_files: int
    appended_sample_count: int
    total_sample_count: int
    scaffold_count: int
    stored_rows: int
    ignored_genome_count: int
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


@dataclass(frozen=True)
class MatrixBuildMemoryPlan:
    duckdb_memory_limit_bytes: int
    commit_batch_bytes: int
    estimated_python_peak_bytes: int
    limiting_scaffold: str


class TorchBackendMissingError(ImportError):
    pass


class Hdf5BackendMissingError(ImportError):
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


def _import_h5py():
    try:
        return importlib.import_module("h5py")
    except ImportError as exc:
        raise Hdf5BackendMissingError(
            "HDF5 matrix input requested but 'h5py' is not installed. "
            'Install with: pip install "zipstrain[matrix]".'
        ) from exc


def _is_hdf5_matrix_store(matrix_file: Path) -> bool:
    if matrix_file.exists():
        try:
            with matrix_file.open("rb") as fh:
                return fh.read(len(HDF5_FILE_SIGNATURE)) == HDF5_FILE_SIGNATURE
        except OSError:
            return False
    return matrix_file.suffix.lower() in MATRIX_HDF5_SUFFIXES


def _hdf5_string_dtype(h5py_module):
    return h5py_module.string_dtype(encoding="utf-8")


def _write_hdf5_string_dataset(
    group,
    name: str,
    values: list[str],
    *,
    h5py_module,
    maxshape=None,
    chunks=None,
) -> None:
    kwargs = {
        "data": np.asarray(values, dtype=object),
        "dtype": _hdf5_string_dtype(h5py_module),
    }
    if maxshape is not None:
        kwargs["maxshape"] = maxshape
    if chunks is not None:
        kwargs["chunks"] = chunks
    group.create_dataset(name, **kwargs)


def _read_hdf5_string_dataset(dataset) -> list[str]:
    if hasattr(dataset, "asstr"):
        return [str(value) for value in dataset.asstr()[...].tolist()]
    raw = dataset[...]
    return [
        value.decode("utf-8") if isinstance(value, (bytes, bytearray, np.bytes_)) else str(value)
        for value in raw.tolist()
    ]


def _write_hdf5_genomes_group(
    h5_file,
    group_name: str,
    genomes: list[GenomeSpec],
    *,
    h5py_module,
) -> None:
    genomes_group = h5_file.create_group(group_name)
    genomes_group.create_dataset(
        "genome_idx",
        data=np.asarray([spec.genome_idx for spec in genomes], dtype=np.int64),
    )
    _write_hdf5_string_dataset(
        genomes_group,
        "genome",
        [spec.genome for spec in genomes],
        h5py_module=h5py_module,
    )
    genomes_group.create_dataset(
        "matrix_length",
        data=np.asarray([spec.matrix_length for spec in genomes], dtype=np.int64),
    )
    genomes_group.create_dataset(
        "true_length",
        data=np.asarray([spec.true_length for spec in genomes], dtype=np.int64),
    )
    genomes_group.create_dataset(
        "scaffold_count",
        data=np.asarray([spec.scaffold_count for spec in genomes], dtype=np.int64),
    )


def _write_hdf5_genome_scaffolds_group(
    h5_file,
    group_name: str,
    genome_scaffolds: list[GenomeScaffoldOffset],
    *,
    h5py_module,
) -> None:
    scaffold_group = h5_file.create_group(group_name)
    scaffold_group.create_dataset(
        "genome_idx",
        data=np.asarray([spec.genome_idx for spec in genome_scaffolds], dtype=np.int64),
    )
    scaffold_group.create_dataset(
        "scaffold_ordinal",
        data=np.asarray([spec.scaffold_ordinal for spec in genome_scaffolds], dtype=np.int64),
    )
    _write_hdf5_string_dataset(
        scaffold_group,
        "genome",
        [spec.genome for spec in genome_scaffolds],
        h5py_module=h5py_module,
    )
    _write_hdf5_string_dataset(
        scaffold_group,
        "chrom",
        [spec.chrom for spec in genome_scaffolds],
        h5py_module=h5py_module,
    )
    for field_name in ("axis_start", "axis_end", "index_base", "vector_length", "min_pos", "max_pos"):
        scaffold_group.create_dataset(
            field_name,
            data=np.asarray([getattr(spec, field_name) for spec in genome_scaffolds], dtype=np.int64),
        )


def _write_hdf5_gene_ranges_group(
    h5_file,
    group_name: str,
    gene_ranges: list[GeneRangeSpec],
    *,
    h5py_module,
) -> None:
    genes_group = h5_file.create_group(group_name)
    genes_group.create_dataset(
        "gene_idx",
        data=np.asarray([spec.gene_idx for spec in gene_ranges], dtype=np.int64),
    )
    genes_group.create_dataset(
        "genome_idx",
        data=np.asarray([spec.genome_idx for spec in gene_ranges], dtype=np.int64),
    )
    _write_hdf5_string_dataset(
        genes_group,
        "genome",
        [spec.genome for spec in gene_ranges],
        h5py_module=h5py_module,
    )
    _write_hdf5_string_dataset(
        genes_group,
        "chrom",
        [spec.chrom for spec in gene_ranges],
        h5py_module=h5py_module,
    )
    _write_hdf5_string_dataset(
        genes_group,
        "gene",
        [spec.gene for spec in gene_ranges],
        h5py_module=h5py_module,
    )
    genes_group.create_dataset(
        "axis_start",
        data=np.asarray([spec.axis_start for spec in gene_ranges], dtype=np.int64),
    )
    genes_group.create_dataset(
        "axis_end",
        data=np.asarray([spec.axis_end for spec in gene_ranges], dtype=np.int64),
    )


def _hdf5_matrix_dataset_path(genome_idx: int) -> str:
    return f"matrices/{int(genome_idx)}"


def _matrix_hdf5_sample_axis_chunk_length(sample_count: int) -> int:
    return max(1, min(int(sample_count), 1024))


def _matrix_hdf5_chunk_sample_count(
    *,
    sample_count: int,
    matrix_length: int,
    dtype_name: str,
    target_batch_mb: float,
) -> int:
    if target_batch_mb <= 0:
        raise ValueError("target_batch_mb must be > 0")
    dtype = COUNT_DTYPES[dtype_name]
    per_sample_bytes = max(1, matrix_length * 4 * np.dtype(dtype).itemsize)
    target_batch_bytes = int(target_batch_mb * (1024 ** 2))
    return max(1, min(sample_count, int(target_batch_bytes // per_sample_bytes) or 1))


def _matrix_hdf5_layout_is_sparse(layout: str) -> bool:
    return str(layout) == CURRENT_MATRIX_HDF5_SPARSE_LAYOUT


def _matrix_hdf5_supported_layouts() -> tuple[str, ...]:
    return (CURRENT_MATRIX_HDF5_LAYOUT, CURRENT_MATRIX_HDF5_SPARSE_LAYOUT)


def _matrix_hdf5_sparse_indices_chunk_length(matrix_length: int) -> int:
    return max(1024, min(max(matrix_length * 4, 1), 1_048_576))


def _dense_matrix_to_sparse_flat_indices(matrix: np.ndarray) -> np.ndarray:
    return np.flatnonzero(matrix.reshape(-1)).astype(np.int64, copy=False)


def _create_sparse_hdf5_genome_store(
    matrices_group,
    *,
    genome_idx: int,
    sample_count: int,
    matrix_length: int,
):
    group = matrices_group.create_group(str(genome_idx))
    group.create_dataset(
        "indptr",
        shape=(sample_count + 1,),
        dtype=np.int64,
        chunks=(_matrix_hdf5_sample_axis_chunk_length(sample_count + 1),),
        maxshape=(None,),
        fillvalue=0,
    )
    group.create_dataset(
        "indices",
        shape=(0,),
        dtype=np.int64,
        chunks=(_matrix_hdf5_sparse_indices_chunk_length(matrix_length),),
        maxshape=(None,),
        fillvalue=0,
    )
    return group


def _append_sparse_hdf5_matrix_row(
    *,
    indptr_dataset,
    indices_dataset,
    sample_row: int,
    flat_indices: np.ndarray,
    current_nnz: int,
) -> int:
    flat_indices = np.asarray(flat_indices, dtype=np.int64)
    next_nnz = current_nnz + int(flat_indices.size)
    if flat_indices.size > 0:
        indices_dataset.resize((next_nnz,))
        indices_dataset[current_nnz:next_nnz] = flat_indices
    indptr_dataset[int(sample_row) + 1] = next_nnz
    return next_nnz


def _load_dense_rows_from_sparse_hdf5(
    *,
    indptr_dataset,
    indices_dataset,
    sample_rows: np.ndarray,
    matrix_length: int,
    numpy_dtype,
) -> np.ndarray:
    sample_rows = np.asarray(sample_rows, dtype=np.int64)
    flat_width = int(matrix_length) * 4
    dense = np.zeros((len(sample_rows), flat_width), dtype=numpy_dtype)
    for out_idx, sample_row in enumerate(sample_rows.tolist()):
        start = int(indptr_dataset[int(sample_row)])
        stop = int(indptr_dataset[int(sample_row) + 1])
        if stop <= start:
            continue
        row_indices = np.asarray(indices_dataset[start:stop], dtype=np.int64)
        if row_indices.size > 0:
            dense[out_idx, row_indices] = 1
    return dense.reshape(len(sample_rows), int(matrix_length), 4)


def _matrix_hdf5_store_is_append_resizable(matrix_hdf5_file: Path) -> bool:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        try:
            sample_idx_ds = h5_file["samples"]["sample_idx"]
            sample_name_ds = h5_file["samples"]["sample_name"]
            if sample_idx_ds.maxshape is None or sample_name_ds.maxshape is None:
                return False
            if sample_idx_ds.maxshape[0] is not None:
                return False
            if sample_name_ds.maxshape[0] is not None:
                return False
            matrices_group = h5_file["matrices"]
            layout = str(h5_file["metadata"].attrs.get("layout", CURRENT_MATRIX_HDF5_LAYOUT))
            for node in matrices_group.values():
                if _matrix_hdf5_layout_is_sparse(layout):
                    indptr_ds = node["indptr"]
                    indices_ds = node["indices"]
                    if indptr_ds.maxshape is None or indptr_ds.maxshape[0] is not None:
                        return False
                    if indices_ds.maxshape is None or indices_ds.maxshape[0] is not None:
                        return False
                else:
                    if node.maxshape is None or node.maxshape[0] is not None:
                        return False
        except Exception:
            return False
    return True

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


def _read_stb_mapping(stb_file: Path) -> pl.DataFrame:
    stb_path = Path(stb_file)
    if not stb_path.exists():
        raise FileNotFoundError(f"STB file does not exist: {stb_path}")
    if not stb_path.is_file():
        raise FileNotFoundError(f"STB path is not a file: {stb_path}")

    mapping = (
        pl.scan_csv(stb_path, separator="\t", has_header=False)
        .select(
            pl.col("column_1").cast(pl.Utf8).alias("scaffold"),
            pl.col("column_2").cast(pl.Utf8).alias("genome"),
        )
        .collect(engine="streaming")
    )
    duplicate_scaffolds = (
        mapping.group_by("scaffold")
        .agg(pl.col("genome").n_unique().alias("genome_count"))
        .filter(pl.col("genome_count") > 1)
        .get_column("scaffold")
        .to_list()
    )
    if duplicate_scaffolds:
        raise ValueError(
            "STB file must map each scaffold to exactly one genome. "
            "Duplicate scaffold mappings found for: "
            + ", ".join(sorted(str(name) for name in duplicate_scaffolds))
        )
    return mapping.unique(["scaffold", "genome"])


def _collect_scaffold_specs_from_bed_and_stb(
    *,
    bed_file: Path,
    stb_file: Path,
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

    stb_mapping = _read_stb_mapping(Path(stb_file))
    if genome is not None:
        stb_mapping = stb_mapping.filter(pl.col("genome") == genome)
    if stb_mapping.height == 0:
        detail = genome if genome is not None else "all"
        raise ValueError(f"No STB rows found for genome scope: {detail}")

    contract_df = stb_mapping.join(
        bed_spans,
        left_on="scaffold",
        right_on="chrom",
        how="left",
    )
    missing_scaffolds = (
        contract_df
        .filter(pl.col("max_end").is_null() | pl.col("min_start").is_null())
        .get_column("scaffold")
        .to_list()
    )
    if missing_scaffolds:
        raise ValueError(
            "BED file is missing scaffold intervals for the following STB scaffolds: "
            + ", ".join(sorted(str(name) for name in missing_scaffolds))
        )

    specs: list[ScaffoldSpec] = []
    for scaffold_idx, row in enumerate(contract_df.sort(["genome", "scaffold"]).iter_rows(named=True)):
        min_start = int(row["min_start"])
        max_end = int(row["max_end"])
        index_base = min_start + 1
        vector_length = max_end - min_start
        if vector_length <= 0:
            raise ValueError(
                f"Invalid BED span for scaffold {row['scaffold']}: start={min_start}, end={max_end}"
            )
        specs.append(
            ScaffoldSpec(
                scaffold_idx=scaffold_idx,
                genome=str(row["genome"]),
                chrom=str(row["scaffold"]),
                index_base=index_base,
                vector_length=vector_length,
                min_pos=index_base,
                max_pos=max_end,
            )
        )
    return specs


def _collect_profile_genome_names(
    profile_paths: list[Path],
    genome: Optional[str] = None,
) -> set[str]:
    scope = pl.lit(True)
    if genome is not None:
        scope = scope & (pl.col("genome") == genome)
    genomes = (
        pl.scan_parquet([str(path) for path in profile_paths])
        .filter(scope)
        .select(pl.col("genome").cast(pl.Utf8))
        .unique()
        .collect(engine="streaming")
        .get_column("genome")
        .to_list()
    )
    return {str(genome_name) for genome_name in genomes}


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


def _subset_contract_to_genomes(
    *,
    contract_genomes: list[GenomeSpec],
    contract_genome_scaffolds: list[GenomeScaffoldOffset],
    genome_names: set[str],
) -> tuple[list[GenomeSpec], list[GenomeScaffoldOffset]]:
    return (
        [spec for spec in contract_genomes if spec.genome in genome_names],
        [offset for offset in contract_genome_scaffolds if offset.genome in genome_names],
    )


def _collect_gene_range_specs(
    *,
    gene_range_table: Path,
    genome_scaffolds: list[GenomeScaffoldOffset],
) -> list[GeneRangeSpec]:
    gene_range_path = Path(gene_range_table)
    if not gene_range_path.exists():
        raise FileNotFoundError(f"Gene range table does not exist: {gene_range_path}")
    if not gene_range_path.is_file():
        raise FileNotFoundError(f"Gene range table path is not a file: {gene_range_path}")

    scaffold_offsets_by_chrom: dict[str, GenomeScaffoldOffset] = {}
    duplicate_scaffolds: set[str] = set()
    for offset in genome_scaffolds:
        if offset.chrom in scaffold_offsets_by_chrom:
            duplicate_scaffolds.add(offset.chrom)
            continue
        scaffold_offsets_by_chrom[offset.chrom] = offset
    if duplicate_scaffolds:
        raise ValueError(
            "Gene range table requires unique scaffold names across the selected matrix scope. "
            "Duplicate scaffold names found: " + ", ".join(sorted(duplicate_scaffolds))
        )

    gene_ranges = (
        pl.scan_csv(
            gene_range_path,
            has_header=False,
            separator="\t",
        )
        .rename(
            {
                "column_1": "gene",
                "column_2": "scaffold",
                "column_3": "start",
                "column_4": "end",
            }
        )
        .select(
            pl.col("gene").cast(pl.Utf8),
            pl.col("scaffold").cast(pl.Utf8),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64),
        )
        .collect(engine="streaming")
    )

    specs: list[GeneRangeSpec] = []
    for row in gene_ranges.iter_rows(named=True):
        chrom = str(row["scaffold"])
        offset = scaffold_offsets_by_chrom.get(chrom)
        if offset is None:
            continue
        start = int(row["start"])
        end = int(row["end"])
        if end < start:
            raise ValueError(
                f"Gene {row['gene']} on scaffold {chrom} has end < start ({end} < {start})."
            )
        if start < offset.min_pos or end > offset.max_pos:
            raise ValueError(
                f"Gene {row['gene']} on scaffold {chrom} falls outside the matrix coordinate range: "
                f"gene={start}-{end}, matrix={offset.min_pos}-{offset.max_pos}"
            )
        axis_start = offset.axis_start + (start - offset.index_base)
        axis_end = offset.axis_start + (end - offset.index_base)
        specs.append(
            GeneRangeSpec(
                gene_idx=-1,
                gene=str(row["gene"]),
                genome_idx=offset.genome_idx,
                genome=offset.genome,
                chrom=chrom,
                axis_start=int(axis_start),
                axis_end=int(axis_end),
            )
        )
    sorted_specs = sorted(
        specs,
        key=lambda item: (item.genome_idx, item.axis_start, item.axis_end, item.gene),
    )
    return [
        GeneRangeSpec(
            gene_idx=idx,
            gene=spec.gene,
            genome_idx=spec.genome_idx,
            genome=spec.genome,
            chrom=spec.chrom,
            axis_start=spec.axis_start,
            axis_end=spec.axis_end,
        )
        for idx, spec in enumerate(sorted_specs)
    ]


def _expand_scaffold_specs_with_gene_ranges(
    *,
    scaffolds: list[ScaffoldSpec],
    gene_range_table: Path,
) -> list[ScaffoldSpec]:
    gene_range_path = Path(gene_range_table)
    if not gene_range_path.exists():
        raise FileNotFoundError(f"Gene range table does not exist: {gene_range_path}")
    if not gene_range_path.is_file():
        raise FileNotFoundError(f"Gene range table path is not a file: {gene_range_path}")

    scaffolds_by_chrom: dict[str, ScaffoldSpec] = {}
    duplicate_scaffolds: set[str] = set()
    for spec in scaffolds:
        if spec.chrom in scaffolds_by_chrom:
            duplicate_scaffolds.add(spec.chrom)
            continue
        scaffolds_by_chrom[spec.chrom] = spec
    if duplicate_scaffolds:
        raise ValueError(
            "Gene range table requires unique scaffold names across the selected matrix scope. "
            "Duplicate scaffold names found: " + ", ".join(sorted(duplicate_scaffolds))
        )

    gene_bounds = (
        pl.scan_csv(
            gene_range_path,
            has_header=False,
            separator="\t",
        )
        .rename(
            {
                "column_1": "gene",
                "column_2": "scaffold",
                "column_3": "start",
                "column_4": "end",
            }
        )
        .select(
            pl.col("scaffold").cast(pl.Utf8),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64),
        )
        .group_by("scaffold")
        .agg(
            pl.col("start").min().alias("min_start"),
            pl.col("end").max().alias("max_end"),
        )
        .collect(engine="streaming")
    )

    expanded_by_chrom: dict[str, ScaffoldSpec] = {spec.chrom: spec for spec in scaffolds}
    for row in gene_bounds.iter_rows(named=True):
        chrom = str(row["scaffold"])
        spec = scaffolds_by_chrom.get(chrom)
        if spec is None:
            continue
        min_start = int(row["min_start"])
        max_end = int(row["max_end"])
        if max_end < min_start:
            raise ValueError(
                f"Gene range table has end < start on scaffold {chrom} ({max_end} < {min_start})."
            )
        new_min = min(spec.min_pos, min_start)
        new_max = max(spec.max_pos, max_end)
        expanded_by_chrom[chrom] = ScaffoldSpec(
            scaffold_idx=spec.scaffold_idx,
            genome=spec.genome,
            chrom=spec.chrom,
            index_base=new_min,
            vector_length=new_max - new_min + 1,
            min_pos=new_min,
            max_pos=new_max,
        )
    return [expanded_by_chrom[spec.chrom] for spec in scaffolds]


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
    conn.execute(
        """
        CREATE TABLE matrix_compare_gene_results (
          sample_idx_1 INTEGER NOT NULL,
          sample_idx_2 INTEGER NOT NULL,
          sample_1 VARCHAR NOT NULL,
          sample_2 VARCHAR NOT NULL,
          genome_idx INTEGER NOT NULL,
          genome VARCHAR NOT NULL,
          gene VARCHAR NOT NULL,
          gene_pop_ani DOUBLE NOT NULL
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
        ("matrix_input_format", "hdf5"),
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


def _ensure_matrix_compare_gene_results_table(
    compare_conn: duckdb.DuckDBPyConnection,
) -> None:
    if _matrix_compare_table_exists(compare_conn, "matrix_compare_gene_results"):
        return
    compare_conn.execute(
        """
        CREATE TABLE matrix_compare_gene_results (
          sample_idx_1 INTEGER NOT NULL,
          sample_idx_2 INTEGER NOT NULL,
          sample_1 VARCHAR NOT NULL,
          sample_2 VARCHAR NOT NULL,
          genome_idx INTEGER NOT NULL,
          genome VARCHAR NOT NULL,
          gene VARCHAR NOT NULL,
          gene_pop_ani DOUBLE NOT NULL
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
        if "matrix_input_format" not in compare_metadata:
            compare_conn.execute(
                "INSERT OR REPLACE INTO matrix_compare_metadata VALUES (?, ?)",
                ["matrix_input_format", "hdf5"],
            )
            compare_metadata["matrix_input_format"] = "hdf5"
        _validate_matrix_compare_db_metadata(
            compare_metadata=compare_metadata,
            matrix_metadata=matrix_metadata,
            genome_scope=genome_scope,
            calculations=calculations,
        )
        _ensure_matrix_compare_completed_table(compare_conn)
        _ensure_matrix_compare_gene_results_table(compare_conn)

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


def _insert_matrix_compare_gene_result_table(
    compare_conn: duckdb.DuckDBPyConnection,
    table: pa.Table,
) -> None:
    if table.num_rows == 0:
        return
    compare_conn.register("_matrix_compare_gene_result_batch", table)
    try:
        compare_conn.execute(
            """
            INSERT INTO matrix_compare_gene_results
            SELECT
              sample_idx_1,
              sample_idx_2,
              sample_1,
              sample_2,
              genome_idx,
              genome,
              gene,
              gene_pop_ani
            FROM _matrix_compare_gene_result_batch
            """
        )
    finally:
        compare_conn.unregister("_matrix_compare_gene_result_batch")


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
                if table is not None:
                    _insert_matrix_compare_result_table(write_conn, table)
                    batch_rows += table.num_rows
                gene_table = _make_gene_arrow_table_from_compare_payload(payload)
                if gene_table is not None:
                    _insert_matrix_compare_gene_result_table(write_conn, gene_table)
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


def _validate_matrix_hdf5_appendable(metadata: dict[str, str]) -> tuple[str, str, str]:
    layout = metadata.get("layout", "")
    if layout not in _matrix_hdf5_supported_layouts():
        raise ValueError(
            "Matrix store layout is not append-compatible with this builder. "
            "Rebuild the matrix store with the current builder."
        )
    semantics = metadata.get("matrix_value_semantics", "")
    if semantics != FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
        raise ValueError(
            "Matrix store value semantics are incompatible with append. "
            "Rebuild the matrix store with the current builder."
        )
    min_cov = metadata.get("coverage_filter_min_cov")
    if min_cov != str(MATRIX_BUILD_MIN_COV):
        raise ValueError(
            "Matrix store coverage filter is incompatible with append. "
            "Rebuild the matrix store with the current builder."
        )
    count_dtype = metadata.get("count_dtype", "")
    if count_dtype not in COUNT_DTYPES:
        raise ValueError(
            f"Matrix store count dtype '{count_dtype}' is not supported for append."
        )
    genome_scope = metadata.get("genome_scope", "all")
    return count_dtype, genome_scope, str(layout)


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
    genome_scope: Optional[str] = None,
    accepted_genomes: Optional[set[str]] = None,
) -> None:
    genome_names = {spec.genome for spec in genomes}
    offsets_by_key = {(offset.genome, offset.chrom): offset for offset in genome_scaffolds}
    scope = pl.lit(True)
    if genome_scope is not None:
        scope = scope & (pl.col("genome") == genome_scope)
    if accepted_genomes is not None:
        scope = scope & pl.col("genome").is_in(sorted(accepted_genomes))
    grouped = (
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
            f"Profile {profile_path.name} does not match the existing matrix store contract: "
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
    stb_file: Optional[Path] = None,
    progress_callback: Optional[BuildProgressCallback] = None,
    sparse: bool = False,
) -> MatrixDbSummary:
    summary = build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=output_file,
        genome=genome,
        count_dtype=count_dtype,
        memory_limit_gb=memory_limit_gb,
        bed_file=bed_file,
        stb_file=stb_file,
        progress_callback=progress_callback,
        sparse=sparse,
        export_batch_mb=max(
            MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
            (commit_batch_gb or 0.0) * 1024.0,
        ),
    )
    return MatrixDbSummary(
        output_file=summary.output_file,
        profile_files=summary.profile_files,
        sample_count=summary.sample_count,
        scaffold_count=summary.scaffold_count,
        stored_rows=summary.stored_rows,
        count_dtype=summary.count_dtype,
        genome_scope=summary.genome_scope,
        memory_limit_gb=summary.memory_limit_gb,
    )


def build_matrix_hdf5(
    profile_dir: Path,
    output_file: Path,
    genome: str = "all",
    count_dtype: str = "uint16",
    memory_limit_gb: float = 16.0,
    bed_file: Optional[Path] = None,
    stb_file: Optional[Path] = None,
    gene_range_table: Optional[Path] = None,
    progress_callback: Optional[BuildProgressCallback] = None,
    export_batch_mb: float = MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
    sparse: bool = False,
) -> MatrixHdf5Summary:
    """Build the current HDF5-backed matrix store from classic profile parquets.

    The resulting matrix store keeps one whole-genome matrix per genome and is
    intended for repeated, resumable matrix comparisons.
    """
    if count_dtype not in COUNT_DTYPES:
        raise ValueError(f"Unsupported count dtype '{count_dtype}'. Choose one of {', '.join(COUNT_DTYPES)}.")
    if export_batch_mb <= 0:
        raise ValueError("export_batch_mb must be > 0")
    if bed_file is None or stb_file is None:
        raise ValueError("build_matrix_hdf5 requires both bed_file and stb_file.")
    profile_paths = discover_profile_parquets(profile_dir)
    genome_scope = None if genome == "all" else genome
    contract_scaffolds = _collect_scaffold_specs_from_bed_and_stb(
        bed_file=Path(bed_file),
        stb_file=Path(stb_file),
        genome=genome_scope,
    )
    if gene_range_table is not None:
        contract_scaffolds = _expand_scaffold_specs_with_gene_ranges(
            scaffolds=contract_scaffolds,
            gene_range_table=Path(gene_range_table),
        )
    contract_genomes, contract_genome_scaffolds = _build_genome_specs(contract_scaffolds)
    observed_genome_names = _collect_profile_genome_names(profile_paths=profile_paths, genome=genome_scope)
    if not observed_genome_names:
        detail = genome_scope if genome_scope is not None else "all"
        raise ValueError(f"No profile rows found for genome scope: {detail}")
    contract_genome_names = {spec.genome for spec in contract_genomes}
    unknown_observed = observed_genome_names - contract_genome_names
    if unknown_observed:
        raise ValueError(
            "Profiles contain genomes that are missing from the provided BED/STB contract: "
            + ", ".join(sorted(unknown_observed))
        )
    genomes, genome_scaffolds = _subset_contract_to_genomes(
        contract_genomes=contract_genomes,
        contract_genome_scaffolds=contract_genome_scaffolds,
        genome_names=observed_genome_names,
    )
    gene_ranges = (
        _collect_gene_range_specs(
            gene_range_table=Path(gene_range_table),
            genome_scaffolds=genome_scaffolds,
        )
        if gene_range_table is not None
        else []
    )
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)

    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    total_work = len(profile_paths) * len(genomes)
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
                "stored_rows": stored_rows,
            }
        )

    h5py_module = _import_h5py()
    tmp_output = output_file.with_suffix(output_file.suffix + ".tmp")
    if tmp_output.exists():
        tmp_output.unlink()
    build_succeeded = False
    try:
        sample_rows = [(idx, path.stem, str(path.resolve()), path) for idx, path in enumerate(profile_paths)]
        metadata_rows = {
            "profiles_dir": str(profile_dir.resolve()),
            "profile_format": "classic_zipstrain_profile_parquet",
            "genome_scope": genome_scope or "all",
            "count_dtype": count_dtype,
            "layout": CURRENT_MATRIX_HDF5_SPARSE_LAYOUT if sparse else CURRENT_MATRIX_HDF5_LAYOUT,
            "matrix_value_semantics": FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS,
            "coverage_filter_min_cov": str(MATRIX_BUILD_MIN_COV),
            "memory_limit_gb": str(memory_limit_gb),
            "export_batch_mb": str(export_batch_mb),
            "separator_rows_between_scaffolds": "1",
            "input_format": "hdf5",
            "has_gene_ranges": "1" if gene_ranges else "0",
        }
        if gene_range_table is not None:
            metadata_rows["gene_range_table"] = str(Path(gene_range_table).resolve())
        if bed_file is not None:
            metadata_rows["bed_file"] = str(Path(bed_file).resolve())
        if stb_file is not None:
            metadata_rows["stb_file"] = str(Path(stb_file).resolve())
        scaffolds_by_genome_idx: dict[int, list[GenomeScaffoldOffset]] = {}
        for offset in genome_scaffolds:
            scaffolds_by_genome_idx.setdefault(offset.genome_idx, []).append(offset)

        with h5py_module.File(str(tmp_output), "w") as h5_file:
            h5_file.attrs["zipstrain_hdf5_version"] = MATRIX_HDF5_FILE_VERSION
            metadata_group = h5_file.create_group("metadata")
            for key, value in metadata_rows.items():
                metadata_group.attrs[str(key)] = str(value)

            sample_count = len(sample_rows)
            sample_chunk_len = _matrix_hdf5_sample_axis_chunk_length(sample_count)
            samples_group = h5_file.create_group("samples")
            samples_group.create_dataset(
                "sample_idx",
                data=np.asarray([sample_idx for sample_idx, _sample_name, _profile_path_str, _profile_path in sample_rows], dtype=np.int64),
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )
            _write_hdf5_string_dataset(
                samples_group,
                "sample_name",
                [sample_name for _sample_idx, sample_name, _profile_path_str, _profile_path in sample_rows],
                h5py_module=h5py_module,
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )

            _write_hdf5_genomes_group(h5_file, "genomes", genomes, h5py_module=h5py_module)
            _write_hdf5_genome_scaffolds_group(
                h5_file,
                "genome_scaffolds",
                genome_scaffolds,
                h5py_module=h5py_module,
            )
            if contract_genomes is not None and contract_genome_scaffolds is not None:
                _write_hdf5_genomes_group(
                    h5_file,
                    MATRIX_HDF5_CONTRACT_GENOMES_GROUP,
                    contract_genomes,
                    h5py_module=h5py_module,
                )
                _write_hdf5_genome_scaffolds_group(
                    h5_file,
                    MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP,
                    contract_genome_scaffolds,
                    h5py_module=h5py_module,
                )

            if gene_ranges:
                _write_hdf5_gene_ranges_group(
                    h5_file,
                    "genes",
                    gene_ranges,
                    h5py_module=h5py_module,
                )

            matrices_group = h5_file.create_group("matrices")
            matrix_datasets: dict[int, object] = {}
            sparse_states: dict[int, tuple[object, object, int]] = {}
            for spec in genomes:
                estimated_python_peak = _estimate_builder_python_peak_bytes(spec.matrix_length, count_dtype)
                if estimated_python_peak > memory_limit_bytes:
                    raise MemoryError(
                        f"Genome {spec.genome} is estimated to require about "
                        f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory, "
                        f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                    )
                if sparse:
                    matrix_group = _create_sparse_hdf5_genome_store(
                        matrices_group,
                        genome_idx=spec.genome_idx,
                        sample_count=len(sample_rows),
                        matrix_length=spec.matrix_length,
                    )
                    matrix_datasets[spec.genome_idx] = matrix_group
                    sparse_states[spec.genome_idx] = (matrix_group["indptr"], matrix_group["indices"], 0)
                else:
                    chunk_samples = _matrix_hdf5_chunk_sample_count(
                        sample_count=len(sample_rows),
                        matrix_length=spec.matrix_length,
                        dtype_name=count_dtype,
                        target_batch_mb=export_batch_mb,
                    )
                    matrix_datasets[spec.genome_idx] = matrices_group.create_dataset(
                        str(spec.genome_idx),
                        shape=(len(sample_rows), spec.matrix_length, 4),
                        dtype=COUNT_DTYPES[count_dtype],
                        chunks=(chunk_samples, spec.matrix_length, 4),
                        maxshape=(None, spec.matrix_length, 4),
                        fillvalue=0,
                    )

            for sample_idx, sample_name, _profile_path_str, profile_path in sample_rows:
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
                    matrix = _load_profile_genome_matrix(
                        profile_path=profile_path,
                        genome_spec=genome_spec,
                        genome_offsets=scaffolds_by_genome_idx[genome_spec.genome_idx],
                        count_dtype=count_dtype,
                        min_cov=MATRIX_BUILD_MIN_COV,
                    )
                    if sparse:
                        indptr_ds, indices_ds, current_nnz = sparse_states[genome_spec.genome_idx]
                        current_nnz = _append_sparse_hdf5_matrix_row(
                            indptr_dataset=indptr_ds,
                            indices_dataset=indices_ds,
                            sample_row=sample_idx,
                            flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                            current_nnz=current_nnz,
                        )
                        sparse_states[genome_spec.genome_idx] = (indptr_ds, indices_ds, current_nnz)
                    else:
                        matrix_datasets[genome_spec.genome_idx][sample_idx, :, :] = matrix
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
        tmp_output.replace(output_file)
        build_succeeded = True
    finally:
        if not build_succeeded and tmp_output.exists():
            tmp_output.unlink(missing_ok=True)
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

    return MatrixHdf5Summary(
        output_file=output_file,
        profile_files=len(profile_paths),
        sample_count=len(profile_paths),
        scaffold_count=len(genome_scaffolds),
        stored_rows=stored_rows,
        count_dtype=count_dtype,
        genome_scope=genome_scope or "all",
        memory_limit_gb=memory_limit_gb,
    )


def _append_matrix_hdf5_in_place(
    *,
    matrix_hdf5_file: Path,
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
    existing_samples: list[tuple[int, str]],
    appended_samples: list[tuple[int, str, str, Path]],
    count_dtype: str,
    layout: str,
    genome_scope: str,
    ignored_genome_count: int,
    memory_limit_gb: float,
    memory_limit_bytes: int,
    progress_callback: Optional[BuildProgressCallback],
) -> MatrixHdf5AppendSummary:
    total_work = len(appended_samples) * len(genomes)
    completed_work = 0
    stored_rows = 0
    total_sample_count = len(existing_samples) + len(appended_samples)
    existing_sample_count = len(existing_samples)

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

    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r+") as h5_file:
        sample_idx_ds = h5_file["samples"]["sample_idx"]
        sample_name_ds = h5_file["samples"]["sample_name"]
        sample_idx_ds.resize((total_sample_count,))
        sample_name_ds.resize((total_sample_count,))
        sample_idx_ds[existing_sample_count:total_sample_count] = np.asarray(
            [sample_idx for sample_idx, _sample_name, _profile_path_str, _profile_path in appended_samples],
            dtype=np.int64,
        )
        sample_name_ds[existing_sample_count:total_sample_count] = np.asarray(
            [sample_name for _sample_idx, sample_name, _profile_path_str, _profile_path in appended_samples],
            dtype=object,
        )

        sparse_states: dict[int, tuple[object, object, int]] = {}
        for spec in genomes:
            estimated_python_peak = _estimate_builder_python_peak_bytes(spec.matrix_length, count_dtype)
            if estimated_python_peak > memory_limit_bytes:
                raise MemoryError(
                    f"Genome {spec.genome} is estimated to require about "
                    f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory, "
                    f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                )
            node = h5_file[_hdf5_matrix_dataset_path(spec.genome_idx)]
            if _matrix_hdf5_layout_is_sparse(layout):
                indptr_ds = node["indptr"]
                indices_ds = node["indices"]
                indptr_ds.resize((total_sample_count + 1,))
                sparse_states[spec.genome_idx] = (
                    indptr_ds,
                    indices_ds,
                    int(indptr_ds[existing_sample_count]),
                )
            else:
                dataset = node
                dataset.resize((total_sample_count, spec.matrix_length, 4))
            for appended_offset, (_sample_idx, sample_name, _profile_path_str, profile_path) in enumerate(appended_samples):
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "processing",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": sample_name,
                            "genome": spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )
                matrix = _load_profile_genome_matrix(
                    profile_path=profile_path,
                    genome_spec=spec,
                    genome_offsets=scaffolds_by_genome_idx[spec.genome_idx],
                    count_dtype=count_dtype,
                    min_cov=MATRIX_BUILD_MIN_COV,
                )
                if _matrix_hdf5_layout_is_sparse(layout):
                    indptr_ds, indices_ds, current_nnz = sparse_states[spec.genome_idx]
                    current_nnz = _append_sparse_hdf5_matrix_row(
                        indptr_dataset=indptr_ds,
                        indices_dataset=indices_ds,
                        sample_row=existing_sample_count + appended_offset,
                        flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                        current_nnz=current_nnz,
                    )
                    sparse_states[spec.genome_idx] = (indptr_ds, indices_ds, current_nnz)
                else:
                    dataset[existing_sample_count + appended_offset, :, :] = matrix
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
                            "scaffold": "",
                            "stored_rows": stored_rows,
                        }
                    )

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

    return MatrixHdf5AppendSummary(
        output_file=matrix_hdf5_file,
        appended_profile_files=len(appended_samples),
        appended_sample_count=len(appended_samples),
        total_sample_count=total_sample_count,
        scaffold_count=len(genome_scaffolds),
        stored_rows=stored_rows,
        ignored_genome_count=ignored_genome_count,
        count_dtype=count_dtype,
        genome_scope=genome_scope or "all",
        memory_limit_gb=memory_limit_gb,
    )


def _append_matrix_hdf5_via_rewrite(
    *,
    matrix_hdf5_file: Path,
    metadata: dict[str, str],
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
    contract_genomes: Optional[list[GenomeSpec]],
    contract_genome_scaffolds: Optional[list[GenomeScaffoldOffset]],
    gene_ranges: Optional[list[GeneRangeSpec]],
    existing_samples: list[tuple[int, str]],
    appended_samples: list[tuple[int, str, str, Path]],
    count_dtype: str,
    layout: str,
    genome_scope: str,
    ignored_genome_count: int,
    memory_limit_gb: float,
    memory_limit_bytes: int,
    export_batch_mb: float,
    progress_callback: Optional[BuildProgressCallback],
) -> MatrixHdf5AppendSummary:
    total_work = len(appended_samples) * len(genomes)
    completed_work = 0
    stored_rows = 0
    source_gene_ranges = None if gene_ranges is not None else _load_matrix_hdf5_gene_ranges(matrix_hdf5_file)

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

    h5py_module = _import_h5py()
    tmp_output = matrix_hdf5_file.with_suffix(matrix_hdf5_file.suffix + ".tmp")
    if tmp_output.exists():
        tmp_output.unlink()
    append_succeeded = False
    try:
        with h5py_module.File(str(matrix_hdf5_file), "r") as src_h5, h5py_module.File(str(tmp_output), "w") as dst_h5:
            dst_h5.attrs["zipstrain_hdf5_version"] = str(
                src_h5.attrs.get("zipstrain_hdf5_version", MATRIX_HDF5_FILE_VERSION)
            )

            metadata_group = dst_h5.create_group("metadata")
            for key, value in metadata.items():
                metadata_group.attrs[str(key)] = str(value)

            all_sample_rows = list(existing_samples) + [
                (sample_idx, sample_name) for sample_idx, sample_name, _profile_path_str, _profile_path in appended_samples
            ]
            total_sample_count = len(all_sample_rows)
            existing_sample_count = len(existing_samples)
            sample_chunk_len = _matrix_hdf5_sample_axis_chunk_length(total_sample_count)
            samples_group = dst_h5.create_group("samples")
            samples_group.create_dataset(
                "sample_idx",
                data=np.asarray([sample_idx for sample_idx, _sample_name in all_sample_rows], dtype=np.int64),
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )
            _write_hdf5_string_dataset(
                samples_group,
                "sample_name",
                [sample_name for _sample_idx, sample_name in all_sample_rows],
                h5py_module=h5py_module,
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )

            _write_hdf5_genomes_group(dst_h5, "genomes", genomes, h5py_module=h5py_module)
            _write_hdf5_genome_scaffolds_group(
                dst_h5,
                "genome_scaffolds",
                genome_scaffolds,
                h5py_module=h5py_module,
            )
            if contract_genomes is not None and contract_genome_scaffolds is not None:
                _write_hdf5_genomes_group(
                    dst_h5,
                    MATRIX_HDF5_CONTRACT_GENOMES_GROUP,
                    contract_genomes,
                    h5py_module=h5py_module,
                )
                _write_hdf5_genome_scaffolds_group(
                    dst_h5,
                    MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP,
                    contract_genome_scaffolds,
                    h5py_module=h5py_module,
                )

            if gene_ranges is not None:
                if gene_ranges:
                    _write_hdf5_gene_ranges_group(
                        dst_h5,
                        "genes",
                        gene_ranges,
                        h5py_module=h5py_module,
                    )
            elif source_gene_ranges:
                _write_hdf5_gene_ranges_group(
                    dst_h5,
                    "genes",
                    source_gene_ranges,
                    h5py_module=h5py_module,
                )

            matrices_group = dst_h5.create_group("matrices")
            for spec in genomes:
                estimated_python_peak = _estimate_builder_python_peak_bytes(spec.matrix_length, count_dtype)
                if estimated_python_peak > memory_limit_bytes:
                    raise MemoryError(
                        f"Genome {spec.genome} is estimated to require about "
                        f"{_format_memory_bytes(estimated_python_peak)} of Python-side working memory, "
                        f"which exceeds the configured total limit of {_format_memory_bytes(memory_limit_bytes)}."
                    )
                src_node = src_h5[_hdf5_matrix_dataset_path(spec.genome_idx)] if _hdf5_matrix_dataset_path(spec.genome_idx) in src_h5 else None
                if _matrix_hdf5_layout_is_sparse(layout):
                    dst_group = _create_sparse_hdf5_genome_store(
                        matrices_group,
                        genome_idx=spec.genome_idx,
                        sample_count=total_sample_count,
                        matrix_length=spec.matrix_length,
                    )
                    dst_indptr = dst_group["indptr"]
                    dst_indices = dst_group["indices"]
                    if src_node is not None:
                        src_indptr = src_node["indptr"]
                        src_indices = src_node["indices"]
                        existing_indices_len = int(src_indices.shape[0])
                        if existing_indices_len > 0:
                            dst_indices.resize((existing_indices_len,))
                            dst_indices[:] = np.asarray(src_indices[...], dtype=np.int64)
                        dst_indptr[: existing_sample_count + 1] = np.asarray(
                            src_indptr[: existing_sample_count + 1],
                            dtype=np.int64,
                        )
                        current_nnz = int(dst_indptr[existing_sample_count])
                    else:
                        current_nnz = 0
                else:
                    src_chunk_samples = (
                        int(src_node.chunks[0])
                        if src_node is not None and src_node.chunks is not None and len(src_node.chunks) > 0
                        else _matrix_hdf5_chunk_sample_count(
                            sample_count=total_sample_count,
                            matrix_length=spec.matrix_length,
                            dtype_name=count_dtype,
                            target_batch_mb=export_batch_mb,
                        )
                    )
                    chunk_samples = max(1, min(total_sample_count, src_chunk_samples))
                    dst_dataset = matrices_group.create_dataset(
                        str(spec.genome_idx),
                        shape=(total_sample_count, spec.matrix_length, 4),
                        dtype=COUNT_DTYPES[count_dtype],
                        chunks=(chunk_samples, spec.matrix_length, 4),
                        maxshape=(None, spec.matrix_length, 4),
                        fillvalue=0,
                    )
                    if src_node is not None:
                        for batch_start in range(0, existing_sample_count, chunk_samples):
                            batch_stop = min(existing_sample_count, batch_start + chunk_samples)
                            dst_dataset[batch_start:batch_stop, :, :] = src_node[batch_start:batch_stop, :, :]

                for appended_offset, (_sample_idx, sample_name, _profile_path_str, profile_path) in enumerate(appended_samples):
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "processing",
                                "completed": completed_work,
                                "total": total_work,
                                "sample_name": sample_name,
                                "genome": spec.genome,
                                "scaffold": "",
                                "stored_rows": stored_rows,
                            }
                        )
                    matrix = _load_profile_genome_matrix(
                        profile_path=profile_path,
                        genome_spec=spec,
                        genome_offsets=scaffolds_by_genome_idx[spec.genome_idx],
                        count_dtype=count_dtype,
                        min_cov=MATRIX_BUILD_MIN_COV,
                    )
                    if _matrix_hdf5_layout_is_sparse(layout):
                        current_nnz = _append_sparse_hdf5_matrix_row(
                            indptr_dataset=dst_indptr,
                            indices_dataset=dst_indices,
                            sample_row=existing_sample_count + appended_offset,
                            flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                            current_nnz=current_nnz,
                        )
                    else:
                        dst_dataset[existing_sample_count + appended_offset, :, :] = matrix
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
                                "scaffold": "",
                                "stored_rows": stored_rows,
                            }
                        )
        tmp_output.replace(matrix_hdf5_file)
        append_succeeded = True
    finally:
        if not append_succeeded and tmp_output.exists():
            tmp_output.unlink(missing_ok=True)

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

    return MatrixHdf5AppendSummary(
        output_file=matrix_hdf5_file,
        appended_profile_files=len(appended_samples),
        appended_sample_count=len(appended_samples),
        total_sample_count=len(existing_samples) + len(appended_samples),
        scaffold_count=len(genome_scaffolds),
        stored_rows=stored_rows,
        ignored_genome_count=ignored_genome_count,
        count_dtype=count_dtype,
        genome_scope=genome_scope or "all",
        memory_limit_gb=memory_limit_gb,
    )


def append_matrix_hdf5(
    profile_dir: Path,
    matrix_hdf5_file: Path,
    memory_limit_gb: float = 16.0,
    progress_callback: Optional[BuildProgressCallback] = None,
    export_batch_mb: float = MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
) -> MatrixHdf5AppendSummary:
    """Append new profile parquets into an existing HDF5 matrix store.

    Genomes outside the persisted BED/STB contract are ignored. Known genomes
    must still match the stored scaffold and coordinate contract.
    """
    if export_batch_mb <= 0:
        raise ValueError("export_batch_mb must be > 0")
    profile_paths = discover_profile_parquets(profile_dir)
    matrix_hdf5_file = Path(matrix_hdf5_file).resolve()
    if not matrix_hdf5_file.exists():
        raise FileNotFoundError(f"Matrix store file does not exist: {matrix_hdf5_file}")
    if not _is_hdf5_matrix_store(matrix_hdf5_file):
        raise ValueError(
            "append-matrix-db expects the current HDF5-backed matrix store format. "
            "If you have a legacy DuckDB matrix database, convert it first with "
            "'zipstrain utilities matrix-db-to-hdf5'."
        )

    metadata = _load_matrix_hdf5_metadata(matrix_hdf5_file)
    count_dtype, genome_scope, layout = _validate_matrix_hdf5_appendable(metadata)
    genomes = _load_matrix_hdf5_genomes(matrix_hdf5_file)
    genome_scaffolds = _load_matrix_hdf5_genome_scaffolds(matrix_hdf5_file)
    contract_genomes, contract_genome_scaffolds = _load_matrix_hdf5_reference_contract(matrix_hdf5_file)
    if contract_genomes is None or contract_genome_scaffolds is None:
        raise ValueError(
            "Matrix store is missing the persisted BED/STB reference contract required for append. "
            "Rebuild the matrix store with both --bed-file and --stb-file."
        )
    contract_genome_names = {spec.genome for spec in contract_genomes}
    allowed_genome_names = (
        contract_genome_names
        if genome_scope == "all"
        else {genome_scope} & contract_genome_names
    )
    existing_samples = _load_matrix_hdf5_samples(matrix_hdf5_file)

    existing_sample_names = {sample_name for _sample_idx, sample_name in existing_samples}
    duplicate_sample_names = sorted(path.stem for path in profile_paths if path.stem in existing_sample_names)
    if duplicate_sample_names:
        raise ValueError(
            "Cannot append profiles whose sample names already exist in the matrix store: "
            + ", ".join(duplicate_sample_names)
        )

    incoming_genome_names = _collect_profile_genome_names(profile_paths=profile_paths, genome=None)
    ignored_genome_names = incoming_genome_names - contract_genome_names
    for profile_path in profile_paths:
        _validate_profile_against_matrix_contract(
            profile_path=profile_path,
            genomes=contract_genomes,
            genome_scaffolds=contract_genome_scaffolds,
            genome_scope=None if genome_scope == "all" else genome_scope,
            accepted_genomes=allowed_genome_names,
        )

    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    next_sample_idx = max((sample_idx for sample_idx, _sample_name in existing_samples), default=-1) + 1
    appended_samples = [
        (next_sample_idx + idx, profile_path.stem, str(profile_path.resolve()), profile_path)
        for idx, profile_path in enumerate(profile_paths)
    ]
    materialized_genome_names = {spec.genome for spec in genomes}
    relevant_incoming_genome_names = incoming_genome_names & allowed_genome_names
    new_genome_names = sorted(relevant_incoming_genome_names - materialized_genome_names)

    rewritten_genomes = genomes
    rewritten_genome_scaffolds = genome_scaffolds
    rewritten_gene_ranges: Optional[list[GeneRangeSpec]] = None
    needs_rewrite = not _matrix_hdf5_store_is_append_resizable(matrix_hdf5_file)

    if new_genome_names:
        new_genome_specs, new_genome_scaffolds = _subset_contract_to_genomes(
            contract_genomes=contract_genomes,
            contract_genome_scaffolds=contract_genome_scaffolds,
            genome_names=set(new_genome_names),
        )
        rewritten_genomes = sorted(genomes + new_genome_specs, key=lambda spec: spec.genome_idx)
        rewritten_genome_scaffolds = sorted(
            genome_scaffolds + new_genome_scaffolds,
            key=lambda spec: (spec.genome_idx, spec.scaffold_ordinal),
        )
        if metadata.get("has_gene_ranges") == "1":
            gene_range_path = metadata.get("gene_range_table")
            if not gene_range_path:
                raise ValueError(
                    "Matrix store contains gene ranges but is missing the source gene_range_table metadata, "
                    "so newly materialized genomes cannot be annotated safely during append."
                )
            rewritten_gene_ranges = _collect_gene_range_specs(
                gene_range_table=Path(gene_range_path),
                genome_scaffolds=rewritten_genome_scaffolds,
            )
        needs_rewrite = True

    if not needs_rewrite:
        return _append_matrix_hdf5_in_place(
            matrix_hdf5_file=matrix_hdf5_file,
            genomes=genomes,
            genome_scaffolds=genome_scaffolds,
            existing_samples=existing_samples,
            appended_samples=appended_samples,
            count_dtype=count_dtype,
            layout=layout,
            genome_scope=genome_scope,
            ignored_genome_count=len(ignored_genome_names),
            memory_limit_gb=memory_limit_gb,
            memory_limit_bytes=memory_limit_bytes,
            progress_callback=progress_callback,
        )
    return _append_matrix_hdf5_via_rewrite(
        matrix_hdf5_file=matrix_hdf5_file,
        metadata=metadata,
        genomes=rewritten_genomes,
        genome_scaffolds=rewritten_genome_scaffolds,
        contract_genomes=contract_genomes,
        contract_genome_scaffolds=contract_genome_scaffolds,
        gene_ranges=rewritten_gene_ranges,
        existing_samples=existing_samples,
        appended_samples=appended_samples,
        count_dtype=count_dtype,
        layout=layout,
        genome_scope=genome_scope,
        ignored_genome_count=len(ignored_genome_names),
        memory_limit_gb=memory_limit_gb,
        memory_limit_bytes=memory_limit_bytes,
        export_batch_mb=export_batch_mb,
        progress_callback=progress_callback,
    )


def append_matrix_db(
    profile_dir: Path,
    matrix_db_file: Path,
    memory_limit_gb: float = 16.0,
    progress_callback: Optional[BuildProgressCallback] = None,
) -> MatrixDbAppendSummary:
    summary = append_matrix_hdf5(
        profile_dir=profile_dir,
        matrix_hdf5_file=matrix_db_file,
        memory_limit_gb=memory_limit_gb,
        progress_callback=progress_callback,
    )
    return MatrixDbAppendSummary(
        output_file=summary.output_file,
        appended_profile_files=summary.appended_profile_files,
        appended_sample_count=summary.appended_sample_count,
        total_sample_count=summary.total_sample_count,
        scaffold_count=summary.scaffold_count,
        stored_rows=summary.stored_rows,
        ignored_genome_count=summary.ignored_genome_count,
        count_dtype=summary.count_dtype,
        genome_scope=summary.genome_scope,
        memory_limit_gb=summary.memory_limit_gb,
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


def _load_matrix_hdf5_metadata(matrix_hdf5_file: Path) -> dict[str, str]:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        if "metadata" not in h5_file:
            raise ValueError(f"HDF5 matrix store is missing /metadata: {matrix_hdf5_file}")
        return {str(key): str(value) for key, value in h5_file["metadata"].attrs.items()}


def _load_matrix_hdf5_samples(matrix_hdf5_file: Path) -> list[tuple[int, str]]:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        sample_group = h5_file["samples"]
        sample_ids = np.asarray(sample_group["sample_idx"][...], dtype=np.int64)
        sample_names = _read_hdf5_string_dataset(sample_group["sample_name"])
    return [(int(sample_idx), str(sample_name)) for sample_idx, sample_name in zip(sample_ids.tolist(), sample_names)]


def _load_matrix_hdf5_genomes_from_group(
    matrix_hdf5_file: Path,
    group_name: str,
    genome: Optional[str] = None,
) -> list[GenomeSpec]:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        if group_name not in h5_file:
            return []
        genome_group = h5_file[group_name]
        genome_ids = np.asarray(genome_group["genome_idx"][...], dtype=np.int64)
        genome_names = _read_hdf5_string_dataset(genome_group["genome"])
        matrix_lengths = np.asarray(genome_group["matrix_length"][...], dtype=np.int64)
        true_lengths = np.asarray(genome_group["true_length"][...], dtype=np.int64)
        scaffold_counts = np.asarray(genome_group["scaffold_count"][...], dtype=np.int64)
    specs = [
        GenomeSpec(
            genome_idx=int(genome_idx),
            genome=str(genome_name),
            matrix_length=int(matrix_length),
            true_length=int(true_length),
            scaffold_count=int(scaffold_count),
        )
        for genome_idx, genome_name, matrix_length, true_length, scaffold_count in zip(
            genome_ids.tolist(),
            genome_names,
            matrix_lengths.tolist(),
            true_lengths.tolist(),
            scaffold_counts.tolist(),
        )
    ]
    if genome is None:
        return specs
    return [spec for spec in specs if spec.genome == genome]


def _load_matrix_hdf5_genomes(
    matrix_hdf5_file: Path,
    genome: Optional[str] = None,
) -> list[GenomeSpec]:
    return _load_matrix_hdf5_genomes_from_group(matrix_hdf5_file, "genomes", genome=genome)


def _load_matrix_hdf5_genome_scaffolds_from_group(
    matrix_hdf5_file: Path,
    group_name: str,
    genome_idx: Optional[int] = None,
    genome: Optional[str] = None,
) -> list[GenomeScaffoldOffset]:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        if group_name not in h5_file:
            return []
        scaffold_group = h5_file[group_name]
        rows = list(
            zip(
                np.asarray(scaffold_group["genome_idx"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["scaffold_ordinal"][...], dtype=np.int64).tolist(),
                _read_hdf5_string_dataset(scaffold_group["genome"]),
                _read_hdf5_string_dataset(scaffold_group["chrom"]),
                np.asarray(scaffold_group["axis_start"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["axis_end"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["index_base"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["vector_length"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["min_pos"][...], dtype=np.int64).tolist(),
                np.asarray(scaffold_group["max_pos"][...], dtype=np.int64).tolist(),
            )
        )
    offsets = [
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
        if (genome_idx is None or int(genome_idx_val) == int(genome_idx))
        and (genome is None or str(genome_name) == genome)
    ]
    if genome_idx is None and genome is None:
        return offsets
    return sorted(offsets, key=lambda item: (item.genome_idx, item.scaffold_ordinal))


def _load_matrix_hdf5_genome_scaffolds(
    matrix_hdf5_file: Path,
    genome_idx: Optional[int] = None,
    genome: Optional[str] = None,
) -> list[GenomeScaffoldOffset]:
    return _load_matrix_hdf5_genome_scaffolds_from_group(
        matrix_hdf5_file,
        "genome_scaffolds",
        genome_idx=genome_idx,
        genome=genome,
    )


def _load_matrix_hdf5_reference_contract(
    matrix_hdf5_file: Path,
) -> tuple[Optional[list[GenomeSpec]], Optional[list[GenomeScaffoldOffset]]]:
    contract_genomes = _load_matrix_hdf5_genomes_from_group(
        matrix_hdf5_file,
        MATRIX_HDF5_CONTRACT_GENOMES_GROUP,
    )
    contract_genome_scaffolds = _load_matrix_hdf5_genome_scaffolds_from_group(
        matrix_hdf5_file,
        MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP,
    )
    if not contract_genomes or not contract_genome_scaffolds:
        return None, None
    return contract_genomes, contract_genome_scaffolds


def _load_matrix_hdf5_gene_ranges(
    matrix_hdf5_file: Path,
    genome_idx: Optional[int] = None,
    genome: Optional[str] = None,
) -> list[GeneRangeSpec]:
    h5py_module = _import_h5py()
    with h5py_module.File(str(matrix_hdf5_file), "r") as h5_file:
        if "genes" not in h5_file:
            return []
        gene_group = h5_file["genes"]
        rows = list(
            zip(
                np.asarray(gene_group["gene_idx"][...], dtype=np.int64).tolist(),
                np.asarray(gene_group["genome_idx"][...], dtype=np.int64).tolist(),
                _read_hdf5_string_dataset(gene_group["genome"]),
                _read_hdf5_string_dataset(gene_group["chrom"]),
                _read_hdf5_string_dataset(gene_group["gene"]),
                np.asarray(gene_group["axis_start"][...], dtype=np.int64).tolist(),
                np.asarray(gene_group["axis_end"][...], dtype=np.int64).tolist(),
            )
        )
    specs = [
        GeneRangeSpec(
            gene_idx=int(gene_idx_val),
            gene=str(gene_name),
            genome_idx=int(genome_idx_val),
            genome=str(genome_name),
            chrom=str(chrom),
            axis_start=int(axis_start),
            axis_end=int(axis_end),
        )
        for gene_idx_val, genome_idx_val, genome_name, chrom, gene_name, axis_start, axis_end in rows
        if (genome_idx is None or int(genome_idx_val) == int(genome_idx))
        and (genome is None or str(genome_name) == genome)
    ]
    return sorted(specs, key=lambda item: (item.genome_idx, item.axis_start, item.axis_end, item.gene))


def _group_gene_ranges_by_genome(
    gene_ranges: list[GeneRangeSpec],
) -> dict[int, list[GeneRangeSpec]]:
    grouped: dict[int, list[GeneRangeSpec]] = {}
    for spec in gene_ranges:
        grouped.setdefault(spec.genome_idx, []).append(spec)
    return grouped


def export_matrix_db_hdf5(
    matrix_db_file: Path,
    output_file: Path,
    progress_callback: Optional[BuildProgressCallback] = None,
    export_batch_mb: float = MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
    sparse: bool = False,
) -> Path:
    """Convert a legacy DuckDB matrix database into the current HDF5 store."""
    matrix_db_file = matrix_db_file.resolve()
    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if export_batch_mb <= 0:
        raise ValueError("export_batch_mb must be > 0")

    h5py_module = _import_h5py()
    read_conn = duckdb.connect(str(matrix_db_file), read_only=True)
    try:
        read_conn.execute("SET preserve_insertion_order=false")
        metadata = _load_matrix_db_metadata(read_conn)
        samples = _load_matrix_db_samples(read_conn)
        genomes = _load_matrix_db_genomes(read_conn)
        genome_scaffolds = _load_matrix_db_genome_scaffolds(read_conn)
        dtype_name = str(metadata.get("count_dtype", ""))
        if dtype_name not in COUNT_DTYPES:
            raise ValueError(
                f"Matrix DB count dtype '{dtype_name}' is not supported for HDF5 export."
            )
        np_dtype = COUNT_DTYPES[dtype_name]
        row_pos_by_sample_idx = {sample_idx: pos for pos, (sample_idx, _sample_name) in enumerate(samples)}
        converted_metadata = {str(key): str(value) for key, value in metadata.items()}
        converted_metadata["source_layout"] = str(metadata.get("layout", ""))
        converted_metadata["layout"] = CURRENT_MATRIX_HDF5_SPARSE_LAYOUT if sparse else CURRENT_MATRIX_HDF5_LAYOUT
        converted_metadata["input_format"] = "hdf5"
        converted_metadata["has_gene_ranges"] = "0"
        tmp_output = output_file.with_suffix(output_file.suffix + ".tmp")
        if tmp_output.exists():
            tmp_output.unlink()
        sample_count = len(samples)
        total_work = sample_count * len(genomes)
        completed_work = 0
        stored_rows = 0
        genome_scope = metadata.get("genome_scope", "all")
        if progress_callback is not None:
            progress_callback(
                {
                    "phase": "start",
                    "completed": completed_work,
                    "total": total_work,
                    "sample_name": "",
                    "genome": genome_scope,
                    "scaffold": "",
                    "stored_rows": stored_rows,
                }
            )
        with h5py_module.File(str(tmp_output), "w") as h5_file:
            h5_file.attrs["zipstrain_hdf5_version"] = MATRIX_HDF5_FILE_VERSION
            metadata_group = h5_file.create_group("metadata")
            for key, value in converted_metadata.items():
                metadata_group.attrs[str(key)] = str(value)

            sample_chunk_len = _matrix_hdf5_sample_axis_chunk_length(sample_count)
            samples_group = h5_file.create_group("samples")
            samples_group.create_dataset(
                "sample_idx",
                data=np.asarray([sample_idx for sample_idx, _sample_name in samples], dtype=np.int64),
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )
            _write_hdf5_string_dataset(
                samples_group,
                "sample_name",
                [sample_name for _sample_idx, sample_name in samples],
                h5py_module=h5py_module,
                chunks=(sample_chunk_len,),
                maxshape=(None,),
            )

            genomes_group = h5_file.create_group("genomes")
            genomes_group.create_dataset(
                "genome_idx",
                data=np.asarray([spec.genome_idx for spec in genomes], dtype=np.int64),
            )
            _write_hdf5_string_dataset(
                genomes_group,
                "genome",
                [spec.genome for spec in genomes],
                h5py_module=h5py_module,
            )
            genomes_group.create_dataset(
                "matrix_length",
                data=np.asarray([spec.matrix_length for spec in genomes], dtype=np.int64),
            )
            genomes_group.create_dataset(
                "true_length",
                data=np.asarray([spec.true_length for spec in genomes], dtype=np.int64),
            )
            genomes_group.create_dataset(
                "scaffold_count",
                data=np.asarray([spec.scaffold_count for spec in genomes], dtype=np.int64),
            )

            scaffold_group = h5_file.create_group("genome_scaffolds")
            scaffold_group.create_dataset(
                "genome_idx",
                data=np.asarray([spec.genome_idx for spec in genome_scaffolds], dtype=np.int64),
            )
            scaffold_group.create_dataset(
                "scaffold_ordinal",
                data=np.asarray([spec.scaffold_ordinal for spec in genome_scaffolds], dtype=np.int64),
            )
            _write_hdf5_string_dataset(
                scaffold_group,
                "genome",
                [spec.genome for spec in genome_scaffolds],
                h5py_module=h5py_module,
            )
            _write_hdf5_string_dataset(
                scaffold_group,
                "chrom",
                [spec.chrom for spec in genome_scaffolds],
                h5py_module=h5py_module,
            )
            for field_name in ("axis_start", "axis_end", "index_base", "vector_length", "min_pos", "max_pos"):
                scaffold_group.create_dataset(
                    field_name,
                    data=np.asarray([getattr(spec, field_name) for spec in genome_scaffolds], dtype=np.int64),
                )

            matrices_group = h5_file.create_group("matrices")
            for spec in genomes:
                if progress_callback is not None:
                    progress_callback(
                        {
                            "phase": "processing",
                            "completed": completed_work,
                            "total": total_work,
                            "sample_name": "",
                            "genome": spec.genome,
                            "scaffold": "",
                            "stored_rows": stored_rows,
                            "detail": f"export_batch_mb={export_batch_mb:g}",
                        }
                    )
                chunk_samples = _matrix_hdf5_chunk_sample_count(
                    sample_count=sample_count,
                    matrix_length=spec.matrix_length,
                    dtype_name=dtype_name,
                    target_batch_mb=export_batch_mb,
                )
                if sparse:
                    matrix_group = _create_sparse_hdf5_genome_store(
                        matrices_group,
                        genome_idx=spec.genome_idx,
                        sample_count=sample_count,
                        matrix_length=spec.matrix_length,
                    )
                    indptr_ds = matrix_group["indptr"]
                    indices_ds = matrix_group["indices"]
                    current_nnz = 0
                else:
                    matrix_dataset = matrices_group.create_dataset(
                        str(spec.genome_idx),
                        shape=(sample_count, spec.matrix_length, 4),
                        dtype=np_dtype,
                        chunks=(chunk_samples, spec.matrix_length, 4),
                        maxshape=(None, spec.matrix_length, 4),
                        fillvalue=0,
                    )
                ordered_sample_ids = [sample_idx for sample_idx, _sample_name in samples]
                for batch_start in range(0, sample_count, chunk_samples):
                    batch_sample_ids = ordered_sample_ids[batch_start:batch_start + chunk_samples]
                    if not batch_sample_ids:
                        continue
                    rows = read_conn.execute(
                        f"""
                        SELECT sample_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob
                        FROM matrix_db_sample_genome_matrices
                        WHERE genome_idx = ? AND sample_idx IN ({','.join(['?'] * len(batch_sample_ids))})
                        ORDER BY sample_idx
                        """,
                        [spec.genome_idx, *batch_sample_ids],
                    ).fetchall()
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "processing",
                                "completed": completed_work,
                                "total": total_work,
                                "sample_name": "",
                                "genome": spec.genome,
                                "scaffold": "",
                                "stored_rows": stored_rows,
                                "detail": f"batch_samples={len(rows)}",
                            }
                        )
                    for sample_idx, stored_dtype, matrix_rows, matrix_cols, matrix_blob in rows:
                        rows_int = int(matrix_rows)
                        cols_int = int(matrix_cols)
                        if rows_int != spec.matrix_length or cols_int != 4:
                            raise ValueError(
                                f"Stored genome matrix shape mismatch for genome_idx={spec.genome_idx}, sample_idx={sample_idx}: "
                                f"expected ({spec.matrix_length}, 4), found ({rows_int}, {cols_int})"
                            )
                        if str(stored_dtype) != dtype_name:
                            raise ValueError(
                                f"Stored genome matrix dtype mismatch for genome_idx={spec.genome_idx}, sample_idx={sample_idx}: "
                                f"expected {dtype_name}, found {stored_dtype}"
                            )
                        matrix = _unpack_matrix(
                            memoryview(matrix_blob),
                            dtype_name,
                            (rows_int, cols_int),
                        )
                        row_pos = row_pos_by_sample_idx[int(sample_idx)]
                        if sparse:
                            current_nnz = _append_sparse_hdf5_matrix_row(
                                indptr_dataset=indptr_ds,
                                indices_dataset=indices_ds,
                                sample_row=row_pos,
                                flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                                current_nnz=current_nnz,
                            )
                        else:
                            matrix_dataset[row_pos, :, :] = matrix
                        completed_work += 1
                        stored_rows += 1
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "sample_name": "",
                                    "genome": spec.genome,
                                    "scaffold": "",
                                    "stored_rows": stored_rows,
                                }
                            )
        tmp_output.replace(output_file)
    finally:
        read_conn.close()
    if progress_callback is not None:
        progress_callback(
            {
                "phase": "done",
                "completed": completed_work,
                "total": total_work,
                "sample_name": "",
                "genome": genome_scope,
                "scaffold": "",
                "stored_rows": stored_rows,
            }
        )
    return output_file


def _plan_chunk_sizes(
    vector_length: int,
    remaining_targets: int,
    dtype_name: str,
    memory_limit_bytes: int,
    backend_kind: str,
) -> tuple[int, int]:
    dtype_bytes = np.dtype(COUNT_DTYPES[dtype_name]).itemsize
    reserve = int(memory_limit_bytes * 0.15)
    budget = max(memory_limit_bytes - reserve, 64 * 1024 * 1024)

    per_position_anchor = 4 * dtype_bytes + 8
    per_position_target = 4 * dtype_bytes + 8
    if backend_kind == "torch":
        per_position_anchor += 4 * dtype_bytes
        per_position_target += 4 * dtype_bytes

    max_targets_full = max(
        1,
        int((budget / max(vector_length, 1) - per_position_anchor) // max(per_position_target, 1)),
    )
    if max_targets_full >= 1:
        return min(remaining_targets, max_targets_full), vector_length

    return 1, vector_length


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
    total_mask = (anchor_cov.unsqueeze(1) > 0) & (target_cov > 0)
    shared_mask = _shared_mask_presence_torch_tensors(
        torch_module=torch_module,
        anchor_t=anchor_t,
        targets_t=targets_t,
    )
    shared = shared_mask.sum(dim=0, dtype=torch_module.int64)
    return totals, shared, total_mask, shared_mask


def _compare_tile_presence_torch_tensors_with_shared_mask(
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


def _accumulate_gene_counts_from_full_torch_masks(
    *,
    torch_module,
    total_mask,
    shared_mask,
    gene_ranges: list[GeneRangeSpec],
) -> tuple[np.ndarray, np.ndarray]:
    if not gene_ranges:
        target_count = int(total_mask.shape[1]) if int(total_mask.ndim) > 1 else 0
        empty = np.zeros((0, target_count), dtype=np.int64)
        return empty, empty

    device = total_mask.device
    prefix_dtype = getattr(torch_module, "int32", torch_module.int64)
    gene_starts = torch_module.tensor(
        [int(gene.axis_start) for gene in gene_ranges],
        dtype=torch_module.int64,
        device=device,
    )
    gene_stops = torch_module.tensor(
        [int(gene.axis_end) for gene in gene_ranges],
        dtype=torch_module.int64,
        device=device,
    )
    total_prefix = total_mask.to(prefix_dtype).cumsum(dim=0)
    shared_prefix = shared_mask.to(prefix_dtype).cumsum(dim=0)
    total_stop = total_prefix.index_select(0, gene_stops)
    shared_stop = shared_prefix.index_select(0, gene_stops)
    start_positions = torch_module.clamp(gene_starts - 1, min=0)
    has_start = (gene_starts > 0).unsqueeze(1).to(prefix_dtype)
    total_start = total_prefix.index_select(0, start_positions) * has_start
    shared_start = shared_prefix.index_select(0, start_positions) * has_start

    gene_total_positions = (total_stop - total_start).detach().cpu().numpy().astype(np.int64, copy=False)
    gene_share_allele_pos = (shared_stop - shared_start).detach().cpu().numpy().astype(np.int64, copy=False)
    return gene_total_positions, gene_share_allele_pos


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


def _max_ibs_from_shared_mask_numpy(shared_mask) -> np.ndarray:
    shared_mask_np = _shared_mask_to_numpy(shared_mask)
    target_count = int(shared_mask_np.shape[1]) if int(shared_mask_np.ndim) > 1 else 0
    current_runs = np.zeros(target_count, dtype=np.int64)
    max_runs = np.zeros(target_count, dtype=np.int64)
    _, max_runs = _update_ibs_numpy(
        shared_mask=shared_mask_np,
        current_runs=current_runs,
        max_runs=max_runs,
    )
    return max_runs


def _max_ibs_from_shared_mask_torch(
    *,
    torch_module,
    shared_mask,
):
    """Return longest consecutive shared-position runs per target on the current torch device."""
    if int(shared_mask.ndim) != 2:
        raise ValueError("shared_mask must be a 2D positions-by-target tensor.")
    position_count = int(shared_mask.shape[0])
    target_count = int(shared_mask.shape[1])
    if target_count == 0:
        return torch_module.zeros(0, dtype=torch_module.int64, device=shared_mask.device)
    if position_count == 0:
        return torch_module.zeros(target_count, dtype=torch_module.int64, device=shared_mask.device)
    nz = shared_mask.to(torch_module.int32)
    csum = nz.cumsum(dim=0, dtype=torch_module.int32)
    zero_mask = (~shared_mask).to(torch_module.int32)
    reset = (csum * zero_mask).cummax(dim=0).values
    runs = csum - reset
    return runs.amax(dim=0).to(torch_module.int64)


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


def _make_gene_arrow_table_from_compare_payload(
    payload: dict[str, object],
) -> Optional[pa.Table]:
    gene_names = list(payload.get("gene_names") or [])
    if not gene_names:
        return None
    gene_total_all = np.asarray(payload.get("gene_total_positions"), dtype=np.int64)
    gene_shared_all = np.asarray(payload.get("gene_share_allele_pos"), dtype=np.int64)
    if gene_total_all.size == 0 or gene_shared_all.size == 0:
        return None
    valid_mask = gene_total_all > 0
    if not bool(valid_mask.any()):
        return None

    gene_idx_flat, target_idx_flat = np.nonzero(valid_mask)
    total_positions = gene_total_all[gene_idx_flat, target_idx_flat].astype(np.float64, copy=False)
    shared_positions = gene_shared_all[gene_idx_flat, target_idx_flat].astype(np.float64, copy=False)
    ani_values = (shared_positions / total_positions) * 100.0

    sample_2_idx_all = np.asarray(payload["sample_2_idx"], dtype=np.int64)
    sample_2_all = list(payload["sample_2"])
    sample_1_idx = int(payload["sample_1_idx"])
    sample_1 = str(payload["sample_1"])
    genome_idx = int(payload["genome_idx"])
    genome = str(payload["genome"])

    row_count = len(gene_idx_flat)
    return pa.Table.from_arrays(
        [
            pa.array([sample_1_idx] * row_count, type=pa.int64()),
            pa.array(sample_2_idx_all[target_idx_flat].astype(np.int64, copy=False), type=pa.int64()),
            pa.array([sample_1] * row_count, type=pa.string()),
            pa.array([sample_2_all[idx] for idx in target_idx_flat.tolist()], type=pa.string()),
            pa.array([genome_idx] * row_count, type=pa.int64()),
            pa.array([genome] * row_count, type=pa.string()),
            pa.array([gene_names[idx] for idx in gene_idx_flat.tolist()], type=pa.string()),
            pa.array(ani_values.tolist(), type=pa.float64()),
        ],
        schema=matrix_compare_gene_result_db_schema(),
    )


def _prepare_torch_matrix(
    compute_backend: MatrixPairComputeBackend,
    matrix,
    matrix_value_semantics: str,
):
    torch_module = compute_backend.torch
    if torch_module is not None and callable(getattr(torch_module, "is_tensor", None)) and torch_module.is_tensor(matrix):
        tensor = matrix.contiguous()
    else:
        tensor = torch_module.from_numpy(np.ascontiguousarray(matrix))
    if compute_backend.device == "cuda" and hasattr(tensor, "pin_memory"):
        tensor = tensor.pin_memory()
        kwargs: dict[str, object] = {"device": compute_backend.device, "non_blocking": True}
    else:
        kwargs = {"device": compute_backend.device}
    if matrix_value_semantics == FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
        kwargs["dtype"] = torch_module.float32
    return tensor.to(**kwargs)


def _torch_host_matrix_dtype(
    torch_module,
    dtype_name: str,
    matrix_value_semantics: str,
):
    if matrix_value_semantics == FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
        return torch_module.float32
    attr_name = dtype_name
    if hasattr(torch_module, attr_name):
        return getattr(torch_module, attr_name)
    return torch_module.int32


def _torch_count_matrix_dtype(
    torch_module,
    dtype_name: str,
):
    if hasattr(torch_module, dtype_name):
        return getattr(torch_module, dtype_name)
    return None


class _Hdf5GenomeMatrixTorchDataset:
    def __init__(
        self,
        matrix_hdf5_file: Path,
        genome_idx: int,
        matrix_length: int,
        dtype_name: str,
        matrix_value_semantics: str,
        torch_module=None,
    ) -> None:
        self.matrix_hdf5_file = str(matrix_hdf5_file)
        self.genome_idx = int(genome_idx)
        self.matrix_length = int(matrix_length)
        self._h5_file = None
        self._matrix_dataset = None
        self._matrix_sparse_indptr = None
        self._matrix_sparse_indices = None
        self.torch_module = torch_module or importlib.import_module("torch")
        self.h5py_module = _import_h5py()
        self.numpy_dtype = COUNT_DTYPES[dtype_name]
        self.host_dtype = _torch_host_matrix_dtype(
            torch_module=self.torch_module,
            dtype_name=dtype_name,
            matrix_value_semantics=matrix_value_semantics,
        )

    def _ensure_open(self) -> None:
        if self._h5_file is not None:
            return
        self._h5_file = self.h5py_module.File(self.matrix_hdf5_file, "r")
        node = self._h5_file[_hdf5_matrix_dataset_path(self.genome_idx)]
        if hasattr(node, "shape"):
            self._matrix_dataset = node
        else:
            self._matrix_sparse_indptr = node["indptr"]
            self._matrix_sparse_indices = node["indices"]

    def _load_batch_numpy(self, sample_rows: np.ndarray) -> np.ndarray:
        self._ensure_open()
        if self._matrix_dataset is not None:
            return np.asarray(self._matrix_dataset[sample_rows.tolist(), :, :])
        return _load_dense_rows_from_sparse_hdf5(
            indptr_dataset=self._matrix_sparse_indptr,
            indices_dataset=self._matrix_sparse_indices,
            sample_rows=sample_rows,
            matrix_length=self.matrix_length,
            numpy_dtype=self.numpy_dtype,
        )

    def close(self) -> None:
        if self._h5_file is not None:
            self._h5_file.close()
        self._h5_file = None
        self._matrix_dataset = None
        self._matrix_sparse_indptr = None
        self._matrix_sparse_indices = None

    def __del__(self) -> None:
        self.close()

    def __len__(self) -> int:
        self._ensure_open()
        if self._matrix_dataset is not None:
            return int(self._matrix_dataset.shape[0])
        return int(self._matrix_sparse_indptr.shape[0] - 1)

    def __getitem__(self, sample_row: int):
        matrix_np = self._load_batch_numpy(np.asarray([int(sample_row)], dtype=np.int64))[0]
        matrix_tensor = self.torch_module.from_numpy(matrix_np)
        if matrix_tensor.dtype != self.host_dtype:
            matrix_tensor = matrix_tensor.to(dtype=self.host_dtype)
        return matrix_tensor

    def load_range(self, start: int, stop: int):
        batch_np = self._load_batch_numpy(np.arange(int(start), int(stop), dtype=np.int64))
        batch_tensor = self.torch_module.from_numpy(batch_np)
        if batch_tensor.dtype != self.host_dtype:
            batch_tensor = batch_tensor.to(dtype=self.host_dtype)
        return batch_tensor.permute(1, 2, 0).contiguous()


class _Hdf5GenomeMatrixNumpyDataset:
    def __init__(
        self,
        matrix_hdf5_file: Path,
        genome_idx: int,
        matrix_length: int,
        dtype_name: str,
    ) -> None:
        self.matrix_hdf5_file = str(matrix_hdf5_file)
        self.genome_idx = int(genome_idx)
        self.matrix_length = int(matrix_length)
        self.h5py_module = _import_h5py()
        self.numpy_dtype = COUNT_DTYPES[dtype_name]
        self._h5_file = None
        self._matrix_dataset = None
        self._matrix_sparse_indptr = None
        self._matrix_sparse_indices = None

    def _ensure_open(self) -> None:
        if self._h5_file is not None:
            return
        self._h5_file = self.h5py_module.File(self.matrix_hdf5_file, "r")
        node = self._h5_file[_hdf5_matrix_dataset_path(self.genome_idx)]
        if hasattr(node, "shape"):
            self._matrix_dataset = node
        else:
            self._matrix_sparse_indptr = node["indptr"]
            self._matrix_sparse_indices = node["indices"]

    def _load_batch_numpy(self, sample_rows: np.ndarray) -> np.ndarray:
        self._ensure_open()
        if self._matrix_dataset is not None:
            return np.asarray(self._matrix_dataset[sample_rows.tolist(), :, :], dtype=self.numpy_dtype)
        return _load_dense_rows_from_sparse_hdf5(
            indptr_dataset=self._matrix_sparse_indptr,
            indices_dataset=self._matrix_sparse_indices,
            sample_rows=sample_rows,
            matrix_length=self.matrix_length,
            numpy_dtype=self.numpy_dtype,
        )

    def close(self) -> None:
        if self._h5_file is not None:
            self._h5_file.close()
        self._h5_file = None
        self._matrix_dataset = None
        self._matrix_sparse_indptr = None
        self._matrix_sparse_indices = None

    def __del__(self) -> None:
        self.close()

    def get_row(self, sample_row: int) -> np.ndarray:
        return self._load_batch_numpy(np.asarray([int(sample_row)], dtype=np.int64))[0]

    def load_indices(self, sample_rows: np.ndarray) -> np.ndarray:
        batch_np = self._load_batch_numpy(sample_rows)
        return np.transpose(batch_np, (1, 2, 0))


def _load_target_queue_block_for_hdf5_torch(
    matrix_hdf5_file: Path,
    genome_idx: int,
    block_rows: list[tuple[int, str]],
    block_start: int,
    matrix_length: int,
    dtype_name: str,
    matrix_value_semantics: str,
    torch_module=None,
) -> tuple[np.ndarray, list[str], np.ndarray, object]:
    if torch_module is None:
        torch_module = importlib.import_module("torch")
    block_ids = np.array([sample_idx for sample_idx, _sample_name in block_rows], dtype=np.int64)
    block_names = [sample_name for _sample_idx, sample_name in block_rows]
    zero_matrix = np.zeros((matrix_length, 4), dtype=COUNT_DTYPES[dtype_name])
    dataset = _Hdf5GenomeMatrixTorchDataset(
        matrix_hdf5_file=matrix_hdf5_file,
        genome_idx=genome_idx,
        matrix_length=matrix_length,
        dtype_name=dtype_name,
        matrix_value_semantics=matrix_value_semantics,
        torch_module=torch_module,
    )
    try:
        target_tensor = dataset.load_range(block_start, block_start + len(block_rows))
    finally:
        dataset.close()
    return block_ids, block_names, zero_matrix, target_tensor


def _load_anchor_queue_batch_for_hdf5_torch(
    matrix_hdf5_file: Path,
    genome_idx: int,
    batch_rows: list[tuple[int, str]],
    batch_start: int,
    matrix_length: int,
    dtype_name: str,
    zero_matrix: np.ndarray,
    matrix_value_semantics: str,
    torch_module=None,
) -> list[tuple[int, str, object]]:
    if not batch_rows:
        return []
    if torch_module is None:
        torch_module = importlib.import_module("torch")
    dataset = _Hdf5GenomeMatrixTorchDataset(
        matrix_hdf5_file=matrix_hdf5_file,
        genome_idx=genome_idx,
        matrix_length=matrix_length,
        dtype_name=dtype_name,
        matrix_value_semantics=matrix_value_semantics,
        torch_module=torch_module,
    )
    try:
        anchor_tensor = dataset.load_range(batch_start, batch_start + len(batch_rows))
    finally:
        dataset.close()
    return [
        (
            sample_idx,
            sample_name,
            anchor_tensor[:, :, pos],
        )
        for pos, (sample_idx, sample_name) in enumerate(batch_rows)
    ]


def _load_target_prefetch_unit_for_hdf5_torch(
    unit_index: int,
    matrix_hdf5_file: Path,
    genome_idx: int,
    block_rows: list[tuple[int, str]],
    block_start: int,
    matrix_length: int,
    dtype_name: str,
    matrix_value_semantics: str,
) -> tuple[int, tuple[np.ndarray, list[str], np.ndarray, object]]:
    return (
        unit_index,
        _load_target_queue_block_for_hdf5_torch(
            matrix_hdf5_file=matrix_hdf5_file,
            genome_idx=genome_idx,
            block_rows=block_rows,
            block_start=block_start,
            matrix_length=matrix_length,
            dtype_name=dtype_name,
            matrix_value_semantics=matrix_value_semantics,
        ),
    )


def _build_matrix_io_executor(kind: str):
    if kind == "thread":
        return ThreadPoolExecutor(max_workers=1)
    if kind == "process":
        return ProcessPoolExecutor(max_workers=1)
    raise ValueError(f"Unsupported matrix I/O executor kind '{kind}'.")


def _compare_anchor_against_target_chunk_torch_device(
    compute_backend: MatrixPairComputeBackend,
    anchor_torch,
    target_torch,
    vector_length: int,
    matrix_value_semantics: str,
    need_ibs: bool = False,
    gene_ranges: Optional[list[GeneRangeSpec]] = None,
):
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
    ibs_shared_mask = None
    gene_total_positions = None
    gene_share_allele_pos = None
    if gene_ranges:
        total_inc, shared_inc, total_mask, shared_mask = _compare_tile_presence_torch_tensors_with_mask(
            torch_module=compute_backend.torch,
            anchor_t=anchor_torch[:vector_length, :],
            targets_t=target_torch[:vector_length, :, :],
        )
        chunk_totals_torch += total_inc
        chunk_shared_torch += shared_inc
        gene_total_positions, gene_share_allele_pos = _accumulate_gene_counts_from_full_torch_masks(
            torch_module=compute_backend.torch,
            total_mask=total_mask,
            shared_mask=shared_mask,
            gene_ranges=gene_ranges,
        )
        if need_ibs:
            ibs_shared_mask = shared_mask
    elif need_ibs:
        total_inc, shared_inc, shared_mask = _compare_tile_presence_torch_tensors_with_shared_mask(
            torch_module=compute_backend.torch,
            anchor_t=anchor_torch[:vector_length, :],
            targets_t=target_torch[:vector_length, :, :],
        )
        chunk_totals_torch += total_inc
        chunk_shared_torch += shared_inc
        ibs_shared_mask = shared_mask
    else:
        total_inc, shared_inc = _compare_tile_presence_torch_tensors(
            torch_module=compute_backend.torch,
            anchor_t=anchor_torch[:vector_length, :],
            targets_t=target_torch[:vector_length, :, :],
        )
        chunk_totals_torch += total_inc
        chunk_shared_torch += shared_inc
    return chunk_totals_torch, chunk_shared_torch, ibs_shared_mask, gene_total_positions, gene_share_allele_pos


def _download_torch_result_tensor_batch(
    compute_backend: MatrixPairComputeBackend,
    totals_tensors: list,
    shared_tensors: list,
    max_run_tensors: Optional[list] = None,
) -> np.ndarray:
    if not totals_tensors:
        channel_count = 3 if max_run_tensors else 2
        return np.zeros((0, channel_count, 0), dtype=np.int64)
    if len(totals_tensors) != len(shared_tensors):
        raise ValueError("totals_tensors and shared_tensors must have the same length")
    if max_run_tensors is not None and len(totals_tensors) != len(max_run_tensors):
        raise ValueError("totals_tensors and max_run_tensors must have the same length")
    torch_module = compute_backend.torch
    lengths = [int(tensor.shape[0]) for tensor in totals_tensors]
    max_length = max(lengths)
    padded_totals = []
    padded_shared = []
    padded_max_runs = [] if max_run_tensors is not None else None
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
    if max_run_tensors is not None and padded_max_runs is not None:
        for max_run_tensor, length in zip(max_run_tensors, lengths):
            if int(max_run_tensor.shape[0]) != length:
                raise ValueError("totals and max-run tensors must have matching lengths")
            pad_width = max_length - length
            if pad_width > 0:
                padded_max_runs.append(torch_module.nn.functional.pad(max_run_tensor, (0, pad_width)))
            else:
                padded_max_runs.append(max_run_tensor)
    stacked_totals = torch_module.stack(padded_totals, dim=0)
    stacked_shared = torch_module.stack(padded_shared, dim=0)
    components = [stacked_totals, stacked_shared]
    if max_run_tensors is not None and padded_max_runs is not None:
        components.append(torch_module.stack(padded_max_runs, dim=0))
    combined = torch_module.stack(components, dim=1)
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
    compute_backend: MatrixPairComputeBackend,
    matrix_value_semantics: str,
    calculations: tuple[str, ...],
    gene_ranges: list[GeneRangeSpec],
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
    max_vector_length = max(spec.matrix_length for spec in genomes)
    global_block_size, _ = _plan_chunk_sizes(
        vector_length=max_vector_length,
        remaining_targets=len(samples),
        dtype_name=dtype_name,
        memory_limit_bytes=memory_limit_bytes,
        backend_kind=compute_backend.kind,
    )
    global_block_size = max(1, global_block_size)
    samples_by_id = {sample_idx: sample_name for sample_idx, sample_name in samples}
    gene_ranges_by_genome = _group_gene_ranges_by_genome(gene_ranges)
    loop = asyncio.get_running_loop()
    target_loader_executor = _build_matrix_io_executor(loader_executor_kind)
    anchor_loader_executor = _build_matrix_io_executor(loader_executor_kind)
    writer_executor = _build_matrix_io_executor(writer_executor_kind)
    pending_write_futures: deque[tuple[asyncio.Future[int], list[dict[str, object]]]] = deque()

    def load_target_block_sync(
        spec: GenomeSpec,
        block_start: int,
        block_rows: list[tuple[int, str]],
    ) -> tuple[np.ndarray, list[str], np.ndarray, object]:
        return _load_target_queue_block_for_hdf5_torch(
            matrix_hdf5_file=matrix_db_file,
            genome_idx=spec.genome_idx,
            block_rows=block_rows,
            block_start=block_start,
            matrix_length=spec.matrix_length,
            dtype_name=dtype_name,
            matrix_value_semantics=matrix_value_semantics,
            torch_module=compute_backend.torch,
        )

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
    target_queue: deque[tuple[int, tuple[np.ndarray, list[str], np.ndarray, object]]] = deque()
    next_target_prefetch_idx = 1
    pending_target_future: Optional[asyncio.Future[tuple[int, tuple[np.ndarray, list[str], np.ndarray, object]]]] = None

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
        prefetch_block_start, prefetch_block_rows, prefetch_spec = work_units[unit_index]
        pending_target_future = loop.run_in_executor(
            target_loader_executor,
            _load_target_prefetch_unit_for_hdf5_torch,
            unit_index,
            matrix_db_file,
            prefetch_spec.genome_idx,
            prefetch_block_rows,
            prefetch_block_start,
            prefetch_spec.matrix_length,
            dtype_name,
            matrix_value_semantics,
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
                block_ids, block_names, zero_matrix, target_matrices = load_target_block_sync(
                    spec,
                    block_start,
                    block_rows,
                )
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
            genome_gene_ranges = gene_ranges_by_genome.get(spec.genome_idx, [])
            tile_targets, _ = _plan_chunk_sizes(
                vector_length=spec.matrix_length,
                remaining_targets=len(block_ids),
                dtype_name=dtype_name,
                memory_limit_bytes=memory_limit_bytes,
                backend_kind=compute_backend.kind,
            )
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
            del target_matrices
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
                max_run_tensors = None
                if batch_device_results and batch_device_results[0]["ibs_shared_mask"] is not None:
                    max_run_tensors = [
                        _max_ibs_from_shared_mask_torch(
                            torch_module=compute_backend.torch,
                            shared_mask=item["ibs_shared_mask"],
                        )
                        for item in batch_device_results
                    ]
                combined_np = _download_torch_result_tensor_batch(
                    compute_backend=compute_backend,
                    totals_tensors=[item["total_positions"] for item in batch_device_results],
                    shared_tensors=[item["share_allele_pos"] for item in batch_device_results],
                    max_run_tensors=max_run_tensors,
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
                            "max_consecutive_length": None
                            if max_run_tensors is None
                            else combined_np[result_idx, 2, :valid_count].astype(np.int64, copy=True),
                            "gene_names": list(item.get("gene_names") or []),
                            "gene_total_positions": None
                            if item.get("gene_total_positions") is None
                            else np.asarray(item["gene_total_positions"], dtype=np.int64)[:, :valid_count].copy(),
                            "gene_share_allele_pos": None
                            if item.get("gene_share_allele_pos") is None
                            else np.asarray(item["gene_share_allele_pos"], dtype=np.int64)[:, :valid_count].copy(),
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

            anchor_queue: deque[tuple[int, str, object]] = deque()
            next_anchor_offset = 0
            pending_anchor_future: Optional[asyncio.Future[list[tuple[int, str, object]]]] = None

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
                    _load_anchor_queue_batch_for_hdf5_torch,
                    matrix_db_file,
                    spec.genome_idx,
                    batch_rows,
                    next_anchor_offset - len(batch_rows),
                    spec.matrix_length,
                    dtype_name,
                    zero_matrix,
                    matrix_value_semantics,
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
                total_inc_torch, shared_inc_torch, ibs_shared_mask, gene_total_inc, gene_shared_inc = _compare_anchor_against_target_chunk_torch_device(
                    compute_backend=compute_backend,
                    anchor_torch=anchor_torch,
                    target_torch=target_torch,
                    vector_length=spec.matrix_length,
                    matrix_value_semantics=matrix_value_semantics,
                    need_ibs="ibs" in calculations,
                    gene_ranges=genome_gene_ranges if "gene" in calculations else None,
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
                        "ibs_shared_mask": None
                        if ibs_shared_mask is None
                        else ibs_shared_mask[:, missing_positions].contiguous(),
                        "gene_names": [gene.gene for gene in genome_gene_ranges] if "gene" in calculations else [],
                        "gene_total_positions": None
                        if gene_total_inc is None
                        else gene_total_inc[:, missing_positions].astype(np.int64, copy=True),
                        "gene_share_allele_pos": None
                        if gene_shared_inc is None
                        else gene_shared_inc[:, missing_positions].astype(np.int64, copy=True),
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
                total_inc_torch, shared_inc_torch, ibs_shared_mask, gene_total_inc, gene_shared_inc = _compare_anchor_against_target_chunk_torch_device(
                    compute_backend=compute_backend,
                    anchor_torch=target_torch[:, :, local_anchor_pos],
                    target_torch=target_torch[:, :, local_anchor_pos + 1:],
                    vector_length=spec.matrix_length,
                    matrix_value_semantics=matrix_value_semantics,
                    need_ibs="ibs" in calculations,
                    gene_ranges=genome_gene_ranges if "gene" in calculations else None,
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
                        "ibs_shared_mask": None
                        if ibs_shared_mask is None
                        else ibs_shared_mask[:, missing_positions].contiguous(),
                        "gene_names": [gene.gene for gene in genome_gene_ranges] if "gene" in calculations else [],
                        "gene_total_positions": None
                        if gene_total_inc is None
                        else gene_total_inc[:, missing_positions].astype(np.int64, copy=True),
                        "gene_share_allele_pos": None
                        if gene_shared_inc is None
                        else gene_shared_inc[:, missing_positions].astype(np.int64, copy=True),
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
            del target_torch
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
    )


def _matrix_compare_reuse_target_chunks_torch(
    matrix_db_file: Path,
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
    compute_backend: MatrixPairComputeBackend,
    matrix_value_semantics: str,
    calculations: tuple[str, ...],
    gene_ranges: list[GeneRangeSpec],
    emit_writer_logs: bool,
    run_start_time: float,
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    return asyncio.run(
        _matrix_compare_reuse_target_chunks_torch_async(
            matrix_db_file=matrix_db_file,
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
            compute_backend=compute_backend,
            matrix_value_semantics=matrix_value_semantics,
            calculations=calculations,
            gene_ranges=gene_ranges,
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
    backend: str = "numpy",
    calculate: Optional[str] = "all",
    emit_writer_logs: bool = False,
    progress_callback: Optional[CompareProgressCallback] = None,
) -> MatrixCompareSummary:
    """Run resumable all-vs-all matrix comparison on the current matrix store.

    Results are written into a DuckDB compare database. Re-running the same
    command against an existing compare database resumes from the remaining
    uncompleted pairs.
    """
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
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    run_start_time = time.monotonic()
    if not _is_hdf5_matrix_store(matrix_db_file):
        raise ValueError(
            "matrix-compare expects the current HDF5-backed matrix store format. "
            "If you have a legacy DuckDB matrix database, convert it first with "
            "'zipstrain utilities matrix-db-to-hdf5'."
        )
    compare_conn: Optional[duckdb.DuckDBPyConnection] = None
    try:
        genome_scope = None if genome == "all" else genome
        gene_ranges: list[GeneRangeSpec] = []
        metadata = _load_matrix_hdf5_metadata(matrix_db_file)
        metadata["input_format"] = "hdf5"
        matrix_value_semantics = metadata.get("matrix_value_semantics")
        if matrix_value_semantics != FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS:
            raise ValueError(
                "This matrix store uses an unsupported storage layout. "
                "Rebuild it with 'zipstrain utilities build-matrix-db'."
            )
        samples = _load_matrix_hdf5_samples(matrix_db_file)
        genomes = _load_matrix_hdf5_genomes(matrix_db_file, genome=genome_scope)
        if not genomes:
            raise ValueError(f"No scaffolds found for genome scope: {genome}")
        genome_scaffolds = _load_matrix_hdf5_genome_scaffolds(matrix_db_file, genome=genome_scope)
        gene_ranges = _load_matrix_hdf5_gene_ranges(matrix_db_file, genome=genome_scope)
        calculations = parse_matrix_calculations(
            calculate,
            include_gene_from_all=bool(gene_ranges),
        )
        if "gene" in calculations and not gene_ranges:
            raise ValueError(
                "Gene ANI was requested, but this matrix store does not contain gene annotations. "
                "Rebuild it with --gene-range-table."
            )
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
            )

        if compute_backend.kind == "torch":
            compare_conn.close()
            compare_conn = None
            return _matrix_compare_reuse_target_chunks_torch(
                matrix_db_file=matrix_db_file,
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
                compute_backend=compute_backend,
                matrix_value_semantics=matrix_value_semantics,
                calculations=calculations,
                gene_ranges=gene_ranges,
                emit_writer_logs=emit_writer_logs,
                run_start_time=run_start_time,
                progress_callback=progress_callback,
            )

        completed_work = 0
        written_rows = 0
        target_chunks = 0
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
                remaining_target_row_indices = np.array(
                    [
                        row_idx
                        for row_idx in range(anchor_offset + 1, len(samples))
                        if (sample_1_idx, samples[row_idx][0]) not in completed_pairs_for_genome
                    ],
                    dtype=np.int64,
                )
                target_ids_all = np.array([sample_idx for sample_idx, _sample_name in remaining_target_rows], dtype=np.int64)
                target_names_all = [sample_name for _sample_idx, sample_name in remaining_target_rows]
                zero_matrix = np.zeros((spec.matrix_length, 4), dtype=COUNT_DTYPES[dtype_name])
                dataset = _Hdf5GenomeMatrixNumpyDataset(
                    matrix_hdf5_file=matrix_db_file,
                    genome_idx=spec.genome_idx,
                    matrix_length=spec.matrix_length,
                    dtype_name=dtype_name,
                )
                try:
                    anchor_matrix = dataset.get_row(anchor_offset)
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
                        max_targets, _ = _plan_chunk_sizes(
                            vector_length=spec.matrix_length,
                            remaining_targets=len(target_ids_all) - target_offset,
                            dtype_name=dtype_name,
                            memory_limit_bytes=memory_limit_bytes,
                            backend_kind=compute_backend.kind,
                        )
                        chunk_ids = target_ids_all[target_offset: target_offset + max_targets]
                        chunk_names = target_names_all[target_offset: target_offset + max_targets]
                        chunk_rows = remaining_target_row_indices[target_offset: target_offset + max_targets]
                        target_matrices = dataset.load_indices(chunk_rows)
                        if "ibs" in calculations:
                            totals_chunk, shared_chunk, shared_mask = _compare_tile_presence_numpy_with_mask(
                                anchor_matrix=anchor_matrix,
                                target_matrices=target_matrices,
                            )
                        else:
                            totals_chunk, shared_chunk = _compare_tile_presence_numpy(
                                anchor_matrix=anchor_matrix,
                                target_matrices=target_matrices,
                            )
                        ibs_max = np.zeros(len(chunk_ids), dtype=np.int64) if "ibs" in calculations else None
                        current_runs = np.zeros(len(chunk_ids), dtype=np.int64) if "ibs" in calculations else None
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
                finally:
                    dataset.close()

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
            )
    finally:
        if compare_conn is not None:
            compare_conn.close()


def export_matrix_compare_parquet(
    matrix_compare_db_file: Path,
    output_file: Path,
    table: str = "genome",
) -> Path:
    """Export matrix compare results from DuckDB into parquet.

    Parameters
    ----------
    matrix_compare_db_file:
        DuckDB compare database produced by :func:`matrix_compare`.
    output_file:
        Destination parquet file.
    table:
        ``"genome"`` for genome-level comparison rows or ``"gene"`` for
        gene-level ANI rows.
    """
    matrix_compare_db_file = Path(matrix_compare_db_file).resolve()
    output_file = Path(output_file).resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    conn = duckdb.connect(str(matrix_compare_db_file), read_only=True)
    try:
        if table == "genome":
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
        elif table == "gene":
            if not _matrix_compare_table_exists(conn, "matrix_compare_gene_results"):
                raise ValueError(
                    "Gene export was requested, but matrix_compare_gene_results is not present in the compare database."
                )
            gene_row_count = int(conn.execute("SELECT count(*) FROM matrix_compare_gene_results").fetchone()[0])
            if gene_row_count == 0:
                raise ValueError(
                    "Gene export was requested, but the compare database does not contain any gene comparison rows."
                )
            query = """
                SELECT
                  sample_1,
                  sample_2,
                  genome,
                  gene,
                  gene_pop_ani
                FROM matrix_compare_gene_results
                ORDER BY sample_idx_1, sample_idx_2, genome_idx, gene
            """
        else:
            raise ValueError(f"Unsupported export table '{table}'. Expected one of: genome, gene.")
        conn.execute(
            f"COPY ({query}) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
            [str(output_file)],
        )
    finally:
        conn.close()
    return output_file
