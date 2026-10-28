"""Blocked-GEMM genome-level popANI / conANI for HDF5 matrix stores.

This is the ``compare_kernel='block'`` kernel. The ``anchor`` kernel in
:mod:`zipstrain.matrix_pairs` compares one anchor
sample against blocks of targets position by position. That is bandwidth
bound: every anchor/target pair re-reads both rows. For genome-level ANI the
per-pair quantities are plain sums over positions, so a whole block of
sample pairs can be computed with dense matrix multiplications instead.

Every supported matrix layout is reduced to a 4-bit allele mask per position
(bit ``b`` set when allele ``b`` is present after the build-time coverage
filter; for conANI, when allele ``b`` is a consensus allele). With ``X`` the
``samples x positions`` mask matrix and ``C = (X != 0)``:

``total_positions = C @ C.T``
    positions covered in both samples.

``share_allele_pos = total_positions - mismatch``
    a covered-in-both position is a mismatch when the two masks share no bit.
    Mismatches can only happen at *informative* positions, where the covered
    samples of the whole cohort do not share a common allele. On those
    positions ``shared = sum_b A_b @ A_b.T`` corrected by inclusion-exclusion
    over multi-allele masks (``A_m = (X & m) == m``), and
    ``mismatch = C_I @ C_I.T - shared_I``.

All GEMMs run in float32, which is exact for integer sums below ``2**24``;
partial sums are flushed into int64 host accumulators before that bound
can be reached. Results are bit-identical to the anchor kernel.

Only genome-level ANI is produced; IBS (``max_consecutive_length``) and
gene-level ANI depend on position order or per-gene masks and must use the
anchor kernel.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import duckdb
import numpy as np
import pyarrow as pa

from zipstrain import matrix_pairs as mp

GEMM_SUPPORTED_ANI_KINDS = ("popani", "conani")
#: float32 represents every integer up to 2**24 exactly.
GEMM_EXACT_FP32_LIMIT = 1 << 24
#: Default number of positions multiplied per GEMM call.
GEMM_POSITION_CHUNK_DEFAULT = 32_768
#: Largest sample block considered by the planner. Diagonal block pairs are
#: computed in full, so larger blocks add ~n*B/2 redundant pair counts.
GEMM_MAX_BLOCK_SAMPLES = 4_096
#: Bound raw count reads and their mask-conversion temporaries independently
#: of the size of the final, one-byte-per-position mask cache.
GEMM_READ_BUFFER_BYTES_DEFAULT = 128 * 1024**2

_ALLELE_BITS = tuple(int(bit) for bit in mp.BITMASK_BASE_BITS)
# Inclusion-exclusion terms for |{b : b in mask_i and b in mask_j}| > 0:
# + singles - pairs + triples - quad. Each entry is (sign, combined mask).
_COMBO_TERMS = tuple(
    (sign, int(sum(combo)))
    for size, sign in ((2, -1.0), (3, 1.0), (4, -1.0))
    for combo in itertools.combinations(_ALLELE_BITS, size)
)
_POPCOUNT = np.array([bin(value).count("1") for value in range(256)], dtype=np.uint8)
_FULL_MASK = np.uint8(sum(_ALLELE_BITS))


@dataclass(frozen=True)
class GemmPlan:
    """Sample block size and position chunk used by the GEMM kernel."""

    block_samples: int
    position_chunk: int
    cache_all_rows: bool


def rows_to_allele_bitmask(
    rows: np.ndarray,
    *,
    matrix_value_semantics: str,
    ani_kind: str,
) -> np.ndarray:
    """Convert a batch of stored matrix rows to ``(B, L)`` uint8 allele masks.

    ``rows`` is what ``_Hdf5GenomeMatrixNumpyDataset._load_batch_numpy``
    returns: ``(B, L)`` packed bitmasks, or ``(B, L, 4)`` counts / presence.
    """
    if matrix_value_semantics == mp.BITMASK_MATRIX_VALUE_SEMANTICS:
        if ani_kind != "popani":
            raise ValueError("Bitmask matrix stores support only popANI.")
        return np.ascontiguousarray(rows, dtype=np.uint8)
    if rows.ndim != 3 or rows.shape[-1] != len(_ALLELE_BITS):
        raise ValueError(f"Expected (samples, positions, 4) matrix rows, got shape {rows.shape}.")
    if ani_kind == "conani":
        if matrix_value_semantics != mp.COUNT_MATRIX_VALUE_SEMANTICS:
            raise ValueError("conANI requires a count matrix store.")
        maxima = rows.max(axis=2, keepdims=True)
        selected = (rows == maxima) & (maxima > 0)
    elif ani_kind == "popani":
        if matrix_value_semantics not in (
            mp.COUNT_MATRIX_VALUE_SEMANTICS,
            mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS,
        ):
            raise ValueError(f"Unsupported matrix value semantics {matrix_value_semantics!r}.")
        selected = rows > 0
    else:
        raise ValueError(f"compare_kernel='block' does not support ani_method {ani_kind!r}.")
    out = np.zeros(rows.shape[:2], dtype=np.uint8)
    for channel, bit in enumerate(_ALLELE_BITS):
        out |= selected[:, :, channel].astype(np.uint8) * np.uint8(bit)
    return out


def _plan_gemm_blocks(
    *,
    sample_count: int,
    matrix_length: int,
    memory_limit_bytes: int,
    device_memory_bytes: Optional[int] = None,
    informative_count: Optional[int] = None,
) -> GemmPlan:
    """Pick the largest sample block that fits the memory budget.

    The whole genome matrix is cached on the host when it fits in a third of
    the budget; row blocks are then zero-copy views. Otherwise two uint8 row
    blocks are resident. Each block also keeps a copy of its informative
    columns. Compute uses ``B x W`` float32 tiles (uploaded per position
    chunk for torch backends) and three ``B x B`` accumulators, plus int64
    host copies of the accumulators.
    """
    length = int(matrix_length)
    width = max(1, min(length, GEMM_POSITION_CHUNK_DEFAULT))
    informative = length if informative_count is None else int(informative_count)
    budget = int(memory_limit_bytes * 0.7)
    full_matrix = int(sample_count) * length
    cache_all = full_matrix <= memory_limit_bytes // 3
    if cache_all:
        budget -= full_matrix
    device_budget = None if device_memory_bytes is None else int(device_memory_bytes * 0.7)

    def host_bytes(block: int) -> int:
        rows = 0 if cache_all else 2 * block * length
        return rows + 2 * block * informative + 12 * block * width * 4 + 3 * block * block * 12

    def device_bytes(block: int) -> int:
        return 12 * block * width * 4 + 3 * block * block * 4

    block = max(1, min(int(sample_count), GEMM_MAX_BLOCK_SAMPLES))
    while block > 1 and (
        host_bytes(block) > budget
        or (device_budget is not None and device_bytes(block) > device_budget)
    ):
        block = max(1, block // 2) if block > 64 else block - 1
    return GemmPlan(block_samples=block, position_chunk=width, cache_all_rows=cache_all)


class _GemmOps:
    """Tiny numpy/torch abstraction for the handful of tensor ops we need."""

    def __init__(self, compute_backend) -> None:
        self.kind = compute_backend.kind
        self.torch = compute_backend.torch
        self.device = compute_backend.device
        self._saved_precision = None
        self._saved_tf32 = None
        self._product = None

    def __enter__(self) -> "_GemmOps":
        if self.kind == "torch":
            torch = self.torch
            # TF32 / reduced precision would break exactness.
            self._saved_precision = torch.get_float32_matmul_precision()
            torch.set_float32_matmul_precision("highest")
            if hasattr(torch.backends, "cuda"):
                self._saved_tf32 = torch.backends.cuda.matmul.allow_tf32
                torch.backends.cuda.matmul.allow_tf32 = False
        return self

    def __exit__(self, *exc) -> None:
        if self.kind == "torch":
            self.torch.set_float32_matmul_precision(self._saved_precision)
            if self._saved_tf32 is not None:
                self.torch.backends.cuda.matmul.allow_tf32 = self._saved_tf32

    def tile(self, host_block: np.ndarray, start: int, stop: int):
        """Columns ``start:stop`` of a host uint8 block, on the compute device.

        Torch backends upload one contiguous ``B x W`` tile at a time, so the
        device never holds whole-genome row blocks.
        """
        view = host_block[:, start:stop]
        if self.kind == "numpy":
            return view
        return self.torch.from_numpy(np.ascontiguousarray(view)).to(self.device)

    def zeros(self, rows: int, cols: int):
        if self.kind == "numpy":
            return np.zeros((rows, cols), dtype=np.float32)
        return self.torch.zeros((rows, cols), dtype=self.torch.float32, device=self.device)

    def covered(self, tile):
        if self.kind == "numpy":
            return (tile != 0).astype(np.float32)
        return (tile != 0).to(self.torch.float32)

    def has_bit(self, tile, bit: int):
        if self.kind == "numpy":
            return ((tile & np.uint8(bit)) != 0).astype(np.float32)
        return ((tile & bit) != 0).to(self.torch.float32)

    def has_all_bits(self, tile, bits: int):
        if self.kind == "numpy":
            return ((tile & np.uint8(bits)) == np.uint8(bits)).astype(np.float32)
        return ((tile & bits) == bits).to(self.torch.float32)

    def addmm_(self, acc, left, right, alpha: float = 1.0) -> None:
        """``acc += alpha * left @ right.T`` in float32."""
        if self.kind == "numpy":
            if self._product is None or self._product.shape != acc.shape:
                self._product = np.empty(acc.shape, dtype=np.float32)
            np.matmul(left, right.T, out=self._product)
            if alpha == 1.0:
                acc += self._product
            else:
                acc -= self._product if alpha == -1.0 else -np.float32(alpha) * self._product
            return
        acc.addmm_(left, right.T, alpha=alpha)

    def to_int64(self, acc) -> np.ndarray:
        if self.kind == "numpy":
            return np.rint(acc).astype(np.int64)
        return self.torch.round(acc).to(self.torch.int64).cpu().numpy()

    def device_memory_bytes(self) -> Optional[int]:
        if self.kind == "torch" and str(self.device).startswith("cuda"):
            free, _total = self.torch.cuda.mem_get_info(self.torch.device(self.device))
            return int(free)
        return None


class _ExactAccumulator:
    """float32 device accumulator with periodic exact flushes to int64."""

    def __init__(self, ops: _GemmOps, rows: int, cols: int, magnitude_per_position: int) -> None:
        self.ops = ops
        self.rows = rows
        self.cols = cols
        self.magnitude = int(magnitude_per_position)
        self.acc = ops.zeros(rows, cols)
        self.host = np.zeros((rows, cols), dtype=np.int64)
        self.pending_positions = 0

    def reserve(self, positions: int) -> None:
        """Flush first if adding ``positions`` could exceed exact float32 range."""
        if (self.pending_positions + positions) * self.magnitude >= GEMM_EXACT_FP32_LIMIT:
            self.flush()
        self.pending_positions += positions

    def flush(self) -> None:
        if self.pending_positions == 0:
            return
        self.host += self.ops.to_int64(self.acc)
        self.acc = self.ops.zeros(self.rows, self.cols)
        self.pending_positions = 0

    def result(self) -> np.ndarray:
        self.flush()
        return self.host


def _block_pair_counts(
    ops: _GemmOps,
    *,
    block_i,
    block_j,
    inf_i,
    inf_j,
    multi_positions: np.ndarray,
    matrix_length: int,
    position_chunk: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return int64 ``(total_positions, share_allele_pos)`` for two row blocks.

    All blocks are host uint8 arrays. ``inf_i`` / ``inf_j`` are the
    informative-position columns of the blocks; ``multi_positions`` indexes
    the informative columns that carry at least one multi-allele mask anywhere
    in the cohort.
    """
    rows_i = int(block_i.shape[0])
    rows_j = int(block_j.shape[0])
    width = max(1, min(int(position_chunk), GEMM_EXACT_FP32_LIMIT // 8))
    # Diagonal block pairs multiply a block with itself: convert each tile once.
    same = block_i is block_j

    def tiles(host_i, host_j, start, stop):
        tile_i = ops.tile(host_i, start, stop)
        return tile_i, (tile_i if same else ops.tile(host_j, start, stop))

    def pair(convert, tile_i, tile_j, *args):
        left = convert(tile_i, *args)
        return left, (left if same else convert(tile_j, *args))

    total = _ExactAccumulator(ops, rows_i, rows_j, magnitude_per_position=1)
    for start in range(0, matrix_length, width):
        stop = min(start + width, matrix_length)
        total.reserve(stop - start)
        ops.addmm_(total.acc, *pair(ops.covered, *tiles(block_i, block_j, start, stop)))
    totals = total.result()

    informative = int(inf_i.shape[1])
    if informative == 0:
        return totals, totals.copy()
    # mismatch = covered_both_I - shared_I; inclusion-exclusion partial sums
    # are bounded by 4 per position (4 single-bit terms of magnitude 1).
    mismatch = _ExactAccumulator(ops, rows_i, rows_j, magnitude_per_position=5)
    for start in range(0, informative, width):
        stop = min(start + width, informative)
        tile_i, tile_j = tiles(inf_i, inf_j, start, stop)
        mismatch.reserve(stop - start)
        ops.addmm_(mismatch.acc, *pair(ops.covered, tile_i, tile_j))
        for bit in _ALLELE_BITS:
            ops.addmm_(mismatch.acc, *pair(ops.has_bit, tile_i, tile_j, bit), alpha=-1.0)
    if multi_positions.size:
        multi_i = np.ascontiguousarray(inf_i[:, multi_positions])
        multi_j = multi_i if same else np.ascontiguousarray(inf_j[:, multi_positions])
        for start in range(0, int(multi_positions.size), width):
            stop = min(start + width, int(multi_positions.size))
            tile_i, tile_j = tiles(multi_i, multi_j, start, stop)
            # 11 combo terms of magnitude <= 1 each per position.
            mismatch.reserve((stop - start) * 3)
            for sign, bits in _COMBO_TERMS:
                # shared gets ``sign`` -> mismatch gets ``-sign``.
                ops.addmm_(mismatch.acc, *pair(ops.has_all_bits, tile_i, tile_j, bits), alpha=-sign)
    mismatches = mismatch.result()
    return totals, totals - mismatches


class _GenomeRows:
    """Loads (and optionally caches) allele-mask row blocks for one genome."""

    def __init__(
        self,
        *,
        matrix_db_file: Path,
        spec,
        dtype_name: str,
        matrix_value_semantics: str,
        ani_kind: str,
        sample_count: int,
        cache_all_rows: bool,
        read_memory_limit_bytes: int = GEMM_READ_BUFFER_BYTES_DEFAULT,
    ) -> None:
        self.dataset = mp._Hdf5GenomeMatrixNumpyDataset(
            matrix_hdf5_file=matrix_db_file,
            genome_idx=spec.genome_idx,
            matrix_length=spec.matrix_length,
            dtype_name=dtype_name,
            matrix_value_semantics=matrix_value_semantics,
        )
        self.semantics = matrix_value_semantics
        self.ani_kind = ani_kind
        self.sample_count = sample_count
        self.matrix_length = int(spec.matrix_length)
        self.read_memory_limit_bytes = int(read_memory_limit_bytes)
        if self.read_memory_limit_bytes < 1:
            raise ValueError("read_memory_limit_bytes must be >= 1")
        self._cache: Optional[np.ndarray] = None
        if cache_all_rows:
            try:
                self._cache = self._read(0, sample_count)
            except BaseException:
                self.dataset.close()
                raise

    def _read(self, start: int, stop: int) -> np.ndarray:
        dataset = self.dataset
        dataset._ensure_open()
        if self.semantics == mp.BITMASK_MATRIX_VALUE_SEMANTICS:
            if dataset._matrix_dataset is not None:
                return np.asarray(dataset._matrix_dataset[start:stop, ...], dtype=np.uint8)
            return dataset._load_batch_numpy(np.arange(start, stop, dtype=np.int64))

        # The planner budgets for packed masks. Reading all count rows before
        # converting them would temporarily require 8--32 times as much memory.
        # Include raw counts, consensus maxima, selection masks and conversion
        # temporaries in a separate bounded workspace, including uncached reads.
        itemsize = np.dtype(dataset.numpy_dtype).itemsize
        bytes_per_position = 5 * itemsize + 16
        if dataset._matrix_dataset is None:
            # Sparse expansion also reads up to four int64 indices and four
            # values per position alongside the dense count buffer.
            bytes_per_position = max(bytes_per_position, 8 * itemsize + 32)
        workspace_bytes = min(self.read_memory_limit_bytes, GEMM_READ_BUFFER_BYTES_DEFAULT)
        positions_per_read = workspace_bytes // bytes_per_position
        if positions_per_read < 1:
            raise MemoryError("GEMM memory budget cannot fit one allele-count position.")
        if dataset._matrix_dataset is None and positions_per_read < self.matrix_length:
            if self.matrix_length * bytes_per_position > self.read_memory_limit_bytes:
                raise MemoryError(
                    "GEMM memory budget cannot fit one expanded sparse count row plus "
                    "mask-conversion buffers. Increase --memory-limit-gb."
                )
            # Sparse rows expand as a whole. Permit one row above the preferred
            # 128 MiB chunk size when it still fits the reserved workspace.
            positions_per_read = self.matrix_length
        batch_samples = max(1, positions_per_read // max(1, self.matrix_length))
        batch_positions = max(1, min(self.matrix_length, positions_per_read))
        masks = np.empty((stop - start, self.matrix_length), dtype=np.uint8)
        for row_start in range(start, stop, batch_samples):
            row_stop = min(row_start + batch_samples, stop)
            for col_start in range(0, self.matrix_length, batch_positions):
                col_stop = min(col_start + batch_positions, self.matrix_length)
                if dataset._matrix_dataset is not None:
                    raw = np.asarray(
                        dataset._matrix_dataset[row_start:row_stop, col_start:col_stop, :],
                        dtype=dataset.numpy_dtype,
                    )
                else:
                    raw = dataset._load_batch_numpy(np.arange(row_start, row_stop, dtype=np.int64))
                masks[row_start - start:row_stop - start, col_start:col_stop] = rows_to_allele_bitmask(
                    raw, matrix_value_semantics=self.semantics, ani_kind=self.ani_kind,
                )
                del raw
        return masks

    def block(self, start: int, stop: int) -> np.ndarray:
        if self._cache is not None:
            return self._cache[start:stop]
        return self._read(start, stop)

    def close(self) -> None:
        self._cache = None
        self.dataset.close()


def _informative_positions(rows: _GenomeRows, block_samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Cohort-wide informative positions and, among them, multi-allele positions.

    A position is informative when at least two samples cover it and the
    covered masks share no common allele bit. Non-informative positions can
    never be mismatches for any pair, so restricting the mismatch GEMMs to the
    informative set is exact.
    """
    length = rows.matrix_length
    common = np.full(length, _FULL_MASK, dtype=np.uint8)
    covered_count = np.zeros(length, dtype=np.int32)
    any_multi = np.zeros(length, dtype=bool)
    for start in range(0, rows.sample_count, block_samples):
        block = rows.block(start, min(start + block_samples, rows.sample_count))
        # Column chunks keep the boolean/uint8 temporaries to ~256 MB.
        width = max(1, (1 << 28) // max(1, block.shape[0]))
        for col in range(0, length, width):
            tile = block[:, col:col + width]
            covered = tile != 0
            common[col:col + width] &= np.bitwise_and.reduce(np.where(covered, tile, _FULL_MASK), axis=0)
            covered_count[col:col + width] += covered.sum(axis=0, dtype=np.int32)
            any_multi[col:col + width] |= (_POPCOUNT[tile] >= 2).any(axis=0)
    informative = np.flatnonzero((common == 0) & (covered_count >= 2))
    multi = np.flatnonzero(any_multi[informative])
    return informative.astype(np.int64), multi.astype(np.int64)


def _completed_by_block_pair(
    completed_pairs: set[tuple[int, int]],
    row_of_sample: dict[int, int],
    block_samples: int,
) -> dict[tuple[int, int], list[tuple[int, int]]]:
    buckets: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for sample_a, sample_b in completed_pairs:
        row_a = row_of_sample.get(int(sample_a))
        row_b = row_of_sample.get(int(sample_b))
        if row_a is None or row_b is None or row_a == row_b:
            continue
        if row_a > row_b:
            row_a, row_b = row_b, row_a
        key = (row_a // block_samples, row_b // block_samples)
        buckets.setdefault(key, []).append((row_a, row_b))
    return buckets


def _insert_completed_arrays(
    compare_conn: duckdb.DuckDBPyConnection,
    sample_idx_1: np.ndarray,
    sample_idx_2: np.ndarray,
    genome_idx: int,
) -> None:
    if sample_idx_1.size == 0:
        return
    table = pa.table(
        {
            "sample_idx_1": pa.array(sample_idx_1, type=pa.int64()),
            "sample_idx_2": pa.array(sample_idx_2, type=pa.int64()),
            "genome_idx": pa.array(np.full(sample_idx_1.size, genome_idx, dtype=np.int64), type=pa.int64()),
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


def _result_table(
    *,
    compare_id_1: np.ndarray,
    compare_id_2: np.ndarray,
    names_1: np.ndarray,
    names_2: np.ndarray,
    genome_id: int,
    genome: str,
    totals: np.ndarray,
    shared: np.ndarray,
) -> pa.Table:
    """Arrow rows in ``matrix_compare_result_db_schema`` order.

    Pairs are canonicalized like the pairwise kernel: ``sample_idx_1`` is the
    smaller compare id and ``sample_1`` its name.
    """
    swap = compare_id_1 > compare_id_2
    first_id = np.where(swap, compare_id_2, compare_id_1)
    second_id = np.where(swap, compare_id_1, compare_id_2)
    first_name = np.where(swap, names_2, names_1)
    second_name = np.where(swap, names_1, names_2)
    rows = int(totals.size)
    return pa.Table.from_arrays(
        [
            pa.array(first_id, type=pa.int64()),
            pa.array(second_id, type=pa.int64()),
            pa.array(first_name.tolist(), type=pa.string()),
            pa.array(second_name.tolist(), type=pa.string()),
            pa.array(np.full(rows, genome_id, dtype=np.int64), type=pa.int64()),
            pa.array([genome] * rows, type=pa.string()),
            pa.array(totals, type=pa.int64()),
            pa.array(shared, type=pa.int64()),
            pa.array(shared / totals * 100.0, type=pa.float64()),
            pa.nulls(rows, type=pa.int64()),
        ],
        schema=mp.matrix_compare_result_db_schema(),
    )


def matrix_compare_gemm(
    *,
    matrix_db_file: Path,
    output_file: Path,
    compare_conn: duckdb.DuckDBPyConnection,
    compare_ids,
    samples: list[tuple[int, str]],
    genomes: list,
    genome_scaffolds: list,
    completed_pairs_by_genome: dict[int, set[tuple[int, int]]],
    requested_pairs: int,
    total_work: int,
    metadata: dict[str, str],
    matrix_value_semantics: str,
    ani_method: str,
    compute_backend,
    memory_limit_bytes: int,
    memory_limit_gb: float,
    min_cov: int,
    genome_scope: Optional[str],
    emit_writer_logs: bool,
    run_start_time: float,
    progress_callback: Optional[Callable[[dict[str, object]], None]],
):
    """Run the GEMM kernel on a prepared compare DB; see module docstring.

    Called by :func:`zipstrain.matrix_pairs.matrix_compare` after the compare
    DB, id catalog and resume state are prepared. Each sample-block pair is
    written in one transaction (result rows with ``total_positions > 0`` plus
    completion keys for every pair in the block pair), so interrupted runs
    resume exactly like the pairwise kernel.
    """
    ani_kind, _cos = mp.parse_matrix_ani_method(ani_method)
    if ani_kind not in GEMM_SUPPORTED_ANI_KINDS:
        raise ValueError(
            f"compare_kernel='block' supports ani_method popani or conani, not {ani_method!r}."
        )
    dtype_name = mp._matrix_dtype_name_from_metadata(metadata)
    sample_count = len(samples)
    sample_matrix_ids = np.array([int(sample_idx) for sample_idx, _name in samples], dtype=np.int64)
    sample_names = np.array([str(name) for _idx, name in samples], dtype=object)
    compare_sample_ids = np.asarray(compare_ids.samples, dtype=np.int64)[sample_matrix_ids]
    row_of_sample = {int(sample_idx): row for row, sample_idx in enumerate(sample_matrix_ids.tolist())}

    completed_work = 0
    written_rows = 0
    block_pairs_done = 0
    sample_blocks_total = 0

    with _GemmOps(compute_backend) as ops:
        device_memory = ops.device_memory_bytes()
        for spec in genomes:
            completed = completed_pairs_by_genome.get(spec.genome_idx, set())
            if len(completed) >= sample_count * (sample_count - 1) // 2:
                continue
            # Conservative plan (informative count unknown) for the cache decision
            # and the informative-position pass, then the final plan.
            plan = _plan_gemm_blocks(
                sample_count=sample_count,
                matrix_length=spec.matrix_length,
                memory_limit_bytes=memory_limit_bytes,
                device_memory_bytes=device_memory,
            )
            compare_genome_id = int(compare_ids.genomes[spec.genome_idx])
            rows = _GenomeRows(
                matrix_db_file=matrix_db_file,
                spec=spec,
                dtype_name=dtype_name,
                matrix_value_semantics=matrix_value_semantics,
                ani_kind=ani_kind,
                sample_count=sample_count,
                cache_all_rows=plan.cache_all_rows,
                # The planner leaves 30% headroom; use at most 20% for raw
                # count reads and conversion, normally in <=128 MiB chunks.
                read_memory_limit_bytes=max(1, memory_limit_bytes // 5),
            )
            try:
                informative, multi = _informative_positions(rows, plan.block_samples)
                plan = _plan_gemm_blocks(
                    sample_count=sample_count,
                    matrix_length=spec.matrix_length,
                    memory_limit_bytes=memory_limit_bytes,
                    device_memory_bytes=device_memory,
                    informative_count=int(informative.size),
                )
                block = plan.block_samples
                starts = list(range(0, sample_count, block))
                sample_blocks_total += len(starts)
                done_buckets = _completed_by_block_pair(completed, row_of_sample, block)
                for bi, start_i in enumerate(starts):
                    stop_i = min(start_i + block, sample_count)
                    host_i = rows.block(start_i, stop_i)
                    inf_i = np.ascontiguousarray(host_i[:, informative])
                    for bj in range(bi, len(starts)):
                        start_j = starts[bj]
                        stop_j = min(start_j + block, sample_count)
                        pending = np.ones((stop_i - start_i, stop_j - start_j), dtype=bool)
                        if bi == bj:
                            pending = np.triu(pending, k=1)
                        for row_a, row_b in done_buckets.get((bi, bj), ()):
                            pending[row_a - start_i, row_b - start_j] = False
                        if not pending.any():
                            continue
                        if bi == bj:
                            host_j, inf_j = host_i, inf_i
                        else:
                            host_j = rows.block(start_j, stop_j)
                            inf_j = np.ascontiguousarray(host_j[:, informative])
                        totals, shared = _block_pair_counts(
                            ops,
                            block_i=host_i,
                            block_j=host_j,
                            inf_i=inf_i,
                            inf_j=inf_j,
                            multi_positions=multi,
                            matrix_length=int(spec.matrix_length),
                            position_chunk=plan.position_chunk,
                        )
                        local_i, local_j = np.nonzero(pending)
                        global_i = local_i + start_i
                        global_j = local_j + start_j
                        pair_totals = totals[local_i, local_j]
                        pair_shared = shared[local_i, local_j]
                        ids_i = compare_sample_ids[global_i]
                        ids_j = compare_sample_ids[global_j]
                        keep = pair_totals > 0
                        table = _result_table(
                            compare_id_1=ids_i[keep],
                            compare_id_2=ids_j[keep],
                            names_1=sample_names[global_i[keep]],
                            names_2=sample_names[global_j[keep]],
                            genome_id=compare_genome_id,
                            genome=spec.genome,
                            totals=pair_totals[keep],
                            shared=pair_shared[keep],
                        )
                        compare_conn.commit()
                        compare_conn.execute("BEGIN")
                        try:
                            mp._insert_matrix_compare_result_table(compare_conn, table)
                            _insert_completed_arrays(
                                compare_conn,
                                np.minimum(ids_i, ids_j),
                                np.maximum(ids_i, ids_j),
                                compare_genome_id,
                            )
                            compare_conn.execute("COMMIT")
                        except Exception:
                            compare_conn.execute("ROLLBACK")
                            raise
                        batch_pairs = int(local_i.size)
                        completed_work += batch_pairs
                        written_rows += table.num_rows
                        block_pairs_done += 1
                        anchor_name = str(sample_names[start_i])
                        if emit_writer_logs:
                            mp._emit_matrix_compare_writer_log(
                                start_time=run_start_time,
                                completed=completed_work,
                                total=total_work,
                                batch_pairs=batch_pairs,
                                batch_rows=table.num_rows,
                                anchor_name=anchor_name,
                                genome=spec.genome,
                                targets_completed=stop_j,
                                targets_total=sample_count,
                                target_chunks=block_pairs_done,
                            )
                        if progress_callback is not None:
                            progress_callback(
                                {
                                    "phase": "advance",
                                    "completed": completed_work,
                                    "total": total_work,
                                    "anchor_name": anchor_name,
                                    "genome": spec.genome,
                                    "scaffold": "",
                                    "targets_completed": stop_j,
                                    "targets_total": sample_count,
                                    "target_chunks": block_pairs_done,
                                }
                            )
            finally:
                rows.close()

    return mp.MatrixCompareSummary(
        output_file=output_file,
        requested_pairs=requested_pairs,
        written_rows=written_rows,
        scaffold_count=len(genome_scaffolds),
        genome_count=len(genomes),
        anchor_groups=sample_blocks_total,
        target_chunks=block_pairs_done,
        min_cov=min_cov,
        genome_scope=genome_scope or metadata.get("genome_scope", "all"),
        backend=compute_backend.kind,
        device=compute_backend.device,
        memory_limit_gb=memory_limit_gb,
        ani_method=ani_method,
    )
