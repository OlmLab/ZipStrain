"""Parity tests: compare_kernel='block' must match the anchor numpy kernel exactly."""

from pathlib import Path

from click.testing import CliRunner
import duckdb
import numpy as np
import polars as pl
import pytest

from zipstrain import cli
from zipstrain import matrix_gemm as mg
from zipstrain import matrix_pairs as mp

SCAFFOLDS = (("chr1", "genome1", 240), ("chr2", "genome1", 130), ("chr3", "genome2", 170))


def _write_random_profiles(profile_dir: Path, sample_count: int = 23, seed: int = 0) -> None:
    """Random multi-strain profiles with mixed alleles, ties and coverage gaps.

    Samples are drawn from three strains that share a backbone and differ at
    ~15% of sites, so most positions are non-informative and many pairs share most alleles
    (near-clonal) while others differ. Some positions carry two alleles, a
    fraction of those as exact count ties (exercises conANI consensus ties),
    and one sample does not cover genome2, so its genome2 pairs have
    total_positions == 0 (no result row, but still marked completed).
    """
    rng = np.random.default_rng(seed)
    profile_dir.mkdir(parents=True, exist_ok=True)
    strains = {}
    for chrom, _genome, length in SCAFFOLDS:
        backbone = rng.integers(0, 4, size=length)
        variants = np.tile(backbone, (3, 1))
        snps = rng.random((3, length)) < 0.15
        variants[snps] = rng.integers(0, 4, size=int(snps.sum()))
        strains[chrom] = variants
    for sample in range(sample_count):
        strain = sample % 3
        coverage = 0.25 + 0.7 * rng.random()
        frames = []
        for chrom, genome, length in SCAFFOLDS:
            pos = np.arange(1, length + 1)
            keep = rng.random(length) < (0.0 if (sample == 5 and genome == "genome2") else coverage)
            if sample == 0:
                keep[[0, length - 1]] = True  # make the contract span the scaffold
            counts = np.zeros((length, 4), dtype=np.int64)
            major = strains[chrom][strain].copy()
            mutate = rng.random(length) < 0.01
            major[mutate] = rng.integers(0, 4, size=int(mutate.sum()))
            counts[np.arange(length), major] = rng.integers(6, 30, size=length)
            minor_pos = np.flatnonzero(rng.random(length) < 0.04)
            minor = rng.integers(0, 4, size=minor_pos.size)
            tie = rng.random(minor_pos.size) < 0.3
            counts[minor_pos, minor] = np.where(
                tie, counts[minor_pos, major[minor_pos]], rng.integers(1, 6, size=minor_pos.size)
            )
            sel = np.flatnonzero(keep)
            frames.append(
                pl.DataFrame(
                    {
                        "chrom": [chrom] * sel.size,
                        "pos": pos[sel].astype(np.int64),
                        "gene": [f"{chrom}_g{p // 60}" for p in pos[sel]],
                        "genome": [genome] * sel.size,
                        "A": counts[sel, 0],
                        "T": counts[sel, 1],
                        "C": counts[sel, 2],
                        "G": counts[sel, 3],
                    }
                )
            )
        pl.concat(frames).write_parquet(profile_dir / f"sample_{sample:02d}.parquet")


def _build_store(tmp_path: Path, *, storage_mode: str, sparse: bool = False, sample_count: int = 23) -> Path:
    profiles = tmp_path / "profiles"
    _write_random_profiles(profiles, sample_count=sample_count)
    frame = (
        pl.scan_parquet(str(profiles / "*.parquet"))
        .group_by(["genome", "chrom"])
        .agg(pl.col("pos").min().alias("lo"), pl.col("pos").max().alias("hi"))
        .collect()
        .sort("chrom")
    )
    bed, stb = tmp_path / "ref.bed", tmp_path / "ref.stb"
    bed.write_text("".join(f"{r['chrom']}\t{r['lo'] - 1}\t{r['hi']}\n" for r in frame.iter_rows(named=True)))
    stb.write_text("".join(f"{r['chrom']}\t{r['genome']}\n" for r in frame.iter_rows(named=True)))
    matrix = tmp_path / f"matrix_{storage_mode}{'_sparse' if sparse else ''}.h5"
    mp.build_matrix_hdf5(
        profile_dir=profiles,
        output_file=matrix,
        bed_file=bed,
        stb_file=stb,
        storage_mode=storage_mode,
        sparse=sparse,
        memory_limit_gb=1,
    )
    return matrix


def _load(compare_db: Path):
    conn = duckdb.connect(str(compare_db), read_only=True)
    try:
        results = pl.from_arrow(
            conn.execute(
                """
                SELECT sample_idx_1, sample_idx_2, sample_1, sample_2, genome_idx, genome,
                       total_positions, share_allele_pos, genome_ani, max_consecutive_length
                FROM matrix_compare_results
                ORDER BY sample_idx_1, sample_idx_2, genome_idx
                """
            ).fetch_arrow_table()
        )
        completed = conn.execute(
            "SELECT sample_idx_1, sample_idx_2, genome_idx FROM matrix_compare_completed_pair_genomes "
            "ORDER BY ALL"
        ).fetchall()
    finally:
        conn.close()
    return results, completed


def _small_plan(block: int = 4, chunk: int = 7, cache_all: bool = False):
    def plan(**_kwargs):
        return mg.GemmPlan(block_samples=block, position_chunk=chunk, cache_all_rows=cache_all)

    return plan


def _assert_same(actual_db: Path, expected_db: Path) -> None:
    actual, actual_completed = _load(actual_db)
    expected, expected_completed = _load(expected_db)
    assert actual_completed == expected_completed
    assert actual.height == expected.height > 0
    assert actual.equals(expected)


def _reference(matrix: Path, tmp_path: Path, ani_method: str = "popani") -> Path:
    out = tmp_path / f"pairwise_{ani_method}.duckdb"
    mp.matrix_compare(matrix, out, backend="numpy", ani_method=ani_method, calculate="ani", memory_limit_gb=1)
    return out


@pytest.mark.parametrize(
    "storage_mode,ani_method,sparse",
    [
        ("bitmask", "popani", False),
        ("bitmask", "popani", True),
        ("counts", "popani", False),
        ("counts", "conani", False),
        ("counts", "conani", True),
    ],
)
@pytest.mark.parametrize("backend", ["numpy", "torch-cpu"])
@pytest.mark.parametrize("plan", [(4, 7, False), (23, 64, True), (5, 1000, True)])
def test_gemm_matches_pairwise(tmp_path, monkeypatch, storage_mode, ani_method, sparse, backend, plan):
    if backend.startswith("torch"):
        pytest.importorskip("torch")
    matrix = _build_store(tmp_path, storage_mode=storage_mode, sparse=sparse)
    expected = _reference(matrix, tmp_path, ani_method)
    monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(*plan))
    out = tmp_path / "gemm.duckdb"
    summary = mp.matrix_compare(
        matrix, out, backend=backend, ani_method=ani_method, calculate="ani",
        memory_limit_gb=1, compare_kernel="block",
    )
    _assert_same(out, expected)
    assert summary.requested_pairs == 23 * 22 // 2
    assert summary.written_rows == _load(expected)[0].height


@pytest.mark.parametrize("backend", ["torch-mps", "torch-cuda"])
def test_gemm_matches_pairwise_on_accelerators(tmp_path, monkeypatch, backend):
    torch = pytest.importorskip("torch")
    if backend == "torch-mps" and not torch.backends.mps.is_available():
        pytest.skip("MPS not available")
    if backend == "torch-cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    matrix = _build_store(tmp_path, storage_mode="counts")
    for ani_method in ("popani", "conani"):
        expected = _reference(matrix, tmp_path, ani_method)
        monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(6, 50))
        out = tmp_path / f"gemm_{ani_method}.duckdb"
        mp.matrix_compare(matrix, out, backend=backend, ani_method=ani_method, calculate="ani",
                          memory_limit_gb=1, compare_kernel="block")
        _assert_same(out, expected)


def test_gemm_default_planner_matches_pairwise(tmp_path):
    matrix = _build_store(tmp_path, storage_mode="bitmask")
    expected = _reference(matrix, tmp_path)
    out = tmp_path / "gemm.duckdb"
    mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel="block")
    _assert_same(out, expected)


def test_gemm_exact_accumulator_flushes(tmp_path, monkeypatch):
    """A tiny exactness bound forces many float32 -> int64 flushes."""
    matrix = _build_store(tmp_path, storage_mode="counts")
    expected = _reference(matrix, tmp_path)
    monkeypatch.setattr(mg, "GEMM_EXACT_FP32_LIMIT", 64)
    monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(7, 1000))
    flushes = 0
    original = mg._ExactAccumulator.flush

    def counting_flush(self):
        nonlocal flushes
        flushes += int(self.pending_positions > 0)
        original(self)

    monkeypatch.setattr(mg._ExactAccumulator, "flush", counting_flush)
    out = tmp_path / "gemm.duckdb"
    mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel="block")
    _assert_same(out, expected)
    assert flushes > 100


class _Interrupt(RuntimeError):
    pass


def _interrupt_after(count: int):
    seen = 0

    def callback(event):
        nonlocal seen
        if event.get("phase") == "advance":
            seen += 1
            if seen >= count:
                raise _Interrupt()

    return callback


def test_gemm_resumes_after_interruption(tmp_path, monkeypatch):
    matrix = _build_store(tmp_path, storage_mode="bitmask")
    expected = _reference(matrix, tmp_path)
    monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(4, 50))
    out = tmp_path / "gemm.duckdb"
    with pytest.raises(_Interrupt):
        mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel="block",
                          progress_callback=_interrupt_after(9))
    partial, partial_completed = _load(out)
    assert 0 < len(partial_completed) < len(_load(expected)[1])
    summary = mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel="block")
    assert summary.written_rows + partial.height == _load(expected)[0].height
    _assert_same(out, expected)
    # A third run has nothing left to do.
    assert mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1,
                             compare_kernel="block").requested_pairs == 0


def test_gemm_resumes_partial_pairwise_run(tmp_path, monkeypatch):
    """Kernels share the completion table, so either can finish the other's run."""
    matrix = _build_store(tmp_path, storage_mode="counts")
    expected = _reference(matrix, tmp_path, "conani")
    out = tmp_path / "mixed.duckdb"
    with pytest.raises(_Interrupt):
        mp.matrix_compare(matrix, out, ani_method="conani", calculate="ani", memory_limit_gb=1,
                          progress_callback=_interrupt_after(7))
    assert 0 < len(_load(out)[1])
    monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(5, 30))
    mp.matrix_compare(matrix, out, ani_method="conani", calculate="ani", memory_limit_gb=1, compare_kernel="block")
    _assert_same(out, expected)


def test_gemm_genome_scope(tmp_path):
    matrix = _build_store(tmp_path, storage_mode="bitmask")
    expected = tmp_path / "pairwise_g2.duckdb"
    mp.matrix_compare(matrix, expected, genome="genome2", calculate="ani", memory_limit_gb=1)
    out = tmp_path / "gemm_g2.duckdb"
    mp.matrix_compare(matrix, out, genome="genome2", calculate="ani", memory_limit_gb=1, compare_kernel="block")
    _assert_same(out, expected)
    assert set(_load(out)[0]["genome"].to_list()) == {"genome2"}


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"calculate": "all"}, "calculate='ani'"),
        ({"calculate": "ani+ibs"}, "calculate='ani'"),
        ({"calculate": "ani", "ani_method": "cosani_0.4"}, "popani or conani"),
        ({"calculate": "ani", "compare_kernel": "fast"}, "compare_kernel must be one of"),
    ],
)
def test_gemm_rejects_unsupported_requests(tmp_path, kwargs, match):
    matrix = _build_store(tmp_path, storage_mode="counts", sample_count=4)
    kwargs = {"compare_kernel": "block", **kwargs}
    out = tmp_path / "bad.duckdb"
    with pytest.raises(ValueError, match=match):
        mp.matrix_compare(matrix, out, memory_limit_gb=1, **kwargs)


def test_gemm_cli(tmp_path):
    matrix = _build_store(tmp_path, storage_mode="bitmask")
    expected = _reference(matrix, tmp_path)
    out = tmp_path / "cli.duckdb"
    result = CliRunner().invoke(
        cli.cli,
        [
            "utilities", "matrix-compare",
            "--matrix-db-file", str(matrix),
            "--output-file", str(out),
            "--memory-limit-gb", "1",
            "--calculate", "ani",
            "--compare-kernel", "block",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "requested_pairs=253" in result.output
    _assert_same(out, expected)


def test_rows_to_allele_bitmask_semantics():
    rng = np.random.default_rng(3)
    counts = rng.integers(0, 3, size=(6, 40, 4)).astype(np.uint16)
    counts[0, 0] = [0, 0, 0, 0]
    counts[0, 1] = [2, 0, 2, 1]
    pop = mg.rows_to_allele_bitmask(counts, matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind="popani")
    con = mg.rows_to_allele_bitmask(counts, matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind="conani")
    presence = mg.rows_to_allele_bitmask(
        (counts > 0).astype(np.uint8),
        matrix_value_semantics=mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS,
        ani_kind="popani",
    )
    bits = np.array(mp.BITMASK_BASE_BITS)
    assert np.array_equal(pop, ((counts > 0) * bits).sum(axis=2))
    assert np.array_equal(presence, pop)
    maxima = counts.max(axis=2, keepdims=True)
    assert np.array_equal(con, (((counts == maxima) & (maxima > 0)) * bits).sum(axis=2))
    assert pop[0, 0] == 0 and con[0, 0] == 0
    assert con[0, 1] == bits[0] | bits[2]
    with pytest.raises(ValueError):
        mg.rows_to_allele_bitmask(pop, matrix_value_semantics=mp.BITMASK_MATRIX_VALUE_SEMANTICS, ani_kind="conani")


def test_planner_respects_memory_budget():
    small = mg._plan_gemm_blocks(sample_count=10_000, matrix_length=3_000_000, memory_limit_bytes=8 * 1024**3)
    large = mg._plan_gemm_blocks(sample_count=10_000, matrix_length=3_000_000, memory_limit_bytes=128 * 1024**3)
    assert 1 <= small.block_samples < large.block_samples <= mg.GEMM_MAX_BLOCK_SAMPLES
    assert 2 * small.block_samples * 3_000_000 < 8 * 1024**3
    tiny = mg._plan_gemm_blocks(sample_count=50, matrix_length=1000, memory_limit_bytes=1024**3)
    assert tiny.block_samples == 50 and tiny.cache_all_rows


@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("cache_all", [False, True])
@pytest.mark.parametrize("ani_method", ["popani", "conani"])
def test_count_mask_reads_are_bounded(tmp_path, monkeypatch, sparse, cache_all, ani_method):
    """The mask cache must not first materialize the entire raw count matrix."""
    h5py = pytest.importorskip("h5py")
    matrix = _build_store(tmp_path, storage_mode="counts", sparse=sparse, sample_count=7)
    metadata = mp._load_matrix_hdf5_metadata(matrix)
    spec = mp._load_matrix_hdf5_genomes(matrix)[0]
    dtype_name = mp._matrix_dtype_name_from_metadata(metadata)
    dataset = mp._Hdf5GenomeMatrixNumpyDataset(
        matrix, spec.genome_idx, spec.matrix_length, dtype_name,
        mp.COUNT_MATRIX_VALUE_SEMANTICS,
    )
    try:
        raw = dataset._load_batch_numpy(np.arange(7))
    finally:
        dataset.close()
    expected = mg.rows_to_allele_bitmask(
        raw, matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind=ani_method,
    )
    # Allow two rows of counts plus conversion workspace, but not all seven.
    bytes_per_position = (8 * raw.dtype.itemsize + 32) if sparse else (5 * raw.dtype.itemsize + 16)
    read_limit = 2 * spec.matrix_length * bytes_per_position
    reads = []
    if sparse:
        original_load = mp._Hdf5GenomeMatrixNumpyDataset._load_batch_numpy

        def tracked_load(self, sample_rows):
            reads.append(len(sample_rows))
            assert len(sample_rows) <= 2
            return original_load(self, sample_rows)

        monkeypatch.setattr(mp._Hdf5GenomeMatrixNumpyDataset, "_load_batch_numpy", tracked_load)
    else:
        original_getitem = h5py.Dataset.__getitem__

        def tracked_getitem(self, key):
            if self.ndim == 3 and self.name.startswith("/matrices/"):
                reads.append(key[0].stop - key[0].start)
                assert reads[-1] <= 2
            return original_getitem(self, key)

        monkeypatch.setattr(h5py.Dataset, "__getitem__", tracked_getitem)
    rows = mg._GenomeRows(
        matrix_db_file=matrix, spec=spec, dtype_name=dtype_name,
        matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind=ani_method,
        sample_count=7, cache_all_rows=cache_all, read_memory_limit_bytes=read_limit,
    )
    try:
        assert np.array_equal(rows.block(0, 7), expected)
        assert reads == [2, 2, 2, 1]
        assert np.array_equal(rows.block(1, 4), expected[1:4])
    finally:
        rows.close()


def test_dense_count_mask_reads_tile_long_rows(tmp_path, monkeypatch):
    h5py = pytest.importorskip("h5py")
    matrix = _build_store(tmp_path, storage_mode="counts", sample_count=4)
    spec = mp._load_matrix_hdf5_genomes(matrix)[0]
    dataset = mp._Hdf5GenomeMatrixNumpyDataset(
        matrix, spec.genome_idx, spec.matrix_length, "uint16", mp.COUNT_MATRIX_VALUE_SEMANTICS,
    )
    try:
        raw = dataset._load_batch_numpy(np.arange(4))
    finally:
        dataset.close()
    expected = mg.rows_to_allele_bitmask(
        raw, matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind="conani",
    )
    original_getitem = h5py.Dataset.__getitem__
    reads = []

    def tracked_getitem(self, key):
        if self.ndim == 3 and self.name.startswith("/matrices/"):
            cells = (key[0].stop - key[0].start) * (key[1].stop - key[1].start)
            reads.append(cells)
            assert cells <= 17
        return original_getitem(self, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", tracked_getitem)
    rows = mg._GenomeRows(
        matrix_db_file=matrix, spec=spec, dtype_name="uint16",
        matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind="conani",
        sample_count=4, cache_all_rows=True, read_memory_limit_bytes=17 * 26,
    )
    try:
        assert np.array_equal(rows.block(0, 4), expected)
        assert len(reads) > 4
    finally:
        rows.close()


def test_sparse_count_read_fails_before_oversized_allocation(tmp_path, monkeypatch):
    matrix = _build_store(tmp_path, storage_mode="counts", sparse=True, sample_count=4)
    spec = mp._load_matrix_hdf5_genomes(matrix)[0]
    closed = []
    original_close = mp._Hdf5GenomeMatrixNumpyDataset.close

    def tracked_close(self):
        closed.append(True)
        original_close(self)

    def oversized_read(self, sample_rows):
        pytest.fail("Sparse counts were loaded before checking the conversion workspace")

    monkeypatch.setattr(mp._Hdf5GenomeMatrixNumpyDataset, "close", tracked_close)
    monkeypatch.setattr(mp._Hdf5GenomeMatrixNumpyDataset, "_load_batch_numpy", oversized_read)
    with pytest.raises(MemoryError, match="expanded sparse count row"):
        mg._GenomeRows(
            matrix_db_file=matrix, spec=spec, dtype_name="uint16",
            matrix_value_semantics=mp.COUNT_MATRIX_VALUE_SEMANTICS, ani_kind="conani",
            sample_count=4, cache_all_rows=True, read_memory_limit_bytes=17 * 26,
        )
    assert closed


@pytest.mark.parametrize("storage,ani_method", [("bitmask", "popani"), ("counts", "conani")])
def test_gemm_appends_only_new_comparisons(tmp_path, storage, ani_method):
    matrix = _build_store(tmp_path, storage_mode=storage, sample_count=4)
    output = tmp_path / "incremental.duckdb"
    mp.matrix_compare(
        matrix, output, ani_method=ani_method, calculate="ani", memory_limit_gb=1,
        compare_kernel="block",
    )
    before, completed_before = _load(output)
    extra = tmp_path / "extra"
    _write_random_profiles(extra, sample_count=6)
    for index in range(4):
        (extra / f"sample_{index:02d}.parquet").unlink()
    mp.append_matrix_hdf5(extra, matrix, memory_limit_gb=1)
    summary = mp.matrix_compare(
        matrix, output, ani_method=ani_method, calculate="ani", memory_limit_gb=1,
        compare_kernel="block",
    )
    after, completed_after = _load(output)
    assert summary.requested_pairs == 9
    assert set(completed_before) <= set(completed_after)
    assert len(completed_after) == 30
    assert after.filter(pl.col("sample_idx_2") < 4).equals(before)
    assert summary.written_rows == after.height - before.height
    _assert_same(output, _reference(matrix, tmp_path, ani_method))


@pytest.mark.parametrize("first_kernel", ["block", "anchor"])
def test_gemm_resumes_after_matrix_reorder(tmp_path, monkeypatch, first_kernel):
    """Compare IDs differ from matrix rows after a reorder; GEMM must remap and canonicalize."""
    from .test_matrix_compare_identity import _reorder_matrix

    matrix = _build_store(tmp_path, storage_mode="bitmask", sparse=True)
    monkeypatch.setattr(mg, "_plan_gemm_blocks", _small_plan(5, 40))
    out = tmp_path / "compare.duckdb"
    with pytest.raises(_Interrupt):
        mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel=first_kernel,
                          progress_callback=_interrupt_after(6))
    before = _load(out)[1]
    assert 0 < len(before) < 506
    order = list(np.random.default_rng(7).permutation(23))
    _reorder_matrix(matrix, order, swap_genomes=True)
    mp.matrix_compare(matrix, out, calculate="ani", memory_limit_gb=1, compare_kernel="block")
    fresh = tmp_path / "fresh.duckdb"
    mp.matrix_compare(matrix, fresh, calculate="ani", memory_limit_gb=1)
    resumed, completed = _load(out)
    expected, expected_completed = _load(fresh)
    assert set(before) <= set(completed)
    assert len(completed) == len(set(completed)) == len(expected_completed) == 506
    assert all(first < second for first, second, _genome in completed)

    def by_name(frame):
        return frame.drop("sample_idx_1", "sample_idx_2", "genome_idx").with_columns(
            s1=pl.min_horizontal("sample_1", "sample_2"), s2=pl.max_horizontal("sample_1", "sample_2")
        ).drop("sample_1", "sample_2").sort(["s1", "s2", "genome"])

    assert by_name(resumed).equals(by_name(expected))
    with duckdb.connect(str(out), read_only=True) as conn:
        assert conn.execute("""
            SELECT count(*) FROM matrix_compare_results r
            LEFT JOIN matrix_compare_samples a ON r.sample_idx_1 = a.sample_idx
            LEFT JOIN matrix_compare_samples b ON r.sample_idx_2 = b.sample_idx
            LEFT JOIN matrix_compare_genomes g ON r.genome_idx = g.genome_idx
            WHERE a.sample_name IS DISTINCT FROM r.sample_1 OR b.sample_name IS DISTINCT FROM r.sample_2
               OR g.genome IS DISTINCT FROM r.genome OR r.sample_idx_1 >= r.sample_idx_2
        """).fetchone()[0] == 0


@pytest.mark.parametrize(("alias", "canonical"), [("gemm", "block"), ("pairwise", "anchor")])
def test_deprecated_kernel_names_still_work(alias, canonical):
    with pytest.warns(FutureWarning, match=f"use '{canonical}'"):
        assert mp.normalize_matrix_compare_kernel(alias) == canonical


def test_gemm_cli_accepts_deprecated_name(tmp_path):
    matrix = _build_store(tmp_path, storage_mode="bitmask")
    expected = _reference(matrix, tmp_path)
    out = tmp_path / "cli.duckdb"
    result = CliRunner().invoke(
        cli.cli,
        [
            "utilities", "matrix-compare",
            "--matrix-db-file", str(matrix),
            "--output-file", str(out),
            "--memory-limit-gb", "1",
            "--calculate", "ani",
            "--compare-kernel", "gemm",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "deprecated; use --compare-kernel block" in result.output
    _assert_same(out, expected)
