"""Comparison identities must survive rebuilding or reordering HDF5 matrices."""

from itertools import combinations

import duckdb
import numpy as np
import polars as pl
import pytest

from zipstrain import matrix_pairs as mp


def _write_inputs(root):
    profiles = root / "profiles"
    profiles.mkdir()
    bed, stb, genes = root / "reference.bed", root / "reference.stb", root / "genes.tsv"
    bed.write_text("chr1\t0\t6\nchr2\t0\t4\nchr3\t0\t4\n")
    stb.write_text("chr1\tgenome1\nchr2\tgenome1\nchr3\tgenome2\n")
    genes.write_text("gene1\tchr1\t1\t6\ngene2\tchr2\t1\t4\ngene3\tchr3\t1\t4\n")
    frames = {}
    for index, name in enumerate(("sample_a", "sample_b", "sample_c", "sample_d")):
        counts = np.random.default_rng(19 + index).integers(0, 12, (14, 4))
        counts[counts < 8] = 0
        counts[0] = [10 + index, 0, 0, 0]
        # One gene has no overlap for B, while its other genes still overlap.
        if name == "sample_b":
            counts[6:10] = 0
        frames[name] = pl.DataFrame({
            "chrom": ["chr1"] * 6 + ["chr2"] * 4 + ["chr3"] * 4,
            "pos": list(range(1, 7)) + list(range(1, 5)) * 2,
            "genome": ["genome1"] * 10 + ["genome2"] * 4,
            "gene": ["gene1"] * 6 + ["gene2"] * 4 + ["gene3"] * 4,
            **{base: counts[:, column] for column, base in enumerate("ATCG")},
        })
    return profiles, (bed, stb, genes), frames


def _build(matrix, profiles, references, frames, names, storage="bitmask", sparse=False):
    for path in profiles.glob("*.parquet"):
        path.unlink()
    for name in names:
        frames[name].write_parquet(profiles / f"{name}.parquet")
    matrix.unlink(missing_ok=True)
    bed, stb, genes = references
    mp.build_matrix_hdf5(
        profile_dir=profiles, output_file=matrix, bed_file=bed, stb_file=stb,
        gene_range_table=genes, storage_mode=storage, sparse=sparse,
        memory_limit_gb=1,
    )


def _reorder_matrix(matrix, order, *, swap_genomes=False):
    h5py = pytest.importorskip("h5py")
    with h5py.File(matrix, "r+") as handle:
        for key, dataset in handle["samples"].items():
            if key != "sample_idx":
                dataset[...] = dataset[...][order]
        for node in handle["matrices"].values():
            if hasattr(node, "shape"):
                node[...] = node[...][order]
            else:
                pointers = node["indptr"][...]
                lengths = np.diff(pointers)[order]
                for key in ("indices", "values"):
                    if key in node:
                        old = node[key][...]
                        node[key][...] = np.concatenate([
                            old[pointers[row]:pointers[row + 1]] for row in order
                        ])
                node["indptr"][...] = np.r_[0, np.cumsum(lengths)]
        if swap_genomes:
            for group in ("genomes", "genome_scaffolds", "genes", "contract_genomes", "contract_genome_scaffolds"):
                if group in handle and "genome_idx" in handle[group]:
                    dataset = handle[group]["genome_idx"]
                    dataset[...] = 1 - dataset[...]
            for key in ("0", "1"):
                handle["matrices"].move(key, f"temporary_{key}")
            handle["matrices"].move("temporary_0", "1")
            handle["matrices"].move("temporary_1", "0")


def _compare(matrix, output, backend, *, ani_method="popani"):
    if backend.startswith("torch"):
        pytest.importorskip("torch")
    return mp.matrix_compare(
        matrix_db_file=matrix, output_file=output, backend=backend,
        ani_method=ani_method, calculate="all", memory_limit_gb=1,
        anchor_queue_size=2, target_queue_size=2, result_transfer_batch_size=2,
    )


def _snapshot(output):
    with duckdb.connect(str(output), read_only=True) as conn:
        samples = conn.execute("SELECT * FROM matrix_compare_samples ORDER BY sample_idx").fetchall()
        genomes = conn.execute("SELECT * FROM matrix_compare_genomes ORDER BY genome_idx").fetchall()
        completed = conn.execute(
            "SELECT * FROM matrix_compare_completed_pair_genomes ORDER BY ALL"
        ).fetchall()
        results = conn.execute("SELECT * FROM matrix_compare_results ORDER BY ALL").pl()
        genes = conn.execute("SELECT * FROM matrix_compare_gene_results ORDER BY ALL").pl()
    return samples, genomes, completed, results, genes


def _normalized(frame):
    first, second = pl.col("sample_1"), pl.col("sample_2")
    frame = frame.with_columns(
        sample_1=pl.when(first < second).then(first).otherwise(second),
        sample_2=pl.when(first < second).then(second).otherwise(first),
    ).drop("sample_idx_1", "sample_idx_2", "genome_idx")
    keys = ["sample_1", "sample_2", "genome"]
    if "gene" in frame.columns:
        keys.append("gene")
    return frame.sort(keys)


def _assert_results_match(resumed, fresh):
    _, _, completed, results, genes = _snapshot(resumed)
    _, _, expected_completed, expected, expected_genes = _snapshot(fresh)
    assert len(completed) == len(expected_completed)
    assert len(completed) == len(set(completed))
    assert all(first < second for first, second, _genome in completed)
    assert _normalized(results).equals(_normalized(expected))
    assert _normalized(genes).equals(_normalized(expected_genes))
    with duckdb.connect(str(resumed), read_only=True) as conn:
        for table in ("matrix_compare_results", "matrix_compare_gene_results"):
            assert conn.execute(f"""
                SELECT count(*) FROM {table} r
                LEFT JOIN matrix_compare_samples a ON r.sample_idx_1 = a.sample_idx
                LEFT JOIN matrix_compare_samples b ON r.sample_idx_2 = b.sample_idx
                LEFT JOIN matrix_compare_genomes g ON r.genome_idx = g.genome_idx
                WHERE a.sample_name IS DISTINCT FROM r.sample_1
                   OR b.sample_name IS DISTINCT FROM r.sample_2
                   OR g.genome IS DISTINCT FROM r.genome
                   OR r.sample_idx_1 >= r.sample_idx_2
            """).fetchone()[0] == 0


@pytest.mark.parametrize("storage,sparse,ani_method", [
    ("bitmask", False, "popani"), ("bitmask", True, "popani"),
    ("counts", False, "popani"), ("counts", True, "popani"),
    ("counts", False, "conani"), ("counts", True, "cosani_0.95"),
])
@pytest.mark.parametrize("backend", ["numpy", "torch-cpu"])
def test_rebuilt_matrix_reuses_existing_comparisons(tmp_path, monkeypatch, storage, sparse, ani_method, backend):
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output, fresh = tmp_path / "matrix.h5", tmp_path / "compare.duckdb", tmp_path / "fresh.duckdb"
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c"], storage, sparse)
    _compare(matrix, output, backend, ani_method=ani_method)
    old_samples, old_genomes, old_completed, old_results, old_genes = _snapshot(output)

    # A now precedes B/C physically. Also exercise nontrivial sample/genome orders.
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c", "sample_d"], storage, sparse)
    _reorder_matrix(matrix, [3, 0, 2, 1], swap_genomes=True)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", lambda vector_length, remaining_targets, **kwargs: (min(2, remaining_targets), vector_length))
    summary = _compare(matrix, output, backend, ani_method=ani_method)
    assert summary.requested_pairs == 5
    samples, genomes, completed, results, genes = _snapshot(output)
    assert samples[:2] == old_samples
    assert genomes == old_genomes
    assert {name for _index, name in samples[2:]} == {"sample_a", "sample_d"}
    assert set(old_completed) <= set(completed)
    assert results.filter((pl.col("sample_idx_1") < 2) & (pl.col("sample_idx_2") < 2)).equals(old_results)
    assert genes.filter((pl.col("sample_idx_1") < 2) & (pl.col("sample_idx_2") < 2)).equals(old_genes)
    _compare(matrix, fresh, "numpy", ani_method=ani_method)
    _assert_results_match(output, fresh)
    assert _compare(matrix, output, backend, ani_method=ani_method).requested_pairs == 0
    assert _snapshot(output)[2] == completed


@pytest.mark.parametrize("backend", ["numpy", "torch-cpu"])
@pytest.mark.parametrize("ani_method", ["popani", "conani", "cosani_0.95"])
def test_reordered_completed_matrix_does_no_computation(tmp_path, monkeypatch, backend, ani_method):
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output = tmp_path / "matrix.h5", tmp_path / "compare.duckdb"
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c"], "counts", True)
    _compare(matrix, output, backend, ani_method=ani_method)
    before = _snapshot(output)
    _reorder_matrix(matrix, [2, 0, 1], swap_genomes=True)

    def unexpected(*args, **kwargs):
        pytest.fail("A complete reordered matrix must not load or compare matrices.")

    monkeypatch.setattr(mp._Hdf5GenomeMatrixNumpyDataset, "get_row", unexpected)
    monkeypatch.setattr(mp, "_matrix_compare_reuse_target_chunks_torch", unexpected)
    summary = _compare(matrix, output, backend, ani_method=ani_method)
    assert summary.requested_pairs == summary.written_rows == 0
    after = _snapshot(output)
    assert after[:3] == before[:3]
    assert after[3].equals(before[3]) and after[4].equals(before[4])


@pytest.mark.parametrize("backend", ["numpy", "torch-cpu"])
def test_subset_rebuild_preserves_ids_and_reintroduced_pairs(tmp_path, backend):
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output, fresh = tmp_path / "matrix.h5", tmp_path / "compare.duckdb", tmp_path / "fresh.duckdb"
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c"])
    _compare(matrix, output, backend)
    catalog = _snapshot(output)[0]
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c"])
    assert _compare(matrix, output, backend).requested_pairs == 0
    assert _snapshot(output)[0] == catalog
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c", "sample_d"])
    assert _compare(matrix, output, backend).requested_pairs == 2
    assert _snapshot(output)[0] == [*catalog, (3, "sample_d")]
    _build(matrix, profiles, references, frames, list(frames))
    assert _compare(matrix, output, backend).requested_pairs == 1
    _compare(matrix, fresh, "numpy")
    _assert_results_match(output, fresh)


@pytest.mark.parametrize("backend", ["numpy", "torch-cpu"])
def test_zero_overlap_completion_survives_rebuild(tmp_path, backend):
    profiles, references, frames = _write_inputs(tmp_path)
    frames["sample_b"] = frames["sample_b"].filter(pl.col("genome") == "genome1")
    frames["sample_c"] = frames["sample_c"].filter(pl.col("genome") == "genome2")
    matrix, output, fresh = tmp_path / "matrix.h5", tmp_path / "compare.duckdb", tmp_path / "fresh.duckdb"
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c"], sparse=True)
    _compare(matrix, output, backend)
    assert _snapshot(output)[3].is_empty()
    assert len(_snapshot(output)[2]) == 2
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c"], sparse=True)
    assert _compare(matrix, output, backend).requested_pairs == 2
    _compare(matrix, fresh, "numpy")
    _assert_results_match(output, fresh)


def test_identity_catalog_rejects_ambiguous_existing_names():
    with duckdb.connect() as conn:
        mp._init_matrix_compare_db_schema(conn)
        conn.execute("INSERT INTO matrix_compare_samples VALUES (0, 'B'), (1, 'B')")
        with pytest.raises(ValueError, match="invalid identity catalog"):
            mp._sync_matrix_compare_catalog(
                conn, table="matrix_compare_samples", index_column="sample_idx",
                name_column="sample_name", rows=[(0, "B"), (1, "C")],
            )
        assert conn.execute("SELECT * FROM matrix_compare_samples ORDER BY sample_idx").fetchall() == [(0, "B"), (1, "B")]


@pytest.mark.parametrize("backend,resume_backend", [("numpy", "torch-cpu"), ("torch-cpu", "numpy")])
def test_interrupted_comparison_resumes_after_reorder(tmp_path, monkeypatch, backend, resume_backend):
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output, fresh = tmp_path / "matrix.h5", tmp_path / "compare.duckdb", tmp_path / "fresh.duckdb"
    _build(matrix, profiles, references, frames, list(frames), sparse=True)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", lambda vector_length, remaining_targets, **kwargs: (min(2, remaining_targets), vector_length))
    monkeypatch.setattr(mp, "MATRIX_COMPARE_CHECKPOINT_BATCH_UNITS", 1)
    monkeypatch.setattr(mp, "MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS", 1)
    monkeypatch.setattr(mp, "MATRIX_COMPARE_WRITE_MAX_ROWS", 1)
    original_mark = mp._mark_completed_pair_genomes
    calls = 0

    def interrupt(conn, completed):
        nonlocal calls
        calls += 1
        original_mark(conn, completed)
        if calls == 2:
            raise RuntimeError("simulated interruption")

    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", interrupt)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        _compare(matrix, output, backend)
    before = _snapshot(output)
    assert 0 < len(before[2]) < 12
    _reorder_matrix(matrix, [3, 1, 0, 2], swap_genomes=True)
    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", original_mark)
    _compare(matrix, output, resume_backend)
    after = _snapshot(output)
    assert after[:2] == before[:2]
    assert set(before[2]) <= set(after[2])
    assert len(after[2]) == len(list(combinations(frames, 2))) * 2
    _compare(matrix, fresh, "numpy")
    _assert_results_match(output, fresh)


def test_rebuilt_matrix_resume_with_process_executors(tmp_path):
    pytest.importorskip("torch")
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output, fresh = tmp_path / "matrix.h5", tmp_path / "compare.duckdb", tmp_path / "fresh.duckdb"
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c"], sparse=True)
    _compare(matrix, output, "numpy")
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c"], sparse=True)
    summary = mp.matrix_compare(
        matrix, output, backend="torch-cpu", calculate="all", memory_limit_gb=1,
        loader_executor_kind="process", writer_executor_kind="process", target_queue_size=2,
    )
    assert summary.requested_pairs == 2
    _compare(matrix, fresh, "numpy")
    _assert_results_match(output, fresh)


def test_identity_reconciliation_rolls_back_on_invalid_genome_catalog(tmp_path):
    profiles, references, frames = _write_inputs(tmp_path)
    matrix, output = tmp_path / "matrix.h5", tmp_path / "compare.duckdb"
    _build(matrix, profiles, references, frames, ["sample_b", "sample_c"])
    _compare(matrix, output, "numpy")
    with duckdb.connect(str(output)) as conn:
        conn.execute("UPDATE matrix_compare_genomes SET genome = 'genome1'")
    _build(matrix, profiles, references, frames, ["sample_a", "sample_b", "sample_c"])
    with pytest.raises(ValueError, match="invalid identity catalog"):
        _compare(matrix, output, "numpy")
    assert _snapshot(output)[0] == [(0, "sample_b"), (1, "sample_c")]


def test_pair_columns_canonicalize_mixed_target_ids():
    columns = mp._compare_pair_arrow_columns(2, "anchor", [4, 0, 3, 1], ["four", "zero", "three", "one"])
    assert [column.to_pylist() for column in columns] == [
        [2, 0, 2, 1], [4, 2, 3, 2],
        ["anchor", "zero", "anchor", "one"], ["four", "anchor", "three", "anchor"],
    ]


def test_resume_remaps_existing_keys_without_scanning_result_tables():
    class CompletionOnlyConnection:
        def execute(self, sql, parameters):
            assert "matrix_compare_completed_pair_genomes" in sql
            assert "matrix_compare_results" not in sql and "JOIN" not in sql
            assert parameters == [9, 9]
            return self

        def fetchall(self):
            return [(4, 8, 7), (1, 4, 7), (4, 8, 99)]

    identities = mp._MatrixCompareIds(samples=np.array([8, 4]), genomes={0: 7})
    completed, pairs, work = mp._load_matrix_compare_resume_state(
        CompletionOnlyConnection(), sample_count=2, genome_ids=[0], compare_ids=identities,
    )
    assert completed == {0: {(0, 1)}}
    assert pairs == work == 0
