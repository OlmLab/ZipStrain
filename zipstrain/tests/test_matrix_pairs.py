from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner
import duckdb
import numpy as np
import polars as pl
import pytest

from zipstrain import cli
from zipstrain import compare as cp
from zipstrain import matrix_pairs as mp


def _write_profiles(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [0, 2, 5],
            "gene": ["gene1", "gene1", "gene2"],
            "genome": ["genome1", "genome1", "genome2"],
            "A": [6, 0, 0],
            "T": [0, 0, 8],
            "C": [0, 3, 0],
            "G": [0, 2, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr2"],
            "pos": [1, 5, 7],
            "gene": ["gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome2", "genome2"],
            "A": [0, 1, 0],
            "T": [9, 1, 0],
            "C": [0, 0, 0],
            "G": [0, 0, 6],
        }
    )
    sample_c = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 5, 7],
            "gene": ["gene1", "gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome1", "genome2", "genome2"],
            "A": [5, 0, 0, 0],
            "T": [0, 7, 4, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 5],
        }
    )
    ignored = pl.DataFrame({"cov": [1]})

    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")
    sample_c.write_parquet(profile_dir / "sample_c.parquet")
    ignored.write_parquet(profile_dir / "ignored.parquet")


def _write_profiles_one_based(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [1, 3, 6],
            "gene": ["gene1", "gene1", "gene2"],
            "genome": ["genome1", "genome1", "genome2"],
            "A": [6, 0, 0],
            "T": [0, 0, 8],
            "C": [0, 3, 0],
            "G": [0, 2, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr2"],
            "pos": [2, 6, 8],
            "gene": ["gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome2", "genome2"],
            "A": [0, 1, 0],
            "T": [9, 1, 0],
            "C": [0, 0, 0],
            "G": [0, 0, 6],
        }
    )
    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")


def _write_many_profiles_same_genome(profile_dir: Path, sample_count: int = 5) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    for idx in range(sample_count):
        frame = pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "pos": [0, 1, 2],
                "gene": ["gene1", "gene1", "gene1"],
                "genome": ["genome1", "genome1", "genome1"],
                "A": [6 if idx % 2 == 0 else 0, 0, 6],
                "T": [0, 6 if idx % 2 else 0, 0],
                "C": [0, 0, 0],
                "G": [0, 0, 0],
            }
        )
        frame.write_parquet(profile_dir / f"sample_{idx}.parquet")


def _load_matrix_db(matrix_db: Path):
    conn = duckdb.connect(str(matrix_db), read_only=True)
    try:
        samples = conn.execute(
            "SELECT sample_idx, sample_name FROM matrix_db_samples ORDER BY sample_idx"
        ).fetchall()
        scaffolds = conn.execute(
            "SELECT scaffold_idx, genome, chrom, vector_length FROM matrix_db_scaffolds ORDER BY scaffold_idx"
        ).fetchall()
        matrices = conn.execute(
            "SELECT sample_idx, scaffold_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob "
            "FROM matrix_db_sample_scaffold_matrices ORDER BY sample_idx, scaffold_idx"
        ).fetchall()
    finally:
        conn.close()
    return samples, scaffolds, matrices


def test_build_matrix_db(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    _write_profiles(profile_dir)
    progress_events: list[dict[str, object]] = []

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert summary.profile_files == 3
    assert summary.sample_count == 3
    assert summary.scaffold_count == 2
    assert summary.count_dtype == "uint16"
    assert summary.stored_rows == 6
    assert progress_events[0]["phase"] == "start"
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 6
    assert progress_events[-1]["total"] == 6
    assert any(event["phase"] == "processing" for event in progress_events)
    assert any(event["sample_name"] == "sample_a" for event in progress_events)
    assert any(event["scaffold"] == "chr1" for event in progress_events)

    samples, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert scaffolds == [
        (0, "genome1", "chr1", 3),
        (1, "genome2", "chr2", 3),
    ]
    conn = duckdb.connect(str(matrix_db), read_only=True)
    try:
        metadata = dict(conn.execute("SELECT key, value FROM matrix_db_metadata").fetchall())
    finally:
        conn.close()
    assert metadata["matrix_value_semantics"] == mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS
    assert metadata["coverage_filter_min_cov"] == str(mp.MATRIX_BUILD_MIN_COV)

    unpacked = {}
    for sample_idx, scaffold_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in matrices:
        unpacked[(sample_idx, scaffold_idx)] = mp._unpack_matrix(
            bytes(matrix_blob),
            count_dtype,
            (matrix_rows, matrix_cols),
        ).tolist()

    assert unpacked[(0, 0)] == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]]
    assert unpacked[(1, 0)] == [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert unpacked[(2, 0)] == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert unpacked[(0, 1)] == [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    assert unpacked[(1, 1)] == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
    assert unpacked[(2, 1)] == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]


def test_build_matrix_db_with_small_commit_batches(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_small_batches.duckdb"
    _write_profiles(profile_dir)

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        commit_batch_gb=1e-9,
    )

    assert summary.stored_rows == 6
    samples, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert len(samples) == 3
    assert len(scaffolds) == 2
    assert len(matrices) == 6


def test_build_matrix_db_with_optional_bed_file(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_bed.duckdb"
    bed_file = tmp_path / "genomes.bed"
    _write_profiles_one_based(profile_dir)
    bed_file.write_text("chr1\t0\t2\nchr1\t2\t5\nchr2\t5\t7\nchr2\t7\t10\n")

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        bed_file=bed_file,
        count_dtype="uint16",
        memory_limit_gb=1.0,
    )

    assert summary.profile_files == 2
    assert summary.scaffold_count == 2

    _samples, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert scaffolds == [
        (0, "genome1", "chr1", 5),
        (1, "genome2", "chr2", 5),
    ]

    unpacked = {}
    for sample_idx, scaffold_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in matrices:
        unpacked[(sample_idx, scaffold_idx)] = mp._unpack_matrix(
            bytes(matrix_blob),
            count_dtype,
            (matrix_rows, matrix_cols),
        ).tolist()

    assert unpacked[(0, 0)] == [
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 1, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    assert unpacked[(1, 0)] == [
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]


def test_matrix_compare_matches_pairwise_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)
    progress_events: list[dict[str, object]] = []

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
        progress_callback=progress_events.append,
    )

    assert summary.requested_pairs == 3
    assert summary.backend == "numpy"
    assert summary.anchor_groups == 2
    assert progress_events[0]["phase"] == "start"
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 6
    assert progress_events[-1]["total"] == 6
    assert any(event["phase"] == "processing" for event in progress_events)
    assert any(event["anchor_name"] == "sample_a" for event in progress_events)
    assert any(event["scaffold"] == "chr1" for event in progress_events)
    actual = pl.read_parquet(output_file).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select([
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
            ])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])

    assert actual.equals(expected)


def test_matrix_compare_with_ibs_matches_pairwise_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare_ibs.parquet"
    _write_profiles(profile_dir)

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani+ibs",
    )

    assert summary.requested_pairs == 3
    actual = pl.read_parquet(output_file).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani+ibs",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select([
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
                "max_consecutive_length",
            ])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])

    assert actual.equals(expected)


def test_matrix_compare_loads_targets_in_batches(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    original_load = mp._load_sample_scaffold_matrices
    call_sizes: list[int] = []

    def tracking_load(conn, scaffold_idx, sample_ids, vector_length):
        call_sizes.append(len(sample_ids))
        return original_load(conn, scaffold_idx, sample_ids, vector_length)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
        position_tile_size=None,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    monkeypatch.setattr(mp, "_load_sample_scaffold_matrices", tracking_load)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)

    mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert call_sizes
    assert 2 in call_sizes
    assert max(call_sizes) == 2


def test_cli_matrix_build_and_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--commit-batch-gb",
            "0.01",
        ],
    )
    assert build_result.exit_code == 0
    assert "scaffolds=2" in build_result.output

    compare_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare",
            "--matrix-db-file",
            str(matrix_db),
            "--output-file",
            str(output_file),
            "--min-cov",
            str(mp.MATRIX_BUILD_MIN_COV),
            "--memory-limit-gb",
            "1",
            "--calculate",
            "ani",
        ],
    )
    assert compare_result.exit_code == 0
    assert output_file.exists()
    assert "requested_pairs=3" in compare_result.output
    assert "anchor_groups=2" in compare_result.output

    out = pl.read_parquet(output_file)
    assert "max_consecutive_length" not in out.columns


def test_matrix_compare_torch_backend_requires_torch(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    def _raise_import_error(_name):
        raise ImportError("torch missing")

    monkeypatch.setattr(mp.importlib, "import_module", _raise_import_error)

    with pytest.raises(mp.TorchBackendMissingError, match="Torch backend requested"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            backend="torch",
            memory_limit_gb=1.0,
            calculate="ani",
        )


def test_matrix_compare_requires_matching_min_cov(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    with pytest.raises(ValueError, match="fixed coverage filter min_cov"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=1,
            memory_limit_gb=1.0,
            backend="numpy",
            calculate="ani",
        )


def test_matrix_compare_torch_reuses_target_chunks_across_anchors(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.parquet"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = SimpleNamespace()

    original_load = mp._load_sample_scaffold_matrices
    call_sizes: list[int] = []

    def tracking_load(conn, scaffold_idx, sample_ids, vector_length):
        call_sizes.append(len(sample_ids))
        return original_load(conn, scaffold_idx, sample_ids, vector_length)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
        position_tile_size=None,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    def identity_prepare(compute_backend, matrix, matrix_value_semantics):
        return matrix

    def zero_compare(
        compute_backend,
        anchor_torch,
        target_torch,
        vector_length: int,
        tile_size: int,
        matrix_value_semantics: str,
        min_cov: int,
        need_ibs: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
        )

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_load_sample_scaffold_matrices", tracking_load)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch", zero_compare)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch",
        calculate="ani",
    )

    assert summary.target_chunks == 3
    assert call_sizes.count(2) == 2
    assert call_sizes.count(1) == 7
