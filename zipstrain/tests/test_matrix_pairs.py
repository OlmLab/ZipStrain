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


def _write_profiles_multiscaffold_same_genome(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 0, 1],
            "gene": ["gene1", "gene1", "gene2", "gene2"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            "A": [6, 6, 6, 6],
            "T": [0, 0, 0, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 0, 1],
            "gene": ["gene1", "gene1", "gene2", "gene2"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            "A": [7, 7, 7, 7],
            "T": [0, 0, 0, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 0],
        }
    )
    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")


def _write_invalid_append_profile(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    invalid = pl.DataFrame(
        {
            "chrom": ["chr9"],
            "pos": [0],
            "gene": ["gene_bad"],
            "genome": ["genome9"],
            "A": [6],
            "T": [0],
            "C": [0],
            "G": [0],
        }
    )
    invalid.write_parquet(profile_dir / "sample_bad.parquet")


def _load_matrix_db(matrix_db: Path):
    conn = duckdb.connect(str(matrix_db), read_only=True)
    try:
        samples = conn.execute(
            "SELECT sample_idx, sample_name FROM matrix_db_samples ORDER BY sample_idx"
        ).fetchall()
        genomes = conn.execute(
            "SELECT genome_idx, genome, matrix_length, true_length, scaffold_count "
            "FROM matrix_db_genomes ORDER BY genome_idx"
        ).fetchall()
        scaffolds = conn.execute(
            "SELECT genome_idx, scaffold_ordinal, genome, chrom, axis_start, axis_end, vector_length "
            "FROM matrix_db_genome_scaffolds ORDER BY genome_idx, scaffold_ordinal"
        ).fetchall()
        matrices = conn.execute(
            "SELECT sample_idx, genome_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob "
            "FROM matrix_db_sample_genome_matrices ORDER BY sample_idx, genome_idx"
        ).fetchall()
    finally:
        conn.close()
    return samples, genomes, scaffolds, matrices


def _load_matrix_compare_db(compare_db: Path):
    conn = duckdb.connect(str(compare_db), read_only=True)
    try:
        metadata = dict(conn.execute("SELECT key, value FROM matrix_compare_metadata ORDER BY key").fetchall())
        completed_pairs = conn.execute(
            """
            SELECT DISTINCT sample_idx_1, sample_idx_2
            FROM matrix_compare_completed_pair_genomes
            ORDER BY sample_idx_1, sample_idx_2
            """
        ).fetchall()
        results = (
            pl.from_arrow(
                conn.execute(
                    """
                    SELECT sample_idx_1, sample_idx_2, sample_1, sample_2, genome_idx, genome,
                           total_positions, share_allele_pos, genome_pop_ani, max_consecutive_length
                    FROM matrix_compare_results
                    ORDER BY sample_idx_1, sample_idx_2, genome_idx
                    """
                ).fetch_arrow_table()
            )
        )
    finally:
        conn.close()
    return metadata, completed_pairs, results


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
    assert all(event["scaffold"] == "" for event in progress_events if "scaffold" in event)

    samples, genomes, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert genomes == [
        (0, "genome1", 3, 3, 1),
        (1, "genome2", 3, 3, 1),
    ]
    assert scaffolds == [
        (0, 0, "genome1", "chr1", 0, 2, 3),
        (1, 0, "genome2", "chr2", 0, 2, 3),
    ]
    conn = duckdb.connect(str(matrix_db), read_only=True)
    try:
        metadata = dict(conn.execute("SELECT key, value FROM matrix_db_metadata").fetchall())
    finally:
        conn.close()
    assert metadata["matrix_value_semantics"] == mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS
    assert metadata["coverage_filter_min_cov"] == str(mp.MATRIX_BUILD_MIN_COV)
    assert metadata["layout"] == "per_sample_per_genome_dense_matrix"
    assert metadata["separator_rows_between_scaffolds"] == "1"

    unpacked = {}
    for sample_idx, genome_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in matrices:
        unpacked[(sample_idx, genome_idx)] = mp._unpack_matrix(
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
    samples, genomes, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert len(samples) == 3
    assert len(genomes) == 2
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

    _samples, genomes, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert genomes == [
        (0, "genome1", 5, 5, 1),
        (1, "genome2", 5, 5, 1),
    ]
    assert scaffolds == [
        (0, 0, "genome1", "chr1", 0, 4, 5),
        (1, 0, "genome2", "chr2", 0, 4, 5),
    ]
    unpacked = {}
    for sample_idx, genome_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in matrices:
        unpacked[(sample_idx, genome_idx)] = mp._unpack_matrix(
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


def test_append_matrix_db_success(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.duckdb"

    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    initial_summary = mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    assert initial_summary.sample_count == 2

    progress_events: list[dict[str, object]] = []
    append_summary = mp.append_matrix_db(
        profile_dir=append_profile_dir,
        matrix_db_file=matrix_db,
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert append_summary.appended_sample_count == 1
    assert append_summary.total_sample_count == 3
    assert append_summary.stored_rows == 2
    assert progress_events[0]["phase"] == "start"
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 2
    assert progress_events[-1]["total"] == 2

    samples, genomes, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert len(genomes) == 2
    assert len(scaffolds) == 2
    assert len(matrices) == 6


def test_append_matrix_db_rejects_incompatible_profile_without_mutation(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_bad"
    matrix_db = tmp_path / "matrix.duckdb"

    _write_profiles(initial_profile_dir)
    (initial_profile_dir / "sample_c.parquet").unlink()
    _write_invalid_append_profile(append_profile_dir)

    mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    samples_before, genomes_before, scaffolds_before, matrices_before = _load_matrix_db(matrix_db)

    with pytest.raises(ValueError, match="does not match the existing matrix DB contract"):
        mp.append_matrix_db(
            profile_dir=append_profile_dir,
            matrix_db_file=matrix_db,
            memory_limit_gb=1.0,
        )

    samples_after, genomes_after, scaffolds_after, matrices_after = _load_matrix_db(matrix_db)
    assert samples_after == samples_before
    assert genomes_after == genomes_before
    assert scaffolds_after == scaffolds_before
    assert matrices_after == matrices_before


def test_matrix_compare_matches_pairwise_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
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
    assert all(event["scaffold"] == "" for event in progress_events if "scaffold" in event)
    metadata, completed_pairs, results = _load_matrix_compare_db(output_file)
    assert metadata["calculate"] == "ani"
    assert len(completed_pairs) == 3
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

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
    output_file = tmp_path / "matrix_compare_ibs.duckdb"
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
    metadata, completed_pairs, results = _load_matrix_compare_db(output_file)
    assert metadata["calculate"] == "ani+ibs"
    assert len(completed_pairs) == 3
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])

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
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    original_load = mp._load_sample_genome_matrices
    call_sizes: list[int] = []

    def tracking_load(conn, genome_idx, sample_ids, matrix_length, dtype_name):
        call_sizes.append(len(sample_ids))
        return original_load(conn, genome_idx, sample_ids, matrix_length, dtype_name)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
        position_tile_size=None,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    monkeypatch.setattr(mp, "_load_sample_genome_matrices", tracking_load)
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
    output_file = tmp_path / "matrix_compare.duckdb"
    export_file = tmp_path / "matrix_compare.parquet"
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

    export_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare-export",
            "--matrix-compare-db-file",
            str(output_file),
            "--output-file",
            str(export_file),
        ],
    )
    assert export_result.exit_code == 0
    out = pl.read_parquet(export_file)
    assert "max_consecutive_length" not in out.columns


def test_cli_append_matrix_db(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.duckdb"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(initial_profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert build_result.exit_code == 0

    append_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "append-matrix-db",
            "--profile-dir",
            str(append_profile_dir),
            "--matrix-db-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert append_result.exit_code == 0
    assert "appended_samples=1" in append_result.output
    assert "total_samples=3" in append_result.output

    samples, _genomes, _scaffolds, matrices = _load_matrix_db(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert len(matrices) == 6


def test_matrix_compare_resumes_after_matrix_append(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.duckdb"
    compare_db = tmp_path / "matrix_compare.duckdb"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    first_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )
    assert first_summary.requested_pairs == 1

    mp.append_matrix_db(
        profile_dir=append_profile_dir,
        matrix_db_file=matrix_db,
        memory_limit_gb=1.0,
    )
    second_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )
    assert second_summary.requested_pairs == 2

    _metadata, completed_pairs, results = _load_matrix_compare_db(compare_db)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=(append_profile_dir if sample_1 == "sample_c" else initial_profile_dir) / f"{sample_1}.parquet",
                mpile_contig_2=(append_profile_dir if sample_2 == "sample_c" else initial_profile_dir) / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(sample_1=pl.lit(sample_1), sample_2=pl.lit(sample_2))
            .select(["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


def test_matrix_compare_resumes_after_interruption(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    compare_db = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )

    original_mark = mp._mark_completed_pair_genomes
    mark_calls = {"count": 0}

    def interrupt_after_first_commit(compare_conn, completed_pair_genomes):
        mark_calls["count"] += 1
        original_mark(compare_conn, completed_pair_genomes)
        if mark_calls["count"] == 2:
            raise RuntimeError("simulated interruption")

    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", interrupt_after_first_commit)

    with pytest.raises(RuntimeError, match="simulated interruption"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=compare_db,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            backend="numpy",
            calculate="ani",
        )

    metadata, completed_pairs_partial, partial_results = _load_matrix_compare_db(compare_db)
    assert metadata["calculate"] == "ani"
    assert completed_pairs_partial == [(0, 1), (0, 2)]
    assert partial_results.height > 0

    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", original_mark)
    resume_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert resume_summary.requested_pairs == 3
    _metadata, completed_pairs, results = _load_matrix_compare_db(compare_db)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

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
            .with_columns(sample_1=pl.lit(sample_1), sample_2=pl.lit(sample_2))
            .select(["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


def test_cli_append_matrix_db_failure_keeps_existing_db(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    bad_profile_dir = tmp_path / "profiles_bad"
    matrix_db = tmp_path / "matrix.duckdb"
    _write_profiles(initial_profile_dir)
    (initial_profile_dir / "sample_c.parquet").unlink()
    _write_invalid_append_profile(bad_profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(initial_profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert build_result.exit_code == 0
    samples_before, _genomes_before, _scaffolds_before, matrices_before = _load_matrix_db(matrix_db)

    append_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "append-matrix-db",
            "--profile-dir",
            str(bad_profile_dir),
            "--matrix-db-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert append_result.exit_code != 0
    assert append_result.exception is not None
    assert "does not match the existing matrix DB contract" in str(append_result.exception)

    samples_after, _genomes_after, _scaffolds_after, matrices_after = _load_matrix_db(matrix_db)
    assert samples_after == samples_before
    assert matrices_after == matrices_before


def test_matrix_compare_torch_backend_requires_torch(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
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


def test_matrix_compare_ignores_requested_min_cov(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=1,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert summary.requested_pairs == 3
    _metadata, completed_pairs, _results = _load_matrix_compare_db(output_file)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]


def test_matrix_compare_torch_reuses_target_chunks_across_anchors(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = SimpleNamespace()

    original_load = mp._load_sample_genome_matrices
    call_sizes: list[int] = []

    def tracking_load(conn, genome_idx, sample_ids, matrix_length, dtype_name):
        call_sizes.append(len(sample_ids))
        return original_load(conn, genome_idx, sample_ids, matrix_length, dtype_name)

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
        need_ibs: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
        )

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_load_sample_genome_matrices", tracking_load)
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


def test_matrix_compare_torch_resumes_after_interruption(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = SimpleNamespace()

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
        need_ibs: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
        )

    original_mark = mp._mark_completed_pair_genomes
    mark_calls = {"count": 0}

    def interrupt_after_first_commit(compare_conn, completed_pair_genomes):
        mark_calls["count"] += 1
        original_mark(compare_conn, completed_pair_genomes)
        if mark_calls["count"] == 2:
            raise KeyboardInterrupt("simulated ctrl-c")

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch", zero_compare)
    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", interrupt_after_first_commit)

    with pytest.raises(KeyboardInterrupt, match="simulated ctrl-c"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            backend="torch",
            calculate="ani",
        )

    _metadata, completed_pairs_partial, _partial_results = _load_matrix_compare_db(output_file)
    assert completed_pairs_partial == [(0, 1)]

    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", original_mark)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch",
        calculate="ani",
    )

    assert summary.requested_pairs == 9
    _metadata, completed_pairs, _results = _load_matrix_compare_db(output_file)
    assert completed_pairs == [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    ]


def test_build_matrix_db_inserts_separator_rows_for_multiscaffold_genome(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_multi.duckdb"
    _write_profiles_multiscaffold_same_genome(profile_dir)

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )

    assert summary.sample_count == 2
    samples, genomes, scaffolds, matrices = _load_matrix_db(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b"]
    assert genomes == [
        (0, "genome1", 5, 4, 2),
    ]
    assert scaffolds == [
        (0, 0, "genome1", "chr1", 0, 1, 2),
        (0, 1, "genome1", "chr2", 3, 4, 2),
    ]
    unpacked = {}
    for sample_idx, genome_idx, count_dtype, matrix_rows, matrix_cols, matrix_blob in matrices:
        unpacked[(sample_idx, genome_idx)] = mp._unpack_matrix(
            bytes(matrix_blob),
            count_dtype,
            (matrix_rows, matrix_cols),
        ).tolist()
    # row 2 is the synthetic separator row between the two scaffolds
    assert unpacked[(0, 0)][2] == [0, 0, 0, 0]
    assert unpacked[(1, 0)][2] == [0, 0, 0, 0]


def test_matrix_compare_ibs_resets_at_separator_rows(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_multi.duckdb"
    output_file = tmp_path / "matrix_multi_compare.duckdb"
    _write_profiles_multiscaffold_same_genome(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani+ibs",
    )

    assert summary.requested_pairs == 1
    _metadata, completed_pairs, results = _load_matrix_compare_db(output_file)
    assert completed_pairs == [(0, 1)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])
    expected = (
        cp.compare_genomes(
            mpile_contig_1=profile_dir / "sample_a.parquet",
            mpile_contig_2=profile_dir / "sample_b.parquet",
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            genome_scope="all",
            ani_method="popani",
            engine="polars",
            calculate="ani+ibs",
            stb_file=None,
        )
        .collect(engine="streaming")
        .with_columns(
            sample_1=pl.lit("sample_a"),
            sample_2=pl.lit("sample_b"),
        )
        .select(
            [
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
                "max_consecutive_length",
            ]
        )
        .sort(["sample_1", "sample_2", "genome"])
    )

    assert actual.equals(expected)
    assert actual.get_column("max_consecutive_length").to_list() == [2]
