from pathlib import Path

import polars as pl
from click.testing import CliRunner
import pytest

from zipstrain import compare as cp
from zipstrain import light
from zipstrain import light_cli


def _write_reference(reference_fasta: Path) -> None:
    reference_fasta.write_text(
        ">chr1\nAAAAAA\n"
    )


def _write_stb(stb_path: Path) -> None:
    stb_path.write_text("chr1\tgenome1\n")


def _write_profiles(profile_1_path: Path, profile_2_path: Path) -> None:
    p1 = pl.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "pos": [1, 2, 3, 4, 5, 6],
            "gene": ["gene1", "gene1", "gene1", "gene2", "gene2", "gene2"],
            "genome": ["genome1"] * 6,
            "A": [10, 10, 0, 10, 0, 10],
            "T": [0, 0, 0, 0, 0, 0],
            "C": [0, 0, 10, 0, 0, 0],
            "G": [0, 0, 0, 0, 10, 0],
        }
    )
    p2 = pl.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "pos": [1, 2, 3, 4, 5, 6],
            "gene": ["gene1", "gene1", "gene1", "gene2", "gene2", "gene2"],
            "genome": ["genome1"] * 6,
            "A": [10, 0, 0, 10, 10, 10],
            "T": [0, 10, 0, 0, 0, 0],
            "C": [0, 0, 10, 0, 0, 0],
            "G": [0, 0, 0, 0, 0, 0],
        }
    )
    p1.write_parquet(profile_1_path)
    p2.write_parquet(profile_2_path)


def _write_light_bam_mode_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t6\n")
    gene_range_table = tmp_path / "gene_ranges.tsv"
    gene_range_table.write_text("gene1\tchr1\t1\t3\ngene2\tchr1\t4\t6\n")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": [0, 5, 10], "max_error_count": [0, 0, 0]}).write_parquet(null_model_file)
    return bam_file, bed_file, gene_range_table, stb_file, null_model_file


def _install_fake_light_bam_profiler(monkeypatch: pytest.MonkeyPatch) -> None:
    def _fake_profile_bam(
        bed_file: str,
        bam_file: str,
        gene_range_table: str,
        stb,
        null_model,
        output_dir: str,
        num_chunks: int = 24,
        max_concurrency: int = 4,
    ) -> None:
        _ = (bed_file, gene_range_table, stb, null_model, num_chunks, max_concurrency)
        out = Path(output_dir) / f"{Path(bam_file).stem}_profile.parquet"
        pl.DataFrame(
            {
                "chrom": ["chr1"] * 6,
                "pos": [1, 2, 3, 4, 5, 6],
                "gene": ["gene1", "gene1", "gene1", "gene2", "gene2", "gene2"],
                "genome": ["genome1"] * 6,
                "A": [10, 10, 0, 10, 0, 10],
                "T": [0, 0, 0, 0, 0, 0],
                "C": [0, 0, 10, 0, 0, 0],
                "G": [0, 0, 0, 0, 10, 0],
            }
        ).write_parquet(out)

    monkeypatch.setattr(light.pf, "profile_bam", _fake_profile_bam)


def test_light_profile_duckdb_folder_tables(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    profile_1 = tmp_path / "p1.parquet"
    _write_reference(reference_fasta)
    _write_profiles(profile_1, tmp_path / "p2.parquet")

    output_dir = tmp_path / "p1_light_profile"
    summary = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=output_dir,
        profile_engine="duckdb",
    )

    output_db = output_dir / "profile.duckdb"
    assert output_db.exists()
    assert summary.coverage_rows == 6
    assert summary.snp_rows == 2

    import duckdb

    con = duckdb.connect(str(output_db), read_only=True)
    try:
        cov_count = con.execute("SELECT COUNT(*) FROM coverage").fetchone()[0]
        snp_count = con.execute("SELECT COUNT(*) FROM snp").fetchone()[0]
        genome_count = con.execute("SELECT COUNT(*) FROM genomes").fetchone()[0]
    finally:
        con.close()

    assert cov_count == 6
    assert snp_count == 2
    assert genome_count == 1


def test_light_profile_bundle_layouts(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    profile_1 = tmp_path / "p1.parquet"
    _write_reference(reference_fasta)
    _write_profiles(profile_1, tmp_path / "p2.parquet")

    duckdb_dir = tmp_path / "p1_duckdb_profile"
    summary_duckdb = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir,
        profile_engine="duckdb",
    )
    assert summary_duckdb.coverage_rows == 6
    assert summary_duckdb.snp_rows == 2
    assert (duckdb_dir / "profile.duckdb").exists()
    assert len(list(duckdb_dir.glob("*.duckdb"))) == 1
    assert not (duckdb_dir / "coverage.parquet").exists()
    assert not (duckdb_dir / "snp.parquet").exists()

    polars_dir = tmp_path / "p1_polars_profile"
    summary_polars = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=polars_dir,
        profile_engine="polars",
    )
    assert summary_polars.coverage_rows == 6
    assert summary_polars.snp_rows == 2
    assert (polars_dir / "coverage.parquet").exists()
    assert (polars_dir / "snp.parquet").exists()
    assert len(list(polars_dir.glob("*.duckdb"))) == 0


def test_light_profile_duckdb_path_is_duckdb_native(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    profile_1 = tmp_path / "p1.parquet"
    _write_reference(reference_fasta)
    _write_profiles(profile_1, tmp_path / "p2.parquet")

    def _fail_split(*args, **kwargs):
        _ = (args, kwargs)
        raise AssertionError("duckdb profile mode should not call _split_light_tables")

    monkeypatch.setattr(light, "_split_light_tables", _fail_split)

    duckdb_dir = tmp_path / "p1_duckdb_native"
    summary = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir,
        profile_engine="duckdb",
    )

    assert summary.coverage_rows == 6
    assert summary.snp_rows == 2
    assert (duckdb_dir / "profile.duckdb").exists()


def test_light_profile_min_cov_filters_coverage_and_snp(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    reference_fasta.write_text(">chr1\nAAAA\n")
    profile = tmp_path / "profile.parquet"
    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "pos": [1, 2, 3, 4],
            "gene": ["gene1", "gene1", "gene1", "gene1"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            # cov = 5 (drop), 6 (keep), 7 (keep SNP), 4 (drop)
            "A": [5, 6, 0, 0],
            "T": [0, 0, 0, 4],
            "C": [0, 0, 7, 0],
            "G": [0, 0, 0, 0],
        }
    ).write_parquet(profile)

    polars_dir = tmp_path / "profile_polars"
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile,
        reference_fasta=reference_fasta,
        output_dir=polars_dir,
        profile_engine="polars",
        min_cov=5,
    )
    cov_polars = pl.read_parquet(polars_dir / "coverage.parquet").sort("pos")
    snp_polars = pl.read_parquet(polars_dir / "snp.parquet").sort("pos")
    assert cov_polars["pos"].to_list() == [2, 3]
    assert snp_polars["pos"].to_list() == [3]

    duckdb_dir = tmp_path / "profile_duckdb"
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir,
        profile_engine="duckdb",
        min_cov=5,
    )

    import duckdb

    con = duckdb.connect(str(duckdb_dir / "profile.duckdb"), read_only=True)
    try:
        cov_positions = [row[0] for row in con.execute("SELECT pos FROM coverage ORDER BY pos").fetchall()]
        snp_positions = [row[0] for row in con.execute("SELECT pos FROM snp ORDER BY pos").fetchall()]
        min_cov_meta = con.execute("SELECT min_cov FROM metadata").fetchone()[0]
    finally:
        con.close()

    assert cov_positions == [2, 3]
    assert snp_positions == [3]
    assert min_cov_meta == 5


def test_light_profile_from_bam_writes_expected_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_light_bam_profiler(monkeypatch)
    reference_fasta = tmp_path / "ref.fasta"
    _write_reference(reference_fasta)
    bam_file, bed_file, gene_range_table, stb_file, null_model_file = _write_light_bam_mode_inputs(tmp_path)

    duckdb_dir = tmp_path / "bam_duckdb_profile"
    duckdb_summary = light.build_light_profile_bundle_from_bam(
        bam_file=bam_file,
        bed_file=bed_file,
        gene_range_table=gene_range_table,
        stb_file=stb_file,
        null_model=null_model_file,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir,
        profile_engine="duckdb",
        min_cov=5,
    )
    assert duckdb_summary.coverage_rows == 6
    assert duckdb_summary.snp_rows == 2
    assert (duckdb_dir / "profile.duckdb").exists()

    polars_dir = tmp_path / "bam_polars_profile"
    polars_summary = light.build_light_profile_bundle_from_bam(
        bam_file=bam_file,
        bed_file=bed_file,
        gene_range_table=gene_range_table,
        stb_file=stb_file,
        null_model=null_model_file,
        reference_fasta=reference_fasta,
        output_dir=polars_dir,
        profile_engine="polars",
        min_cov=5,
        num_chunks=2,
        max_concurrency=2,
    )
    assert polars_summary.coverage_rows == 6
    assert polars_summary.snp_rows == 2
    assert (polars_dir / "coverage.parquet").exists()
    assert (polars_dir / "snp.parquet").exists()


def test_light_cli_profile_from_bam(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_light_bam_profiler(monkeypatch)
    reference_fasta = tmp_path / "ref.fasta"
    _write_reference(reference_fasta)
    bam_file, bed_file, gene_range_table, stb_file, null_model_file = _write_light_bam_mode_inputs(tmp_path)
    out_dir = tmp_path / "cli_bam_profile"

    runner = CliRunner()
    res = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--bam-file",
            str(bam_file),
            "--bed-file",
            str(bed_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model_file),
            "--gene-range-table",
            str(gene_range_table),
            "--reference-fasta",
            str(reference_fasta),
            "--engine",
            "polars",
            "--num-chunks",
            "2",
            "--max-concurrency",
            "2",
            "--min-cov",
            "5",
            "--output-dir",
            str(out_dir),
        ],
    )
    assert res.exit_code == 0, res.output
    assert (out_dir / "coverage.parquet").exists()
    assert (out_dir / "snp.parquet").exists()


def test_light_cli_profile_from_bam_requires_null_model(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    _write_reference(reference_fasta)
    bam_file, bed_file, gene_range_table, stb_file, _ = _write_light_bam_mode_inputs(tmp_path)
    out_dir = tmp_path / "cli_bam_profile_missing_null"

    runner = CliRunner()
    res = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--bam-file",
            str(bam_file),
            "--bed-file",
            str(bed_file),
            "--stb-file",
            str(stb_file),
            "--gene-range-table",
            str(gene_range_table),
            "--reference-fasta",
            str(reference_fasta),
            "--output-dir",
            str(out_dir),
        ],
    )
    assert res.exit_code != 0
    assert "--null-model" in res.output


def test_light_cli_profile_from_bam_duckdb_engine(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_light_bam_profiler(monkeypatch)
    reference_fasta = tmp_path / "ref.fasta"
    _write_reference(reference_fasta)
    bam_file, bed_file, gene_range_table, stb_file, null_model_file = _write_light_bam_mode_inputs(tmp_path)
    out_dir = tmp_path / "cli_bam_profile_duckdb"

    runner = CliRunner()
    res = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--bam-file",
            str(bam_file),
            "--bed-file",
            str(bed_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model_file),
            "--gene-range-table",
            str(gene_range_table),
            "--reference-fasta",
            str(reference_fasta),
            "--engine",
            "duckdb",
            "--min-cov",
            "5",
            "--output-dir",
            str(out_dir),
        ],
    )
    assert res.exit_code == 0, res.output
    assert (out_dir / "profile.duckdb").exists()


def test_light_compare_matches_regular(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1.parquet"
    profile_2 = tmp_path / "p2.parquet"
    _write_reference(reference_fasta)
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)

    duckdb_dir_1 = tmp_path / "p1_profile_duckdb"
    duckdb_dir_2 = tmp_path / "p2_profile_duckdb"
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir_1,
        profile_engine="duckdb",
    )
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_2,
        reference_fasta=reference_fasta,
        output_dir=duckdb_dir_2,
        profile_engine="duckdb",
    )

    regular = (
        cp.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=1,
            min_gene_compare_len=1,
            engine="duckdb",
            stb_file=stb_path,
        )
        .collect()
        .sort("genome")
    )

    light_duckdb = (
        light.compare_light_genomes(
            profile_db_1=duckdb_dir_1,
            profile_db_2=duckdb_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="duckdb",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )

    light_polars = (
        light.compare_light_genomes(
            profile_db_1=duckdb_dir_1,
            profile_db_2=duckdb_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="polars",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )

    profile_dir_1 = tmp_path / "p1_profile_polars"
    profile_dir_2 = tmp_path / "p2_profile_polars"
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=profile_dir_1,
        profile_engine="polars",
    )
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_2,
        reference_fasta=reference_fasta,
        output_dir=profile_dir_2,
        profile_engine="polars",
    )
    light_from_parquet_bundle = (
        light.compare_light_genomes(
            profile_db_1=profile_dir_1,
            profile_db_2=profile_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="duckdb",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )
    light_from_parquet_bundle_polars = (
        light.compare_light_genomes(
            profile_db_1=profile_dir_1,
            profile_db_2=profile_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="polars",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )

    # Normalize floating-point precision for robust equality.
    float_cols = [
        col for col, dtype in regular.schema.items() if dtype in {pl.Float32, pl.Float64}
    ]
    if float_cols:
        regular = regular.with_columns([pl.col(c).round(8) for c in float_cols])
        light_duckdb = light_duckdb.with_columns([pl.col(c).round(8) for c in float_cols])
        light_polars = light_polars.with_columns([pl.col(c).round(8) for c in float_cols])
        light_from_parquet_bundle = light_from_parquet_bundle.with_columns([pl.col(c).round(8) for c in float_cols])
        light_from_parquet_bundle_polars = light_from_parquet_bundle_polars.with_columns(
            [pl.col(c).round(8) for c in float_cols]
        )

    light_from_parquet_bundle = light_from_parquet_bundle.select(
        [
            pl.col(col).cast(dtype)
            for col, dtype in regular.schema.items()
        ]
    )
    light_from_parquet_bundle_polars = light_from_parquet_bundle_polars.select(
        [
            pl.col(col).cast(dtype)
            for col, dtype in regular.schema.items()
        ]
    )

    assert regular.equals(light_duckdb)
    assert regular.equals(light_polars)
    assert regular.equals(light_from_parquet_bundle)
    assert regular.equals(light_from_parquet_bundle_polars)

    ani_cols = ["genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    light_duckdb_ani = (
        light.compare_light_genomes(
            profile_db_1=duckdb_dir_1,
            profile_db_2=duckdb_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="duckdb",
            calculate="ani",
        )
        .collect()
        .sort("genome")
    )
    light_polars_ani = (
        light.compare_light_genomes(
            profile_db_1=duckdb_dir_1,
            profile_db_2=duckdb_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="polars",
            calculate="ani",
        )
        .collect()
        .sort("genome")
    )

    regular_ani = regular.select(ani_cols)
    if "genome_pop_ani" in regular_ani.columns:
        regular_ani = regular_ani.with_columns(pl.col("genome_pop_ani").round(8))
        light_duckdb_ani = light_duckdb_ani.with_columns(pl.col("genome_pop_ani").round(8))
        light_polars_ani = light_polars_ani.with_columns(pl.col("genome_pop_ani").round(8))

    assert light_duckdb_ani.columns == ani_cols
    assert light_polars_ani.columns == ani_cols
    assert regular_ani.equals(light_duckdb_ani)
    assert regular_ani.equals(light_polars_ani)


def test_light_compare_rejects_direct_duckdb_file_input(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1.parquet"
    profile_2 = tmp_path / "p2.parquet"
    _write_reference(reference_fasta)
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)

    profile_dir_1 = tmp_path / "p1_profile_duckdb"
    profile_dir_2 = tmp_path / "p2_profile_duckdb"
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        reference_fasta=reference_fasta,
        output_dir=profile_dir_1,
        profile_engine="duckdb",
    )
    light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_2,
        reference_fasta=reference_fasta,
        output_dir=profile_dir_2,
        profile_engine="duckdb",
    )

    try:
        light.compare_light_genomes(
            profile_db_1=profile_dir_1 / "profile.duckdb",
            profile_db_2=profile_dir_2,
            stb_file=stb_path,
            min_cov=1,
            min_gene_compare_len=1,
            engine="duckdb",
            calculate="ani",
        ).collect()
    except ValueError as exc:
        assert "Provide a profile directory" in str(exc)
    else:
        raise AssertionError("Expected direct .duckdb profile input to be rejected")


def test_light_cli_profile_and_compare(tmp_path: Path) -> None:
    reference_fasta = tmp_path / "ref.fasta"
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1.parquet"
    profile_2 = tmp_path / "p2.parquet"
    _write_reference(reference_fasta)
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)

    profile_dir_1 = tmp_path / "p1_light_profile"
    profile_dir_2 = tmp_path / "p2_light_profile"
    out = tmp_path / "comparison.parquet"

    runner = CliRunner()
    res_1 = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--profile-parquet",
            str(profile_1),
            "--reference-fasta",
            str(reference_fasta),
            "--output-dir",
            str(profile_dir_1),
        ],
    )
    assert res_1.exit_code == 0, res_1.output

    res_2 = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--profile-parquet",
            str(profile_2),
            "--reference-fasta",
            str(reference_fasta),
            "--engine",
            "polars",
            "--output-dir",
            str(profile_dir_2),
        ],
    )
    assert res_2.exit_code == 0, res_2.output

    res_cmp = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(profile_dir_1),
            "--profile-2",
            str(profile_dir_2),
            "--stb-file",
            str(stb_path),
            "--output-file",
            str(out),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
        ],
    )
    assert res_cmp.exit_code == 0, res_cmp.output
    assert out.exists()

    out_df = pl.read_parquet(out)
    assert "sample_1" in out_df.columns
    assert "sample_2" in out_df.columns

    out_polars = tmp_path / "comparison_polars.parquet"
    res_cmp_polars = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(profile_dir_1),
            "--profile-2",
            str(profile_dir_2),
            "--stb-file",
            str(stb_path),
            "--output-file",
            str(out_polars),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--engine",
            "polars",
        ],
    )
    assert res_cmp_polars.exit_code == 0, res_cmp_polars.output
    assert out_polars.exists()

    out_ani = tmp_path / "comparison_ani.parquet"
    res_cmp_ani = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(profile_dir_1),
            "--profile-2",
            str(profile_dir_2),
            "--stb-file",
            str(stb_path),
            "--output-file",
            str(out_ani),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--calculate",
            "ani",
        ],
    )
    assert res_cmp_ani.exit_code == 0, res_cmp_ani.output
    out_ani_df = pl.read_parquet(out_ani)
    assert out_ani_df.columns == ["genome", "total_positions", "share_allele_pos", "genome_pop_ani", "sample_1", "sample_2"]

    out_no_stb = tmp_path / "comparison_no_stb.parquet"
    res_cmp_no_stb = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(profile_dir_1),
            "--profile-2",
            str(profile_dir_2),
            "--output-file",
            str(out_no_stb),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--calculate",
            "ani",
        ],
    )
    assert res_cmp_no_stb.exit_code != 0

    direct_db_input = profile_dir_1 / "profile.duckdb"
    out_direct_db = tmp_path / "comparison_direct_db.parquet"
    res_cmp_direct_db = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(direct_db_input),
            "--profile-2",
            str(profile_dir_2),
            "--stb-file",
            str(stb_path),
            "--output-file",
            str(out_direct_db),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--calculate",
            "ani",
        ],
    )
    assert res_cmp_direct_db.exit_code != 0
