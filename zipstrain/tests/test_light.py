from pathlib import Path

import duckdb
import polars as pl
from click.testing import CliRunner
import pytest

from zipstrain import compare as cp
from zipstrain import light
from zipstrain import light_cli


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



def _write_companion_stats(profile_path: Path) -> None:
    prefix = profile_path.stem[:-len("_profile")] if profile_path.stem.endswith("_profile") else profile_path.stem
    pl.DataFrame(
        {
            "genome": ["genome1", "genome1"],
            "gene": ["gene1", "gene2"],
            "length": [3, 3],
            "breadth": [1.0, 1.0],
            "coverage": [10.0, 10.0],
        }
    ).write_parquet(profile_path.with_name(f"{prefix}_gene_stats.parquet"))
    pl.DataFrame(
        {
            "genome": ["genome1"],
            "coverage": [10.0],
            "breadth": [1.0],
            "genome_length": [6],
            "gap_mean": [0.0],
            "gap_std": [0.0],
            "5x_cov_sites": [6],
            "heterogeneity": [0.0],
            "ber": [1.0],
            "fug": [1.0],
            "reads_mapped": [10],
        }
    ).write_parquet(profile_path.with_name(f"{prefix}_genome_stats.parquet"))



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
        out_dir = Path(output_dir)
        out = out_dir / f"{Path(bam_file).stem}_profile.parquet"
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
        _write_companion_stats(out)

    monkeypatch.setattr(light.pf, "profile_bam", _fake_profile_bam)



def test_light_profile_duckdb_tables_and_stats(tmp_path: Path) -> None:
    profile_1 = tmp_path / "sample_profile.parquet"
    _write_profiles(profile_1, tmp_path / "p2.parquet")
    _write_companion_stats(profile_1)

    output_db = tmp_path / "sample.light_profile.duckdb"
    summary = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile_1,
        output_file=output_db,
        include_gene_stats=True,
        include_genome_stats=True,
    )

    assert output_db.exists()
    assert summary.coverage_rows == 6
    assert summary.gene_stats_rows == 2
    assert summary.genome_stats_rows == 1

    con = duckdb.connect(str(output_db), read_only=True)
    try:
        coverage = con.execute("SELECT chrom, pos, gene, genome, bit_representation FROM mpileup ORDER BY genome, chrom, pos").fetchall()
        tables = {row[0] for row in con.execute("SHOW TABLES").fetchall()}
    finally:
        con.close()

    assert coverage[0] == ("chr1", 1, "gene1", "genome1", 1)
    assert coverage[2] == ("chr1", 3, "gene1", "genome1", 4)
    assert tables == {"mpileup", "gene_stats", "genome_stats"}



def test_light_profile_min_cov_filters_masked_rows(tmp_path: Path) -> None:
    profile = tmp_path / "profile.parquet"
    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "pos": [1, 2, 3, 4],
            "gene": ["gene1", "gene1", "gene1", "gene1"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            "A": [5, 6, 0, 0],
            "T": [0, 0, 0, 4],
            "C": [0, 0, 7, 0],
            "G": [0, 0, 0, 0],
        }
    ).write_parquet(profile)

    output_db = tmp_path / "profile.light.duckdb"
    summary = light.build_light_profile_bundle_from_parquet(
        profile_parquet=profile,
        output_file=output_db,
        min_cov=5,
    )
    assert summary.coverage_rows == 2

    con = duckdb.connect(str(output_db), read_only=True)
    try:
        rows = con.execute("SELECT pos, bit_representation FROM mpileup ORDER BY pos").fetchall()
    finally:
        con.close()

    assert rows == [(2, 1), (3, 4)]



def test_light_profile_from_bam_writes_expected_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_light_bam_profiler(monkeypatch)
    bam_file, bed_file, gene_range_table, stb_file, null_model_file = _write_light_bam_mode_inputs(tmp_path)

    output_db = tmp_path / "bam.light_profile.duckdb"
    summary = light.build_light_profile_bundle_from_bam(
        bam_file=bam_file,
        bed_file=bed_file,
        gene_range_table=gene_range_table,
        stb_file=stb_file,
        null_model=null_model_file,
        output_file=output_db,
        min_cov=5,
        include_gene_stats=True,
        include_genome_stats=True,
        num_chunks=2,
        max_concurrency=2,
    )

    assert summary.coverage_rows == 6
    assert summary.gene_stats_rows == 2
    assert summary.genome_stats_rows == 1
    assert output_db.exists()



def test_light_cli_profile_from_bam_requires_null_model(tmp_path: Path) -> None:
    bam_file, bed_file, gene_range_table, stb_file, _ = _write_light_bam_mode_inputs(tmp_path)
    output_db = tmp_path / "cli_bam_profile_missing_null.duckdb"

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
            "--output-file",
            str(output_db),
        ],
    )
    assert res.exit_code != 0
    assert "--null-model" in res.output



def test_light_compare_matches_regular(tmp_path: Path) -> None:
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1_profile.parquet"
    profile_2 = tmp_path / "p2_profile.parquet"
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)

    profile_db_1 = tmp_path / "p1.light_profile.duckdb"
    profile_db_2 = tmp_path / "p2.light_profile.duckdb"
    light.build_light_profile_bundle_from_parquet(profile_parquet=profile_1, output_file=profile_db_1)
    light.build_light_profile_bundle_from_parquet(profile_parquet=profile_2, output_file=profile_db_2)

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

    light_all = (
        light.compare_light_genomes(
            profile_db_1=profile_db_1,
            profile_db_2=profile_db_2,
            min_gene_compare_len=1,
            calculate="all",
        )
        .collect()
        .sort("genome")
    )

    float_cols = [col for col, dtype in regular.schema.items() if dtype in {pl.Float32, pl.Float64}]
    if float_cols:
        regular = regular.with_columns([pl.col(c).round(8) for c in float_cols])
        light_all = light_all.with_columns([pl.col(c).round(8) for c in float_cols])

    light_all = light_all.select([pl.col(col).cast(dtype) for col, dtype in regular.schema.items()])
    assert regular.equals(light_all)

    ani_cols = ["genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    regular_ani = regular.select(ani_cols)
    light_ani = (
        light.compare_light_genomes(
            profile_db_1=profile_db_1,
            profile_db_2=profile_db_2,
            min_gene_compare_len=1,
            calculate="ani",
        )
        .collect()
        .sort("genome")
    )
    light_ani = light_ani.with_columns(pl.col("genome_pop_ani").round(8))
    regular_ani = regular_ani.with_columns(pl.col("genome_pop_ani").round(8))
    assert light_ani.columns == ani_cols
    assert regular_ani.equals(light_ani)



def test_light_compare_accepts_direct_duckdb_file_input(tmp_path: Path) -> None:
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1_profile.parquet"
    profile_2 = tmp_path / "p2_profile.parquet"
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)

    profile_db_1 = tmp_path / "p1.light_profile.duckdb"
    profile_db_2 = tmp_path / "p2.light_profile.duckdb"
    light.build_light_profile_bundle_from_parquet(profile_parquet=profile_1, output_file=profile_db_1)
    light.build_light_profile_bundle_from_parquet(profile_parquet=profile_2, output_file=profile_db_2)

    out = light.compare_light_genomes(
        profile_db_1=profile_db_1,
        profile_db_2=profile_db_2,
        min_gene_compare_len=1,
        calculate="ani",
    ).collect()
    assert out.height == 1



def test_light_cli_profile_and_compare(tmp_path: Path) -> None:
    stb_path = tmp_path / "stb.tsv"
    profile_1 = tmp_path / "p1_profile.parquet"
    profile_2 = tmp_path / "p2_profile.parquet"
    _write_stb(stb_path)
    _write_profiles(profile_1, profile_2)
    _write_companion_stats(profile_1)
    _write_companion_stats(profile_2)

    profile_db_1 = tmp_path / "p1.light_profile.duckdb"
    profile_db_2 = tmp_path / "p2.light_profile.duckdb"
    out = tmp_path / "comparison.parquet"

    runner = CliRunner()
    res_1 = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--profile-parquet",
            str(profile_1),
            "--output-file",
            str(profile_db_1),
        ],
    )
    assert res_1.exit_code == 0, res_1.output

    res_2 = runner.invoke(
        light_cli.cli,
        [
            "profile",
            "--profile-parquet",
            str(profile_2),
            "--output-file",
            str(profile_db_2),
        ],
    )
    assert res_2.exit_code == 0, res_2.output

    res_cmp = runner.invoke(
        light_cli.cli,
        [
            "compare",
            "--profile-1",
            str(profile_db_1),
            "--profile-2",
            str(profile_db_2),
            "--output-file",
            str(out),
            "--min-gene-compare-len",
            "1",
            "--calculate",
            "ani",
        ],
    )
    assert res_cmp.exit_code == 0, res_cmp.output
    assert out.exists()

    out_df = pl.read_parquet(out)
    assert out_df.columns == ["genome", "total_positions", "share_allele_pos", "genome_pop_ani", "sample_1", "sample_2"]
