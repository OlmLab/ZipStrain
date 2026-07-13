from click.testing import CliRunner
from zipstrain import cli
from zipstrain import compare as cp
from zipstrain import database as db
from pathlib import Path
import pytest
import polars as pl
import tomllib

a_chr1=[ 1, 0, 1, 0, 0, 0, 4, 4, 2, 1]
t_chr1=[ 0, 2, 3, 1, 0, 0, 1, 0, 0, 1]
c_chr1=[10, 1, 2, 0, 0, 0, 1, 0, 0, 0]
g_chr1=[ 0, 2, 0, 2, 3, 1, 0, 2, 7, 0]
a_chr2=[ 0, 0, 4, 3, 2, 4, 6, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0]
t_chr2=[ 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2]
c_chr2=[ 0, 0, 1, 2, 1, 1, 2, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
g_chr2=[ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
a_chr3=[ 0, 0, 1, 0, 0, 0,15, 0, 0, 0, 0, 0, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_chr3=[ 1, 0, 0, 1, 0,17, 1, 0,10,11, 0, 2, 3, 1, 3, 2,12,10, 1, 6, 1, 0, 0, 1, 0, 0, 8, 5, 0, 1]
c_chr3=[ 0, 0,19,21,18, 1, 2, 10, 0, 0,12, 1, 0, 0, 0,13, 0, 0,11, 0, 1, 0, 1, 2, 1,20, 0, 1, 0, 0]
g_chr3=[20,20, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 2, 3, 1, 1, 2, 7, 0]   
gene_locs = (["NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA","NA"])+ \
            (["NA","NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA",
              "NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA","NA"])+ \
            (["NA","NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA",
              "NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA","NA",
              "NA","NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA"])


def _project_version() -> str:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def _write_profile_db_for_compare_config_tests(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir(parents=True, exist_ok=True)
    profile_paths: dict[str, Path] = {}
    for idx, sample_name in enumerate(("sample_a", "sample_b", "sample_c"), start=1):
        profile_path = profile_dir / f"{sample_name}.parquet"
        pl.DataFrame(
            {
                "chrom": ["chr1"],
                "genome": ["genome1"],
                "pos": [idx],
                "gene": ["NA"],
                "A": [1],
                "T": [0],
                "C": [0],
                "G": [0],
            }
        ).write_parquet(profile_path)
        profile_paths[sample_name] = profile_path

    profile_db = tmp_path / "profile_db.parquet"
    pl.DataFrame(
        {
            "profile_name": list(profile_paths.keys()),
            "profile_location": [str(path) for path in profile_paths.values()],
        }
    ).write_parquet(profile_db)
    return profile_db, profile_paths


def _write_genome_compare_table_for_config_tests(compare_path: Path) -> None:
    pl.DataFrame(
        {
            "genome": ["genome1"],
            "total_positions": [10],
            "share_allele_pos": [10],
            "genome_pop_ani": [100.0],
            "max_consecutive_length": [10],
            "shared_genes_count": [0],
            "identical_gene_count": [0],
            "perc_id_genes": [0.0],
            "sample_1": ["sample_a"],
            "sample_2": ["sample_b"],
        }
    ).write_parquet(compare_path)


def _write_gene_compare_table_for_config_tests(compare_path: Path) -> None:
    pl.DataFrame(
        {
            "genome": ["genome1"],
            "gene": ["gene1"],
            "total_positions": [10],
            "share_allele_pos": [10],
            "ani": [100.0],
            "sample_1": ["sample_a"],
            "sample_2": ["sample_b"],
        }
    ).write_parquet(compare_path)


def _reference_profile_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "genome": ["genome1", "genome1", "genome1", "genome2"],
            "gene": ["gene1", "gene1", "gene1", "gene2"],
            "pos": [1, 2, 3, 1],
            "A": [6, 0, 3, 10],
            "C": [0, 6, 2, 0],
            "G": [0, 0, 0, 0],
            "T": [0, 0, 0, 0],
            cli.ut.REF_BASE_BITMASK_COLUMN: [1, 1, 2, 8],
        }
    )

@pytest.fixture
def profile_1()->pl.LazyFrame:
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20,
        "genome": ["genome1"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2))),
        "gene": gene_locs[:30],
        "A": a_chr1 + a_chr2 ,
        "T": t_chr1 + t_chr2 ,
        "C": c_chr1 + c_chr2 ,
        "G": g_chr1 + g_chr2 ,
    }).lazy()

@pytest.fixture
def profile_2()->pl.LazyFrame:
    """Exactly the same as profile_1"""
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "genome": ["genome1"]*30 + ["genome2"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()

@pytest.fixture
def profile_3()->pl.LazyFrame:
    """Exactly the same as profile_1"""
    return pl.DataFrame({
        "chrom": ["chr3"]*30,
        "genome": ["genome2"]*30,
        "pos":list(range(0,len(a_chr3))),
        "gene": gene_locs[30:],
        "A": a_chr3,
        "T": t_chr3,
        "C": c_chr3,
        "G": g_chr3,
    }).lazy()
@pytest.fixture
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()

def test_cli_top_level_layout():
    commands = set(cli.cli.commands)
    assert {"compare", "map", "profile", "test", "utilities"} == commands


def test_cli_top_level_help_order_and_branding():
    import click as _click

    ctx = _click.Context(cli.cli)
    assert cli.cli.list_commands(ctx) == ["test", "map", "profile", "compare", "utilities"]

    result = CliRunner().invoke(cli.cli, ["-h"])
    assert result.exit_code == 0
    out = result.output
    # Collapse whitespace since the help text wraps across terminal-width lines.
    normalized = " ".join(out.split())
    assert "Parsa Ghadermazi and Matt Olm" in normalized
    assert "University of Colorado Boulder" in normalized
    assert "https://github.com/OlmLab/ZipStrain" in normalized
    # Commands appear in the curated order, not alphabetically.
    assert out.index("\n  test") < out.index("\n  map") < out.index("\n  profile") < out.index("\n  compare")


def test_cli_version_flag():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--version"])
    assert result.exit_code == 0
    assert "zipstrain" in result.output.lower()
    assert _project_version() in result.output


def test_sort_profile_cli_sorts_in_place_and_adds_metadata(tmp_path: Path):
    profile_path = tmp_path / "old_profile.parquet"
    tmp_work_dir = tmp_path / "sort_tmp"
    pl.DataFrame(
        {
            "chrom": ["chr2", "chr1", "chr1"],
            "genome": ["genome2", "genome1", "genome1"],
            "pos": [5, 2, 1],
            "gene": ["gene2", "gene1", "gene1"],
            "A": [0, 5, 4],
            "T": [1, 0, 0],
            "C": [0, 0, 0],
            "G": [0, 0, 0],
        }
    ).write_parquet(profile_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "sort-profile",
            "--input-profile",
            str(profile_path),
            "--tmp-dir",
            str(tmp_work_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    sorted_profile = pl.read_parquet(profile_path)
    assert sorted_profile.to_dicts() == sorted_profile.sort(["chrom", "pos"]).to_dicts()
    assert (
        pl.read_parquet_metadata(profile_path)[cli.pf.PROFILE_SORTED_METADATA_KEY]
        == cli.pf.PROFILE_SORTED_METADATA_VALUE
    )
    assert list(tmp_work_dir.iterdir()) == []

def test_cli_profile_compare(profile_1:pl.LazyFrame,
                             profile_2:pl.LazyFrame,
                             profile_3:pl.LazyFrame,
                             stb:pl.LazyFrame,
                             tmp_path):
    profile_1_dir=tmp_path/"profile_1.parquet"
    profile_2_dir=tmp_path/"profile_2.parquet"
    profile_3_dir=tmp_path/"profile_3.parquet"
    stb_path=tmp_path/"stb.parquet"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    profile_3.sink_parquet(profile_3_dir)
    stb.sink_csv(stb_path,separator="\t",include_header=False)
    runner = CliRunner()
    result = runner.invoke(cli.cli, [
        "utilities",
        "single_compare_genome", 
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_2_dir),
        "--stb-file", str(stb_path),
        "--output-file", str(tmp_path/"output.parquet"),
        "--duckdb-memory-limit", "512MB",
        "--duckdb-temp-directory", str(tmp_path),
    ])
    lf1 = pl.read_parquet(tmp_path/"output.parquet")
    assert result.exit_code == 0 
    assert lf1.shape[0] == 2
    result = runner.invoke(cli.cli, [
        "utilities",
        "single_compare_genome",
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_2_dir),
        "--stb-file", str(stb_path),
        "--engine", "duckdb",
        "--output-file", str(tmp_path/"output_duckdb.parquet"),
    ])
    lf1_duckdb = pl.read_parquet(tmp_path/"output_duckdb.parquet")
    assert result.exit_code == 0
    # The two engines agree, but derived percentages (e.g. genome_con_ani) can
    # differ by ~1 ULP due to float non-associativity, so compare non-float
    # columns exactly and float columns with tolerance.
    a = lf1.sort("genome")
    b = lf1_duckdb.sort("genome")
    assert a.columns == b.columns
    for col, dtype in a.schema.items():
        if dtype.is_float():
            assert a[col].to_list() == pytest.approx(b[col].to_list(), abs=1e-9, nan_ok=True)
        else:
            assert a[col].equals(b[col])
    result = runner.invoke(cli.cli, [
        "utilities",
        "single_compare_genome", 
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_3_dir),
        "--stb-file", str(stb_path),
        "--output-file", str(tmp_path/"output.parquet"),
    ])
    lf2 = pl.read_parquet(tmp_path/"output.parquet")
    assert result.exit_code == 0
    assert lf2.shape[0] == 2
    


def test_single_compare_genome_duckdb_scope_skips_prefilter(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, monkeypatch):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    stb_path = tmp_path / "stb.tsv"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    def _should_not_prefilter(*args, **kwargs):
        raise AssertionError("duckdb_prefilter_by_scope should not be called for --engine duckdb")

    monkeypatch.setattr(cli.cp, "duckdb_prefilter_by_scope", _should_not_prefilter)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--stb-file",
            str(stb_path),
            "--engine",
            "duckdb",
            "--genome",
            "genome1",
            "--output-file",
            str(tmp_path / "scoped_duckdb.parquet"),
        ],
    )
    assert result.exit_code == 0
    assert (tmp_path / "scoped_duckdb.parquet").exists()


def test_single_compare_gene_duckdb_scope_skips_prefilter(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, monkeypatch):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    stb_path = tmp_path / "stb.tsv"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    def _should_not_prefilter(*args, **kwargs):
        raise AssertionError("duckdb_prefilter_by_scope should not be called for --engine duckdb")

    monkeypatch.setattr(cli.cp, "duckdb_prefilter_by_scope", _should_not_prefilter)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_gene",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--stb-file",
            str(stb_path),
            "--engine",
            "duckdb",
            "--scope",
            "genome1:gene1",
            "--output-file",
            str(tmp_path / "scoped_gene_duckdb.parquet"),
        ],
    )
    assert result.exit_code == 0
    assert (tmp_path / "scoped_gene_duckdb.parquet").exists()


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_single_compare_gene_without_stb_succeeds(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, tmp_path, engine):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    output_path = tmp_path / f"gene_no_stb_{engine}.parquet"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_gene",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--engine",
            engine,
            "--scope",
            "genome1:gene1",
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert output_path.exists()


def test_single_compare_genome_polars_scope_uses_polars_prefilter(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, monkeypatch):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    stb_path = tmp_path / "stb.tsv"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    calls = {"polars": 0}
    original_polars_prefilter = cp.polars_prefilter_by_scope

    def _should_not_use_duckdb(*args, **kwargs):
        raise AssertionError("duckdb_prefilter_by_scope should not be called for --engine polars")

    def _polars_prefilter(*args, **kwargs):
        calls["polars"] += 1
        return original_polars_prefilter(*args, **kwargs)

    monkeypatch.setattr(cli.cp, "duckdb_prefilter_by_scope", _should_not_use_duckdb)
    monkeypatch.setattr(cli.cp, "polars_prefilter_by_scope", _polars_prefilter)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--stb-file",
            str(stb_path),
            "--engine",
            "polars",
            "--genome",
            "genome1",
            "--output-file",
            str(tmp_path / "scoped_polars.parquet"),
        ],
    )
    assert result.exit_code == 0
    assert calls["polars"] == 1


def test_single_compare_gene_polars_scope_uses_polars_prefilter(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, monkeypatch):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    stb_path = tmp_path / "stb.tsv"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    calls = {"polars": 0}
    original_polars_prefilter = cp.polars_prefilter_by_scope

    def _should_not_use_duckdb(*args, **kwargs):
        raise AssertionError("duckdb_prefilter_by_scope should not be called for --engine polars")

    def _polars_prefilter(*args, **kwargs):
        calls["polars"] += 1
        return original_polars_prefilter(*args, **kwargs)

    monkeypatch.setattr(cli.cp, "duckdb_prefilter_by_scope", _should_not_use_duckdb)
    monkeypatch.setattr(cli.cp, "polars_prefilter_by_scope", _polars_prefilter)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_gene",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--stb-file",
            str(stb_path),
            "--engine",
            "polars",
            "--scope",
            "genome1:gene1",
            "--output-file",
            str(tmp_path / "scoped_gene_polars.parquet"),
        ],
    )
    assert result.exit_code == 0
    assert calls["polars"] == 1


def test_single_compare_genome_calculate_controls_output_columns(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path):
    profile_1_dir = tmp_path / "profile_1.parquet"
    profile_2_dir = tmp_path / "profile_2.parquet"
    stb_path = tmp_path / "stb.tsv"
    out_path = tmp_path / "calc_cols.parquet"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_2_dir),
            "--stb-file",
            str(stb_path),
            "--calculate",
            "ibs",
            "--output-file",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    out = pl.read_parquet(out_path)
    assert set(out.columns) == {"genome", "max_consecutive_length", "sample_1", "sample_2"}


def test_single_compare_genome_without_stb_reports_only_nonzero_genomes_polars(profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, tmp_path):
    profile_1_dir = tmp_path / "profile_2.parquet"
    profile_3_dir = tmp_path / "profile_3.parquet"
    out_path = tmp_path / "no_stb_polars.parquet"
    profile_2.sink_parquet(profile_1_dir)
    profile_3.sink_parquet(profile_3_dir)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_3_dir),
            "--engine",
            "polars",
            "--output-file",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    out = pl.read_parquet(out_path).sort("genome")
    assert out.get_column("genome").to_list() == ["genome2"]
    assert out.get_column("total_positions").to_list() == [22]


def test_single_compare_genome_without_stb_reports_only_nonzero_genomes_duckdb(profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, tmp_path):
    profile_1_dir = tmp_path / "profile_2.parquet"
    profile_3_dir = tmp_path / "profile_3.parquet"
    out_path = tmp_path / "no_stb_duckdb.parquet"
    profile_2.sink_parquet(profile_1_dir)
    profile_3.sink_parquet(profile_3_dir)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_dir),
            "--mpileup-contig-2",
            str(profile_3_dir),
            "--engine",
            "duckdb",
            "--output-file",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    out = pl.read_parquet(out_path).sort("genome")
    assert out.get_column("genome").to_list() == ["genome2"]
    assert out.get_column("total_positions").to_list() == [22]


def test_generate_genome_pairs_command(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, tmp_path):
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir()
    profile_1.sink_parquet(profile_dir / "profile_1_profile.parquet")
    profile_2.sink_parquet(profile_dir / "profile_2_profile.parquet")
    profile_3.sink_parquet(profile_dir / "profile_3_profile.parquet")
    pl.DataFrame({"cov": [1]}).write_parquet(profile_dir / "ignore_me.parquet")

    runner = CliRunner()
    output_file = tmp_path / "pairs.parquet"
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "generate-genome-pairs",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(output_file),
            "--write-batch-size",
            "2",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "profiles=3" in result.output
    assert "pairs=3" in result.output
    pairs = pl.read_parquet(output_file).sort(["sample_name_1", "sample_name_2"])
    assert pairs.columns == ["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"]
    assert pairs.height == 3
    assert set(zip(pairs["sample_name_1"].to_list(), pairs["sample_name_2"].to_list())) == {
        ("profile_1", "profile_2"),
        ("profile_1", "profile_3"),
        ("profile_2", "profile_3"),
    }


def test_reference_ani_is_added_to_gene_and_genome_stats():
    profile_lf = _reference_profile_frame().lazy()
    gene_bed = pl.DataFrame(
        {
            "gene": ["gene1", "gene2"],
            "scaffold": ["chr1", "chr2"],
            "start": [1, 1],
            "end": [3, 1],
        }
    ).lazy()
    stb = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2"],
            "genome": ["genome1", "genome2"],
        }
    ).lazy()
    genome_bed = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [3, 1],
        }
    ).lazy()
    read_locs = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr1", "chr1", "chr2"],
            "loc": [1, 2, 3, 1],
        }
    ).lazy()

    gene_stats = cli.ut.get_gene_stats(
        profile=profile_lf,
        gene_bed=gene_bed,
        stb=stb,
    ).collect().sort(["genome", "gene"])
    assert gene_stats.columns == ["genome", "gene", "length", "breadth", "coverage", "ref_ani"]
    assert gene_stats.to_dicts() == [
        {
            "genome": "genome1",
            "gene": "gene1",
            "length": 3,
            "breadth": 1.0,
            "coverage": pytest.approx((6 + 6 + 5) / 3),
            "ref_ani": pytest.approx(66.66666666666666),
        },
        {
            "genome": "genome2",
            "gene": "gene2",
            "length": 1,
            "breadth": 1.0,
            "coverage": 10.0,
            "ref_ani": pytest.approx(0.0),
        },
    ]

    genome_stats = cli.ut.get_genome_stats(
        profile=profile_lf,
        bed=genome_bed,
        stb=stb,
        read_loc_table=read_locs,
    ).collect().sort("genome")
    assert "ref_ani" in genome_stats.columns
    by_genome = genome_stats.rows_by_key("genome", unique=True, named=True)
    assert by_genome["genome1"]["ref_ani"] == pytest.approx(66.66666666666666)
    assert by_genome["genome2"]["ref_ani"] == pytest.approx(0.0)


def test_reference_ani_is_omitted_from_stats_without_reference_base():
    profile_lf = _reference_profile_frame().drop(cli.ut.REF_BASE_BITMASK_COLUMN).lazy()
    gene_bed = pl.DataFrame(
        {
            "gene": ["gene1", "gene2"],
            "scaffold": ["chr1", "chr2"],
            "start": [1, 1],
            "end": [3, 1],
        }
    ).lazy()
    stb = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2"],
            "genome": ["genome1", "genome2"],
        }
    ).lazy()
    genome_bed = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [3, 1],
        }
    ).lazy()
    read_locs = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr1", "chr1", "chr2"],
            "loc": [1, 2, 3, 1],
        }
    ).lazy()

    gene_stats = cli.ut.get_gene_stats(
        profile=profile_lf,
        gene_bed=gene_bed,
        stb=stb,
    ).collect()
    genome_stats = cli.ut.get_genome_stats(
        profile=profile_lf,
        bed=genome_bed,
        stb=stb,
        read_loc_table=read_locs,
    ).collect()
    assert "ref_ani" not in gene_stats.columns
    assert "ref_ani" not in genome_stats.columns


def test_get_snp_reference_cli_outputs_only_reference_snps(tmp_path):
    profile_path = tmp_path / "profile_with_ref.parquet"
    output_path = tmp_path / "reference_snps.parquet"
    _reference_profile_frame().write_parquet(profile_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "get-snp-reference",
            "--profile-file",
            str(profile_path),
            "--min-cov",
            "5",
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    out = pl.read_parquet(output_path).sort(["chrom", "pos"])
    # Full classified SNV table: monomorphic reference matches (chr1 pos1) are
    # omitted; chr1 pos2 / chr2 pos1 are fixed substitutions (SNS); chr1 pos3 is
    # a consensus SNP that still carries the reference allele (con_SNV).
    assert out.select(["chrom", "pos", "ref_base", "con_base", "class"]).to_dicts() == [
        {"chrom": "chr1", "pos": 2, "ref_base": "A", "con_base": "C", "class": "SNS"},
        {"chrom": "chr1", "pos": 3, "ref_base": "C", "con_base": "A", "class": "con_SNV"},
        {"chrom": "chr2", "pos": 1, "ref_base": "T", "con_base": "A", "class": "SNS"},
    ]
    row = out.filter(pl.col("pos") == 3).row(0, named=True)
    assert row["allele_count"] == 2
    assert row["position_coverage"] == 5
    assert row["ref_freq"] == pytest.approx(0.4)
    assert row["con_freq"] == pytest.approx(0.6)


def test_get_snp_reference_cli_respects_min_cov(tmp_path):
    profile_path = tmp_path / "profile_with_ref.parquet"
    output_path = tmp_path / "reference_snps_min_cov.parquet"
    _reference_profile_frame().write_parquet(profile_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "get-snp-reference",
            "--profile-file",
            str(profile_path),
            "--min-cov",
            "7",
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    out = pl.read_parquet(output_path).sort(["chrom", "pos"])
    # min_cov=7 drops chr1 pos2 (cov 6) and pos3 (cov 5); only chr2 pos1 remains.
    assert out.select(["chrom", "pos", "class"]).to_dicts() == [
        {"chrom": "chr2", "pos": 1, "class": "SNS"},
    ]


def test_get_snp_reference_cli_requires_reference_base_column(tmp_path):
    profile_path = tmp_path / "profile_without_ref.parquet"
    output_path = tmp_path / "reference_snps_missing_ref.parquet"
    _reference_profile_frame().drop(cli.ut.REF_BASE_BITMASK_COLUMN).write_parquet(profile_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "get-snp-reference",
            "--profile-file",
            str(profile_path),
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code != 0
    assert "ref_base_bitmask" in result.output
    assert "--reference-fasta" in result.output


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_single_compare_genome_strips_profile_suffix_from_sample_columns(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, engine):
    profile_1_path = tmp_path / "sample_alpha_profile.parquet"
    profile_2_path = tmp_path / "sample_beta_profile.parquet"
    output_path = tmp_path / f"sample_suffix_{engine}.parquet"
    stb_path = tmp_path / "stb.tsv"

    profile_1.sink_parquet(profile_1_path)
    profile_2.sink_parquet(profile_2_path)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_path),
            "--mpileup-contig-2",
            str(profile_2_path),
            "--stb-file",
            str(stb_path),
            "--engine",
            engine,
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    out = pl.read_parquet(output_path)
    assert set(out.get_column("sample_1").unique().to_list()) == {"sample_alpha"}
    assert set(out.get_column("sample_2").unique().to_list()) == {"sample_beta"}


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_single_compare_genome_writes_mismatch_tolerant_metadata(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, engine):
    profile_1_path = tmp_path / f"metadata_left_{engine}.parquet"
    profile_2_path = tmp_path / f"metadata_right_{engine}.parquet"
    output_path = tmp_path / f"metadata_compare_{engine}.parquet"
    stb_path = tmp_path / "stb.tsv"

    profile_1.sink_parquet(
        profile_1_path,
        metadata={
            cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_hash_shared",
            cli.ut.PROFILE_GENE_HASH_METADATA_KEY: "gene_hash_left",
            cli.ut.PROFILE_NULL_MODEL_HASH_METADATA_KEY: "null_hash_shared",
            cli.ut.PROFILE_STB_HASH_METADATA_KEY: "stb_hash_shared",
        },
    )
    profile_2.sink_parquet(
        profile_2_path,
        metadata={
            cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_hash_shared",
            cli.ut.PROFILE_GENE_HASH_METADATA_KEY: "gene_hash_right",
            cli.ut.PROFILE_NULL_MODEL_HASH_METADATA_KEY: "null_hash_shared",
            cli.ut.PROFILE_STB_HASH_METADATA_KEY: "stb_hash_shared",
        },
    )
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_genome",
            "--mpileup-contig-1",
            str(profile_1_path),
            "--mpileup-contig-2",
            str(profile_2_path),
            "--stb-file",
            str(stb_path),
            "--engine",
            engine,
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    metadata = pl.read_parquet_metadata(output_path)
    assert metadata[cli.ut.COMPARE_KIND_METADATA_KEY] == "genome"
    assert metadata[cli.ut.COMPARE_SCOPE_METADATA_KEY] == "all"
    assert metadata[cli.ut.COMPARE_MIN_COV_METADATA_KEY] == "5"
    assert metadata[cli.ut.COMPARE_MIN_GENE_COMPARE_LEN_METADATA_KEY] == "100"
    assert metadata[cli.ut.COMPARE_ENGINE_METADATA_KEY] == engine
    assert metadata[cli.ut.COMPARE_USES_STB_METADATA_KEY] == "1"
    assert metadata[cli.ut.COMPARE_REFERENCE_HASH_METADATA_KEY] == "ref_hash_shared"
    assert metadata[cli.ut.COMPARE_GENE_HASH_METADATA_KEY] == cli.ut.COMPARE_METADATA_MISSING_VALUE
    assert metadata[cli.ut.COMPARE_NULL_MODEL_HASH_METADATA_KEY] == "null_hash_shared"
    assert metadata[cli.ut.COMPARE_STB_HASH_METADATA_KEY] == "stb_hash_shared"


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_single_compare_gene_writes_scope_metadata(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, tmp_path, engine):
    profile_1_path = tmp_path / f"gene_metadata_left_{engine}.parquet"
    profile_2_path = tmp_path / f"gene_metadata_right_{engine}.parquet"
    output_path = tmp_path / f"gene_metadata_compare_{engine}.parquet"

    profile_1.sink_parquet(
        profile_1_path,
        metadata={
            cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_hash_shared",
            cli.ut.PROFILE_GENE_HASH_METADATA_KEY: "gene_hash_shared",
            cli.ut.PROFILE_NULL_MODEL_HASH_METADATA_KEY: "null_hash_left",
        },
    )
    profile_2.sink_parquet(
        profile_2_path,
        metadata={
            cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_hash_shared",
            cli.ut.PROFILE_GENE_HASH_METADATA_KEY: "gene_hash_shared",
            cli.ut.PROFILE_NULL_MODEL_HASH_METADATA_KEY: "null_hash_right",
        },
    )

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "single_compare_gene",
            "--mpileup-contig-1",
            str(profile_1_path),
            "--mpileup-contig-2",
            str(profile_2_path),
            "--engine",
            engine,
            "--scope",
            "genome1:gene1",
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    metadata = pl.read_parquet_metadata(output_path)
    assert metadata[cli.ut.COMPARE_KIND_METADATA_KEY] == "gene"
    assert metadata[cli.ut.COMPARE_SCOPE_METADATA_KEY] == "genome1:gene1"
    assert metadata[cli.ut.COMPARE_ENGINE_METADATA_KEY] == engine
    assert metadata[cli.ut.COMPARE_USES_STB_METADATA_KEY] == "0"
    assert metadata[cli.ut.COMPARE_REFERENCE_HASH_METADATA_KEY] == "ref_hash_shared"
    assert metadata[cli.ut.COMPARE_GENE_HASH_METADATA_KEY] == "gene_hash_shared"
    assert metadata[cli.ut.COMPARE_NULL_MODEL_HASH_METADATA_KEY] == cli.ut.COMPARE_METADATA_MISSING_VALUE
    assert metadata[cli.ut.COMPARE_STB_HASH_METADATA_KEY] == cli.ut.COMPARE_METADATA_MISSING_VALUE


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_chunk_genome_compare_command(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, stb: pl.LazyFrame, tmp_path, engine):
    profile_1_path = tmp_path / "profile_1.parquet"
    profile_2_path = tmp_path / "profile_2.parquet"
    profile_3_path = tmp_path / "profile_3.parquet"
    stb_path = tmp_path / "stb.tsv"
    output_file = tmp_path / f"chunk_compare_{engine}.parquet"

    profile_1.sink_parquet(profile_1_path)
    profile_2.sink_parquet(profile_2_path)
    profile_3.sink_parquet(profile_3_path)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    pair_table = pl.DataFrame(
        {
            "sample_name_1": ["profile_1", "profile_1"],
            "sample_name_2": ["profile_2", "profile_3"],
            "profile_location_1": [str(profile_1_path), str(profile_1_path)],
            "profile_location_2": [str(profile_2_path), str(profile_3_path)],
        }
    )
    pair_table_path = tmp_path / "pairs.parquet"
    pair_table.write_parquet(pair_table_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "chunk-genome-compare",
            "--pair-table",
            str(pair_table_path),
            "--stb-file",
            str(stb_path),
            "--output-file",
            str(output_file),
            "--workers",
            "1",
            "--calculate",
            "ani",
            "--engine",
            engine,
        ],
    )

    assert result.exit_code == 0, result.output
    assert "pairs=2" in result.output
    assert "avg_wall_s_per_pair=" in result.output

    actual = pl.read_parquet(output_file).sort(["sample_1", "sample_2", "genome"])
    expected_frames = []
    for sample_1, sample_2, left, right in [
        ("profile_1", "profile_2", profile_1_path, profile_2_path),
        ("profile_1", "profile_3", profile_1_path, profile_3_path),
    ]:
        expected_frames.append(
            cp.compare_genomes(
                mpile_contig_1=left,
                mpile_contig_2=right,
                min_cov=5,
                min_gene_compare_len=100,
                genome_scope="all",
                ani_method="popani",
                engine=engine,
                stb_file=stb_path,
                calculate="ani",
            )
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select(["genome", "total_positions", "share_allele_pos", "genome_pop_ani", "sample_1", "sample_2"])
            .collect(engine="streaming")
        )
    expected = pl.concat(expected_frames, how="vertical_relaxed").sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_chunk_genome_compare_without_stb_reports_only_nonzero_genomes(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, tmp_path, engine):
    profile_1_path = tmp_path / "profile_1.parquet"
    profile_2_path = tmp_path / "profile_2.parquet"
    profile_3_path = tmp_path / "profile_3.parquet"
    output_file = tmp_path / f"chunk_compare_no_stb_{engine}.parquet"

    profile_1.sink_parquet(profile_1_path)
    profile_2.sink_parquet(profile_2_path)
    profile_3.sink_parquet(profile_3_path)

    pair_table = pl.DataFrame(
        {
            "sample_name_1": ["profile_1", "profile_1"],
            "sample_name_2": ["profile_2", "profile_3"],
            "profile_location_1": [str(profile_1_path), str(profile_1_path)],
            "profile_location_2": [str(profile_2_path), str(profile_3_path)],
        }
    )
    pair_table_path = tmp_path / "pairs_no_stb.parquet"
    pair_table.write_parquet(pair_table_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "chunk-genome-compare",
            "--pair-table",
            str(pair_table_path),
            "--output-file",
            str(output_file),
            "--workers",
            "1",
            "--calculate",
            "ani",
            "--engine",
            engine,
        ],
    )

    assert result.exit_code == 0, result.output
    actual = pl.read_parquet(output_file).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2, left, right in [
        ("profile_1", "profile_2", profile_1_path, profile_2_path),
        ("profile_1", "profile_3", profile_1_path, profile_3_path),
    ]:
        expected_frames.append(
            cp.compare_genomes(
                mpile_contig_1=left,
                mpile_contig_2=right,
                min_cov=5,
                min_gene_compare_len=100,
                genome_scope="all",
                ani_method="popani",
                engine=engine,
                stb_file=None,
                calculate="ani",
            )
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select(["genome", "total_positions", "share_allele_pos", "genome_pop_ani", "sample_1", "sample_2"])
            .collect(engine="streaming")
        )
    expected = pl.concat(expected_frames, how="vertical_relaxed").sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)
    assert actual.filter(
        (pl.col("sample_1") == "profile_1")
        & (pl.col("sample_2") == "profile_3")
    ).height == 0


def test_compare_genomes_batch_passes_duckdb_threads(tmp_path, monkeypatch):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    captured = {}

    def _fake_lazy_run_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_compares", _fake_lazy_run_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "--profile-db",
            str(profile_db),
            "--run-dir",
            str(tmp_path / "run"),
            "--scope",
            "all",
            "--ani-method",
            "conani",
            "--calculate",
            "ani",
            "--engine",
            "duckdb",
            "--duckdb-threads",
            "7",
        ],
    )
    assert result.exit_code == 0
    assert captured["duckdb_threads"] == 7
    assert captured["compare_engine"] == "duckdb"
    assert captured["calculate"] == "ani"


def test_compare_genomes_uses_default_versioned_docker_image(tmp_path, monkeypatch):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    captured = {}

    def _fake_lazy_run_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_compares", _fake_lazy_run_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "--profile-db",
            str(profile_db),
            "--run-dir",
            str(tmp_path / "run"),
            "--container-engine",
            "docker",
        ],
    )
    assert result.exit_code == 0
    assert isinstance(captured["container_engine"], cli.tm.DockerEngine)
    assert captured["container_engine"].address == f"parsaghadermazi/zipstrain:{_project_version()}"


def test_compare_genomes_honors_container_address_override(tmp_path, monkeypatch):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    captured = {}

    def _fake_lazy_run_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_compares", _fake_lazy_run_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "--profile-db",
            str(profile_db),
            "--run-dir",
            str(tmp_path / "run"),
            "--container-engine",
            "docker",
            "--container-address",
            "parsaghadermazi/zipstrain:custom-tag",
        ],
    )
    assert result.exit_code == 0
    assert isinstance(captured["container_engine"], cli.tm.DockerEngine)
    assert captured["container_engine"].address == "parsaghadermazi/zipstrain:custom-tag"


def test_compare_genes_batch_passes_duckdb_threads(tmp_path, monkeypatch):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    captured = {}

    def _fake_lazy_run_gene_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_gene_compares", _fake_lazy_run_gene_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "--compare-genes",
            "--profile-db",
            str(profile_db),
            "--run-dir",
            str(tmp_path / "run"),
            "--scope",
            "all:all",
            "--engine",
            "duckdb",
            "--duckdb-threads",
            "9",
        ],
    )
    assert result.exit_code == 0
    assert captured["duckdb_threads"] == 9
    assert captured["compare_engine"] == "duckdb"


def test_compare_genes_honors_apptainer_container_address_override(tmp_path, monkeypatch):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    captured = {}

    def _fake_lazy_run_gene_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_gene_compares", _fake_lazy_run_gene_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "--compare-genes",
            "--profile-db",
            str(profile_db),
            "--run-dir",
            str(tmp_path / "run"),
            "--container-engine",
            "apptainer",
            "--container-address",
            "/scratch/containers/zipstrain-0.10.1.img",
        ],
    )
    assert result.exit_code == 0
    assert isinstance(captured["container_engine"], cli.tm.ApptainerEngine)
    assert captured["container_engine"].address == "/scratch/containers/zipstrain-0.10.1.img"


def _write_profiles_csv_for_compare(tmp_path):
    """A profiles CSV (profile_name,profile_location) pointing at real profile parquets."""
    _profile_db_parquet, profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    csv_path = tmp_path / "profiles.csv"
    lines = ["profile_name,profile_location"]
    lines += [f"{name},{path}" for name, path in profile_paths.items()]
    csv_path.write_text("\n".join(lines) + "\n")
    return csv_path


def test_compare_accepts_profiles_csv_without_build_profile_db(tmp_path, monkeypatch):
    csv_path = _write_profiles_csv_for_compare(tmp_path)
    captured = {}
    monkeypatch.setattr(cli.tm, "lazy_run_compares", lambda **kwargs: captured.update(kwargs))

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        ["compare", "--profile-db", str(csv_path), "--run-dir", str(tmp_path / "run")],
    )
    assert result.exit_code == 0, result.output
    # Routed to the genome path by default, with the CSV loaded into a profile db.
    assert isinstance(captured["comps_db"], cli.db.GenomeComparisonDatabase)


def test_compare_genes_flag_routes_to_gene_compare(tmp_path, monkeypatch):
    csv_path = _write_profiles_csv_for_compare(tmp_path)
    captured = {}
    monkeypatch.setattr(cli.tm, "lazy_run_gene_compares", lambda **kwargs: captured.update(kwargs))
    # Genome path must NOT be called in gene mode.
    monkeypatch.setattr(
        cli.tm, "lazy_run_compares",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("genome path called in gene mode")),
    )

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        ["compare", "--compare-genes", "--profile-db", str(csv_path), "--run-dir", str(tmp_path / "run")],
    )
    assert result.exit_code == 0, result.output
    assert isinstance(captured["comps_db"], cli.db.GeneComparisonDatabase)


def test_compare_method_matrix_routes_to_matrix_workflow(tmp_path, monkeypatch):
    csv_path = _write_profiles_csv_for_compare(tmp_path)
    captured = {}

    def _fake_run_matrix_compare(**kwargs):
        captured.update(kwargs)
        out = Path(kwargs["run_dir"]) / "all_comparisons.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame({"genome": ["g"], "genome_pop_ani": [100.0]}).write_parquet(out)
        return out

    monkeypatch.setattr(cli.matrix_workflow, "run_matrix_compare", _fake_run_matrix_compare)
    # Standard path must not run in matrix mode.
    monkeypatch.setattr(
        cli.tm, "lazy_run_compares",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("standard path called in matrix mode")),
    )

    result = CliRunner().invoke(
        cli.cli,
        [
            "compare", "--method", "matrix",
            "--profile-db", str(csv_path),
            "--stb-file", str(tmp_path / "ref.stb"),
            "--backend", "numpy",
            "--run-dir", str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0, result.output
    # Profiles were passed as (name, location) pairs, and matrix knobs propagated.
    assert [name for name, _loc in captured["profiles"]] == ["sample_a", "sample_b", "sample_c"]
    assert captured["backend"] == "numpy"
    assert captured["compare_genes"] is False


def test_compare_default_scope_differs_by_mode(tmp_path, monkeypatch):
    csv_path = _write_profiles_csv_for_compare(tmp_path)

    genome_db = {}
    monkeypatch.setattr(cli.tm, "lazy_run_compares", lambda **kwargs: genome_db.update(comps_db=kwargs["comps_db"]))
    runner = CliRunner()
    assert runner.invoke(cli.cli, ["compare", "--profile-db", str(csv_path), "--run-dir", str(tmp_path / "g")]).exit_code == 0
    assert genome_db["comps_db"].config.scope == "all"

    gene_db = {}
    monkeypatch.setattr(cli.tm, "lazy_run_gene_compares", lambda **kwargs: gene_db.update(comps_db=kwargs["comps_db"]))
    assert runner.invoke(cli.cli, ["compare", "--compare-genes", "--profile-db", str(csv_path), "--run-dir", str(tmp_path / "gene")]).exit_code == 0
    assert gene_db["comps_db"].config.scope == "all:all"


def test_profile_command_calls_lazy_run_profile(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov":list(range(100)),"max_error_count":[0 for _ in range(100)]}).write_parquet(null_model)
    
    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--reference-fasta",
            str(reference_fasta),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--gene-range-table",
            str(tmp_path / "gene_range.tsv"),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0
    assert captured["run_dir"] == tmp_path / "run"
    assert captured["reference_fasta_file"] == reference_fasta
    assert captured["profiling_contract_file"] is None
    assert captured["min_mapq"] == cli.pf.PROFILE_MIN_MAPQ_DEFAULT
    assert captured["min_baseq"] == cli.pf.PROFILE_MIN_BASEQ_DEFAULT
    assert captured["min_read_ani"] == cli.pf.PROFILE_MIN_READ_ANI_DEFAULT
    assert captured["read_inclusion"] == cli.pf.PROFILE_READ_INCLUSION_DEFAULT


def test_profile_command_passes_custom_read_filters(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": list(range(10)), "max_error_count": [0 for _ in range(10)]}).write_parquet(null_model)

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
            "--min-mapq",
            "7",
            "--min-baseq",
            "21",
            "--min-read-ani",
            "0.97",
            "--read-inclusion",
            "proper-pairs",
        ],
    )
    assert result.exit_code == 0
    assert captured["min_mapq"] == 7
    assert captured["min_baseq"] == 21
    assert captured["min_read_ani"] == 0.97
    assert captured["read_inclusion"] == "proper-pairs"


def test_profile_command_auto_generates_assets_from_minimal_inputs(tmp_path, monkeypatch):
    """profile with just input-table + reference + stb builds its own assets."""
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference = tmp_path / "reference.fna"
    reference.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTTGGGGCC\n")
    stb = tmp_path / "reference.stb"
    stb.write_text("chr1\tgenome1\nchr2\tgenome1\n")
    run_dir = tmp_path / "run"

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table", str(input_table),
            "--reference-fasta", str(reference),
            "--stb-file", str(stb),
            "--run-dir", str(run_dir),
            "--max-total-reads", "100",
        ],
    )
    assert result.exit_code == 0, result.output

    assets_dir = run_dir / cli.pf.DEFAULT_PROFILING_ASSETS_DIRNAME
    assert assets_dir.is_dir()
    # The paths handed to profiling point into the generated assets directory
    # and the files actually exist.
    assert captured["null_model_file"] == assets_dir / cli.pf.ASSET_NULL_MODEL_FILENAME
    assert captured["bed_file"] == assets_dir / cli.pf.ASSET_BED_FILENAME
    assert captured["genome_length_file"] == assets_dir / cli.pf.ASSET_GENOME_LENGTH_FILENAME
    assert captured["null_model_file"].exists()
    assert captured["bed_file"].exists()
    assert captured["genome_length_file"].exists()


def test_profile_command_errors_without_reference_when_autogen_needed(tmp_path, monkeypatch):
    """Missing reference AND missing bed/genome-length is a usage error."""
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    stb = tmp_path / "reference.stb"
    stb.write_text("chr1\tgenome1\n")

    monkeypatch.setattr(cli.tm, "lazy_run_profile", lambda **kwargs: None)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table", str(input_table),
            "--stb-file", str(stb),
            "--run-dir", str(tmp_path / "run"),
        ],
    )
    assert result.exit_code != 0
    assert "reference-fasta is required" in str(result.exception)


def test_profile_command_uses_default_versioned_apptainer_image(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": list(range(100)), "max_error_count": [0 for _ in range(100)]}).write_parquet(null_model)

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--reference-fasta",
            str(reference_fasta),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
            "--container-engine",
            "apptainer",
        ],
    )
    assert result.exit_code == 0
    assert isinstance(captured["container_engine"], cli.tm.ApptainerEngine)
    assert captured["container_engine"].address == f"docker://parsaghadermazi/zipstrain:{_project_version()}"


def test_profile_command_honors_container_address_override(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": list(range(100)), "max_error_count": [0 for _ in range(100)]}).write_parquet(null_model)

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--reference-fasta",
            str(reference_fasta),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
            "--container-engine",
            "apptainer",
            "--container-address",
            "/scratch/containers/zipstrain-0.10.1.sif",
        ],
    )
    assert result.exit_code == 0
    assert isinstance(captured["container_engine"], cli.tm.ApptainerEngine)
    assert captured["container_engine"].address == "/scratch/containers/zipstrain-0.10.1.sif"


def test_profile_command_allows_missing_gene_range_table(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": list(range(10)), "max_error_count": [0 for _ in range(10)]}).write_parquet(null_model)

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--reference-fasta",
            str(reference_fasta),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0
    assert captured["gene_range_table"] is None


def test_profile_command_allows_missing_reference_fasta(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    null_model = tmp_path / "null_model.parquet"
    pl.DataFrame({"cov": list(range(10)), "max_error_count": [0 for _ in range(10)]}).write_parquet(null_model)

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0
    assert captured["reference_fasta_file"] is None


def test_prepare_profiling_creates_null_model_and_contract(tmp_path):
    reference_fasta = tmp_path / "reference.fna"
    gene_fasta = tmp_path / "genes.fna"
    stb_file = tmp_path / "mapping.stb"
    output_dir = tmp_path / "profiling_assets"
    reference_fasta.write_text(">chr1\nACGTACGT\n>chr2\nTTTTAAAA\n")
    gene_fasta.write_text(
        ">chr1_1 # 1 # 8 # 1 # ID=1_1;partial=00\nACGTACGT\n"
        ">chr2_1 # 1 # 8 # 1 # ID=2_1;partial=00\nTTTTAAAA\n"
    )
    stb_file.write_text("chr1\tgenome1\nchr2\tgenome2\n")

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "prepare_profiling",
            "--reference-fasta",
            str(reference_fasta),
            "--gene-fasta",
            str(gene_fasta),
            "--stb-file",
            str(stb_file),
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "genomes_bed_file.bed").exists()
    assert (output_dir / "gene_range_table.tsv").exists()
    assert (output_dir / "genome_lengths.parquet").exists()
    assert (output_dir / "null_model.parquet").exists()
    assert (output_dir / "profiling_contract.json").exists()
    assert (output_dir / "reference.fasta").read_text() == reference_fasta.read_text()

    contract = cli.ut.read_profile_contract_file(output_dir / "profiling_contract.json")
    assert contract == {
        "reference_hash": cli.ut.sha256_file(reference_fasta),
        "gene_hash": cli.ut.sha256_file(gene_fasta),
        "null_model_hash": cli.ut.sha256_file(output_dir / "null_model.parquet"),
        "stb_hash": cli.ut.sha256_file(stb_file),
    }

    null_model = pl.read_parquet(output_dir / "null_model.parquet")
    assert null_model.columns == ["cov", "max_error_count"]
    assert null_model.height > 0


def test_profile_command_passes_profiling_contract(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    null_model = tmp_path / "null_model.parquet"
    profiling_contract = tmp_path / "profiling_contract.json"
    pl.DataFrame({"cov": list(range(10)), "max_error_count": [0 for _ in range(10)]}).write_parquet(null_model)
    profiling_contract.write_text("{}\n")

    captured = {}

    def _fake_lazy_run_profile(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_profile", _fake_lazy_run_profile)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "profile",
            "--input-table",
            str(input_table),
            "--reference-fasta",
            str(reference_fasta),
            "--stb-file",
            str(tmp_path / "stb.tsv"),
            "--null-model",
            str(null_model),
            "--profiling-contract",
            str(profiling_contract),
            "--bed-file",
            str(tmp_path / "bed.bed"),
            "--genome-length-file",
            str(tmp_path / "genome_length.parquet"),
            "--run-dir",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0
    assert captured["profiling_contract_file"] == profiling_contract
    assert captured["reference_fasta_file"] == reference_fasta


def test_profile_single_passes_profiling_contract_to_profile_bam(tmp_path, monkeypatch):
    reference_fasta = tmp_path / "reference.fna"
    bed_file = tmp_path / "genomes.bed"
    bam_file = tmp_path / "sample.bam"
    stb_file = tmp_path / "mapping.stb"
    null_model = tmp_path / "null_model.parquet"
    gene_range = tmp_path / "gene_range.tsv"
    profiling_contract = tmp_path / "profiling_contract.json"
    reference_fasta.write_text(">chr1\nACGTACGTAA\n")
    bed_file.write_text("chr1\t0\t10\n")
    bam_file.write_text("")
    stb_file.write_text("chr1\tgenome1\n")
    gene_range.write_text("gene1\tchr1\t1\t10\n")
    pl.DataFrame({"cov": [0, 1], "max_error_count": [0, 0]}).write_parquet(null_model)
    cli.ut.write_profile_contract_file(
        {
            "reference_hash": "refhash",
            "gene_hash": "genehash",
            "null_model_hash": "modelhash",
            "stb_hash": "stbhash",
        },
        profiling_contract,
    )

    captured = {}

    def _fake_profile_bam(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.pf, "profile_bam", _fake_profile_bam)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "profile-single",
            "--reference-fasta",
            str(reference_fasta),
            "--bed-file",
            str(bed_file),
            "--bam-file",
            str(bam_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model),
            "--gene-range-table",
            str(gene_range),
            "--profiling-contract",
            str(profiling_contract),
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["profile_contract"] == {
        "reference_hash": "refhash",
        "gene_hash": "genehash",
        "null_model_hash": "modelhash",
        "stb_hash": "stbhash",
    }
    assert captured["min_mapq"] == cli.pf.PROFILE_MIN_MAPQ_DEFAULT
    assert captured["min_baseq"] == cli.pf.PROFILE_MIN_BASEQ_DEFAULT
    assert captured["min_read_ani"] == cli.pf.PROFILE_MIN_READ_ANI_DEFAULT
    assert captured["read_inclusion"] == cli.pf.PROFILE_READ_INCLUSION_DEFAULT


def test_profile_single_passes_custom_read_filters_to_profile_bam(tmp_path, monkeypatch):
    bed_file = tmp_path / "genomes.bed"
    bam_file = tmp_path / "sample.bam"
    stb_file = tmp_path / "mapping.stb"
    null_model = tmp_path / "null_model.parquet"
    bed_file.write_text("chr1\t0\t10\n")
    bam_file.write_text("")
    stb_file.write_text("chr1\tgenome1\n")
    pl.DataFrame({"cov": [0, 1], "max_error_count": [0, 0]}).write_parquet(null_model)

    captured = {}

    def _fake_profile_bam(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.pf, "profile_bam", _fake_profile_bam)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "profile-single",
            "--bed-file",
            str(bed_file),
            "--bam-file",
            str(bam_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model),
            "--min-mapq",
            "9",
            "--min-baseq",
            "17",
            "--min-read-ani",
            "0.96",
            "--read-inclusion",
            "paired",
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["min_mapq"] == 9
    assert captured["min_baseq"] == 17
    assert captured["min_read_ani"] == 0.96
    assert captured["read_inclusion"] == "paired"


def test_merge_stat_tables_command(tmp_path):
    stat_a = tmp_path / "sample_a_genome_stats.parquet"
    stat_b = tmp_path / "sample_b_genome_stats.parquet"
    output_file = tmp_path / "merged_stats.parquet"
    pl.DataFrame(
        {
            "genome": ["g1"],
            "length": [100],
            "breadth": [0.5],
            "coverage": [10.0],
        }
    ).write_parquet(stat_a)
    pl.DataFrame(
        {
            "genome": ["g2"],
            "length": [200],
            "breadth": [0.7],
            "coverage": [20.0],
        }
    ).write_parquet(stat_b)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "merge-stat-tables",
            "--stat-table",
            str(stat_a),
            "--stat-table",
            str(stat_b),
            "--output-file",
            str(output_file),
        ],
    )
    assert result.exit_code == 0, result.output
    merged = pl.read_parquet(output_file).sort("sample")
    assert merged["sample"].to_list() == ["sample_a", "sample_b"]


def test_get_coverage_stats_command(profile_2: pl.LazyFrame, tmp_path):
    profile_path = tmp_path / "profile.parquet"
    gene_bed_path = tmp_path / "genes.tsv"
    genome_bed_path = tmp_path / "genomes.bed"

    profile_2.sink_parquet(profile_path)
    pl.DataFrame(
        {
            "gene": ["gene1", "gene2", "gene3", "gene1", "gene2", "gene3"],
            "scaffold": ["chr1", "chr2", "chr2", "chr3", "chr3", "chr3"],
            "start": [2, 2, 11, 3, 11, 22],
            "end": [5, 6, 15, 6, 15, 26],
        }
    ).write_csv(gene_bed_path, separator="\t", include_header=False)
    pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2", "chr3"],
            "start": [0, 0, 0],
            "end": [10, 20, 30],
        }
    ).write_csv(genome_bed_path, separator="\t", include_header=False)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "get-coverage-stats",
            "--profile-parquet",
            str(profile_path),
            "--gene-bed",
            str(gene_bed_path),
            "--genome-bed",
            str(genome_bed_path),
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "sample1",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "sample1_gene_stats.parquet" in result.output
    assert "sample1_genome_stats.parquet" in result.output

    gene_stats = pl.read_parquet(tmp_path / "sample1_gene_stats.parquet")
    genome_stats = pl.read_parquet(tmp_path / "sample1_genome_stats.parquet")
    assert gene_stats.columns == ["genome", "gene", "length", "breadth", "coverage", "5x_cov_sites", "ber"]
    assert genome_stats.columns == ["genome", "length", "breadth", "coverage", "5x_cov_sites", "ber"]
    assert gene_stats.height == 6
    assert genome_stats.height == 2


def test_run_group_removed():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["run", "--help"])
    assert result.exit_code != 0


def test_profile_utility_commands_exist_under_utilities():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["utilities", "--help"])
    assert result.exit_code == 0
    assert "prepare_profiling" in result.output
    assert "profile-single" in result.output
    assert "adjust-sequence-errors" in result.output


def test_profile_subcommands_removed():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["profile", "batch"])
    assert result.exit_code != 0


def test_gene_tools_commands_exist_under_utilities():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["utilities", "--help"])
    assert result.exit_code == 0
    assert "gene-range-table" in result.output


def test_gene_tools_group_removed():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["gene_tools", "--help"])
    assert result.exit_code != 0


def test_generate_stb(tmp_path):
    # Create multiple fasta files in a temporary directory
    fasta1 = tmp_path / "genome1.fasta"
    fasta1.write_text(">chr1_1\nACGTACGTAC\n>chr2_1\nACGTACGTACGTACGTACGT\n")
    fasta2 = tmp_path / "genome2.fasta"
    fasta2.write_text(">chr1_2\nACGTACGTACGTACGTACGTACGTACGT\n")
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["utilities","generate_stb", "-g",  str(tmp_path) , "-o", str(tmp_path / "stb_output.tsv"), "--extension", ".fasta"])
    assert result.exit_code == 0
    stb_path = tmp_path / "stb_output.tsv"
    assert stb_path.exists()
    stb=pl.read_csv(stb_path,separator="\t",has_header=False).rename(
        {
            "column_1":"scaffold",
            "column_2":"genome"
        }
    )
    result_dict = stb.rows_by_key("scaffold", unique=True, named=True)
    assert result_dict["chr1_1"]["genome"] == "genome1"
    assert result_dict["chr2_1"]["genome"] == "genome1"
    assert result_dict["chr1_2"]["genome"] == "genome2"


def test_cli_build_genome_db_from_cache(tmp_path):
    abundance = tmp_path / "sylph.csv"
    abundance.write_text(
        "Genome_file,sample1,sample2\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.6,0.0\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.0,0.4\n"
    )
    cache_dir = tmp_path / "cache"
    output_dir = tmp_path / "out"
    genomes_dir = cache_dir / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    fasta_1 = genomes_dir / "GCF_000001405.40.fna"
    fasta_1.write_text(">chr1\nATCG\n")
    fasta_2 = genomes_dir / "GCA_123456.1.fna"
    fasta_2.write_text(">chr2\nGGTT\n")

    local_db = cli.bdb.LocalGenomeDB(cache_dir / cli.bdb.LocalGenomeDB.DEFAULT_DB_NAME)
    local_db.add_or_update_genomes(
        pl.DataFrame(
            {
                "accession": ["GCF_000001405.40", "GCA_123456.1"],
                "genome_name": ["GCF_000001405.40", "GCA_123456.1"],
                "download_url": [None, None],
            }
        ),
        source_tool="sylph",
    )
    local_db.set_location("GCF_000001405.40", fasta_1)
    local_db.set_location("GCA_123456.1", fasta_2)
    local_db.sync()
    local_db.save()

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-genome-db",
            "--tool",
            "sylph",
            "--abundance-table",
            str(abundance),
            "--cache-dir",
            str(cache_dir),
            "--output-dir",
            str(output_dir),
        ],
    )
    assert result.exit_code == 0
    assert "Genome DB Build Report" in result.output
    assert "Selected genomes (non-zero)" in result.output
    assert "Cache directory" in result.output
    assert (output_dir / "reference_genomes.fna").exists()
    assert (output_dir / "reference_genomes.stb").exists()
    report_file = output_dir / "genome_db_build_report.txt"
    assert report_file.exists()
    report_text = report_file.read_text()
    assert "Genome DB Build Report" in report_text
    assert "Failed IDs (0):" in report_text


def test_cli_build_genome_db_passes_retry_options(tmp_path, monkeypatch):
    abundance = tmp_path / "sylph.csv"
    abundance.write_text(
        "Genome_file,Eff_cov\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.6\n"
    )
    cache_dir = tmp_path / "cache"
    output_dir = tmp_path / "out"
    captured = {}

    def _fake_build_reference_from_abundance(**kwargs):
        captured.update(kwargs)
        out_dir = kwargs["output_dir"]
        out_dir.mkdir(parents=True, exist_ok=True)
        out_fasta = out_dir / "reference_genomes.fna"
        out_stb = out_dir / "reference_genomes.stb"
        out_fasta.write_text(">x\nATCG\n")
        out_stb.write_text("x\tg\n")
        extracted = pl.DataFrame({"accession": ["GCF_000001405.40"]})
        report = pl.DataFrame(
            {
                "accession": ["GCF_000001405.40"],
                "status": ["failed"],
                "location": [None],
                "url": ["https://example.test/GCF_000001405.40.fna"],
                "error": ["simulated error"],
            }
        )
        summary = {
            "selected_genomes": 1,
            "cached_before_download": 0,
            "attempted_downloads": 1,
            "downloaded_now": 0,
            "failed_downloads": 1,
            "cached_after_download": 0,
            "missing_after_retries": 1,
            "max_download_attempts": kwargs["max_download_attempts"],
        }
        return out_fasta, out_stb, extracted, report, summary

    monkeypatch.setattr(cli.bdb, "build_reference_from_abundance", _fake_build_reference_from_abundance)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-genome-db",
            "--tool",
            "sylph",
            "--abundance-table",
            str(abundance),
            "--cache-dir",
            str(cache_dir),
            "--output-dir",
            str(output_dir),
            "--download-retries",
            "5",
            "--retry-backoff-seconds",
            "0.25",
            "--download-workers",
            "6",
        ],
    )
    assert result.exit_code == 0
    assert captured["max_download_attempts"] == 5
    assert captured["backoff_base_seconds"] == 0.25
    assert captured["download_workers"] == 6
    report_text = (output_dir / "genome_db_build_report.txt").read_text()
    assert "Failed IDs (1):" in report_text
    assert "GCF_000001405.40" in report_text


def test_cli_build_profile_db_minimal_columns(tmp_path):
    profile_path = tmp_path / "sample_profile.parquet"
    pl.DataFrame({
        "chrom": ["chr2", "chr1", "chr1"],
        "genome": ["genome1", "genome1", "genome1"],
        "pos": [0, 1, 2],
        "gene": ["NA", "NA", "NA"],
        "A": [1, 0, 0],
        "T": [0, 1, 0],
        "C": [0, 0, 1],
        "G": [0, 0, 0],
    }).write_parquet(profile_path)

    input_csv = tmp_path / "profiles.csv"
    pl.DataFrame({
        "profile_name": ["sample_1"],
        "profile_location": [str(profile_path)],
    }).write_csv(input_csv)

    output_db = tmp_path / "profiles.parquet"
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-profile-db",
            "-p",
            str(input_csv),
            "-o",
            str(output_db),
        ],
    )
    assert result.exit_code == 0
    assert output_db.exists()

    built_db = pl.read_parquet(output_db)
    assert set(built_db.columns) == {"profile_name", "profile_location"}


def test_cli_build_profile_db_rejects_metadata_mismatch_without_allow_mismatch(tmp_path):
    profile_a = tmp_path / "a.parquet"
    profile_b = tmp_path / "b.parquet"
    base = pl.DataFrame(
        {"chrom": ["chr1"], "genome": ["g1"], "pos": [1], "gene": ["NA"], "A": [1], "T": [0], "C": [0], "G": [0]}
    )
    base.write_parquet(profile_a, metadata={cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_a"})
    base.write_parquet(profile_b, metadata={cli.ut.PROFILE_REFERENCE_HASH_METADATA_KEY: "ref_b"})

    input_csv = tmp_path / "profiles_mismatch.csv"
    pl.DataFrame(
        {
            "profile_name": ["a", "b"],
            "profile_location": [str(profile_a), str(profile_b)],
        }
    ).write_csv(input_csv)

    output_db = tmp_path / "profiles_mismatch.parquet"
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-profile-db",
            "-p",
            str(input_csv),
            "-o",
            str(output_db),
        ],
    )
    assert result.exit_code != 0
    assert "Profile contract metadata do not match" in str(result.exception)

    allow_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-profile-db",
            "-p",
            str(input_csv),
            "-o",
            str(output_db),
            "--allow-mismatch",
        ],
    )
    assert allow_result.exit_code == 0


def test_cli_to_complete_table_uses_current_comparison_state(tmp_path):
    profile_db, profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    current_compare = tmp_path / "current_genome_compare.parquet"
    output_csv = tmp_path / "remaining_pairs.csv"
    _write_genome_compare_table_for_config_tests(current_compare)

    runner = CliRunner()
    complete_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "to-complete-table",
            "--profile-db",
            str(profile_db),
            "--comp-db-file",
            str(current_compare),
            "--output-file",
            str(output_csv),
        ],
    )
    assert complete_result.exit_code == 0

    remaining = pl.read_csv(output_csv)
    assert remaining.select(["sample_name_1", "sample_name_2"]).rows() == [
        ("sample_a", "sample_c"),
        ("sample_b", "sample_c"),
    ]
    assert remaining["profile_location_1"].to_list() == [
        str(profile_paths["sample_a"]),
        str(profile_paths["sample_b"]),
    ]
    assert remaining["profile_location_2"].to_list() == [
        str(profile_paths["sample_c"]),
        str(profile_paths["sample_c"]),
    ]


def test_cli_merge_parquet_batched_mode(tmp_path):
    input_dir = tmp_path / "parts"
    input_dir.mkdir()
    for idx in range(5):
        pl.DataFrame({"part": [idx], "value": [idx + 1]}).write_parquet(input_dir / f"part_{idx:03d}.parquet")

    output_file = tmp_path / "merged.parquet"
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "merge_parquet",
            "--input-dir",
            str(input_dir),
            "--output-file",
            str(output_file),
            "--batch-size",
            "2",
        ],
    )

    assert result.exit_code == 0
    assert output_file.exists()
    merged = pl.read_parquet(output_file).sort("part")
    assert merged.shape == (5, 2)
    assert "batch 1/3 start" in result.output


def test_cli_merge_parquet_rejects_compare_metadata_mismatch(tmp_path):
    input_dir = tmp_path / "parts"
    input_dir.mkdir()
    common_metadata = {
        cli.ut.COMPARE_KIND_METADATA_KEY: "genome",
        cli.ut.COMPARE_MIN_COV_METADATA_KEY: "5",
    }
    pl.DataFrame({"part": [0]}).write_parquet(
        input_dir / "part_a.parquet",
        metadata={
            **common_metadata,
            cli.ut.COMPARE_SCOPE_METADATA_KEY: "all",
        },
    )
    pl.DataFrame({"part": [1]}).write_parquet(
        input_dir / "part_b.parquet",
        metadata={
            **common_metadata,
            cli.ut.COMPARE_SCOPE_METADATA_KEY: "genome1",
        },
    )

    output_file = tmp_path / "merged.parquet"
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "merge_parquet",
            "--input-dir",
            str(input_dir),
            "--output-file",
            str(output_file),
            "--batch-size",
            "1",
        ],
    )

    assert result.exit_code != 0
    assert "Parquet metadata mismatch across inputs" in str(result.exception)


def test_cli_merge_parquet_allow_mismatch_rewrites_compare_metadata_to_na(tmp_path):
    input_dir = tmp_path / "parts"
    input_dir.mkdir()
    common_metadata = {
        cli.ut.COMPARE_KIND_METADATA_KEY: "genome",
        cli.ut.COMPARE_MIN_COV_METADATA_KEY: "5",
    }
    pl.DataFrame({"part": [0]}).write_parquet(
        input_dir / "part_a.parquet",
        metadata={
            **common_metadata,
            cli.ut.COMPARE_SCOPE_METADATA_KEY: "all",
        },
    )
    pl.DataFrame({"part": [1]}).write_parquet(
        input_dir / "part_b.parquet",
        metadata={
            **common_metadata,
            cli.ut.COMPARE_SCOPE_METADATA_KEY: "genome1",
        },
    )

    output_file = tmp_path / "merged.parquet"
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "merge_parquet",
            "--input-dir",
            str(input_dir),
            "--output-file",
            str(output_file),
            "--batch-size",
            "1",
            "--allow-mismatch",
        ],
    )

    assert result.exit_code == 0
    metadata = pl.read_parquet_metadata(output_file)
    assert metadata[cli.ut.COMPARE_KIND_METADATA_KEY] == "genome"
    assert metadata[cli.ut.COMPARE_MIN_COV_METADATA_KEY] == "5"
    assert metadata[cli.ut.COMPARE_SCOPE_METADATA_KEY] == cli.ut.COMPARE_METADATA_MISSING_VALUE


def test_cli_adjust_sequence_errors(tmp_path):
    profile_path = tmp_path / "profile.parquet"
    null_model_path = tmp_path / "null_model.parquet"
    output_path = tmp_path / "adjusted.parquet"

    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "genome": ["genome1", "genome1"],
            "gene": ["gene1", "gene1"],
            "pos": [1, 2],
            "A": [3, 2],
            "C": [1, 0],
            "G": [0, 0],
            "T": [0, 1],
        }
    ).write_parquet(profile_path)
    pl.DataFrame(
        {
            "cov": [4, 3],
            "max_error_count": [2, 2],
        }
    ).write_parquet(null_model_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "adjust-sequence-errors",
            "--profile-parquet",
            str(profile_path),
            "--null-model",
            str(null_model_path),
            "--output-file",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    adjusted = pl.read_parquet(output_path)
    assert adjusted.columns == ["chrom", "genome", "gene", "pos", "A", "C", "G", "T"]
    assert set(adjusted.select("A", "C", "G", "T").rows()) == {(3, 0, 0, 0), (2, 0, 0, 0)}


def test_cli_adjust_sequence_errors_rejects_same_input_output(tmp_path):
    profile_path = tmp_path / "profile.parquet"
    null_model_path = tmp_path / "null_model.parquet"

    pl.DataFrame(
        {
            "chrom": ["chr1"],
            "genome": ["genome1"],
            "gene": ["gene1"],
            "pos": [1],
            "A": [3],
            "C": [1],
            "G": [0],
            "T": [0],
        }
    ).write_parquet(profile_path)
    pl.DataFrame(
        {
            "cov": [4],
            "max_error_count": [2],
        }
    ).write_parquet(null_model_path)

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "adjust-sequence-errors",
            "--profile-parquet",
            str(profile_path),
            "--null-model",
            str(null_model_path),
            "--output-file",
            str(profile_path),
        ],
    )

    assert result.exit_code != 0
    assert "must be different paths" in result.output


def test_dash_h_shows_help_globally():
    runner = CliRunner()
    for args in (["-h"], ["profile", "-h"], ["utilities", "prepare_profiling", "-h"]):
        result = runner.invoke(cli.cli, args)
        assert result.exit_code == 0, result.output
        assert "Usage:" in result.output


def test_profile_help_is_organized_into_sections_with_defaults():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["profile", "--help"])
    assert result.exit_code == 0, result.output
    for header in (
        "Required inputs:",
        "Optional inputs:",
        "Optional pre-built assets",
        "Profiling parameters:",
        "Running parameters:",
    ):
        assert header in result.output
    # Defaults are surfaced for parameters that have them.
    assert "[default: 0.001]" in result.output   # --error-rate
    assert "[default: 8]" in result.output        # --num-procs
    assert "[default: local]" in result.output    # --execution-mode


def _write_small_parquet(path, rows=3):
    pl.DataFrame({"genome": ["g"] * rows, "genome_pop_ani": [100.0] * rows}).write_parquet(path)


def test_maybe_write_csv_writes_small_and_respects_flags(tmp_path):
    pq = tmp_path / "all_comparisons.parquet"
    _write_small_parquet(pq)

    # Default: small parquet gets a companion csv.
    out = cli._maybe_write_csv(pq, no_csv=False, force_csv=False)
    assert out == tmp_path / "all_comparisons.csv"
    assert out.exists()
    out.unlink()

    # --no-csv: nothing written.
    assert cli._maybe_write_csv(pq, no_csv=True, force_csv=False) is None
    assert not (tmp_path / "all_comparisons.csv").exists()


def test_maybe_write_csv_skips_huge_unless_forced(tmp_path, monkeypatch):
    pq = tmp_path / "big.parquet"
    _write_small_parquet(pq)
    # Pretend the estimate is over the threshold.
    monkeypatch.setattr(cli, "_estimated_csv_mb", lambda p: cli.CSV_SIZE_THRESHOLD_MB + 1)

    assert cli._maybe_write_csv(pq, no_csv=False, force_csv=False) is None
    assert not (tmp_path / "big.csv").exists()

    forced = cli._maybe_write_csv(pq, no_csv=False, force_csv=True)
    assert forced is not None and forced.exists()


def test_finalize_profile_outputs_writes_stat_csvs(tmp_path):
    sample = tmp_path / "sampleA"
    sample.mkdir()
    # genome_stats needs the presence-input columns.
    pl.DataFrame({"genome": ["g"], "coverage": [10.0], "breadth": [1.0], "ber": [1.0], "fug": [0.1]}).write_parquet(
        sample / "sampleA_genome_stats.parquet"
    )
    _write_small_parquet(sample / "sampleA_gene_stats.parquet")
    _write_small_parquet(sample / "sampleA_profile.parquet")  # no ref bitmask -> no SNVs, no csv

    cli._finalize_profile_outputs(
        tmp_path,
        no_csv=False,
        force_csv=False,
        emit_snvs=True,
        snv_min_cov=5,
        presence_ber=0.5,
        presence_fug=2.0,
        presence_min_cov_use_fug=2.0,
        presence_min_coverage=0.1,
    )

    assert (sample / "sampleA_genome_stats.csv").exists()
    assert (sample / "sampleA_gene_stats.csv").exists()
    assert not (sample / "sampleA_profile.csv").exists()
    # Presence column was added to the genome stats table.
    assert "presence" in pl.read_parquet(sample / "sampleA_genome_stats.parquet").columns
    # No reference bitmask in the profile -> SNV calling skipped without error.
    assert not (sample / "sampleA_SNVs.parquet").exists()


def test_add_presence_column_calls_present_and_absent(tmp_path):
    stats = pl.DataFrame(
        {
            "genome": ["present_one", "absent_one"],
            "coverage": [10.0, 0.0],
            "breadth": [1.0, 0.0],
            "ber": [1.0, 0.0],
            "fug": [0.1, 5.0],
        }
    )
    out = cli._add_presence_column(stats, ber=0.5, fug=1.0, min_cov_use_fug=2.0, min_coverage=0.1)
    by_genome = {row["genome"]: row["presence"] for row in out.iter_rows(named=True)}
    assert by_genome == {"present_one": "present", "absent_one": "absent"}


def test_presence_fug_direction_and_coverage():
    # coverage 0.5 is between min_coverage (0.1) and min_cov_use_fug (2), so these
    # rows go through the FUG path; low_cov is below the coverage floor.
    stats = pl.DataFrame(
        {
            "genome": ["hi_fug", "lo_fug", "low_cov"],
            "coverage": [0.5, 0.5, 0.05],
            "ber":      [0.9, 0.9, 0.9],
            "fug":      [0.8, 0.3, 0.8],   # 0.8/0.632=1.27 > 1; 0.3/0.632=0.47 < 1
        }
    )
    out = cli._add_presence_column(stats, ber=0.5, fug=1.0, min_cov_use_fug=2.0, min_coverage=0.1)
    calls = {r["genome"]: r["presence"] for r in out.iter_rows(named=True)}
    assert calls["hi_fug"] == "present"     # fug/0.632 = 1.27 > 1 -> present (correct direction)
    assert calls["lo_fug"] == "absent"      # fug/0.632 = 0.47 < 1 -> absent
    assert calls["low_cov"] == "absent"     # coverage 0.05 < 0.1 floor


def test_discover_taxonomy_file_next_to_reference_and_stb(tmp_path):
    ref = tmp_path / "reference_genomes.fna"; ref.write_text(">x\nA\n")
    stb = tmp_path / "reference_genomes.stb"; stb.write_text("x\tg\n")
    tax = tmp_path / "reference_genomes_taxonomy.tsv"; tax.write_text("genome\tgenome_taxonomy\ng\td__Bacteria\n")
    assert cli._discover_taxonomy_file(str(ref), str(stb), None) == tax
    # explicit wins
    other = tmp_path / "custom.tsv"; other.write_text("genome\tgenome_taxonomy\n")
    assert cli._discover_taxonomy_file(str(ref), str(stb), str(other)) == other
    # none when absent
    tax.unlink()
    assert cli._discover_taxonomy_file(str(ref), str(stb), None) is None


def test_finalize_adds_genome_taxonomy_column(tmp_path):
    sample = tmp_path / "sampleA"; sample.mkdir()
    pl.DataFrame({"genome": ["GCF_1.1", "GCF_2.1"], "coverage": [10.0, 10.0], "ber": [1.0, 1.0], "fug": [1.0, 1.0]}).write_parquet(
        sample / "sampleA_genome_stats.parquet"
    )
    tax = tmp_path / "tax.tsv"; tax.write_text("genome\tgenome_taxonomy\nGCF_1.1\td__Bacteria;s__Foo\n")

    cli._finalize_profile_outputs(
        tmp_path, no_csv=True, force_csv=False, emit_snvs=False, snv_min_cov=5,
        presence_ber=0.5, presence_fug=1.0, presence_min_cov_use_fug=2.0, presence_min_coverage=0.1,
        taxonomy_file=tax,
    )
    d = pl.read_parquet(sample / "sampleA_genome_stats.parquet")
    assert "genome_taxonomy" in d.columns
    by = {r["genome"]: r["genome_taxonomy"] for r in d.iter_rows(named=True)}
    assert by["GCF_1.1"] == "d__Bacteria;s__Foo"
    assert by["GCF_2.1"] is None  # unmatched -> null


def test_add_percent_compared_and_coverage_overlap(tmp_path):
    import pathlib
    # A run dir with two samples' genome_stats + a genome_lengths asset.
    (tmp_path / "profiling_assets").mkdir()
    pl.DataFrame({"genome": ["g1"], "genome_length": [1000]}).write_parquet(
        tmp_path / "profiling_assets" / "genome_lengths.parquet"
    )
    profiles = []
    for name, covered in [("s1", 800), ("s2", 600)]:
        d = tmp_path / name
        d.mkdir()
        prof = d / f"{name}_profile.parquet"
        pl.DataFrame({"genome": ["g1"]}).write_parquet(prof)  # placeholder
        pl.DataFrame({"genome": ["g1"], "5x_cov_sites": [covered]}).write_parquet(
            d / f"{name}_genome_stats.parquet"
        )
        profiles.append((name, str(prof)))

    comp_path = tmp_path / "all_comparisons.parquet"
    pl.DataFrame(
        {"sample_1": ["s1"], "sample_2": ["s2"], "genome": ["g1"], "total_positions": [500]}
    ).write_parquet(comp_path)

    cli._add_percent_compared(comp_path, cli._discover_genome_lengths(tmp_path, [p for _n, p in profiles]))
    cli._add_coverage_overlap(comp_path, profiles)

    out = pl.read_parquet(comp_path)
    # percent_compared = 500 / 1000
    assert out["percent_compared"][0] == pytest.approx(0.5)
    # coverage_overlap = 500 / (800 + 600 - 500) = 500/900
    assert out["coverage_overlap"][0] == pytest.approx(500 / 900)


def test_get_snp_reference_cli_vcf_includes_subconsensus_snv_sites(tmp_path):
    """`--fmt vcf` emits any site with a non-reference allele, including SNV sites
    where the reference is still the consensus (the superset definition)."""
    from zipstrain import utils

    # pos 1: ref=A, counts A=1 C=7 -> consensus C (differs from ref)  [con != ref]
    # pos 2: ref=A, counts A=8 C=2 -> consensus A (== ref), minor C   [SNV: con == ref]
    profile_path = tmp_path / "profile_subconsensus.parquet"
    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "genome": ["g", "g"],
            "gene": ["NA", "NA"],
            "pos": [1, 2],
            "A": [1, 8],
            "C": [7, 2],
            "G": [0, 0],
            "T": [0, 0],
            utils.REF_BASE_BITMASK_COLUMN: [1, 1],
        }
    ).write_parquet(profile_path)
    output_path = tmp_path / "reference_snps.vcf"

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities", "get-snp-reference",
            "--profile-file", str(profile_path),
            "--min-cov", "5",
            "--fmt", "vcf",
            "--output-file", str(output_path),
        ],
    )
    assert result.exit_code == 0, result.output
    positions = [
        line.split("\t")[1]
        for line in output_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    # Superset: both the consensus-differs site (1) AND the sub-consensus SNV site (2).
    assert positions == ["1", "2"], positions
