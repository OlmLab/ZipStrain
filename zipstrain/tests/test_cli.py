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
            "reference_db_id": ["ref_1", "ref_1", "ref_1"],
            "gene_db_id": ["gene_ref_1", "gene_ref_1", "gene_ref_1"],
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
    assert {"compare", "profile", "test", "utilities"} == commands


def test_cli_version_flag():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--version"])
    assert result.exit_code == 0
    assert "zipstrain" in result.output.lower()
    assert _project_version() in result.output

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
    assert lf1.sort("genome").equals(lf1_duckdb.sort("genome"))
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


def test_generate_genome_pairs_command(profile_1: pl.LazyFrame, profile_2: pl.LazyFrame, profile_3: pl.LazyFrame, tmp_path):
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir()
    profile_1.sink_parquet(profile_dir / "profile_1.parquet")
    profile_2.sink_parquet(profile_dir / "profile_2.parquet")
    profile_3.sink_parquet(profile_dir / "profile_3.parquet")
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


def test_compare_genomes_batch_passes_duckdb_threads(tmp_path, monkeypatch):
    comp_obj = tmp_path / "genome_comp.json"
    comp_obj.write_text("{}")
    captured = {}

    monkeypatch.setattr(
        cli.db.GenomeComparisonDatabase,
        "load_obj",
        classmethod(lambda cls, path: object()),
    )

    def _fake_lazy_run_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_compares", _fake_lazy_run_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "genomes",
            "--genome-comparison-object",
            str(comp_obj),
            "--run-dir",
            str(tmp_path / "run"),
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


def test_compare_genes_batch_passes_duckdb_threads(tmp_path, monkeypatch):
    comp_obj = tmp_path / "gene_comp.json"
    comp_obj.write_text("{}")
    captured = {}

    monkeypatch.setattr(
        cli.db.GeneComparisonDatabase,
        "load_obj",
        classmethod(lambda cls, path: object()),
    )

    def _fake_lazy_run_gene_compares(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli.tm, "lazy_run_gene_compares", _fake_lazy_run_gene_compares)
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "compare",
            "genes",
            "--gene-comparison-object",
            str(comp_obj),
            "--run-dir",
            str(tmp_path / "run"),
            "--engine",
            "duckdb",
            "--duckdb-threads",
            "9",
        ],
    )
    assert result.exit_code == 0
    assert captured["duckdb_threads"] == 9
    assert captured["compare_engine"] == "duckdb"


def test_profile_command_calls_lazy_run_profile(tmp_path, monkeypatch):
    input_table = tmp_path / "input.csv"
    input_table.write_text("sample_name,bamfile\nsample1,/tmp/sample1.bam\n")
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


def test_profile_command_allows_missing_gene_range_table(tmp_path, monkeypatch):
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
    assert captured["gene_range_table"] is None


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
        "reference_db_id": ["ref_1"],
        "gene_db_id": ["gene_ref_1"],
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
    assert set(built_db.columns) == {"profile_name", "profile_location", "reference_db_id", "gene_db_id"}


def test_cli_build_genome_comparison_config_with_current_comp_table(tmp_path):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    current_compare = tmp_path / "current_genome_compare.parquet"
    output_json = tmp_path / "genome_compare_config.json"
    stb_file = tmp_path / "reference.stb"
    _write_genome_compare_table_for_config_tests(current_compare)
    stb_file.write_text("chr1\tgenome1\n")

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-genome-comparison-config",
            "--profile-db",
            str(profile_db),
            "--gene-db-id",
            "gene_ref_1",
            "--reference-genome-id",
            "ref_1",
            "--scope",
            "genome1",
            "--min-cov",
            "7",
            "--min-gene-compare-len",
            "150",
            "--stb-file-loc",
            str(stb_file),
            "--current-comp-table",
            str(current_compare),
            "--output-file",
            str(output_json),
        ],
    )

    assert result.exit_code == 0
    loaded = db.GenomeComparisonDatabase.load_obj(output_json)
    assert loaded.comp_db_loc == current_compare.resolve()
    assert loaded.config.reference_id == "ref_1"
    assert loaded.config.gene_db_id == "gene_ref_1"
    assert loaded.config.scope == "genome1"
    assert loaded.config.min_cov == 7
    assert loaded.config.min_gene_compare_len == 150
    assert loaded.config.stb_file_loc == str(stb_file)


def test_cli_build_gene_comparison_config_with_current_comp_table(tmp_path):
    profile_db, _profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    current_compare = tmp_path / "current_gene_compare.parquet"
    output_json = tmp_path / "gene_compare_config.json"
    stb_file = tmp_path / "reference.stb"
    _write_gene_compare_table_for_config_tests(current_compare)
    stb_file.write_text("chr1\tgenome1\n")

    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-gene-comparison-config",
            "--profile-db",
            str(profile_db),
            "--gene-db-id",
            "gene_ref_1",
            "--reference-genome-id",
            "ref_1",
            "--scope",
            "genome1:gene1",
            "--min-cov",
            "6",
            "--min-gene-compare-len",
            "125",
            "--stb-file-loc",
            str(stb_file),
            "--current-comp-table",
            str(current_compare),
            "--output-file",
            str(output_json),
        ],
    )

    assert result.exit_code == 0
    loaded = db.GeneComparisonDatabase.load_obj(output_json)
    assert loaded.comp_db_loc == current_compare.resolve()
    assert loaded.config.reference_genome_id == "ref_1"
    assert loaded.config.gene_db_id == "gene_ref_1"
    assert loaded.config.scope == "genome1:gene1"
    assert loaded.config.min_cov == 6
    assert loaded.config.min_gene_compare_len == 125
    assert loaded.config.stb_file_loc == str(stb_file)


def test_cli_to_complete_table_uses_current_comparison_state(tmp_path):
    profile_db, profile_paths = _write_profile_db_for_compare_config_tests(tmp_path)
    current_compare = tmp_path / "current_genome_compare.parquet"
    config_json = tmp_path / "genome_compare_config.json"
    output_csv = tmp_path / "remaining_pairs.csv"
    stb_file = tmp_path / "reference.stb"
    _write_genome_compare_table_for_config_tests(current_compare)
    stb_file.write_text("chr1\tgenome1\n")

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-genome-comparison-config",
            "--profile-db",
            str(profile_db),
            "--gene-db-id",
            "gene_ref_1",
            "--reference-genome-id",
            "ref_1",
            "--scope",
            "all",
            "--min-cov",
            "5",
            "--min-gene-compare-len",
            "100",
            "--stb-file-loc",
            str(stb_file),
            "--current-comp-table",
            str(current_compare),
            "--output-file",
            str(config_json),
        ],
    )
    assert build_result.exit_code == 0

    complete_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "to-complete-table",
            "--genome-comparison-object",
            str(config_json),
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
