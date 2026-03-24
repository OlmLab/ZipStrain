from click.testing import CliRunner
from zipstrain import cli
import pytest
import polars as pl

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
    assert {"compare_gene", "compare_genome", "profile", "test", "utilities"} <= commands
    assert "run" not in commands
    assert "compare" not in commands
    assert "gene_tools" not in commands

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
        "--mpileup-contig-2", str(profile_2_dir),
        "--stb-file", str(stb_path),
        "--calculate", "ani",
        "--output-file", str(tmp_path/"output_ani.parquet"),
    ])
    lf1_ani = pl.read_parquet(tmp_path/"output_ani.parquet")
    assert result.exit_code == 0
    assert set(lf1_ani.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "sample_1",
        "sample_2",
    }
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


def test_run_compare_genomes_passes_duckdb_threads(tmp_path, monkeypatch):
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
            "compare_genome",
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
    assert captured["ani_method"] == "conani"
    assert captured["calculate"] == "ani"


def test_run_compare_genes_passes_duckdb_threads(tmp_path, monkeypatch):
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
            "compare_gene",
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


def test_profile_command_passes_null_model(tmp_path, monkeypatch):
    input_table = tmp_path / "bams.csv"
    input_table.write_text("sample_name,bamfile\ns1,/tmp/a.bam\n")
    stb_file = tmp_path / "stb.tsv"
    null_model = tmp_path / "null_model.parquet"
    gene_range = tmp_path / "genes.tsv"
    bed_file = tmp_path / "regions.bed"
    genome_lengths = tmp_path / "genome_lengths.parquet"
    for path in [stb_file, null_model, gene_range, bed_file, genome_lengths]:
        path.write_text("x\n")

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
            str(stb_file),
            "--null-model",
            str(null_model),
            "--gene-range-table",
            str(gene_range),
            "--bed-file",
            str(bed_file),
            "--genome-length-file",
            str(genome_lengths),
            "--run-dir",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == 0
    assert captured["null_model_file"] == null_model


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
