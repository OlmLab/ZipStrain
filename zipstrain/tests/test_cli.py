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
        "compare",
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
        "compare",
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
        "compare",
        "single_compare_genome", 
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_3_dir),
        "--stb-file", str(stb_path),
        "--output-file", str(tmp_path/"output.parquet"),
    ])
    lf2 = pl.read_parquet(tmp_path/"output.parquet")
    assert result.exit_code == 0
    assert lf2.shape[0] == 2
    


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


def test_cli_build_genome_db_no_download(tmp_path):
    abundance = tmp_path / "sylph.csv"
    abundance.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.6\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.4\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    report_file = tmp_path / "report.csv"

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
            "--db-file",
            str(db_path),
            "--genomes-dir",
            str(genomes_dir),
            "--no-download",
            "--report-file",
            str(report_file),
        ],
    )
    assert result.exit_code == 0
    assert db_path.exists()
    db_df = pl.read_parquet(db_path)
    assert set(db_df["accession"].to_list()) == {"GCF_000001405.40", "GCA_123456.1"}
    assert report_file.exists()
