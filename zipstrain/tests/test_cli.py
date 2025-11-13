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
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2))),
        "gene": gene_locs[:30],
        "A": a_chr1 + a_chr2 ,
        "T": t_chr1 + t_chr2 ,
        "C": c_chr1 + c_chr2 ,
        "G": g_chr1 + g_chr2 ,
    }).lazy()

@pytest.fixture
def scaffold_1()->str:
    return "chr1\nchr2"

@pytest.fixture
def profile_2()->pl.LazyFrame:
    """Exactly the same as profile_1"""
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
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
        "pos":list(range(0,len(a_chr3))),
        "gene": gene_locs[30:],
        "A": a_chr3,
        "T": t_chr3,
        "C": c_chr3,
        "G": g_chr3,
    }).lazy()
@pytest.fixture
def scaffold_2():
    return "chr1\nchr2\nchr3"

@pytest.fixture
def scaffold_3():
    return "chr3"

@pytest.fixture
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()

@pytest.fixture
def null_model()->pl.LazyFrame:
    return pl.DataFrame({
        "cov":list(range(100)),
        "max_error_count":[int(i*0.1) for i in range(100)],
    }).lazy()
    

def test_cli_profile_compare(profile_1:pl.LazyFrame,
                             profile_2:pl.LazyFrame,
                             profile_3:pl.LazyFrame,
                             scaffold_1:str,
                             scaffold_2:str,
                             scaffold_3:str,
                             stb:pl.LazyFrame,
                             null_model:pl.LazyFrame,
                             tmp_path):
    profile_1_dir=tmp_path/"profile_1.parquet"
    profile_2_dir=tmp_path/"profile_2.parquet"
    profile_3_dir=tmp_path/"profile_3.parquet"
    scaffold_1_path=tmp_path/"scaffold_1.txt"
    scaffold_2_path=tmp_path/"scaffold_2.txt"
    scaffold_3_path=tmp_path/"scaffold_3.txt"
    stb_path=tmp_path/"stb.parquet"
    null_model_path=tmp_path/"null_model.parquet"
    profile_1.sink_parquet(profile_1_dir)
    profile_2.sink_parquet(profile_2_dir)
    profile_3.sink_parquet(profile_3_dir)
    scaffold_1_path.write_text(scaffold_1)
    scaffold_2_path.write_text(scaffold_2)
    scaffold_3_path.write_text(scaffold_3)
    stb.sink_csv(stb_path,separator="\t",include_header=False)
    null_model.sink_parquet(null_model_path)
    runner = CliRunner()
    result = runner.invoke(cli.cli, [
        "compare",
        "single_compare_genome", 
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_2_dir),
        "--scaffolds-1", str(scaffold_1_path),
        "--scaffolds-2", str(scaffold_2_path),
        "--stb-file", str(stb_path),
        "--null-model", str(null_model_path),
        "--output-file", str(tmp_path/"output.parquet"),
    ])
    lf1 = pl.read_parquet(tmp_path/"output.parquet")
    assert result.exit_code == 0 
    assert lf1.shape[0] == 2
    result = runner.invoke(cli.cli, [
        "compare",
        "single_compare_genome", 
        "--mpileup-contig-1", str(profile_1_dir),
        "--mpileup-contig-2", str(profile_3_dir),
        "--scaffolds-1", str(scaffold_1_path),
        "--scaffolds-2", str(scaffold_3_path),
        "--stb-file", str(stb_path),
        "--null-model", str(null_model_path),
        "--output-file", str(tmp_path/"output.parquet"),
    ])
    lf2 = pl.read_parquet(tmp_path/"output.parquet")
    assert result.exit_code == 0
    assert lf2.shape[0] == 2