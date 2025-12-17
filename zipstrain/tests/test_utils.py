from zipstrain import database,compare,utils
import polars as pl
import pytest
import pdb

a_chr1=     [ 1, 0, 1, 0, 0, 0, 4, 4, 2, 1]
t_chr1=     [ 0, 2, 3, 1, 0, 0, 1, 0, 0, 1]
c_chr1=     [10, 1, 2, 0, 0, 0, 1, 0, 0, 0]
g_chr1=     [ 0, 2, 0, 2, 3, 1, 0, 2, 7, 0]
a_chr2=     [ 0, 0, 4, 3, 2, 4, 0, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0]
a_chr2_1mut=[ 0, 0, 4, 3, 2, 4, 8, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0] # A single mutation at index (0->1)
t_chr2=     [ 1, 0, 0, 1, 0, 0, 7, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2]
t_chr2_1mut=[ 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2] # A single mutation at index (19->20)
c_chr2=     [ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
c_chr2_1mut=[ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0] # A single mutation at index (18->19)
g_chr2=     [ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
g_chr2_1mut=[ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
a_chr3=     [ 0, 0, 1, 0, 0, 0,15, 0, 0, 0, 0, 0, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_chr3=     [ 1, 0, 0, 1, 0,17, 1, 0,10,11, 0, 2, 3, 1, 3, 2,12,10, 1, 6, 1, 0, 0, 1, 0, 0, 8, 5, 0, 1]
c_chr3=     [ 0, 0,19,21,18, 1, 2,10, 0, 0,12, 1, 0, 0, 0,13, 0, 0,11, 0, 1, 0, 1, 2, 1,20, 0, 1, 0, 0]
g_chr3=     [20,20, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 2, 3, 1, 1, 2, 7, 0]

gene_locs = (["NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA","NA"])+ \
            (["NA","NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA",
              "NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA","NA"])+ \
            (["NA","NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA",
              "NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA","NA",
              "NA","NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA"])



@pytest.fixture(scope="module")
def profile_1()->pl.LazyFrame:
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    },schema={
        "chrom": pl.String,
        "pos": pl.Int64,
        "gene": pl.String,
        "A": pl.Int32,
        "T": pl.Int32,
        "C": pl.Int32,
        "G": pl.Int32,
    }).lazy()

@pytest.fixture(scope="module")
def scaffold_1()->str:
    return "chr1\nchr2\nchr3"

@pytest.fixture(scope="module")
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()

@pytest.fixture(scope="module")
def bed_table()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "start":[0,0,0],
        "end":[10,20,30]
    }).lazy()

@pytest.fixture(scope="module")
def genome_length_lf()->pl.LazyFrame:
    return pl.DataFrame({
        "genome":["genome1","genome2"],
        "genome_length":[30,30]
    }).lazy()

    
def test_get_genome_stats(stb,bed_table):
    result = utils.extract_genome_length(stb, bed_table).collect().rows_by_key("genome",unique=True,named=True)
    assert result["genome1"]["genome_length"] == 30
    assert result["genome2"]["genome_length"] == 30

def test_estimate_genome_presence_interface(profile_1,
                                  bed_table,
                                  stb,):
    
    
    read_loc_table = pl.LazyFrame(
        {
            "scaffold": ["chr1"]*7 + ["chr2"]*14,
            "loc": [0,1,2,3,6,7,9] + [0,2,3,6,7,11,12,13,14,15,16,17,18,19],
        }
    )
    
    result = utils.get_genome_stats(profile=profile_1,
                                    bed=bed_table,
                                    stb=stb,
                                    read_loc_table=read_loc_table
                                          ).collect()
    
    result_dict = result.rows_by_key("genome", unique=True, named=True)
    assert result_dict["genome1"]["is_present"] 
    assert result_dict["genome2"]["is_present"] 
    assert result_dict["genome1"]["coverage"] > 0
    assert result_dict["genome2"]["coverage"] > 0
    assert result_dict["genome1"]["breadth"] ==1
    assert result_dict["genome2"]["breadth"] ==1
    assert result_dict["genome1"]["fug"]>0
    assert result_dict["genome2"]["fug"]is None