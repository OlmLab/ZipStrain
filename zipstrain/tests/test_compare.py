import pytest
import polars as pl
from zipstrain import compare
from click.testing import CliRunner
import json

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
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()

@pytest.fixture
def scaffold_1()->str:
    return "chr1\nchr2\nchr3"

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
def scaffold_2():
    return "chr1\nchr3"

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

@pytest.mark.parametrize("min_cov,min_gene_compare_len", [(1, 1), (5, 1), (1, 3), (5, 3),(5,5)])
def test_compare_profiles_profile_1_2_mc_mgcl(profile_1,profile_2,stb,null_model,min_cov,min_gene_compare_len):
    res_dict=compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        null_model=null_model,
        scaffold_to_genome=stb,
        min_cov=min_cov ,
        min_gene_compare_len=min_gene_compare_len
    ).collect().fill_null(-1).rows_by_key(key="genome",unique=True,named=True)

    a_genome_1=a_chr1 + a_chr2 
    t_genome_1=t_chr1 + t_chr2 
    c_genome_1=c_chr1 + c_chr2 
    g_genome_1=g_chr1 + g_chr2 

    assert res_dict["genome1"]["total_positions"]==len([i for i in zip(a_genome_1,t_genome_1,c_genome_1,g_genome_1) if sum(i)>=min_cov])

    a_genome_2=a_chr3
    t_genome_2=t_chr3
    c_genome_2=c_chr3
    g_genome_2=g_chr3

    assert res_dict["genome2"]["total_positions"]==len([i for i in zip(a_genome_2,t_genome_2,c_genome_2,g_genome_2) if sum(i)>=min_cov])

    assert res_dict["genome1"]["share_allele_pos"]==len([i for i in zip(a_genome_1,t_genome_1,c_genome_1,g_genome_1) if sum(i)>=min_cov])

    assert res_dict["genome2"]["share_allele_pos"]==len([i for i in zip(a_genome_2,t_genome_2,c_genome_2,g_genome_2) if sum(i)>=min_cov])

    assert res_dict["genome1"]["genome_pop_ani"]==100.0
    
    assert res_dict["genome2"]["genome_pop_ani"]==100.0

    assert res_dict["genome1"]["max_consecutive_length"]==max([len([i for i in zip(a_chr1,t_chr1,c_chr1,g_chr1) if sum(i)>=min_cov]),len([i for i in zip(a_chr2,t_chr2,c_chr2,g_chr2) if sum(i)>=min_cov])])

    assert res_dict["genome2"]["max_consecutive_length"]==max([len([i for i in zip(a_chr3,t_chr3,c_chr3,g_chr3) if sum(i)>=min_cov])])


    covered_gene_counts=profile_1.filter((pl.col("A")+pl.col("T")+pl.col("C")+pl.col("G")>=min_cov) & (pl.col("gene")!="NA")).join(stb,left_on="chrom",right_on="scaffold").group_by(["gene","genome"]).agg(pl.len()).filter(pl.col("len")>=min_gene_compare_len).group_by("genome").agg(pl.len()).collect().rows_by_key(key="genome",unique=True,named=True)
    
    
    if "genome1" in covered_gene_counts and "genome1" in res_dict:
        assert res_dict["genome1"]["shared_genes_count"]==covered_gene_counts["genome1"]["len"]
    
    if "genome2" in covered_gene_counts and "genome2" in res_dict:
        assert res_dict["genome2"]["shared_genes_count"]==covered_gene_counts["genome2"]["len"]

    if "genome1" in res_dict and "genome1" in covered_gene_counts:
        assert res_dict["genome1"]["identical_gene_count"]==covered_gene_counts["genome1"]["len"]

    if "genome2" in res_dict and "genome2" in covered_gene_counts:
        assert res_dict["genome2"]["identical_gene_count"]==covered_gene_counts["genome2"]["len"]

    if "genome1" in res_dict:
        assert res_dict["genome1"]["perc_id_genes"]==100.0 if "genome1" in covered_gene_counts else -1
    
    if "genome2" in res_dict:
        assert res_dict["genome2"]["perc_id_genes"]==100.0 if "genome2" in covered_gene_counts else -1




