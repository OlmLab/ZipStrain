import random

import pytest
import polars as pl
from zipstrain import compare


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
        "genome": ["genome1"]*30 + ["genome2"]*30,
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
        "genome": ["genome1"]*30 + ["genome2"]*30,
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
def large_profile_1()->pl.LazyFrame:

    chr1=pl.DataFrame({
        "chrom": ["chr1"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr2=pl.DataFrame({
        "chrom": ["chr2"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr3=pl.DataFrame({
        "chrom": ["chr3"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    return pl.concat([chr1,chr2,chr3])
@pytest.fixture
def large_profile_2()->pl.LazyFrame:
    """A slightly different version of large_profile_1"""
    chr1=pl.DataFrame({
        "chrom": ["chr1"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr2=pl.DataFrame({
        "chrom": ["chr2"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr3=pl.DataFrame({
        "chrom": ["chr3"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    return pl.concat([chr1,chr2,chr3])
    
@pytest.mark.parametrize("min_cov,min_gene_compare_len", [(1, 1), (5, 1), (1, 3), (5, 3),(5,5)])
def test_compare_profiles_profile_1_2_mc_mgcl(profile_1,profile_2,min_cov,min_gene_compare_len):
    res_dict=compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=min_cov ,
        min_gene_compare_len=min_gene_compare_len,
        calculate="all",
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


    covered_gene_counts=profile_1.filter((pl.col("A")+pl.col("T")+pl.col("C")+pl.col("G")>=min_cov) & (pl.col("gene")!="NA")).group_by(["gene","genome"]).agg(pl.len()).filter(pl.col("len")>=min_gene_compare_len).group_by("genome").agg(pl.len()).collect().rows_by_key(key="genome",unique=True,named=True)
    
    
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

def test_con_ani_expression():
    profile_1=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [10,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [0,0,5],
    }).lazy()
    profile_2=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [5,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [5,0,5],
    }).lazy()
    profile_3=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [0,0,0],
        "T": [10,10,0],
        "C": [0,0,0],
        "G": [0,0,10],
    }).lazy()

    assert compare.get_shared_locs(profile_1, profile_2, ani_method="conani").select(pl.col("surr")).sum().collect()[0,0]==3
    assert compare.get_shared_locs(profile_3, profile_2, ani_method="conani").select(pl.col("surr")).sum().collect()[0,0]==2

@pytest.mark.parametrize("threshold1,threshold2", [(0.1,0.8),(0.2,0.9),(0.3,0.95),(0.4,0.99),(0.1,0.2),(0.5,0.6),(0.7,0.8),(0,1)])
def test_cos_ani_expression(threshold1,threshold2):
    profile_1=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [10,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [0,0,5],
    }).lazy()
    profile_2=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [5,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [5,0,5],
    }).lazy()
    profile_3=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [0,0,0],
        "T": [10,10,0],
        "C": [0,0,0],
        "G": [0,0,10],
    }).lazy()

    assert compare.get_shared_locs(profile_1, profile_2, ani_method=f"cosani_{threshold1}").select(pl.col("surr")).sum().collect()[0,0]>=compare.get_shared_locs(profile_3, profile_2, ani_method=f"cosani_{threshold2}").select(pl.col("surr")).sum().collect()[0,0]


def test_duckdb_compare_genomes_to_parquet(profile_1, profile_2, stb, tmp_path):
    p1 = tmp_path / "p1.parquet"
    p2 = tmp_path / "p2.parquet"
    stb_path = tmp_path / "stb.tsv"
    out_path = tmp_path / "genome_comp.parquet"
    profile_1.sink_parquet(p1)
    profile_2.sink_parquet(p2)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    compare.duckdb_compare_genomes_to_parquet(
        mpile1=p1,
        mpile2=p2,
        output_file=out_path,
        stb_file=stb_path,
        sample_1_name="s1",
        sample_2_name="s2",
        min_cov=5,
        min_gene_compare_len=3,
        genome_scope="all",
        ani_method="popani",
        calculate="all",
        memory_limit="512MB",
        temp_directory=tmp_path,
    )

    out = pl.read_parquet(out_path)
    assert out.shape[0] == 2
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "max_consecutive_length",
        "shared_genes_count",
        "identical_gene_count",
        "perc_id_genes",
        "sample_1",
        "sample_2",
    }
    assert set(out["sample_1"].to_list()) == {"s1"}
    assert set(out["sample_2"].to_list()) == {"s2"}


def test_duckdb_compare_genomes_to_parquet_ani_only(profile_1, profile_2, stb, tmp_path):
    p1 = tmp_path / "p1.parquet"
    p2 = tmp_path / "p2.parquet"
    stb_path = tmp_path / "stb.tsv"
    out_path = tmp_path / "genome_comp_ani.parquet"
    profile_1.sink_parquet(p1)
    profile_2.sink_parquet(p2)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    compare.duckdb_compare_genomes_to_parquet(
        mpile1=p1,
        mpile2=p2,
        output_file=out_path,
        stb_file=stb_path,
        sample_1_name="s1",
        sample_2_name="s2",
        min_cov=5,
        min_gene_compare_len=3,
        genome_scope="all",
        ani_method="popani",
        calculate="ani",
        memory_limit="512MB",
        temp_directory=tmp_path,
    )

    out = pl.read_parquet(out_path)
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "sample_1",
        "sample_2",
    }


def test_compare_genomes_ani_only_matches_full_subset(profile_1, profile_2, stb, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    full = (
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            stb_file=stb_path,
            calculate="all",
        )
        .collect()
        .select("genome", "total_positions", "share_allele_pos", "genome_pop_ani")
        .sort("genome")
    )
    ani_only = (
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            stb_file=stb_path,
            calculate="ani",
        )
        .collect()
        .sort("genome")
    )
    assert ani_only.equals(full)


def test_compare_genomes_polars_with_categorical_keys(tmp_path):
    p1 = (
        pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2"],
                "genome": ["genome1", "genome1", "genome2", "genome2"],
                "pos": [0, 1, 0, 1],
                "gene": ["geneA", "geneA", "geneB", "geneB"],
                "A": [10, 5, 8, 8],
                "T": [0, 0, 0, 0],
                "C": [0, 0, 0, 0],
                "G": [0, 0, 0, 0],
            }
        )
        .with_columns(
            pl.col("chrom").cast(pl.Categorical),
            pl.col("genome").cast(pl.Categorical),
            pl.col("gene").cast(pl.Categorical),
        )
        .lazy()
    )
    p2 = (
        pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2"],
                "genome": ["genome1", "genome1", "genome2", "genome2"],
                "pos": [0, 1, 0, 1],
                "gene": ["geneA", "geneA", "geneB", "geneB"],
                "A": [9, 4, 7, 7],
                "T": [0, 0, 0, 0],
                "C": [0, 0, 0, 0],
                "G": [0, 0, 0, 0],
            }
        )
        .with_columns(
            pl.col("chrom").cast(pl.Categorical),
            pl.col("genome").cast(pl.Categorical),
            pl.col("gene").cast(pl.Categorical),
        )
        .lazy()
    )
    stb_path = tmp_path / "stb.tsv"
    pl.DataFrame({"scaffold": ["chr1", "chr2"], "genome": ["genome1", "genome2"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )

    out = compare.compare_genomes(
        mpile_contig_1=p1,
        mpile_contig_2=p2,
        min_cov=1,
        min_gene_compare_len=1,
        engine="polars",
        stb_file=stb_path,
    ).collect()

    assert out.shape[0] == 2
    assert set(out["genome"].to_list()) == {"genome1", "genome2"}


def test_compare_genomes_default_calculate_is_all_metrics(profile_1, profile_2, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    pl.DataFrame({"scaffold": ["chr1", "chr2", "chr3"], "genome": ["genome1", "genome1", "genome2"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )
    out = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="polars",
        stb_file=stb_path,
    ).collect()
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "max_consecutive_length",
        "shared_genes_count",
        "identical_gene_count",
        "perc_id_genes",
    }


def test_compare_genomes_calculate_all_matches_between_engines(profile_1, profile_2, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    pl.DataFrame({"scaffold": ["chr1", "chr2", "chr3"], "genome": ["genome1", "genome1", "genome2"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )
    out_polars = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="polars",
        calculate="all",
        stb_file=stb_path,
    ).collect().sort("genome")
    out_duckdb = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="duckdb",
        calculate="all",
        stb_file=stb_path,
    ).collect().sort("genome")
    assert out_polars.equals(out_duckdb)
