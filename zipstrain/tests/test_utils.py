from zipstrain import database,compare,utils,cli
import polars as pl
import numpy as np
from click.testing import CliRunner
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
        "genome": ["genome1"]*10+["genome1"]*20+["genome2"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    },schema={
        "chrom": pl.String,
        "pos": pl.Int64,
        "genome": pl.String,
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
def gene_bed()->pl.LazyFrame:
    return pl.DataFrame({
        "gene": ["gene1", "gene2", "gene3", "gene1", "gene2", "gene3"],
        "scaffold": ["chr1", "chr2", "chr2", "chr3", "chr3", "chr3"],
        "start": [2, 2, 11, 3, 11, 22],
        "end": [5, 6, 15, 6, 15, 26],
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


def test_get_gene_stats(profile_1, gene_bed, stb):
    result = utils.get_gene_stats(
        profile=profile_1,
        gene_bed=gene_bed,
        stb=stb,
    ).collect()

    assert result.columns == ["genome", "gene", "length", "breadth", "coverage"]
    g1 = result.filter((pl.col("genome") == "genome1") & (pl.col("gene") == "gene1")).to_dicts()[0]
    g2 = result.filter((pl.col("genome") == "genome2") & (pl.col("gene") == "gene1")).to_dicts()[0]
    assert g1["length"] == 4
    assert g1["breadth"] == pytest.approx(1.0)
    assert g1["coverage"] == pytest.approx(13 / 4)
    assert g2["length"] == 4
    assert g2["breadth"] == pytest.approx(1.0)
    assert g2["coverage"] == pytest.approx(76 / 4)

def test_estimate_genome_presence_interface(profile_1,
                                  bed_table,
                                  stb,):
    
    # This test only checks to see if function provides expected columns. Otherwise, the logic of the test case is wrong because zero reads are mapped to genome2 while
    #breadth and coverage are both > 0. Another test below checks with reasonable read mapping.
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
    assert result_dict["genome1"]["coverage"] > 0
    assert result_dict["genome2"]["coverage"] > 0
    assert result_dict["genome1"]["breadth"] ==1
    assert result_dict["genome2"]["breadth"] ==1
    assert result_dict["genome1"]["fug"]>0
    assert result_dict["genome2"]["fug"]is None
    assert result_dict["genome1"]["reads_mapped"]==21
    assert result_dict["genome2"]["reads_mapped"]==0
    assert result_dict["genome1"]["ber"]>0
    assert result_dict["genome2"]["ber"]>0

def test_estimate_genome_presence_makes_sense():
    # This test case will include two genomes of about real lenght expected from a bacterial genome. Each genome will have a few scaffolds. we will build each case here:
    np.random.seed(10)
    # specs:
    # read length: 150bp
    read_length = 150
    #genome 1:
    genome1_scaffold_lengths = [500000, 300000, 200000]  # Total length: 10000bp

    # We design genome 1 to exist and the reads are well distributed across the genome
    # Number of reads mapped to genome 1
    num_reads_genome1 = 10000
    reads_per_scaffold_genome1 = [5_000, 3_000, 2_000]
    bed_table=pl.LazyFrame({
        "scaffold": [f"genome1_scaf{i+1}" for i in range(len(genome1_scaffold_lengths))],
        "start": [0]*len(genome1_scaffold_lengths),
        "end": genome1_scaffold_lengths
    }) 
    read_loc_table_1=pl.concat([
        pl.DataFrame({
            "scaffold": [f"genome1_scaf{i+1}"]*reads_per_scaffold_genome1[i],
            "loc": np.random.randint(0, genome1_scaffold_lengths[i]-read_length, size=reads_per_scaffold_genome1[i]).tolist()
        }) for i in range(len(genome1_scaffold_lengths))
    ]).lazy()
    # we need a profile table that has chrom,pos,A,T,C,G columns for all positions that are covered by at least one read. We can get the coverage profile since we know where the reads are mapped.
    # Since we do not care about the actual nucleotides here, we can just set A=coverage and T=C=G=0 for simplicity.
    cov=[np.zeros(length,dtype=int) for length in genome1_scaffold_lengths]
    for i in range(len(genome1_scaffold_lengths)):
        scaffold_reads = read_loc_table_1.filter(pl.col("scaffold")==f"genome1_scaf{i+1}").collect()
        for loc in scaffold_reads["loc"]:
            cov[i][loc:loc+read_length] += 1
    
    profile_1=pl.concat([
        pl.DataFrame({
            "chrom": [f"genome1_scaf{i+1}"]*genome1_scaffold_lengths[i],
            "genome": ["genome1"]*genome1_scaffold_lengths[i],
            "pos": list(range(genome1_scaffold_lengths[i])),
            "gene": ["NA"]*genome1_scaffold_lengths[i],
            "A": cov[i].tolist(),
            "T": [0]*genome1_scaffold_lengths[i],
            "C": [0]*genome1_scaffold_lengths[i],
            "G": [0]*genome1_scaffold_lengths[i],
        }) for i in range(len(genome1_scaffold_lengths))
    ]).lazy().filter(pl.col("A")>0)
    
    # Now we want to add genome 2 that does not exist in the sample. However, some reads are mapped to it because of similarity in only a small region with the reads.
    # At the same time we'd like to test for some edge cases. As a result we design genome 2 such that:
    # - it has one scaffold with zero reads mapped to it
    # - it has one scaffold with just one read mapped to it
    # - it has one scaffold with a few reads mapped to it but all concentrated in two small regions leading to low breadth and high unexpected gap.
    genome2_scaffold_lengths = [400000, 250000, 150000]  # Total length: 800000bp
    num_reads_genome2 = 100
    reads_per_scaffold_genome2 = [0, 1, 99]
    bed_table_2=pl.LazyFrame({
        "scaffold": [f"genome2_scaf{i+1}" for i in range(len(genome2_scaffold_lengths))],
        "start": [0]*len(genome2_scaffold_lengths),
        "end": genome2_scaffold_lengths
    })
    read_locs={
        "genome2_scaf1": [],
        "genome2_scaf2": [50000],
        "genome2_scaf3": np.random.choice([1000+np.random.randint(0, 100), 50000+np.random.randint(0, 100), 120000+np.random.randint(0, 100)], size=99, replace=True).tolist()
    }
    reads_2_dict={
        "scaffold": [],
        "loc": []
    }
    for scaf in read_locs:
        if len(read_locs[scaf])>0:
            reads_2_dict["scaffold"] += [scaf]*len(read_locs[scaf])
            reads_2_dict["loc"] += read_locs[scaf]
    read_loc_table_2=pl.DataFrame(reads_2_dict).lazy()
    
    # we need a profile table that has chrom,pos,A,T,C,G columns for all positions that are covered by at least one read. We can get the coverage profile since we know where the reads are mapped.
    # Since we do not care about the actual nucleotides here, we can just set A=coverage and T=C=G=0 for simplicity.
    cov2=[np.zeros(length,dtype=int) for length in genome2_scaffold_lengths]
    for i in range(len(genome2_scaffold_lengths)):
        scaffold_reads = read_loc_table_2.filter(pl.col("scaffold")==f"genome2_scaf{i+1}").collect()
        for loc in scaffold_reads["loc"]:
            cov2[i][loc:loc+read_length] += 1   
    profile_2=pl.concat([
        pl.DataFrame({
            "chrom": [f"genome2_scaf{i+1}"]*genome2_scaffold_lengths[i],
            "genome": ["genome2"]*genome2_scaffold_lengths[i],
            "pos": list(range(genome2_scaffold_lengths[i])),
            "gene": ["NA"]*genome2_scaffold_lengths[i],
            "A": cov2[i].tolist(),
            "T": [0]*genome2_scaffold_lengths[i],
            "C": [0]*genome2_scaffold_lengths[i],
            "G": [0]*genome2_scaffold_lengths[i],
        }) for i in range(len(genome2_scaffold_lengths))
    ]).lazy().filter(pl.col("A")>0)    
    full_profile=pl.concat([profile_1,profile_2])
    full_stb=pl.LazyFrame({
        "scaffold":[f"genome1_scaf{i+1}" for i in range(len(genome1_scaffold_lengths))]+
                   [f"genome2_scaf{i+1}" for i in range(len(genome2_scaffold_lengths))],
        "genome":["genome1"]*len(genome1_scaffold_lengths)+
                 ["genome2"]*len(genome2_scaffold_lengths),
    }).lazy()
    full_bed=pl.concat([bed_table,bed_table_2])
    full_read_loc_table=pl.concat([read_loc_table_1,read_loc_table_2])
    result = utils.get_genome_stats(profile=full_profile,
                                    bed=full_bed,
                                    stb=full_stb,
                                    read_loc_table=full_read_loc_table
                                          ).collect()
    genome_lengths={
        "genome1": sum(genome1_scaffold_lengths),
        "genome2": sum(genome2_scaffold_lengths)
    }
    result_dict = result.rows_by_key("genome", unique=True, named=True)
    assert result_dict["genome1"]["reads_mapped"]==num_reads_genome1   
    assert result_dict["genome2"]["reads_mapped"]==num_reads_genome2
    assert result_dict["genome1"]["coverage"] == (num_reads_genome1 * read_length) / genome_lengths["genome1"]
    assert result_dict["genome2"]["coverage"] == (num_reads_genome2 * read_length) / genome_lengths["genome2"]
    assert result_dict["genome1"]["ber"] > result_dict["genome2"]["ber"] 
