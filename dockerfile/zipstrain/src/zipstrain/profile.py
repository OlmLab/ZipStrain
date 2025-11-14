"""zipstrain.profile
========================
This module provides functions and utilities to profile a bamfile.
By profile we mean generating gene, genome, and nucleotide counts at each position on the reference.
This is a fundamental step for downstream analysis in zipstrain.
"""
import pathlib
import polars as pl
from typing import Generator


def parse_gene_loc_table(fasta_file:pathlib.Path) -> Generator[tuple,None,None]:
    """
    Extract gene locations from a FASTA assuming it is from prodigal yield gene info.

    Parameters:
    fasta_file (pathlib.Path): Path to the FASTA file.

    Returns:
    Tuple: A tuple containing:
        - gene_ID
        - scaffold
        - start
        - end
    """
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split()
                gene_id = parts[0]
                scaffold = "_".join(gene_id.split('_')[:-1])
                start = parts[2]
                end=parts[4]      
                yield gene_id, scaffold,start,end


def build_gene_loc_table(fasta_file:pathlib.Path,scaffold:set)->pl.DataFrame:
    """
    Build a gene location table from a FASTA file.

    Parameters:
    fasta_file (pathlib.Path): Path to the FASTA file.

    Returns:
    pl.DataFrame: A Polars DataFrame containing gene locations.
    """
    scaffolds = []
    gene_ids = []
    pos=[]
    for genes in parse_gene_loc_table(fasta_file):
        if genes[1] in scaffold:
            scaffolds.extend([genes[1]]* (int(genes[3])-int(genes[2])+1))
            gene_ids.extend([genes[0]]* (int(genes[3])-int(genes[2])+1))
            pos.extend(list(range(int(genes[2]), int(genes[3])+1)))
    return pl.DataFrame({
        "scaffold":scaffolds,
        "gene":gene_ids,
        "pos":pos
    })
    
def build_gene_range_table(fasta_file:pathlib.Path)->pl.DataFrame:
    """
    Build a gene location table in the form of <gene scaffold start end> from a FASTA file.
    Parameters:
    fasta_file (pathlib.Path): Path to the FASTA file.

    Returns:
    pl.DataFrame: A Polars DataFrame containing gene locations.
    """
    out=[]
    for parsed_annot in parse_gene_loc_table(fasta_file):
        out.append(parsed_annot)
    return pl.DataFrame(out, schema=["gene", "scaffold", "start", "end"],orient='row')



def add_gene_info_to_mpileup(mpileup_df:pl.LazyFrame, gene_range:pl.DataFrame)->pl.DataFrame:
    mpileup_df=mpileup_df.with_columns(pl.col("gene").fill_null("NA"))
    for gene, scaffold, start, end in gene_range.iter_rows():
        mpileup_df=mpileup_df.with_columns(
            pl.when((pl.col("chrom") == scaffold) & (pl.col("pos") >= start) & (pl.col("pos") <= end))
            .then(gene)
            .otherwise(pl.col("gene"))
            .alias("gene")
        )
    return mpileup_df


def get_strain_hetrogeneity(profile:pl.LazyFrame,
                            stb:pl.LazyFrame, 
                            min_cov=5,
                            freq_threshold=0.8)->pl.LazyFrame:
    """
    Calculate strain heterogeneity for each genome based on nucleotide frequencies.
    The definition of strain heterogeneity here is the fraction of sites that have enough coverage
    (min_cov) and have a dominant nucleotide with frequency less than freq_threshold.

    Args:
        profile (pl.LazyFrame): The profile LazyFrame containing nucleotide counts.
        stb (pl.LazyFrame): The scaffold-to-bin mapping LazyFrame. First column is 'scaffold', second column is 'bin'.
        min_cov (int): The minimum coverage threshold.
        freq_threshold (float): The frequency threshold for dominant nucleotides.

    Returns:
    pl.LazyFrame: A LazyFrame containing strain heterogeneity information grouped by genome.
    """
    # Calculate the total number of sites with sufficient coverage
    profile = profile.with_columns(
        (pl.col("A")+pl.col("T")+pl.col("C")+pl.col("G")).alias("coverage")
    ).filter(pl.col("coverage") >= min_cov)
    
    profile = profile.with_columns(
        (pl.max_horizontal(["A", "T", "C", "G"])/pl.col("coverage") < freq_threshold)
        .cast(pl.Int8)
        .alias("heterogeneous_site")
    )
    
    profile = profile.join(stb, left_on="chrom", right_on="scaffold", how="left").group_by("genome").agg([
        pl.len().alias(f"total_sites_at_{min_cov}_coverage"),
        pl.sum("heterogeneous_site").alias("heterogeneous_sites")
    ])
    
    strain_heterogeneity = profile.with_columns(
        (pl.col("heterogeneous_sites")/pl.col(f"total_sites_at_{min_cov}_coverage")).alias("strain_heterogeneity")
    )
    return strain_heterogeneity