"""zipstrain.compare
========================
This module provides all comparison functions for zipstrain.

"""
import polars as pl

def coverage_filter(mpile_frame:pl.LazyFrame, min_cov:int,engine:str)-> pl.LazyFrame:
    """
    Filter the mpile lazyframe based on minimum coverage at each loci.
    
    Parameters:
    mpile_frame (pl.LazyFrame): The input LazyFrame containing coverage data.
    min_cov (int): The minimum coverage threshold.
    
    Returns:
    pl.LazyFrame: Filtered LazyFrame with positions having coverage >= min_cov.
    """
    mpile_frame = mpile_frame.with_columns(
        (pl.col("A") + pl.col("C") + pl.col("G") + pl.col("T")).alias("cov")
    )
    return mpile_frame.filter(pl.col("cov") >= min_cov).collect(engine=engine).lazy()

def adjust_for_sequence_errors(mpile_frame:pl.LazyFrame, null_model:pl.LazyFrame) -> pl.LazyFrame:
    """
    Adjust the mpile frame for sequence errors based on the null model.
    
    Parameters:
    mpile_frame (pl.LazyFrame): The input LazyFrame containing coverage data.
    null_model (pl.LazyFrame): The null model LazyFrame containing error counts.
    
    Returns:
    pl.LazyFrame: Adjusted LazyFrame with sequence errors accounted for.
    """
    return mpile_frame.join(null_model, on="cov", how="left").with_columns([
        pl.when(pl.col(base) >= pl.col("max_error_count"))
        .then(pl.col(base))
        .otherwise(0)
        .alias(base)
        for base in ["A", "T", "C", "G"]
    ]).drop("max_error_count")

def get_shared_locs(mpile_contig_1:pl.LazyFrame, mpile_contig_2:pl.LazyFrame) -> pl.LazyFrame:
    """
    Retuerns a lazyframe with ATCG information for shared scaffolds and positions between two mpileup files.
    Parameters:
    mpile_contig_1 (pl.LazyFrame): The first mpileup LazyFrame.
    mpile_contig_2 (pl.LazyFrame): The second mpileup LazyFrame.
    Returns:
    pl.LazyFrame: Merged LazyFrame containing shared scaffolds and positions with ATCG information.
    """
    mpile_contig= mpile_contig_1.join(
        mpile_contig_2,
        on=["chrom", "pos"],
        how="inner",
        suffix="_2"  # To distinguish lf2 columns
    ).select(
        surr=pl.col("A") * pl.col("A_2") + pl.col("C") * pl.col("C_2") + pl.col("G") * pl.col("G_2") + pl.col("T") * pl.col("T_2"),
        scaffold=pl.col("chrom"),
        pos=pl.col("pos"),
        gene=pl.col("gene")
    )
    return mpile_contig

def add_contiguity_info(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """ Adds group id information to the lazy frame. If on the same scaffold and not popANI, then they are in the same group.
    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    Returns:
    pl.LazyFrame: Updated LazyFrame with group id information added.
    """

    mpile_contig= mpile_contig.sort(["scaffold", "pos"])
    mpile_contig = mpile_contig.with_columns([
        (pl.col("scaffold").shift(1).fill_null(pl.col("scaffold").first()).alias("prev_scaffold")),
    ])
    mpile_contig = mpile_contig.with_columns([
        (((pl.col("scaffold") != pl.col("prev_scaffold")) | (pl.col("surr") == 0))).cum_sum().alias("group_id")
    ])
    return mpile_contig

def add_genome_info(mpile_contig:pl.LazyFrame, scaffold_to_genome:pl.LazyFrame) -> pl.LazyFrame:
    """
    Adds genome information to the mpileup LazyFrame based on scaffold to genome mapping.
    
    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    scaffold_to_genome (pl.LazyFrame): The LazyFrame mapping scaffolds to genomes.
    
    Returns:
    pl.LazyFrame: Updated LazyFrame with genome information added.
    """
    return mpile_contig.join(
        scaffold_to_genome, on="scaffold", how="left"
    ).fill_null("NA")

def calculate_pop_ani(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """
    Calculates the population ANI (Average Nucleotide Identity) for the given mpileup LazyFrame.
    ###NOTE### Remember that this function should be applied to the merged mpileup using get_shared_locs.

    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
    pl.LazyFrame: Updated LazyFrame with population ANI information added.
    """
    return mpile_contig.group_by("genome").agg(
            total_positions=pl.len(),
            share_allele_pos=(pl.col("surr") > 0 ).sum()
        ).with_columns(
            genome_pop_ani=pl.col("share_allele_pos")/pl.col("total_positions")*100,
        )

def get_longest_consecutive_blocks(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """
    Calculates the longest consecutive blocks for each genome in the mpileup LazyFrame for any genome.
    
    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
    pl.LazyFrame: Updated LazyFrame with longest consecutive blocks information added.
    """
    block_lengths = (
        mpile_contig.group_by(["genome", "scaffold", "group_id"])
        .agg(pl.len().alias("length"))
    ) 
    return block_lengths.group_by("genome").agg(pl.col("length").max().alias("max_consecutive_length"))

def get_gene_ani(mpile_contig:pl.LazyFrame, min_gene_compare_len:int) -> pl.LazyFrame:
    """
    Calculates gene ANI (Average Nucleotide Identity) for each gene in each genome.
    
    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    min_gene_compare_len (int): Minimum length of the gene to consider for comparison.
    
    Returns:
    pl.LazyFrame: Updated LazyFrame with gene ANI information added.
    """
    return mpile_contig.group_by(["genome", "gene"]).agg(
        total_positions=pl.len(),
        share_allele_pos=(pl.col("surr") > 0).sum()
    ).filter(pl.col("total_positions") >= min_gene_compare_len).with_columns(
        identical=(pl.col("share_allele_pos") == pl.col("total_positions")),
    ).filter(pl.col("gene") != "NA").group_by("genome").agg(
        shared_genes_count=pl.len(),
        identical_gene_count=pl.col("identical").sum()
    ).with_columns(perc_id_genes=pl.col("identical_gene_count") / pl.col("shared_genes_count") * 100)

def get_unique_scaffolds(mpile_contig:pl.LazyFrame,batch_size:int=10000) -> set:
    """
    Retrieves unique scaffolds from the mpileup LazyFrame.
    
    Parameters:
    mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
    set: A set of unique scaffold names.
    """
    scaffolds = set()
    start_index = 0
    while True:
        batch = mpile_contig.slice(start_index, batch_size).select("chrom").fetch()
        if batch.height == 0:
            break
        scaffolds.update(batch["chrom"].to_list())
        start_index += batch_size
    return scaffolds 


def compare_genomes(mpile_contig_1:pl.LazyFrame,
              mpile_contig_2:pl.LazyFrame,
              null_model:pl.LazyFrame,
              scaffold_to_genome:pl.LazyFrame,
              min_cov:int=5,
              min_gene_compare_len:int=100,
              memory_mode:str="heavy",
              chrom_batch_size:int=10000,
              shared_scaffolds:list=None,
              scaffold_scope:list=None,
              engine="streaming"
            ):
    if memory_mode == "heavy":
        if scaffold_scope is not None:
            mpile_contig_1 = mpile_contig_1.filter(pl.col("chrom").is_in(scaffold_scope)).collect(engine=engine).lazy()
            mpile_contig_2 = mpile_contig_2.filter(pl.col("chrom").is_in(scaffold_scope)).collect(engine=engine).lazy()
        lf1=coverage_filter(mpile_contig_1, min_cov,engine=engine)
        lf1=adjust_for_sequence_errors(lf1, null_model)
        lf2=coverage_filter(mpile_contig_2, min_cov,engine=engine)
        lf2=adjust_for_sequence_errors(lf2, null_model)
        ### Now we need to only keep (scaffold, pos) that are in both lf1 and lf2
        lf = get_shared_locs(lf1, lf2)
        ## Add Contiguity Information
        lf = add_contiguity_info(lf)
        ## Let's add genome information for all scaffolds and positions
        lf = add_genome_info(lf, scaffold_to_genome)
        ## Let's calculate popANI
        genome_comp= calculate_pop_ani(lf)
        ## Calculate longest consecutive blocks
        max_consecutive_per_genome = get_longest_consecutive_blocks(lf)
        ## Calculate gene ani for each gene in each genome
        gene= get_gene_ani(lf, min_gene_compare_len)
        genome_comp=genome_comp.join(max_consecutive_per_genome, on="genome", how="left")
        genome_comp=genome_comp.join(gene, on="genome", how="left")
    
    elif memory_mode == "light":
        shared_scaffolds_batches = [shared_scaffolds[i:i + chrom_batch_size] for i in range(0, len(shared_scaffolds), chrom_batch_size)]
        lf_list=[]
        for scaffold in shared_scaffolds_batches:
            lf1= coverage_filter(mpile_contig_1.filter(pl.col("chrom").is_in(scaffold)), min_cov)
            lf1=adjust_for_sequence_errors(lf1, null_model)
            lf2= coverage_filter(mpile_contig_2.filter(pl.col("chrom").is_in(scaffold)), min_cov)
            lf2=adjust_for_sequence_errors(lf2, null_model)
            ### Now we need to only keep (scaffold, pos) that are in both lf1 and lf2
            lf = get_shared_locs(lf1, lf2)
            ## Lets add contiguity information
            lf= add_contiguity_info(lf)
            lf_list.append(lf)
        lf= pl.concat(lf_list)
        lf= add_genome_info(lf, scaffold_to_genome)
        genome_comp= calculate_pop_ani(lf)
        max_consecutive_per_genome = get_longest_consecutive_blocks(lf)
        gene= get_gene_ani(lf, min_gene_compare_len)
        genome_comp=genome_comp.join(max_consecutive_per_genome, on="genome", how="left")
        genome_comp=genome_comp.join(gene, on="genome", how="left")
    else:
        raise ValueError("Invalid memory_mode. Choose either 'heavy' or 'light'.")
    return genome_comp








