"""
zipstrain.utils
========================
This module contains the command-line interface (CLI) implementation for the zipstrain application.
"""
import rich_click as click
import zipstrain.utils as ut
import zipstrain.compare as cp
import zipstrain.profile as pf
import zipstrain.task_manager as tm
import zipstrain.database as db
import polars as pl
import pathlib
import asyncio

@click.group()
def cli():
    """ZipStrain CLI"""
    pass

@cli.group()
def utilities():
    """The commands in this group are related to various utility functions that mainly prepare input files for profiling and comparison."""
    pass

@utilities.command("build-null-model")
@click.option('--error-rate', '-e', default=0.001, help="Error rate for the sequencing technology.")
@click.option('--max-total-reads', '-m', default=10000, help="Maximum coverage to consider for a base")
@click.option('--p-threshold', '-p', default=0.05, help="Significance threshold for the Poisson distribution.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
@click.option('--model-type', '-t', default="poisson", type=click.Choice(['poisson']), help="Type of null model to build.")
def build_null_model(error_rate, max_total_reads, p_threshold, output_file, model_type):
    """
    Build a null model for sequencing errors based on the Poisson distribution.

    Args:
    error_rate (float): Error rate for the sequencing technology.
    max_total_reads (int): Maximum total reads to consider.
    p_threshold (float): Significance threshold for the Poisson distribution.
    """
    if model_type == "poisson":
        df_thresh = ut.build_null_poisson(error_rate, max_total_reads, p_threshold)
    else:
        raise ValueError(f"Unsupported model type: {model_type}")
    df_thresh = pl.DataFrame(df_thresh, schema=["cov", "max_error_count"])
    df_thresh.write_parquet(output_file)

@utilities.command("merge_parquet")
@click.option('--input-dir', '-i', required=True, help="Directory containing Parquet files to merge.") 
@click.option('--output-file', '-o', required=True, help="Path to save the merged Parquet file.")
def merge_parquet(input_dir, output_file):
    """
    Merge multiple Parquet files in a directory into a single Parquet file, adding gene information.

    Args:
    input_dir (str): Directory containing Parquet files to merge.
    output_file (str): Path to save the merged Parquet file.
    """
    input_path = pathlib.Path(input_dir)
    parquet_files = list(input_path.glob("*.parquet"))
    if not parquet_files:
        raise ValueError(f"No Parquet files found in directory: {input_dir}")
    
    mpileup_df = pl.concat([pl.scan_parquet(pf) for pf in parquet_files])
    mpileup_df.sink_parquet(pathlib.Path(output_file), compression='zstd')


@utilities.command("process_mpileup")
@click.option('--gene-range-table-loc', '-g', required=True, help="Location of the gene range table in TSV format.")
@click.option('--batch-bed', '-b', required=True, help="Location of the batch BED file.")
@click.option('--batch-size', '-s', default=10000, help="Buffer size for processing stdin from samtools.")
@click.option('--output-file', '-o', required=True, help="Location to save the output Parquet file.")
def process_mpileup(gene_range_table_loc, batch_bed, batch_size, output_file):
    """
    Process mpileup files and save the results in a Parquet file.

    Args:
    gene_range_table_loc (str): Path to the gene range table in TSV format.
    batch_bed (str): Path to the batch BED file.
    output_file (str): Path to save the output Parquet file.
    """
    ut.process_mpileup_function(gene_range_table_loc, batch_bed, batch_size, output_file)
    
@utilities.command("make_bed")
@click.option('--db-fasta-dir', '-d', required=True, help="Path to the database in fasta format.")
@click.option('--max-scaffold-length', '-m', default=500000, help="Maximum scaffold length to split into multiple entries.")
@click.option('--output-file', '-o', required=True, help="Path to save the output BED file.")
def make_bed(db_fasta_dir, max_scaffold_length, output_file):
    """
    Create a BED file from the database in fasta format.

    Args:
    db_fasta_dir (str): Path to the fasta file.
    max_scaffold_length (int): Splits scaffolds longer than this into multiple entries of length <= max_scaffold_length.
    output_file (str): Path to save the output BED file.
    """
    bed_df = ut.make_the_bed(db_fasta_dir, max_scaffold_length)
    bed_df.write_csv(output_file, separator='\t', include_header=False)

@utilities.command("get_genome_lengths")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--bed-file', '-b', required=True, help="Path to the BED file.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
def get_genome_lengths(stb_file, bed_file, output_file):
    """
    Extract the genome length information from the scaffold-to-genome mapping table.

    Args:
    stb_file (str): Path to the scaffold-to-genome mapping file.
    bed_file (str): Path to the BED file containing genomic regions.
    output_file (str): Path to save the output Parquet file.
    """
    stb = pl.scan_csv(stb_file, separator='\t',has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").alias("genome")
    )
    
    bed_table = pl.scan_csv(bed_file, separator='\t',has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end")
    ).select(["scaffold", "start", "end"])
    genome_length = ut.extract_genome_length(stb, bed_table)
    genome_length.sink_parquet(output_file, compression='zstd')

@utilities.command("genome_breadth_matrix")
@click.option('--profile', '-p', type=str, required=True, help="Path to the profile Parquet file.")
@click.option('--genome-length', '-g', type=str, required=True, help="Path to the genome length Parquet file.")
@click.option('--stb', '-s', type=str, required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--min-cov', '-c', default=1, help="Minimum coverage to consider a position.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
def genome_breadth_matrix(profile, genome_length, stb, min_cov, output_file):
    """
    Generate a genome breadth matrix from the given profiles and scaffold-to-genome mapping.

    Args:
    profiles (list): List of profiles in the format 'name:path_to_profile'.
    genome_length (str): Path to the genome length Parquet file.
    stb_file (str): Path to the scaffold-to-genome mapping file.
    min_cov (int): Minimum coverage to consider a position.
    output_file (str): Path to save the output Parquet file.
    """
    genome_length = pl.scan_parquet(genome_length)
    stb = pl.scan_csv(stb, separator='\t', has_header=False).select(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").alias("genome")
    )
    profile_dir= pathlib.Path(profile)
    profile = pl.scan_parquet(profile)
    lf=ut.get_genome_breadth_matrix(profile,profile_dir.name, genome_length,stb, min_cov)
    lf.sink_parquet(output_file, compression='zstd')

@utilities.command("collect_breadth_tables")
@click.option('--breadth-tables-dir', '-d', required=True, help="Directory containing breadth tables in Parquet format.")
@click.option('--extension', '-e', default='parquet', help="File extension of the breadth tables.")
@click.option('--output-file', '-o', required=True, help="Path to save the collected breadth tables.")
def collect_breadth(breadth_tables_dir, extension, output_file):
    """
    Collect multiple genome breadth tables into a single LazyFrame.

    Args:
    breadth_tables_dir (str): Directory containing breadth tables in Parquet format.
    extension (str): File extension of the breadth tables.
    output_file (str): Path to save the collected breadth tables.
    """
    breadth_tables = list(pathlib.Path(breadth_tables_dir).glob(f"*.{extension}"))
    if not breadth_tables:
        raise ValueError(f"No breadth tables found in directory: {breadth_tables_dir}")

    lazy_frames = [pl.scan_parquet(str(pf)) for pf in breadth_tables]
    combined_lf = ut.collect_breadth_tables(lazy_frames)
    combined_lf.sink_parquet(output_file, compression='zstd')
    
@utilities.command("strain_heterogeneity")
@click.option('--profile-file', '-p', required=True, help="Path to the profile Parquet file.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--min-cov', '-c', default=5, help="Minimum coverage to consider a position.")
@click.option('--freq-threshold', '-f', default=0.8, help="Frequency threshold to define dominant nucleotide.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
def strain_heterogeneity(profile_file, stb_file, min_cov, freq_threshold, output_file):
    """
    Calculate strain heterogeneity for each genome based on nucleotide frequencies.

    Args:
    profile_file (str): Path to the profile Parquet file.
    stb_file (str): Path to the scaffold-to-genome mapping file.
    min_cov (int): Minimum coverage to consider a position.
    freq_threshold (float): Frequency threshold to define dominant nucleotide.
    output_file (str): Path to save the output Parquet file.
    """
    profile = pl.scan_parquet(profile_file)
    stb = pl.scan_csv(stb_file, separator="\t", has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").alias("genome")
    ).select(["scaffold", "genome"])
    
    het_profile = pf.get_strain_hetrogeneity(profile, stb, min_cov=min_cov, freq_threshold=freq_threshold)
    het_profile.sink_parquet(output_file, compression='zstd')

@utilities.command("build-profile-db")
@click.option('--profile-db-csv', '-p', required=True, help="Path to the profile database CSV file.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
def build_profile_db(profile_db_csv, output_file):
    """
    Build a profile database from the given CSV file.

    Args:
    profile_db_csv (str): Path to the profile database CSV file.
    """
    profile_db = db.ProfileDatabase.from_csv(pathlib.Path(profile_db_csv))
    profile_db.save_as_new_database(pathlib.Path(output_file))


@utilities.command("build-comparison-config")
@click.option('--profile-db', '-p', required=True, help="Path to the profile database Parquet file.")
@click.option('--gene-db-id', '-g', required=True, help="Gene database ID.")
@click.option('--reference-db-id', '-r', required=True, help="Reference fasta ID.")
@click.option('--scope', '-s', default="all", help="Genome scope for comparison.")
@click.option('--min-cov', '-c', default=5, help="Minimum coverage to consider a position.")
@click.option('--min-gene-compare-len', '-l', default=200, help="Minimum gene length to consider for comparison.")
@click.option('--null-model-p-value', '-n', default=0.05, help="P-value threshold for the null model to detect sequencing error.")
@click.option('--stb-file-loc', '-t', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--null-model-loc', '-m', required=True, help="Path to the null model Parquet file.")
@click.option('--current-comp-table', '-a', default=None, help="Path to the current comparison table in Parquet format.")
@click.option('--output-file', '-o', required=True, help="Path to save the output configuration JSON file.")    
def build_comparison_config(profile_db, gene_db_id, reference_genome_id, scope, min_cov, min_gene_compare_len, null_model_p_value, stb_file_loc, null_model_loc, current_comp_table, output_file):
    """
    Build a comparison configuration JSON file from the given parameters.

    Args:
    profile_db (str): Path to the profile database Parquet file.
    gene_db_id (str): Gene database ID.
    reference_genome_id (str): Reference genome ID.
    scope (str): Genome scope for comparison.
    min_cov (int): Minimum coverage to consider a position.
    min_gene_compare_len (int): Minimum gene length to consider for comparison.
    null_model_p_value (float): P-value threshold for the null model to detect sequencing error.
    stb_file_loc (str): Path to the scaffold-to-genome mapping file.
    null_model_loc (str): Path to the null model Parquet file.
    current_comp_table (str): Path to the current comparison table in Parquet format.
    output_file (str): Path to save the output configuration JSON file.
    """
    conf_obj=db.GenomeComparisonConfig(
        gene_db_id=gene_db_id,
        reference_id=reference_genome_id,
        scope=scope,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        null_model_p_value=null_model_p_value,
        stb_file_loc=stb_file_loc,
        null_model_loc=null_model_loc,
    )
    comp_obj=db.GenomeComparisonDatabase(
        profile_db=profile_db,
        config=conf_obj,
        comp_db_loc=current_comp_table
    )
    comp_obj.dump_obj(pathlib.Path(output_file))


@utilities.command("to-complete-table")
@click.option("--genome-comparison-object", "-g", required=True, help="Path to the genome comparison object in json format.")
@click.option("--output-file", "-o", required=True, help="Path to save the completed pairs csv file.")
def to_complete_table(genome_comparison_object, output_file):
    """
    Generate a table of completed genome comparison pairs and save it to a csv file.

    Parameters:
    genome_comparison_object (str): Path to the genome comparison object in json format.
    output_file (str): Path to save the completed pairs Parquet file.
    """
    genome_comp_db=db.GenomeComparisonDatabase.load_obj(pathlib.Path(genome_comparison_object))
    completed_pairs=genome_comp_db.to_complete_input_table()
    completed_pairs.sink_csv(pathlib.Path(output_file), compression='zstd', engine="streaming")


@cli.group()
def gene_tools():
    """Holds anything related to gene analysis."""
    pass


@gene_tools.command("gene-range-table")
@click.option('--gene-file', '-g', required=True, help="location of gene file. Prodigal's nucleotide fasta output")
@click.option('--output-file', '-o', required=True, help="location to save output tsv file")
def get_gene_range_table(gene_file, output_file):
    """
    Main function to build and save the gene location table.

    Args:
    gene_file (str): Path to the gene FASTA file.
    output_file (str): Path to save the output TSV file.
    """
    gene_locs=pf.build_gene_range_table(pathlib.Path(gene_file))
    gene_locs.sink_csv(pathlib.Path(output_file), separator="\t", include_header=False)


@gene_tools.command("gene-loc-table")
@click.option('--gene-file', '-g', required=True, help="location of gene file. Prodigal's nucleotide fasta output")
@click.option('--scaffold-list', '-s', required=True, help="location of scaffold list. A text file with each line containing a scaffold name.")
@click.option('--output-file', '-o', required=True, help="location to save output parquet file")
def get_gene_loc_table(gene_file, scaffold_list, output_file):
    """
    Main function to build and save the gene location table.

    Args:
    gene_file (str): Path to the gene FASTA file.
    scaffold_list (str): Path to the scaffold list file.
    output_file (str): Path to save the output Parquet file.
    """
    scaffolds=set(pl.read_csv(pathlib.Path(scaffold_list), has_header=False,separator="\t").select(pl.col("column_1")).to_series().to_list())
    gene_locs=pf.build_gene_loc_table(pathlib.Path(gene_file), scaffolds)
    gene_locs.write_parquet(pathlib.Path(output_file))



@cli.group()
def compare():
    """The commands in this group are related to comparing profiled samples."""
    pass

@compare.command("single_compare_genome")
@click.option('--mpileup-contig-1', '-m1', required=True, help="Path to the first mpileup file.")
@click.option('--mpileup-contig-2', '-m2', required=True, help="Path to the second mpileup file.")
@click.option('--scaffolds-1', '-s1', required=True, help="Path to the list of scaffolds for the first mpileup file.")
@click.option('--scaffolds-2', '-s2', required=True, help="Path to the list of scaffolds for the second mpileup file.")
@click.option('--null-model', '-n', required=True, help="Path to the null model Parquet file.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold to genome mapping file.")
@click.option('--min-cov', '-c', default=5, help="Minimum coverage to consider a position.")
@click.option('--min-gene-compare-len', '-l', default=100, help="Minimum gene length to consider for comparison.")
@click.option('--memory-mode', '-m', default="heavy", type=click.Choice(['heavy', 'light'], case_sensitive=False), help="Memory mode for processing: 'heavy' or 'light'.")
@click.option('--chrom-batch-size', '-b', default=10000, help="Batch size for processing chromosomes. Only used in light memory mode.")
@click.option('--output-file', '-o', required=True, help="Path to save the parquet file.")
@click.option('--genome', '-g', default="all", help="If provided, do the comparison only for the specified genome.")
@click.option('--engine', '-e', default="streaming", type=click.Choice(['streaming', 'gpu',"auto"], case_sensitive=False), help="Engine to use for processing: 'streaming', 'gpu' or 'auto'.")
def single_compare_genome(mpileup_contig_1, mpileup_contig_2, scaffolds_1, scaffolds_2, null_model, stb_file, min_cov, min_gene_compare_len, memory_mode, chrom_batch_size, output_file, genome, engine):
    """
    Main function to compare two mpileup files and calculate genome and gene statistics.
    
    Args:
    mpileup_contig_1 (str): Path to the first mpileup file.
    mpileup_contig_2 (str): Path to the second mpileup file.
    scaffolds_1 (str): Path to the list of scaffolds for the first mpileup file.
    scaffolds_2 (str): Path to the list of scaffolds for the second mpileup file.
    null_model (str): Path to the null model Parquet file.
    gene_locs (str): Path to the gene locations Parquet file.
    min_cov (int): Minimum coverage to consider a position.
    min_gene_compare_len (int): Minimum gene length to consider for comparison.
    memory_mode (str): Memory mode for processing: 'heavy' or 'light'.
    chrom_batch_size (int): Batch size for processing chromosomes. Only used in light memory mode.
    output_file (str): Path to save the parquet file.
    genome (str): If provided, do the comparison only for the specified genome.
    stb_file (str): Path to the scaffold to genome mapping file.
    """
    with pl.StringCache():
        mpile_contig_1 = pl.scan_parquet(mpileup_contig_1).with_columns(
            (pl.col("chrom").cast(pl.Categorical).alias("chrom"),
             pl.col("gene").cast(pl.Categorical).alias("gene"))
        )
        mpile_contig_2 = pl.scan_parquet(mpileup_contig_2).with_columns(
            (pl.col("chrom").cast(pl.Categorical).alias("chrom"),
             pl.col("gene").cast(pl.Categorical).alias("gene"))
        )

        stb = pl.scan_csv(stb_file, separator="\t", has_header=False).with_columns(
            pl.col("column_1").alias("scaffold").cast(pl.Categorical),
            pl.col("column_2").alias("genome").cast(pl.Categorical)
        ).select(["scaffold", "genome"])
        if genome != "all":
            stb = stb.filter(pl.col("genome") == genome)

    null_model = pl.scan_parquet(null_model)
    mpile_contig_1_name = pathlib.Path(mpileup_contig_1).name
    mpile_contig_2_name = pathlib.Path(mpileup_contig_2).name
    if genome != "all":
        scaffold_scope = stb.filter(pl.col("genome") == genome).collect()["scaffold"].to_list()
    else:
        scaffold_scope = None

    if memory_mode == "light":
        scaffolds_1 = pl.scan_csv(scaffolds_1, separator="\t", has_header=False).select(pl.col("column_1").alias("scaffold"))
        scaffolds_2 = pl.scan_csv(scaffolds_2, separator="\t", has_header=False).select(pl.col("column_1").alias("scaffold"))
        shared_scaffolds = list(set(scaffolds_1["scaffold"].to_list()).intersection(set(scaffolds_2["scaffold"].to_list())))
        mpile_contig_1 = mpile_contig_1.filter(pl.col("chrom").is_in(shared_scaffolds))
        mpile_contig_2 = mpile_contig_2.filter(pl.col("chrom").is_in(shared_scaffolds))
    else:
        shared_scaffolds=None
        
    
    comp = cp.compare_genomes(mpile_contig_1=mpile_contig_1, 
                     mpile_contig_2=mpile_contig_2, 
                     null_model=null_model,
                     scaffold_to_genome=stb, 
                     min_cov=min_cov,
                     min_gene_compare_len=min_gene_compare_len, 
                     memory_mode=memory_mode, 
                     chrom_batch_size=chrom_batch_size, 
                     shared_scaffolds=shared_scaffolds, 
                     scaffold_scope=scaffold_scope, 
                     engine=engine)
    comp=comp.join(
        stb.select("genome").unique(),
        left_on=pl.col("genome"),
        right_on=pl.col("genome"),
        how="full",
        coalesce=True
    ).fill_null(0)

    comp=comp.with_columns(pl.lit(mpile_contig_1_name).alias("sample_1"), pl.lit(mpile_contig_2_name).alias("sample_2")).fill_null(0)
    
    comp.sink_parquet(output_file,engine=engine)


@cli.group()
def profile():
    """The commands in this group are related to profiling bam files."""
    pass


@profile.command("prepare_profiling",help="Prepare the files needed for profiling bam files and save them in the specified output directory.")
@click.option('--reference-fasta', '-r', required=True, help="Path to the reference genome in FASTA format.")
@click.option('--gene-fasta', '-g', required=True, help="Path to the gene annotations in FASTA format.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--output-dir', '-o', required=True, help="Directory to save the profiling database.")
def prepare_profiling(reference_fasta, gene_fasta, stb_file, output_dir):
    """
    Prepare the files needed for profiling bam files and save them in the specified output directory.
    """
    output_dir=pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bed_df = ut.make_the_bed(reference_fasta)
    bed_df.write_csv(output_dir / "genomes_bed_file.bed", separator='\t', include_header=False)
    gene_range_table = pf.build_gene_range_table(pathlib.Path(gene_fasta))
    gene_range_table.write_csv(output_dir / "gene_range_table.tsv", separator='\t', include_header=False)
    
    stb = pl.scan_csv(stb_file, separator='\t',has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").alias("genome")
    )

    bed_df = bed_df.lazy()
    genome_length = ut.extract_genome_length(stb, bed_df)
    genome_length.sink_parquet(output_dir / "genome_lengths.parquet", compression='zstd')


@profile.command("profile-single")
@click.option('--bed-file', '-b', required=True, help="Path to the BED file describing regions to be profiled.")
@click.option('--bam-file', '-a', required=True, help="Path to the BAM file to be profiled.")
@click.option('--gene-range-table', '-g', required=True, help="Path to the gene range table.")
@click.option('--num-workers', '-n', default=1, help="Number of workers to use for profiling.")
@click.option('--output-dir', '-o', required=True, help="Directory to save the profiling output.")
def profile_single(bed_file, bam_file, gene_range_table, num_workers, output_dir):
    """
    Profile a single BAM file using the provided BED file and gene range table.
    
    """
    output_dir=pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pf.profile_bam(
        bed_file=bed_file,
        bam_file=bam_file,
        gene_range_table=gene_range_table,
        output_dir=output_dir,
        num_workers=num_workers
    )

@cli.group()
def run():
    """The commands in this group are related to running zipstrain workflows."""
    pass


@run.command("profile")
@click.option('--input-table', '-i', required=True, help="Path to the input table in TSV format containing sample names and paths to bam files.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--gene-range-table', '-g', required=True, help="Path to the gene range table file.")
@click.option('--bed-file', '-b', required=True, help="Path to the BED file for profiling regions.")
@click.option('--genome-length-file', '-l', required=True, help="Path to the genome length file.")
@click.option('--run-dir', '-r', required=True, help="Directory to save the run data.")
@click.option('--num-procs', '-n', default=8, help="Number of processors to use for each profiling task.")
@click.option('--max-concurrent-batches', '-m', default=5, help="Maximum number of concurrent batches to run.")
@click.option('--poll-interval', '-p', default=1, help="Polling interval in seconds to check the status of batches.")
@click.option('--execution-mode', '-e', default="local", help="Execution mode: 'local' or 'slurm'.")
@click.option('--slurm-config', '-c', default=None, help="Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.")
@click.option('--container-engine', '-o', default="local", help="Container engine to use: 'local', 'docker' or 'apptainer'.")
@click.option('--task-per-batch', '-t', default=10, help="Number of tasks to include in each batch.")
def profile(input_table, stb_file, gene_range_table, bed_file, genome_length_file, run_dir, num_procs, max_concurrent_batches, poll_interval, execution_mode, slurm_config, container_engine, task_per_batch):
    """
    Run BAM file profiling in batches using the specified execution mode and container engine.

    Args:
    input_table (str): Path to the input table in TSV format containing sample names and BAM file paths.
    stb_file (str): Path to the scaffold-to-genome mapping file.
    gene_range_table (str): Path to the gene range table file.
    bed_file (str): Path to the BED file for profiling regions.
    genome_length_file (str): Path to the genome length file.
    run_dir (str): Directory to save the run data.
    num_procs (int): Number of processors to use for each profiling task.
    max_concurrent_batches (int): Maximum number of concurrent batches to run.
    poll_interval (int): Polling interval in seconds to check the status of batches.
    execution_mode (str): Execution mode: 'local' or 'slurm'.
    slurm_config (str): Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.
    container_engine (str): Container engine to use: 'local', 'docker' or 'apptainer'.
    task_per_batch (int): Number of tasks to include in each batch.
    """
    # Load the BAM files table
    bams_lf = pl.scan_csv(input_table)
    
    # Validate required columns exist
    required_columns = {"sample_name", "bamfile"}
    actual_columns = set(bams_lf.collect_schema().names())
    if not required_columns.issubset(actual_columns):
        missing = required_columns - actual_columns
        raise ValueError(f"Input table missing required columns: {missing}")
    
    run_dir = pathlib.Path(run_dir)
    slurm_conf = None
    if execution_mode == "slurm":
        if slurm_config is None:
            raise ValueError("SLURM configuration file must be provided when execution mode is 'slurm'.")
        slurm_conf = tm.SlurmConfig.from_json(slurm_config)
    
    if container_engine == "local":
        container_engine_obj = tm.LocalEngine(address="")
    elif container_engine == "docker":
        container_engine_obj = tm.DockerEngine(address="parsaghadermazi/zipstrain:amd64")
    elif container_engine == "apptainer":
        container_engine_obj = tm.ApptainerEngine(address="docker://parsaghadermazi/zipstrain:amd64")
    else:
        raise ValueError("Invalid container engine. Choose from 'local', 'docker', or 'apptainer'.")
    
    tm.lazy_run_profile(
        run_dir=run_dir,
        container_engine=container_engine_obj,
        bams_lf=bams_lf,
        stb_file=pathlib.Path(stb_file),
        gene_range_table=pathlib.Path(gene_range_table),
        bed_file=pathlib.Path(bed_file),
        genome_length_file=pathlib.Path(genome_length_file),
        num_procs=num_procs,
        tasks_per_batch=task_per_batch,
        max_concurrent_batches=max_concurrent_batches,
        poll_interval=poll_interval,
        execution_mode=execution_mode,
        slurm_config=slurm_conf,
    )


@run.command("compare_genomes")
@click.option("--genome-comparison-object", "-g", required=True, help="Path to the genome comparison object in json format.")
@click.option("--run-dir", "-r", required=True, help="Directory to save the run data.")
@click.option("--max-concurrent-batches", "-m", default=5, help="Maximum number of concurrent batches to run.")
@click.option("--poll-interval", "-p", default=1, help="Polling interval in seconds to check the status of batches.")
@click.option("--execution-mode", "-e", default="local", help="Execution mode: 'local' or 'slurm'.")
@click.option("--slurm-config", "-s", default=None, help="Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.")
@click.option("--container-engine", "-c", default="local", help="Container engine to use: 'local', 'docker' or 'apptainer'.")
@click.option("--task-per-batch", "-t", default=10, help="Number of tasks to include in each batch.")
@click.option("--polars-engine", "-a", default="streaming", type=click.Choice(['streaming', 'gpu', 'auto'], case_sensitive=False), help="Polars engine to use: 'streaming', 'gpu' or 'auto'.")  
@click.option("--chrom-batch-size", "-b", default=10000, help="Batch size for processing chromosomes. Only used in light memory mode.")
@click.option("--memory-mode", "-h", default="heavy", type=click.Choice(['heavy', 'light'], case_sensitive=False), help="Memory mode for processing: 'heavy' or 'light'.")
def compare_genomes(genome_comparison_object, run_dir, max_concurrent_batches, poll_interval, execution_mode, slurm_config, container_engine, task_per_batch, polars_engine, chrom_batch_size, memory_mode):
    """
    Run genome comparisons in batches using the specified execution mode and container engine.

    Args:
    genome_comparison_object (str): Path to the genome comparison object in json format.
    run_dir (str): Directory to save the run data.
    max_concurrent_batches (int): Maximum number of concurrent batches to run.
    poll_interval (int): Polling interval in seconds to check the status of batches.
    execution_mode (str): Execution mode: 'local' or 'slurm'.
    slurm_config (str): Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.
    container_engine (str): Container engine to use: 'local', 'docker' or 'apptainer'.
    task_per_batch (int): Number of tasks to include in each batch.
    """
    genome_comp_db=db.GenomeComparisonDatabase.load_obj(pathlib.Path(genome_comparison_object))
    run_dir=pathlib.Path(run_dir)
    slurm_conf=None
    if execution_mode == "slurm":
        if slurm_config is None:
            raise ValueError("SLURM configuration file must be provided when execution mode is 'slurm'.")
        slurm_conf = tm.SlurmConfig.from_json(slurm_config)
    
    if container_engine == "local":
        container_engine_obj = tm.LocalEngine(address="")
    elif container_engine == "docker":
        container_engine_obj = tm.DockerEngine(address="parsaghadermazi/zipstrain:amd64") #could go to a toml or json config file
    elif container_engine == "apptainer":
        container_engine_obj = tm.ApptainerEngine(address="docker://parsaghadermazi/zipstrain:amd64") #could go to a toml or json config file
    else:
        raise ValueError("Invalid container engine. Choose from 'local', 'docker', or 'apptainer'.")
    tm.lazy_run_compares(
        comps_db=genome_comp_db,
        container_engine=container_engine_obj,
        run_dir=run_dir,
        max_concurrent_batches=max_concurrent_batches,
        polars_engine=polars_engine,
        execution_mode=execution_mode,
        slurm_config=slurm_conf,
        memory_mode=memory_mode,
        chrom_batch_size=chrom_batch_size,
        tasks_per_batch=task_per_batch,
        poll_interval=poll_interval,
    )



@run.command("build-comp-database")
@click.option("--profile-db-dir", "-p", required=True, help="Directory containing profile either in parquet format.")
@click.option("--config-file", "-c", required=True, help="Path to the  genome comparsion database config file in json format.")
@click.option("--output-dir", "-o", required=True, help="Directory to genome comparison database object.")
@click.option("--comp-db-file", "-f", required=False, help="The initial database file. If provided only additional comparisons will be added to this database.")
def build_comp_database(profile_db_dir, config_file, output_dir, comp_db_file):
    """
    Build a genome comparison database from the given profiles and configuration.

    Parameters:
    profile_db_dir (str): Directory containing profile either in parquet format.
    config_file (str): Path to the genome comparison database config file in json format.
    """
    profile_db_dir=pathlib.Path(profile_db_dir)
    profile_db=db.ProfileDatabase(
        db_loc=profile_db_dir,
    )
    existing_db_loc=pathlib.Path(comp_db_file) if comp_db_file is not None else None
    if existing_db_loc is not None and not existing_db_loc.exists():
        raise FileNotFoundError(f"{existing_db_loc} does not exist.")
    obj=db.GenomeComparisonDatabase(
        profile_db=profile_db,
        config=db.GenomeComparisonConfig.from_json(pathlib.Path(config_file)),
        comp_db_loc=existing_db_loc,
    )
    obj.dump_obj(pathlib.Path(output_dir))


        
        
@cli.command("test")
def test():
    """Run basic tests to ensure ZipStrain is setup correctly."""
    ### Check samtools installation
    if all([ut.check_samtools()]):
        click.echo("ZipStrain setup looks good!")
    else:
        click.echo("There are issues with the ZipStrain setup. Please check the above messages.")

if __name__ == "__main__":
    cli()