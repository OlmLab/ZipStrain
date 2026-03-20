"""zipstrain.profile
========================
This module provides functions and utilities to profile a bamfile.
By profile we mean generating gene, genome, and nucleotide counts at each position on the reference.
This is a fundamental step for downstream analysis in zipstrain.
"""
import pathlib
import polars as pl
from typing import Generator
from zipstrain import utils
import asyncio
import os
import duckdb

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

def adjust_for_sequence_errors(mpile_frame:pl.LazyFrame, null_model:pl.LazyFrame) -> pl.LazyFrame:
    """
    Adjust the mpile frame for sequence errors based on the null model.
    
    Args:
        mpile_frame (pl.LazyFrame): The input LazyFrame containing coverage data.
        null_model (pl.LazyFrame): The null model LazyFrame containing error counts.
    
    Returns:
        pl.LazyFrame: Adjusted LazyFrame with sequence errors accounted for.
    """
    mpile_frame = mpile_frame.with_columns(
        pl.sum_horizontal(["A", "T", "C", "G"]).alias("cov")
    )
    return mpile_frame.join(null_model, on="cov", how="left").with_columns([
        pl.when(pl.col(base) >= pl.col("max_error_count"))
        .then(pl.col(base))
        .otherwise(0)
        .alias(base)
        for base in ["A", "T", "C", "G"]
    ]).drop("max_error_count")
    
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

def add_genome_info_to_mpileup(mpileup_df:pl.LazyFrame, scaffold_to_genome:pl.LazyFrame)->pl.LazyFrame:
    mpileup_df=mpileup_df.join(scaffold_to_genome,
                               left_on="chrom", right_on="scaffold", how="left").with_columns(
        pl.col("genome").fill_null("NA")    )
    return mpileup_df

def add_gene_info_to_mpileup(mpileup_df:pl.LazyFrame, gene_range:pl.LazyFrame)->pl.LazyFrame:
    mpileup_df=mpileup_df.sort(["chrom", "pos"])
    gene_range=gene_range.sort(["scaffold", "start"])
    annotated_mpileup=mpileup_df.join_asof(
        gene_range,
        left_on="pos",
        right_on="start",
        by_left="chrom",
        by_right="scaffold",
        strategy="backward",
        check_sortedness=False,
    ).with_columns(
            pl.when(pl.col("pos") <= pl.col("end"))
            .then(pl.col("gene"))
            .otherwise(pl.lit("NA"))
            .alias("gene")
        )
    return annotated_mpileup


def _duckdb_quote_sql_string(value: str) -> str:
    """Quote a string literal for embedding in DuckDB SQL."""
    return value.replace("'", "''")


def _annotate_mpileup_chunk_with_duckdb(
    adjusted_mpileup_parquet: pathlib.Path,
    output_parquet: pathlib.Path,
    scaffold_to_genome: pl.LazyFrame,
    gene_range: pl.LazyFrame,
) -> None:
    """Annotate one adjusted mpileup chunk with genome+gene in DuckDB."""
    conn = duckdb.connect()
    try:
        conn.register(
            "stb_src",
            scaffold_to_genome.select(["scaffold", "genome"]).collect().to_arrow(),
        )
        conn.register(
            "gene_src",
            gene_range.select(["gene", "scaffold", "start", "end"]).collect().to_arrow(),
        )
        in_sql = _duckdb_quote_sql_string(str(adjusted_mpileup_parquet))
        out_sql = _duckdb_quote_sql_string(str(output_parquet))
        conn.execute(
            f"""
            COPY (
              WITH mpileup AS (
                SELECT
                  CAST(chrom AS VARCHAR) AS chrom,
                  CAST(pos AS INTEGER) AS pos,
                  CAST(A AS INTEGER) AS A,
                  CAST(C AS INTEGER) AS C,
                  CAST(G AS INTEGER) AS G,
                  CAST(T AS INTEGER) AS T
                FROM read_parquet('{in_sql}')
              ),
              stb AS (
                SELECT
                  CAST(scaffold AS VARCHAR) AS scaffold,
                  CAST(genome AS VARCHAR) AS genome
                FROM stb_src
              ),
              gene_range AS (
                SELECT
                  CAST(gene AS VARCHAR) AS gene,
                  CAST(scaffold AS VARCHAR) AS scaffold,
                  CAST(start AS INTEGER) AS start,
                  CAST("end" AS INTEGER) AS "end"
                FROM gene_src
                ORDER BY scaffold, start
              ),
              mp_with_genome AS (
                SELECT
                  mp.chrom,
                  COALESCE(stb.genome, 'NA') AS genome,
                  mp.pos,
                  mp.A,
                  mp.C,
                  mp.G,
                  mp.T
                FROM mpileup AS mp
                LEFT JOIN stb ON mp.chrom = stb.scaffold
              ),
              annotated AS (
                SELECT
                  mg.chrom,
                  mg.genome,
                  CASE
                    WHEN mg.pos <= gr."end" THEN gr.gene
                    ELSE 'NA'
                  END AS gene,
                  mg.pos,
                  mg.A,
                  mg.C,
                  mg.G,
                  mg.T
                FROM mp_with_genome AS mg
                ASOF LEFT JOIN gene_range AS gr
                  ON mg.chrom = gr.scaffold
                 AND mg.pos >= gr.start
              )
              SELECT
                chrom,
                genome,
                gene,
                pos,
                A,
                C,
                G,
                T
              FROM annotated
            ) TO '{out_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        conn.close()


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



async def _profile_chunk_task(
    bed_file:pathlib.Path,
    bam_file:pathlib.Path,
    gene_range_table:pathlib.Path,
    stb:pl.LazyFrame,
    null_model:pl.LazyFrame,
    output_dir:pathlib.Path,
    chunk_id:int
)->None:
    cmd=["samtools", "mpileup", "-A", "-l", str(bed_file.absolute()), str(bam_file.absolute())]
    cmd += ["|", "zipstrain", "utilities", "process_mpileup", "--output-file", f"{bam_file.stem}_{chunk_id}_pre.parquet"]
    proc = await asyncio.create_subprocess_shell(
                " ".join(cmd),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=output_dir
            )
    stdout, stderr = await proc.communicate()
    if proc.returncode != 0:
        raise Exception(f"Command failed with error: {stderr.decode().strip()}")
    else:
        mpileup=pl.scan_parquet(output_dir/f"{bam_file.stem}_{chunk_id}_pre.parquet")
        scaffolds=pl.scan_csv(bed_file,has_header=False,separator="\t").select("column_1").unique().collect().to_series().to_list()
        if pathlib.Path(output_dir/f"{bam_file.stem}_{chunk_id}_pre.parquet").exists():
            mpileup=adjust_for_sequence_errors(
                mpile_frame=mpileup,
                null_model=null_model
            )
            mpileup = add_genome_info_to_mpileup(
                mpileup_df=mpileup.select(["chrom", "pos", "A", "C", "G", "T"]),
                scaffold_to_genome=stb.filter(pl.col("scaffold").is_in(scaffolds)).select(["scaffold", "genome"]),
            )
            mpileup = add_gene_info_to_mpileup(
                mpileup_df=mpileup,
                gene_range=pl.scan_csv(
                    gene_range_table,
                    has_header=False,
                    separator="\t",
                ).rename({
                    "column_1": "gene",
                    "column_2": "scaffold",
                    "column_3": "start",
                    "column_4": "end",
                }).filter(pl.col("scaffold").is_in(scaffolds)),
            )
            mpileup.select(["chrom", "genome", "gene", "pos", "A", "C", "G", "T"]).sink_parquet(
                output_dir / f"{bam_file.stem}_{chunk_id}.parquet",
                compression="zstd",
                engine="streaming",
            )
    cmd=["samtools", "view", "-F", "132", "-L", str(bed_file.absolute()), str(bam_file.absolute()), "|", "zipstrain", "utilities", "process-read-locs", "--output-file", f"{bam_file.stem}_read_locs_{chunk_id}.parquet"]
    proc = await asyncio.create_subprocess_shell(
                " ".join(cmd),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=output_dir
            )
    stdout, stderr = await proc.communicate()
    if proc.returncode != 0:
        raise Exception(f"Command failed with error: {stderr.decode().strip()}")

async def _run_profile_chunk_with_semaphore(
    semaphore: asyncio.Semaphore,
    *,
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    gene_range_table: pathlib.Path,
    stb: pl.LazyFrame,
    null_model: pl.LazyFrame,
    output_dir: pathlib.Path,
    chunk_id: int,
) -> None:
    """Bound chunk-level profiling concurrency."""
    async with semaphore:
        await _profile_chunk_task(
            bed_file=bed_file,
            bam_file=bam_file,
            gene_range_table=gene_range_table,
            stb=stb,
            null_model=null_model,
            output_dir=output_dir,
            chunk_id=chunk_id,
        )

async def profile_bam_in_chunks(
    bed_file:str,
    bam_file:str,
    gene_range_table:str,
    stb:pl.LazyFrame,
    null_model:pl.LazyFrame,
    output_dir:str,
    num_chunks:int=24,
    max_concurrency:int=4,
)->None:
    """
    Profile a BAM file in chunks using provided BED files.

    Parameters:
    bed_file (list[pathlib.Path]): A bed file describing all regions to be profiled.
    bam_file (pathlib.Path): Path to the BAM file.
    gene_range_table (pathlib.Path): Path to the gene range table.
    stb (pl.LazyFrame): The scaffold-to-genome mapping table.
    null_model (pl.LazyFrame): The null model to be used for adjusting for sequence errors.
    output_dir (pathlib.Path): Directory to save output files.
    num_chunks (int): Number of BED chunks to create.
    max_concurrency (int): Maximum number of chunks to process concurrently.
    """
    
    output_dir=pathlib.Path(output_dir)
    bam_file=pathlib.Path(bam_file)
    bed_file=pathlib.Path(bed_file)
    gene_range_table=pathlib.Path(gene_range_table)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"tmp").mkdir(exist_ok=True)
    bed_lf=pl.scan_csv(bed_file,has_header=False,separator="\t")
    gene_range_lf = pl.scan_csv(
        gene_range_table,
        has_header=False,
        separator="\t",
    ).rename({
        "column_1": "gene",
        "column_2": "scaffold",
        "column_3": "start",
        "column_4": "end",
    })
    bed_chunks=utils.split_lf_to_chunks(bed_lf, num_chunks)
    bed_chunk_files=[]
    for chunk_id, bed_file in enumerate(bed_chunks):
        bed_file.sink_csv(output_dir/"tmp"/f"bed_chunk_{chunk_id}.bed",include_header=False,separator="\t")
        bed_chunk_files.append(output_dir/"tmp"/f"bed_chunk_{chunk_id}.bed")
    tasks = []
    semaphore = asyncio.Semaphore(max(1, int(max_concurrency)))
    for chunk_id, bed_chunk_file in enumerate(bed_chunk_files):
        tasks.append(_run_profile_chunk_with_semaphore(
            semaphore,
            bed_file=bed_chunk_file,
            bam_file=bam_file,
            gene_range_table=gene_range_table,
            stb=stb,
            null_model=null_model,
            output_dir=output_dir/"tmp",
            chunk_id=chunk_id
        ))
    await asyncio.gather(*tasks) 
    pfs=[(output_dir/"tmp"/f"{bam_file.stem}_{chunk_id}.parquet", output_dir/"tmp"/f"{bam_file.stem}_read_locs_{chunk_id}.parquet" ) for chunk_id in range(len(bed_chunk_files)) if (output_dir/"tmp"/f"{bam_file.stem}_{chunk_id}.parquet").exists()]

    mpile_container: list[pl.LazyFrame] = []
    read_loc_pfs: list[pl.LazyFrame] = []
    
    for pf, read_loc_pf in pfs:
        
        if pf.exists():
            mpile_container.append(pl.scan_parquet(pf).lazy())
        
        if read_loc_pf.exists():
            read_loc_pfs.append(pl.scan_parquet(read_loc_pf).lazy())
    if mpile_container:
        mpileup_df = pl.concat(mpile_container)
        mpileup_df.sink_parquet(output_dir/f"{bam_file.stem}_profile.parquet", compression='zstd', engine='streaming')
        utils.get_gene_stats(
            profile=mpileup_df,
            gene_bed=gene_range_lf,
            stb=stb,
        ).sink_parquet(
            output_dir/f"{bam_file.stem}_gene_stats.parquet",
            compression='zstd',
            engine='streaming',
        )
    
    if read_loc_pfs:
        read_loc_df = pl.concat(read_loc_pfs).rename(
            {
                "chrom":"scaffold",
            "pos":"loc",
        }
    )

    if mpile_container and read_loc_pfs:
        utils.get_genome_stats(
            profile=mpileup_df,
            read_loc_table=read_loc_df,
            stb=stb,
            bed=bed_lf.rename({"column_1":"scaffold","column_2":"start","column_3":"end"}),
        ).sink_parquet(output_dir/f"{bam_file.stem}_genome_stats.parquet", compression='zstd', engine='streaming')
    
    
    os.system(f"rm -r {output_dir}/tmp")

def profile_bam(
    bed_file:str,
    bam_file:str,
    gene_range_table:str,
    stb:pl.LazyFrame,
    null_model:pl.LazyFrame,
    output_dir:str,
    num_chunks:int=24,
    max_concurrency:int=4,
)->None:
    """
    Profile a BAM file in chunks using provided BED files.

    Parameters:
    bed_file (list[pathlib.Path]): A bed file describing all regions to be profiled.
    bam_file (pathlib.Path): Path to the BAM file.
    gene_range_table (pathlib.Path): Path to the gene range table.
    stb (pl.LazyFrame): Scaffold-to-genome mapping table.
    null_model (pl.LazyFrame): The null model to be used for adjusting for sequence errors.
    output_dir (pathlib.Path): Directory to save output files.
    num_chunks (int): Number of BED chunks to create.
    max_concurrency (int): Maximum number of chunks to process concurrently.
    """
    asyncio.run(profile_bam_in_chunks(
        bed_file=bed_file,
        bam_file=bam_file,
        gene_range_table=gene_range_table,
        stb=stb,
        null_model=null_model,
        output_dir=output_dir,
        num_chunks=num_chunks,
        max_concurrency=max_concurrency,
    ))
