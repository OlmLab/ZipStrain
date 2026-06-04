"""zipstrain.profile
========================
This module provides functions and utilities to profile a bamfile.
By profile we mean generating gene, genome, and nucleotide counts at each position on the reference.
This is a fundamental step for downstream analysis in zipstrain.
"""
import pathlib
import polars as pl
from bisect import bisect_right
from typing import Generator, Optional
from zipstrain import utils
import asyncio
import os
import re
import shutil
import duckdb
import tempfile
import pyarrow as pa
import pyarrow.parquet as pq


PROFILE_SORTED_METADATA_KEY = "zipstrain_sorted_by"
PROFILE_SORTED_METADATA_VALUE = "chrom,pos"
PROFILE_WRITE_BATCH_SIZE = 10_000
PROFILE_PARQUET_SCHEMA = pa.schema(
    [
        ("chrom", pa.string()),
        ("genome", pa.string()),
        ("gene", pa.string()),
        ("pos", pa.int32()),
        ("A", pa.int32()),
        ("C", pa.int32()),
        ("G", pa.int32()),
        ("T", pa.int32()),
    ]
)
MPILEUP_INDEL_RE = re.compile(r"\^.|[\$]|[+-](\d+)")

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


def adjust_profile_parquet_for_sequence_errors(
    profile_parquet: pathlib.Path,
    null_model_parquet: pathlib.Path,
    output_file: pathlib.Path,
) -> None:
    """
    Apply sequence-error adjustment to an existing profile parquet.

    The output preserves the input column order/schema, except that temporary
    columns introduced during adjustment are not written.
    """
    input_path = pathlib.Path(profile_parquet).expanduser().resolve()
    output_path = pathlib.Path(output_file).expanduser().resolve()
    if input_path == output_path:
        raise ValueError("Input profile parquet and output file must be different paths.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile_scan = pl.scan_parquet(input_path)
    input_columns = profile_scan.collect_schema().names()
    adjusted = adjust_for_sequence_errors(
        mpile_frame=profile_scan,
        null_model=pl.scan_parquet(null_model_parquet),
    )
    adjusted.select(input_columns).sink_parquet(
        output_path,
        compression="zstd",
        engine="streaming",
    )
    
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


def empty_gene_range_table() -> pl.LazyFrame:
    """Return an empty gene-range LazyFrame with the expected schema."""
    return pl.DataFrame(
        schema={
            "gene": pl.Utf8,
            "scaffold": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
        }
    ).lazy()


def empty_gene_stats_table() -> pl.LazyFrame:
    """Return an empty gene-stats LazyFrame with the expected schema."""
    return pl.DataFrame(
        schema={
            "genome": pl.Utf8,
            "gene": pl.Utf8,
            "length": pl.Int64,
            "breadth": pl.Float64,
            "coverage": pl.Float64,
        }
    ).lazy()


def normalize_gene_range_table_path(
    gene_range_table: Optional[str | pathlib.Path],
) -> Optional[pathlib.Path]:
    """Treat missing or empty gene-range files as absent annotations."""
    if gene_range_table in (None, ""):
        return None
    path = pathlib.Path(gene_range_table)
    if path.exists() and path.stat().st_size == 0:
        return None
    return path

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


def _write_sorted_profile_with_metadata(
    profile_lf: pl.LazyFrame,
    output_file: pathlib.Path,
    tmp_dir: pathlib.Path,
) -> None:
    """Write the final profile parquet and only fall back to sorting when needed."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = tmp_dir / f"{output_file.stem}.candidate.parquet"
    candidate_path.unlink(missing_ok=True)

    profile_lf.sink_parquet(
        candidate_path,
        compression="zstd",
        engine="streaming",
    )

    if _parquet_rows_are_coordinate_sorted(candidate_path):
        pl.scan_parquet(candidate_path).sink_parquet(
            output_file,
            compression="zstd",
            engine="streaming",
            metadata={PROFILE_SORTED_METADATA_KEY: PROFILE_SORTED_METADATA_VALUE},
        )
    else:
        _sort_existing_profile_parquet(
            input_file=candidate_path,
            output_file=output_file,
            tmp_dir=tmp_dir,
        )


def _sort_existing_profile_parquet(
    input_file: pathlib.Path,
    output_file: pathlib.Path,
    tmp_dir: pathlib.Path,
) -> None:
    """Sort an existing classic profile parquet and attach sortedness metadata."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    sorted_path = tmp_dir / f"{output_file.stem}.sorted.parquet"
    sorted_path.unlink(missing_ok=True)

    conn = duckdb.connect()
    try:
        in_sql = _duckdb_quote_sql_string(str(input_file))
        out_sql = _duckdb_quote_sql_string(str(sorted_path))
        conn.execute(
            f"""
            COPY (
              SELECT
                chrom,
                genome,
                gene,
                pos,
                A,
                C,
                G,
                T
              FROM read_parquet('{in_sql}')
              ORDER BY chrom, pos
            ) TO '{out_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        conn.close()

    pl.scan_parquet(sorted_path).sink_parquet(
        output_file,
        compression="zstd",
        engine="streaming",
        metadata={PROFILE_SORTED_METADATA_KEY: PROFILE_SORTED_METADATA_VALUE},
    )


def sort_profile_parquet_in_place(
    profile_parquet: pathlib.Path,
    tmp_dir: Optional[pathlib.Path] = None,
) -> None:
    """Sort a classic profile parquet in place and attach sortedness metadata."""
    profile_path = pathlib.Path(profile_parquet).expanduser().resolve()
    if not profile_path.exists():
        raise ValueError(f"Profile parquet does not exist: {profile_path}")

    tmp_root = (
        pathlib.Path(tmp_dir).expanduser().resolve()
        if tmp_dir is not None
        else profile_path.parent
    )
    tmp_root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"{profile_path.stem}.sort.", dir=tmp_root) as work_dir_str:
        work_dir = pathlib.Path(work_dir_str)
        rewritten_path = work_dir / profile_path.name
        _sort_existing_profile_parquet(
            input_file=profile_path,
            output_file=rewritten_path,
            tmp_dir=work_dir,
        )
        metadata = pl.read_parquet_metadata(rewritten_path)
        if metadata.get(PROFILE_SORTED_METADATA_KEY) != PROFILE_SORTED_METADATA_VALUE:
            raise RuntimeError("Failed to attach sortedness metadata to rewritten profile parquet.")
        os.replace(rewritten_path, profile_path)


def _profile_parquet_is_coordinate_sorted(profile_parquet: pathlib.Path) -> bool:
    try:
        metadata = pl.read_parquet_metadata(profile_parquet)
    except Exception:
        return False
    return metadata.get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE


def _parquet_rows_are_coordinate_sorted(profile_parquet: pathlib.Path) -> bool:
    """Check row order without materializing the full parquet in memory."""
    check_expr = (
        pl.when(pl.col("chrom").shift(1).is_null())
        .then(True)
        .otherwise(
            (pl.col("chrom") > pl.col("chrom").shift(1))
            | (
                (pl.col("chrom") == pl.col("chrom").shift(1))
                & (pl.col("pos") >= pl.col("pos").shift(1))
            )
        )
        .all()
        .alias("is_sorted")
    )
    return bool(
        pl.scan_parquet(profile_parquet)
        .select(check_expr)
        .collect(engine="streaming")
        .item()
    )


def _build_null_model_threshold_lookup(null_model: pl.LazyFrame) -> dict[int, int]:
    return {
        int(cov): int(max_error_count)
        for cov, max_error_count in (
            null_model.select(["cov", "max_error_count"]).collect().iter_rows()
        )
    }


def _build_scaffold_to_genome_lookup(stb: pl.LazyFrame) -> dict[str, str]:
    return {
        str(scaffold): str(genome)
        for scaffold, genome in (
            stb.select(["scaffold", "genome"]).collect().iter_rows()
        )
    }


def _build_gene_lookup(
    gene_range: pl.LazyFrame,
) -> dict[str, tuple[list[int], list[tuple[int, int, str]]]]:
    gene_df = (
        gene_range.select(["gene", "scaffold", "start", "end"])
        .collect()
        .sort(["scaffold", "start"])
    )
    gene_lookup: dict[str, tuple[list[int], list[tuple[int, int, str]]]] = {}
    for gene, scaffold, start, end in gene_df.iter_rows():
        chrom = str(scaffold)
        if chrom not in gene_lookup:
            gene_lookup[chrom] = ([], [])
        starts, intervals = gene_lookup[chrom]
        starts.append(int(start))
        intervals.append((int(start), int(end), str(gene)))
    return gene_lookup


def _lookup_gene_name(
    gene_lookup: dict[str, tuple[list[int], list[tuple[int, int, str]]]],
    chrom: str,
    pos: int,
) -> str:
    starts_and_intervals = gene_lookup.get(chrom)
    if starts_and_intervals is None:
        return "NA"
    starts, intervals = starts_and_intervals
    interval_index = bisect_right(starts, pos) - 1
    if interval_index < 0:
        return "NA"
    _, end, gene = intervals[interval_index]
    return gene if pos <= end else "NA"


async def _stream_profile_mpileup_chunk_to_parquet(
    *,
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    output_parquet: pathlib.Path,
    null_model_thresholds: dict[int, int],
    scaffold_to_genome_lookup: dict[str, str],
    gene_lookup: dict[str, tuple[list[int], list[tuple[int, int, str]]]],
    batch_size: int = PROFILE_WRITE_BATCH_SIZE,
) -> None:
    """Stream one mpileup chunk directly into a final annotated parquet file."""
    cmd = ["samtools", "mpileup", "-A", "-l", str(bed_file.absolute()), str(bam_file.absolute())]
    proc = await asyncio.create_subprocess_shell(
        " ".join(cmd),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=output_parquet.parent,
    )

    chroms: list[str] = []
    genomes: list[str] = []
    genes: list[str] = []
    positions: list[int] = []
    As: list[int] = []
    Cs: list[int] = []
    Gs: list[int] = []
    Ts: list[int] = []
    writer: Optional[pq.ParquetWriter] = None
    success = False

    def flush_batch() -> None:
        nonlocal writer
        if not chroms:
            return
        batch = pa.RecordBatch.from_arrays(
            [
                pa.array(chroms, type=pa.string()),
                pa.array(genomes, type=pa.string()),
                pa.array(genes, type=pa.string()),
                pa.array(positions, type=pa.int32()),
                pa.array(As, type=pa.int32()),
                pa.array(Cs, type=pa.int32()),
                pa.array(Gs, type=pa.int32()),
                pa.array(Ts, type=pa.int32()),
            ],
            schema=PROFILE_PARQUET_SCHEMA,
        )
        if writer is None:
            writer = pq.ParquetWriter(str(output_parquet), PROFILE_PARQUET_SCHEMA, compression="zstd")
        writer.write_table(pa.Table.from_batches([batch]))
        chroms.clear()
        genomes.clear()
        genes.clear()
        positions.clear()
        As.clear()
        Cs.clear()
        Gs.clear()
        Ts.clear()

    try:
        while True:
            raw_line = await proc.stdout.readline()
            if not raw_line:
                break
            line = raw_line.decode().strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom, pos_str, _, _, bases = fields[:5]
            counts = utils.count_bases(utils.clean_bases(bases, MPILEUP_INDEL_RE))
            coverage = counts["A"] + counts["C"] + counts["G"] + counts["T"]
            threshold = null_model_thresholds.get(coverage)
            if threshold is None:
                adjusted_a = adjusted_c = adjusted_g = adjusted_t = 0
            else:
                adjusted_a = counts["A"] if counts["A"] >= threshold else 0
                adjusted_c = counts["C"] if counts["C"] >= threshold else 0
                adjusted_g = counts["G"] if counts["G"] >= threshold else 0
                adjusted_t = counts["T"] if counts["T"] >= threshold else 0
            pos = int(pos_str)
            chroms.append(chrom)
            genomes.append(scaffold_to_genome_lookup.get(chrom, "NA"))
            genes.append(_lookup_gene_name(gene_lookup, chrom, pos))
            positions.append(pos)
            As.append(adjusted_a)
            Cs.append(adjusted_c)
            Gs.append(adjusted_g)
            Ts.append(adjusted_t)
            if len(chroms) >= batch_size:
                flush_batch()
        stderr = await proc.stderr.read()
        if await proc.wait() != 0:
            raise Exception(f"Command failed with error: {stderr.decode().strip()}")
        flush_batch()
        success = True
    finally:
        if writer is not None:
            writer.close()
        if not success:
            output_parquet.unlink(missing_ok=True)


def _annotate_mpileup_chunk_with_duckdb(
    adjusted_mpileup_parquet: pathlib.Path,
    output_parquet: pathlib.Path,
    scaffold_to_genome: pl.LazyFrame,
    gene_range: pl.LazyFrame,
    mpileup_sorted: Optional[bool] = None,
    gene_range_sorted: bool = False,
) -> None:
    """Annotate one adjusted mpileup chunk with genome+gene in DuckDB."""
    if mpileup_sorted is None:
        mpileup_sorted = _profile_parquet_is_coordinate_sorted(adjusted_mpileup_parquet)
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
        mpileup_order_sql = "" if mpileup_sorted else "ORDER BY chrom, pos"
        gene_range_order_sql = "" if gene_range_sorted else "ORDER BY scaffold, start"
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
                {mpileup_order_sql}
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
                {gene_range_order_sql}
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
              ORDER BY chrom, pos
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
    scaffold_to_genome_lookup: dict[str, str],
    gene_lookup: dict[str, tuple[list[int], list[tuple[int, int, str]]]],
    null_model_thresholds: dict[int, int],
    output_dir:pathlib.Path,
    chunk_id:int
)->None:
    await _stream_profile_mpileup_chunk_to_parquet(
        bed_file=bed_file,
        bam_file=bam_file,
        output_parquet=output_dir / f"{bam_file.stem}_{chunk_id}.parquet",
        null_model_thresholds=null_model_thresholds,
        scaffold_to_genome_lookup=scaffold_to_genome_lookup,
        gene_lookup=gene_lookup,
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
    scaffold_to_genome_lookup: dict[str, str],
    gene_lookup: dict[str, tuple[list[int], list[tuple[int, int, str]]]],
    null_model_thresholds: dict[int, int],
    output_dir: pathlib.Path,
    chunk_id: int,
) -> None:
    """Bound chunk-level profiling concurrency."""
    async with semaphore:
        await _profile_chunk_task(
            bed_file=bed_file,
            bam_file=bam_file,
            scaffold_to_genome_lookup=scaffold_to_genome_lookup,
            gene_lookup=gene_lookup,
            null_model_thresholds=null_model_thresholds,
            output_dir=output_dir,
            chunk_id=chunk_id,
        )

async def profile_bam_in_chunks(
    bed_file:str,
    bam_file:str,
    gene_range_table: Optional[str],
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
    gene_range_table (pathlib.Path | None): Optional path to the gene range table.
    stb (pl.LazyFrame): The scaffold-to-genome mapping table.
    null_model (pl.LazyFrame): The null model to be used for adjusting for sequence errors.
    output_dir (pathlib.Path): Directory to save output files.
    num_chunks (int): Number of BED chunks to create.
    max_concurrency (int): Maximum number of chunks to process concurrently.
    """
    
    output_dir=pathlib.Path(output_dir)
    bam_file=pathlib.Path(bam_file)
    bed_file=pathlib.Path(bed_file)
    gene_range_table_path = normalize_gene_range_table_path(gene_range_table)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"tmp").mkdir(exist_ok=True)
    bed_lf=pl.scan_csv(bed_file,has_header=False,separator="\t")
    if gene_range_table_path is None:
        gene_range_lf = empty_gene_range_table()
    else:
        gene_range_lf = pl.scan_csv(
            gene_range_table_path,
            has_header=False,
            separator="\t",
        ).rename({
            "column_1": "gene",
            "column_2": "scaffold",
            "column_3": "start",
            "column_4": "end",
        })
    null_model_thresholds = _build_null_model_threshold_lookup(null_model)
    scaffold_to_genome_lookup = _build_scaffold_to_genome_lookup(stb)
    gene_lookup = _build_gene_lookup(gene_range_lf)
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
            scaffold_to_genome_lookup=scaffold_to_genome_lookup,
            gene_lookup=gene_lookup,
            null_model_thresholds=null_model_thresholds,
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
        _write_sorted_profile_with_metadata(
            profile_lf=mpileup_df,
            output_file=output_dir / f"{bam_file.stem}_profile.parquet",
            tmp_dir=output_dir / "tmp",
        )
        if gene_range_table_path is None:
            empty_gene_stats_table().sink_parquet(
                output_dir/f"{bam_file.stem}_gene_stats.parquet",
                compression='zstd',
                engine='streaming',
            )
        else:
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
    
    
    shutil.rmtree(output_dir / "tmp", ignore_errors=True)

def profile_bam(
    bed_file:str,
    bam_file:str,
    gene_range_table: Optional[str],
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
    gene_range_table (pathlib.Path | None): Optional path to the gene range table.
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
