"""
zipstrain.utils
========================
This module provides utility functions for profiling and compare operations.
"""
import pathlib
import polars as pl
import sys
import re
import pyarrow as pa
import pyarrow.parquet as pq
from intervaltree import IntervalTree
from collections import defaultdict,Counter
from functools import reduce
from scipy.stats import poisson
import subprocess
import pdb

def build_null_poisson(error_rate:float=0.001,
                       max_total_reads:int=10000,
                       p_threshold:float=0.05)->list[float]:
    """
    Build a null model to correct for sequencing errors based on the Poisson distribution.

    Parameters:
    error_rate (float): Error rate for the sequencing technology.
    max_total_reads (int): Maximum total reads to consider.
    p_threshold (float): Significance threshold for the Poisson distribution.

    Returns:
    pl.DataFrame: DataFrame containing total reads and maximum error count thresholds.
    """ 
    records = []
    for n in range(1, max_total_reads + 1):
        lam = n * (error_rate / 3)
        k = 0
        while poisson.sf(k - 1, lam) > p_threshold:
            k += 1
        records.append((n, k - 1))
    return records



def clean_bases(bases: str, indel_re: re.Pattern) -> str:
    """
    Remove read start/end markers and indels from bases string using regex.
    Returns cleaned uppercase string of bases only.
    Args:
        bases (str): The bases string from mpileup.
        indel_re (re.Pattern): Compiled regex pattern to match indels and markers.
    
    """
    cleaned = []
    i = 0
    while i < len(bases):
        m = indel_re.match(bases, i)
        if m:
            if m.group(0).startswith('+') or m.group(0).startswith('-'):
                # indel length
                indel_len = int(m.group(1))
                i = m.end() + indel_len
            else:
                i = m.end()
        else:
            cleaned.append(bases[i].upper())
            i += 1
    return ''.join(cleaned)

def count_bases(bases: str):
    """
    Count occurrences of A, C, G, T in the cleaned bases string.
    Args:
        bases (str): Cleaned bases string.
    Returns:
        dict: Dictionary with counts of A, C, G, T.
    """
    counts = Counter(bases)
    return {
        'A': counts.get('A', 0),
        'C': counts.get('C', 0),
        'G': counts.get('G', 0),
        'T': counts.get('T', 0),
    }

def process_mpileup_function(gene_range_table_loc, batch_bed, batch_size, output_file):
    """
    Process mpileup files and save the results in a Parquet file.

    Parameters:
    gene_range_table_loc (str): Path to the gene range table in TSV format.
    batch_bed (str): Path to the batch BED file.
    batch_size (int): Buffer size for processing stdin from samtools.
    output_file (str): Path to save the output Parquet file.
    """
    indel_re = re.compile(r'\^.|[\$]|[+-](\d+)')
    gene_ranges_pl = pl.scan_csv(gene_range_table_loc,separator='\t', has_header=False).rename({
        "column_1": "scaffold",
        "column_2": "start",
        "column_3": "end",
        "column_4": "gene"
    })
    scaffolds = pl.read_csv(batch_bed, separator='\t', has_header=False)["column_1"].unique().to_list()
    gene_ranges_pl = gene_ranges_pl.filter(pl.col("scaffold").is_in(scaffolds)).collect()
    gene_ranges = defaultdict(IntervalTree)
    for row in gene_ranges_pl.iter_rows(named=True):
        gene_ranges[row["scaffold"]].addi(row["start"], row["end"] + 1, row["gene"])

    schema = pa.schema([
        ('chrom', pa.string()),
        ('pos', pa.int32()),
        ('gene', pa.string()),
        ('A', pa.uint16()),
        ('C', pa.uint16()),
        ('G', pa.uint16()),
        ('T', pa.uint16()),
    ])

    chroms = []
    positions = []
    genes = []
    As = []
    Cs = []
    Gs = []
    Ts = []

    writer = None
    def flush_batch():
        nonlocal writer
        if not chroms:
            return
        batch = pa.RecordBatch.from_arrays([
            pa.array(chroms, type=pa.string()),
            pa.array(positions, type=pa.int32()),
            pa.array(genes, type=pa.string()),
            pa.array(As, type=pa.uint16()),
            pa.array(Cs, type=pa.uint16()),
            pa.array(Gs, type=pa.uint16()),
            pa.array(Ts, type=pa.uint16()),
        ], schema=schema)

        if writer is None:
            # Open writer for the first time
            writer = pq.ParquetWriter(output_file, schema, compression='zstd')
        writer.write_table(pa.Table.from_batches([batch]))

        # Clear buffers
        chroms.clear()
        positions.clear()
        genes.clear()
        As.clear()
        Cs.clear()
        Gs.clear()
        Ts.clear()
    for line in sys.stdin:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        chrom, pos, _, _, bases = fields[:5]

        cleaned = clean_bases(bases, indel_re)
        counts = count_bases(cleaned)

        chroms.append(chrom)
        positions.append(int(pos))
        matches = gene_ranges[chrom][int(pos)]
        genes.append(next(iter(matches)).data if matches else "NA")
        As.append(counts['A'])
        Cs.append(counts['C'])
        Gs.append(counts['G'])
        Ts.append(counts['T'])

        if len(chroms) >= batch_size:
            flush_batch()

    # Flush remaining data
    flush_batch()

    if writer:
        writer.close()

def process_read_location(output_file:str, batch_size:int=10000)->None:
    """
    This function takes the output of samtools view -F 132 and processes it to extract read locations in a parquet file.
    """
    schema = pa.schema([
        ('chrom', pa.string()),
        ('pos', pa.int32()),
    ])
    writer = None
    chroms = []
    positions = []
    def flush_batch():
        nonlocal writer
        if not chroms:
            return
        batch = pa.RecordBatch.from_arrays([
            pa.array(chroms, type=pa.string()),
            pa.array(positions, type=pa.int32()),
        ], schema=schema)

        if writer is None:
            # Open writer for the first time
            writer = pq.ParquetWriter(output_file, schema, compression='zstd')
        writer.write_table(pa.Table.from_batches([batch]))

        # Clear buffers
        chroms.clear()
        positions.clear()
    for line in sys.stdin:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 4:
            continue
        chrom, pos = fields[2], fields[3]
        chroms.append(chrom)
        positions.append(int(pos))
        if len(chroms) >= batch_size:
            flush_batch()
    # Flush remaining data
    flush_batch()
    if writer:
        writer.close()


def extract_genome_length(stb: pl.LazyFrame, bed_table: pl.LazyFrame) -> pl.LazyFrame:
    """
    Extract the genome length information from the scaffold-to-genome mapping table.

    Parameters:
    stb (pl.LazyFrame): Scaffold-to-bin mapping table.
    bed_table (pl.LazyFrame): BED table containing genomic regions.

    Returns:
    pl.LazyFrame: A LazyFrame containing the genome lengths.
    """
    lf= bed_table.select(
        pl.col("scaffold"),
        (pl.col("end") - pl.col("start")).alias("scaffold_length")
    ).group_by("scaffold").agg(
        scaffold_length=pl.sum("scaffold_length")
    ).select(
        pl.col("scaffold").alias("scaffold"),
        pl.col("scaffold_length")
    ).join(
        stb.select(
            pl.col("scaffold").alias("scaffold"),
            pl.col("genome").alias("genome")
        ),
        on="scaffold",
        how="left"
    ).group_by("genome").agg(
        genome_length=pl.sum("scaffold_length")
    ).select(
        pl.col("genome"),
        pl.col("genome_length")
    )
    return lf

def make_the_bed(db_fasta_dir: str | pathlib.Path, max_scaffold_length: int = 500_000) -> pl.DataFrame:
    """
    Create a BED file from the database in fasta format.

    Parameters:
    db_fasta_dir (Union[str, pathlib.Path]): Path to the fasta file.
    max_scaffold_length (int): Splits scaffolds longer than this into multiple entries of length <= max_scaffold_length.

    Returns:
    pl.LazyFrame: A LazyFrame containing the BED data.
    """
    db_fasta_dir = pathlib.Path(db_fasta_dir)
    if not db_fasta_dir.is_file():
        raise FileNotFoundError(f"{db_fasta_dir} is not a valid fasta file.")

    records = []
    with db_fasta_dir.open() as f:
        scaffold = None
        seq_chunks = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Process the previous scaffold
                if scaffold is not None:
                    seq = ''.join(seq_chunks)
                    for start in range(0, len(seq), max_scaffold_length):
                        end = min(start + max_scaffold_length, len(seq))
                        records.append((scaffold, start, end))
                # Start new scaffold
                scaffold = line[1:].split()[0]  # ID only (up to first whitespace)
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # Don't forget the last scaffold
        if scaffold is not None:
            seq = ''.join(seq_chunks)
            for start in range(0, len(seq), max_scaffold_length):
                end = min(start + max_scaffold_length, len(seq))
                records.append((scaffold, start, end))

    return pl.DataFrame(records, schema=["scaffold", "start", "end"], orient="row")


def get_genome_breadth_matrix(
                              profile:pl.LazyFrame,
                              name:str,
                              genome_length: pl.LazyFrame,
                              stb: pl.LazyFrame,
                              min_cov: int = 1)-> pl.LazyFrame:
    """
    Get the genome breadth matrix from the provided profiles and scaffold-to-genome mapping.
    Parameters:
    profiles (list): List of tuples containing profile names and their corresponding LazyFrames.
    stb (pl.LazyFrame): Scaffold-to-genome mapping table.
    min_cov (int): Minimum coverage to consider a position. 
    Returns:
    pl.LazyFrame: A LazyFrame containing the genome breadth matrix.
    """
    profile = profile.filter((pl.col("A") + pl.col("C") + pl.col("G") + pl.col("T")) >= min_cov)
    profile=profile.group_by("chrom").agg(
        breadth=pl.count()
    ).select(
        pl.col("chrom").alias("scaffold"),
        pl.col("breadth")
    ).join(
        stb,
        on="scaffold",
        how="left"
    )
    profile=profile.join(genome_length, on="genome", how="left")
    
    profile=profile.group_by("genome").agg(
        genome_length=pl.first("genome_length"),
        breadth=pl.col("breadth").sum())
    profile = profile.with_columns(
        (pl.col("breadth")/ pl.col("genome_length")).alias("breadth")
    )
    return profile.select(
            pl.col("genome"),
            pl.col("breadth").alias(name)
        )
        
def collect_breadth_tables(
    breadth_tables: list[pl.LazyFrame],
) -> pl.LazyFrame:
    """
    Collect multiple genome breadth tables into a single LazyFrame.
    
    Parameters:
    breadth_tables (list[pl.LazyFrame]): List of LazyFrames containing genome breadth data.
    
    Returns:
    pl.LazyFrame: A LazyFrame containing the combined genome breadth data.
    """
    if not breadth_tables:
        raise ValueError("No breadth tables provided.")

    return reduce(lambda x, y: x.join(y, on="genome", how="outer", coalesce=True), breadth_tables)

def check_samtools():
    try:
        result = subprocess.run(
            ["samtools", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        return True
    except:
        print("Samtools is not installed or not found in PATH. Please install samtools to use all of the ZipStrain's functionalities.")
        return False

def split_lf_to_chunks(lf:pl.LazyFrame,num_chunks:int)->list[pl.LazyFrame]:
    """
    Split a Polars LazyFrame into smaller chunks.

    Parameters:
    lf (pl.LazyFrame): The input LazyFrame to be split.
    num_chunks (int): The number of chunks to split the LazyFrame into.

    Returns:
    list[pl.LazyFrame]: A list of smaller LazyFrames.
    """
    total_rows = lf.select(pl.count()).collect().item()
    chunk_size = total_rows // num_chunks
    chunks = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_chunks - 1 else total_rows
        chunk = lf.slice(start, end - start)
        chunks.append(chunk)
    return chunks


def get_genome_gaps(
    read_loc_table: pl.LazyFrame,
    stb: pl.LazyFrame,
    genome_length: pl.LazyFrame,
                    )-> pl.LazyFrame:
    read_loc_table=read_loc_table.sort(["scaffold",'loc'])
    read_loc_table=read_loc_table.with_columns(
        (pl.col("loc") - pl.col("loc").shift(1).over("scaffold")).alias("gap_length")
    ).join(
        stb,
        on="scaffold",
        how="left"
    )
    delta=read_loc_table.group_by("genome").agg(
        rn=pl.len()).join(
        genome_length,
        on="genome",
        how="left"
    ).with_columns(
        delta=(pl.col("genome_length")/pl.col("rn")).round().alias("delta")).select(
        pl.col("genome"),
        pl.col("delta"),
        pl.col("rn")
        )
    read_loc_table=read_loc_table.join(
        delta,
        on="genome",
        how="left"
    )
    read_loc_table=read_loc_table.filter(
        pl.col("gap_length") > pl.col("delta")
    ).group_by(["genome","gap_length"]).agg(
        pd=(pl.len()/(pl.col("rn").first()-1)),
        delta=pl.col("delta").first()
    ).with_columns(
        pd= pl.col("pd") * (pl.col("gap_length")-pl.col("delta"))
    ).group_by("genome").agg(
        fug=(pl.col("delta").first()-pl.col("pd").sum())/pl.col("delta").first()
    )
    return read_loc_table.select(
        pl.col("genome"),
        pl.col("fug")
    )

def get_genome_stats(
    profile:pl.LazyFrame,
    bed: pl.LazyFrame,
    stb: pl.LazyFrame,
    read_loc_table: pl.LazyFrame,
    ber:float=0.5,
    fug:float=2,
    min_cov_use_fug:int=0.1
)->pl.LazyFrame:

    genome_lengths=extract_genome_length(stb, bed)
    genome_gap_stats= get_genome_gaps(read_loc_table, stb, genome_lengths)
    profile=profile.join(
        stb,
        left_on="chrom",
        right_on="scaffold",
        how="left"
    ).select(
        pl.col("chrom"),
        pl.col("genome"),
        (pl.col("A")+pl.col("C")+pl.col("G")+pl.col("T")).alias("coverage")
    )
    profile=profile.group_by("genome").agg(
        total_covered_sites=pl.len(),
        coverage=pl.col("coverage").sum()
    ).join(
        genome_lengths,
        on="genome",
        how="left"
    ).join(
        genome_gap_stats,
        on="genome",
        how="left"
    ).with_columns(
        coverage=(pl.col("coverage")/pl.col("genome_length")),
        breadth=(pl.col("total_covered_sites")/pl.col("genome_length")),
    ).with_columns(
        ber=pl.col("breadth")/(1-(-0.883*pl.col("coverage")).exp()),
        fug=pl.col("fug")
    ).with_columns(
    
    pl.when(pl.col("coverage") > min_cov_use_fug)
    .then(
        pl.col("ber") > ber
        ).otherwise(
            (pl.col("fug")/0.632 < fug) &
            (pl.col("ber") > ber)
        ).fill_null(False).alias("is_present"))

    return profile.select(
        pl.col("genome"),
        pl.col("coverage"),
        pl.col("breadth"),
        pl.col("ber"),
        pl.col("fug"),
        pl.col("is_present")
    )


