"""
zipstrain.utils
========================
This module provides utility functions for profiling and compare operations.
"""
import pathlib
from tempfile import TemporaryDirectory
import polars as pl
import sys
import re
import pyarrow as pa
import pyarrow.parquet as pq
from collections import Counter
from functools import reduce
from scipy.stats import poisson
import subprocess
import duckdb


class CallPresence:
    """This class provides methods to use the information """
    def validate_input(self,lf:pl.LazyFrame)->pl.LazyFrame:
        required_columns = {"genome", "coverage", "breadth", "ber", "fug"}
        missing_columns = required_columns - set(lf.collect_schema().names)
        if missing_columns:
            raise ValueError(f"Input LazyFrame is missing required columns: {missing_columns}")
        return lf
    
    def metapresence(self,
                       lf:pl.LazyFrame,
                       ber:float=0.5,
                       fug:float=0.5,
                       min_cov_use_fug:int=0.1
                       )->pl.LazyFrame:
        """
        Call presence/absence of genomes based on breadth, coverage, ber, and fug.
        Parameters:
        lf (pl.LazyFrame): Input LazyFrame with genome statistics.
        ber (float): Breadth error rate threshold.
        fug (float): Fragmented unassembled genome threshold.
        min_cov_use_fug (int): Minimum coverage to use fug for presence call.
        Returns:
        pl.LazyFrame: LazyFrame with presence/absence calls.
        """
        lf=lf.with_columns(
            pl.when(pl.col("coverage") > min_cov_use_fug)
            .then(
                pl.col("ber") > ber
                ).otherwise(
                    (pl.col("fug")/0.632 < fug) &
                    (pl.col("ber") > ber)
                ).fill_null(False).alias("is_present"))
        return lf.select(
            pl.col("genome"),
            pl.col("coverage"),
            pl.col("breadth"),
            pl.col("ber"),
            pl.col("fug"),
            pl.col("is_present")
        )
        
    def breadth_only(
        self,
        lf:pl.LazyFrame,
        breadth:float=0.5
        )->pl.LazyFrame:
        """
        Call presence/absence of genomes based on breadth only.
        Parameters:
        lf (pl.LazyFrame): Input LazyFrame with genome statistics.
        breadth (float): Breadth threshold.
        Returns:
        pl.LazyFrame: LazyFrame with presence/absence calls.
        """

        lf=lf.with_columns(
            (pl.col("breadth") > breadth).fill_null(False).alias("is_present"))
        return lf.select(
            pl.col("genome"),
            pl.col("coverage"),
            pl.col("breadth"),
            pl.col("ber"),
            pl.col("fug"),
            pl.col("is_present")
        )
    
    def coverage_only(
        self,
        lf:pl.LazyFrame,
        coverage:float=0.1
        )->pl.LazyFrame:
        """
        Call presence/absence of genomes based on coverage only.
        Parameters:
        lf (pl.LazyFrame): Input LazyFrame with genome statistics.
        coverage (float): Coverage threshold.
        Returns:
        pl.LazyFrame: LazyFrame
        """
        lf=lf.with_columns(
            (pl.col("coverage") > coverage).fill_null(False).alias("is_present"))
        return lf.select(
            pl.col("genome"),
            pl.col("coverage"),
            pl.col("breadth"),
            pl.col("ber"),
            pl.col("fug"),
            pl.col("is_present")
        )           
    
    def __call__(self, method: str, lf:pl.LazyFrame, **kwargs) -> pl.LazyFrame:
        self.validate_input(lf)
        return self.__getattribute__(method)(lf, **kwargs)




class EstimateAbundance:
    """This class provides methods to estimate abundance of genomes based on coverage."""
    def validate_input(self,lf:pl.LazyFrame)->pl.LazyFrame:
        required_columns = {"genome", "coverage","is_present","Rn"}
        missing_columns = required_columns - set(lf.collect_schema().names)
        if missing_columns:
            raise ValueError(f"Input LazyFrame is missing required columns: {missing_columns}")
        return lf

    def coverage_ratio(
        self,
        lf:pl.LazyFrame
        )->pl.LazyFrame:
        """
        Estimate abundance based on coverage ratio.
        Parameters:
        lf (pl.LazyFrame): Input LazyFrame with genome statistics.
        Returns:
        pl.LazyFrame: LazyFrame with estimated abundance.
        """
        lf=lf.with_columns(
            abundance=pl.when(pl.col("is_present"))
            .then(
                pl.col("coverage") /pl.col("coverage").sum()
            ).otherwise(pl.lit(0.0))
        )
        return lf
    
    def reads_ratio(
        self,
        lf:pl.LazyFrame
        )->pl.LazyFrame:
        """
        Estimate abundance based on reads ratio.
        Parameters:
        lf (pl.LazyFrame): Input LazyFrame with genome statistics.
        Returns:
        pl.LazyFrame: LazyFrame with estimated abundance.
        """
        lf=lf.with_columns(
            abundance=pl.when(pl.col("is_present"))
            .then(
                pl.col("Rn") /pl.col("total_reads").sum()
            ).otherwise(pl.lit(0.0))
        )
        return lf
    
    



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

def process_mpileup_function(batch_size, output_file):
    """
    Process mpileup files and save the results in a Parquet file.

    Parameters:
    gene_range_table_loc (str): Path to the gene range table in TSV format.
    batch_bed (str): Path to the batch BED file.
    batch_size (int): Buffer size for processing stdin from samtools.
    output_file (str): Path to save the output Parquet file.
    """
    indel_re = re.compile(r'\^.|[\$]|[+-](\d+)')


    schema = pa.schema([
        ('chrom', pa.string()),
        ('pos', pa.int32()),
        ('A', pa.uint16()),
        ('C', pa.uint16()),
        ('G', pa.uint16()),
        ('T', pa.uint16()),
    ])

    chroms = []
    positions = []
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
        breadth=pl.len()
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
    total_rows = lf.select(pl.len()).collect().item()
    if total_rows == 0:
        return []
    num_chunks = max(1, min(int(num_chunks), total_rows))
    chunk_size = (total_rows + num_chunks - 1) // num_chunks
    chunks = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, total_rows)
        if start >= total_rows:
            break
        chunk = lf.slice(start, end - start)
        chunks.append(chunk)
    return chunks


def _duckdb_quote_sql_string(value: str) -> str:
    return value.replace("'", "''")


def _duckdb_connect_with_temp_dir(temp_dir: str) -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect()
    conn.execute(f"SET temp_directory = '{_duckdb_quote_sql_string(temp_dir)}'")
    conn.execute("SET preserve_insertion_order = false")
    return conn


def get_genome_gaps(
    read_loc_table: pl.LazyFrame,
    stb: pl.LazyFrame,
    genome_length: pl.LazyFrame,
                    )-> pl.LazyFrame:
    # The read-location table can be much larger than the profile itself. Stream it to a
    # temporary parquet and let DuckDB handle the global sort/window with disk spill.
    with TemporaryDirectory(prefix="zipstrain_genome_gaps_") as tmp_dir:
        read_loc_path = pathlib.Path(tmp_dir) / "read_locs.parquet"
        read_loc_table.sink_parquet(read_loc_path, compression="zstd", engine="streaming")

        conn = _duckdb_connect_with_temp_dir(tmp_dir)
        try:
            conn.register(
                "stb_src",
                stb.select("scaffold", "genome").collect().to_arrow(),
            )
            conn.register(
                "genome_length_src",
                genome_length.select("genome", "genome_length").collect().to_arrow(),
            )
            read_loc_sql = _duckdb_quote_sql_string(str(read_loc_path))
            table = conn.execute(
                f"""
                WITH read_loc AS (
                  SELECT
                    CAST(scaffold AS VARCHAR) AS scaffold,
                    CAST(loc AS BIGINT) AS loc
                  FROM read_parquet('{read_loc_sql}')
                ),
                stb AS (
                  SELECT
                    CAST(scaffold AS VARCHAR) AS scaffold,
                    CAST(genome AS VARCHAR) AS genome
                  FROM stb_src
                ),
                genome_length AS (
                  SELECT
                    CAST(genome AS VARCHAR) AS genome,
                    CAST(genome_length AS DOUBLE) AS genome_length
                  FROM genome_length_src
                ),
                mapped AS (
                  SELECT
                    rl.scaffold,
                    rl.loc,
                    rl.loc - LAG(rl.loc) OVER (PARTITION BY rl.scaffold ORDER BY rl.loc) AS gap_length,
                    stb.genome
                  FROM read_loc rl
                  LEFT JOIN stb ON rl.scaffold = stb.scaffold
                ),
                delta AS (
                  SELECT
                    m.genome,
                    COUNT(*)::BIGINT AS rn,
                    AVG(m.gap_length)::DOUBLE AS gap_mean,
                    STDDEV_SAMP(m.gap_length)::DOUBLE AS gap_std,
                    ROUND(gl.genome_length / COUNT(*))::DOUBLE AS delta
                  FROM mapped m
                  LEFT JOIN genome_length gl ON m.genome = gl.genome
                  GROUP BY m.genome, gl.genome_length
                ),
                gap_excess AS (
                  SELECT
                    m.genome,
                    m.gap_length,
                    d.delta,
                    d.rn,
                    d.gap_mean,
                    d.gap_std,
                    (COUNT(*)::DOUBLE / NULLIF(d.rn - 1, 0)) * (m.gap_length - d.delta) AS pd
                  FROM mapped m
                  INNER JOIN delta d ON m.genome = d.genome
                  WHERE m.gap_length > d.delta
                  GROUP BY m.genome, m.gap_length, d.delta, d.rn, d.gap_mean, d.gap_std
                )
                SELECT
                  CAST(genome AS VARCHAR) AS genome,
                  (MIN(delta) - SUM(pd)) / MIN(delta) AS fug,
                  MIN(rn)::BIGINT AS rn,
                  MIN(gap_mean)::DOUBLE AS gap_mean,
                  MIN(gap_std)::DOUBLE AS gap_std
                FROM gap_excess
                GROUP BY genome
                """
            ).fetch_arrow_table()
        finally:
            conn.close()

    return pl.from_arrow(table).lazy()

def get_genome_stats(
    profile:pl.LazyFrame,
    bed: pl.LazyFrame,
    stb: pl.LazyFrame,
    read_loc_table: pl.LazyFrame,
    comp_min_cov_breadth: int = 5,
    hetro_min_freq: float = 0.8,
    hetro_min_cov: int = 5
    
)->pl.LazyFrame:
    genome_lengths = extract_genome_length(stb, bed)
    genome_gap_stats = get_genome_gaps(read_loc_table, stb, genome_lengths)
    with TemporaryDirectory(prefix="zipstrain_genome_stats_") as tmp_dir:
        profile_path = pathlib.Path(tmp_dir) / "profile.parquet"
        profile.select(
            pl.col("genome").cast(pl.Utf8).alias("genome"),
            pl.col("A"),
            pl.col("C"),
            pl.col("G"),
            pl.col("T"),
        ).sink_parquet(profile_path, compression="zstd", engine="streaming")

        conn = _duckdb_connect_with_temp_dir(tmp_dir)
        try:
            conn.register(
                "genome_length_src",
                genome_lengths.select("genome", "genome_length").collect().to_arrow(),
            )
            conn.register(
                "genome_gap_src",
                genome_gap_stats.select("genome", "fug", "rn", "gap_mean", "gap_std").collect().to_arrow(),
            )
            profile_sql = _duckdb_quote_sql_string(str(profile_path))
            table = conn.execute(
                f"""
                WITH profile_rows AS (
                  SELECT
                    CAST(genome AS VARCHAR) AS genome,
                    CAST(A AS DOUBLE) AS A,
                    CAST(C AS DOUBLE) AS C,
                    CAST(G AS DOUBLE) AS G,
                    CAST(T AS DOUBLE) AS T,
                    CAST(A AS DOUBLE) + CAST(C AS DOUBLE) + CAST(G AS DOUBLE) + CAST(T AS DOUBLE) AS coverage,
                    GREATEST(
                      CAST(A AS DOUBLE),
                      CAST(C AS DOUBLE),
                      CAST(G AS DOUBLE),
                      CAST(T AS DOUBLE)
                    ) AS max_base_count
                  FROM read_parquet('{profile_sql}')
                ),
                profile_stats AS (
                  SELECT
                    genome,
                    COUNT(*)::BIGINT AS total_covered_sites,
                    SUM(coverage)::DOUBLE AS covered_bases,
                    SUM(CASE WHEN coverage >= {int(comp_min_cov_breadth)} THEN 1 ELSE 0 END)::BIGINT AS cov_sites,
                    SUM(
                      CASE
                        WHEN coverage > {int(hetro_min_cov)}
                         AND (max_base_count / NULLIF(coverage, 0)) <= {float(hetro_min_freq)}
                        THEN 1
                        ELSE 0
                      END
                    )::DOUBLE
                    / NULLIF(
                        SUM(CASE WHEN coverage > {int(hetro_min_cov)} THEN 1 ELSE 0 END),
                        0
                      ) AS heterogeneity
                  FROM profile_rows
                  GROUP BY genome
                )
                SELECT
                  CAST(p.genome AS VARCHAR) AS genome,
                  p.covered_bases / gl.genome_length AS coverage,
                  p.total_covered_sites::DOUBLE / gl.genome_length AS breadth,
                  CAST(gl.genome_length AS BIGINT) AS genome_length,
                  CAST(gg.gap_mean AS DOUBLE) AS gap_mean,
                  CAST(gg.gap_std AS DOUBLE) AS gap_std,
                  CAST(p.cov_sites AS BIGINT) AS "{int(comp_min_cov_breadth)}x_cov_sites",
                  CAST(p.heterogeneity AS DOUBLE) AS heterogeneity,
                  (p.total_covered_sites::DOUBLE / gl.genome_length)
                    / (1 - EXP(-0.883 * (p.covered_bases / gl.genome_length))) AS ber,
                  CAST(gg.fug AS DOUBLE) AS fug,
                  COALESCE(CAST(gg.rn AS BIGINT), 0) AS reads_mapped
                FROM profile_stats AS p
                LEFT JOIN genome_length_src AS gl
                  ON p.genome = CAST(gl.genome AS VARCHAR)
                LEFT JOIN genome_gap_src AS gg
                  ON p.genome = CAST(gg.genome AS VARCHAR)
                """
            ).fetch_arrow_table()
        finally:
            conn.close()

    return pl.from_arrow(table).lazy()


def get_gene_stats(
    profile:pl.LazyFrame,
    gene_bed: pl.LazyFrame,
    stb: pl.LazyFrame,
)->pl.LazyFrame:
    with TemporaryDirectory(prefix="zipstrain_gene_stats_") as tmp_dir:
        profile_path = pathlib.Path(tmp_dir) / "profile.parquet"
        profile.select(
            pl.col("genome").cast(pl.Utf8).alias("genome"),
            pl.col("gene").cast(pl.Utf8).alias("gene"),
            pl.col("A"),
            pl.col("C"),
            pl.col("G"),
            pl.col("T"),
        ).sink_parquet(profile_path, compression="zstd", engine="streaming")

        conn = _duckdb_connect_with_temp_dir(tmp_dir)
        try:
            conn.register(
                "gene_bed_src",
                gene_bed.select("gene", "scaffold", "start", "end").collect().to_arrow(),
            )
            conn.register(
                "stb_src",
                stb.select("scaffold", "genome").collect().to_arrow(),
            )
            profile_sql = _duckdb_quote_sql_string(str(profile_path))
            table = conn.execute(
                f"""
                WITH gene_lengths AS (
                  SELECT DISTINCT
                    CAST(s.genome AS VARCHAR) AS genome,
                    CAST(g.gene AS VARCHAR) AS gene,
                    CASE
                      WHEN CAST(g."end" AS BIGINT) >= CAST(g.start AS BIGINT)
                      THEN CAST(g."end" AS BIGINT) - CAST(g.start AS BIGINT) + 1
                      ELSE 0
                    END AS length
                  FROM gene_bed_src AS g
                  LEFT JOIN stb_src AS s
                    ON CAST(g.scaffold AS VARCHAR) = CAST(s.scaffold AS VARCHAR)
                ),
                covered_stats AS (
                  SELECT
                    CAST(genome AS VARCHAR) AS genome,
                    CAST(gene AS VARCHAR) AS gene,
                    COUNT(*)::BIGINT AS total_covered_sites,
                    SUM(
                      CAST(A AS DOUBLE) +
                      CAST(C AS DOUBLE) +
                      CAST(G AS DOUBLE) +
                      CAST(T AS DOUBLE)
                    )::DOUBLE AS covered_bases
                  FROM read_parquet('{profile_sql}')
                  WHERE gene IS NOT NULL
                    AND CAST(gene AS VARCHAR) <> 'NA'
                  GROUP BY 1, 2
                )
                SELECT
                  CAST(gl.genome AS VARCHAR) AS genome,
                  CAST(gl.gene AS VARCHAR) AS gene,
                  CAST(gl.length AS BIGINT) AS length,
                  CASE
                    WHEN gl.length > 0
                    THEN COALESCE(cs.total_covered_sites, 0)::DOUBLE / gl.length
                    ELSE 0.0
                  END AS breadth,
                  CASE
                    WHEN gl.length > 0
                    THEN COALESCE(cs.covered_bases, 0)::DOUBLE / gl.length
                    ELSE 0.0
                  END AS coverage
                FROM gene_lengths AS gl
                LEFT JOIN covered_stats AS cs
                  ON gl.genome = cs.genome
                 AND gl.gene = cs.gene
                """
            ).fetch_arrow_table()
        finally:
            conn.close()

    return pl.from_arrow(table).lazy()
