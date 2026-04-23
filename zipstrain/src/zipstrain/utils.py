"""
zipstrain.utils
========================
This module provides utility functions for profiling and compare operations.
"""
import os
import pathlib
import time
from tempfile import TemporaryDirectory
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import polars as pl
import sys
import re
import pyarrow as pa
import pyarrow.parquet as pq
from collections import Counter
from scipy.stats import poisson
import subprocess
import duckdb

CLASSIC_PROFILE_REQUIRED_COLUMNS = {"chrom", "pos", "gene", "genome", "A", "T", "C", "G"}
PROFILE_COVERAGE_STATS_REQUIRED_COLUMNS = {"chrom", "pos", "gene", "genome", "A", "T", "C", "G"}
DEFAULT_COVERAGE_STATS_SITE_THRESHOLD = 5
GENOME_PAIR_TABLE_SCHEMA = pa.schema(
    [
        pa.field("sample_name_1", pa.string()),
        pa.field("sample_name_2", pa.string()),
        pa.field("profile_location_1", pa.string()),
        pa.field("profile_location_2", pa.string()),
    ]
)

GENOME_COMPARE_OUTPUT_DTYPES = {
    "genome": pl.Utf8,
    "total_positions": pl.Int64,
    "share_allele_pos": pl.Int64,
    "genome_pop_ani": pl.Float64,
    "max_consecutive_length": pl.Int64,
    "shared_genes_count": pl.Int64,
    "identical_gene_count": pl.Int64,
    "perc_id_genes": pl.Float64,
    "sample_1": pl.Utf8,
    "sample_2": pl.Utf8,
}


def _merge_parquet_chunk(parquet_files: list[pathlib.Path], output_file: pathlib.Path) -> None:
    """Merge a bounded list of parquet files into one parquet file."""
    pl.concat([pl.scan_parquet(path) for path in parquet_files]).sink_parquet(
        output_file,
        compression="zstd",
    )


def _iter_parquet_batches(parquet_files: list[pathlib.Path], batch_size: int):
    """Yield parquet paths in fixed-size batches."""
    for start in range(0, len(parquet_files), batch_size):
        yield parquet_files[start:start + batch_size]


def merge_parquet_files(
    input_dir: str | pathlib.Path,
    output_file: str | pathlib.Path,
    batch_size: int = -1,
    progress_callback: Callable[[str], None] | None = None,
) -> pathlib.Path:
    """
    Merge parquet files from a directory into a single parquet file.

    When ``batch_size`` is ``-1`` the merge is performed in a single pass, which
    preserves the current behavior. For positive ``batch_size`` values, inputs
    are merged in staged batches to avoid building an excessively large lazy
    concatenation plan at once.
    """
    if batch_size == 0 or batch_size < -1:
        raise ValueError("batch_size must be -1 or a positive integer")

    input_path = pathlib.Path(input_dir)
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_resolved = output_path.resolve()

    parquet_files = sorted(
        path for path in input_path.glob("*.parquet")
        if path.resolve() != output_resolved
    )
    if not parquet_files:
        raise ValueError(f"No Parquet files found in directory: {input_dir}")

    def emit(message: str) -> None:
        if progress_callback is not None:
            progress_callback(message)

    total_files = len(parquet_files)
    emit(f"merge_parquet: found {total_files} parquet files")

    if batch_size == -1 or total_files <= batch_size:
        emit(f"merge_parquet: single-pass merge start ({total_files} files)")
        _merge_parquet_chunk(parquet_files, output_path)
        emit(f"merge_parquet: single-pass merge done -> {output_path}")
        return output_path

    with TemporaryDirectory(prefix="zipstrain_merge_parquet_", dir=output_path.parent) as tmp_dir:
        tmp_root = pathlib.Path(tmp_dir)
        batch_outputs: list[pathlib.Path] = []
        batch_count = (total_files + batch_size - 1) // batch_size
        processed_files = 0

        for batch_idx, batch_files in enumerate(_iter_parquet_batches(parquet_files, batch_size), start=1):
            batch_output = tmp_root / f"batch_{batch_idx:05d}.parquet"
            processed_files += len(batch_files)
            emit(
                f"merge_parquet: batch {batch_idx}/{batch_count} start "
                f"({len(batch_files)} files; processed {processed_files}/{total_files})"
            )
            _merge_parquet_chunk(batch_files, batch_output)
            batch_outputs.append(batch_output)
            emit(
                f"merge_parquet: batch {batch_idx}/{batch_count} done "
                f"-> {batch_output.name}"
            )

        emit(
            f"merge_parquet: final merge start "
            f"({len(batch_outputs)} batch outputs)"
        )
        _merge_parquet_chunk(batch_outputs, output_path)
        emit(f"merge_parquet: final merge done -> {output_path}")
        return output_path


def _read_gene_coverage_bed(gene_bed_file: str | pathlib.Path) -> pl.DataFrame:
    gene_bed = pl.read_csv(
        gene_bed_file,
        separator="\t",
        has_header=False,
        infer_schema_length=10000,
    )
    if gene_bed.width == 4:
        gene_bed = gene_bed.rename(
            {
                "column_1": "gene",
                "column_2": "scaffold",
                "column_3": "start",
                "column_4": "end",
            }
        ).with_columns(pl.lit(None, dtype=pl.Utf8).alias("genome"))
    elif gene_bed.width == 5:
        gene_bed = gene_bed.rename(
            {
                "column_1": "gene",
                "column_2": "scaffold",
                "column_3": "start",
                "column_4": "end",
                "column_5": "genome",
            }
        )
    else:
        raise ValueError("Gene BED must have 4 columns (gene, scaffold, start, end) or 5 columns (+ genome).")
    return gene_bed.select(
        pl.col("gene").cast(pl.Utf8),
        pl.col("scaffold").cast(pl.Utf8),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("genome").cast(pl.Utf8),
    )


def _read_genome_coverage_bed(genome_bed_file: str | pathlib.Path) -> pl.DataFrame:
    genome_bed = pl.read_csv(
        genome_bed_file,
        separator="\t",
        has_header=False,
        infer_schema_length=10000,
    )
    if genome_bed.width == 3:
        genome_bed = genome_bed.rename(
            {
                "column_1": "scaffold",
                "column_2": "start",
                "column_3": "end",
            }
        ).with_columns(pl.lit(None, dtype=pl.Utf8).alias("genome"))
    elif genome_bed.width == 4:
        genome_bed = genome_bed.rename(
            {
                "column_1": "scaffold",
                "column_2": "start",
                "column_3": "end",
                "column_4": "genome",
            }
        )
    else:
        raise ValueError("Genome BED must have 3 columns (scaffold, start, end) or 4 columns (+ genome).")
    return genome_bed.select(
        pl.col("scaffold").cast(pl.Utf8),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("genome").cast(pl.Utf8),
    )


def _validate_unique_scaffold_genome_mapping(scaffold_to_genome: pl.DataFrame) -> pl.DataFrame:
    conflicts = (
        scaffold_to_genome.drop_nulls(["scaffold", "genome"])
        .group_by("scaffold")
        .agg(pl.col("genome").n_unique().alias("genome_count"))
        .filter(pl.col("genome_count") > 1)
    )
    if conflicts.height > 0:
        raise ValueError("Scaffold-to-genome mapping must be unique per scaffold.")
    return scaffold_to_genome.drop_nulls(["scaffold", "genome"]).unique(["scaffold", "genome"])


def _resolve_scaffold_to_genome_mapping(
    profile_parquet: str | pathlib.Path,
) -> pl.DataFrame:
    mapping = pl.scan_parquet(profile_parquet).select(
        pl.col("chrom").cast(pl.Utf8).alias("scaffold"),
        pl.col("genome").cast(pl.Utf8).alias("genome"),
    ).unique().collect()
    return _validate_unique_scaffold_genome_mapping(mapping)


def _profile_rows_cte_sql(profile_sql: str) -> str:
    return f"""
        profile_rows AS (
          SELECT
            CAST(p.chrom AS VARCHAR) AS chrom,
            CAST(p.pos AS BIGINT) AS pos,
            CAST(p.genome AS VARCHAR) AS genome,
            CAST(p.gene AS VARCHAR) AS gene,
            CAST(p.A AS DOUBLE) AS A,
            CAST(p.C AS DOUBLE) AS C,
            CAST(p.G AS DOUBLE) AS G,
            CAST(p.T AS DOUBLE) AS T,
            CAST(p.A AS DOUBLE) + CAST(p.C AS DOUBLE) + CAST(p.G AS DOUBLE) + CAST(p.T AS DOUBLE) AS site_coverage
          FROM read_parquet('{profile_sql}') AS p
        )
    """


def get_coverage_stats(
    profile_parquet: str | pathlib.Path,
    gene_bed_file: str | pathlib.Path,
    genome_bed_file: str | pathlib.Path,
    output_dir: str | pathlib.Path,
    prefix: str,
    cov_site_threshold: int = DEFAULT_COVERAGE_STATS_SITE_THRESHOLD,
) -> dict[str, object]:
    profile_path = pathlib.Path(profile_parquet).expanduser().resolve()
    output_path = pathlib.Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    gene_output = output_path / f"{prefix}_gene_stats.parquet"
    genome_output = output_path / f"{prefix}_genome_stats.parquet"
    if gene_output.exists():
        raise FileExistsError(f"Output file already exists: {gene_output}")
    if genome_output.exists():
        raise FileExistsError(f"Output file already exists: {genome_output}")
    if cov_site_threshold < 1:
        raise ValueError("cov_site_threshold must be >= 1")

    profile_schema = set(pq.read_schema(profile_path).names)
    missing_columns = PROFILE_COVERAGE_STATS_REQUIRED_COLUMNS - profile_schema
    if missing_columns:
        raise ValueError(f"Profile parquet is missing required columns: {sorted(missing_columns)}")

    gene_bed = _read_gene_coverage_bed(gene_bed_file)
    genome_bed = _read_genome_coverage_bed(genome_bed_file)
    scaffold_to_genome = _resolve_scaffold_to_genome_mapping(
        profile_parquet=profile_path,
    )
    profile_sql = _duckdb_quote_sql_string(str(profile_path))
    cov_sites_column = f"{int(cov_site_threshold)}x_cov_sites"
    profile_rows_cte = _profile_rows_cte_sql(profile_sql=profile_sql)

    with TemporaryDirectory(prefix="zipstrain_coverage_stats_") as tmp_dir:
        conn = _duckdb_connect_with_temp_dir(tmp_dir)
        try:
            conn.register("gene_bed_src", gene_bed.to_arrow())
            conn.register("genome_bed_src", genome_bed.to_arrow())
            conn.register("scaffold_map_src", scaffold_to_genome.to_arrow())

            gene_table = conn.execute(
                f"""
                WITH
                {profile_rows_cte},
                gene_lengths AS (
                  SELECT
                    sm.genome AS genome,
                    CAST(g.gene AS VARCHAR) AS gene,
                    SUM(
                      CASE
                        WHEN CAST(g."end" AS BIGINT) >= CAST(g.start AS BIGINT)
                        THEN CAST(g."end" AS BIGINT) - CAST(g.start AS BIGINT) + 1
                        ELSE 0
                      END
                    )::BIGINT AS length
                  FROM gene_bed_src AS g
                  INNER JOIN scaffold_map_src AS sm
                    ON CAST(g.scaffold AS VARCHAR) = CAST(sm.scaffold AS VARCHAR)
                  GROUP BY 1, 2
                ),
                covered_stats AS (
                  SELECT
                    CAST(genome AS VARCHAR) AS genome,
                    CAST(gene AS VARCHAR) AS gene,
                    COUNT(*)::BIGINT AS total_covered_sites,
                    SUM(site_coverage)::DOUBLE AS covered_bases,
                    SUM(CASE WHEN site_coverage >= {int(cov_site_threshold)} THEN 1 ELSE 0 END)::BIGINT AS cov_sites
                  FROM profile_rows
                  WHERE gene IS NOT NULL
                    AND CAST(gene AS VARCHAR) <> 'NA'
                    AND genome IS NOT NULL
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
                  END AS coverage,
                  COALESCE(CAST(cs.cov_sites AS BIGINT), 0) AS "{cov_sites_column}",
                  CASE
                    WHEN gl.length > 0 AND COALESCE(cs.covered_bases, 0) > 0
                    THEN
                      (COALESCE(cs.total_covered_sites, 0)::DOUBLE / gl.length) /
                      (1 - EXP(-0.883 * (COALESCE(cs.covered_bases, 0)::DOUBLE / gl.length)))
                    ELSE 0.0
                  END AS ber
                FROM gene_lengths AS gl
                LEFT JOIN covered_stats AS cs
                  ON gl.genome = cs.genome
                 AND gl.gene = cs.gene
                ORDER BY genome, gene
                """
            ).fetch_arrow_table()

            genome_table = conn.execute(
                f"""
                WITH
                {profile_rows_cte},
                genome_lengths AS (
                  SELECT
                    sm.genome AS genome,
                    SUM(
                      CASE
                        WHEN CAST(gb."end" AS BIGINT) >= CAST(gb.start AS BIGINT)
                        THEN CAST(gb."end" AS BIGINT) - CAST(gb.start AS BIGINT)
                        ELSE 0
                      END
                    )::BIGINT AS length
                  FROM genome_bed_src AS gb
                  INNER JOIN scaffold_map_src AS sm
                    ON CAST(gb.scaffold AS VARCHAR) = CAST(sm.scaffold AS VARCHAR)
                  GROUP BY 1
                ),
                covered_stats AS (
                  SELECT
                    CAST(genome AS VARCHAR) AS genome,
                    COUNT(*)::BIGINT AS total_covered_sites,
                    SUM(site_coverage)::DOUBLE AS covered_bases,
                    SUM(CASE WHEN site_coverage >= {int(cov_site_threshold)} THEN 1 ELSE 0 END)::BIGINT AS cov_sites
                  FROM profile_rows
                  WHERE genome IS NOT NULL
                  GROUP BY 1
                )
                SELECT
                  CAST(gl.genome AS VARCHAR) AS genome,
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
                  END AS coverage,
                  COALESCE(CAST(cs.cov_sites AS BIGINT), 0) AS "{cov_sites_column}",
                  CASE
                    WHEN gl.length > 0 AND COALESCE(cs.covered_bases, 0) > 0
                    THEN
                      (COALESCE(cs.total_covered_sites, 0)::DOUBLE / gl.length) /
                      (1 - EXP(-0.883 * (COALESCE(cs.covered_bases, 0)::DOUBLE / gl.length)))
                    ELSE 0.0
                  END AS ber
                FROM genome_lengths AS gl
                LEFT JOIN covered_stats AS cs
                  ON gl.genome = cs.genome
                ORDER BY genome
                """
            ).fetch_arrow_table()
        finally:
            conn.close()

    pq.write_table(gene_table, gene_output, compression="zstd")
    pq.write_table(genome_table, genome_output, compression="zstd")
    return {
        "gene_stats_file": gene_output,
        "genome_stats_file": genome_output,
        "gene_rows": gene_table.num_rows,
        "genome_rows": genome_table.num_rows,
        "cov_sites_column": cov_sites_column,
    }


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


def infer_sample_name_from_stat_table(stat_table: str | pathlib.Path) -> str:
    """Infer sample name from a gene/genome stats parquet file name."""
    stem = pathlib.Path(stat_table).stem
    for suffix in ("_gene_stats", "_genome_stats", "_stats"):
        if stem.endswith(suffix):
            return stem[:-len(suffix)]
    return stem


def infer_sample_name_from_profile(profile_path: str | pathlib.Path) -> str:
    """Infer sample name from a classic profile parquet path."""
    return pathlib.Path(profile_path).stem


def _looks_like_classic_profile_parquet(path: pathlib.Path) -> bool:
    try:
        names = set(pq.read_schema(path).names)
    except Exception:
        return False
    return CLASSIC_PROFILE_REQUIRED_COLUMNS.issubset(names)


def discover_classic_profile_parquets(profile_dir: str | pathlib.Path) -> list[pathlib.Path]:
    """Return classic ZipStrain profile parquets from a directory."""
    input_dir = pathlib.Path(profile_dir)
    if not input_dir.exists():
        raise FileNotFoundError(f"Profile directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Profile path is not a directory: {input_dir}")
    profiles = sorted(
        path for path in input_dir.glob("*.parquet")
        if path.is_file() and _looks_like_classic_profile_parquet(path)
    )
    sample_names = [path.stem for path in profiles]
    if len(sample_names) != len(set(sample_names)):
        raise ValueError("Profile file stems must be unique to derive unique sample names.")
    return profiles


def _scan_pairs_table(table_path: pathlib.Path) -> pl.LazyFrame:
    suffix = table_path.suffix.lower()
    if suffix == ".parquet":
        return pl.scan_parquet(table_path)
    if suffix == ".csv":
        return pl.scan_csv(table_path, separator=",", has_header=True)
    if suffix in {".tsv", ".tab", ".txt"}:
        return pl.scan_csv(table_path, separator="\t", has_header=True)
    raise ValueError(f"Unsupported pair table format for {table_path}. Use parquet, csv, or tsv.")


def _empty_genome_compare_frame(calculate: str | tuple[str, ...]) -> pl.DataFrame:
    from zipstrain import compare as cp

    columns = cp.genome_metric_output_columns(calculate) + ["sample_1", "sample_2"]
    return pl.DataFrame(
        {
            column: pl.Series(name=column, values=[], dtype=GENOME_COMPARE_OUTPUT_DTYPES[column])
            for column in columns
        }
    )


def _normalize_pair_batch(batch: pl.DataFrame) -> pl.DataFrame:
    columns = set(batch.columns)
    if {"sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"}.issubset(columns):
        normalized = batch.select(["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"])
    elif {"sample_name_1", "sample_name_2", "profile_1", "profile_2"}.issubset(columns):
        normalized = (
            batch.rename({"profile_1": "profile_location_1", "profile_2": "profile_location_2"})
            .select(["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"])
        )
    elif {"sample_1", "sample_2", "profile_1", "profile_2"}.issubset(columns):
        normalized = (
            batch.rename(
                {
                    "sample_1": "sample_name_1",
                    "sample_2": "sample_name_2",
                    "profile_1": "profile_location_1",
                    "profile_2": "profile_location_2",
                }
            ).select(["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"])
        )
    elif {"profile_location_1", "profile_location_2"}.issubset(columns):
        normalized = batch.with_columns(
            pl.col("profile_location_1")
            .map_elements(lambda value: infer_sample_name_from_profile(value), return_dtype=pl.Utf8)
            .alias("sample_name_1"),
            pl.col("profile_location_2")
            .map_elements(lambda value: infer_sample_name_from_profile(value), return_dtype=pl.Utf8)
            .alias("sample_name_2"),
        ).select(["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"])
    elif {"profile_1", "profile_2"}.issubset(columns):
        normalized = (
            batch.rename({"profile_1": "profile_location_1", "profile_2": "profile_location_2"})
            .with_columns(
                pl.col("profile_location_1")
                .map_elements(lambda value: infer_sample_name_from_profile(value), return_dtype=pl.Utf8)
                .alias("sample_name_1"),
                pl.col("profile_location_2")
                .map_elements(lambda value: infer_sample_name_from_profile(value), return_dtype=pl.Utf8)
                .alias("sample_name_2"),
            )
            .select(["sample_name_1", "sample_name_2", "profile_location_1", "profile_location_2"])
        )
    else:
        raise ValueError(
            "Pair table must include one of these column sets: "
            "[sample_name_1, sample_name_2, profile_location_1, profile_location_2], "
            "[sample_name_1, sample_name_2, profile_1, profile_2], "
            "[sample_1, sample_2, profile_1, profile_2], "
            "[profile_location_1, profile_location_2], or [profile_1, profile_2]."
        )

    normalized = normalized.with_columns(
        pl.col("sample_name_1").cast(pl.Utf8),
        pl.col("sample_name_2").cast(pl.Utf8),
        pl.col("profile_location_1").cast(pl.Utf8),
        pl.col("profile_location_2").cast(pl.Utf8),
    )
    if normalized.height and any(count > 0 for count in normalized.null_count().row(0)):
        raise ValueError("Pair table contains null values in required columns.")
    return normalized


def generate_genome_pairs(
    profile_dir: str | pathlib.Path,
    output_file: str | pathlib.Path,
    write_batch_size: int = 100_000,
    progress_callback: Callable[[str], None] | None = None,
) -> dict[str, object]:
    """Generate all non-redundant classic profile pairs and write them to parquet."""
    if write_batch_size < 1:
        raise ValueError("write_batch_size must be >= 1")

    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles = discover_classic_profile_parquets(profile_dir)
    sample_rows = [(path.stem, str(path.resolve())) for path in profiles]
    total_profiles = len(sample_rows)
    total_pairs = total_profiles * (total_profiles - 1) // 2

    def emit(message: str) -> None:
        if progress_callback is not None:
            progress_callback(message)

    emit(f"generate_genome_pairs: found {total_profiles} profiles")
    written_pairs = 0
    buffer = {
        "sample_name_1": [],
        "sample_name_2": [],
        "profile_location_1": [],
        "profile_location_2": [],
    }

    with pq.ParquetWriter(output_path, GENOME_PAIR_TABLE_SCHEMA, compression="zstd") as writer:
        if total_pairs == 0:
            writer.write_table(pa.Table.from_arrays([pa.array([], type=field.type) for field in GENOME_PAIR_TABLE_SCHEMA], schema=GENOME_PAIR_TABLE_SCHEMA))
        else:
            for left_idx, (sample_1, profile_1) in enumerate(sample_rows[:-1]):
                for sample_2, profile_2 in sample_rows[left_idx + 1:]:
                    buffer["sample_name_1"].append(sample_1)
                    buffer["sample_name_2"].append(sample_2)
                    buffer["profile_location_1"].append(profile_1)
                    buffer["profile_location_2"].append(profile_2)
                    if len(buffer["sample_name_1"]) >= write_batch_size:
                        writer.write_table(pa.Table.from_pydict(buffer, schema=GENOME_PAIR_TABLE_SCHEMA))
                        written_pairs += len(buffer["sample_name_1"])
                        emit(
                            f"generate_genome_pairs: wrote {written_pairs}/{total_pairs} pairs"
                        )
                        for key in buffer:
                            buffer[key].clear()
            if buffer["sample_name_1"]:
                writer.write_table(pa.Table.from_pydict(buffer, schema=GENOME_PAIR_TABLE_SCHEMA))
                written_pairs += len(buffer["sample_name_1"])
                emit(f"generate_genome_pairs: wrote {written_pairs}/{total_pairs} pairs")

    emit(f"generate_genome_pairs: done -> {output_path}")
    return {
        "output_file": output_path,
        "profiles": total_profiles,
        "pairs": total_pairs,
        "write_batch_size": write_batch_size,
    }


def _schema_signature(schema: pl.Schema) -> tuple[tuple[str, pl.DataType], ...]:
    return tuple((name, schema[name]) for name in schema.names())


def merge_stat_tables(
    stat_tables: list[str | pathlib.Path],
    output_file: str | pathlib.Path,
    progress_callback: Callable[[str], None] | None = None,
) -> pathlib.Path:
    """Concatenate gene/genome stat parquets and add a sample column from file names."""
    if not stat_tables:
        raise ValueError("No stat tables provided.")

    stat_paths = [pathlib.Path(path) for path in stat_tables]
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def emit(message: str) -> None:
        if progress_callback is not None:
            progress_callback(message)

    base_schema: pl.Schema | None = None
    base_signature: tuple[tuple[str, pl.DataType], ...] | None = None
    for path in stat_paths:
        schema = pl.scan_parquet(path).collect_schema()
        signature = _schema_signature(schema)
        if base_signature is None:
            base_schema = schema
            base_signature = signature
            continue
        if signature != base_signature:
            raise ValueError(
                "Stat table schema mismatch. "
                f"Expected {list(base_signature)}, found {list(signature)} in {path}."
            )

    assert base_schema is not None
    output_columns = ["sample", *base_schema.names()]
    emit(f"merge_stat_tables: found {len(stat_paths)} stat tables")
    emit(f"merge_stat_tables: streaming concat start -> {output_path}")
    pl.concat(
        [
            pl.scan_parquet(path)
            .with_columns(pl.lit(infer_sample_name_from_stat_table(path)).alias("sample"))
            .select(output_columns)
            for path in stat_paths
        ],
        how="vertical",
    ).sink_parquet(
        output_path,
        compression="zstd",
        engine="streaming",
    )
    emit(f"merge_stat_tables: streaming concat done -> {output_path}")
    return output_path


def _run_genome_compare_pair(
    row: dict[str, str],
    *,
    stb_file: str | pathlib.Path,
    min_cov: int,
    min_gene_compare_len: int,
    genome_scope: str,
    ani_method: str,
    calculate: str | tuple[str, ...],
    engine: str,
    duckdb_memory_limit: str | None,
    duckdb_temp_directory: str | pathlib.Path | None,
    duckdb_threads: int | None,
) -> tuple[pl.DataFrame, float]:
    from zipstrain import compare as cp

    start = time.perf_counter()
    profile_1 = row["profile_location_1"]
    profile_2 = row["profile_location_2"]
    sample_1 = row["sample_name_1"]
    sample_2 = row["sample_name_2"]
    output_columns = cp.genome_metric_output_columns(calculate) + ["sample_1", "sample_2"]

    profile_1_for_compare = profile_1
    profile_2_for_compare = profile_2
    if engine == "polars" and genome_scope != "all":
        profile_1_for_compare, profile_2_for_compare = cp.duckdb_prefilter_by_scope(
            mpile1=profile_1,
            mpile2=profile_2,
            genome_scope=genome_scope,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )

    frame = (
        cp.compare_genomes(
            mpile_contig_1=profile_1_for_compare,
            mpile_contig_2=profile_2_for_compare,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            ani_method=ani_method,
            duckdb_memory_limit=duckdb_memory_limit,
            duckdb_temp_directory=duckdb_temp_directory,
            duckdb_threads=duckdb_threads,
            engine=engine,
            stb_file=stb_file,
            calculate=calculate,
        )
        .with_columns(
            sample_1=pl.lit(sample_1),
            sample_2=pl.lit(sample_2),
        )
        .select(output_columns)
        .collect(engine="streaming")
    )
    return frame, time.perf_counter() - start


def chunk_genome_compare(
    pair_table: str | pathlib.Path,
    output_file: str | pathlib.Path,
    stb_file: str | pathlib.Path,
    workers: int | None = None,
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    calculate: str | tuple[str, ...] = "all",
    engine: str = "polars",
    duckdb_memory_limit: str | None = None,
    duckdb_temp_directory: str | pathlib.Path | None = None,
    duckdb_threads: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
) -> dict[str, object]:
    """Run classic genome compare on one pair table using in-process parallel workers."""
    from zipstrain import compare as cp

    if workers is not None and workers < 1:
        raise ValueError("workers must be >= 1")
    if engine not in {"polars", "duckdb"}:
        raise ValueError("engine must be one of: polars, duckdb")

    calculations = cp.parse_genome_calculations(calculate)
    pair_table_path = pathlib.Path(pair_table)
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def emit(message: str) -> None:
        if progress_callback is not None:
            progress_callback(message)

    pair_frame = _normalize_pair_batch(_scan_pairs_table(pair_table_path).collect())
    total_pairs = pair_frame.height
    resolved_workers = max(1, min(workers or (os.cpu_count() or 1), max(total_pairs, 1)))
    start_time = time.perf_counter()
    emit(
        f"chunk_genome_compare: start pairs={total_pairs} "
        f"workers={resolved_workers} engine={engine} genome={genome_scope}"
    )

    if total_pairs == 0:
        _empty_genome_compare_frame(calculations).write_parquet(output_path, compression="zstd")
        elapsed = time.perf_counter() - start_time
        emit(f"chunk_genome_compare: done pairs=0 rows=0 elapsed={elapsed:.2f}s")
        return {
            "output_file": output_path,
            "pairs": 0,
            "completed_pairs": 0,
            "rows": 0,
            "workers": resolved_workers,
            "elapsed_seconds": elapsed,
            "avg_wall_seconds_per_pair": 0.0,
            "avg_compute_seconds_per_pair": 0.0,
            "avg_seconds_per_genome_row": 0.0,
            "avg_genome_rows_per_pair": 0.0,
            "engine": engine,
            "genome_scope": genome_scope,
        }

    completed_pairs = 0
    total_rows = 0
    cumulative_pair_compute_seconds = 0.0
    emit("chunk_genome_compare: dispatching pair comparisons")
    ordered_results: list[tuple[pl.DataFrame, float] | None] = [None] * total_pairs
    with ThreadPoolExecutor(max_workers=resolved_workers) as executor:
        future_to_index = {
            executor.submit(
                _run_genome_compare_pair,
                row,
                stb_file=stb_file,
                min_cov=min_cov,
                min_gene_compare_len=min_gene_compare_len,
                genome_scope=genome_scope,
                ani_method=ani_method,
                calculate=calculations,
                engine=engine,
                duckdb_memory_limit=duckdb_memory_limit,
                duckdb_temp_directory=duckdb_temp_directory,
                duckdb_threads=duckdb_threads,
            ): idx
            for idx, row in enumerate(pair_frame.iter_rows(named=True))
        }
        for future in as_completed(future_to_index):
            idx = future_to_index[future]
            frame, pair_elapsed = future.result()
            ordered_results[idx] = (frame, pair_elapsed)
            completed_pairs += 1
            total_rows += frame.height
            cumulative_pair_compute_seconds += pair_elapsed
            emit(
                f"chunk_genome_compare: progress pairs={completed_pairs}/{total_pairs} "
                f"rows={total_rows} elapsed={time.perf_counter() - start_time:.2f}s"
            )

    result_frames: list[pl.DataFrame] = []
    for result in ordered_results:
        assert result is not None
        frame, _pair_elapsed = result
        result_frames.append(frame)

    if result_frames:
        pl.concat(result_frames, how="vertical_relaxed").write_parquet(output_path, compression="zstd")
    else:
        _empty_genome_compare_frame(calculations).write_parquet(output_path, compression="zstd")

    elapsed_seconds = time.perf_counter() - start_time
    summary = {
        "output_file": output_path,
        "pairs": total_pairs,
        "completed_pairs": completed_pairs,
        "rows": total_rows,
        "workers": resolved_workers,
        "elapsed_seconds": elapsed_seconds,
        "avg_wall_seconds_per_pair": elapsed_seconds / max(completed_pairs, 1),
        "avg_compute_seconds_per_pair": cumulative_pair_compute_seconds / max(completed_pairs, 1),
        "avg_seconds_per_genome_row": elapsed_seconds / max(total_rows, 1),
        "avg_genome_rows_per_pair": total_rows / max(completed_pairs, 1),
        "engine": engine,
        "genome_scope": genome_scope,
    }
    emit(
        f"chunk_genome_compare: done pairs={completed_pairs}/{total_pairs} rows={total_rows} "
        f"elapsed={elapsed_seconds:.2f}s"
    )
    return summary

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
