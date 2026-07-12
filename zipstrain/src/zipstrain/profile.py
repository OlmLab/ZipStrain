"""zipstrain.profile
========================
This module provides functions and utilities to profile a bamfile.
By profile we mean generating gene, genome, and nucleotide counts at each position on the reference.
This is a fundamental step for downstream analysis in zipstrain.
"""
import pathlib
import polars as pl
from typing import Generator, Optional
from zipstrain import utils
import asyncio
import os
import re
import shutil
import shlex
import duckdb
import tempfile
import subprocess
import pyarrow as pa
import pyarrow.parquet as pq


PROFILE_SORTED_METADATA_KEY = "zipstrain_sorted_by"
PROFILE_SORTED_METADATA_VALUE = "chrom,pos"
PROFILE_WRITE_BATCH_SIZE = 10_000
MPILEUP_ASYNCIO_STREAM_LIMIT_BYTES = 10 * 1024 * 1024
PROFILE_SUBPROCESS_STDERR_TAIL_BYTES = 64 * 1024
PROFILE_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
PROFILE_MIN_MAPQ_DEFAULT = 0
PROFILE_MIN_BASEQ_DEFAULT = 13
READ_INCLUSION_ALL_MAPPED = "all-mapped"
READ_INCLUSION_PAIRED = "paired"
READ_INCLUSION_PROPER_PAIRS = "proper-pairs"
PROFILE_READ_INCLUSION_CHOICES = (
    READ_INCLUSION_PROPER_PAIRS,
    READ_INCLUSION_PAIRED,
    READ_INCLUSION_ALL_MAPPED,
)
SAM_FLAG_PAIRED = 0x1
SAM_FLAG_PROPER_PAIR = 0x2
SAM_FLAG_UNMAP = 0x4
SAM_FLAG_MUNMAP = 0x8
RAW_PROFILE_PARQUET_FIELDS = [
    ("chrom", pa.string()),
    ("pos", pa.int32()),
    ("A", pa.int32()),
    ("C", pa.int32()),
    ("G", pa.int32()),
    ("T", pa.int32()),
]
PROFILE_PARQUET_BASE_FIELDS = [
    ("chrom", pa.string()),
    ("genome", pa.string()),
    ("gene", pa.string()),
    ("pos", pa.int32()),
    ("A", pa.int32()),
    ("C", pa.int32()),
    ("G", pa.int32()),
    ("T", pa.int32()),
]
MPILEUP_INDEL_RE = re.compile(r"\^.|[\$]|[+-](\d+)")
VCF_SNP_BASE_ORDER = ("A", "C", "G", "T")


async def _read_stream_bounded_tail(
    stream,
    *,
    tail_bytes: int = PROFILE_SUBPROCESS_STDERR_TAIL_BYTES,
    chunk_size: int = 8192,
) -> bytes:
    """Drain an async subprocess stream while retaining only a bounded tail.

    Some tools, notably ``samtools mpileup``, can write warnings to stderr while
    streaming large stdout. If stderr is piped but not drained concurrently, the
    child process can block forever once the OS pipe buffer fills.
    """
    if stream is None:
        return b""
    tail = bytearray()
    while True:
        chunk = await stream.read(chunk_size)
        if not chunk:
            break
        tail.extend(chunk)
        if len(tail) > tail_bytes:
            del tail[: len(tail) - tail_bytes]
    return bytes(tail)


def _shell_join(args: list[str | pathlib.Path]) -> str:
    return " ".join(shlex.quote(str(arg)) for arg in args)


async def _run_profile_shell_pipeline(
    command: str,
    *,
    cwd: pathlib.Path,
    output_file: pathlib.Path | None = None,
) -> None:
    proc = await asyncio.create_subprocess_shell(
        f"bash -o pipefail -c {shlex.quote(command)}",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        limit=MPILEUP_ASYNCIO_STREAM_LIMIT_BYTES,
        cwd=cwd,
    )
    stdout_task = asyncio.create_task(_read_stream_bounded_tail(proc.stdout))
    stderr_task = asyncio.create_task(_read_stream_bounded_tail(proc.stderr))
    try:
        returncode = await proc.wait()
        stdout_tail, stderr_tail = await asyncio.gather(stdout_task, stderr_task)
        if returncode != 0:
            if output_file is not None:
                output_file.unlink(missing_ok=True)
            stderr_text = stderr_tail.decode(errors="replace").strip()
            stdout_text = stdout_tail.decode(errors="replace").strip()
            details = stderr_text or stdout_text or f"exit status {returncode}"
            raise Exception(f"Command failed with error: {details}")
    finally:
        for task in (stdout_task, stderr_task):
            if not task.done():
                task.cancel()
                try:
                    await task
                except asyncio.CancelledError:
                    pass


def _validate_profile_filter_settings(
    *,
    min_mapq: int,
    min_baseq: int,
    min_read_ani: float | None,
    read_inclusion: str,
) -> None:
    if min_mapq < 0:
        raise ValueError("min_mapq must be >= 0")
    if min_baseq < 0:
        raise ValueError("min_baseq must be >= 0")
    if min_read_ani is not None and not (0.0 <= min_read_ani <= 1.0):
        raise ValueError("min_read_ani must be between 0 and 1")
    if read_inclusion not in PROFILE_READ_INCLUSION_CHOICES:
        raise ValueError(
            f"read_inclusion must be one of: {', '.join(PROFILE_READ_INCLUSION_CHOICES)}"
        )


def _aligned_query_bases_from_cigar(cigar: str) -> int:
    if cigar == "*" or not cigar:
        return 0
    total = 0
    for length_str, op in PROFILE_CIGAR_RE.findall(cigar):
        if op in {"M", "I", "=", "X"}:
            total += int(length_str)
    return total


def _extract_nm_from_optional_fields(optional_fields: list[str]) -> int | None:
    for field in optional_fields:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return None


def _sam_alignment_passes_profile_filters(
    *,
    flag: int,
    mapq: int,
    cigar: str,
    optional_fields: list[str],
    min_mapq: int,
    min_read_ani: float | None,
    read_inclusion: str,
) -> bool:
    if flag & SAM_FLAG_UNMAP:
        return False

    if read_inclusion == READ_INCLUSION_PAIRED:
        if not (flag & SAM_FLAG_PAIRED):
            return False
        if flag & SAM_FLAG_MUNMAP:
            return False
    elif read_inclusion == READ_INCLUSION_PROPER_PAIRS:
        if not (flag & SAM_FLAG_PAIRED):
            return False
        if not (flag & SAM_FLAG_PROPER_PAIR):
            return False
        if flag & SAM_FLAG_MUNMAP:
            return False

    if mapq < min_mapq:
        return False

    if min_read_ani is not None:
        nm = _extract_nm_from_optional_fields(optional_fields)
        if nm is None:
            raise ValueError(
                "min_read_ani filtering requires NM tags in the BAM alignments."
            )
        aligned_query_bases = _aligned_query_bases_from_cigar(cigar)
        if aligned_query_bases <= 0:
            return False
        read_ani = 1.0 - (float(nm) / float(aligned_query_bases))
        if read_ani < min_read_ani:
            return False

    return True


def _profile_alignment_filters_require_prefilter(
    *,
    min_mapq: int,
    min_read_ani: float | None,
    read_inclusion: str,
) -> bool:
    return (
        min_mapq != PROFILE_MIN_MAPQ_DEFAULT
        or min_read_ani is not None
        or read_inclusion != READ_INCLUSION_ALL_MAPPED
    )


def _filter_bam_for_profiling(
    *,
    input_bam: pathlib.Path,
    output_bam: pathlib.Path,
    min_mapq: int,
    min_read_ani: float | None,
    read_inclusion: str,
) -> None:
    _validate_profile_filter_settings(
        min_mapq=min_mapq,
        min_baseq=PROFILE_MIN_BASEQ_DEFAULT,
        min_read_ani=min_read_ani,
        read_inclusion=read_inclusion,
    )
    with tempfile.TemporaryFile(mode="w+t") as in_stderr_file, tempfile.TemporaryFile(mode="w+t") as out_stderr_file:
        samtools_view_in = subprocess.Popen(
            ["samtools", "view", "-h", str(input_bam)],
            stdout=subprocess.PIPE,
            stderr=in_stderr_file,
            text=True,
        )
        samtools_view_out = subprocess.Popen(
            ["samtools", "view", "-b", "-o", str(output_bam), "-"],
            stdin=subprocess.PIPE,
            stderr=out_stderr_file,
            text=True,
        )
        assert samtools_view_in.stdout is not None
        assert samtools_view_out.stdin is not None

        try:
            for line in samtools_view_in.stdout:
                if line.startswith("@"):
                    samtools_view_out.stdin.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                if _sam_alignment_passes_profile_filters(
                    flag=int(fields[1]),
                    mapq=int(fields[4]),
                    cigar=fields[5],
                    optional_fields=fields[11:],
                    min_mapq=min_mapq,
                    min_read_ani=min_read_ani,
                    read_inclusion=read_inclusion,
                ):
                    samtools_view_out.stdin.write(line)
        finally:
            samtools_view_in.stdout.close()
            samtools_view_out.stdin.close()

        in_rc = samtools_view_in.wait()
        out_rc = samtools_view_out.wait()
        in_stderr_file.seek(0)
        out_stderr_file.seek(0)
        in_stderr = in_stderr_file.read()
        out_stderr = out_stderr_file.read()
    if in_rc != 0:
        raise RuntimeError(f"samtools view failed while reading BAM: {in_stderr.strip()}")
    if out_rc != 0:
        raise RuntimeError(f"samtools view failed while writing filtered BAM: {out_stderr.strip()}")

    index_proc = subprocess.run(
        ["samtools", "index", str(output_bam)],
        check=False,
        capture_output=True,
        text=True,
    )
    if index_proc.returncode != 0:
        raise RuntimeError(
            f"samtools index failed for filtered BAM: {index_proc.stderr.strip()}"
        )


def raw_profile_parquet_schema(include_reference_base: bool = True) -> pa.Schema:
    fields = list(RAW_PROFILE_PARQUET_FIELDS)
    if include_reference_base:
        fields.append((utils.REF_BASE_BITMASK_COLUMN, pa.uint8()))
    return pa.schema(fields)


def profile_parquet_schema(include_reference_base: bool = True) -> pa.Schema:
    fields = list(PROFILE_PARQUET_BASE_FIELDS)
    if include_reference_base:
        fields.append((utils.REF_BASE_BITMASK_COLUMN, pa.uint8()))
    return pa.schema(fields)

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
    return mpile_frame.join(
        null_model,
        on="cov",
        how="left",
        maintain_order="left",
    ).with_columns([
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


def empty_gene_stats_table(include_ref_ani: bool = False) -> pl.LazyFrame:
    """Return an empty gene-stats LazyFrame with the expected schema."""
    schema = {
        "genome": pl.Utf8,
        "gene": pl.Utf8,
        "length": pl.Int64,
        "breadth": pl.Float64,
        "coverage": pl.Float64,
    }
    if include_ref_ani:
        schema["ref_ani"] = pl.Float64
    return pl.DataFrame(schema=schema).lazy()


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
    mpileup_df=mpileup_df.join(
        scaffold_to_genome,
        left_on="chrom",
        right_on="scaffold",
        how="left",
        maintain_order="left",
    ).with_columns(
        pl.col("genome").fill_null("NA")
    )
    return mpileup_df.set_sorted(["chrom", "pos"])

def add_gene_info_to_mpileup(mpileup_df:pl.LazyFrame, gene_range:pl.LazyFrame)->pl.LazyFrame:
    mpileup_df=mpileup_df.sort(["chrom", "pos"]).set_sorted(["chrom", "pos"])
    gene_range=gene_range.sort(["scaffold", "start"]).set_sorted(["scaffold", "start"])
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
    return annotated_mpileup.set_sorted(["chrom", "pos"])


def _duckdb_quote_sql_string(value: str) -> str:
    """Quote a string literal for embedding in DuckDB SQL."""
    return value.replace("'", "''")


def _write_sorted_profile_with_metadata(
    profile_lf: pl.LazyFrame,
    output_file: pathlib.Path,
    tmp_dir: pathlib.Path,
    metadata: Optional[dict[str, str]] = None,
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

    output_metadata = dict(metadata or {})
    output_metadata[PROFILE_SORTED_METADATA_KEY] = PROFILE_SORTED_METADATA_VALUE

    if _parquet_rows_are_coordinate_sorted(candidate_path):
        pl.scan_parquet(candidate_path).sink_parquet(
            output_file,
            compression="zstd",
            engine="streaming",
            metadata=output_metadata,
        )
    else:
        _sort_existing_profile_parquet(
            input_file=candidate_path,
            output_file=output_file,
            tmp_dir=tmp_dir,
            metadata=output_metadata,
        )


def _sort_existing_profile_parquet(
    input_file: pathlib.Path,
    output_file: pathlib.Path,
    tmp_dir: pathlib.Path,
    metadata: Optional[dict[str, str]] = None,
) -> None:
    """Sort an existing classic profile parquet and attach sortedness metadata."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    sorted_path = tmp_dir / f"{output_file.stem}.sorted.parquet"
    sorted_path.unlink(missing_ok=True)

    conn = duckdb.connect()
    try:
        in_sql = _duckdb_quote_sql_string(str(input_file))
        out_sql = _duckdb_quote_sql_string(str(sorted_path))
        input_columns = pq.read_schema(input_file).names
        select_columns_sql = ",\n                ".join(f'"{column}"' for column in input_columns)
        conn.execute(
            f"""
            COPY (
              SELECT
                {select_columns_sql}
              FROM read_parquet('{in_sql}')
              ORDER BY chrom, pos
            ) TO '{out_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        conn.close()

    output_metadata = utils._read_custom_parquet_metadata(input_file)
    if metadata:
        output_metadata.update(metadata)
    output_metadata[PROFILE_SORTED_METADATA_KEY] = PROFILE_SORTED_METADATA_VALUE

    pl.scan_parquet(sorted_path).sink_parquet(
        output_file,
        compression="zstd",
        engine="streaming",
        metadata=output_metadata,
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


async def _stream_profile_mpileup_chunk_to_parquet(
    *,
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    reference_fasta: pathlib.Path | None,
    output_parquet: pathlib.Path,
    batch_size: int = PROFILE_WRITE_BATCH_SIZE,
    min_mapq: int = PROFILE_MIN_MAPQ_DEFAULT,
    min_baseq: int = PROFILE_MIN_BASEQ_DEFAULT,
) -> None:
    """Stream one mpileup chunk into a raw counts parquet file."""
    include_reference_base = reference_fasta is not None
    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-A",
        "-q",
        str(min_mapq),
        "-Q",
        str(min_baseq),
    ]
    if reference_fasta is not None:
        mpileup_cmd.extend(["-f", str(reference_fasta.absolute())])
    mpileup_cmd.extend([
        "-l",
        str(bed_file.absolute()),
        str(bam_file.absolute()),
    ])
    converter_cmd = [
        "zipstrain",
        "utilities",
        "process_mpileup",
        "--batch-size",
        str(batch_size),
        "--output-file",
        output_parquet.name,
    ]
    if include_reference_base:
        converter_cmd.append("--include-reference-base")

    await _run_profile_shell_pipeline(
        f"{_shell_join(mpileup_cmd)} | {_shell_join(converter_cmd)}",
        cwd=output_parquet.parent,
        output_file=output_parquet,
    )


def _process_raw_profile_chunk_with_polars(
    *,
    raw_parquet: pathlib.Path,
    output_parquet: pathlib.Path,
    null_model: pl.LazyFrame,
    scaffold_to_genome: pl.LazyFrame,
    gene_range: pl.LazyFrame,
    include_reference_base: bool = True,
) -> None:
    """Adjust and annotate one coordinate-sorted raw chunk using Polars."""
    raw_columns = ["chrom", "pos", "A", "C", "G", "T"]
    if include_reference_base:
        raw_columns.append(utils.REF_BASE_BITMASK_COLUMN)

    raw_lf = pl.scan_parquet(raw_parquet).set_sorted(["chrom", "pos"])
    adjusted = adjust_for_sequence_errors(raw_lf, null_model).select(raw_columns)
    with_genome = add_genome_info_to_mpileup(
        mpileup_df=adjusted,
        scaffold_to_genome=scaffold_to_genome,
    )

    gene_schema_names = gene_range.collect_schema().names()
    if gene_schema_names:
        annotated = (
            with_genome.set_sorted(["chrom", "pos"])
            .join_asof(
                gene_range,
                left_on="pos",
                right_on="start",
                by_left="chrom",
                by_right="scaffold",
                strategy="backward",
                check_sortedness=False,
            )
            .with_columns(
                pl.when(pl.col("pos") <= pl.col("end"))
                .then(pl.col("gene"))
                .otherwise(pl.lit("NA"))
                .alias("gene")
            )
            .drop(["start", "end"])
            .set_sorted(["chrom", "pos"])
        )
    else:
        annotated = with_genome.with_columns(pl.lit("NA").alias("gene")).set_sorted(["chrom", "pos"])

    final_columns = ["chrom", "genome", "gene", "pos", "A", "C", "G", "T"]
    if include_reference_base:
        final_columns.append(utils.REF_BASE_BITMASK_COLUMN)

    annotated.select(final_columns).sink_parquet(
        output_parquet,
        compression="zstd",
        engine="streaming",
    )


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


def _reference_sharing_prepared_profile(
    profile: pl.LazyFrame,
    *,
    min_cov: int,
    consumer_name: str,
) -> tuple[pl.LazyFrame, list[str]]:
    """Prepare a reference-aware profile with coverage and shared-reference flags."""
    schema_names = profile.collect_schema().names()
    if utils.REF_BASE_BITMASK_COLUMN not in schema_names:
        raise ValueError(
            "Profile parquet is missing the required ref_base_bitmask column. "
            f"Re-run profiling with --reference-fasta to use {consumer_name}."
        )

    coverage = (pl.col("A") + pl.col("C") + pl.col("G") + pl.col("T")).alias("coverage")
    shares_reference = (
        (
            ((pl.col(utils.REF_BASE_BITMASK_COLUMN) == 1) & (pl.col("A") > 0))
            | ((pl.col(utils.REF_BASE_BITMASK_COLUMN) == 2) & (pl.col("C") > 0))
            | ((pl.col(utils.REF_BASE_BITMASK_COLUMN) == 4) & (pl.col("G") > 0))
            | ((pl.col(utils.REF_BASE_BITMASK_COLUMN) == 8) & (pl.col("T") > 0))
        )
        .cast(pl.Int64)
        .alias("share_ref_site")
    )

    prepared = (
        profile.with_columns(coverage)
        .filter(
            (pl.col("coverage") >= min_cov)
            & (pl.col(utils.REF_BASE_BITMASK_COLUMN) > 0)
        )
        .with_columns(shares_reference)
    )
    return prepared, schema_names


def get_reference_ani(
    profile: pl.LazyFrame,
    *,
    agg_level: str = "genome",
    min_cov: int = 5,
) -> pl.LazyFrame:
    """
    Calculate reference-relative ANI from a profile carrying ``ref_base_bitmask``.

    The reference base is encoded as a one-hot bitmask:
    - A -> 1
    - C -> 2
    - G -> 4
    - T -> 8
    - other / unknown -> 0

    A covered position is treated as sharing the reference allele when the
    observed nucleotide counts contain the reference allele after profile
    sequence-error adjustment.
    """
    prepared, _ = _reference_sharing_prepared_profile(
        profile,
        min_cov=min_cov,
        consumer_name="reference ANI",
    )

    if agg_level == "scaffold":
        group_cols = ["chrom"]
        renamed_cols = {"chrom": "scaffold"}
        sort_cols = ["scaffold"]
    elif agg_level == "genome":
        group_cols = ["genome"]
        renamed_cols = {}
        sort_cols = ["genome"]
    elif agg_level == "gene":
        prepared = prepared.filter(pl.col("gene") != "NA")
        group_cols = ["genome", "chrom", "gene"]
        renamed_cols = {"chrom": "scaffold"}
        sort_cols = ["genome", "scaffold", "gene"]
    else:
        raise ValueError("agg_level must be one of: scaffold, genome, gene")

    result = (
        prepared.group_by(group_cols)
        .agg(
            total_positions=pl.len().cast(pl.Int64),
            share_allele_pos=pl.sum("share_ref_site").cast(pl.Int64),
        )
        .with_columns(
            (pl.col("share_allele_pos") / pl.col("total_positions") * 100.0)
            .cast(pl.Float64)
            .alias("ani_reference")
        )
    )
    if renamed_cols:
        result = result.rename(renamed_cols)
    return result.sort(sort_cols)


def get_reference_snps(
    profile: pl.LazyFrame,
    *,
    min_cov: int = 5,
) -> pl.LazyFrame:
    """
    Return profile-like rows that are SNPs relative to the reference.

    A row is emitted when it passes ``min_cov``, has a known reference base,
    and does not retain the reference allele after profile sequence-error
    adjustment. This matches the reference-sharing logic used by
    ``get_reference_ani``.
    """
    prepared, input_columns = _reference_sharing_prepared_profile(
        profile,
        min_cov=min_cov,
        consumer_name="get-snp-reference",
    )
    return (
        prepared.filter(pl.col("share_ref_site") == 0)
        .select(input_columns)
        .sort(["chrom", "pos"])
    )


def _reference_snp_alt_bases(row: dict[str, object], ref_base: str) -> list[str]:
    return [
        base
        for base in VCF_SNP_BASE_ORDER
        if base != ref_base and int(row[base]) > 0
    ]


def _write_reference_snps_vcf_from_parquet(
    *,
    snp_parquet: pathlib.Path,
    output_file: pathlib.Path,
    min_cov: int,
) -> None:
    parquet_file = pq.ParquetFile(snp_parquet)
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##source=zipstrain get-snp-reference\n")
        handle.write(f"##zipstrain_min_cov={min_cov}\n")
        handle.write(
            '##FILTER=<ID=PASS,Description="Site passes ZipStrain reference SNP filters">\n'
        )
        handle.write(
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total adjusted coverage at the site">\n'
        )
        handle.write(
            '##INFO=<ID=ACGT,Number=4,Type=Integer,Description="Adjusted A,C,G,T counts in the profile row">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for batch in parquet_file.iter_batches():
            for row in pl.from_arrow(batch).iter_rows(named=True):
                ref_base = utils.decode_reference_base_bitmask(row[utils.REF_BASE_BITMASK_COLUMN])
                if ref_base is None:
                    raise ValueError(
                        "Reference SNP VCF export encountered an unknown ref_base_bitmask value."
                    )
                alt_bases = _reference_snp_alt_bases(row, ref_base)
                if not alt_bases:
                    raise ValueError(
                        "Reference SNP VCF export encountered a row without alternate alleles."
                    )
                dp = sum(int(row[base]) for base in VCF_SNP_BASE_ORDER)
                info = (
                    f"DP={dp};ACGT="
                    f"{int(row['A'])},{int(row['C'])},{int(row['G'])},{int(row['T'])}"
                )
                handle.write(
                    "\t".join(
                        [
                            str(row["chrom"]),
                            str(int(row["pos"])),
                            ".",
                            ref_base,
                            ",".join(alt_bases),
                            ".",
                            "PASS",
                            info,
                        ]
                    )
                    + "\n"
                )


def write_reference_snps(
    profile: pl.LazyFrame,
    *,
    output_file: str | pathlib.Path,
    min_cov: int = 5,
    fmt: str = "parquet",
) -> pathlib.Path:
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    snps = get_reference_snps(profile, min_cov=min_cov)
    if fmt == "parquet":
        snps.sink_parquet(output_path, compression="zstd")
        return output_path
    if fmt != "vcf":
        raise ValueError("fmt must be one of: parquet, vcf")

    with tempfile.TemporaryDirectory(prefix=f"{output_path.stem}.snps.", dir=output_path.parent) as work_dir_str:
        work_dir = pathlib.Path(work_dir_str)
        snp_parquet = work_dir / "reference_snps.parquet"
        rendered_vcf = work_dir / output_path.name
        snps.sink_parquet(snp_parquet, compression="zstd")
        _write_reference_snps_vcf_from_parquet(
            snp_parquet=snp_parquet,
            output_file=rendered_vcf,
            min_cov=min_cov,
        )
        os.replace(rendered_vcf, output_path)
    return output_path

def _profile_chunk_paths(
    *,
    bam_file: pathlib.Path,
    output_dir: pathlib.Path,
    chunk_id: int,
) -> tuple[pathlib.Path, pathlib.Path, pathlib.Path, pathlib.Path, pathlib.Path]:
    raw_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.raw.parquet"
    raw_sorted_path = output_dir / f"{bam_file.stem}_{chunk_id}.raw.sorted.parquet"
    candidate_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.candidate.parquet"
    final_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.parquet"
    read_locs_path = output_dir / f"{bam_file.stem}_read_locs_{chunk_id}.parquet"
    return raw_chunk_path, raw_sorted_path, candidate_chunk_path, final_chunk_path, read_locs_path


async def _generate_profile_raw_chunk_task(
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    reference_fasta: pathlib.Path | None,
    output_dir: pathlib.Path,
    chunk_id: int,
    min_mapq: int,
    min_baseq: int,
) -> tuple[pathlib.Path, pathlib.Path]:
    raw_chunk_path, _, _, _, read_locs_path = _profile_chunk_paths(
        bam_file=bam_file,
        output_dir=output_dir,
        chunk_id=chunk_id,
    )
    await _stream_profile_mpileup_chunk_to_parquet(
        bed_file=bed_file,
        bam_file=bam_file,
        reference_fasta=reference_fasta,
        output_parquet=raw_chunk_path,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
    )

    view_cmd = [
        "samtools",
        "view",
        "-F",
        "132",
        "-L",
        str(bed_file.absolute()),
        str(bam_file.absolute()),
    ]
    read_locs_cmd = [
        "zipstrain",
        "utilities",
        "process-read-locs",
        "--output-file",
        read_locs_path.name,
    ]
    await _run_profile_shell_pipeline(
        f"{_shell_join(view_cmd)} | {_shell_join(read_locs_cmd)}",
        cwd=output_dir,
        output_file=read_locs_path,
    )
    return raw_chunk_path, read_locs_path


def _postprocess_profile_raw_chunk(
    *,
    bam_file: pathlib.Path,
    reference_fasta: pathlib.Path | None,
    scaffold_to_genome: pl.LazyFrame,
    gene_range: pl.LazyFrame,
    null_model: pl.LazyFrame,
    output_dir: pathlib.Path,
    chunk_id: int,
) -> pathlib.Path:
    include_reference_base = reference_fasta is not None
    raw_chunk_path, raw_sorted_path, candidate_chunk_path, final_chunk_path, _ = _profile_chunk_paths(
        bam_file=bam_file,
        output_dir=output_dir,
        chunk_id=chunk_id,
    )
    raw_for_processing = raw_chunk_path
    if not _parquet_rows_are_coordinate_sorted(raw_chunk_path):
        _sort_existing_profile_parquet(
            input_file=raw_chunk_path,
            output_file=raw_sorted_path,
            tmp_dir=output_dir,
        )
        raw_for_processing = raw_sorted_path

    _process_raw_profile_chunk_with_polars(
        raw_parquet=raw_for_processing,
        output_parquet=candidate_chunk_path,
        null_model=null_model,
        scaffold_to_genome=scaffold_to_genome,
        gene_range=gene_range,
        include_reference_base=include_reference_base,
    )
    if _parquet_rows_are_coordinate_sorted(candidate_chunk_path):
        os.replace(candidate_chunk_path, final_chunk_path)
    else:
        _sort_existing_profile_parquet(
            input_file=candidate_chunk_path,
            output_file=final_chunk_path,
            tmp_dir=output_dir,
        )
    return final_chunk_path


async def _run_profile_raw_chunk_with_semaphore(
    semaphore: asyncio.Semaphore,
    *,
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    reference_fasta: pathlib.Path | None,
    output_dir: pathlib.Path,
    chunk_id: int,
    min_mapq: int,
    min_baseq: int,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Bound chunk-level subprocess work."""
    async with semaphore:
        return await _generate_profile_raw_chunk_task(
            bed_file=bed_file,
            bam_file=bam_file,
            reference_fasta=reference_fasta,
            output_dir=output_dir,
            chunk_id=chunk_id,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
        )

async def profile_bam_in_chunks(
    bed_file:str,
    bam_file:str,
    reference_fasta: Optional[str],
    gene_range_table: Optional[str],
    stb:pl.LazyFrame,
    null_model:pl.LazyFrame,
    output_dir:str,
    num_chunks:int=24,
    max_concurrency:int=4,
    profile_contract: Optional[dict[str, str]] = None,
    min_mapq: int = PROFILE_MIN_MAPQ_DEFAULT,
    min_baseq: int = PROFILE_MIN_BASEQ_DEFAULT,
    min_read_ani: float | None = None,
    read_inclusion: str = READ_INCLUSION_ALL_MAPPED,
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
    reference_fasta = pathlib.Path(reference_fasta) if reference_fasta is not None else None
    include_reference_base = reference_fasta is not None
    gene_range_table_path = normalize_gene_range_table_path(gene_range_table)
    _validate_profile_filter_settings(
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        min_read_ani=min_read_ani,
        read_inclusion=read_inclusion,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"tmp").mkdir(exist_ok=True)
    filtered_bam_file: pathlib.Path | None = None
    effective_bam_file = bam_file
    if _profile_alignment_filters_require_prefilter(
        min_mapq=min_mapq,
        min_read_ani=min_read_ani,
        read_inclusion=read_inclusion,
    ):
        filtered_bam_file = output_dir / "tmp" / f"{bam_file.stem}.filtered.bam"
        _filter_bam_for_profiling(
            input_bam=bam_file,
            output_bam=filtered_bam_file,
            min_mapq=min_mapq,
            min_read_ani=min_read_ani,
            read_inclusion=read_inclusion,
        )
        effective_bam_file = filtered_bam_file
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
        gene_range_lf = gene_range_lf.sort(["scaffold", "start"]).set_sorted(["scaffold", "start"])
    bed_chunks=utils.split_lf_to_chunks(bed_lf, num_chunks)
    bed_chunk_files=[]
    for chunk_id, bed_file in enumerate(bed_chunks):
        bed_file.sink_csv(output_dir/"tmp"/f"bed_chunk_{chunk_id}.bed",include_header=False,separator="\t")
        bed_chunk_files.append(output_dir/"tmp"/f"bed_chunk_{chunk_id}.bed")
    raw_tasks = []
    semaphore = asyncio.Semaphore(max(1, int(max_concurrency)))
    for chunk_id, bed_chunk_file in enumerate(bed_chunk_files):
        raw_tasks.append(
            _run_profile_raw_chunk_with_semaphore(
                semaphore,
                bed_file=bed_chunk_file,
                bam_file=effective_bam_file,
                reference_fasta=reference_fasta,
                output_dir=output_dir / "tmp",
                chunk_id=chunk_id,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
            )
        )
    await asyncio.gather(*raw_tasks)

    for chunk_id in range(len(bed_chunk_files)):
        _postprocess_profile_raw_chunk(
            bam_file=effective_bam_file,
            reference_fasta=reference_fasta,
            scaffold_to_genome=stb,
            gene_range=gene_range_lf,
            null_model=null_model,
            output_dir=output_dir / "tmp",
            chunk_id=chunk_id,
        )
    pfs = [
        (
            output_dir / "tmp" / f"{effective_bam_file.stem}_{chunk_id}.parquet",
            output_dir / "tmp" / f"{effective_bam_file.stem}_read_locs_{chunk_id}.parquet",
        )
        for chunk_id in range(len(bed_chunk_files))
        if (output_dir / "tmp" / f"{effective_bam_file.stem}_{chunk_id}.parquet").exists()
    ]

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
            metadata=utils.profile_contract_metadata_from_values(profile_contract),
        )
        if gene_range_table_path is None:
            empty_gene_stats_table(include_ref_ani=include_reference_base).sink_parquet(
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
    reference_fasta: Optional[str],
    gene_range_table: Optional[str],
    stb:pl.LazyFrame,
    null_model:pl.LazyFrame,
    output_dir:str,
    num_chunks:int=24,
    max_concurrency:int=4,
    profile_contract: Optional[dict[str, str]] = None,
    min_mapq: int = PROFILE_MIN_MAPQ_DEFAULT,
    min_baseq: int = PROFILE_MIN_BASEQ_DEFAULT,
    min_read_ani: float | None = None,
    read_inclusion: str = READ_INCLUSION_ALL_MAPPED,
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
        reference_fasta=reference_fasta,
        gene_range_table=gene_range_table,
        stb=stb,
        null_model=null_model,
        output_dir=output_dir,
        num_chunks=num_chunks,
        max_concurrency=max_concurrency,
        profile_contract=profile_contract,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        min_read_ani=min_read_ani,
        read_inclusion=read_inclusion,
    ))
