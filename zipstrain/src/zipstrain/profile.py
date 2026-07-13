"""zipstrain.profile
========================
This module provides functions and utilities to profile a bamfile.
By profile we mean generating gene, genome, and nucleotide counts at each position on the reference.
This is a fundamental step for downstream analysis in zipstrain.
"""
import pathlib
import polars as pl
from typing import Generator, Optional
from dataclasses import dataclass
from zipstrain import utils
import asyncio
import json
import os
import re
import shutil
import duckdb
from concurrent.futures import ThreadPoolExecutor
import tempfile
import subprocess
import pyarrow as pa
import pyarrow.parquet as pq


PROFILE_SORTED_METADATA_KEY = "zipstrain_sorted_by"
PROFILE_SORTED_METADATA_VALUE = "chrom,pos"
PROFILE_WRITE_BATCH_SIZE = 10_000
MPILEUP_ASYNCIO_STREAM_LIMIT_BYTES = 10 * 1024 * 1024
PROFILE_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
PROFILE_MIN_MAPQ_DEFAULT = 0
PROFILE_MIN_BASEQ_DEFAULT = 13
# Read-ANI floor applied by default so low-identity (typically mis-mapped) reads
# do not inflate SNV calls; matches inStrain's default of 0.95.
PROFILE_MIN_READ_ANI_DEFAULT = 0.95
READ_INCLUSION_ALL_MAPPED = "all-mapped"
READ_INCLUSION_PAIRED = "paired"
READ_INCLUSION_PROPER_PAIRS = "proper-pairs"
PROFILE_READ_INCLUSION_CHOICES = (
    READ_INCLUSION_PROPER_PAIRS,
    READ_INCLUSION_PAIRED,
    READ_INCLUSION_ALL_MAPPED,
)
# Default eligibility: for paired reads keep only those whose mate maps to the same
# scaffold (matching inStrain's default "paired_only": drops half-mapped orphans and
# cross-scaffold pairs), while keeping genuinely single-end reads.
PROFILE_READ_INCLUSION_DEFAULT = READ_INCLUSION_PAIRED
SAM_FLAG_PAIRED = 0x1
SAM_FLAG_PROPER_PAIR = 0x2
SAM_FLAG_UNMAP = 0x4
SAM_FLAG_MUNMAP = 0x8
VCF_SNP_BASE_ORDER = ("A", "C", "G", "T")


def read_stb(stb_file) -> pl.LazyFrame:
    """Scan a scaffold-to-genome (STB) TSV as ``scaffold, genome``.

    Leading/trailing whitespace around the columns is stripped so scaffold names
    match the BAM/FASTA exactly. Some STB files (which inStrain tolerates) carry
    stray spaces around the tab; without stripping, those scaffolds fail to join
    and are silently dropped to genome ``NA``.
    """
    return pl.scan_csv(stb_file, separator="\t", has_header=False).select(
        pl.col("column_1").cast(pl.Utf8).str.strip_chars().alias("scaffold"),
        pl.col("column_2").cast(pl.Utf8).str.strip_chars().alias("genome"),
    )
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
    rname: str = "*",
    rnext: str = "*",
) -> bool:
    if flag & SAM_FLAG_UNMAP:
        return False

    if read_inclusion == READ_INCLUSION_PAIRED:
        # For paired reads, match inStrain's default 'paired_only': keep a read only if
        # its mate is mapped to the SAME scaffold. This drops half-mapped orphans (mate
        # unmapped) and pairs split across scaffolds, both of which otherwise inflate
        # coverage relative to inStrain. Mate-same-scaffold is encoded in SAM as
        # RNEXT == "=" (or RNEXT == RNAME when the mapper writes the name out).
        # Genuinely single-end reads (not flagged paired) can't be part of a pair, so
        # they are kept, matching this mode's single-end-safe behavior.
        if flag & SAM_FLAG_PAIRED:
            if flag & SAM_FLAG_MUNMAP:
                return False
            if rnext != "=" and rnext != rname:
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
            # No NM tag on this alignment: we can't estimate read ANI, so leave
            # the read unfiltered rather than failing the run. (BAMs from a
            # mapper that omits NM simply don't get ANI filtering.)
            return True
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
    samtools_view_in = subprocess.Popen(
        ["samtools", "view", "-h", str(input_bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    samtools_view_out = subprocess.Popen(
        ["samtools", "view", "-b", "-o", str(output_bam), "-"],
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
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
                rname=fields[2],
                rnext=fields[6],
            ):
                samtools_view_out.stdin.write(line)
    finally:
        samtools_view_in.stdout.close()
        samtools_view_out.stdin.close()

    in_stderr = samtools_view_in.stderr.read() if samtools_view_in.stderr is not None else ""
    out_stderr = samtools_view_out.stderr.read() if samtools_view_out.stderr is not None else ""
    in_rc = samtools_view_in.wait()
    out_rc = samtools_view_out.wait()
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


def _stream_profile_mpileup_chunk_to_parquet(
    *,
    bed_file: pathlib.Path,
    bam_file: pathlib.Path,
    reference_fasta: pathlib.Path | None,
    output_parquet: pathlib.Path,
    batch_size: int = PROFILE_WRITE_BATCH_SIZE,
    min_mapq: int = PROFILE_MIN_MAPQ_DEFAULT,
    min_baseq: int = PROFILE_MIN_BASEQ_DEFAULT,
) -> None:
    """Stream one mpileup chunk into a raw counts parquet file.

    Runs synchronously with ``subprocess`` (not asyncio subprocess): asyncio's default
    ThreadedChildWatcher spawns one waitpid thread per child, and this profiler launches
    many short-lived samtools processes per chunk, so under concurrent invocation those
    threads accumulate (100s of them) and asyncio.run()'s shutdown deadlocks. ``subprocess``
    reaps its own children via waitpid, so no child watcher is involved. Callers run this
    in a ThreadPoolExecutor worker to keep chunk-level concurrency.
    """
    include_reference_base = reference_fasta is not None
    cmd = [
        "samtools",
        "mpileup",
        "-A",
        "-q",
        str(min_mapq),
        "-Q",
        str(min_baseq),
    ]
    if reference_fasta is not None:
        cmd.extend(["-f", str(reference_fasta.absolute())])
    cmd.extend([
        "-l",
        str(bed_file.absolute()),
        str(bam_file.absolute()),
    ])

    chroms: list[str] = []
    positions: list[int] = []
    As: list[int] = []
    Cs: list[int] = []
    Gs: list[int] = []
    Ts: list[int] = []
    ref_base_bitmasks: list[int] = []
    parquet_schema = raw_profile_parquet_schema(include_reference_base=include_reference_base)
    writer = pq.ParquetWriter(str(output_parquet), parquet_schema, compression="zstd")
    success = False

    def flush_batch() -> None:
        if not chroms:
            return
        arrays = [
            pa.array(chroms, type=pa.string()),
            pa.array(positions, type=pa.int32()),
            pa.array(As, type=pa.int32()),
            pa.array(Cs, type=pa.int32()),
            pa.array(Gs, type=pa.int32()),
            pa.array(Ts, type=pa.int32()),
        ]
        if include_reference_base:
            arrays.append(pa.array(ref_base_bitmasks, type=pa.uint8()))
        batch = pa.RecordBatch.from_arrays(arrays, schema=parquet_schema)
        writer.write_table(pa.Table.from_batches([batch]))
        chroms.clear()
        positions.clear()
        As.clear()
        Cs.clear()
        Gs.clear()
        Ts.clear()
        ref_base_bitmasks.clear()

    # stderr goes to a temp file (not a pipe) so it can never fill a pipe buffer and
    # deadlock the stdout stream; it is only read back if samtools exits non-zero.
    try:
        with tempfile.TemporaryFile() as errf:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf)
            try:
                for raw_line in proc.stdout:
                    line = raw_line.decode().strip()
                    if not line:
                        continue
                    fields = line.split("\t")
                    if len(fields) < 5:
                        continue
                    chrom, pos_str, ref_base, _, bases = fields[:5]
                    counts = utils.count_mpileup_bases(bases, ref_base, MPILEUP_INDEL_RE)
                    chroms.append(chrom)
                    positions.append(int(pos_str))
                    As.append(counts["A"])
                    Cs.append(counts["C"])
                    Gs.append(counts["G"])
                    Ts.append(counts["T"])
                    if include_reference_base:
                        ref_base_bitmasks.append(utils.encode_reference_base_bitmask(ref_base))
                    if len(chroms) >= batch_size:
                        flush_batch()
            finally:
                proc.stdout.close()
            if proc.wait() != 0:
                errf.seek(0)
                err_text = errf.read().decode(errors="replace").strip()
                raise Exception(f"Command failed with error: {err_text}")
            flush_batch()
            success = True
    finally:
        writer.close()
        if not success:
            output_parquet.unlink(missing_ok=True)


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


SNV_TABLE_COLUMNS = [
    "chrom", "genome", "gene", "pos", "position_coverage", "allele_count",
    "ref_base", "con_base", "var_base", "ref_freq", "con_freq", "var_freq",
    "A", "C", "G", "T", "cryptic", "class", "ref_base_bitmask",
]


def get_reference_snps(
    profile: pl.LazyFrame,
    *,
    min_cov: int = 5,
) -> pl.LazyFrame:
    """Classify variant sites relative to the reference (inStrain-parity table).

    One row per covered position that is *divergent* — i.e. not a monomorphic
    match to the reference. Because the profile counts are already sequence-error
    adjusted, a base with ``count > 0`` is a "passing" allele. With
    ``alleles`` = the passing bases, ``con`` = the consensus (most common) base
    and ``ref`` = the reference base, each site is labelled:

    - ``SNS``     : one allele, ``con != ref`` (fixed substitution)
    - ``SNV``     : ≥2 alleles, ``con == ref`` (reference is majority; minor variant)
    - ``con_SNV`` : ≥2 alleles, ``con != ref``, reference still among the alleles
    - ``pop_SNV`` : ≥2 alleles, ``con != ref``, reference absent from the alleles

    Monomorphic reference sites (one allele equal to the reference) are omitted.
    Frequencies are over the (error-adjusted) position coverage. ``cryptic`` is
    always ``False`` (ZipStrain profiles a single mismatch level). Requires a
    profile carrying ``ref_base_bitmask`` (i.e. ``--reference-fasta``).
    """
    schema_names = profile.collect_schema().names()
    if utils.REF_BASE_BITMASK_COLUMN not in schema_names:
        raise ValueError(
            "Profile parquet is missing the required ref_base_bitmask column. "
            "Re-run profiling with --reference-fasta to call SNVs."
        )
    bm = pl.col(utils.REF_BASE_BITMASK_COLUMN)
    cov = (pl.col("A") + pl.col("C") + pl.col("G") + pl.col("T"))
    base_order = ["A", "C", "G", "T"]
    max_count = pl.max_horizontal("A", "C", "G", "T")

    def _first_base(cond_counts):
        # Pick the first base (A,C,G,T order) whose count equals the target.
        expr = pl.lit(None, dtype=pl.Utf8)
        for base, val in reversed(cond_counts):
            expr = pl.when(val).then(pl.lit(base)).otherwise(expr)
        return expr

    prepared = (
        profile.with_columns(cov.alias("position_coverage"))
        .filter((pl.col("position_coverage") >= min_cov) & (bm > 0))
        .with_columns(
            ref_base=pl.when(bm == 1).then(pl.lit("A")).when(bm == 2).then(pl.lit("C"))
            .when(bm == 4).then(pl.lit("G")).when(bm == 8).then(pl.lit("T")).otherwise(pl.lit("N")),
            ref_count=pl.when(bm == 1).then(pl.col("A")).when(bm == 2).then(pl.col("C"))
            .when(bm == 4).then(pl.col("G")).when(bm == 8).then(pl.col("T")).otherwise(pl.lit(0)),
            con_base=_first_base([(b, pl.col(b) == max_count) for b in base_order]),
            con_count=max_count,
            allele_count=((pl.col("A") > 0).cast(pl.Int64) + (pl.col("C") > 0).cast(pl.Int64)
                          + (pl.col("G") > 0).cast(pl.Int64) + (pl.col("T") > 0).cast(pl.Int64)),
        )
        # var_base = highest-count base other than the consensus.
        .with_columns(
            _va=pl.when(pl.col("con_base") != "A").then(pl.col("A")).otherwise(pl.lit(-1)),
            _vc=pl.when(pl.col("con_base") != "C").then(pl.col("C")).otherwise(pl.lit(-1)),
            _vg=pl.when(pl.col("con_base") != "G").then(pl.col("G")).otherwise(pl.lit(-1)),
            _vt=pl.when(pl.col("con_base") != "T").then(pl.col("T")).otherwise(pl.lit(-1)),
        )
        .with_columns(_var_max=pl.max_horizontal("_va", "_vc", "_vg", "_vt"))
        .with_columns(
            var_base=_first_base([
                ("A", pl.col("_va") == pl.col("_var_max")),
                ("C", pl.col("_vc") == pl.col("_var_max")),
                ("G", pl.col("_vg") == pl.col("_var_max")),
                ("T", pl.col("_vt") == pl.col("_var_max")),
            ]),
            var_count=pl.max_horizontal(pl.col("_var_max"), pl.lit(0)),
        )
        .with_columns(
            ref_freq=pl.col("ref_count") / pl.col("position_coverage"),
            con_freq=pl.col("con_count") / pl.col("position_coverage"),
            var_freq=pl.col("var_count") / pl.col("position_coverage"),
            cryptic=pl.lit(False),
            **{"class": pl.when((pl.col("allele_count") == 1) & (pl.col("con_base") != pl.col("ref_base")))
                .then(pl.lit("SNS"))
                .when((pl.col("allele_count") >= 2) & (pl.col("con_base") == pl.col("ref_base")))
                .then(pl.lit("SNV"))
                .when((pl.col("allele_count") >= 2) & (pl.col("con_base") != pl.col("ref_base")) & (pl.col("ref_count") > 0))
                .then(pl.lit("con_SNV"))
                .when((pl.col("allele_count") >= 2) & (pl.col("con_base") != pl.col("ref_base")) & (pl.col("ref_count") == 0))
                .then(pl.lit("pop_SNV"))
                .otherwise(pl.lit(None))},
        )
        .filter(pl.col("class").is_not_null())
    )
    return prepared.select(SNV_TABLE_COLUMNS).sort(["chrom", "pos"])



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


def _profile_chunk_task(
    bed_file:pathlib.Path,
    bam_file:pathlib.Path,
    reference_fasta: pathlib.Path | None,
    scaffold_to_genome: pl.LazyFrame,
    gene_range: pl.LazyFrame,
    null_model: pl.LazyFrame,
    output_dir:pathlib.Path,
    chunk_id:int,
    min_mapq: int,
    min_baseq: int,
)->None:
    include_reference_base = reference_fasta is not None
    raw_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.raw.parquet"
    raw_sorted_path = output_dir / f"{bam_file.stem}_{chunk_id}.raw.sorted.parquet"
    candidate_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.candidate.parquet"
    final_chunk_path = output_dir / f"{bam_file.stem}_{chunk_id}.parquet"

    _stream_profile_mpileup_chunk_to_parquet(
        bed_file=bed_file,
        bam_file=bam_file,
        reference_fasta=reference_fasta,
        output_parquet=raw_chunk_path,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
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

    # Synchronous shell pipe (samtools view | process-read-locs); absolute output path
    # so cwd is not needed. Run via subprocess (not asyncio) to avoid the child-watcher
    # thread accumulation described in _stream_profile_mpileup_chunk_to_parquet.
    read_locs_out = (output_dir / f"{bam_file.stem}_read_locs_{chunk_id}.parquet").absolute()
    cmd=["samtools", "view", "-F", "132", "-L", str(bed_file.absolute()), str(bam_file.absolute()), "|", "zipstrain", "utilities", "process-read-locs", "--output-file", str(read_locs_out)]
    result = subprocess.run(
        " ".join(cmd),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise Exception(f"Command failed with error: {result.stderr.decode().strip()}")

def _ensure_faidx(reference_fasta: pathlib.Path) -> None:
    """Create <reference_fasta>.fai once (idempotent) so concurrent chunk mpileups
    don't race to build it. No cwd is passed so subprocess can use posix_spawn."""
    fai = reference_fasta.with_name(reference_fasta.name + ".fai")
    if fai.exists():
        return
    result = subprocess.run(
        ["samtools", "faidx", str(reference_fasta.absolute())],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0 and not fai.exists():
        raise Exception(f"samtools faidx failed: {result.stderr.decode().strip()}")


def profile_bam_in_chunks(
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

    # Build the FASTA index once up front. Otherwise the per-chunk samtools mpileup
    # processes race to create <reference>.fai concurrently and corrupt it.
    if reference_fasta is not None:
        _ensure_faidx(reference_fasta)

    # Run chunks in a thread pool bounded by max_concurrency. Each chunk task is
    # synchronous and manages its own subprocesses via `subprocess` (which reaps its
    # own children); this deliberately avoids asyncio subprocess, whose default
    # ThreadedChildWatcher spawns a waitpid thread per child and deadlocks
    # asyncio.run()'s shutdown once hundreds accumulate under concurrent invocation.
    with ThreadPoolExecutor(max_workers=max(1, int(max_concurrency))) as executor:
        futures = [
            executor.submit(
                _profile_chunk_task,
                bed_file=bed_chunk_file,
                bam_file=effective_bam_file,
                reference_fasta=reference_fasta,
                scaffold_to_genome=stb,
                gene_range=gene_range_lf,
                null_model=null_model,
                output_dir=output_dir/"tmp",
                chunk_id=chunk_id,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
            )
            for chunk_id, bed_chunk_file in enumerate(bed_chunk_files)
        ]
        for future in futures:
            future.result()
    pfs=[(output_dir/"tmp"/f"{effective_bam_file.stem}_{chunk_id}.parquet", output_dir/"tmp"/f"{effective_bam_file.stem}_read_locs_{chunk_id}.parquet" ) for chunk_id in range(len(bed_chunk_files)) if (output_dir/"tmp"/f"{effective_bam_file.stem}_{chunk_id}.parquet").exists()]

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
    profile_bam_in_chunks(
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
    )


# ---------------------------------------------------------------------------
# Profiling-asset preparation and caching
#
# These helpers let `zipstrain profile` build the intermediate files it needs
# (bed, gene range table, genome lengths, null model, profiling contract)
# on the fly, instead of forcing the user to run `prepare_profiling` first.
# ---------------------------------------------------------------------------

# File names used inside the profiling-assets directory. These match the names
# produced by the `prepare_profiling` utility so the two paths are compatible.
ASSET_BED_FILENAME = "genomes_bed_file.bed"
ASSET_GENE_RANGE_FILENAME = "gene_range_table.tsv"
ASSET_GENOME_LENGTH_FILENAME = "genome_lengths.parquet"
ASSET_NULL_MODEL_FILENAME = "null_model.parquet"
ASSET_CONTRACT_FILENAME = "profiling_contract.json"
ASSET_CACHE_MANIFEST_FILENAME = "cache_manifest.json"

DEFAULT_PROFILING_ASSETS_DIRNAME = "profiling_assets"


@dataclass
class ProfilingAssets:
    """Concrete paths to the intermediate files needed to profile a BAM.

    ``gene_range_table`` and ``profiling_contract_file`` may be ``None`` when
    the caller supplied explicit assets and did not provide those optional
    files; the auto-preparation path always populates every field.
    """

    bed_file: pathlib.Path
    gene_range_table: pathlib.Path | None
    genome_length_file: pathlib.Path
    null_model_file: pathlib.Path
    profiling_contract_file: pathlib.Path | None


def _build_null_model_frame(
    error_rate: float, max_total_reads: int, p_threshold: float, model_type: str
) -> pl.DataFrame:
    """Build the null-model DataFrame for the requested model type."""
    if model_type == "poisson":
        rows = utils.build_null_poisson(error_rate, max_total_reads, p_threshold)
    else:
        raise ValueError(f"Unsupported model type: {model_type}")
    return pl.DataFrame(rows, schema=["cov", "max_error_count"], orient="row")


def prepare_profiling_assets(
    *,
    reference_fasta: str | pathlib.Path,
    stb_file: str | pathlib.Path,
    output_dir: str | pathlib.Path,
    gene_fasta: str | pathlib.Path | None = None,
    error_rate: float = 0.001,
    max_total_reads: int = 10000,
    p_threshold: float = 0.05,
    model_type: str = "poisson",
) -> ProfilingAssets:
    """Build every intermediate profiling asset into ``output_dir``.

    This is the shared implementation behind the ``prepare_profiling`` utility
    and the auto-preparation performed by ``zipstrain profile``. It writes the
    bed file, gene range table, genome length table, null model, and profiling
    contract, and returns their paths.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bed_path = output_dir / ASSET_BED_FILENAME
    bed_df = utils.make_the_bed(reference_fasta)
    bed_df.write_csv(bed_path, separator="\t", include_header=False)

    gene_range_path = output_dir / ASSET_GENE_RANGE_FILENAME
    if gene_fasta is None:
        empty_gene_range_table().sink_csv(
            gene_range_path, separator="\t", include_header=False
        )
    else:
        build_gene_range_table(pathlib.Path(gene_fasta)).write_csv(
            gene_range_path, separator="\t", include_header=False
        )

    stb = read_stb(stb_file)
    genome_length_path = output_dir / ASSET_GENOME_LENGTH_FILENAME
    utils.extract_genome_length(stb, bed_df.lazy()).sink_parquet(
        genome_length_path, compression="zstd"
    )

    null_model_path = output_dir / ASSET_NULL_MODEL_FILENAME
    _build_null_model_frame(
        error_rate, max_total_reads, p_threshold, model_type
    ).write_parquet(null_model_path)

    contract = {
        "reference_hash": utils.sha256_file(reference_fasta),
        "gene_hash": utils.sha256_file(gene_fasta)
        if gene_fasta is not None
        else utils.PROFILE_CONTRACT_MISSING_VALUE,
        "stb_hash": utils.sha256_file(stb_file),
        "null_model_hash": utils.sha256_file(null_model_path),
    }
    contract_path = output_dir / ASSET_CONTRACT_FILENAME
    utils.write_profile_contract_file(contract, contract_path)

    return ProfilingAssets(
        bed_file=bed_path,
        gene_range_table=gene_range_path,
        genome_length_file=genome_length_path,
        null_model_file=null_model_path,
        profiling_contract_file=contract_path,
    )


def _file_signature(path: str | pathlib.Path) -> dict[str, object]:
    """Return a lightweight (path, size, mtime) signature for cache checks."""
    stat = os.stat(path)
    return {
        "path": str(pathlib.Path(path).resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _build_cache_manifest(
    *,
    reference_fasta: str | pathlib.Path,
    stb_file: str | pathlib.Path,
    gene_fasta: str | pathlib.Path | None,
    error_rate: float,
    max_total_reads: int,
    p_threshold: float,
    model_type: str,
) -> dict[str, object]:
    """Build the cache manifest describing the inputs used to generate assets."""
    return {
        "reference": _file_signature(reference_fasta),
        "stb": _file_signature(stb_file),
        "gene": _file_signature(gene_fasta) if gene_fasta is not None else None,
        "null_model_params": {
            "error_rate": error_rate,
            "max_total_reads": max_total_reads,
            "p_threshold": p_threshold,
            "model_type": model_type,
        },
    }


def _assets_from_dir(assets_dir: pathlib.Path) -> ProfilingAssets:
    return ProfilingAssets(
        bed_file=assets_dir / ASSET_BED_FILENAME,
        gene_range_table=assets_dir / ASSET_GENE_RANGE_FILENAME,
        genome_length_file=assets_dir / ASSET_GENOME_LENGTH_FILENAME,
        null_model_file=assets_dir / ASSET_NULL_MODEL_FILENAME,
        profiling_contract_file=assets_dir / ASSET_CONTRACT_FILENAME,
    )


def _cached_assets_are_valid(
    assets_dir: pathlib.Path, expected_manifest: dict[str, object]
) -> bool:
    """True if the cached assets exist and match the expected input manifest."""
    manifest_path = assets_dir / ASSET_CACHE_MANIFEST_FILENAME
    if not manifest_path.exists():
        return False
    try:
        stored = json.loads(manifest_path.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    if stored != expected_manifest:
        return False
    assets = _assets_from_dir(assets_dir)
    return all(
        path.exists()
        for path in (
            assets.bed_file,
            assets.gene_range_table,
            assets.genome_length_file,
            assets.null_model_file,
            assets.profiling_contract_file,
        )
    )


def resolve_profiling_assets(
    *,
    run_dir: str | pathlib.Path,
    reference_fasta: str | pathlib.Path | None,
    stb_file: str | pathlib.Path,
    gene_fasta: str | pathlib.Path | None = None,
    null_model_file: str | pathlib.Path | None = None,
    bed_file: str | pathlib.Path | None = None,
    genome_length_file: str | pathlib.Path | None = None,
    gene_range_table: str | pathlib.Path | None = None,
    profiling_contract_file: str | pathlib.Path | None = None,
    error_rate: float = 0.001,
    max_total_reads: int = 10000,
    p_threshold: float = 0.05,
    model_type: str = "poisson",
    force_prepare: bool = False,
) -> ProfilingAssets:
    """Resolve every profiling asset, auto-generating any that weren't supplied.

    Explicitly-provided paths always win. When any of the required assets
    (bed, genome length, null model) is missing, the full asset set is built
    into ``run_dir/profiling_assets`` (reusing a valid cached copy when the
    inputs and null-model parameters are unchanged and ``force_prepare`` is
    False), and any explicitly-provided paths override the generated ones.
    """
    needs_generation = force_prepare or any(
        value is None
        for value in (null_model_file, bed_file, genome_length_file)
    )

    generated: ProfilingAssets | None = None
    if needs_generation:
        if reference_fasta is None:
            raise ValueError(
                "--reference-fasta is required to auto-generate profiling assets "
                "(the bed and genome-length files are derived from it). Provide it, "
                "or pass pre-built --bed-file and --genome-length-file."
            )
        assets_dir = pathlib.Path(run_dir) / DEFAULT_PROFILING_ASSETS_DIRNAME
        expected_manifest = _build_cache_manifest(
            reference_fasta=reference_fasta,
            stb_file=stb_file,
            gene_fasta=gene_fasta,
            error_rate=error_rate,
            max_total_reads=max_total_reads,
            p_threshold=p_threshold,
            model_type=model_type,
        )
        if not force_prepare and _cached_assets_are_valid(assets_dir, expected_manifest):
            generated = _assets_from_dir(assets_dir)
        else:
            generated = prepare_profiling_assets(
                reference_fasta=reference_fasta,
                stb_file=stb_file,
                output_dir=assets_dir,
                gene_fasta=gene_fasta,
                error_rate=error_rate,
                max_total_reads=max_total_reads,
                p_threshold=p_threshold,
                model_type=model_type,
            )
            (assets_dir / ASSET_CACHE_MANIFEST_FILENAME).write_text(
                json.dumps(expected_manifest, indent=2, sort_keys=True) + "\n"
            )

    def _resolve(explicit, generated_path):
        if explicit is not None:
            return pathlib.Path(explicit)
        return generated_path

    generated_bed = generated.bed_file if generated else None
    generated_gene_range = generated.gene_range_table if generated else None
    generated_genome_length = generated.genome_length_file if generated else None
    generated_null_model = generated.null_model_file if generated else None
    generated_contract = generated.profiling_contract_file if generated else None

    return ProfilingAssets(
        bed_file=_resolve(bed_file, generated_bed),
        gene_range_table=_resolve(gene_range_table, generated_gene_range),
        genome_length_file=_resolve(genome_length_file, generated_genome_length),
        null_model_file=_resolve(null_model_file, generated_null_model),
        profiling_contract_file=_resolve(profiling_contract_file, generated_contract),
    )
