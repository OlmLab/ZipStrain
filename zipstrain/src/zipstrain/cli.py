"""
zipstrain.utils
========================
This module contains the command-line interface (CLI) implementation for the zipstrain application.
"""
import click as click
from zipstrain import __version__
import zipstrain.utils as ut
import zipstrain.compare as cp
import zipstrain.profile as pf
import zipstrain.task_manager as tm
import zipstrain.database as db
import zipstrain.build_db as bdb
import zipstrain.matrix_pairs as mp
import zipstrain.healthcheck as hc
import zipstrain.mapping as mapping
import zipstrain.matrix_workflow as matrix_workflow
from zipstrain.run_logger import RunLogger
import polars as pl
import pathlib
import shutil
import sys
import time
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import box
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn

DEFAULT_CONTAINER_REPOSITORY = "parsaghadermazi/zipstrain"


class SectionedCommand(click.Command):
    """A click.Command that renders its options grouped under section headers.

    Set ``option_sections`` to an ordered mapping of ``section title -> [param
    names]``. Options are listed under their section in the given order; any
    option not assigned to a section (e.g. ``--help``) is shown last under
    "Other options".
    """

    option_sections: dict[str, list[str]] = {}

    def format_options(self, ctx, formatter):
        options_by_name = {
            param.name: param
            for param in self.get_params(ctx)
            if isinstance(param, click.Option)
        }
        assigned: set[str] = set()

        def _records(names):
            records = []
            for name in names:
                param = options_by_name.get(name)
                if param is None:
                    continue
                record = param.get_help_record(ctx)
                if record is not None:
                    records.append(record)
                    assigned.add(name)
            return records

        for section_title, names in self.option_sections.items():
            records = _records(names)
            if records:
                with formatter.section(section_title):
                    formatter.write_dl(records)

        leftover = [
            record
            for name, param in options_by_name.items()
            if name not in assigned
            and (record := param.get_help_record(ctx)) is not None
        ]
        if leftover:
            with formatter.section("Other options"):
                formatter.write_dl(leftover)


class OrderedGroup(click.Group):
    """A click.Group that lists its commands in a fixed, curated order."""

    COMMAND_ORDER = ["test", "map", "profile", "compare", "utilities"]

    def list_commands(self, ctx):
        commands = list(super().list_commands(ctx))
        ordered = [name for name in self.COMMAND_ORDER if name in commands]
        ordered += [name for name in commands if name not in self.COMMAND_ORDER]
        return ordered


@click.group(
    cls=OrderedGroup,
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog="Source & docs: https://github.com/OlmLab/ZipStrain",
)
@click.version_option(version=__version__, prog_name="zipstrain")
def cli():
    """ZipStrain — fast strain-level metagenomic profiling and comparison.

    A typical run goes: map reads to BAMs, profile them at nucleotide
    resolution, then compare samples by ANI.

    Developed by Parsa Ghadermazi and Matt Olm in the Olm Lab at the
    University of Colorado Boulder.
    """
    pass

@cli.group(short_help="Lower-level helper commands.")
def utilities():
    """The commands in this group are related to various utility functions that mainly prepare input files for profiling and comparison."""
    pass


def _emit_stderr_log(prefix: str, **fields: object) -> None:
    payload = " ".join(f"{key}={value}" for key, value in fields.items())
    click.echo(f"{prefix} {payload}".rstrip(), err=True)
    sys.stderr.flush()


# Companion-CSV emission. Outputs are written as parquet; a matching .csv is
# also written when its estimated size is below this threshold, unless the user
# opts out with --no-csv or forces it with --force-csv.
CSV_SIZE_THRESHOLD_MB = 100
_CSV_BYTES_PER_CELL_ESTIMATE = 16


def _estimated_csv_mb(parquet_path: pathlib.Path) -> float:
    """Cheaply estimate the CSV size of a parquet file from its row/column counts."""
    lazy = pl.scan_parquet(parquet_path)
    n_cols = len(lazy.collect_schema().names())
    n_rows = lazy.select(pl.len()).collect().item()
    return (n_rows * n_cols * _CSV_BYTES_PER_CELL_ESTIMATE) / (1024 * 1024)


def _maybe_write_csv(
    parquet_path: pathlib.Path,
    *,
    no_csv: bool,
    force_csv: bool,
    console: "Console | None" = None,
) -> pathlib.Path | None:
    """Write a ``.csv`` next to ``parquet_path`` unless disabled or too large.

    Returns the CSV path when written, otherwise ``None``. ``--no-csv`` always
    wins; ``--force-csv`` writes regardless of the estimated size.
    """
    parquet_path = pathlib.Path(parquet_path)
    if no_csv or not parquet_path.exists():
        return None
    if not force_csv and _estimated_csv_mb(parquet_path) >= CSV_SIZE_THRESHOLD_MB:
        if console is not None:
            console.print(
                f"[yellow]Skipping CSV for {parquet_path.name}[/] "
                f"(estimated > {CSV_SIZE_THRESHOLD_MB} MB; use --force-csv to write it)."
            )
        return None
    csv_path = parquet_path.with_suffix(".csv")
    pl.read_parquet(parquet_path).write_csv(csv_path)
    return csv_path


def _discover_taxonomy_file(reference_fasta, stb_file, explicit) -> pathlib.Path | None:
    """Locate a genome-taxonomy TSV (explicit, or next to the reference/STB)."""
    if explicit is not None:
        return pathlib.Path(explicit)
    candidates = []
    if reference_fasta is not None:
        ref = pathlib.Path(reference_fasta)
        candidates.append(ref.with_name(mapping.REFERENCE_TAXONOMY_NAME))
    if stb_file is not None:
        candidates.append(pathlib.Path(stb_file).with_name(mapping.REFERENCE_TAXONOMY_NAME))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _discover_genome_lengths(run_dir, profile_locations) -> pathlib.Path | None:
    """Locate a genome_lengths.parquet asset for percent_compared.

    Checks the run's own profiling_assets, then the profiling_assets next to each
    input profile (profiles live at <run>/<sample>/<sample>_profile.parquet, so
    the assets sit two levels up).
    """
    candidates = [pathlib.Path(run_dir) / "profiling_assets" / pf.ASSET_GENOME_LENGTH_FILENAME]
    for loc in profile_locations:
        loc = pathlib.Path(loc)
        candidates.append(loc.parent.parent / "profiling_assets" / pf.ASSET_GENOME_LENGTH_FILENAME)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _add_percent_compared(comparison_parquet, genome_lengths_file) -> None:
    """Add ``percent_compared`` (= total_positions / genome_length) in place.

    A no-op when the genome-length asset is missing or the comparison has no
    ``total_positions`` column. Works on the final parquet regardless of engine
    (standard or matrix), so it never touches the compare compute paths.
    """
    comparison_parquet = pathlib.Path(comparison_parquet)
    if genome_lengths_file is None or not comparison_parquet.exists():
        return
    comp = pl.read_parquet(comparison_parquet)
    if "total_positions" not in comp.columns or "percent_compared" in comp.columns:
        return
    lengths = (
        pl.read_parquet(genome_lengths_file)
        .select(["genome", "genome_length"])
        .unique(subset=["genome"])
    )
    comp = comp.join(lengths, on="genome", how="left").with_columns(
        percent_compared=(
            pl.col("total_positions") / pl.col("genome_length")
        ).cast(pl.Float64)
    ).drop("genome_length")
    comp.write_parquet(comparison_parquet, compression="zstd")


def _add_coverage_overlap(comparison_parquet, profiles) -> None:
    """Add ``coverage_overlap`` in place: compared / positions-covered-in-either.

    Per-sample covered-position counts come from each sample's genome_stats
    (the ``Nx_cov_sites`` column, a sibling of its profile). A no-op when those
    tables or the ``sample_1``/``sample_2`` columns are unavailable. Exact when
    the compare ``--min-cov`` matches the profiling coverage cutoff (both 5 by
    default); otherwise approximate.
    """
    comparison_parquet = pathlib.Path(comparison_parquet)
    if not comparison_parquet.exists():
        return
    comp = pl.read_parquet(comparison_parquet)
    needed = {"sample_1", "sample_2", "genome", "total_positions"}
    if not needed.issubset(comp.columns) or "coverage_overlap" in comp.columns:
        return
    covered_frames = []
    for name, location in profiles:
        gs_path = pathlib.Path(str(location).replace("_profile.parquet", "_genome_stats.parquet"))
        if not gs_path.exists():
            return  # need every sample's stats to be faithful; skip otherwise
        gs = pl.read_parquet(gs_path)
        cov_cols = [c for c in gs.columns if c.endswith("x_cov_sites")]
        if not cov_cols:
            return
        covered_frames.append(
            gs.select(
                pl.lit(name).alias("sample"),
                pl.col("genome"),
                pl.col(cov_cols[0]).cast(pl.Int64).alias("covered"),
            )
        )
    covered = pl.concat(covered_frames)
    comp = (
        comp.join(
            covered.rename({"sample": "sample_1", "covered": "covered_1"}),
            on=["sample_1", "genome"], how="left",
        )
        .join(
            covered.rename({"sample": "sample_2", "covered": "covered_2"}),
            on=["sample_2", "genome"], how="left",
        )
        .with_columns(
            coverage_overlap=(
                pl.col("total_positions")
                / (pl.col("covered_1") + pl.col("covered_2") - pl.col("total_positions"))
            ).cast(pl.Float64)
        )
        .drop(["covered_1", "covered_2"])
    )
    comp.write_parquet(comparison_parquet, compression="zstd")


def _add_presence_column(
    genome_stats: pl.DataFrame,
    *,
    ber: float,
    fug: float,
    min_cov_use_fug: float,
    min_coverage: float,
) -> pl.DataFrame:
    """Add a 'presence' column ('present'/'absent') using the Metapresence logic.

    A genome is called present when both of the following hold:

    - **BER / FUG.** Above ``min_cov_use_fug`` coverage, its breadth-error ratio
      ``ber`` exceeds the threshold; at or below that coverage the ``fug``
      criterion is also required. ``fug`` (fraction of non-excess gaps) is ~0.632
      for uniformly-distributed reads, so present requires ``fug / 0.632 > fug``
      (higher = more uniform = present).
    - **Minimum coverage.** At least ``min_coverage`` mean coverage.
    """
    ber_fug_call = (
        pl.when(pl.col("coverage") > min_cov_use_fug)
        .then(pl.col("ber") > ber)
        .otherwise((pl.col("fug") / 0.632 > fug) & (pl.col("ber") > ber))
        .fill_null(False)
    )
    is_present = ber_fug_call & (pl.col("coverage") >= min_coverage).fill_null(False)
    return genome_stats.with_columns(
        pl.when(is_present).then(pl.lit("present")).otherwise(pl.lit("absent")).alias("presence")
    )


def _finalize_profile_outputs(
    run_dir: pathlib.Path,
    *,
    emit_snvs: bool,
    snv_min_cov: int,
    presence_ber: float,
    presence_fug: float,
    presence_min_cov_use_fug: float,
    presence_min_coverage: float,
    taxonomy_file: pathlib.Path | None = None,
    console=None,
) -> None:
    """Post-process a completed profile run: presence calls and SNV tables."""
    run_dir = pathlib.Path(run_dir)

    taxonomy = None
    if taxonomy_file is not None and pathlib.Path(taxonomy_file).exists():
        taxonomy = pl.read_csv(taxonomy_file, separator="\t").select("genome", "genome_taxonomy")

    # Add presence and (when available) taxonomy columns to each genome_stats table.
    for genome_stats_path in sorted(run_dir.glob("*/*_genome_stats.parquet")):
        stats = pl.read_parquet(genome_stats_path)
        if {"coverage", "ber", "fug"}.issubset(stats.columns):
            stats = _add_presence_column(
                stats, ber=presence_ber, fug=presence_fug, min_cov_use_fug=presence_min_cov_use_fug, min_coverage=presence_min_coverage
            )
        if taxonomy is not None and "genome_taxonomy" not in stats.columns:
            stats = stats.join(taxonomy, on="genome", how="left")
        stats.write_parquet(genome_stats_path)

    # Call SNPs/SNVs relative to the reference for each profile (needs a reference).
    if emit_snvs:
        for profile_path in sorted(run_dir.glob("*/*_profile.parquet")):
            profile_lf = pl.scan_parquet(profile_path)
            if ut.REF_BASE_BITMASK_COLUMN not in profile_lf.collect_schema().names():
                if console is not None:
                    console.print(
                        f"[yellow]Skipping SNV calls for {profile_path.name}[/] "
                        "(needs --reference-fasta so reference bases are recorded)."
                    )
                continue
            snv_path = profile_path.with_name(profile_path.name.replace("_profile.parquet", "_SNVs.parquet"))
            pf.get_reference_snps(profile_lf, min_cov=snv_min_cov).sink_parquet(snv_path, compression="zstd")


def _default_container_address(engine_kind: str) -> str:
    if __version__ in {"", "unknown"}:
        raise ValueError(
            "Cannot infer a default container image because the ZipStrain version is unknown. "
            "Pass --container-address explicitly."
        )
    if engine_kind == "docker":
        return f"{DEFAULT_CONTAINER_REPOSITORY}:{__version__}"
    if engine_kind == "apptainer":
        return f"docker://{DEFAULT_CONTAINER_REPOSITORY}:{__version__}"
    raise ValueError(f"Container defaults are only defined for docker/apptainer, not {engine_kind!r}.")


def _build_container_engine(container_engine: str, container_address: str | None) -> tm.Engine:
    if container_engine == "local":
        return tm.LocalEngine(address="")
    if container_engine == "docker":
        return tm.DockerEngine(address=container_address or _default_container_address("docker"))
    if container_engine == "apptainer":
        return tm.ApptainerEngine(address=container_address or _default_container_address("apptainer"))
    raise ValueError("Invalid container engine. Choose from 'local', 'docker', or 'apptainer'.")


def _load_profile_database(profile_db: str, allow_mismatch: bool = False) -> "db.ProfileDatabase":
    """Load a ProfileDatabase from either a profiles CSV or a pre-built parquet.

    Passing a ``.csv`` (with ``profile_name,profile_location`` columns) builds
    the profile database in memory, so ``zipstrain compare`` can be run without
    first calling ``zipstrain utilities build-profile-db``. Any other extension
    is treated as a pre-built profile-database parquet.
    """
    path = pathlib.Path(profile_db)
    if path.suffix.lower() == ".csv":
        return db.ProfileDatabase.from_csv(path, allow_mismatch=allow_mismatch)
    return db.ProfileDatabase(path)


class _ThrottledMatrixLogger:
    def __init__(self, prefix: str, detail_formatter):
        self.prefix = prefix
        self.detail_formatter = detail_formatter
        self.last_percent_bucket = -1
        self.last_log_time = 0.0
        self.logged_advance = False
        self.start_time = time.monotonic()
        self.last_processing_detail: str | None = None

    def _elapsed_field(self) -> str:
        return f"{time.monotonic() - self.start_time:.1f}s"

    def __call__(self, event: dict[str, object]) -> None:
        phase = str(event.get("phase", ""))
        completed = int(event.get("completed", 0))
        total = int(event.get("total", 0))
        now = time.monotonic()
        detail = event.get("detail")

        if phase in {"start", "done"}:
            self.last_percent_bucket = -1 if total <= 0 else int((completed / max(total, 1)) * 100) // 5
            self.last_log_time = now
            if phase == "start":
                self.logged_advance = False
                self.last_processing_detail = None
            _emit_stderr_log(
                f"{self.prefix} {phase.upper()}",
                elapsed=self._elapsed_field(),
                completed=completed,
                total=total,
                **self.detail_formatter(event),
            )
            return

        if phase == "processing":
            if detail is None:
                return
            detail_str = str(detail)
            should_log = (
                detail_str != self.last_processing_detail
                or (now - self.last_log_time) >= 5.0
            )
            if not should_log:
                return
            self.last_processing_detail = detail_str
            self.last_log_time = now
            _emit_stderr_log(
                f"{self.prefix} PROCESSING",
                elapsed=self._elapsed_field(),
                detail=detail_str,
                completed=completed,
                total=total,
                **self.detail_formatter(event),
            )
            return

        if phase != "advance":
            return

        percent_bucket = int((completed / max(total, 1)) * 100) // 5 if total > 0 else 0
        should_log = (
            (not self.logged_advance)
            or
            total <= 20
            or completed == total
            or percent_bucket > self.last_percent_bucket
            or (now - self.last_log_time) >= 5.0
        )
        if not should_log:
            return

        self.logged_advance = True
        self.last_percent_bucket = percent_bucket
        self.last_log_time = now
        _emit_stderr_log(
            f"{self.prefix} PROGRESS",
            elapsed=self._elapsed_field(),
            completed=completed,
            total=total,
            percent=f"{(completed / max(total, 1)) * 100:.1f}" if total > 0 else "0.0",
            **self.detail_formatter(event),
        )


def _build_null_model_frame(error_rate: float, max_total_reads: int, p_threshold: float, model_type: str) -> pl.DataFrame:
    return pf._build_null_model_frame(error_rate, max_total_reads, p_threshold, model_type)

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
    df_thresh = _build_null_model_frame(error_rate, max_total_reads, p_threshold, model_type)
    df_thresh.write_parquet(output_file)

@utilities.command("merge_parquet")
@click.option('--input-dir', '-i', required=True, help="Directory containing Parquet files to merge.") 
@click.option('--output-file', '-o', required=True, help="Path to save the merged Parquet file.")
@click.option(
    '--batch-size',
    type=int,
    default=-1,
    show_default=True,
    help="Number of parquet files to merge per stage. Use -1 for the current single-pass behavior.",
)
@click.option(
    '--allow-mismatch',
    is_flag=True,
    help="Allow merge even when ZipStrain comparison metadata differ across input parquet files.",
)
def merge_parquet(input_dir, output_file, batch_size, allow_mismatch):
    """
    Merge multiple Parquet files in a directory into a single Parquet file.

    Args:
    input_dir (str): Directory containing Parquet files to merge.
    output_file (str): Path to save the merged Parquet file.
    batch_size (int): Maximum number of parquet files to merge per staged batch.
    """
    def _emit_progress(message: str) -> None:
        click.echo(message, err=True)
        sys.stderr.flush()

    ut.merge_parquet_files(
        input_dir=input_dir,
        output_file=output_file,
        batch_size=batch_size,
        allow_mismatch=allow_mismatch,
        progress_callback=_emit_progress,
    )


@utilities.command("parquet-to-csv")
@click.option('--input-file', '-i', required=True, help="Input parquet file.")
@click.option('--output-file', '-o', default=None, help="Output CSV file. Defaults to the input path with a .csv suffix.")
@click.option('--separator', default=",", show_default=True, help="CSV field separator.")
@click.option('--no-header', is_flag=True, default=False, show_default=True, help="Do not write a header row.")
def parquet_to_csv(input_file, output_file, separator, no_header):
    """Convert a parquet table to CSV."""
    try:
        output_path = ut.parquet_to_csv(
            input_file=input_file,
            output_file=output_file,
            separator=separator,
            include_header=not no_header,
        )
    except Exception as exc:
        raise click.UsageError(str(exc)) from exc
    click.echo(f"wrote={output_path}")


@utilities.command("adjust-sequence-errors")
@click.option('--profile-parquet', '-p', required=True, help="Input profile parquet file.")
@click.option('--null-model', '-n', required=True, help="Null model parquet file.")
@click.option('--output-file', '-o', required=True, help="Path to save the sequence-adjusted profile parquet.")
def adjust_sequence_errors(profile_parquet, null_model, output_file):
    """
    Apply ZipStrain's sequence-error adjustment to an existing profile parquet.

    The output preserves the input column order and drops temporary columns used
    during the adjustment.
    """
    try:
        pf.adjust_profile_parquet_for_sequence_errors(
            profile_parquet=pathlib.Path(profile_parquet),
            null_model_parquet=pathlib.Path(null_model),
            output_file=pathlib.Path(output_file),
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


@utilities.command("sort-profile")
@click.option('--input-profile', '-i', required=True, help="Classic profile parquet to sort in place and mark as sorted.")
@click.option('--tmp-dir', required=False, default=None, help="Temporary working directory used while rewriting the profile.")
def sort_profile(input_profile, tmp_dir):
    """
    Sort a classic profile parquet in place and attach the sortedness metadata marker.
    """
    try:
        pf.sort_profile_parquet_in_place(
            profile_parquet=pathlib.Path(input_profile),
            tmp_dir=pathlib.Path(tmp_dir) if tmp_dir is not None else None,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


@utilities.command("process_mpileup")
@click.option('--batch-size', '-s', default=10000, help="Buffer size for processing stdin from samtools.")
@click.option('--output-file', '-o', required=True, help="Location to save the output Parquet file.")
def process_mpileup(batch_size, output_file):
    """
    Process mpileup files and save the results in a Parquet file.

    Args:
    gene_range_table_loc (str): Path to the gene range table in TSV format.
    batch_bed (str): Path to the batch BED file.
    output_file (str): Path to save the output Parquet file.
    """
    ut.process_mpileup_function(batch_size, output_file)
    
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
    stb = pf.read_stb(stb_file)
    
    bed_table = pl.scan_csv(bed_file, separator='\t',has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end")
    ).select(["scaffold", "start", "end"])
    genome_length = ut.extract_genome_length(stb, bed_table)
    genome_length.sink_parquet(output_file, compression='zstd')

@utilities.command("merge-stat-tables")
@click.option('--stat-table', '-s', multiple=True, required=True, help="Stat parquet file to include. Repeat for multiple files.")
@click.option('--output-file', '-o', required=True, help="Path to save the merged stat table.")
def merge_stat_tables(stat_table, output_file):
    """
    Concatenate stat tables and add a sample column inferred from each file name.
    """
    def _emit_progress(message: str) -> None:
        click.echo(message, err=True)
        sys.stderr.flush()

    ut.merge_stat_tables(
        stat_tables=list(stat_table),
        output_file=output_file,
        progress_callback=_emit_progress,
    )


@utilities.command("get-coverage-stats")
@click.option('--profile-parquet', '-p', required=True, help="Classic profile parquet file.")
@click.option('--gene-bed', '-g', required=True, help="Gene BED/range file. Supports 4 columns (gene, scaffold, start, end) or 5 columns (+ genome).")
@click.option('--genome-bed', '-b', required=True, help="Genome BED file. Supports 3 columns (scaffold, start, end) or 4 columns (+ genome).")
@click.option('--output-dir', '-o', required=True, help="Directory to write the gene/genome stats parquet files.")
@click.option('--prefix', required=True, help="Prefix for the output files.")
def get_coverage_stats(profile_parquet, gene_bed, genome_bed, output_dir, prefix):
    """
    Build coverage-only gene and genome stats from an existing profile parquet.
    """
    summary = ut.get_coverage_stats(
        profile_parquet=profile_parquet,
        gene_bed_file=gene_bed,
        genome_bed_file=genome_bed,
        output_dir=output_dir,
        prefix=prefix,
    )
    click.echo(
        f"gene_stats={summary['gene_stats_file']} "
        f"genome_stats={summary['genome_stats_file']} "
        f"gene_rows={summary['gene_rows']} "
        f"genome_rows={summary['genome_rows']} "
        f"cov_sites_column={summary['cov_sites_column']}"
    )


@utilities.command("generate-genome-pairs")
@click.option('--profile-dir', '-p', required=True, help="Directory containing classic ZipStrain profile parquets.")
@click.option('--output-file', '-o', required=True, help="Output parquet file with all non-redundant profile pairs.")
@click.option('--write-batch-size', type=int, default=100000, show_default=True, help="How many pairs to buffer before writing a parquet row group.")
def generate_genome_pairs(profile_dir, output_file, write_batch_size):
    """Generate a pair table ready for chunk-genome-compare or other compare utilities."""
    def _emit_progress(message: str) -> None:
        click.echo(message, err=True)
        sys.stderr.flush()

    summary = ut.generate_genome_pairs(
        profile_dir=profile_dir,
        output_file=output_file,
        write_batch_size=write_batch_size,
        progress_callback=_emit_progress,
    )
    click.echo(
        f"wrote={summary['output_file']} "
        f"profiles={summary['profiles']} "
        f"pairs={summary['pairs']}"
    )


@utilities.command("chunk-genome-compare")
@click.option('--pair-table', '-p', required=True, help="Pair table parquet/csv/tsv with profile locations.")
@click.option('--stb-file', '-s', required=False, default=None, help="Optional scaffold-to-genome mapping file. When provided, all genomes from the mapping appear in the output; otherwise only genomes with comparable loci are reported.")
@click.option('--output-file', '-o', required=True, help="Path to save the merged genome comparison parquet.")
@click.option('--workers', '-w', type=int, default=None, help="Parallel workers inside this compare chunk. Defaults to CPU count capped by pair count.")
@click.option('--min-cov', '-c', default=5, show_default=True, help="Minimum coverage to consider a position.")
@click.option('--min-gene-compare-len', '-l', default=100, show_default=True, help="Minimum gene length to consider for comparison.")
@click.option('--genome', '-g', default="all", show_default=True, help="Optional genome scope.")
@click.option('--ani-method', '-a', default="popani", show_default=True, help="ANI calculation method to use.")
@click.option('--calculate', default="all", show_default=True, help="Genome metrics to compute: ani, conani, ibs, identical_genes. Combine with '+', or use all.")
@click.option('--engine', type=click.Choice(["polars", "duckdb"]), default="polars", show_default=True, help="Execution engine for compare.")
@click.option('--duckdb-memory-limit', default=None, help="DuckDB memory limit (for scoped filtering or duckdb engine).")
@click.option('--duckdb-temp-directory', default=None, help="Directory DuckDB can use for spill files.")
@click.option('--duckdb-threads', type=int, default=None, help="Number of DuckDB worker threads.")
def chunk_genome_compare(pair_table, stb_file, output_file, workers, min_cov, min_gene_compare_len, genome, ani_method, calculate, engine, duckdb_memory_limit, duckdb_temp_directory, duckdb_threads):
    """Run classic genome compare over one pair-table chunk using Python-side parallel workers."""
    def _emit_progress(message: str) -> None:
        click.echo(message, err=True)
        sys.stderr.flush()

    summary = ut.chunk_genome_compare(
        pair_table=pair_table,
        output_file=output_file,
        stb_file=stb_file,
        workers=workers,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome,
        ani_method=ani_method,
        calculate=calculate,
        engine=engine,
        duckdb_memory_limit=duckdb_memory_limit,
        duckdb_temp_directory=duckdb_temp_directory,
        duckdb_threads=duckdb_threads,
        progress_callback=_emit_progress,
    )
    click.echo(
        f"wrote={summary['output_file']} "
        f"pairs={summary['pairs']} "
        f"rows={summary['rows']} "
        f"workers={summary['workers']} "
        f"elapsed_s={summary['elapsed_seconds']:.2f} "
        f"avg_wall_s_per_pair={summary['avg_wall_seconds_per_pair']:.4f} "
        f"avg_compute_s_per_pair={summary['avg_compute_seconds_per_pair']:.4f} "
        f"avg_s_per_genome_row={summary['avg_seconds_per_genome_row']:.4f} "
        f"avg_genome_rows_per_pair={summary['avg_genome_rows_per_pair']:.2f}"
    )
    
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
    stb = pf.read_stb(stb_file)
    
    het_profile = pf.get_strain_hetrogeneity(profile, stb, min_cov=min_cov, freq_threshold=freq_threshold)
    het_profile.sink_parquet(output_file, compression='zstd')

@utilities.command("get-snp-reference")
@click.option('--profile-file', '-p', required=True, help="Path to the profile Parquet file.")
@click.option('--min-cov', '-c', default=5, show_default=True, help="Minimum coverage required for a site to contribute.")
@click.option('--fmt', default="parquet", show_default=True, type=click.Choice(["parquet", "vcf"]), help="Output format for the SNP-only result.")
@click.option('--output-file', '-o', required=True, help="Path to save the SNP-only output file.")
def get_snp_reference(profile_file, min_cov, fmt, output_file):
    """
    Emit profile-like rows that are SNPs relative to the reference.
    """
    profile = pl.scan_parquet(profile_file)
    try:
        pf.write_reference_snps(
            profile,
            output_file=output_file,
            min_cov=min_cov,
            fmt=fmt,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

@utilities.command("build-profile-db")
@click.option('--profile-db-csv', '-p', required=True, help="Path to the profile database CSV file.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
@click.option('--allow-mismatch', is_flag=True, help="Skip profile contract metadata validation when building the profile database.")
def build_profile_db(profile_db_csv, output_file, allow_mismatch):
    """
    Build a profile database from the given CSV file.

    Args:
    profile_db_csv (str): Path to the profile database CSV file.
    """
    profile_db = db.ProfileDatabase.from_csv(pathlib.Path(profile_db_csv), allow_mismatch=allow_mismatch)
    profile_db.save_as_new_database(pathlib.Path(output_file))


@utilities.command("build-matrix-db")
@click.option('--profile-dir', '-p', required=True, help="Directory containing classic ZipStrain profile parquets.")
@click.option('--output-file', '-o', required=True, help="Output matrix store.")
@click.option('--genome', '-g', default="all", show_default=True, help="Optional genome scope.")
@click.option('--bed-file', '-b', required=True, help="BED file defining scaffold extents for the matrix contract.")
@click.option('--stb-file', '-s', required=True, help="STB file defining scaffold-to-genome mapping for the matrix contract.")
@click.option('--gene-range-table', default=None, help="Optional headerless TSV of gene, scaffold, start, end to store gene-coordinate ranges for gene ANI.")
@click.option(
    '--count-dtype',
    type=click.Choice(sorted(mp.COUNT_DTYPES.keys())),
    default="uint16",
    show_default=True,
    help="Stored dtype for whole-genome dense matrices.",
)
@click.option(
    '--memory-limit-gb',
    type=float,
    default=16.0,
    show_default=True,
    help="Approximate maximum memory budget for the entire matrix build process.",
)
@click.option(
    '--export-batch-mb',
    type=float,
    default=mp.MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
    show_default=True,
    help="Approximate matrix-store sample-axis chunk target size, in MiB.",
)
@click.option(
    '--sparse',
    is_flag=True,
    default=False,
    help="Store genome matrices sparsely in HDF5. Compare currently materializes them back to dense on load.",
)
def build_matrix_db(profile_dir, output_file, genome, bed_file, stb_file, gene_range_table, count_dtype, memory_limit_gb, export_batch_mb, sparse):
    """
    Build a matrix store directly from classic profile parquets.

    Each genome is stored as one sample-major matrix with shape
    samples x positions x 4, while preserving the same whole-genome axis and
    separator-row rules used by matrix compare.
    """
    progress_console = Console(stderr=True)
    use_progress_bar = sys.stderr.isatty()

    if use_progress_bar:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("stored={task.fields[stored_rows]}"),
            TimeElapsedColumn(),
            console=progress_console,
            transient=True,
        ) as progress:
            task_id = progress.add_task(
                "Building matrix store",
                total=1,
                completed=0,
                stored_rows="0",
            )

            def _progress_callback(event: dict[str, object]) -> None:
                total = int(event.get("total", 0)) or 1
                completed = int(event.get("completed", 0))
                stored_rows = str(event.get("stored_rows", 0))
                progress.update(
                    task_id,
                    total=total,
                    completed=completed,
                    stored_rows=stored_rows,
                )

            summary = mp.build_matrix_hdf5(
                profile_dir=pathlib.Path(profile_dir),
                output_file=pathlib.Path(output_file),
                genome=genome,
                bed_file=pathlib.Path(bed_file),
                stb_file=pathlib.Path(stb_file),
                gene_range_table=pathlib.Path(gene_range_table) if gene_range_table is not None else None,
                count_dtype=count_dtype,
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
                sparse=sparse,
                progress_callback=_progress_callback,
            )
    else:
        progress_logger = _ThrottledMatrixLogger(
            "MATRIX-BUILD",
            lambda event: {"stored_rows": event.get("stored_rows", 0)},
        )
        summary = mp.build_matrix_hdf5(
            profile_dir=pathlib.Path(profile_dir),
            output_file=pathlib.Path(output_file),
            genome=genome,
            bed_file=pathlib.Path(bed_file),
            stb_file=pathlib.Path(stb_file),
            gene_range_table=pathlib.Path(gene_range_table) if gene_range_table is not None else None,
            count_dtype=count_dtype,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
            sparse=sparse,
            progress_callback=progress_logger,
        )
    click.echo(
        f"wrote={summary.output_file} "
        f"samples={summary.sample_count} "
        f"scaffolds={summary.scaffold_count} "
        f"stored_rows={summary.stored_rows}"
    )


@utilities.command("append-matrix-db")
@click.option('--profile-dir', '-p', required=True, help="Directory containing new classic ZipStrain profile parquets to append.")
@click.option('--matrix-db-file', '-m', required=True, help="Existing matrix store to append to.")
@click.option(
    '--memory-limit-gb',
    type=float,
    default=16.0,
    show_default=True,
    help="Approximate maximum memory budget for the append process.",
)
@click.option(
    '--export-batch-mb',
    type=float,
    default=mp.MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
    show_default=True,
    help="Approximate matrix-store sample-axis chunk target size, in MiB, used when rewriting the store.",
)
def append_matrix_db(profile_dir, matrix_db_file, memory_limit_gb, export_batch_mb):
    """
    Append new profiles to an existing matrix store.

    The append uses the existing genome/scaffold layout already stored in the
    matrix store as the contract. Contract-valid genomes are appended or
    materialized as needed. Genomes outside the stored contract are ignored.
    Known genomes with incompatible scaffolds or coordinate ranges still fail.
    Sample names must also be new.
    """
    progress_console = Console(stderr=True)
    use_progress_bar = sys.stderr.isatty()

    if use_progress_bar:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("stored={task.fields[stored_rows]}"),
            TimeElapsedColumn(),
            console=progress_console,
            transient=True,
        ) as progress:
            task_id = progress.add_task(
                "Appending matrix store",
                total=1,
                completed=0,
                stored_rows="0",
            )

            def _progress_callback(event: dict[str, object]) -> None:
                total = int(event.get("total", 0)) or 1
                completed = int(event.get("completed", 0))
                stored_rows = str(event.get("stored_rows", 0))
                progress.update(
                    task_id,
                    total=total,
                    completed=completed,
                    stored_rows=stored_rows,
                )

            summary = mp.append_matrix_hdf5(
                profile_dir=pathlib.Path(profile_dir),
                matrix_hdf5_file=pathlib.Path(matrix_db_file),
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
                progress_callback=_progress_callback,
            )
    else:
        progress_logger = _ThrottledMatrixLogger(
            "MATRIX-APPEND",
            lambda event: {"stored_rows": event.get("stored_rows", 0)},
        )
        summary = mp.append_matrix_hdf5(
            profile_dir=pathlib.Path(profile_dir),
            matrix_hdf5_file=pathlib.Path(matrix_db_file),
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
            progress_callback=progress_logger,
        )
    click.echo(
        f"wrote={summary.output_file} "
        f"appended_samples={summary.appended_sample_count} "
        f"total_samples={summary.total_sample_count} "
        f"scaffolds={summary.scaffold_count} "
        f"stored_rows={summary.stored_rows} "
        f"ignored_genomes={summary.ignored_genome_count}"
    )


@utilities.command("matrix-compare")
@click.option('--matrix-db-file', '-m', required=True, help="Input matrix store from build-matrix-db or matrix-db-to-hdf5.")
@click.option('--output-file', '-o', required=True, help="Output DuckDB compare database. If it already exists, only remaining pairs are added.")
@click.option('--genome', '-g', default="all", show_default=True, help="Optional genome scope.")
@click.option('--memory-limit-gb', type=float, default=16.0, show_default=True, help="Approximate memory budget for compare.")
@click.option('--anchor-queue-size', type=int, default=1, show_default=True, help="Host-side torch anchor queue size. Only one anchor is transferred to the GPU at a time.")
@click.option('--target-queue-size', type=int, default=1, show_default=True, help="Host-side torch target queue size. `1` keeps the current synchronous target-load behavior.")
@click.option('--result-transfer-batch-size', type=int, default=mp.MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE_DEFAULT, show_default=True, help="Number of torch compare units to batch before transferring result vectors back to CPU.")
@click.option(
    '--loader-executor',
    type=click.Choice(mp.MATRIX_IO_EXECUTOR_KINDS),
    default="thread",
    show_default=True,
    help="Executor used for torch loader prefetch work.",
)
@click.option(
    '--writer-executor',
    type=click.Choice(mp.MATRIX_IO_EXECUTOR_KINDS),
    default="thread",
    show_default=True,
    help="Executor used for torch result writing/checkpoint work.",
)
@click.option('--calculate', default="all", show_default=True, help="Matrix metrics to compute: ani, ibs, gene. Combine with '+'. 'all' means ani+ibs, and also gene when the matrix store has gene annotations.")
@click.option(
    '--backend',
    type=click.Choice(mp.MATRIX_PAIR_BACKENDS),
    default="numpy",
    show_default=True,
    help="Compute backend. Torch backends use CPU/CUDA/MPS depending on selection.",
)
def matrix_compare(
    matrix_db_file,
    output_file,
    genome,
    memory_limit_gb,
    anchor_queue_size,
    target_queue_size,
    result_transfer_batch_size,
    loader_executor,
    writer_executor,
    calculate,
    backend,
):
    """
    Run resumable matrix compare on all non-redundant, non-self sample pairs.

    This command writes results into a DuckDB compare database. If that compare
    DB already exists, only pairs that are not yet marked completed are
    processed and appended.
    """
    progress_console = Console(stderr=True)
    use_progress_bar = sys.stderr.isatty()

    if use_progress_bar:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("chunks={task.fields[target_chunks]}"),
            TimeElapsedColumn(),
            console=progress_console,
            transient=True,
        ) as progress:
            task_id = progress.add_task(
                "Comparing matrix store",
                total=1,
                completed=0,
                target_chunks="0",
            )

            def _progress_callback(event: dict[str, object]) -> None:
                total = int(event.get("total", 0)) or 1
                completed = int(event.get("completed", 0))
                target_chunks = str(event.get("target_chunks", 0))
                description = "Comparing matrix store"
                if event.get("phase") == "processing" and event.get("detail"):
                    description = f"Comparing matrix store: {event['detail']}"
                progress.update(
                    task_id,
                    description=description,
                    total=total,
                    completed=completed,
                    target_chunks=target_chunks,
                )

            summary = mp.matrix_compare(
                matrix_db_file=pathlib.Path(matrix_db_file),
                output_file=pathlib.Path(output_file),
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome=genome,
                memory_limit_gb=memory_limit_gb,
                anchor_queue_size=anchor_queue_size,
                target_queue_size=target_queue_size,
                result_transfer_batch_size=result_transfer_batch_size,
                loader_executor_kind=loader_executor,
                writer_executor_kind=writer_executor,
                backend=backend,
                calculate=calculate,
                emit_writer_logs=False,
                progress_callback=_progress_callback,
            )
    else:
        summary = mp.matrix_compare(
            matrix_db_file=pathlib.Path(matrix_db_file),
            output_file=pathlib.Path(output_file),
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            genome=genome,
            memory_limit_gb=memory_limit_gb,
            anchor_queue_size=anchor_queue_size,
            target_queue_size=target_queue_size,
            result_transfer_batch_size=result_transfer_batch_size,
            loader_executor_kind=loader_executor,
            writer_executor_kind=writer_executor,
            backend=backend,
            calculate=calculate,
            emit_writer_logs=True,
            progress_callback=None,
        )
    click.echo(
        f"wrote={summary.output_file} "
        f"requested_pairs={summary.requested_pairs} "
        f"rows={summary.written_rows} "
        f"scaffolds={summary.scaffold_count} "
        f"anchor_groups={summary.anchor_groups} "
        f"target_chunks={summary.target_chunks}"
    )


@utilities.command("matrix-db-to-hdf5")
@click.option('--matrix-db-file', '-m', required=True, help="Input legacy DuckDB matrix database from older matrix-store builds.")
@click.option('--output-file', '-o', required=False, default=None, help="Output HDF5 matrix file. Defaults to the same path with a .h5 suffix.")
@click.option(
    '--export-batch-mb',
    type=float,
    default=mp.MATRIX_HDF5_EXPORT_TARGET_BATCH_MB_DEFAULT,
    show_default=True,
    help="Approximate target batch size, in MiB, for streaming sample matrices from DuckDB into HDF5.",
)
@click.option(
    '--sparse',
    is_flag=True,
    default=False,
    help="Store genome matrices sparsely in HDF5. Compare currently materializes them back to dense on load.",
)
def matrix_db_to_hdf5(matrix_db_file, output_file, export_batch_mb, sparse):
    """
    Convert a legacy DuckDB matrix database into the current matrix-store format.
    """
    matrix_db_path = pathlib.Path(matrix_db_file)
    output_path = pathlib.Path(output_file) if output_file is not None else matrix_db_path.with_suffix(".h5")
    progress_console = Console(stderr=True)
    use_progress_bar = sys.stderr.isatty()

    if use_progress_bar:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("stored={task.fields[stored_rows]}"),
            TimeElapsedColumn(),
            console=progress_console,
            transient=True,
        ) as progress:
            task_id = progress.add_task(
                "Converting matrix DB to HDF5",
                total=1,
                completed=0,
                stored_rows="0",
            )

            def _progress_callback(event: dict[str, object]) -> None:
                total = int(event.get("total", 0)) or 1
                completed = int(event.get("completed", 0))
                stored_rows = str(event.get("stored_rows", 0))
                description = "Converting matrix DB to HDF5"
                if event.get("phase") == "processing" and event.get("genome"):
                    description = f"Converting matrix DB to HDF5: {event['genome']}"
                progress.update(
                    task_id,
                    description=description,
                    total=total,
                    completed=completed,
                    stored_rows=stored_rows,
                )

            written = mp.export_matrix_db_hdf5(
                matrix_db_file=matrix_db_path,
                output_file=output_path,
                progress_callback=_progress_callback,
                export_batch_mb=export_batch_mb,
                sparse=sparse,
            )
    else:
        progress_logger = _ThrottledMatrixLogger(
            "MATRIX-HDF5",
            lambda event: {"stored_rows": event.get("stored_rows", 0), "genome": event.get("genome", "")},
        )
        written = mp.export_matrix_db_hdf5(
            matrix_db_file=matrix_db_path,
            output_file=output_path,
            progress_callback=progress_logger,
            export_batch_mb=export_batch_mb,
            sparse=sparse,
        )
    click.echo(f"wrote={written}")


@utilities.command("matrix-compare-export")
@click.option('--matrix-compare-db-file', '-m', required=True, help="Input DuckDB matrix compare database from matrix-compare.")
@click.option('--output-file', '-o', required=True, help="Output parquet file.")
@click.option(
    '--table',
    'result_table',
    type=click.Choice(["genome", "gene"]),
    default="genome",
    show_default=True,
    help="Which compare results table to export.",
)
def matrix_compare_export(matrix_compare_db_file, output_file, result_table):
    """
    Export a matrix compare DuckDB database to parquet.
    """
    try:
        exported = mp.export_matrix_compare_parquet(
            matrix_compare_db_file=pathlib.Path(matrix_compare_db_file),
            output_file=pathlib.Path(output_file),
            table=result_table,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"wrote={exported}")


@utilities.command("build-genome-db")
@click.option('--tool', '-t', required=True, type=click.Choice(sorted(bdb.ADAPTERS.keys())), help="Abundance tool format to parse (for example, sylph).")
@click.option('--abundance-table', '-a', required=True, help="Path to abundance table (csv/tsv/parquet).")
@click.option('--cache-dir', '-c', required=True, help="Genome cache directory. Reuses already-downloaded genomes.")
@click.option('--output-dir', '-o', required=True, help="Directory where concatenated reference FASTA and STB are written.")
@click.option('--download-retries', type=int, default=8, show_default=True, help="Max download attempts per genome before skipping.")
@click.option('--retry-backoff-seconds', type=float, default=10.0, show_default=True, help="Base seconds for exponential retry backoff.")
@click.option('--download-workers', type=int, default=1, show_default=True, help="Parallel genome download workers.")
def build_genome_db(tool, abundance_table, cache_dir, output_dir, download_retries, retry_backoff_seconds, download_workers):
    """
    Build a reference bundle from an abundance table.

    This command extracts genomes with non-zero abundance in at least one sample,
    reuses/downloads them via a local cache directory, and writes:
      - reference_genomes.fna
      - reference_genomes.stb
    into the output directory.
    """
    if download_retries < 1:
        raise ValueError("--download-retries must be >= 1")
    if retry_backoff_seconds < 0:
        raise ValueError("--retry-backoff-seconds must be >= 0")
    if download_workers < 1:
        raise ValueError("--download-workers must be >= 1")

    console = Console()
    progress_callback = None
    if console.is_terminal:
        progress = Progress(
            TextColumn("[bold cyan]Genome Cache"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total} processed"),
            TextColumn("remaining: {task.fields[remaining]}"),
            TextColumn("{task.fields[last_event]}"),
            TimeElapsedColumn(),
            console=console,
            transient=True,
        )
        progress_state = {"task_id": None}

        def _progress_callback(event: str, fetched: int, remaining: int, total: int, accession: str | None) -> None:
            total_display = max(total, 1)
            if event == "downloaded":
                last_event = f"downloaded {accession}"
            elif event == "failed":
                last_event = f"failed {accession}"
            elif event == "already_present":
                last_event = f"cached {accession}"
            elif event == "retry":
                last_event = f"retrying {accession}"
            elif event == "done":
                last_event = "downloads complete"
            elif event == "assembling_reference":
                last_event = "building concatenated reference"
            elif event == "completed":
                last_event = "done"
            else:
                last_event = "starting"

            if progress_state["task_id"] is None:
                progress_state["task_id"] = progress.add_task(
                    "genomes",
                    total=total_display,
                    completed=min(fetched, total_display),
                    remaining=max(remaining, 0),
                    last_event=last_event,
                )
            else:
                progress.update(
                    progress_state["task_id"],
                    total=total_display,
                    completed=min(fetched, total_display),
                    remaining=max(remaining, 0),
                    last_event=last_event,
                )

        progress_callback = _progress_callback
        with progress:
            out_fasta, out_stb, extracted, report, summary = bdb.build_reference_from_abundance(
                tool_name=tool,
                abundance_table=pathlib.Path(abundance_table),
                cache_dir=pathlib.Path(cache_dir),
                output_dir=pathlib.Path(output_dir),
                progress_callback=progress_callback,
                max_download_attempts=download_retries,
                backoff_base_seconds=retry_backoff_seconds,
                download_workers=download_workers,
            )
    else:
        out_fasta, out_stb, extracted, report, summary = bdb.build_reference_from_abundance(
            tool_name=tool,
            abundance_table=pathlib.Path(abundance_table),
            cache_dir=pathlib.Path(cache_dir),
            output_dir=pathlib.Path(output_dir),
            max_download_attempts=download_retries,
            backoff_base_seconds=retry_backoff_seconds,
            download_workers=download_workers,
        )
    failed_rows = report.filter(pl.col("status") == "failed") if report.height else report
    rate_limited_failures = (
        failed_rows.filter(pl.col("error").fill_null("").str.contains("Too Many Requests", literal=True)).height
        if failed_rows.height
        else 0
    )
    report_file = pathlib.Path(output_dir) / "genome_db_build_report.txt"
    report_lines = [
        "Genome DB Build Report",
        f"Tool: {tool}",
        f"Selected genomes (non-zero): {summary['selected_genomes']}",
        f"Cached before run: {summary['cached_before_download']}",
        f"Download attempts (new): {summary['attempted_downloads']}",
        f"Downloaded now: {summary['downloaded_now']}",
        f"Failed downloads: {summary['failed_downloads']}",
        f"Skipped after retries: {summary['missing_after_retries']}",
        f"Retry attempts/genome: {download_retries}",
        f"Retry backoff base (s): {retry_backoff_seconds}",
        f"Download workers: {download_workers}",
        f"Available in cache after run: {summary['cached_after_download']}",
        f"Concatenated FASTA: {out_fasta}",
        f"STB file: {out_stb}",
        f"Cache directory: {pathlib.Path(cache_dir)}",
        "",
        f"Failed IDs ({failed_rows.height}):",
    ]
    if failed_rows.height == 0:
        report_lines.append("- none")
    else:
        for row in failed_rows.select("accession", "error").iter_rows(named=True):
            error = row["error"] if row["error"] not in (None, "") else "no error message"
            report_lines.append(f"- {row['accession']}: {error}")
    if rate_limited_failures > 0:
        report_lines.extend(
            [
                "",
                "Rate-limit hint:",
                "- Reduce --download-workers (for example, 1-2)",
                "- Increase --download-retries (for example, 8)",
                "- Increase --retry-backoff-seconds (for example, 10-20)",
            ]
        )
    report_file.write_text("\n".join(report_lines) + "\n")
    report_table = Table(box=box.SIMPLE_HEAVY, show_header=False, pad_edge=False)
    report_table.add_column("metric", style="bold cyan")
    report_table.add_column("value", style="white")
    report_table.add_row("Tool", tool)
    report_table.add_row("Selected genomes (non-zero)", str(summary["selected_genomes"]))
    report_table.add_row("Cached before run", str(summary["cached_before_download"]))
    report_table.add_row("Download attempts (new)", str(summary["attempted_downloads"]))
    report_table.add_row("Downloaded now", str(summary["downloaded_now"]))
    report_table.add_row("Failed downloads", str(summary["failed_downloads"]))
    report_table.add_row("Skipped after retries", str(summary["missing_after_retries"]))
    report_table.add_row("Retry attempts/genome", str(download_retries))
    report_table.add_row("Retry backoff base (s)", str(retry_backoff_seconds))
    report_table.add_row("Download workers", str(download_workers))
    report_table.add_row("Available in cache after run", str(summary["cached_after_download"]))
    if rate_limited_failures > 0:
        report_table.add_row("Rate-limited failures", str(rate_limited_failures))
    report_table.add_row("Concatenated FASTA", str(out_fasta))
    report_table.add_row("STB file", str(out_stb))
    report_table.add_row("Cache directory", str(pathlib.Path(cache_dir)))
    report_table.add_row("Report file", str(report_file))
    Console().print(Panel(report_table, title="Genome DB Build Report", border_style="green"))


@utilities.command("to-complete-table")
@click.option("--profile-db", "-p", required=True, help="Path to the profile database parquet file.")
@click.option("--comp-db-file", "-c", required=False, default=None, help="Optional current genome comparison parquet.")
@click.option("--output-file", "-o", required=True, help="Path to save the completed pairs csv file.")
def to_complete_table(profile_db, comp_db_file, output_file):
    """
    Generate the not-yet-completed genome-comparison sample pairs and save them to a csv file.

    Parameters:
    profile_db (str): Path to the profile database parquet file.
    output_file (str): Path to save the remaining-pairs CSV file.
    """
    genome_comp_db=db.GenomeComparisonDatabase(
        profile_db=db.ProfileDatabase(pathlib.Path(profile_db)),
        config=db.GenomeComparisonConfig(
            scope="all",
            min_cov=5,
            min_gene_compare_len=100,
            stb_file_loc=None,
        ),
        comp_db_loc=comp_db_file,
    )
    completed_pairs=genome_comp_db.to_complete_input_table()
    completed_pairs.sink_csv(pathlib.Path(output_file), engine="streaming")

@utilities.command("presence-profile")
@click.option('--profile-file', '-p', required=True, help="Path to the profile Parquet file.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--bed-file', '-b', required=True, help="Path to the BED file.")
@click.option('--read-loc-file', '-r', required=True, help="Path to the read location table.")
@click.option('--min-cov-fug', '-c', default=0.1, help="Minimum coverage to use fug.")
@click.option('--fug-threshold', '-f', default=2, help="FUG threshold.")
@click.option('--ber', '-e', default=0.5, help="Minimum ratio of breadth over expected breadth to consider presence.")
@click.option('--output-file', '-o', required=True, help="Path to save the output Parquet file.")
def presence_profile(profile_file, stb_file, bed_file, read_loc_file, min_cov_fug, fug_threshold, ber, output_file):
    """
    Generate a presence profile for genomes based on the given profile and read location data.

    Args:
    profile_file (str): Path to the profile Parquet file.
    stb_file (str): Path to the scaffold-to-genome mapping file.
    bed_file (str): Path to the BED file.
    read_loc_file (str): Path to the read location table.
    min_cov_fug (float): Minimum coverage to use fug.
    fug_threshold (float): FUG threshold.
    ber (float): Minimum ratio of breadth over expected breadth to consider presence.
    output_file (str): Path to save the output Parquet file.
    """
    profile = pl.scan_parquet(profile_file)
    stb = pf.read_stb(stb_file)
    bed = pl.scan_csv(bed_file, separator="\t", has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end")
    ).select(["scaffold", "start", "end"])
    read_loc_table = pl.scan_parquet(read_loc_file).rename({
        "chrom":"scaffold",
        "pos":"loc"
    })
    presence_df = ut.get_genome_stats(
        profile=profile,
        stb=stb,
        bed=bed,
        read_loc_table=read_loc_table,
        min_cov_use_fug=min_cov_fug,
        fug=fug_threshold,
        ber=ber
    )
    presence_df.sink_parquet(output_file, compression='zstd')

@utilities.command("process-read-locs")
@click.option("--output-file", "-o", required=True, help="Path to save the processed read locations Parquet file.")
def process_read_locs(output_file):
    """
    Process read locations and save them to a Parquet file.

    Args:
    output_file (str): Path to save the output Parquet file.
    """
    ut.process_read_location(output_file=pathlib.Path(output_file))

@utilities.command("generate_stb")
@click.option('--genomes-dir-file', '-g', required=True, help="Path to the genomes directory file. A text file with each line containing a genome fasta file path.")
@click.option('--output-file', '-o', required=True, help="Path to save the output scaffold-to-genome mapping file.")
@click.option('--extension', '-e', default=".fasta", help="File extension of the genome fasta files.")
def generate_stb(genomes_dir_file, output_file, extension):
    """
    Generate a scaffold-to-genome mapping file from the given genomes directory file.

    Args:
    genomes_dir_file (str): Path to the genomes directory file.
    output_file (str): Path to save the output scaffold-to-genome mapping file.
    extension (str): File extension of the genome fasta files.
    """
    with open(output_file, 'w') as out_f:
        for genome in pathlib.Path(genomes_dir_file).glob(f"*{extension}"):
            genome_name = genome.stem
            with open(genome, 'r') as gf:
                for line in gf:
                    if line.startswith('>'):
                        scaffold_name = line[1:].strip().split()[0]
                        out_f.write(f"{scaffold_name}\t{genome_name}\n")
    
        
    


@utilities.command("gene-range-table")
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
    gene_locs.write_csv(pathlib.Path(output_file), separator="\t", include_header=False)


@utilities.command("single_compare_genome")
@click.option('--mpileup-contig-1', '-m1', required=True, help="Path to the first mpileup file.")
@click.option('--mpileup-contig-2', '-m2', required=True, help="Path to the second mpileup file.")
@click.option('--stb-file', '-s', required=False, default=None, help="Optional scaffold-to-genome mapping file. When provided, all genomes from the mapping appear in the output; otherwise only genomes with comparable loci are reported.")
@click.option('--min-cov', '-c', default=5, help="Minimum coverage to consider a position.")
@click.option('--min-gene-compare-len', '-l', default=100, help="Minimum gene length to consider for comparison.")
@click.option('--output-file', '-o', required=True, help="Path to save the parquet file.")
@click.option('--genome', '-g', default="all", help="If provided, do the comparison only for the specified genome.")
@click.option('--ani-method', '-a', default="popani", help="ANI calculation method to use (e.g., 'popani', 'conani', 'cosani_0.4').")
@click.option('--calculate', default="all", show_default=True, help="Genome metrics to compute: ani, conani, ibs, identical_genes. Combine with '+', or use all.")
@click.option('--engine', type=click.Choice(["polars", "duckdb"]), default="polars", show_default=True, help="Execution engine for compare.")
@click.option('--duckdb-memory-limit', default=None, help="DuckDB memory limit (e.g., 2GB, 1024MB).")
@click.option('--duckdb-temp-directory', default=None, help="Directory DuckDB can use for spill files.")
@click.option('--duckdb-threads', type=int, default=None, help="Number of DuckDB worker threads.")
def single_compare_genome(mpileup_contig_1, mpileup_contig_2, stb_file, min_cov, min_gene_compare_len, output_file, genome, ani_method, calculate, engine, duckdb_memory_limit, duckdb_temp_directory, duckdb_threads):
    """
    Main function to compare two mpileup files and calculate genome and gene statistics.
    
    Args:
    mpileup_contig_1 (str): Path to the first mpileup file.
    mpileup_contig_2 (str): Path to the second mpileup file.
    min_cov (int): Minimum coverage to consider a position.
    min_gene_compare_len (int): Minimum gene length to consider for comparison.
    output_file (str): Path to save the parquet file.
    genome (str): If provided, do the comparison only for the specified genome.
    stb_file (str): Optional path to the scaffold to genome mapping file.
    """
    mpile_contig_1_name = ut.infer_sample_name_from_profile(mpileup_contig_1)
    mpile_contig_2_name = ut.infer_sample_name_from_profile(mpileup_contig_2)
    if duckdb_threads is not None and duckdb_threads < 1:
        raise ValueError("--duckdb-threads must be >= 1")
    calculations = cp.parse_genome_calculations(calculate)
    output_cols = cp.genome_metric_output_columns(calculations)
    compare_metadata = ut.build_single_compare_metadata(
        mpileup_contig_1,
        mpileup_contig_2,
        compare_kind="genome",
        scope=genome,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        engine=engine,
        uses_stb=stb_file is not None,
    )

    mpile_1_for_compare = mpileup_contig_1
    mpile_2_for_compare = mpileup_contig_2
    if engine == "polars" and genome != "all":
        mpile_1_for_compare, mpile_2_for_compare = cp.polars_prefilter_by_scope(
            mpile1=mpileup_contig_1,
            mpile2=mpileup_contig_2,
            genome_scope=genome,
        )

    if engine == "duckdb":
        cp.duckdb_compare_genomes_to_parquet(
            mpile1=mpileup_contig_1,
            mpile2=mpileup_contig_2,
            output_file=output_file,
            stb_file=stb_file,
            sample_1_name=mpile_contig_1_name,
            sample_2_name=mpile_contig_2_name,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome,
            ani_method=ani_method,
            calculate=calculations,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )
        ut.rewrite_parquet_with_metadata(output_file, compare_metadata)
        return

    comp = cp.compare_genomes(
        mpile_contig_1=mpile_1_for_compare,
        mpile_contig_2=mpile_2_for_compare,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome,
        ani_method=ani_method,
        duckdb_memory_limit=duckdb_memory_limit,
        duckdb_temp_directory=duckdb_temp_directory,
        duckdb_threads=duckdb_threads,
        engine="polars",
        stb_file=stb_file,
        calculate=calculations,
    ).with_columns(
        sample_1=pl.lit(mpile_contig_1_name),
        sample_2=pl.lit(mpile_contig_2_name),
    ).select(output_cols + ["sample_1", "sample_2"])
    comp.sink_parquet(output_file, compression='zstd', metadata=compare_metadata)

@utilities.command("single_compare_gene")
@click.option('--mpileup-contig-1', '-m1', required=True, help="Path to the first mpileup file.")
@click.option('--mpileup-contig-2', '-m2', required=True, help="Path to the second mpileup file.")
@click.option('--stb-file', '-s', required=False, default=None, help="Optional scaffold-to-genome mapping file. Currently unused for gene compare, and accepted for workflow consistency.")
@click.option('--min-cov', '-c', default=5, help="Minimum coverage to consider a position.")
@click.option('--min-gene-compare-len', '-l', default=100, help="Minimum gene length to consider for comparison.")
@click.option('--output-file', '-o', required=True, help="Path to save the parquet file.")
@click.option('--scope', '-g', default="all:all", help="If provided, do the comparison only for the specified genome-gene pair.")
@click.option('--ani-method', '-a', default="popani", help="ANI calculation method to use (e.g., 'popani', 'conani', 'cosani_0.4').")
@click.option('--engine', type=click.Choice(["polars", "duckdb"]), default="polars", show_default=True, help="Execution engine for compare.")
@click.option('--duckdb-memory-limit', default=None, help="DuckDB memory limit (e.g., 2GB, 1024MB).")
@click.option('--duckdb-temp-directory', default=None, help="Directory DuckDB can use for spill files.")
@click.option('--duckdb-threads', type=int, default=None, help="Number of DuckDB worker threads.")
def single_compare_gene(mpileup_contig_1, mpileup_contig_2, stb_file, min_cov, min_gene_compare_len, output_file, scope, ani_method, engine, duckdb_memory_limit, duckdb_temp_directory, duckdb_threads):
    """
    Compare two mpileup files and calculate gene-level comparison statistics.
    
    Args:
    mpileup_contig_1 (str): Path to the first mpileup file.
    mpileup_contig_2 (str): Path to the second mpileup file.
    stb_file (str | None): Optional path to the scaffold to genome mapping file.
    min_cov (int): Minimum coverage to consider a position.
    min_gene_compare_len (int): Minimum gene length to consider for comparison.
    output_file (str): Path to save the parquet file.
    scope (str): If provided, do the comparison only for the specified genome-gene pair.
    ani_method (str): ANI calculation method to use.
    """

    mpile_contig_1_name = ut.infer_sample_name_from_profile(mpileup_contig_1)
    mpile_contig_2_name = ut.infer_sample_name_from_profile(mpileup_contig_2)
    if duckdb_threads is not None and duckdb_threads < 1:
        raise ValueError("--duckdb-threads must be >= 1")

    _ = stb_file
    if ":" not in scope:
        raise ValueError("scope must be in 'GENOME:GENE' format (e.g., all:all).")
    genome_scope, gene_scope = scope.split(":", 1)
    compare_metadata = ut.build_single_compare_metadata(
        mpileup_contig_1,
        mpileup_contig_2,
        compare_kind="gene",
        scope=scope,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        engine=engine,
        uses_stb=stb_file is not None,
    )

    mpile_1_for_compare = mpileup_contig_1
    mpile_2_for_compare = mpileup_contig_2
    if engine == "polars" and (genome_scope != "all" or gene_scope != "all"):
        mpile_1_for_compare, mpile_2_for_compare = cp.polars_prefilter_by_scope(
            mpile1=mpileup_contig_1,
            mpile2=mpileup_contig_2,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
        )

    if engine == "duckdb":
        cp.duckdb_compare_genes_to_parquet(
            mpile1=mpileup_contig_1,
            mpile2=mpileup_contig_2,
            output_file=output_file,
            sample_1_name=mpile_contig_1_name,
            sample_2_name=mpile_contig_2_name,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
            ani_method=ani_method,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )
        ut.rewrite_parquet_with_metadata(output_file, compare_metadata)
        return

    gene_comp = cp.compare_genes(
        mpile_contig_1=mpile_1_for_compare,
        mpile_contig_2=mpile_2_for_compare,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
        ani_method=ani_method,
        duckdb_memory_limit=duckdb_memory_limit,
        duckdb_temp_directory=duckdb_temp_directory,
        duckdb_threads=duckdb_threads,
        engine="polars",
    ).with_columns(
        sample_1=pl.lit(mpile_contig_1_name),
        sample_2=pl.lit(mpile_contig_2_name),
    ).select(
        "genome",
        "gene",
        "total_positions",
        "share_allele_pos",
        "ani",
        "sample_1",
        "sample_2",
    )
    gene_comp.sink_parquet(output_file, compression='zstd', metadata=compare_metadata)

@utilities.command("prepare_profiling", help="Prepare the files needed for profiling bam files and save them in the specified output directory.")
@click.option('--reference-fasta', '-r', required=True, help="Path to the reference genome in FASTA format.")
@click.option('--gene-fasta', '-g', default=None, help="Optional path to the gene annotations in FASTA format.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--error-rate', '-e', default=0.001, show_default=True, help="Error rate for the sequencing technology when building the null model.")
@click.option('--max-total-reads', '-m', default=10000, show_default=True, help="Maximum coverage to consider when building the null model.")
@click.option('--p-threshold', '-p', default=0.05, show_default=True, help="Significance threshold for the null model.")
@click.option('--model-type', '-t', default="poisson", show_default=True, type=click.Choice(['poisson']), help="Type of null model to build.")
@click.option('--output-dir', '-o', required=True, help="Directory to save the profiling database.")
def prepare_profiling(reference_fasta, gene_fasta, stb_file, error_rate, max_total_reads, p_threshold, model_type, output_dir):
    """
    Prepare the files needed for profiling bam files and save them in the specified output directory.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(reference_fasta, output_dir / "reference.fasta")
    pf.prepare_profiling_assets(
        reference_fasta=reference_fasta,
        stb_file=stb_file,
        output_dir=output_dir,
        gene_fasta=gene_fasta,
        error_rate=error_rate,
        max_total_reads=max_total_reads,
        p_threshold=p_threshold,
        model_type=model_type,
    )


@utilities.command("profile-single")
@click.option('--reference-fasta', '-r', default=None, help="Optional reference FASTA used for mpileup. When provided, ref_base_bitmask is added to the profile.")
@click.option('--bed-file', '-b', required=True, help="Path to the BED file describing regions to be profiled.")
@click.option('--bam-file', '-a', required=True, help="Path to the BAM file to be profiled.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--null-model', '-m', required=True, help="Path to the null model parquet file.") 
@click.option('--gene-range-table', '-g', default=None, help="Optional path to the gene range table.")
@click.option('--profiling-contract', default=None, help="Optional profiling_contract.json from prepare_profiling. When provided, its hashes are written into the profile parquet metadata.")
@click.option('--num-chunks', '-n', default=24, show_default=True, help="Number of BED chunks to create for profiling.")
@click.option('--max-concurrency', '-c', default=4, show_default=True, help="Maximum number of profiling chunks to run concurrently.")
@click.option('--min-mapq', default=pf.PROFILE_MIN_MAPQ_DEFAULT, show_default=True, type=int, help="Minimum mapping quality for a read to be used during profiling.")
@click.option('--min-baseq', default=pf.PROFILE_MIN_BASEQ_DEFAULT, show_default=True, type=int, help="Minimum base quality for a base to be counted during profiling.")
@click.option('--min-read-ani', default=pf.PROFILE_MIN_READ_ANI_DEFAULT, show_default=True, type=float, help="Minimum read ANI (from the NM tag / aligned span) to use a read; filters low-identity/mis-mapped reads. Reads lacking an NM tag are kept. Pass 0 to disable.")
@click.option('--read-inclusion', default=pf.PROFILE_READ_INCLUSION_DEFAULT, show_default=True, type=click.Choice(pf.PROFILE_READ_INCLUSION_CHOICES), help="Which mapped reads are eligible: 'paired' (inStrain-style paired_only) keeps a paired read only if its mate maps to the same scaffold, dropping half-mapped orphans and cross-scaffold pairs, while keeping genuinely single-end reads; 'proper-pairs' keeps only proper pairs; 'all-mapped' keeps every mapped read.")
@click.option('--output-dir', '-o', required=True, help="Directory to save the profiling output.")
def profile_single(reference_fasta, bed_file, bam_file, stb_file, null_model, gene_range_table, profiling_contract, num_chunks, max_concurrency, min_mapq, min_baseq, min_read_ani, read_inclusion, output_dir):
    """
    Profile a single BAM file using the provided BED file and optional gene range table.
    
    """
    output_dir=pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stb= pf.read_stb(stb_file)
    null_model=pl.scan_parquet(null_model)
    profile_contract_values = (
        ut.read_profile_contract_file(profiling_contract)
        if profiling_contract is not None
        else None
    )
    pf.profile_bam(
        bed_file=bed_file,
        bam_file=bam_file,
        reference_fasta=reference_fasta,
        gene_range_table=gene_range_table,
        stb=stb,
        null_model=null_model,
        output_dir=output_dir,
        num_chunks=num_chunks,
        max_concurrency=max_concurrency,
        profile_contract=profile_contract_values,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        min_read_ani=min_read_ani,
        read_inclusion=read_inclusion,
    )

@cli.command("map", cls=SectionedCommand, short_help="Map reads to sorted BAMs.")
@click.option('--reads-table', '-i', required=True, help="CSV of reads to map, with columns 'sample_name,reads1[,reads2]' (reads2 blank/absent for single-end).")
@click.option('--output-dir', '-o', required=True, help="Directory to write BAMs, the reference FASTA/STB, and a samples.txt ready for `zipstrain profile`.")
@click.option('--reference-fasta', '-f', default=None, help="Reference FASTA to map against. If omitted, Sylph automatically picks and builds a reference from the reads.")
@click.option('--stb-file', '-s', default=None, help="Scaffold-to-genome mapping file. Required when --reference-fasta is provided.")
@click.option('--sylph-db', default=None, help="Path to the Sylph database. Used when no --reference-fasta is given; downloaded from --sylph-db-url if the path does not exist.")
@click.option('--sylph-db-url', default=mapping.DEFAULT_SYLPH_DB_URL, show_default=True, help="URL to download the Sylph database from when --sylph-db is missing.")
@click.option('--genome-cache-dir', default=None, help="Directory that caches genome FASTAs downloaded during Sylph-based reference building. Required when no --reference-fasta is given.")
@click.option('--predict-genes', is_flag=True, default=False, show_default=True, help="Also run prodigal to emit a gene FASTA (for gene-level profiling via `profile --gene-fasta`).")
@click.option('--non-competitive', is_flag=True, default=False, show_default=True, help="Pass -a to Bowtie2 for non-competitive mapping (report all alignments).")
@click.option('--force', is_flag=True, default=False, show_default=True, help="Redo every step from scratch, ignoring cached outputs. By default `map` resumes: completed Sylph tables, reference, index, and BAMs are reused.")
@click.option('--threads', '-t', default=4, show_default=True, help="Threads for Sylph, Bowtie2, and samtools.")
def map_command(reads_table, output_dir, reference_fasta, stb_file, sylph_db, sylph_db_url, genome_cache_dir, predict_genes, non_competitive, force, threads):
    """
    Map sequencing reads to BAM files, ready for `zipstrain profile`.

    Provide a reads table (``sample_name,reads1[,reads2]``). If you do not pass
    ``--reference-fasta``, ZipStrain runs Sylph to pick reference genomes from
    the reads automatically, downloading and caching them, then maps against the
    built reference. Outputs sorted, indexed BAMs, the reference FASTA + STB, and
    a ``samples.txt`` you can hand straight to ``zipstrain profile``.
    """
    console = Console()
    start = time.monotonic()
    step_number = {"n": 0}

    def _on_step(message: str) -> None:
        step_number["n"] += 1
        console.print(f"[bold cyan]›[/] [{step_number['n']}] {message}")

    console.print(Panel.fit("[bold magenta]ZipStrain map[/]", border_style="magenta"))
    with RunLogger(output_dir, command="map", argv=sys.argv) as run_log:
        try:
            results = mapping.run_map(
                reads_table=reads_table,
                output_dir=output_dir,
                reference_fasta=reference_fasta,
                stb_file=stb_file,
                sylph_db=sylph_db,
                sylph_db_url=sylph_db_url,
                genome_cache_dir=genome_cache_dir,
                threads=threads,
                predict_genes_flag=predict_genes,
                non_competitive=non_competitive,
                force=force,
                progress_callback=run_log.as_callback(_on_step),
            )
        except (ValueError, RuntimeError) as exc:
            raise click.UsageError(str(exc)) from exc

    elapsed = time.monotonic() - start
    elapsed_str = f"{int(elapsed // 3600)}:{int(elapsed % 3600 // 60):02d}:{int(elapsed % 60):02d}"
    out_dir = pathlib.Path(output_dir).absolute()

    lines = [
        "[bold green]Mapping complete![/]",
        "",
        f"Elapsed: {elapsed_str}",
        "",
        f"Output:    {out_dir}",
        f"Reference: {pathlib.Path(results['reference_fasta']).absolute()}",
        f"STB:       {pathlib.Path(results['stb_file']).absolute()}",
    ]
    if "gene_fasta" in results:
        lines.append(f"Genes:     {pathlib.Path(results['gene_fasta']).absolute()}")
    lines += [
        f"Samples:   {pathlib.Path(results['samples_table']).absolute()}",
        "",
        "Next:      zipstrain profile "
        f"--input-table {results['samples_table']} "
        f"--reference-fasta {results['reference_fasta']} "
        f"--stb-file {results['stb_file']} --run-dir <run_dir>",
    ]
    console.print(Panel("\n".join(lines), title="Summary", border_style="green", expand=True))


map_command.option_sections = {
    "Required inputs": [
        "reads_table",
        "output_dir",
    ],
    "Reference (Sylph auto-picks one if omitted)": [
        "reference_fasta",
        "stb_file",
        "sylph_db",
        "sylph_db_url",
        "genome_cache_dir",
    ],
    "Options": [
        "predict_genes",
        "non_competitive",
        "force",
        "threads",
    ],
}


@cli.command("profile", cls=SectionedCommand, short_help="Profile BAMs at nucleotide resolution.")
@click.option('--input-table', '-i', required=True, help="Path to the input table in TSV format containing sample names and paths to bam files.")
@click.option('--reference-fasta', '-f', default=None, help="Reference FASTA. Used for mpileup (adds ref_base_bitmask) and required to auto-generate the bed/genome-length assets when they are not supplied.")
@click.option('--stb-file', '-s', required=True, help="Path to the scaffold-to-genome mapping file.")
@click.option('--run-dir', '-r', required=True, help="Directory to save the run data (sample outputs, profiling_assets, and logs).")
@click.option('--null-model', '-u', default=None, help="Pre-built null model parquet file. Auto-generated into <run-dir>/profiling_assets if not provided.")
@click.option('--gene-fasta', default=None, help="Gene FASTA. When provided, a gene range table is auto-generated from it for gene-level profiling.")
@click.option('--gene-range-table', '-g', default=None, help="Pre-built gene range table file. Overrides --gene-fasta auto-generation.")
@click.option('--profiling-contract', default=None, help="Pre-built profiling_contract.json. When provided, its hashes are written into each profile parquet metadata. Auto-generated otherwise.")
@click.option('--bed-file', '-b', default=None, help="Pre-built BED file for profiling regions. Auto-generated into <run-dir>/profiling_assets if not provided.")
@click.option('--genome-length-file', '-l', default=None, help="Pre-built genome length file. Auto-generated into <run-dir>/profiling_assets if not provided.")
@click.option('--error-rate', default=0.001, show_default=True, help="Error rate used when auto-generating the null model.")
@click.option('--max-total-reads', default=10000, show_default=True, help="Maximum coverage considered when auto-generating the null model.")
@click.option('--p-threshold', default=0.05, show_default=True, help="Significance threshold used when auto-generating the null model.")
@click.option('--model-type', default="poisson", show_default=True, type=click.Choice(['poisson']), help="Null model type used when auto-generating the null model.")
@click.option('--force-prepare', is_flag=True, default=False, show_default=True, help="Regenerate all auto-generated profiling assets even if valid cached copies exist.")
@click.option('--min-mapq', default=pf.PROFILE_MIN_MAPQ_DEFAULT, show_default=True, type=int, help="Minimum mapping quality for a read to be used during profiling.")
@click.option('--min-baseq', default=pf.PROFILE_MIN_BASEQ_DEFAULT, show_default=True, type=int, help="Minimum base quality for a base to be counted during profiling.")
@click.option('--min-read-ani', default=pf.PROFILE_MIN_READ_ANI_DEFAULT, show_default=True, type=float, help="Minimum read ANI (from the NM tag / aligned span) to use a read; filters low-identity/mis-mapped reads. Reads lacking an NM tag are kept. Pass 0 to disable.")
@click.option('--read-inclusion', default=pf.PROFILE_READ_INCLUSION_DEFAULT, show_default=True, type=click.Choice(pf.PROFILE_READ_INCLUSION_CHOICES), help="Which mapped reads are eligible: 'paired' (inStrain-style paired_only) keeps a paired read only if its mate maps to the same scaffold, dropping half-mapped orphans and cross-scaffold pairs, while keeping genuinely single-end reads; 'proper-pairs' keeps only proper pairs; 'all-mapped' keeps every mapped read.")
@click.option('--num-procs', '-n', default=8, show_default=True, help="Number of processors to use for each profiling task.")
@click.option('--max-concurrent-batches', '-m', default=5, show_default=True, help="Maximum number of concurrent batches to run.")
@click.option('--poll-interval', '-p', default=1, show_default=True, help="Polling interval in seconds to check the status of batches.")
@click.option('--task-per-batch', '-t', default=10, show_default=True, help="Number of tasks to include in each batch.")
@click.option('--execution-mode', '-e', default="local", show_default=True, help="Execution mode: 'local' or 'slurm'.")
@click.option('--slurm-config', '-c', default=None, help="Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.")
@click.option('--container-engine', '-o', default="local", show_default=True, help="Container engine to use: 'local', 'docker' or 'apptainer'.")
@click.option('--container-address', default=None, help="Optional container image/address override. Defaults to the current ZipStrain version tag for docker/apptainer.")
@click.option('--no-snvs', is_flag=True, default=False, show_default=True, help="Do not call SNVs/SNPs (per-sample <sample>_SNVs.parquet). SNV calling needs --reference-fasta.")
@click.option('--snv-min-cov', default=5, show_default=True, help="Minimum coverage for a site to be eligible as an SNV/SNP call.")
@click.option('--presence-ber', default=0.5, show_default=True, help="Breadth-error-ratio threshold for the genome present/absent call (the Metapresence paper recommends ~0.8).")
@click.option('--presence-fug', default=1.0, show_default=True, help="FUG threshold for the present/absent call at low coverage. A genome is present when fug/0.632 exceeds this (fug ~ 0.632 under uniform coverage, so 1.0 means at least as uniform as random).")
@click.option('--presence-min-cov-use-fug', default=2.0, show_default=True, help="Coverage above which the present/absent call uses BER alone (below it, FUG is also required).")
@click.option('--presence-min-coverage', default=0.1, show_default=True, help="Minimum mean coverage required to call a genome present.")
@click.option('--genome-taxonomy', default=None, help="Optional genome->taxonomy TSV to add a genome_taxonomy column to genome_stats. Auto-discovered next to the reference/STB when produced by `zipstrain map` (Sylph route).")
def profile(input_table, reference_fasta, stb_file, null_model, gene_fasta, gene_range_table, profiling_contract, bed_file, genome_length_file, error_rate, max_total_reads, p_threshold, model_type, force_prepare, run_dir, num_procs, max_concurrent_batches, poll_interval, execution_mode, slurm_config, container_engine, container_address, task_per_batch, min_mapq, min_baseq, min_read_ani, read_inclusion, no_snvs, snv_min_cov, presence_ber, presence_fug, presence_min_cov_use_fug, presence_min_coverage, genome_taxonomy):
    """
    Run BAM file profiling in batches using the specified execution mode and container engine.

    Any profiling assets (null model, bed file, genome length table, gene range
    table, profiling contract) that are not supplied explicitly are generated
    automatically into a ``profiling_assets`` directory inside ``run-dir`` and
    reused on subsequent runs when the inputs are unchanged. This means a
    minimal run needs only ``--input-table``, ``--reference-fasta``, and
    ``--stb-file``.
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
    with RunLogger(run_dir, command="profile", argv=sys.argv) as run_log:
        slurm_conf = None
        if execution_mode == "slurm":
            if slurm_config is None:
                raise ValueError("SLURM configuration file must be provided when execution mode is 'slurm'.")
            slurm_conf = tm.SlurmConfig.from_json(slurm_config)

        run_log.step("Resolving profiling assets")
        assets = pf.resolve_profiling_assets(
            run_dir=run_dir,
            reference_fasta=reference_fasta,
            stb_file=stb_file,
            gene_fasta=gene_fasta,
            null_model_file=null_model,
            bed_file=bed_file,
            genome_length_file=genome_length_file,
            gene_range_table=gene_range_table,
            profiling_contract_file=profiling_contract,
            error_rate=error_rate,
            max_total_reads=max_total_reads,
            p_threshold=p_threshold,
            model_type=model_type,
            force_prepare=force_prepare,
        )

        container_engine_obj = _build_container_engine(container_engine, container_address)

        run_log.step("Profiling BAM files")
        tm.lazy_run_profile(
            run_dir=run_dir,
            container_engine=container_engine_obj,
            bams_lf=bams_lf,
            reference_fasta_file=pathlib.Path(reference_fasta) if reference_fasta is not None else None,
            stb_file=pathlib.Path(stb_file),
            null_model_file=assets.null_model_file,
            gene_range_table=assets.gene_range_table,
            profiling_contract_file=assets.profiling_contract_file,
            bed_file=assets.bed_file,
            genome_length_file=assets.genome_length_file,
            num_procs=num_procs,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            min_read_ani=min_read_ani,
            read_inclusion=read_inclusion,
            tasks_per_batch=task_per_batch,
            max_concurrent_batches=max_concurrent_batches,
            poll_interval=poll_interval,
            execution_mode=execution_mode,
            slurm_config=slurm_conf,
        )

        run_log.step("Finalizing outputs (presence, SNVs, taxonomy)")
        _finalize_profile_outputs(
            run_dir,
            emit_snvs=not no_snvs,
            snv_min_cov=snv_min_cov,
            presence_ber=presence_ber,
            presence_fug=presence_fug,
            presence_min_cov_use_fug=presence_min_cov_use_fug,
            presence_min_coverage=presence_min_coverage,
            taxonomy_file=_discover_taxonomy_file(reference_fasta, stb_file, genome_taxonomy),
            console=Console(),
        )


profile.option_sections = {
    "Required inputs": [
        "input_table",
        "reference_fasta",
        "stb_file",
        "run_dir",
    ],
    "Optional inputs": [
        "gene_fasta",
    ],
    "Optional pre-built assets (auto-generated if omitted)": [
        "null_model",
        "gene_range_table",
        "profiling_contract",
        "bed_file",
        "genome_length_file",
    ],
    "Profiling parameters": [
        "error_rate",
        "max_total_reads",
        "p_threshold",
        "model_type",
        "force_prepare",
        "min_mapq",
        "min_baseq",
        "min_read_ani",
        "read_inclusion",
    ],
    "SNV calling and presence": [
        "no_snvs",
        "snv_min_cov",
        "presence_ber",
        "presence_fug",
        "presence_min_cov_use_fug",
        "presence_min_coverage",
        "genome_taxonomy",
    ],
    "Running parameters": [
        "num_procs",
        "max_concurrent_batches",
        "poll_interval",
        "task_per_batch",
        "execution_mode",
        "slurm_config",
        "container_engine",
        "container_address",
    ],
}


def _run_matrix_compare_method(
    *,
    profile_database,
    run_dir,
    stb_file,
    bed_file,
    gene_range_table,
    scope,
    backend,
    memory_limit_gb,
    compare_genes,
    no_csv=False,
    force_csv=False,
    run_log=None,
):
    """Drive the matrix-store comparison route with progress + a summary panel."""
    profiles = [
        (row["profile_name"], row["profile_location"])
        for row in profile_database.db.collect().iter_rows(named=True)
    ]

    console = Console()
    start = time.monotonic()
    step_number = {"n": 0}

    def _on_step(message: str) -> None:
        step_number["n"] += 1
        console.print(f"[bold cyan]›[/] [{step_number['n']}] {message}")
        if run_log is not None:
            run_log.step(message)

    console.print(Panel.fit("[bold magenta]ZipStrain compare (matrix)[/]", border_style="magenta"))
    try:
        output = matrix_workflow.run_matrix_compare(
            profiles=profiles,
            run_dir=run_dir,
            stb_file=stb_file,
            bed_file=bed_file,
            gene_range_table=gene_range_table,
            scope=scope if scope is not None else "all",
            backend=backend,
            memory_limit_gb=memory_limit_gb,
            compare_genes=compare_genes,
            progress_callback=_on_step,
        )
    except ModuleNotFoundError as exc:
        raise click.UsageError(
            "The matrix method needs the matrix extra (h5py, torch). Install it with "
            f'`pip install "zipstrain[matrix]"`. Missing: {exc.name}'
        ) from exc
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    if not compare_genes:
        _add_percent_compared(
            output,
            _discover_genome_lengths(run_dir, [loc for _name, loc in profiles]),
        )
        _add_coverage_overlap(output, profiles)
    _maybe_write_csv(output, no_csv=no_csv, force_csv=force_csv, console=console)

    elapsed = time.monotonic() - start
    elapsed_str = f"{int(elapsed // 3600)}:{int(elapsed % 3600 // 60):02d}:{int(elapsed % 60):02d}"
    console.print(
        Panel(
            "\n".join(
                [
                    "[bold green]Comparison complete![/]",
                    "",
                    f"Elapsed: {elapsed_str}",
                    "",
                    f"Output:  {pathlib.Path(output).absolute()}",
                    f"Store:   {(pathlib.Path(run_dir) / 'intermediate_files' / matrix_workflow.MATRIX_STORE_NAME).absolute()}",
                ]
            ),
            title="Summary",
            border_style="green",
            expand=True,
        )
    )


@cli.command("compare", cls=SectionedCommand, short_help="Compare samples (genomes or genes).")
@click.option("--profile-db", required=True, help="Profiles to compare: either a CSV with 'profile_name,profile_location' columns (built in memory, no build-profile-db needed) or a pre-built profile-database parquet.")
@click.option("--run-dir", "-r", required=True, help="Directory to save the run data.")
@click.option("--method", type=click.Choice(["standard", "matrix"]), default="standard", show_default=True, help="Comparison engine: 'standard' (direct pairwise) or 'matrix' (reusable matrix store, good for repeated all-vs-all).")
@click.option("--compare-genes", is_flag=True, default=False, show_default=True, help="Compare genes instead of genomes.")
@click.option("--scope", default=None, help="Comparison scope. Defaults to 'all' for genomes and 'all:all' for genes.")
@click.option("--min-cov", default=5, show_default=True, help="Minimum coverage to consider a position.")
@click.option("--min-gene-compare-len", default=100, show_default=True, help="Minimum gene length to consider for comparison.")
@click.option("--stb-file", default=None, help="Scaffold-to-genome mapping file. Required for --method matrix.")
@click.option("--comp-db-file", default=None, help="Optional existing comparison parquet to resume/extend (standard method). Auto-detected from --run-dir if omitted.")
@click.option("--allow-mismatch", is_flag=True, default=False, show_default=True, help="Skip profile contract validation when building the profile database from a CSV.")
@click.option("--ani-method", "-a", default="popani", show_default=True, help="ANI calculation method (e.g., 'popani', 'conani', 'cosani_0.4').")
@click.option("--engine", type=click.Choice(["polars", "duckdb"]), default="polars", show_default=True, help="Comparison engine for standard compare tasks.")
@click.option("--calculate", default="all", show_default=True, help="Genome metrics to compute (genome mode only): ani, conani, ibs, identical_genes. Combine with '+', or use all.")
@click.option("--bed-file", default=None, help="BED file for the matrix store (--method matrix). Auto-discovered from profiling_assets if omitted.")
@click.option("--gene-range-table", default=None, help="Gene range table for gene ANI (--method matrix). Auto-discovered from profiling_assets if omitted.")
@click.option("--backend", type=click.Choice(mp.MATRIX_PAIR_BACKENDS), default="numpy", show_default=True, help="Compute backend for --method matrix (numpy, or torch on CPU/CUDA/MPS).")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True, help="Approximate memory budget for --method matrix.")
@click.option("--duckdb-memory-limit", "-d", default=None, help="DuckDB memory limit for compare tasks (e.g., 2GB).")
@click.option("--duckdb-threads", type=int, default=None, help="Number of DuckDB worker threads for compare tasks.")
@click.option("--max-concurrent-batches", "-m", default=5, show_default=True, help="Maximum number of concurrent batches to run.")
@click.option("--poll-interval", "-p", default=1, show_default=True, help="Polling interval in seconds to check the status of batches.")
@click.option("--task-per-batch", "-t", default=10, show_default=True, help="Number of tasks to include in each batch.")
@click.option("--execution-mode", "-e", default="local", show_default=True, help="Execution mode: 'local' or 'slurm'.")
@click.option("--slurm-config", "-s", default=None, help="Path to the SLURM configuration file in json format. Required if execution mode is 'slurm'.")
@click.option("--container-engine", "-c", default="local", show_default=True, help="Container engine to use: 'local', 'docker' or 'apptainer'.")
@click.option("--container-address", default=None, help="Optional container image/address override. Defaults to the current ZipStrain version tag for docker/apptainer.")
@click.option("--no-csv", is_flag=True, default=False, show_default=True, help="Do not write a companion .csv next to the comparison parquet.")
@click.option("--force-csv", is_flag=True, default=False, show_default=True, help="Write the companion .csv even when the estimated size exceeds 100 MB.")
def compare(profile_db, run_dir, method, compare_genes, scope, min_cov, min_gene_compare_len, stb_file, comp_db_file, allow_mismatch, ani_method, engine, calculate, bed_file, gene_range_table, backend, memory_limit_gb, duckdb_memory_limit, duckdb_threads, max_concurrent_batches, poll_interval, task_per_batch, execution_mode, slurm_config, container_engine, container_address, no_csv, force_csv):
    """
    Compare profiled samples at the genome level (default) or gene level (--compare-genes).

    ``--profile-db`` may be a CSV of ``profile_name,profile_location`` rows, so
    there is no need to run ``zipstrain utilities build-profile-db`` first; a
    pre-built profile-database parquet is also accepted.

    Both methods write ``<run-dir>/all_comparisons.parquet``. Re-running with the
    same ``--run-dir`` and a profiles table that includes new samples extends the
    existing comparison, computing only the new pairs.
    """
    run_dir = pathlib.Path(run_dir)
    with RunLogger(run_dir, command="compare", argv=sys.argv) as run_log:
        run_log.step("Loading profile database")
        profile_database = _load_profile_database(profile_db, allow_mismatch=allow_mismatch)

        if method == "matrix":
            _run_matrix_compare_method(
                profile_database=profile_database,
                run_dir=run_dir,
                stb_file=stb_file,
                bed_file=bed_file,
                gene_range_table=gene_range_table,
                scope=scope,
                backend=backend,
                memory_limit_gb=memory_limit_gb,
                compare_genes=compare_genes,
                no_csv=no_csv,
                force_csv=force_csv,
                run_log=run_log,
            )
            return

        if duckdb_threads is not None and duckdb_threads < 1:
            raise ValueError("--duckdb-threads must be >= 1")
        slurm_conf = None
        if execution_mode == "slurm":
            if slurm_config is None:
                raise ValueError("SLURM configuration file must be provided when execution mode is 'slurm'.")
            slurm_conf = tm.SlurmConfig.from_json(slurm_config)

        # Standard method: auto-resume from a prior run in the same run-dir when the
        # user did not pass an explicit comparison to extend.
        if comp_db_file is None:
            existing_output = run_dir / ("all_gene_comparisons.parquet" if compare_genes else "all_comparisons.parquet")
            if existing_output.exists():
                comp_db_file = str(existing_output)

        container_engine_obj = _build_container_engine(container_engine, container_address)

        if compare_genes:
            resolved_scope = scope if scope is not None else "all:all"
            comps_db = db.GeneComparisonDatabase(
                profile_db=profile_database,
                config=db.GeneComparisonConfig(
                    scope=resolved_scope,
                    min_cov=min_cov,
                    min_gene_compare_len=min_gene_compare_len,
                    stb_file_loc=stb_file,
                ),
                comp_db_loc=comp_db_file,
            )
            run_log.step("Comparing genes across sample pairs")
            tm.lazy_run_gene_compares(
                comps_db=comps_db,
                container_engine=container_engine_obj,
                run_dir=run_dir,
                max_concurrent_batches=max_concurrent_batches,
                execution_mode=execution_mode,
                slurm_config=slurm_conf,
                compare_engine=engine,
                tasks_per_batch=task_per_batch,
                poll_interval=poll_interval,
                ani_method=ani_method,
                duckdb_memory_limit=duckdb_memory_limit,
                duckdb_threads=duckdb_threads,
            )
        else:
            resolved_scope = scope if scope is not None else "all"
            cp.parse_genome_calculations(calculate)
            comps_db = db.GenomeComparisonDatabase(
                profile_db=profile_database,
                config=db.GenomeComparisonConfig(
                    scope=resolved_scope,
                    min_cov=min_cov,
                    min_gene_compare_len=min_gene_compare_len,
                    stb_file_loc=stb_file,
                ),
                comp_db_loc=comp_db_file,
            )
            run_log.step("Comparing genomes across sample pairs")
            tm.lazy_run_compares(
                comps_db=comps_db,
                container_engine=container_engine_obj,
                run_dir=run_dir,
                max_concurrent_batches=max_concurrent_batches,
                execution_mode=execution_mode,
                slurm_config=slurm_conf,
                ani_method=ani_method,
                compare_engine=engine,
                calculate=calculate,
                duckdb_memory_limit=duckdb_memory_limit,
                duckdb_threads=duckdb_threads,
                tasks_per_batch=task_per_batch,
                poll_interval=poll_interval,
            )

        run_log.step("Writing comparison outputs")
        output_name = "all_gene_comparisons.parquet" if compare_genes else "all_comparisons.parquet"
        if not compare_genes:
            db_rows = profile_database.db.select(["profile_name", "profile_location"]).collect()
            profile_pairs = list(zip(db_rows["profile_name"].to_list(), db_rows["profile_location"].to_list()))
            _add_percent_compared(
                run_dir / output_name,
                _discover_genome_lengths(run_dir, [loc for _n, loc in profile_pairs]),
            )
            _add_coverage_overlap(run_dir / output_name, profile_pairs)
        _maybe_write_csv(run_dir / output_name, no_csv=no_csv, force_csv=force_csv, console=Console())


compare.option_sections = {
    "Required inputs": [
        "profile_db",
        "run_dir",
    ],
    "Comparison parameters": [
        "method",
        "compare_genes",
        "scope",
        "min_cov",
        "min_gene_compare_len",
        "stb_file",
        "ani_method",
        "calculate",
        "comp_db_file",
        "allow_mismatch",
    ],
    "Matrix method (--method matrix)": [
        "bed_file",
        "gene_range_table",
        "backend",
        "memory_limit_gb",
    ],
    "Output": [
        "no_csv",
        "force_csv",
    ],
    "Standard method / engine": [
        "engine",
        "duckdb_memory_limit",
        "duckdb_threads",
        "max_concurrent_batches",
        "poll_interval",
        "task_per_batch",
        "execution_mode",
        "slurm_config",
        "container_engine",
        "container_address",
    ],
}

@cli.command("test", short_help="Check your environment is ready.")
def test():
    """Run a lightweight ZipStrain health check."""
    hc.render_health_report(hc.collect_health_report(), console=Console())

if __name__ == "__main__":
    cli()
