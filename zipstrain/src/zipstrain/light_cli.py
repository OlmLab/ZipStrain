"""CLI entrypoint for zipstrain-light."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click

from zipstrain import light


@click.group()
def cli() -> None:
    """ZipStrain-Light CLI."""


@cli.command("profile")
@click.option("--output-file", "-o", required=True, help="Output DuckDB file.")
@click.option("--profile-parquet", default=None, help="Existing ZipStrain profile parquet input.")
@click.option("--bam-file", default=None, help="BAM file to profile directly.")
@click.option("--bed-file", default=None, help="BED file for profiling (required with --bam-file).")
@click.option("--stb-file", default=None, help="STB file for profiling (required with --bam-file).")
@click.option("--null-model", default=None, help="Null model parquet (required with --bam-file).")
@click.option("--gene-range-table", default=None, help="Optional gene range table TSV.")
@click.option("--num-chunks", type=int, default=24, show_default=True, help="Number of BED chunks to create in BAM profiling mode.")
@click.option("--max-concurrency", type=int, default=4, show_default=True, help="Maximum number of BAM profiling chunks to run concurrently.")
@click.option("--min-cov", type=int, default=5, show_default=True, help="Retain loci with corrected coverage strictly greater than this value.")
@click.option("--gene-stats/--no-gene-stats", default=True, show_default=True, help="Write the companion gene stats parquet into the output DuckDB file.")
@click.option("--genome-stats/--no-genome-stats", default=True, show_default=True, help="Write the companion genome stats parquet into the output DuckDB file.")
@click.option("--duckdb-memory-limit", default=None, help="DuckDB memory limit for writing profile DB (for example 2GB).")
@click.option("--duckdb-temp-directory", default=None, help="DuckDB temp directory for spill files.")
@click.option("--duckdb-threads", type=int, default=None, help="DuckDB worker thread count.")
def profile(
    output_file: str,
    profile_parquet: Optional[str],
    bam_file: Optional[str],
    bed_file: Optional[str],
    stb_file: Optional[str],
    null_model: Optional[str],
    gene_range_table: Optional[str],
    num_chunks: int,
    max_concurrency: int,
    min_cov: int,
    gene_stats: bool,
    genome_stats: bool,
    duckdb_memory_limit: Optional[str],
    duckdb_temp_directory: Optional[str],
    duckdb_threads: Optional[int],
) -> None:
    """Build a lightweight profile DuckDB file."""
    if profile_parquet is None and bam_file is None:
        raise click.UsageError("Provide one of --profile-parquet or --bam-file.")
    if profile_parquet is not None and bam_file is not None:
        raise click.UsageError("Use only one of --profile-parquet or --bam-file.")

    if bam_file is not None:
        missing = [
            name
            for name, value in (
                ("--bed-file", bed_file),
                ("--stb-file", stb_file),
                ("--null-model", null_model),
            )
            if value is None
        ]
        if missing:
            raise click.UsageError(f"Missing required options with --bam-file: {', '.join(missing)}")

        summary = light.build_light_profile_bundle_from_bam(
            bam_file=bam_file,
            bed_file=bed_file,
            gene_range_table=gene_range_table,
            stb_file=stb_file,
            null_model=null_model,
            output_file=output_file,
            min_cov=min_cov,
            include_gene_stats=gene_stats,
            include_genome_stats=genome_stats,
            num_chunks=num_chunks,
            max_concurrency=max_concurrency,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )
    else:
        summary = light.build_light_profile_bundle_from_parquet(
            profile_parquet=profile_parquet,
            output_file=output_file,
            min_cov=min_cov,
            include_gene_stats=gene_stats,
            include_genome_stats=genome_stats,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )

    click.echo(
        f"wrote={Path(output_file)} min_cov={min_cov} coverage_rows={summary.coverage_rows} "
        f"gene_stats_rows={summary.gene_stats_rows} genome_stats_rows={summary.genome_stats_rows}"
    )


@cli.command("compare")
@click.option("--profile-1", "-p1", required=True, help="First light profile DuckDB file.")
@click.option("--profile-2", "-p2", required=True, help="Second light profile DuckDB file.")
@click.option(
    "--stb-file",
    "-s",
    default=None,
    help="Unused in the current schema; kept for backward compatibility.",
)
@click.option("--output-file", "-o", required=True, help="Output comparison parquet.")
@click.option("--genome", "-g", default="all", show_default=True, help="Genome scope.")
@click.option(
    "--min-gene-compare-len",
    "-l",
    type=int,
    default=100,
    show_default=True,
    help="Minimum shared positions per gene for identical gene metrics.",
)
@click.option("--ani-method", "-a", default="popani", show_default=True, help="ANI method (currently only popani).")
@click.option(
    "--calculate",
    default="all",
    show_default=True,
    help="Metric groups to compute: ani, ibs, genes, or '+' combinations (for example ani+ibs).",
)
@click.option("--duckdb-memory-limit", default=None, help="DuckDB memory limit (for example 2GB).")
@click.option("--duckdb-temp-directory", default=None, help="DuckDB temp directory for spill files.")
@click.option("--duckdb-threads", type=int, default=None, help="DuckDB worker thread count.")
@click.option("--sample-1", default=None, help="Sample label for first profile (defaults to the input file stem).")
@click.option("--sample-2", default=None, help="Sample label for second profile (defaults to the input file stem).")
def compare(
    profile_1: str,
    profile_2: str,
    stb_file: str,
    output_file: str,
    genome: str,
    min_gene_compare_len: int,
    ani_method: str,
    calculate: str,
    duckdb_memory_limit: Optional[str],
    duckdb_temp_directory: Optional[str],
    duckdb_threads: Optional[int],
    sample_1: Optional[str],
    sample_2: Optional[str],
) -> None:
    """Compare two light profiles and write genome comparison parquet."""
    profile_1_path = Path(profile_1)
    profile_2_path = Path(profile_2)
    if not profile_1_path.exists():
        raise click.UsageError(f"--profile-1 does not exist: {profile_1_path}")
    if not profile_2_path.exists():
        raise click.UsageError(f"--profile-2 does not exist: {profile_2_path}")

    sample_1_name = sample_1 or profile_1_path.stem
    sample_2_name = sample_2 or profile_2_path.stem
    light.compare_light_genomes_to_parquet(
        profile_db_1=profile_1,
        profile_db_2=profile_2,
        stb_file=stb_file,
        output_file=output_file,
        sample_1_name=sample_1_name,
        sample_2_name=sample_2_name,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome,
        ani_method=ani_method,
        memory_limit=duckdb_memory_limit,
        temp_directory=duckdb_temp_directory,
        threads=duckdb_threads,
        calculate=calculate,
    )
    click.echo(f"wrote={Path(output_file)}")


if __name__ == "__main__":
    cli()
