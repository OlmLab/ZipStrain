"""zipstrain.matrix_workflow
========================
High-level orchestration for the matrix route of ``zipstrain compare``.

``zipstrain compare --method matrix`` funnels through here: it builds (or
incrementally appends to) a matrix store from a profile database, runs the
resumable all-vs-all matrix compare, and exports a single
``all_comparisons.parquet`` — the same top-level output the standard method
produces. All the heavy intermediates (the HDF5 store, the resumable DuckDB,
staged profile symlinks) live under ``<run_dir>/intermediate_files``.
"""
from __future__ import annotations

import pathlib
import shutil
from typing import Callable

import polars as pl

from zipstrain import matrix_pairs as mp


# Sentinel genome value for profile positions that are not assigned to any
# binned genome. These cannot be compared at the genome level and are dropped
# before building the matrix store.
UNBINNED_GENOME = "NA"


ProgressCallback = Callable[[str], None]

INTERMEDIATE_DIRNAME = "intermediate_files"
MATRIX_STORE_NAME = "matrix_db.h5"
MATRIX_COMPARE_DB_NAME = "matrix_compare.duckdb"
MATRIX_STAGE_DIRNAME = "matrix_profiles"
GENOME_OUTPUT_NAME = "all_comparisons.parquet"
GENE_OUTPUT_NAME = "all_gene_comparisons.parquet"


def _notify(callback: ProgressCallback | None, message: str) -> None:
    if callback is not None:
        callback(message)


def discover_bed_file(profile_locations: list[str | pathlib.Path]) -> pathlib.Path | None:
    """Find ``profiling_assets/genomes_bed_file.bed`` near any profile."""
    return _discover_asset(profile_locations, "genomes_bed_file.bed")


def discover_gene_range_table(profile_locations: list[str | pathlib.Path]) -> pathlib.Path | None:
    """Find a non-empty ``profiling_assets/gene_range_table.tsv`` near any profile."""
    return _discover_asset(profile_locations, "gene_range_table.tsv", require_nonempty=True)


def _discover_asset(
    profile_locations: list[str | pathlib.Path],
    asset_name: str,
    *,
    require_nonempty: bool = False,
) -> pathlib.Path | None:
    for location in profile_locations:
        profile_path = pathlib.Path(location).resolve()
        for parent in profile_path.parents:
            candidate = parent / "profiling_assets" / asset_name
            if candidate.exists() and (not require_nonempty or candidate.stat().st_size > 0):
                return candidate
    return None


def _existing_store_sample_names(store: pathlib.Path) -> set[str]:
    if not store.exists():
        return set()
    return {name for _idx, name in mp._load_matrix_hdf5_samples(store)}


def _clean_stb(stb_file: pathlib.Path, dest: pathlib.Path) -> pathlib.Path:
    """Write a whitespace-stripped copy of an STB file.

    Some STB files (e.g. from MaxBin) pad their tab-separated columns with
    spaces. The matrix build cross-checks scaffold names against the bed file
    exactly, so we strip surrounding whitespace from both columns first.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    cleaned_lines = []
    for line in pathlib.Path(stb_file).read_text().splitlines():
        if not line.strip():
            continue
        fields = [field.strip() for field in line.split("\t")]
        cleaned_lines.append("\t".join(fields))
    dest.write_text("\n".join(cleaned_lines) + "\n")
    return dest


def _stage_profiles(
    profiles: list[tuple[str, str | pathlib.Path]], dest_dir: pathlib.Path
) -> pathlib.Path:
    """Stage each ``(name, path)`` profile into ``dest_dir`` as ``<name>.parquet``.

    Rows with the unbinned genome sentinel are dropped so the staged profiles
    match the STB/BED contract the matrix store requires. When a profile has no
    unbinned rows it is symlinked (no copy); otherwise a filtered copy is written.
    """
    shutil.rmtree(dest_dir, ignore_errors=True)
    dest_dir.mkdir(parents=True, exist_ok=True)
    for name, location in profiles:
        source = pathlib.Path(location).resolve()
        dest = dest_dir / f"{name}.parquet"
        lazy = pl.scan_parquet(source)
        has_unbinned = (
            "genome" in lazy.collect_schema().names()
            and lazy.filter(pl.col("genome") == UNBINNED_GENOME).limit(1).collect().height > 0
        )
        if has_unbinned:
            lazy.filter(pl.col("genome") != UNBINNED_GENOME).sink_parquet(dest)
        else:
            dest.symlink_to(source)
    return dest_dir


def run_matrix_compare(
    *,
    profiles: list[tuple[str, str | pathlib.Path]],
    run_dir: str | pathlib.Path,
    stb_file: str | pathlib.Path | None,
    bed_file: str | pathlib.Path | None = None,
    gene_range_table: str | pathlib.Path | None = None,
    scope: str = "all",
    backend: str = "numpy",
    memory_limit_gb: float = 16.0,
    compare_genes: bool = False,
    progress_callback: ProgressCallback | None = None,
) -> pathlib.Path:
    """Build/extend a matrix store, run resumable compare, and export to parquet.

    ``profiles`` is a list of ``(sample_name, profile_parquet_path)``. Re-running
    with a superset of profiles appends only the new samples to the existing
    store and computes only the new pairs, then re-exports the full table.
    Returns the path to the exported comparison parquet.
    """
    run_dir = pathlib.Path(run_dir)
    intermediate = run_dir / INTERMEDIATE_DIRNAME
    intermediate.mkdir(parents=True, exist_ok=True)
    store = intermediate / MATRIX_STORE_NAME
    compare_db = intermediate / MATRIX_COMPARE_DB_NAME

    if stb_file is None:
        raise ValueError("--stb-file is required for --method matrix.")

    profile_locations = [location for _name, location in profiles]
    if bed_file is None:
        bed_file = discover_bed_file(profile_locations)
        if bed_file is None:
            raise ValueError(
                "Could not auto-discover a bed file "
                "(profiling_assets/genomes_bed_file.bed near the profiles). "
                "Pass --bed-file explicitly for --method matrix."
            )
    if gene_range_table is None:
        gene_range_table = discover_gene_range_table(profile_locations)

    cleaned_stb = _clean_stb(pathlib.Path(stb_file), intermediate / "reference_cleaned.stb")

    existing_names = _existing_store_sample_names(store)
    new_profiles = [(name, location) for name, location in profiles if name not in existing_names]
    stage_dir = intermediate / MATRIX_STAGE_DIRNAME

    if not store.exists():
        _notify(progress_callback, f"Building matrix store from {len(new_profiles)} profile(s)")
        _stage_profiles(new_profiles, stage_dir)
        mp.build_matrix_hdf5(
            profile_dir=stage_dir,
            output_file=store,
            genome=scope,
            bed_file=pathlib.Path(bed_file),
            stb_file=cleaned_stb,
            gene_range_table=pathlib.Path(gene_range_table) if gene_range_table is not None else None,
            memory_limit_gb=memory_limit_gb,
        )
    elif new_profiles:
        _notify(progress_callback, f"Appending {len(new_profiles)} new profile(s) to the matrix store")
        _stage_profiles(new_profiles, stage_dir)
        mp.append_matrix_db(
            profile_dir=stage_dir,
            matrix_db_file=store,
            memory_limit_gb=memory_limit_gb,
        )
    else:
        _notify(progress_callback, "Matrix store already contains all profiles; reusing it")

    _notify(progress_callback, f"Comparing all-vs-all (backend={backend})")
    mp.matrix_compare(
        matrix_db_file=store,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        genome=scope,
        memory_limit_gb=memory_limit_gb,
        backend=backend,
        calculate="all",
        emit_writer_logs=False,
        progress_callback=None,
    )

    table = "gene" if compare_genes else "genome"
    output = run_dir / (GENE_OUTPUT_NAME if compare_genes else GENOME_OUTPUT_NAME)
    if output.exists():
        output.unlink()
    _notify(progress_callback, f"Exporting {table}-level results")
    mp.export_matrix_compare_parquet(
        matrix_compare_db_file=compare_db,
        output_file=output,
        table=table,
    )
    return output
