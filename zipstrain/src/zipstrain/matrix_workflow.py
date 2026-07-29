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

import duckdb
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


def _validate_existing_matrix_run(
    *,
    store: pathlib.Path,
    compare_db: pathlib.Path,
    scope: str,
    min_cov: int,
    ani_method: str,
) -> None:
    metadata = mp._load_matrix_hdf5_metadata(store)
    normalized_ani_method, _threshold = mp.parse_matrix_ani_method(ani_method)
    storage_mode = mp._matrix_storage_mode_from_metadata(metadata)
    if storage_mode != mp.MATRIX_STORAGE_COUNTS and normalized_ani_method != "popani":
        raise ValueError(
            f"{normalized_ani_method} requires a count matrix store. "
            "Use a different --run-dir."
        )
    try:
        stored_min_cov = int(metadata["coverage_filter_min_cov"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Existing matrix store has no valid build-time min_cov.") from exc
    if stored_min_cov != min_cov:
        raise ValueError(
            f"This run directory contains a matrix built with min_cov={stored_min_cov}; "
            f"requested min_cov={min_cov}. Use a different --run-dir."
        )

    if not compare_db.exists() or compare_db.stat().st_size == 0:
        return
    conn = duckdb.connect(str(compare_db), read_only=True)
    try:
        if not mp._matrix_compare_table_exists(conn, "matrix_compare_metadata"):
            raise ValueError(
                "Existing matrix compare DB is missing metadata and cannot be resumed safely."
            )
        compare_metadata = {
            str(key): str(value)
            for key, value in conn.execute(
                "SELECT key, value FROM matrix_compare_metadata"
            ).fetchall()
        }
    finally:
        conn.close()
    compare_metadata.setdefault("matrix_input_format", "hdf5")
    gene_ranges = mp._load_matrix_hdf5_gene_ranges(store)
    calculations = mp.parse_matrix_calculations(
        "all",
        include_gene_from_all=bool(gene_ranges),
    )
    mp._validate_matrix_compare_db_metadata(
        compare_metadata=compare_metadata,
        matrix_metadata=metadata,
        genome_scope=None if scope == "all" else scope,
        calculations=calculations,
        min_cov=min_cov,
        ani_method=normalized_ani_method,
    )


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
    min_cov: int = mp.MATRIX_BUILD_MIN_COV,
    ani_method: str = "popani",
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

    if store.exists():
        _validate_existing_matrix_run(
            store=store,
            compare_db=compare_db,
            scope=scope,
            min_cov=min_cov,
            ani_method=ani_method,
        )
    existing_names = _existing_store_sample_names(store)
    new_profiles = [(name, location) for name, location in profiles if name not in existing_names]
    stage_dir = intermediate / MATRIX_STAGE_DIRNAME

    if not store.exists():
        _notify(progress_callback, f"Building matrix store from {len(new_profiles)} profile(s)")
        _stage_profiles(new_profiles, stage_dir)
        ani_kind, _cos_threshold = mp.parse_matrix_ani_method(ani_method)
        storage_mode = (
            mp.MATRIX_STORAGE_BITMASK
            if ani_kind == "popani"
            else mp.MATRIX_STORAGE_COUNTS
        )
        mp.build_matrix_hdf5(
            profile_dir=stage_dir,
            output_file=store,
            genome=scope,
            bed_file=pathlib.Path(bed_file),
            stb_file=cleaned_stb,
            gene_range_table=pathlib.Path(gene_range_table) if gene_range_table is not None else None,
            storage_mode=storage_mode,
            count_dtype=None,
            min_cov=min_cov,
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
        min_cov=min_cov,
        genome=scope,
        memory_limit_gb=memory_limit_gb,
        backend=backend,
        ani_method=ani_method,
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
