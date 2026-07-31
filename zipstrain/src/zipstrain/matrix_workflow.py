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

from zipstrain import compare as cp
from zipstrain import matrix_pairs as mp
from zipstrain import utils as ut


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


def _matrix_compare_db_path(
    intermediate: pathlib.Path,
    ani_method: str,
    ani_methods: tuple[str, ...],
) -> pathlib.Path:
    if len(ani_methods) == 1:
        return intermediate / MATRIX_COMPARE_DB_NAME
    return intermediate / f"matrix_compare_{cp.ani_method_suffix(ani_method)}.duckdb"


def _matrix_method_export_path(
    intermediate: pathlib.Path,
    ani_method: str,
) -> pathlib.Path:
    return intermediate / f"matrix_compare_{cp.ani_method_suffix(ani_method)}.parquet"


def _combine_matrix_method_exports(
    *,
    method_exports: list[tuple[str, pathlib.Path]],
    output: pathlib.Path,
    table: str,
    calculate: str,
) -> pathlib.Path:
    """Join compact per-method matrix results into one public comparison table."""
    if not method_exports:
        raise ValueError("At least one matrix ANI result is required.")
    ani_methods = tuple(method for method, _path in method_exports)
    keys = ["sample_1", "sample_2", "genome"]
    if table == "gene":
        keys.append("gene")

    combined: pl.LazyFrame | None = None
    for index, (method, path) in enumerate(method_exports):
        frame = pl.scan_parquet(path)
        schema_names = frame.collect_schema().names()
        selects = [pl.col(column) for column in keys]
        if index == 0 and "total_positions" in schema_names:
            selects.append(pl.col("total_positions"))

        if table == "genome":
            metric_bases = (
                "share_allele_pos",
                "genome_ani",
                "max_consecutive_length",
            )
        else:
            metric_bases = ("gene_pop_ani",)
        for base in metric_bases:
            if base not in schema_names:
                continue
            public_base = "gene_ani" if base == "gene_pop_ani" else base
            selects.append(
                pl.col(base).alias(
                    cp.ani_metric_column(public_base, method, ani_methods)
                )
            )
        method_frame = frame.select(selects)
        combined = (
            method_frame
            if combined is None
            else combined.join(method_frame, on=keys, how="inner", validate="1:1")
        )

    metadata = {
        key: value
        for key, value in pl.read_parquet_metadata(method_exports[0][1]).items()
        if key != "ARROW:schema"
    }
    metadata[ut.COMPARE_ANI_METHOD_METADATA_KEY] = ",".join(ani_methods)
    metadata[ut.COMPARE_CALCULATE_METADATA_KEY] = calculate
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()
    assert combined is not None
    combined_columns = combined.collect_schema().names()
    ordered_columns = [*keys]
    if "total_positions" in combined_columns:
        ordered_columns.append("total_positions")
    if table == "genome":
        for method in ani_methods:
            ordered_columns.extend(
                column
                for column in (
                    cp.ani_metric_column(
                        "share_allele_pos", method, ani_methods
                    ),
                    cp.ani_metric_column("genome_ani", method, ani_methods),
                )
                if column in combined_columns
            )
        ordered_columns.extend(
            cp.ani_metric_column("max_consecutive_length", method, ani_methods)
            for method in ani_methods
            if cp.ani_metric_column(
                "max_consecutive_length", method, ani_methods
            )
            in combined_columns
        )
    else:
        ordered_columns.extend(
            cp.ani_metric_column("gene_ani", method, ani_methods)
            for method in ani_methods
            if cp.ani_metric_column("gene_ani", method, ani_methods)
            in combined_columns
        )
    combined = combined.select(ordered_columns)
    combined.sink_parquet(output, compression="zstd", metadata=metadata)
    return output


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
    calculations: tuple[str, ...],
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
    calculate: str = "all",
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
    ani_methods = cp.parse_ani_methods(ani_method)

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
    calculations = mp.parse_matrix_calculations(
        calculate,
        include_gene_from_all=gene_range_table is not None,
    )
    if "gene" in calculations and gene_range_table is None:
        raise ValueError("Gene comparison requires --gene-range-table.")

    cleaned_stb = _clean_stb(pathlib.Path(stb_file), intermediate / "reference_cleaned.stb")

    if store.exists():
        for method in ani_methods:
            _validate_existing_matrix_run(
                store=store,
                compare_db=_matrix_compare_db_path(intermediate, method, ani_methods),
                scope=scope,
                min_cov=min_cov,
                ani_method=method,
                calculations=calculations,
            )
    existing_names = _existing_store_sample_names(store)
    new_profiles = [(name, location) for name, location in profiles if name not in existing_names]
    stage_dir = intermediate / MATRIX_STAGE_DIRNAME

    if not store.exists():
        _notify(progress_callback, f"Building matrix store from {len(new_profiles)} profile(s)")
        _stage_profiles(new_profiles, stage_dir)
        parsed_methods = [mp.parse_matrix_ani_method(method)[0] for method in ani_methods]
        storage_mode = (
            mp.MATRIX_STORAGE_BITMASK
            if parsed_methods == ["popani"]
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

    table = "gene" if "gene" in calculations else "genome"
    output = run_dir / GENOME_OUTPUT_NAME
    method_exports: list[tuple[str, pathlib.Path]] = []
    for method in ani_methods:
        compare_db = _matrix_compare_db_path(intermediate, method, ani_methods)
        _notify(
            progress_callback,
            f"Comparing all-vs-all with {method} (backend={backend})",
        )
        mp.matrix_compare(
            matrix_db_file=store,
            output_file=compare_db,
            min_cov=min_cov,
            genome=scope,
            memory_limit_gb=memory_limit_gb,
            backend=backend,
            ani_method=method,
            calculate="+".join(calculations),
            emit_writer_logs=False,
            progress_callback=None,
        )

        if len(ani_methods) == 1:
            if output.exists():
                output.unlink()
            _notify(progress_callback, f"Exporting {table}-level results")
            mp.export_matrix_compare_parquet(
                matrix_compare_db_file=compare_db,
                output_file=output,
                table=table,
            )
            return output

        method_output = _matrix_method_export_path(intermediate, method)
        if method_output.exists():
            method_output.unlink()
        mp.export_matrix_compare_parquet(
            matrix_compare_db_file=compare_db,
            output_file=method_output,
            table=table,
        )
        method_exports.append((method, method_output))

    _notify(
        progress_callback,
        f"Combining {len(method_exports)} ANI methods into {table}-level results",
    )
    _combine_matrix_method_exports(
        method_exports=method_exports,
        output=output,
        table=table,
        calculate="+".join(calculations),
    )
    return output
