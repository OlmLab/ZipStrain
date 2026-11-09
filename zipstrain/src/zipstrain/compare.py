"""zipstrain.compare
========================
This module provides all comparison functions for zipstrain.

"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable, Optional, Union

import numpy as np
import polars as pl
from zipstrain import options
from zipstrain.gene_info import _index_genes, _scaffold_runs
from zipstrain.options import PROFILE_SORTED_METADATA_KEY, PROFILE_SORTED_METADATA_VALUE


def _read_custom_parquet_metadata(source) -> dict[str, str]:
    """Read a Parquet file's ZipStrain metadata without loading :mod:`zipstrain.utils`."""
    from zipstrain.utils import _read_custom_parquet_metadata as read
    return read(source)

_PAIR_BYTES_PER_ROW = 128
_AUTO_WORKERS_MIN_ROWS = 50_000_000


def parse_genome_calculations(calculate: Optional[Union[str, Iterable[str]]] = None) -> tuple[str, ...]:
    """Return the genome metrics in a selection for routes that only calculate genome rows."""
    explicit = options.explicit_metrics(calculate)
    for metric in explicit:
        if metric in options.GENE_METRICS:
            raise ValueError(
                f"{metric} is calculated per gene, and this route only calculates genome metrics. "
                "Use --preprofile-dir with V2 profiles, or request gene_ani on its own with 1.x profiles."
            )
    return options.GENOME_METRICS if options.requests_all(calculate) else explicit


def genome_metric_output_columns(calculate: Optional[Union[str, Iterable[str]]] = None) -> list[str]:
    """Return ordered output columns for selected genome-level calculations."""
    calculations = parse_genome_calculations(calculate)
    cols = ["genome"]
    if "genome_ani" in calculations:
        cols.extend(["total_positions", "share_allele_pos", "genome_ani"])
    if "max_identical_stretch" in calculations:
        cols.append("max_identical_stretch")
    return cols


def metrics_in_frame(frame: pl.LazyFrame) -> tuple[str, ...]:
    """Return the metrics a comparison frame holds, for recording in its metadata."""
    columns = set(frame.collect_schema().names())
    markers = {"genome_ani": "genome_ani", "max_identical_stretch": "max_identical_stretch",
               "gene_ani": "gene_ani", "gene_selection": "dN"}
    return tuple(metric for metric in options.METRICS if markers[metric] in columns)


class PolarsANIExpressions:
    """ 
    Any kind of ANI calculation based on two profiles should be implemented as a method of this class.
    In defining this method, the following rules should be followed:
    
    -   The method returns a Polars expression (pl.Expr).
    
    -   When applied to a row, the method returns a zero if that position is a SNV. Otherwise it should return a number greater than zero.
    
    -   A, T, C, G columns in the first profile are named "A", "T", "C", "G" and in the second profile they are named "A_2", "T_2", "C_2", "G_2".
    
    1. popani: Population ANI based on the shared alleles between two profiles.
    2. conani: Consensus ANI based on the consensus alleles between two profiles.
    3. cosani_<threshold>: Generalized cosine similarity ANI where threshold is a float value between 0 and 1. Once the similarity is below the threshold, it is considered a SNV.
    """
    BASES = ["A", "C", "G", "T"]
    MPILE_1_BASES = ["A", "T", "C", "G"]
    MPILE_2_BASES = ["A_2", "T_2", "C_2", "G_2"]

    def _base_counts(self):
        """Return the two samples' A/C/G/T counts as Int64 Polars expressions.

        The base counts arrive as smaller integer types, so the per-position
        products below can overflow at high coverage. We cast to Int64 once
        here, up front, and the ANI formulas that follow read as plain
        arithmetic on ``sample_1``/``sample_2`` instead of repeating
        ``.cast(pl.Int64)`` on every column.

        Returns:
            tuple[dict[str, pl.Expr], dict[str, pl.Expr]]: two dicts keyed by
            base letter (``"A"``, ``"C"``, ``"G"``, ``"T"``); the first holds
            sample 1's counts and the second sample 2's.
        """
        sample_1 = {base: pl.col(base).cast(pl.Int64) for base in self.BASES}
        sample_2 = {base: pl.col(f"{base}_2").cast(pl.Int64) for base in self.BASES}
        return sample_1, sample_2

    def popani(self):
        # popANI treats a position as shared when the two samples have at least
        # one allele in common. The dot product of their per-base counts is
        # positive exactly in that case (and zero only when they share no base,
        # i.e. an SNV), which is the >0 / ==0 signal every caller expects.
        sample_1, sample_2 = self._base_counts()
        return (sample_1["A"] * sample_2["A"]
                + sample_1["C"] * sample_2["C"]
                + sample_1["G"] * sample_2["G"]
                + sample_1["T"] * sample_2["T"])

    def conani(self):
        # conANI compares only the consensus (most abundant) base in each
        # sample: the position counts as shared when the consensus base is the
        # same in both samples. `max_horizontal` gives each sample's consensus
        # count; the position matches when that consensus falls on the same base.
        consensus_1 = pl.max_horizontal(*[pl.col(base) for base in self.MPILE_1_BASES])
        consensus_2 = pl.max_horizontal(*[pl.col(base) for base in self.MPILE_2_BASES])
        return pl.when((pl.col("A") == consensus_1) & (pl.col("A_2") == consensus_2)
                       | (pl.col("T") == consensus_1) & (pl.col("T_2") == consensus_2)
                       | (pl.col("C") == consensus_1) & (pl.col("C_2") == consensus_2)
                       | (pl.col("G") == consensus_1) & (pl.col("G_2") == consensus_2)).then(1).otherwise(0)

    def generalized_cos_ani(self, threshold: float = 0.4):
        # Treat each position as a 4-vector of A/C/G/T counts per sample and
        # take the cosine similarity between the two samples' vectors; a
        # position counts as shared when that similarity clears `threshold`,
        # and is otherwise treated as an SNV.
        sample_1, sample_2 = self._base_counts()

        dot_product = (sample_1["A"] * sample_2["A"]
                       + sample_1["C"] * sample_2["C"]
                       + sample_1["G"] * sample_2["G"]
                       + sample_1["T"] * sample_2["T"])
        magnitude_1 = (sample_1["A"]**2 + sample_1["C"]**2 + sample_1["G"]**2 + sample_1["T"]**2)**0.5
        magnitude_2 = (sample_2["A"]**2 + sample_2["C"]**2 + sample_2["G"]**2 + sample_2["T"]**2)**0.5
        cosine_similarity = dot_product / (magnitude_1 * magnitude_2)

        return pl.when(cosine_similarity >= threshold).then(1).otherwise(0)

    def __getattribute__(self, name):
        if name.startswith("cosani_"):
            try:
                threshold = float(name.split("_")[1])
            except ValueError:
                raise AttributeError(f"Invalid threshold in method name: {name}")
            return lambda: self.generalized_cos_ani(threshold)
        else:
            return super().__getattribute__(name)

def _join_stb_requested_genomes(
    genome_comp: pl.LazyFrame,
    stb_file: Union[str, Path],
    genome_scope: str,
) -> pl.LazyFrame:
    genomes_utf8 = (
        pl.scan_csv(
            stb_file, separator="\t", has_header=False, new_columns=["scaffold", "genome"]
        )
        .select(pl.col("genome").cast(pl.Utf8).str.strip_chars())
        .unique()
    )
    genome_dtype = genome_comp.collect_schema().get("genome")
    if genome_dtype == pl.Categorical:
        categories = sorted(
            set(
                genomes_utf8.select("genome")
                .collect(engine="streaming")["genome"]
                .to_list()
            )
        )
        enum_dtype = pl.Enum(categories)
        genomes = genomes_utf8.with_columns(pl.col("genome").cast(enum_dtype))
        genome_comp = genome_comp.with_columns(pl.col("genome").cast(enum_dtype))
    else:
        genomes = genomes_utf8.with_columns(pl.col("genome").cast(genome_dtype))
    if genome_scope != "all":
        genomes = genomes.filter(pl.col("genome") == genome_scope)
    return genomes.join(genome_comp, on="genome", how="left")


def _as_lazy_profile(source: Union[str, Path, pl.LazyFrame]) -> pl.LazyFrame:
    if isinstance(source, pl.LazyFrame):
        frame = source
    if isinstance(source, (str, Path)):
        frame = pl.scan_parquet(source)
    elif not isinstance(source, pl.LazyFrame):
        raise TypeError(f"Unsupported source type for profile input: {type(source)}")
    return frame


def profile_is_gene_free(source: Union[str, Path, pl.LazyFrame]) -> bool:
    """Return whether a profile is V2 shaped, with genes kept in a shared catalog."""
    return "gene" not in _as_lazy_profile(source).collect_schema().names()


def require_same_profile_generation(sources: Iterable[Union[str, Path, pl.LazyFrame]]) -> None:
    """Reject mixed 1.x and V2 inputs, whose gene results depend on argument order."""
    if len({profile_is_gene_free(source) for source in sources}) > 1:
        raise ValueError(
            "Cannot compare 1.x profiles with V2 profiles. Profile all samples with the same version."
        )


def _profile_scope_predicate(genome_scope: str = "all", gene_scope: str = "all") -> pl.Expr:
    expr = pl.lit(True)
    if genome_scope != "all":
        expr = expr & (pl.col("genome") == genome_scope)
    if gene_scope != "all":
        expr = expr & (pl.col("gene") == gene_scope)
    return expr


def _validate_unmarked_profile(frame: pl.LazyFrame) -> list[dict]:
    """Check coordinate and count invariants for input without a profile order marker."""
    invalid_counts = pl.any_horizontal(pl.col(base).is_null() | (pl.col(base) < 0) for base in "ACGT")
    summary = frame.select(
        (pl.any_horizontal(pl.col("chrom", "genome", "pos").is_null()) | (pl.col("pos") < 1) | invalid_counts)
        .any().alias("invalid"),
        pl.struct("chrom", "pos").is_duplicated().any().alias("duplicate"),
        pl.struct("chrom", "genome").unique().implode().alias("scaffolds"),
    ).collect(engine="streaming").row(0, named=True)
    if summary["invalid"]:
        raise ValueError("Profiles require nonnegative counts and valid coordinates")
    if summary["duplicate"]:
        raise ValueError("Duplicate profile coordinates")
    return summary["scaffolds"]


def _join_position_frames(frames: list[pl.LazyFrame], *, sorted_inputs: bool) -> pl.LazyFrame:
    """Pass the profiler's composite sort order directly to Polars."""
    left, right = frames
    keys = ["chrom", "pos"]
    if sorted_inputs:
        return left.set_sorted(keys).join(right.set_sorted(keys), on=keys, maintain_order="left", suffix="_2")
    return left.join(right, on=keys, suffix="_2")


def _profiles_marked_sorted(*sources) -> bool:
    """Return whether every source is a profile file carrying the profiler's sort marker."""
    return all(
        isinstance(source, (str, Path))
        and _read_custom_parquet_metadata(source).get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE
        for source in sources
    )


def _shared_loci_polars(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
) -> pl.LazyFrame:
    ani_expr = getattr(PolarsANIExpressions(), ani_method)()
    con_expr = PolarsANIExpressions().conani()
    require_same_profile_generation((mpile1, mpile2))
    inputs = [_as_lazy_profile(source) for source in (mpile1, mpile2)]
    scope = _profile_scope_predicate(genome_scope, gene_scope)
    coverage = pl.sum_horizontal("A", "C", "G", "T") >= min_cov
    frames = [frame.filter(scope & coverage) for frame in inputs]
    joined = _join_position_frames(frames, sorted_inputs=_profiles_marked_sorted(mpile1, mpile2))
    gene = [] if profile_is_gene_free(mpile1) else ["gene"]
    return (
        joined
        .with_columns(ani_expr.alias("surr"), con_expr.alias("con_surr"))
        .select(
            pl.col("surr"),
            pl.col("con_surr"),
            scaffold=pl.col("chrom"),
            pos=pl.col("pos"),
            genome=pl.col("genome"),
            *gene,
        )
    )


def polars_prefilter_by_scope(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Scope-filter both profiles while retaining their lazy Parquet scans."""
    if genome_scope == gene_scope == "all":
        return _as_lazy_profile(mpile1), _as_lazy_profile(mpile2)
    scope = _profile_scope_predicate(genome_scope, gene_scope)

    def filtered(source):
        frame = _as_lazy_profile(source).filter(scope)
        if isinstance(source, (str, Path)) and (
            _read_custom_parquet_metadata(source).get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE
        ):
            frame = frame.set_sorted(["chrom", "pos"])
        return frame

    return filtered(mpile1), filtered(mpile2)


def add_genome_info(mpile_contig:pl.LazyFrame, scaffold_to_genome:pl.LazyFrame) -> pl.LazyFrame:
    """
    Adds genome information to the mpileup LazyFrame based on scaffold to genome mapping.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
        scaffold_to_genome (pl.LazyFrame): The LazyFrame mapping scaffolds to genomes.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with genome information added.
    """
    return mpile_contig.join(
        scaffold_to_genome, on="scaffold", how="left"
    ).fill_null("NA")

def get_unique_scaffolds(mpile_contig:pl.LazyFrame,batch_size:int=10000) -> set:
    """
    Retrieves unique scaffolds from the mpileup LazyFrame.

    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
        batch_size (int): The number of rows to process in each batch. Default is 10000.
    Returns:
        set: A set of unique scaffold names.
    """
    scaffolds = set()
    start_index = 0
    while True:
        batch = mpile_contig.slice(start_index, batch_size).select("chrom").collect()
        if batch.height == 0:
            break
        scaffolds.update(batch["chrom"].to_list())
        start_index += batch_size
    return scaffolds 


def compare_genomes_polars(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = options.COMPARE_MIN_COV,
    min_gene_compare_len: int = options.MIN_GENE_COMPARE_LEN,
    genome_scope: str = "all",
    ani_method: str = options.ANI_METHOD,
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
) -> pl.LazyFrame:
    """Compare two profiles fully in Polars and return genome-level statistics.

    V2 profiles use :func:`compare_profiles`, so they get the same values on every
    route. ``min_gene_compare_len`` has no effect here because genome rows do not
    depend on genes. It remains for callers that pass one setting to both routes.
    """
    require_same_profile_generation((mpile_contig_1, mpile_contig_2))
    if profile_is_gene_free(mpile_contig_1):
        selected = options.resolve_metrics(calculate, genes_available=False)
        genome_comp = compare_profiles(
            mpile_contig_1, mpile_contig_2, calculate=selected, ani_method=ani_method, min_cov=min_cov,
        )
        if genome_scope != "all":
            genome_comp = genome_comp.filter(pl.col("genome") == genome_scope)
        if stb_file is not None:
            genome_comp = _join_stb_requested_genomes(genome_comp, stb_file, genome_scope)
            counts = [name for name in ("total_positions", "share_allele_pos", "max_identical_stretch")
                      if name in genome_comp.collect_schema().names()]
            genome_comp = genome_comp.with_columns(pl.col(counts).fill_null(0))
        return genome_comp.select(genome_metric_output_columns(selected))
    calculations = parse_genome_calculations(calculate)
    shared = _shared_loci_polars(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        ani_method=ani_method,
    ).select(pl.col("scaffold").alias("chrom"), pl.col("genome").cast(pl.String), "pos",
             (pl.col("surr") > 0).alias("match"))
    if not _profiles_marked_sorted(mpile_contig_1, mpile_contig_2):
        shared = shared.sort("chrom", "pos")
    wants_stretch = "max_identical_stretch" in calculations

    def counts() -> pl.DataFrame:
        # The V2 streaming scan, keeping the 1.x stretch definition (a run opened by a mismatch counts it).
        genome_counts, _ = _scan_shared_positions(shared, None, stretches=wants_stretch,
                                                  count_opening_mismatch=True)
        return genome_counts.with_columns(
            (pl.col("share_allele_pos") / pl.col("total_positions") * 100).alias("genome_ani"))

    genome_comp = pl.defer(counts, schema={**_GENOME_COUNT_SCHEMA, "genome_ani": pl.Float64}).select(
        "genome",
        *(["total_positions", "share_allele_pos", "genome_ani"] if "genome_ani" in calculations else []),
        *(["max_identical_stretch"] if wants_stretch else []),
    )

    if stb_file is not None:
        genome_comp = _join_stb_requested_genomes(
            genome_comp=genome_comp,
            stb_file=stb_file,
            genome_scope=genome_scope,
        )

    casts: list[pl.Expr] = []
    if "genome_ani" in calculations:
        casts.extend(
            [
                pl.col("total_positions").fill_null(0).cast(pl.Int64),
                pl.col("share_allele_pos").fill_null(0).cast(pl.Int64),
                pl.col("genome_ani").fill_null(0.0).cast(pl.Float64),
            ]
        )
    if "max_identical_stretch" in calculations:
        casts.append(pl.col("max_identical_stretch").fill_null(0).cast(pl.Int64))
    genome_comp = genome_comp.with_columns(casts)

    return genome_comp.select(genome_metric_output_columns(calculations))


def compare_genes_polars(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = options.COMPARE_MIN_COV,
    min_gene_compare_len: int = options.MIN_GENE_COMPARE_LEN,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = options.ANI_METHOD,
) -> pl.LazyFrame:
    """Compare two 1.x profiles fully in Polars and return gene-level ANI statistics."""
    require_same_profile_generation((mpile_contig_1, mpile_contig_2))
    if profile_is_gene_free(mpile_contig_1):
        raise ValueError("V2 gene comparisons require --preprofile-dir (compare_profiles)")
    shared = _shared_loci_polars(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
        ani_method=ani_method,
    )
    return (
        shared.select(["genome", "gene", "surr"])
        .group_by(["genome", "gene"])
        .agg(
            total_positions=pl.len(),
            share_allele_pos=(pl.col("surr") > 0).sum(),
        )
        .filter(pl.col("total_positions") >= min_gene_compare_len)
        .with_columns(
            ani=pl.col("share_allele_pos") / pl.col("total_positions") * 100,
        )
    )


@dataclass
class _AlignmentInput:
    """One profile's checked frame and identity, reusable across every pair it joins."""

    frame: pl.LazyFrame
    marked: bool
    reference_id: Optional[str] = None
    filters: Optional[str] = None
    scaffolds: Optional[list[dict]] = None


def _alignment_input(source) -> _AlignmentInput:
    """Check one profile's column types and read the identity used to align it."""
    from zipstrain import gene_info as gi
    profile = _as_lazy_profile(source)
    schema = profile.collect_schema()
    if any(not schema[column].is_integer() for column in ("pos", *"ACGT")):
        raise ValueError("Profile counts and positions must be integer columns")
    frame = profile.select(pl.col("chrom", "genome").cast(pl.String), pl.col("pos", *"ACGT").cast(pl.Int64))
    metadata = pl.read_parquet_metadata(source) if isinstance(source, (str, Path)) else {}
    return _AlignmentInput(
        frame, metadata.get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE,
        metadata.get(gi.REFERENCE_ID), metadata.get("zipstrain_profile_filters"),
    )


def _scaffold_assignments(item: _AlignmentInput) -> list[dict]:
    """Return a profile's scaffold-to-genome rows, validating an unmarked profile once."""
    if item.scaffolds is None:
        item.scaffolds = (
            item.frame.select("chrom", "genome").unique().collect(engine="streaming").to_dicts()
            if item.marked else _validate_unmarked_profile(item.frame)
        )
    return item.scaffolds


def _align_inputs(first: _AlignmentInput, second: _AlignmentInput, min_cov: int) -> pl.LazyFrame:
    """Inner-join two prepared profiles on shared coordinates after depth filtering."""
    if any((first.reference_id, second.reference_id)) and first.reference_id != second.reference_id:
        raise ValueError("Profiles use incompatible reference identities")
    if any((first.filters, second.filters)) and first.filters != second.filters:
        raise ValueError("Profiles use incompatible profiling filters")
    sorted_inputs = first.marked and second.marked
    if not sorted_inputs:
        assignments = {}
        for item in (first, second):
            for row in _scaffold_assignments(item):
                if assignments.setdefault(row["chrom"], row["genome"]) != row["genome"]:
                    raise ValueError("Conflicting genome assignments between profiles")
    coverage = pl.sum_horizontal("A", "C", "G", "T") >= min_cov
    joined = _join_position_frames(
        [item.frame.filter(coverage) for item in (first, second)], sorted_inputs=sorted_inputs,
    ).drop("genome_2")
    return joined if sorted_inputs else joined.sort(["chrom", "pos"])


def align_profiles(profile_1, profile_2, *, min_cov=options.COMPARE_MIN_COV) -> pl.LazyFrame:
    """Validate and align two position profiles, retaining both samples' counts.

    Apply retained-depth filtering independently and inner-join shared coordinates
    in scaffold and position order. Check reference and profiling-filter metadata.
    The marked path trusts the order and row invariants established by profiling.
    The fallback validates unmarked counts and coordinates and checks scaffold
    assignments before joining. No gene labels are required.
    See ``docs/v2.md`` for the exact input/output contract.
    This is :attr:`PairComparison.aligned` for a pair used once.
    """
    return PairComparison(profile_1, profile_2, min_cov=min_cov).aligned


_STRETCH_SCHEMA = {"genome": pl.String, "chrom": pl.String, "stretch_id": pl.Int64, "start": pl.Int64,
                   "end": pl.Int64, "length": pl.Int64}
_STRETCH_BACKGROUND_SCHEMA = {"genome": pl.String, "_positions": pl.Int64, "_mismatches": pl.Int64,
                              "_segments": pl.Int64}


def _scan_identical_stretches(matches: pl.LazyFrame, min_length: int = 1, *,
                              chunk_size: int = 1 << 20) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Find identical stretches and each genome's background in one streaming NumPy pass.

    ``matches`` holds ``chrom``, ``genome``, ``pos`` and boolean ``match`` in scaffold and
    position order. A stretch is a run of consecutive matching positions and breaks at a
    mismatch, a coordinate gap, or a new scaffold. Return the stretches of at least
    ``min_length`` (``_STRETCH_SCHEMA``, in stream order) and, per genome, the shared positions,
    mismatches and unbroken segments (scaffold starts plus coordinate gaps), which together
    give the number of runs the significance test expects. ``stretch_id`` numbers every
    stretch in stream order from 1, before the length filter. Only run boundaries are kept,
    so memory follows the number of reported stretches, not the number of positions.
    """
    found: list[tuple] = []
    background: dict[str, list[int]] = {}
    pending = None  # [genome, chrom, stretch_id, start, end, length] of a run reaching a batch end
    next_id, previous_chrom, previous_pos, previous_match = 1, None, 0, False

    def emit(run):
        if run[5] >= min_length:
            found.append(tuple(run))

    frames = matches.select("chrom", "genome", "pos", "match")
    for batch in frames.collect_batches(engine="streaming", chunk_size=chunk_size):
        if batch.is_empty():
            continue
        positions, flags, scaffolds = batch["pos"].to_numpy(), batch["match"].to_numpy(), batch["chrom"]
        for first, last in _scaffold_runs(scaffolds):
            chrom, genome = scaffolds[first], batch["genome"][first]
            pos, match = positions[first:last], flags[first:last]
            continues = chrom == previous_chrom and int(pos[0]) == previous_pos + 1
            adjacent = np.diff(pos) == 1
            counts = background.setdefault(genome, [0, 0, 0])
            counts[0] += pos.size
            counts[1] += int(pos.size - match.sum())
            counts[2] += int(not continues) + int((~adjacent).sum())
            # A run starts at a match whose predecessor is a mismatch or not adjacent, and
            # ends at a match whose successor is; the last row of a slice always closes one.
            opened = match & ~np.r_[continues and previous_match, adjacent & match[:-1]]
            closed = match & ~np.r_[adjacent & match[1:], False]
            merge = bool(continues and previous_match and match[0])
            if merge:
                opened[0] = True
            run_starts, run_ends = np.flatnonzero(opened), np.flatnonzero(closed)
            if pending is not None and not merge:
                emit(pending)
                pending = None
            fresh = len(run_starts) - merge
            ids = np.arange(next_id - merge, next_id + fresh)
            next_id += fresh
            lengths = run_ends - run_starts + 1
            starts, ends = pos[run_starts], pos[run_ends]
            for k in range(len(run_starts)):
                run = [genome, chrom, int(ids[k]), int(starts[k]), int(ends[k]), int(lengths[k])]
                if k == 0 and merge:
                    run[2], run[3], run[5] = pending[2], pending[3], pending[5] + run[5]
                    pending = None
                if k == len(run_starts) - 1 and run_ends[k] == pos.size - 1:
                    pending = run
                else:
                    emit(run)
            previous_chrom, previous_pos, previous_match = chrom, int(pos[-1]), bool(match[-1])
    if pending is not None:
        emit(pending)
    stretches = pl.DataFrame(found, schema=_STRETCH_SCHEMA, orient="row")
    totals = pl.DataFrame([(name, *values) for name, values in background.items()],
                          schema=_STRETCH_BACKGROUND_SCHEMA, orient="row")
    return stretches, totals


def _stretch_significance(stretches: pl.LazyFrame, background: pl.LazyFrame) -> pl.LazyFrame:
    """Add how surprising each stretch is if mismatches fall at random at the genome's popANI rate.

    With ``M`` mismatches among ``n`` shared positions, ``d = M / n``. Every mismatch, scaffold
    start and coordinate gap begins a run, so a genome with ``S`` such segments holds
    ``M + S`` runs of matches. A run reaches ``length`` with probability
    ``q = (1 - d) ** length``, so ``expected_runs = (M + S) q`` (an E-value) and ``p_value``
    is the chance that any run does. ``background`` has ``genome``, ``_positions``,
    ``_mismatches`` and ``_segments``, as :func:`_scan_identical_stretches` returns.
    """
    rate = pl.col("_mismatches") / pl.col("_positions")
    runs = pl.col("_mismatches") + pl.col("_segments")
    q = (pl.col("length") * (-rate).log1p()).exp()
    # log1p(-q) stays accurate for tiny q; 1 - exp(x) loses digits near 0, where -x - x²/2 does not.
    x = runs * (-q).log1p()
    p_value = pl.when(x > -1e-3).then(-x - x * x / 2).otherwise(1 - x.exp())
    return stretches.join(background, on="genome", how="left").with_columns(
        rate.alias("background_divergence"), (runs * q).alias("expected_runs"), p_value.alias("p_value"),
    ).drop("_positions", "_mismatches", "_segments")


_MOSAIC_RATE_FLOOR = 1e-12
_MOSAIC_STATES = ("low_divergence", "background", "high_divergence")
_MOSAIC_SCHEMA = {
    "genome": pl.String, "chrom": pl.String, "state": pl.String, "start": pl.Int64, "end": pl.Int64,
    "callable_positions": pl.Int64, "mismatches": pl.Int64, "divergence": pl.Float64,
    "background_divergence": pl.Float64, "log_likelihood_ratio": pl.Float64,
}


def _mosaic_events(matches: pl.LazyFrame) -> pl.LazyFrame:
    """Return each mismatch and each scaffold's start and end, ranked along shared positions.

    ``rank`` counts shared positions in stream order, so the number of matches between two
    events is their rank difference, whatever coordinate gaps lie between. ``kind`` is 0 for a
    scaffold start (ranked just before its first position), 1 for a mismatch and 2 for a
    scaffold end. The result has one row per mismatch and two per scaffold.
    """
    ranked = matches.with_row_index("rank").with_columns(pl.col("rank").cast(pl.Int64),
                                                         pl.col("genome").cast(pl.String), pl.col("chrom").cast(pl.String))
    scaffolds = ranked.group_by("genome", "chrom").agg(
        pl.col("rank").min().alias("first_rank"), pl.col("rank").max().alias("last_rank"),
        pl.col("pos").min().alias("first_pos"), pl.col("pos").max().alias("last_pos"))
    columns = ["genome", "chrom", "rank", "pos", "kind"]
    return pl.concat([
        scaffolds.select("genome", "chrom", (pl.col("first_rank") - 1).alias("rank"),
                         pl.col("first_pos").cast(pl.Int64).alias("pos"), pl.lit(0, pl.Int8).alias("kind")),
        ranked.filter(~pl.col("match")).select("genome", "chrom", "rank", pl.col("pos").cast(pl.Int64),
                                               pl.lit(1, pl.Int8).alias("kind")),
        scaffolds.select("genome", "chrom", pl.col("last_rank").alias("rank"),
                         pl.col("last_pos").cast(pl.Int64).alias("pos"), pl.lit(2, pl.Int8).alias("kind")),
    ]).select(columns)


def _viterbi(emissions: np.ndarray, scaffold_starts: np.ndarray, switch_probability: float) -> np.ndarray:
    """Return the most likely state of each segment, restarting at every scaffold.

    ``emissions`` is ``(segments, 3)`` log-likelihoods. A segment keeps the previous state with
    probability ``1 - switch_probability`` and moves to each other state with half the rest.
    Each scaffold starts in the background state with that same probability.
    """
    stay, move = float(np.log1p(-switch_probability)), float(np.log(switch_probability / 2))
    count = len(emissions)
    ends = np.r_[scaffold_starts[1:], count] - 1
    is_start, is_end = np.zeros(count, dtype=bool), np.zeros(count, dtype=bool)
    is_start[scaffold_starts], is_end[ends] = True, True
    argmax = lambda score: 0 if score[0] >= score[1] and score[0] >= score[2] else (1 if score[1] >= score[2] else 2)
    came_from, final, score = [None] * count, {}, [0.0, 0.0, 0.0]
    # Plain Python floats: three states make NumPy's per-call overhead the dominant cost.
    for index, (row, first, last) in enumerate(zip(emissions.tolist(), is_start.tolist(), is_end.tolist())):
        if first:
            score = [move + row[0], stay + row[1], move + row[2]]
        else:
            best = argmax(score)
            top = score[best] + move
            back, new = [], []
            for state in range(3):
                kept = score[state] + stay
                back.append(state if kept >= top else best)
                new.append(max(kept, top) + row[state])
            came_from[index], score = back, new
        if last:
            final[index] = argmax(score)
    path = np.empty(count, dtype=np.int8)
    for first, last in zip(scaffold_starts.tolist(), ends.tolist()):
        state = final[last]
        for index in range(last, first - 1, -1):
            path[index] = state
            if index > first:
                state = came_from[index][state]
    return path


def _decode_mosaics(events: pl.DataFrame, *, factor: float, switch_probability: float,
                    max_iterations: int = 10) -> pl.DataFrame:
    """Segment one genome's scaffolds into low-divergence, background and high-divergence blocks.

    Each segment is the run of matches up to and including the next mismatch, or up to a
    scaffold end. In a state with mismatch rate ``r``, a segment of ``g`` matches closed by a
    mismatch has likelihood ``(1 - r) ** g * r``, and a final segment ``(1 - r) ** g``. The
    background rate starts at the genome's popANI rate and is refitted from background segments,
    so large mosaics do not set their own background. The other states have ``factor`` times
    fewer or more expected mutations: ``1 - (1 - d) ** (1 / factor)`` and ``1 - (1 - d) ** factor``.
    Return the low and high blocks, merged across consecutive segments in one state.
    """
    empty = pl.DataFrame(schema=_MOSAIC_SCHEMA)
    events = events.sort("chrom", "rank", "kind")
    chrom, rank = events["chrom"].to_numpy(), events["rank"].to_numpy()
    pos, kind = events["pos"].to_numpy(), events["kind"].to_numpy()
    closes = np.flatnonzero(kind != 0)
    opens = closes - 1
    hit = (kind[closes] == 1).astype(np.int64)
    gaps = rank[closes] - rank[opens] - hit
    positions = gaps + hit
    if not hit.sum():
        return empty
    starts = np.where(kind[opens] == 0, pos[opens], pos[opens] + 1)
    ends = pos[closes]
    scaffold_starts = np.flatnonzero(kind[opens] == 0)
    rate = hit.sum() / positions.sum()
    def emission_scores(rate):
        # Clip so a genome at (or refitted to) rate 0 or 1 gives finite scores, not 0 * -inf.
        rates = np.clip([1 - (1 - rate) ** (1 / factor), rate, 1 - (1 - rate) ** factor],
                        _MOSAIC_RATE_FLOOR, 1 - _MOSAIC_RATE_FLOOR)
        return gaps[:, None] * np.log1p(-rates) + hit[:, None] * np.log(rates)

    for _ in range(max_iterations):
        emissions = emission_scores(rate)
        path = _viterbi(emissions, scaffold_starts, switch_probability)
        background = path == 1
        if not hit[background].sum():
            break
        refitted = hit[background].sum() / positions[background].sum()
        if abs(refitted - rate) <= 1e-3 * rate:
            rate = refitted
            break
        rate = refitted
    emissions = emission_scores(rate)
    path = _viterbi(emissions, scaffold_starts, switch_probability)
    new_block = np.r_[True, (path[1:] != path[:-1]) | (kind[opens][1:] == 0)]
    block = np.cumsum(new_block) - 1
    segments = pl.DataFrame({
        "block": block, "chrom": chrom[closes], "state": path, "start": starts, "end": ends,
        "positions": positions, "hit": hit,
        "gain": emissions[np.arange(len(path)), path] - emissions[:, 1],
    }).filter(pl.col("state") != 1)
    if segments.is_empty():
        return empty
    genome = events["genome"][0]
    return segments.group_by("block", maintain_order=True).agg(
        pl.col("chrom").first(), pl.col("state").first(), pl.col("start").first(), pl.col("end").last(),
        pl.col("positions").sum().alias("callable_positions"), pl.col("hit").sum().alias("mismatches"),
        pl.col("gain").sum().alias("log_likelihood_ratio"),
    ).select(
        pl.lit(genome).alias("genome"), "chrom",
        pl.col("state").replace_strict(dict(enumerate(_MOSAIC_STATES)), return_dtype=pl.String),
        "start", "end", "callable_positions", "mismatches",
        (pl.col("mismatches") / pl.col("callable_positions")).alias("divergence"),
        pl.lit(rate).alias("background_divergence"), "log_likelihood_ratio",
    ).cast(_MOSAIC_SCHEMA)


def _scan_shared_positions(shared: pl.LazyFrame, genes: Optional[pl.DataFrame], *, stretches: bool,
                           chunk_size: int = 1 << 20, gene_index: Optional[dict] = None,
                           count_opening_mismatch: bool = False,
                           ) -> tuple[pl.DataFrame, Optional[pl.DataFrame]]:
    """Calculate genome and gene counts from one streaming pass over aligned positions.

    ``shared`` holds ``chrom``, ``genome``, ``pos`` and boolean ``match`` columns, with rows
    grouped by scaffold and in position order. Return one frame per genome with
    ``total_positions``, ``share_allele_pos`` and ``max_identical_stretch``, and, when
    ``genes`` is given, one per gene with ``gene_total_positions`` and
    ``gene_share_allele_pos``. Overlapping genes each count all their own positions.
    A stretch breaks at a mismatch, a coordinate gap, or a new scaffold. 1.x profiles keep
    their definition with ``count_opening_mismatch``, which counts the mismatch that opens a run.
    """
    gene_slices, gene_total, gene_share = {}, None, None
    if genes is not None:
        gene_total = np.zeros(genes.height, dtype=np.int64)
        gene_share = np.zeros(genes.height, dtype=np.int64)
        gene_slices = _index_genes(genes) if gene_index is None else gene_index
    genome_counts: dict[str, list[int]] = {}
    previous_chrom, previous_pos, previous_run = None, 0, 0
    for batch in shared.collect_batches(engine="streaming", chunk_size=chunk_size):
        if batch.is_empty():
            continue
        positions = batch["pos"].to_numpy()
        matches = batch["match"].to_numpy()
        scaffolds = batch["chrom"]
        for first, last in _scaffold_runs(scaffolds):
            chrom, genome = scaffolds[first], batch["genome"][first]
            pos, match = positions[first:last], matches[first:last]
            counts = genome_counts.setdefault(genome, [0, 0, 0])
            counts[0] += pos.size
            counts[1] += int(match.sum())
            if stretches:
                edges = np.r_[np.flatnonzero(~match | np.r_[False, np.diff(pos) != 1]), pos.size]
                if edges[0] != 0:
                    edges = np.r_[0, edges]
                # A run opened by a mismatch does not count that position.
                runs = np.diff(edges) - (0 if count_opening_mismatch else ~match[edges[:-1]])
                if chrom == previous_chrom and pos[0] == previous_pos + 1 and match[0]:
                    runs[0] += previous_run
                counts[2] = max(counts[2], int(runs.max()))
                previous_chrom, previous_pos, previous_run = chrom, int(pos[-1]), int(runs[-1])
            if chrom in gene_slices:
                index, gene_starts, gene_ends = gene_slices[chrom]
                prefix = np.r_[0, np.cumsum(match, dtype=np.int64)]
                low = np.searchsorted(pos, gene_starts, side="left")
                high = np.searchsorted(pos, gene_ends, side="right")
                gene_total[index] += high - low
                gene_share[index] += prefix[high] - prefix[low]
    genome_frame = pl.DataFrame(
        [(name, *counts) for name, counts in genome_counts.items()],
        schema={"genome": pl.String, "total_positions": pl.Int64, "share_allele_pos": pl.Int64,
                "max_identical_stretch": pl.Int64},
        orient="row",
    )
    if genes is None:
        return genome_frame, None
    return genome_frame, genes.select("gene_id", "genome").with_columns(
        pl.Series("gene_total_positions", gene_total), pl.Series("gene_share_allele_pos", gene_share),
    )


def compare_profiles(profile_1, profile_2, *, preprofile_dir=None, calculate=options.CALCULATE,
                     ani_method=options.ANI_METHOD, min_cov=options.COMPARE_MIN_COV,
                     min_gene_compare_len=options.MIN_GENE_COMPARE_LEN,
                     codon_counts_1=None, codon_counts_2=None,
                     sample_1=None, sample_2=None, allele_integration=options.ALLELE_INTEGRATION) -> pl.LazyFrame:
    """Compare one pair using ordinary Polars functions and explicit shared assets.

    Return one row per genome, or one row per qualifying gene if any gene
    metric was requested. Genome metrics repeat across genes in combined runs.
    ``calculate`` follows :func:`zipstrain.options.resolve_metrics`.
    Gene metrics require an authenticated ``preprofile_dir/gene_info.parquet``;
    overlapping genes are aggregated independently from coordinate intervals.
    ``gene_selection`` additionally requires both observed codon tables.
    SNS/SNV integration follows the gene-overlap method (consensus or weighted). Unknown/uncallable ANI and rates are
    null, not zero. Input count filters, profile identities, and annotation IDs
    are validated before using gene IDs. See ``docs/v2.md`` for examples.
    This is :meth:`PairComparison.compare` for a pair used once.
    """
    if min_gene_compare_len < 1:
        raise ValueError("min_gene_compare_len must be positive")
    pair = PairComparison(
        profile_1, profile_2, preprofile_dir=preprofile_dir, ani_method=ani_method, min_cov=min_cov,
        sample_1=sample_1, sample_2=sample_2, codon_counts_1=codon_counts_1, codon_counts_2=codon_counts_2,
        allele_integration=allele_integration,
    )
    return pair.compare(calculate, min_gene_compare_len=min_gene_compare_len)


def _file_digest(path) -> str:
    """Return a file's SHA-256, hashing each unchanged file version once per process."""
    path = Path(path).resolve()
    stat = path.stat()
    return _cached_digest(str(path), stat.st_mtime_ns, stat.st_size)


@lru_cache(maxsize=64)
def _cached_digest(path: str, _mtime_ns: int, _size: int) -> str:
    from zipstrain.utils import sha256_file
    return sha256_file(path)


def _pair_names(profile_1, profile_2, sample_1=None, sample_2=None) -> tuple[str, str]:
    """Name a pair's samples, defaulting to the profile file names without ``_profile``."""
    def name(source, given, default):
        if given:
            return given
        return Path(source).stem.removesuffix("_profile") if isinstance(source, (str, Path)) else default
    return name(profile_1, sample_1, "sample_1"), name(profile_2, sample_2, "sample_2")


class PairComparison:
    """Two validated profiles and their shared-position join, reused by every pair analysis.

    Construction does once what each analysis needs: the profile-generation, reference and
    filter checks, optional preprofile assets, sample names and the ANI method. Every
    analysis is lazy and reads no profile rows until collected. Polars plans built on
    ``aligned`` and ``matches`` share one join and one streaming pass when collected
    together with :meth:`collect`::

        pair = PairComparison("A_profile.parquet", "B_profile.parquet", preprofile_dir="assets")
        positions, blocks = pair.collect(
            pair.matches.group_by("genome").agg(pl.len()),
            pair.mosaic_regions(),
        )

    :meth:`compare` and :meth:`identical_stretches` count in NumPy batches and run their
    own pass when collected. Codon tables are only read when ``compare(calculate=...)``
    with ``gene_selection`` is collected.
    """

    def __init__(self, profile_1, profile_2, *, preprofile_dir=None, ani_method=options.ANI_METHOD,
                 min_cov=options.COMPARE_MIN_COV, sample_1=None, sample_2=None,
                 codon_counts_1=None, codon_counts_2=None, allele_integration=options.ALLELE_INTEGRATION):
        from zipstrain import gene_info as gi, matrix_pairs as mp
        if min_cov < 1:
            raise ValueError("min_cov must be positive")
        self.profiles = (profile_1, profile_2)
        self.preprofile_dir = None if preprofile_dir is None else Path(preprofile_dir)
        self.gene_file, self.genes = None, None
        if self.preprofile_dir is not None:
            self.gene_file = self.preprofile_dir/"gene_info.parquet"
            gi.validate_profile_assets(list(self.profiles), self.gene_file)
            self.genes = gi.read_gene_info(self.gene_file)
        self.names = _pair_names(profile_1, profile_2, sample_1, sample_2)
        self.codon_counts = (codon_counts_1, codon_counts_2)
        self.allele_integration = allele_integration
        self.min_cov = min_cov
        self.ani_method, _threshold = mp.parse_matrix_ani_method(ani_method)
        require_same_profile_generation(self.profiles)
        # Reject incompatible samples before any analysis spends time on them.
        self.aligned: pl.LazyFrame = _align_inputs(
            _alignment_input(profile_1), _alignment_input(profile_2), min_cov)

    @property
    def matches(self) -> pl.LazyFrame:
        """Shared positions as ``chrom``, ``genome``, ``pos`` and whether they match under ``ani_method``."""
        match = getattr(PolarsANIExpressions(), self.ani_method)() > 0
        return self.aligned.select("chrom", "genome", "pos", match.alias("match"))

    @staticmethod
    def collect(*frames: pl.LazyFrame, engine: str = "streaming") -> list[pl.DataFrame]:
        """Collect several analyses together, so Polars runs their shared join once."""
        return pl.collect_all(list(frames), engine=engine)

    def compare(self, calculate=options.CALCULATE, *,
                min_gene_compare_len=options.MIN_GENE_COMPARE_LEN) -> pl.LazyFrame:
        """Return genome or gene comparison rows, as described in :func:`compare_profiles`."""
        from zipstrain import codons
        if min_gene_compare_len < 1:
            raise ValueError("min_gene_compare_len must be positive")
        codon_tables = dict(zip(self.names, self.codon_counts, strict=True))
        calculations = options.resolve_metrics(
            calculate, genes_available=self.genes is not None and self.genes.height > 0,
            codon_tables_present=[name for name, counts in codon_tables.items() if counts is not None],
            codon_tables_missing=[name for name, counts in codon_tables.items() if counts is None],
        )
        rates = None
        if "gene_selection" in calculations:
            def rates():
                # At collect time, so building the plan reads no codons and hashes no profile.
                for counts, positions in zip(self.codon_counts, self.profiles, strict=True):
                    if pl.read_parquet_metadata(counts).get("zipstrain_source_profile_sha256") != _file_digest(positions):
                        raise ValueError("Codon evidence does not match its sample position profile")
                return codons.pairwise_selection(*self.codon_counts, self.gene_file, min_cov=self.min_cov,
                                                 allele_integration=self.allele_integration)
        universe = [_as_lazy_profile(source).select(pl.col("genome").cast(pl.String)).unique()
                    for source in self.profiles]
        return _compare_aligned(
            self.aligned, self.names, calculations=calculations, ani_method=self.ani_method,
            min_gene_compare_len=min_gene_compare_len, genes=self.genes, rates=rates,
            genome_universe=[*universe, *_catalog_genomes(self.preprofile_dir, self.genes)],
        )

    def identical_stretches(self, min_length: int) -> pl.LazyFrame:
        """Return identical stretches of at least ``min_length``, as in :func:`identical_stretches`."""
        if min_length < 1:
            raise ValueError("min_length must be positive")
        matches, names = self.matches, self.names

        def run() -> pl.DataFrame:
            stretches, background = _scan_identical_stretches(matches, min_length)
            return _stretch_significance(stretches.lazy(), background.lazy()).select(
                pl.lit(names[0]).alias("sample_1"), pl.lit(names[1]).alias("sample_2"), *_STRETCH_SCHEMA,
                "background_divergence", "expected_runs", "p_value").collect()

        schema = {"sample_1": pl.String, "sample_2": pl.String, **_STRETCH_SCHEMA,
                  "background_divergence": pl.Float64, "expected_runs": pl.Float64, "p_value": pl.Float64}
        # Its own NumPy pass when collected, like :meth:`compare`.
        return pl.defer(run, schema=schema)

    def mosaic_regions(self, *, factor: float = 20.0, switch_probability: float = 1e-3) -> pl.LazyFrame:
        """Return blocks whose mismatches are much rarer or much denser than the genome's background.

        Under the background model, popANI mismatches fall at random along each genome's shared
        positions at one rate. A three-state hidden Markov model over the runs of matches between
        mismatches finds blocks that fit ``factor`` times fewer (``low_divergence``, such as recent
        recombination) or ``factor`` times more (``high_divergence``, such as DNA imported from a
        more distant lineage) expected mutations better than the background. The background
        rate starts at the genome's popANI and is refitted without the blocks.
        ``switch_probability`` is the prior chance that the state changes at a mismatch.

        Return one row per block: ``sample_1``, ``sample_2``, ``genome``, ``chrom``, ``state``,
        ``start``, ``end``, ``callable_positions``, ``mismatches``, ``divergence``,
        ``background_divergence`` and ``log_likelihood_ratio`` (the block's log-likelihood in its
        state minus in the background state; larger is stronger). A block runs from the first
        position after a mismatch, or the scaffold's first shared position, to the mismatch that
        ends it, or the scaffold's last shared position. Lengths count shared positions, so
        positions not covered in both samples are skipped rather than counted as matches.
        Coverage artefacts, such as reads recruited from another species or multi-copy genes,
        also look like ``high_divergence`` blocks, so check depth before reading them as imports.
        """
        from functools import partial
        if factor <= 1:
            raise ValueError("factor must be greater than 1")
        if not 0 < switch_probability < 1:
            raise ValueError("switch_probability must be between 0 and 1")
        decode = partial(_decode_mosaics, factor=factor, switch_probability=switch_probability)
        return (
            _mosaic_events(self.matches)
            .group_by("genome").map_groups(decode, schema=_MOSAIC_SCHEMA)
            .sort("genome", "chrom", "start")
            .select(pl.lit(self.names[0]).alias("sample_1"), pl.lit(self.names[1]).alias("sample_2"),
                    *_MOSAIC_SCHEMA)
        )


def identical_stretches(profile_1, profile_2, *, min_length: int, ani_method=options.ANI_METHOD,
                        min_cov=options.COMPARE_MIN_COV, sample_1=None, sample_2=None) -> pl.LazyFrame:
    """Return every identical stretch of at least ``min_length`` positions between two profiles.

    A stretch is a run of consecutive shared positions that all match under ``ani_method``.
    It breaks at a mismatch, at a position not covered by ``min_cov`` in both samples, and at
    a new scaffold, exactly as ``max_identical_stretch`` does, so the longest stretch here
    has that length. Return one row per stretch with ``sample_1``, ``sample_2``, ``genome``,
    ``chrom``, ``stretch_id``, ``start``, ``end``, ``length``, ``background_divergence``,
    ``expected_runs`` and ``p_value``. The last three test the stretch against mismatches
    falling at random at the genome's popANI rate ``d``: with ``M`` mismatches and ``S`` unbroken
    segments (scaffold starts and coordinate gaps) the genome has ``M + S`` runs,
    ``expected_runs = (M + S)(1 - d) ** length`` is the E-value, and
    ``p_value`` is the chance that some run is at least this long. Large mosaics change ``d``
    itself; :meth:`PairComparison.mosaic_regions` refits it. ``start`` and ``end`` are
    one-based, inclusive profile coordinates. ``stretch_id`` numbers every stretch of the
    pair in scaffold and position order before the length filter, so it is unique within
    the pair but has gaps. This uses the table method and needs no preprofile assets.
    Profiles carrying the profiler's sort marker are merge-joined, as in
    :func:`compare_profiles`. Collect with ``engine="streaming"`` for large genomes.
    This is :meth:`PairComparison.identical_stretches` for a pair used once.
    """
    if min_length < 1:
        raise ValueError("min_length must be positive")
    pair = PairComparison(profile_1, profile_2, ani_method=ani_method, min_cov=min_cov,
                          sample_1=sample_1, sample_2=sample_2)
    return pair.identical_stretches(min_length)


def mosaic_regions(profile_1, profile_2, *, factor: float = 20.0, switch_probability: float = 1e-3,
                   ani_method=options.ANI_METHOD, min_cov=options.COMPARE_MIN_COV,
                   sample_1=None, sample_2=None) -> pl.LazyFrame:
    """Return blocks of unusually low or high divergence between two profiles.

    This is :meth:`PairComparison.mosaic_regions` for a pair used once. It uses the table
    method and needs no preprofile assets.
    """
    pair = PairComparison(profile_1, profile_2, ani_method=ani_method, min_cov=min_cov,
                          sample_1=sample_1, sample_2=sample_2)
    return pair.mosaic_regions(factor=factor, switch_probability=switch_probability)


def _catalog_genomes(preprofile_dir, genes: Optional[pl.DataFrame]) -> list[pl.LazyFrame]:
    """Return the genomes a comparison must report even when no position is shared."""
    lengths = Path(preprofile_dir)/"genome_lengths.parquet" if preprofile_dir is not None else None
    if lengths is not None and lengths.exists():
        return [pl.scan_parquet(lengths).select(pl.col("genome").cast(pl.String)).unique()]
    return [genes.lazy().select(pl.col("genome").cast(pl.String)).unique()] if genes is not None else []


_GENOME_COUNT_SCHEMA = {"genome": pl.String, "total_positions": pl.Int64, "share_allele_pos": pl.Int64,
                        "max_identical_stretch": pl.Int64}
_GENE_COUNT_SCHEMA = {"gene_id": pl.UInt64, "genome": pl.String, "gene_total_positions": pl.Int64,
                      "gene_share_allele_pos": pl.Int64}


def _compare_aligned(aligned: pl.LazyFrame, names: tuple[str, str], *, calculations: tuple[str, ...],
                     ani_method: str, min_gene_compare_len: int, genes: Optional[pl.DataFrame],
                     genome_universe: list, rates=None,
                     gene_index: Optional[dict] = None) -> pl.LazyFrame:
    """Build one pair's comparison rows from already validated inputs and resolved metrics.

    The result is lazy: the NumPy counting pass runs when it is collected. ``rates`` is a
    selection frame, or a function returning one that is called at collect time.
    """
    from zipstrain import codons
    shared = aligned.select(
        "chrom", "genome", "pos", (getattr(PolarsANIExpressions(), ani_method)() > 0).alias("match"))
    gene_rows = "gene_ani" in calculations or rates is not None
    wants_stretch = "max_identical_stretch" in calculations

    def rows(genome_counts: pl.DataFrame, gene_counts: Optional[pl.DataFrame], rate_frame) -> pl.LazyFrame:
        return _comparison_rows(genome_counts, gene_counts, rate_frame, names, calculations=calculations,
                                min_gene_compare_len=min_gene_compare_len, genome_universe=genome_universe)

    def run() -> pl.DataFrame:
        genome_counts, gene_counts = _scan_shared_positions(
            shared, genes if gene_rows else None, stretches=wants_stretch, gene_index=gene_index)
        rate_frame = rates() if callable(rates) else rates
        return rows(genome_counts, gene_counts, rate_frame).collect(engine="streaming")

    empty_rates = pl.LazyFrame(schema=codons._SELECTION_SCHEMA) if rates is not None else None
    schema = rows(pl.DataFrame(schema=_GENOME_COUNT_SCHEMA),
                  pl.DataFrame(schema=_GENE_COUNT_SCHEMA) if gene_rows else None, empty_rates).collect_schema()
    return pl.defer(run, schema=schema)


def _comparison_rows(genome_counts: pl.DataFrame, gene_counts: Optional[pl.DataFrame], rates,
                     names: tuple[str, str], *, calculations: tuple[str, ...], min_gene_compare_len: int,
                     genome_universe: list) -> pl.LazyFrame:
    """Turn one pair's genome and gene counts into its output rows."""
    from zipstrain import gene_info as gi
    keys = ["sample_1", "sample_2", "genome"]
    sample_columns = [pl.lit(name).alias(key) for key, name in zip(keys[:2], names, strict=True)]
    wants_ani, wants_stretch = "genome_ani" in calculations, "max_identical_stretch" in calculations
    genome = None
    if wants_ani or wants_stretch:
        columns = ["total_positions", "share_allele_pos", "max_identical_stretch"]
        genome = pl.concat(genome_universe).unique().join(genome_counts.lazy(), on="genome", how="left").with_columns(
            pl.col(columns).fill_null(0),
        ).with_columns(
            pl.when(pl.col("total_positions") > 0)
            .then(pl.col("share_allele_pos")/pl.col("total_positions")*100).alias("genome_ani"),
        ).select(
            "genome",
            *(["total_positions", "share_allele_pos", "genome_ani"] if wants_ani else []),
            *(["max_identical_stretch"] if wants_stretch else []),
        )
    if gene_counts is not None:
        gene = gene_counts.lazy().filter(pl.col("gene_total_positions") >= min_gene_compare_len)
        if "gene_ani" in calculations:
            gene = gene.with_columns(
                (100 * pl.col("gene_share_allele_pos") / pl.col("gene_total_positions")).alias("gene_ani")
            )
        if rates is not None:
            gene = gene.join(rates.select("gene_id", *gi.PAIR_SELECTION_COLUMNS), on="gene_id", how="left")
        if genome is not None:
            gene = gene.join(genome, on="genome", how="left")
        return gene.with_columns(sample_columns).select(*keys, "gene_id", pl.exclude(*keys, "gene_id"))
    assert genome is not None
    return genome.with_columns(sample_columns).select(*keys, pl.exclude(keys))


def resolve_collection_metrics(profiles, preprofile_dir, calculate, codon_counts=None) -> tuple:
    """Decide once, for a whole set of samples, which metrics every pair will calculate.

    Return the gene catalog, the resolved metrics, and each sample's codon table, given by name
    in ``codon_counts`` or found beside its profile. Scheduled and local runs share this, so
    one sample's missing table can never give pairs with different columns.
    """
    from zipstrain import codons, gene_info as gi
    genes = gi.read_gene_info(Path(preprofile_dir)/"gene_info.parquet") if preprofile_dir is not None else None
    codon_paths = {}
    for name, path in profiles:
        candidate = (codon_counts or {}).get(name) or codons.sibling_codon_table(path)
        if candidate is not None and Path(candidate).exists():
            codon_paths[name] = Path(candidate)
    calculations = options.resolve_metrics(
        calculate, genes_available=genes is not None and genes.height > 0,
        codon_tables_present=sorted(codon_paths),
        codon_tables_missing=[name for name, _ in profiles if name not in codon_paths],
    )
    return genes, calculations, codon_paths


def _pair_workers(requested, pair_count: int, memory_limit_gb: float, largest_rows: int) -> tuple[int, int]:
    """Choose concurrent pair processes and each one's Polars threads from CPU and memory.

    Pairs barely benefit from more threads, but scale with more processes that each have
    few, so ``auto`` uses half the allocated CPUs at two threads each, unless the work is too
    small to repay starting processes. It stays within the memory one pair holds: roughly 128
    bytes per profile row in the join and gene counters.
    """
    from zipstrain.resource_limits import duckdb_thread_budget
    if requested != "auto" and (not isinstance(requested, int) or requested < 1):
        raise ValueError("Pair workers must be a positive integer or 'auto'")
    cpus = duckdb_thread_budget()
    by_memory = max(1, int(memory_limit_gb * 1024**3 // max(1, largest_rows * _PAIR_BYTES_PER_ROW)))
    wanted = requested
    if requested == "auto":
        wanted = max(1, cpus // 2) if pair_count * largest_rows >= _AUTO_WORKERS_MIN_ROWS else 1
    workers = max(1, min(pair_count, wanted, by_memory))
    return workers, max(1, cpus // workers)


@dataclass(frozen=True)
class _PairJob:
    """One pair's inputs and settings, in a form a worker process can receive."""

    names: tuple[str, str]
    paths: tuple[Path, Path]
    target: Path
    identity: str
    metadata: dict
    calculations: tuple[str, ...]
    ani_method: str
    min_cov: int
    min_gene_compare_len: int
    preprofile_dir: Optional[Path]
    rates: Optional[Path]


@lru_cache(maxsize=None)
def _shared_input(path: Path) -> _AlignmentInput:
    """Read each profile's identity once per process, however many pairs it joins."""
    return _alignment_input(path)


@lru_cache(maxsize=None)
def _shared_sample_genomes(path: Path) -> pl.DataFrame:
    """Read each profile's genomes once per process, however many pairs it joins."""
    return _shared_input(path).frame.select("genome").unique().collect(engine="streaming")


@lru_cache(maxsize=1)
def _shared_catalog(preprofile_dir: Optional[Path]) -> tuple:
    """Authenticate and index the gene catalog once per process."""
    from zipstrain import gene_info as gi
    genes = gi.read_gene_info(preprofile_dir/"gene_info.parquet") if preprofile_dir is not None else None
    return (genes, _index_genes(genes) if genes is not None else None,
            [frame.collect() for frame in _catalog_genomes(preprofile_dir, genes)])


def _clear_shared_inputs() -> None:
    """Forget per-process inputs so a later run sees files that have changed."""
    for function in (_shared_input, _shared_sample_genomes, _shared_catalog):
        function.cache_clear()


def _compare_pair_job(job: _PairJob) -> None:
    """Calculate and cache one pair's rows, atomically, in a worker or in this process."""
    genes, gene_index, catalog_genomes = _shared_catalog(job.preprofile_dir)
    gene_rows = any(metric in options.GENE_METRICS for metric in job.calculations)
    genome_rows = any(metric in options.GENOME_METRICS for metric in job.calculations)
    frame = _compare_aligned(
        _align_inputs(*map(_shared_input, job.paths), job.min_cov), job.names,
        calculations=job.calculations, ani_method=job.ani_method,
        min_gene_compare_len=job.min_gene_compare_len, genes=genes if gene_rows else None,
        rates=pl.scan_parquet(job.rates) if job.rates is not None else None,
        gene_index=gene_index if gene_rows else None,
        genome_universe=[*(_shared_sample_genomes(path).lazy() for path in job.paths if genome_rows),
                         *(genomes.lazy() for genomes in catalog_genomes)],
    )
    with TemporaryDirectory(dir=job.target.parent) as staging:
        temporary = Path(staging) / job.target.name
        frame.sink_parquet(temporary, engine="streaming", compression="zstd", metadata={
            **job.metadata, "zipstrain_comparison_id": job.identity,
        })
        temporary.replace(job.target)


def compare_profile_collection(profiles, run_dir, *, preprofile_dir=None,
                               calculate=options.CALCULATE, ani_method=options.ANI_METHOD,
                               min_cov=options.COMPARE_MIN_COV,
                               min_gene_compare_len=options.MIN_GENE_COMPARE_LEN,
                               codon_counts=None, allele_integration=options.ALLELE_INTEGRATION,
                               memory_limit_gb=16.0, workers=1):
    """Run/resume V2 Polars comparisons for named profile paths and export Parquet.

    ``profiles`` is an iterable of ``(sample_name, profile_path)``. Codon files
    may be supplied by sample name in ``codon_counts``; otherwise use sibling
    ``<sample>_codon_counts.parquet`` files. Cache each pair under
    ``intermediate_files/v2_pairs`` with a content-addressed input/settings
    contract in Parquet metadata. Only complete, authenticated pair outputs are
    reused. Added samples compute only new pairs; changed inputs/settings never
    reuse stale results. Final exports are atomically replaced and contain only
    the requested collection's pairs. Selection metrics use the explicitly documented gene-overlap SNS/SNV definitions.

    Inputs are validated and each sample's codons packed once, however many pairs it joins.
    ``workers`` runs that many pairs in separate processes, within the CPU allocation and
    ``memory_limit_gb``; ``"auto"`` chooses the count. Results do not depend on it. More than
    one worker re-imports the calling script, so a script needs ``if __name__ == "__main__"``.
    """
    from zipstrain import codons, gene_info as gi, matrix_pairs as mp, pair_selection, utils
    from zipstrain.resource_limits import job_pool, map_jobs
    if min_cov < 1:
        raise ValueError("min_cov must be positive")
    if min_gene_compare_len < 1:
        raise ValueError("min_gene_compare_len must be positive")
    profiles = sorted((str(name), Path(path)) for name, path in profiles)
    if len(profiles) < 2 or len({name for name, _ in profiles}) != len(profiles):
        raise ValueError("At least two uniquely named samples are required")
    require_same_profile_generation([path for _, path in profiles])
    genes, calculations, codon_paths = resolve_collection_metrics(profiles, preprofile_dir, calculate, codon_counts)
    selection = "gene_selection" in calculations
    gene_rows = any(metric in options.GENE_METRICS for metric in calculations)
    # Genome rows do not depend on the threshold, so it must not invalidate their cache.
    gene_threshold = min_gene_compare_len if gene_rows else None
    if preprofile_dir is not None:
        metadata = gi.validate_profile_assets([p for _, p in profiles], Path(preprofile_dir)/"gene_info.parquet")
    else:
        metadata = {gi.SCHEMA_KEY: "2"}
    metadata.update({"zipstrain_compare_reference_id": metadata.get(gi.REFERENCE_ID, "NA"),
                     "zipstrain_compare_annotation_id": metadata.get(gi.ANNOTATION_ID, "NA"),
                     "zipstrain_compare_min_cov": str(min_cov),
                     "zipstrain_compare_min_gene_compare_len": "NA" if gene_threshold is None else str(gene_threshold),
                     "zipstrain_compare_allele_integration": allele_integration,
                     "zipstrain_compare_selection_method": codons.METHOD if selection else "not_calculated",
                     "zipstrain_compare_kind": "gene" if gene_rows else "genome",
                     utils.COMPARE_ENGINE_METADATA_KEY: "table",
                     "zipstrain_compare_ani_method": ani_method,
                     "zipstrain_compare_calculate": "+".join(calculations),
                     "zipstrain_allele_integration": allele_integration,
                     "zipstrain_selection_method": codons.METHOD if selection else "not_calculated"})
    run_dir = Path(run_dir)
    cache = run_dir/"intermediate_files"/"v2_pairs"
    cache.mkdir(parents=True, exist_ok=True)
    hashes = {name: utils.sha256_file(path) for name, path in profiles}
    codon_hashes = {name: utils.sha256_file(path) for name, path in codon_paths.items()} if selection else {}
    pairs = {}
    for i, (name_1, _) in enumerate(profiles):
        for name_2, _ in profiles[i+1:]:
            contract = {"version": "v2-pair-5", "names": [name_1, name_2],
                "profiles": [hashes[name_1], hashes[name_2]], "calculate": calculations,
                "ani_method": ani_method, "min_cov": min_cov,
                "min_gene_compare_len": gene_threshold, "allele_integration": allele_integration,
                "codons": [codon_hashes.get(name_1), codon_hashes.get(name_2)], **metadata}
            identity = gi.content_digest(contract)
            pairs[name_1, name_2] = identity, cache / f"{identity}.parquet"
    missing = {pair: value for pair, value in pairs.items()
               if utils._read_custom_parquet_metadata(value[1]).get("zipstrain_comparison_id") != value[0]}
    if missing:
        method, _threshold = mp.parse_matrix_ani_method(ani_method)
        directory = Path(preprofile_dir) if preprofile_dir is not None else None
        paths = dict(profiles)
        _clear_shared_inputs()
        # Fail on any unusable profile before starting worker processes.
        for path in paths.values():
            _shared_input(path)
        rows = max(pl.scan_parquet(path).select(pl.len()).collect().item() for path in paths.values())
        worker_count, threads = _pair_workers(workers, len(missing), memory_limit_gb, rows)
        with job_pool(worker_count, threads) as pool:
            rates = {}
            if selection:
                rates = dict(zip(pairs, pair_selection.stage_selection_pairs(
                    profiles, hashes, directory/"gene_info.parquet", cache.parent/"pair_selection",
                    min_cov=min_cov, genome="all", allele_integration=allele_integration,
                    memory_limit_gb=memory_limit_gb, codon_tables=codon_paths, codon_digests=codon_hashes,
                    pool=pool, workers=worker_count,
                ), strict=True))
            map_jobs(pool, _compare_pair_job, [
                _PairJob(pair, (paths[pair[0]], paths[pair[1]]), target, identity, metadata, calculations,
                         method, min_cov, min_gene_compare_len, directory, rates.get(pair))
                for pair, (identity, target) in missing.items()
            ])
        _clear_shared_inputs()
    files = [target for _, target in pairs.values()]
    output = run_dir / "all_comparisons.parquet"
    with TemporaryDirectory(dir=run_dir) as staging:
        temporary = Path(staging) / output.name
        pl.scan_parquet(files).sink_parquet(temporary, engine="streaming", compression="zstd", metadata=metadata)
        temporary.replace(output)
    return output


def compare_genomes(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = options.COMPARE_MIN_COV,
    min_gene_compare_len: int = options.MIN_GENE_COMPARE_LEN,
    genome_scope: str = "all",
    ani_method: str = options.ANI_METHOD,
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
) -> pl.LazyFrame:
    """Compare two profiles at the genome level using Polars."""
    return compare_genomes_polars(
        mpile_contig_1=mpile_contig_1,
        mpile_contig_2=mpile_contig_2,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome_scope,
        ani_method=ani_method,
        stb_file=stb_file,
        calculate=calculate,
    )


def matrix_surr_numpy(
    counts_np: "np.ndarray",
    N: int,
    L: int,
    min_cov: int,
    ani_method: str,
    chunk_size: int,
) -> "tuple[np.ndarray, np.ndarray]":
    """CPU all-pairs match/shared counts via numpy, chunked over positions.

    Parameters
    ----------
    counts_np:
        int32 array of shape (N, L, 4) — samples × positions × bases (A,C,G,T).
    N, L:
        Number of samples and positions.
    min_cov:
        Minimum total coverage to include a position for a sample.
    ani_method:
        ``"popani"``, ``"conani"``, or ``"cosani_<threshold>"``.
    chunk_size:
        Positions processed per iteration. Peak memory is O(N²×chunk_size).

    Returns
    -------
    matches, shared : np.ndarray
        Both float64 arrays of shape (N, N).
    """
    import numpy as np

    cov = counts_np.sum(axis=2)
    mask = (cov >= min_cov)

    matches = np.zeros((N, N), dtype=np.float64)
    shared  = np.zeros((N, N), dtype=np.float64)

    for start in range(0, L, chunk_size):
        c = counts_np[:, start:start + chunk_size, :]
        m = mask[:, start:start + chunk_size]

        shared += m.astype(np.float64) @ m.astype(np.float64).T

        if ani_method == "conani":
            argmax_allele = c.argmax(axis=2)
            H = np.zeros((N, c.shape[1], 4), dtype=np.int8)
            rows, cols = np.nonzero(m)
            H[rows, cols, argmax_allele[rows, cols]] = 1
            H_flat = H.reshape(N, -1).astype(np.int32)
            matches += (H_flat @ H_flat.T).astype(np.float64)
        elif ani_method == "popani":
            c_i64 = c.astype(np.int64)
            surr = np.einsum('ilk,jlk->ijl', c_i64, c_i64)
            mask_pair = m[:, np.newaxis, :] & m[np.newaxis, :, :]
            surr[~mask_pair] = 0
            matches += (surr > 0).sum(axis=2).astype(np.float64)
        elif ani_method.startswith("cosani_"):
            threshold = float(ani_method.split("_", 1)[1])
            c_f64 = c.astype(np.float64)
            norms = np.sqrt((c_f64 ** 2).sum(axis=2, keepdims=True))
            norms = np.where(norms == 0, 1.0, norms)
            normed = c_f64 / norms
            cos_sim = np.einsum('ilk,jlk->ijl', normed, normed)
            mask_pair = m[:, np.newaxis, :] & m[np.newaxis, :, :]
            cos_sim[~mask_pair] = 0.0
            matches += (cos_sim >= threshold).sum(axis=2).astype(np.float64)
        else:
            raise ValueError(f"Unknown ani_method: {ani_method}")

    return matches, shared


def matrix_surr_torch(
    counts_np: "np.ndarray",
    N: int,
    L: int,
    min_cov: int,
    ani_method: str,
    device: "torch.device",
    chunk_size: int,
) -> "tuple[np.ndarray, np.ndarray]":
    """GPU all-pairs match/shared counts via PyTorch, chunked over positions.

    Same semantics as :func:`matrix_surr_numpy` but runs on the given torch
    device (``"cuda"`` or ``"mps"``).  MPS uses float32 throughout since it
    does not support float64/int64.

    Returns
    -------
    matches, shared : np.ndarray
        Both float64 arrays of shape (N, N), moved to CPU.
    """
    try:
        import torch
    except ImportError:
        raise ImportError(
            'PyTorch is required for GPU comparison. Install with: pip install "zipstrain[matrix]". '
            "For CUDA, install the matching Torch wheel from pytorch.org."
        )

    is_mps = (device.type == "mps")

    matches = torch.zeros((N, N), dtype=torch.float32, device=device)
    shared  = torch.zeros((N, N), dtype=torch.float32, device=device)

    counts_t = torch.from_numpy(counts_np).to(device)
    cov = counts_t.sum(dim=2)
    mask = (cov >= min_cov)

    for start in range(0, L, chunk_size):
        c = counts_t[:, start:start + chunk_size, :]
        m = mask[:, start:start + chunk_size]
        m_f = m.float()
        shared += m_f @ m_f.T

        if ani_method == "conani":
            argmax_allele = c.argmax(dim=2)
            C = c.shape[1]
            H = torch.zeros(N, C, 4, dtype=torch.float32, device=device)
            rows, cols = torch.nonzero(m, as_tuple=True)
            H[rows, cols, argmax_allele[rows, cols]] = 1.0
            H_flat = H.reshape(N, -1)
            matches += H_flat @ H_flat.T
        elif ani_method == "popani":
            c_f = c.float() if is_mps else c.to(torch.int64)
            surr = torch.einsum('ilk,jlk->ijl', c_f, c_f)
            mask_pair = m.unsqueeze(1) & m.unsqueeze(0)
            surr[~mask_pair] = 0
            matches += (surr > 0).sum(dim=2).float()
        elif ani_method.startswith("cosani_"):
            threshold = float(ani_method.split("_", 1)[1])
            c_f = c.float()
            norms = c_f.pow(2).sum(dim=2, keepdim=True).sqrt()
            norms = torch.where(norms == 0, torch.ones_like(norms), norms)
            normed = c_f / norms
            cos_sim = torch.einsum('ilk,jlk->ijl', normed, normed)
            mask_pair = m.unsqueeze(1) & m.unsqueeze(0)
            cos_sim[~mask_pair] = 0.0
            matches += (cos_sim >= threshold).sum(dim=2).float()
        else:
            raise ValueError(f"Unknown ani_method: {ani_method}")

    return matches.cpu().to(torch.float64).numpy(), shared.cpu().to(torch.float64).numpy()


def compare_genes(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = options.COMPARE_MIN_COV,
    min_gene_compare_len: int = options.MIN_GENE_COMPARE_LEN,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
) -> pl.LazyFrame:
    """Compare two profiles at the gene level using Polars."""
    return compare_genes_polars(
        mpile_contig_1=mpile_contig_1,
        mpile_contig_2=mpile_contig_2,
        min_cov=min_cov,
        min_gene_compare_len=min_gene_compare_len,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
        ani_method=ani_method,
    )
