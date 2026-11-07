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
from zipstrain.profile import PROFILE_SORTED_METADATA_KEY, PROFILE_SORTED_METADATA_VALUE
from zipstrain.utils import _read_custom_parquet_metadata

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
    sorted_inputs = all(
        isinstance(source, (str, Path))
        and _read_custom_parquet_metadata(source).get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE
        for source in (mpile1, mpile2)
    )
    joined = _join_position_frames(frames, sorted_inputs=sorted_inputs)
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


def add_contiguity_info(mpile_contig: pl.LazyFrame, *, coordinates_sorted: bool = False) -> pl.LazyFrame:
    """ Adds group id information to the lazy frame. If on the same scaffold and not popANI, then they are in the same group.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
        coordinates_sorted: Rows are already grouped by scaffold and in position order.

    Returns:
        pl.LazyFrame: Updated LazyFrame with group id information added.
    """

    has_genome = "genome" in mpile_contig.collect_schema().names()
    sort_cols = (["genome"] if has_genome else []) + ["scaffold", "pos"]

    def changed(column: str) -> pl.Expr:
        """True on the first row of each new run of `column` (once sorted)."""
        return pl.col(column) != pl.col(column).shift(1).fill_null(pl.col(column))

    # A contiguous block is a run of consecutive shared positions. It breaks at
    # the first row of a new scaffold (or new genome, when that column exists),
    # at a gap in genomic coordinates, or at a position that is itself an SNV
    # (surr == 0). Marking those break points and taking a cumulative sum
    # numbers the blocks, so `group_id` is constant within each block and
    # increments at every break.
    position_gap = (pl.col("pos") - pl.col("pos").shift(1) != 1).fill_null(True)
    break_expr = changed("scaffold") | position_gap | (pl.col("surr") == 0)
    if has_genome:
        break_expr = (
            changed("genome")
            | changed("scaffold")
            | position_gap
            | (pl.col("surr") == 0)
        )

    ordered = mpile_contig if coordinates_sorted else mpile_contig.sort(sort_cols)
    return ordered.with_columns(
        break_expr.cast(pl.Int64).cum_sum().alias("group_id")
    )

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

def calculate_pop_ani(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """
    Calculates the population ANI (Average Nucleotide Identity) for the given mpileup LazyFrame.
    NOTE: Remember that this function should be applied to the merged mpileup using `_shared_loci_polars`
    or ``align_profiles`` followed by an ANI match expression.

    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with population ANI information added.
    """
    return mpile_contig.group_by("genome").agg(
            total_positions=pl.len(),
            share_allele_pos=(pl.col("surr") > 0 ).sum()
        ).with_columns(
            genome_ani=pl.col("share_allele_pos")/pl.col("total_positions")*100,
        )

def get_longest_consecutive_blocks(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """
    Calculates the longest consecutive blocks for each genome in the mpileup LazyFrame for any genome.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with longest consecutive blocks information added.
    """
    block_lengths = (
        mpile_contig.group_by(["genome", "scaffold", "group_id"])
        .agg(pl.len().alias("length"))
    ) 
    return block_lengths.group_by("genome").agg(pl.col("length").max().alias("max_identical_stretch"))

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
    )
    genome_comp_parts: list[pl.LazyFrame] = []
    if "genome_ani" in calculations:
        genome_comp_parts.append(calculate_pop_ani(shared))
    if "max_identical_stretch" in calculations:
        genome_comp_parts.append(get_longest_consecutive_blocks(add_contiguity_info(shared)))

    genome_comp = genome_comp_parts[0]
    for part in genome_comp_parts[1:]:
        genome_comp = genome_comp.join(part, on="genome", how="left")

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
    """
    if min_cov < 1:
        raise ValueError("min_cov must be positive")
    require_same_profile_generation((profile_1, profile_2))
    return _align_inputs(_alignment_input(profile_1), _alignment_input(profile_2), min_cov)


def _index_genes(genes: pl.DataFrame) -> dict:
    """Group a catalog's gene rows and coordinates by scaffold for interval counting."""
    starts, ends = genes["start"].to_numpy(), genes["end"].to_numpy()
    index = {}
    for (chrom,), rows in genes.with_row_index("row").group_by("chrom", maintain_order=True):
        row = rows["row"].to_numpy()
        index[chrom] = (row, starts[row], ends[row])
    return index


def _scan_shared_positions(shared: pl.LazyFrame, genes: Optional[pl.DataFrame], *, stretches: bool,
                           chunk_size: int = 1 << 20, gene_index: Optional[dict] = None,
                           ) -> tuple[pl.DataFrame, Optional[pl.DataFrame]]:
    """Calculate genome and gene counts from one streaming pass over aligned positions.

    ``shared`` holds ``chrom``, ``genome``, ``pos`` and boolean ``match`` columns, with rows
    grouped by scaffold and in position order. Return one frame per genome with
    ``total_positions``, ``share_allele_pos`` and ``max_identical_stretch``, and, when
    ``genes`` is given, one per gene with ``gene_total_positions`` and
    ``gene_share_allele_pos``. Overlapping genes each count all their own positions.
    A stretch breaks at a mismatch, a coordinate gap, or a new scaffold.
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
        boundaries = np.flatnonzero(scaffolds.ne_missing(scaffolds.shift(1)).to_numpy())
        for first, last in zip(boundaries, np.r_[boundaries[1:], batch.height], strict=True):
            chrom, genome = scaffolds[int(first)], batch["genome"][int(first)]
            pos, match = positions[first:last], matches[first:last]
            counts = genome_counts.setdefault(genome, [0, 0, 0])
            counts[0] += pos.size
            counts[1] += int(match.sum())
            if stretches:
                edges = np.r_[np.flatnonzero(~match | np.r_[False, np.diff(pos) != 1]), pos.size]
                if edges[0] != 0:
                    edges = np.r_[0, edges]
                # A run opened by a mismatch does not count that position.
                runs = np.diff(edges) - (~match[edges[:-1]])
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
    """
    from zipstrain import codons, gene_info as gi, matrix_pairs as mp
    if min_gene_compare_len < 1:
        raise ValueError("min_gene_compare_len must be positive")
    if min_cov < 1:
        raise ValueError("min_cov must be positive")
    genes, gene_file = None, None
    if preprofile_dir is not None:
        gene_file = Path(preprofile_dir)/"gene_info.parquet"
        gi.validate_profile_assets([profile_1, profile_2], gene_file)
        genes = gi.read_gene_info(gene_file)
    names = (sample_1 or (Path(profile_1).stem.removesuffix("_profile") if isinstance(profile_1, (str, Path)) else "sample_1"),
             sample_2 or (Path(profile_2).stem.removesuffix("_profile") if isinstance(profile_2, (str, Path)) else "sample_2"))
    codon_tables = {name: counts for name, counts in zip(names, (codon_counts_1, codon_counts_2), strict=True)}
    calculations = options.resolve_metrics(
        calculate, genes_available=genes is not None and genes.height > 0,
        codon_tables_present=[name for name, counts in codon_tables.items() if counts is not None],
        codon_tables_missing=[name for name, counts in codon_tables.items() if counts is None],
    )
    rates = None
    if "gene_selection" in calculations:
        from zipstrain.utils import sha256_file
        for counts, positions in ((codon_counts_1, profile_1), (codon_counts_2, profile_2)):
            if pl.read_parquet_metadata(counts).get("zipstrain_source_profile_sha256") != sha256_file(positions):
                raise ValueError("Codon evidence does not match its sample position profile")
    ani_method, _threshold = mp.parse_matrix_ani_method(ani_method)
    require_same_profile_generation((profile_1, profile_2))
    first, second = _alignment_input(profile_1), _alignment_input(profile_2)
    # Reject incompatible samples before spending time on selection.
    aligned = _align_inputs(first, second, min_cov)
    if "gene_selection" in calculations:
        rates = codons.pairwise_selection(codon_counts_1, codon_counts_2, gene_file, min_cov=min_cov,
                                          allele_integration=allele_integration)
    universe = [_as_lazy_profile(source).select(pl.col("genome").cast(pl.String)).unique()
                for source in (profile_1, profile_2)]
    return _compare_aligned(
        aligned, names, calculations=calculations, ani_method=ani_method,
        min_gene_compare_len=min_gene_compare_len, genes=genes, rates=rates,
        genome_universe=[*universe, *_catalog_genomes(preprofile_dir, genes)],
    )


def _catalog_genomes(preprofile_dir, genes: Optional[pl.DataFrame]) -> list[pl.LazyFrame]:
    """Return the genomes a comparison must report even when no position is shared."""
    lengths = Path(preprofile_dir)/"genome_lengths.parquet" if preprofile_dir is not None else None
    if lengths is not None and lengths.exists():
        return [pl.scan_parquet(lengths).select(pl.col("genome").cast(pl.String)).unique()]
    return [genes.lazy().select(pl.col("genome").cast(pl.String)).unique()] if genes is not None else []


def _compare_aligned(aligned: pl.LazyFrame, names: tuple[str, str], *, calculations: tuple[str, ...],
                     ani_method: str, min_gene_compare_len: int, genes: Optional[pl.DataFrame],
                     genome_universe: list, rates: Optional[pl.LazyFrame] = None,
                     gene_index: Optional[dict] = None) -> pl.LazyFrame:
    """Build one pair's comparison rows from already validated inputs and resolved metrics."""
    from zipstrain import gene_info as gi
    shared = aligned.select(
        "chrom", "genome", "pos", (getattr(PolarsANIExpressions(), ani_method)() > 0).alias("match"))
    keys = ["sample_1", "sample_2", "genome"]
    sample_columns = [pl.lit(name).alias(key) for key, name in zip(keys[:2], names, strict=True)]
    gene_rows = "gene_ani" in calculations or rates is not None
    wants_ani, wants_stretch = "genome_ani" in calculations, "max_identical_stretch" in calculations
    genome_counts, gene_counts = _scan_shared_positions(
        shared, genes if gene_rows else None, stretches=wants_stretch, gene_index=gene_index)
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
    if gene_rows:
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
