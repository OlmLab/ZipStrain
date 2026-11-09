"""Shared comparison and profiling defaults, and the one ``--calculate`` parser.

Every signature, CLI option and workflow reads its defaults and metric names
from here, so entry points cannot disagree. This module imports nothing from
ZipStrain and is safe to import anywhere.
"""
from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from typing import Optional, Union

# Comparison defaults.
COMPARE_MIN_COV = 5
# Comparison tables keep every gene with a shared position. Analyses choose their
# own stricter threshold (see ``ANALYSIS_MIN_GENE_POSITIONS``), and a table filtered
# at this value can be filtered further later, while the reverse is not possible.
MIN_GENE_COMPARE_LEN = 1
ANI_METHOD = "popani"
ALLELE_INTEGRATION = "consensus"
CALCULATE = "all"
ANALYSIS_MIN_GENE_POSITIONS = 100

# Profiling read filters. Low-identity reads are typically mis-mapped, and half-mapped
# or cross-scaffold pairs are dropped while single-end reads stay (inStrain's defaults).
PROFILE_MIN_READ_ANI = 0.95
PROFILE_READ_INCLUSION = "paired"
PROFILE_READ_INCLUSION_CHOICES = ("proper-pairs", "paired", "all-mapped")
PROFILE_MIN_MAPQ = 0
PROFILE_MIN_BASEQ = 13
PROFILE_MIN_FREQ = 0.01
# Sequencing-error null model that profiling uses to drop error-level allele counts.
NULL_MODEL_ERROR_RATE = 0.001
NULL_MODEL_MAX_COVERAGE = 50_000
NULL_MODEL_P_THRESHOLD = 1e-6
# Profiles carry this Parquet metadata when rows are sorted by scaffold and position.
PROFILE_SORTED_METADATA_KEY = "zipstrain_sorted_by"
PROFILE_SORTED_METADATA_VALUE = "chrom,pos"

# Choices the CLI shows while it loads. They live here so building the CLI does not
# import the modules that implement them; those modules re-export these names.
MATRIX_PAIR_BACKENDS = ("numpy", "torch", "torch-cpu", "torch-cuda", "torch-mps")
MATRIX_COMPARE_KERNEL_ANCHOR = "anchor"
MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE = 1
MATRIX_STORAGE_MODES = ("bitmask", "counts")
MATRIX_COUNT_DTYPE_CHOICES = ("auto", "uint16", "uint32")
MATRIX_IO_EXECUTOR_KINDS = ("thread", "process")
MATRIX_HDF5_EXPORT_TARGET_BATCH_MB = 128.0
ABUNDANCE_TOOLS = ("sylph",)  # must match build_db.ADAPTERS
SYLPH_DB_URL = "http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"

GENOME_METRICS = ("genome_ani", "max_identical_stretch")
GENE_METRICS = ("gene_ani", "gene_selection")
METRICS = (*GENOME_METRICS, *GENE_METRICS)
CALCULATE_HELP = (
    "Metrics to compute, combined with '+': genome_ani; max_identical_stretch; gene_ani; "
    "gene_selection: dN, dS, dN/dS, pN, pS, pN/pS per gene (requires --collect-codons); "
    "or all, meaning every metric the inputs support."
)
MIN_GENE_COMPARE_LEN_HELP = "Minimum shared positions for a gene to get a comparison row."

# Removed spellings, kept only to point people at the current name.
_RENAMED_METRICS = {
    "ani": "genome_ani",
    "popani": "genome_ani",
    "ibs": "max_identical_stretch",
    "max_block": "max_identical_stretch",
    "max_consecutive_length": "max_identical_stretch",
    "gene": "gene_ani",
    "genes": "gene_ani",
    "identical_genes": "gene_ani",
    "gene_identity": "gene_ani",
    "id_genes": "gene_ani",
    "dnds": "gene_selection",
}


def _tokens(calculate: Optional[Union[str, Iterable[str]]]) -> list[str]:
    if calculate is None:
        return [CALCULATE]
    parts = [calculate] if isinstance(calculate, str) else calculate
    tokens = [
        token.strip().lower()
        for part in parts
        for token in str(part).replace(",", "+").split("+")
        if token.strip()
    ]
    if not tokens:
        raise ValueError(f"No metric requested. Choose from: all, {', '.join(METRICS)}")
    for token in tokens:
        if token != "all" and token not in METRICS:
            suggestion = _RENAMED_METRICS.get(token)
            hint = f" Did you mean '{suggestion}'?" if suggestion else f" Choose from: all, {', '.join(METRICS)}"
            raise ValueError(f"Unknown metric '{token}'.{hint}")
    return tokens


def explicit_metrics(calculate: Optional[Union[str, Iterable[str]]]) -> tuple[str, ...]:
    """Return metrics named directly, in canonical order, leaving ``all`` out."""
    tokens = _tokens(calculate)
    return tuple(metric for metric in METRICS if metric in tokens)


def requests_all(calculate: Optional[Union[str, Iterable[str]]]) -> bool:
    """Return whether the selection contains ``all``."""
    return "all" in _tokens(calculate)


def parse_calculate(calculate: Optional[Union[str, Iterable[str]]] = None) -> tuple[str, ...]:
    """Validate a selection and return its metrics, with ``all`` naming every metric."""
    return METRICS if requests_all(calculate) else explicit_metrics(calculate)


def resolve_metrics(
    calculate: Optional[Union[str, Iterable[str]]],
    *,
    genes_available: bool,
    codon_tables_present: Sequence[str] = (),
    codon_tables_missing: Sequence[str] = (),
) -> tuple[str, ...]:
    """Return the metrics a run will calculate from what its inputs support.

    A metric named explicitly must be satisfiable or this raises ``ValueError``.
    ``all`` keeps every supported metric and warns about each one it skips. Codon
    tables present for only some samples are an error either way, because the
    result would silently depend on which samples were profiled with codons.
    """
    explicit = explicit_metrics(calculate)
    selected = list(parse_calculate(calculate))
    if not genes_available:
        for metric in GENE_METRICS:
            if metric in explicit:
                raise ValueError(
                    f"{metric} requires a gene catalog: pass --preprofile-dir with a "
                    "gene_info.parquet that contains genes."
                )
        if "gene_ani" in selected:
            warnings.warn(
                "Skipping gene_ani and gene_selection: no gene catalog with genes was provided "
                "(--preprofile-dir). Calculating genome metrics only.",
                stacklevel=2,
            )
        return tuple(metric for metric in selected if metric in GENOME_METRICS)
    if "gene_selection" not in selected:
        return tuple(selected)
    missing = sorted(codon_tables_missing)
    if missing and ("gene_selection" in explicit or codon_tables_present):
        raise ValueError(
            "gene_selection requires a codon table for every sample (profile with "
            f"--collect-codons). Missing for: {', '.join(missing)}"
        )
    if missing:
        warnings.warn(
            "Skipping gene_selection: no sample has a codon table (profile with --collect-codons).",
            stacklevel=2,
        )
        selected.remove("gene_selection")
    return tuple(selected)
