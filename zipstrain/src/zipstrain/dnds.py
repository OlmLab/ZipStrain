"""Codon summaries and pairwise dN/dS calculations.

The database annotation remains one row per gene. Sample codons are materialized
only when requested and only for covered coding positions, keeping both the
default profiling path and large reference databases lightweight.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import permutations
from pathlib import Path
from typing import Optional, Union

import duckdb
import polars as pl


BASES = "ACGT"
CODON_PROFILE_SCHEMA_VERSION = "2"
CODON_SCHEMA_METADATA_KEY = "zipstrain_codon_profile_schema"
GENE_INFO_HASH_METADATA_KEY = "zipstrain_gene_info_hash"
ALLELE_INTEGRATION_METADATA_KEY = "zipstrain_dnds_allele_integration"
DEFAULT_MEMORY_LIMIT = "1GB"
ALLELE_INTEGRATION_MODES = ("consensus", "weighted")
# Polars common-subplan elimination can corrupt the repeated lazy branches used
# for SNS/SNV summaries. Only the compact per-gene result is materialized with
# this optimization disabled; the multi-million-row codon input stays lazy.
DNDS_QUERY_OPTIMIZATIONS = pl.QueryOptFlags(comm_subplan_elim=False)

CODON_PROFILE_COLUMNS = (
    "gene_id",
    "codon_id",
    "reference_codon_code",
    "min_coverage",
    *(f"{base}_{offset}" for offset in range(3) for base in BASES),
)

DNDS_RESULT_COLUMNS = (
    "callable_codons",
    "sns_count",
    "snv_count",
    "sns_synonymous_changes",
    "sns_nonsynonymous_changes",
    "snv_synonymous_changes",
    "snv_nonsynonymous_changes",
    "synonymous_sites",
    "nonsynonymous_sites",
    "pS",
    "pN",
    "pN_pS",
    "dS",
    "dN",
    "dN_dS",
    "sns_stop_changes",
    "snv_stop_changes",
    "allele_tie_sites",
)

REFERENCE_DNDS_RESULT_COLUMNS = tuple(f"ref_{name}" for name in DNDS_RESULT_COLUMNS)


# Tables 1 and 11 use the same internal codon translations. Alternative start
# codons are relevant only at initiation and do not change substitutions within
# a coding sequence.
STANDARD_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def codon_profile_path(profile_path: Union[str, Path]) -> Path:
    """Return the deterministic codon-sidecar path for a profile parquet."""
    path = Path(profile_path)
    suffix = "_profile.parquet"
    if not path.name.endswith(suffix):
        raise ValueError(
            f"Profile filename must end with {suffix!r} to locate its codon sidecar: {path}"
        )
    return path.with_name(path.name[: -len(suffix)] + "_codon_profile.parquet")


def _decode_codon(code: int) -> str:
    return "".join(
        (
            BASES[(code // 16) % 4],
            BASES[(code // 4) % 4],
            BASES[code % 4],
        )
    )


def _encode_codon(codon: str) -> int:
    values = {base: index for index, base in enumerate(BASES)}
    return values[codon[0]] * 16 + values[codon[1]] * 4 + values[codon[2]]


def _synonymous_sites(codon: str) -> float:
    amino_acid = STANDARD_GENETIC_CODE[codon]
    if amino_acid == "*":
        return 0.0
    synonymous = 0
    for position in range(3):
        for alt in BASES:
            if alt == codon[position]:
                continue
            mutated = codon[:position] + alt + codon[position + 1 :]
            if STANDARD_GENETIC_CODE[mutated] == amino_acid:
                synonymous += 1
    return synonymous / 3.0


def _path_changes(codon_1: str, codon_2: str) -> tuple[float, float]:
    differing = [index for index in range(3) if codon_1[index] != codon_2[index]]
    if not differing:
        return 0.0, 0.0

    valid_paths: list[tuple[int, int]] = []
    fallback_paths: list[tuple[int, int]] = []
    for order in permutations(differing):
        current = list(codon_1)
        synonymous = 0
        nonsynonymous = 0
        has_intermediate_stop = False
        for step_index, position in enumerate(order):
            before = "".join(current)
            current[position] = codon_2[position]
            after = "".join(current)
            if STANDARD_GENETIC_CODE[before] == STANDARD_GENETIC_CODE[after]:
                synonymous += 1
            else:
                nonsynonymous += 1
            if step_index < len(order) - 1 and STANDARD_GENETIC_CODE[after] == "*":
                has_intermediate_stop = True
        fallback_paths.append((synonymous, nonsynonymous))
        if not has_intermediate_stop:
            valid_paths.append((synonymous, nonsynonymous))

    paths = valid_paths or fallback_paths
    return (
        sum(path[0] for path in paths) / len(paths),
        sum(path[1] for path in paths) / len(paths),
    )


@lru_cache(maxsize=1)
def codon_pair_lookup() -> pl.DataFrame:
    """Return the 4096-row lookup used by all dN/dS calculations."""
    rows: list[tuple[int, float, float, float, float, bool, int]] = []
    for code_1 in range(64):
        codon_1 = _decode_codon(code_1)
        aa_1 = STANDARD_GENETIC_CODE[codon_1]
        for code_2 in range(64):
            codon_2 = _decode_codon(code_2)
            aa_2 = STANDARD_GENETIC_CODE[codon_2]
            stop_change = int((aa_1 == "*") != (aa_2 == "*"))
            valid = aa_1 != "*" and aa_2 != "*"
            if valid:
                synonymous, nonsynonymous = _path_changes(codon_1, codon_2)
                synonymous_sites = (
                    _synonymous_sites(codon_1) + _synonymous_sites(codon_2)
                ) / 2.0
                nonsynonymous_sites = 3.0 - synonymous_sites
            else:
                synonymous = 0.0
                nonsynonymous = 0.0
                synonymous_sites = 0.0
                nonsynonymous_sites = 0.0
            rows.append(
                (
                    code_1 * 64 + code_2,
                    synonymous,
                    nonsynonymous,
                    synonymous_sites,
                    nonsynonymous_sites,
                    valid,
                    stop_change,
                )
            )
    return pl.DataFrame(
        rows,
        schema={
            "pair_code": pl.UInt16,
            "synonymous_changes": pl.Float64,
            "nonsynonymous_changes": pl.Float64,
            "synonymous_sites": pl.Float64,
            "nonsynonymous_sites": pl.Float64,
            "valid_dnds": pl.Boolean,
            "stop_change": pl.UInt8,
        },
        orient="row",
    )


def codon_profile_exists(profile_path: Union[str, Path]) -> bool:
    """True when a profile has a usable codon sidecar beside it.

    Used to decide whether ``--calculate all`` can include dN/dS, so a profile
    whose name does not follow the sidecar convention is reported as simply not
    having one rather than raising.
    """
    if profile_path is None:
        return False
    try:
        return codon_profile_path(profile_path).exists()
    except (ValueError, TypeError):
        return False


def _quote(value: Union[str, Path]) -> str:
    return str(value).replace("'", "''")


def validate_allele_integration(value: str) -> str:
    """Return a canonical allele-integration mode or raise a useful error."""
    mode = str(value).strip().lower()
    if mode not in ALLELE_INTEGRATION_MODES:
        supported = ", ".join(ALLELE_INTEGRATION_MODES)
        raise ValueError(f"allele_integration must be one of: {supported}")
    return mode


def write_codon_profile(
    *,
    profile_path: Union[str, Path],
    gene_info_path: Union[str, Path],
    output_path: Union[str, Path],
    temp_directory: Union[str, Path],
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
    metadata: Optional[dict[str, str]] = None,
) -> None:
    """Build a covered-codon, coordinate-independent allele summary out of core."""
    from zipstrain import utils

    profile_path = Path(profile_path)
    gene_info_path = Path(gene_info_path)
    output_path = Path(output_path)
    temp_directory = Path(temp_directory)
    temp_directory.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    profile_columns = set(pl.scan_parquet(profile_path).collect_schema().names())
    if utils.REF_BASE_BITMASK_COLUMN not in profile_columns:
        raise ValueError(
            "--prepare-dnds requires a reference FASTA so the profile contains "
            "ref_base_bitmask. Re-run profiling with --reference-fasta."
        )

    gene_columns = set(pl.scan_parquet(gene_info_path).collect_schema().names())
    required_gene_columns = {
        "gene_id", "scaffold", "start", "end", "strand", "phase",
        "genetic_code", "first_codon_id", "n_codons",
    }
    missing = required_gene_columns - gene_columns
    if missing:
        raise ValueError(
            "Gene information table is missing required dN/dS columns: "
            + ", ".join(sorted(missing))
        )
    unsupported_codes = (
        pl.scan_parquet(gene_info_path)
        .filter(~pl.col("genetic_code").is_in([1, 11]))
        .select("genetic_code")
        .unique()
        .collect(engine="streaming")
        .get_column("genetic_code")
        .to_list()
    )
    if unsupported_codes:
        raise ValueError(
            "dN/dS currently supports genetic codes 1 and 11; unsupported "
            f"codes found: {', '.join(map(str, sorted(unsupported_codes)))}"
        )

    raw_output = temp_directory / f"{output_path.stem}.duckdb.parquet"
    raw_output.unlink(missing_ok=True)
    con = duckdb.connect()
    try:
        con.execute(f"SET memory_limit='{_quote(memory_limit)}'")
        con.execute(f"SET temp_directory='{_quote(temp_directory)}'")
        con.execute(
            f"""
            COPY (
              WITH mapped AS (
                SELECT
                  g.gene_id::UINTEGER AS gene_id,
                  g.first_codon_id::UBIGINT AS first_codon_id,
                  g.n_codons::UBIGINT AS n_codons,
                  g.strand::TINYINT AS strand,
                  g.phase::INTEGER AS phase,
                  p.pos::BIGINT AS pos,
                  g.start::BIGINT AS gene_start,
                  g.end::BIGINT AS gene_end,
                  (p.A::BIGINT + p.C::BIGINT + p.G::BIGINT + p.T::BIGINT) AS coverage,
                  p.A::UINTEGER AS genomic_A,
                  p.C::UINTEGER AS genomic_C,
                  p.G::UINTEGER AS genomic_G,
                  p.T::UINTEGER AS genomic_T,
                  CASE p.ref_base_bitmask
                    WHEN 1 THEN 0 WHEN 2 THEN 1 WHEN 4 THEN 2 WHEN 8 THEN 3
                    ELSE NULL
                  END::INTEGER AS genomic_reference_base
                FROM read_parquet('{_quote(profile_path)}') p
                JOIN read_parquet('{_quote(gene_info_path)}') g
                  ON p.chrom = g.scaffold
                 AND p.pos BETWEEN g.start AND g.end
              ), oriented AS (
                SELECT *,
                  CASE WHEN strand = 1
                    THEN pos - gene_start - phase
                    ELSE gene_end - pos - phase
                  END AS coding_index,
                  CASE WHEN strand = 1
                    THEN genomic_reference_base
                    ELSE 3 - genomic_reference_base
                  END AS coding_reference_base,
                  CASE WHEN strand = 1 THEN genomic_A ELSE genomic_T END AS coding_A,
                  CASE WHEN strand = 1 THEN genomic_C ELSE genomic_G END AS coding_C,
                  CASE WHEN strand = 1 THEN genomic_G ELSE genomic_C END AS coding_G,
                  CASE WHEN strand = 1 THEN genomic_T ELSE genomic_A END AS coding_T
                FROM mapped
                WHERE coverage > 0
                  AND genomic_reference_base IS NOT NULL
              ), coded AS (
                SELECT
                  gene_id,
                  first_codon_id + FLOOR(coding_index / 3)::UBIGINT AS codon_id,
                  (coding_index % 3)::INTEGER AS codon_offset,
                  coding_reference_base,
                  coverage,
                  coding_A,
                  coding_C,
                  coding_G,
                  coding_T
                FROM oriented
                WHERE coding_index >= 0
                  AND coding_index < n_codons * 3
              )
              SELECT
                gene_id,
                codon_id,
                SUM(coding_reference_base * CASE codon_offset WHEN 0 THEN 16 WHEN 1 THEN 4 ELSE 1 END)::UTINYINT AS reference_codon_code,
                MIN(coverage)::UINTEGER AS min_coverage,
                MAX(CASE WHEN codon_offset = 0 THEN coding_A END)::UINTEGER AS A_0,
                MAX(CASE WHEN codon_offset = 0 THEN coding_C END)::UINTEGER AS C_0,
                MAX(CASE WHEN codon_offset = 0 THEN coding_G END)::UINTEGER AS G_0,
                MAX(CASE WHEN codon_offset = 0 THEN coding_T END)::UINTEGER AS T_0,
                MAX(CASE WHEN codon_offset = 1 THEN coding_A END)::UINTEGER AS A_1,
                MAX(CASE WHEN codon_offset = 1 THEN coding_C END)::UINTEGER AS C_1,
                MAX(CASE WHEN codon_offset = 1 THEN coding_G END)::UINTEGER AS G_1,
                MAX(CASE WHEN codon_offset = 1 THEN coding_T END)::UINTEGER AS T_1,
                MAX(CASE WHEN codon_offset = 2 THEN coding_A END)::UINTEGER AS A_2,
                MAX(CASE WHEN codon_offset = 2 THEN coding_C END)::UINTEGER AS C_2,
                MAX(CASE WHEN codon_offset = 2 THEN coding_G END)::UINTEGER AS G_2,
                MAX(CASE WHEN codon_offset = 2 THEN coding_T END)::UINTEGER AS T_2
              FROM coded
              GROUP BY gene_id, codon_id
              HAVING COUNT(*) = 3 AND COUNT(DISTINCT codon_offset) = 3
              ORDER BY codon_id
            ) TO '{_quote(raw_output)}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        con.close()

    output_metadata = dict(metadata or {})
    output_metadata[CODON_SCHEMA_METADATA_KEY] = CODON_PROFILE_SCHEMA_VERSION
    output_metadata[GENE_INFO_HASH_METADATA_KEY] = utils.sha256_file(gene_info_path)
    pl.scan_parquet(raw_output).select(CODON_PROFILE_COLUMNS).sink_parquet(
        output_path,
        compression="zstd",
        engine="streaming",
        metadata=output_metadata,
    )
    raw_output.unlink(missing_ok=True)


@lru_cache(maxsize=1)
def codon_site_lookup() -> pl.DataFrame:
    """Reference-codon synonymous and nonsynonymous opportunities."""
    rows = []
    for code in range(64):
        codon = _decode_codon(code)
        valid = STANDARD_GENETIC_CODE[codon] != "*"
        synonymous = _synonymous_sites(codon) if valid else 0.0
        rows.append((code, synonymous, 3.0 - synonymous if valid else 0.0, valid))
    return pl.DataFrame(
        rows,
        schema={
            "reference_codon_code": pl.UInt8,
            "synonymous_sites": pl.Float64,
            "nonsynonymous_sites": pl.Float64,
            "valid_reference": pl.Boolean,
        },
        orient="row",
    )


def _reference_base_expr() -> pl.Expr:
    multiplier = pl.when(pl.col("codon_offset") == 0).then(16).when(
        pl.col("codon_offset") == 1
    ).then(4).otherwise(1)
    return (
        (pl.col("reference_codon_code").cast(pl.UInt16) // multiplier) % 4
    ).cast(pl.UInt8)


def _replace_codon_base(code: pl.Expr, offset: pl.Expr, base: pl.Expr) -> pl.Expr:
    multiplier = pl.when(offset == 0).then(16).when(offset == 1).then(4).otherwise(1)
    old_base = (code.cast(pl.Int32) // multiplier) % 4
    return (code.cast(pl.Int32) + (base.cast(pl.Int32) - old_base) * multiplier).cast(
        pl.UInt8
    )


def _codon_positions(codons: pl.LazyFrame) -> pl.LazyFrame:
    """Expand one codon row into three rows without expanding allele choices."""
    frames = []
    for offset in range(3):
        frames.append(
            codons.select(
                "gene_id",
                "codon_id",
                "reference_codon_code",
                pl.lit(offset, dtype=pl.UInt8).alias("codon_offset"),
                *(pl.col(f"{base}_{offset}").alias(base) for base in BASES),
            )
        )
    return (
        pl.concat(frames)
        .with_columns(
            coverage=pl.sum_horizontal(*(pl.col(base) for base in BASES)).cast(pl.UInt64),
            allele_count=pl.sum_horizontal(
                *((pl.col(base) > 0).cast(pl.UInt8) for base in BASES)
            ),
            reference_base=_reference_base_expr(),
        )
    )


def _singleton_base_expr(columns: tuple[str, ...] = tuple(BASES)) -> pl.Expr:
    expr = pl.when(pl.col(columns[0]) > 0).then(pl.lit(0, dtype=pl.UInt8))
    for index, column in enumerate(columns[1:], start=1):
        expr = expr.when(pl.col(column) > 0).then(pl.lit(index, dtype=pl.UInt8))
    return expr.otherwise(None)


def _event_frame(
    positions: pl.LazyFrame,
    *,
    source_base: pl.Expr,
    target_base: pl.Expr,
    weight: pl.Expr,
) -> pl.LazyFrame:
    source_code = _replace_codon_base(
        pl.col("reference_codon_code"), pl.col("codon_offset"), source_base
    )
    target_code = _replace_codon_base(
        pl.col("reference_codon_code"), pl.col("codon_offset"), target_base
    )
    return positions.select(
        "gene_id",
        "codon_id",
        "codon_offset",
        (source_code.cast(pl.UInt16) * 64 + target_code.cast(pl.UInt16)).alias(
            "pair_code"
        ),
        weight.cast(pl.Float64).alias("weight"),
    )


def _reference_event_frames(
    codons: pl.LazyFrame,
    allele_integration: str,
) -> tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]:
    positions = _codon_positions(codons)
    sample_base = _singleton_base_expr()
    sns_sites = positions.filter(
        (pl.col("allele_count") == 1) & (sample_base != pl.col("reference_base"))
    )
    sns_events = _event_frame(
        sns_sites,
        source_base=pl.col("reference_base"),
        target_base=sample_base,
        weight=pl.lit(1.0),
    )

    snv_sites = positions.filter(pl.col("allele_count") > 1)
    alternative_counts = [
        pl.when(pl.col("reference_base") != index)
        .then(pl.col(base).cast(pl.Float64))
        .otherwise(0.0)
        for index, base in enumerate(BASES)
    ]
    snv_sites = snv_sites.with_columns(
        max_support=pl.max_horizontal(*alternative_counts)
    ).with_columns(
        tied_candidates=pl.sum_horizontal(
            *(
                (
                    (pl.col("reference_base") != index)
                    & (pl.col(base).cast(pl.Float64) == pl.col("max_support"))
                    & (pl.col(base) > 0)
                ).cast(pl.UInt8)
                for index, base in enumerate(BASES)
            )
        )
    )
    tie_sites = snv_sites.filter(pl.col("tied_candidates") > 1).select(
        "gene_id", "codon_id", "codon_offset"
    )

    candidates = []
    for index, base in enumerate(BASES):
        candidate = snv_sites.filter(
            (pl.col(base) > 0) & (pl.col("reference_base") != index)
        )
        if allele_integration == "consensus":
            candidate = candidate.filter(
                pl.col(base).cast(pl.Float64) == pl.col("max_support")
            )
            weight = 1.0 / pl.col("tied_candidates")
        else:
            weight = pl.col(base).cast(pl.Float64) / pl.col("coverage")
        candidates.append(
            _event_frame(
                candidate,
                source_base=pl.col("reference_base"),
                target_base=pl.lit(index, dtype=pl.UInt8),
                weight=weight,
            )
        )
    return sns_events, pl.concat(candidates), sns_sites, snv_sites, tie_sites


def _pairwise_event_frames(
    codons: pl.LazyFrame,
    allele_integration: str,
) -> tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]:
    left = _codon_positions(codons)
    right = _codon_positions(
        codons.select(
            "gene_id",
            "codon_id",
            "reference_codon_code",
            *(
                pl.col(f"{column}_right").alias(column)
                for column in CODON_PROFILE_COLUMNS[4:]
            ),
        )
    )
    # The caller supplies an already joined codon frame. Its right-hand counts
    # carry a `_right` suffix; normalize those after position expansion.
    right = right.select(
        "codon_id",
        "codon_offset",
        *(pl.col(base).alias(f"{base}_right") for base in BASES),
        pl.col("coverage").alias("coverage_right"),
        pl.col("allele_count").alias("allele_count_right"),
    )
    positions = left.join(
        right, on=["codon_id", "codon_offset"], how="inner", maintain_order="left"
    )
    left_base = _singleton_base_expr()
    right_base = _singleton_base_expr(tuple(f"{base}_right" for base in BASES))
    sns_sites = positions.filter(
        (pl.col("allele_count") == 1)
        & (pl.col("allele_count_right") == 1)
        & (left_base != right_base)
    )
    sns_events = _event_frame(
        sns_sites,
        source_base=left_base,
        target_base=right_base,
        weight=pl.lit(1.0),
    )

    snv_sites = positions.filter(
        (pl.col("allele_count") > 1) | (pl.col("allele_count_right") > 1)
    )
    supports = [
        pl.col(left_base_name).cast(pl.Float64) * pl.col(f"{right_base_name}_right")
        for left_index, left_base_name in enumerate(BASES)
        for right_index, right_base_name in enumerate(BASES)
        if left_index != right_index
    ]
    snv_sites = snv_sites.with_columns(max_support=pl.max_horizontal(*supports))
    tie_expressions = []
    for left_index, left_base_name in enumerate(BASES):
        for right_index, right_base_name in enumerate(BASES):
            if left_index == right_index:
                continue
            support = pl.col(left_base_name).cast(pl.Float64) * pl.col(
                f"{right_base_name}_right"
            )
            tie_expressions.append(
                ((support == pl.col("max_support")) & (support > 0)).cast(pl.UInt8)
            )
    snv_sites = snv_sites.with_columns(
        tied_candidates=pl.sum_horizontal(*tie_expressions)
    )
    tie_sites = snv_sites.filter(pl.col("tied_candidates") > 1).select(
        "gene_id", "codon_id", "codon_offset"
    )

    candidates = []
    for left_index, left_base_name in enumerate(BASES):
        for right_index, right_base_name in enumerate(BASES):
            if left_index == right_index:
                continue
            support = pl.col(left_base_name).cast(pl.Float64) * pl.col(
                f"{right_base_name}_right"
            )
            candidate = snv_sites.filter(support > 0)
            if allele_integration == "consensus":
                candidate = candidate.filter(support == pl.col("max_support"))
                weight = 1.0 / pl.col("tied_candidates")
            else:
                weight = (
                    pl.col(left_base_name).cast(pl.Float64) / pl.col("coverage")
                ) * (
                    pl.col(f"{right_base_name}_right").cast(pl.Float64)
                    / pl.col("coverage_right")
                )
            candidates.append(
                _event_frame(
                    candidate,
                    source_base=pl.lit(left_index, dtype=pl.UInt8),
                    target_base=pl.lit(right_index, dtype=pl.UInt8),
                    weight=weight,
                )
            )
    return sns_events, pl.concat(candidates), sns_sites, snv_sites, tie_sites


def _event_totals(events: pl.LazyFrame, category: str) -> pl.LazyFrame:
    annotated = events.join(
        codon_pair_lookup().lazy().select(
            "pair_code", "synonymous_changes", "nonsynonymous_changes",
            "valid_dnds", "stop_change",
        ),
        on="pair_code",
        how="left",
    )
    valid = annotated.filter(pl.col("valid_dnds"))
    changes = valid.group_by("gene_id").agg(
        (pl.col("synonymous_changes") * pl.col("weight"))
        .sum()
        .alias(f"{category}_synonymous_changes"),
        (pl.col("nonsynonymous_changes") * pl.col("weight"))
        .sum()
        .alias(f"{category}_nonsynonymous_changes"),
    )
    stops = annotated.group_by("gene_id").agg(
        (pl.col("stop_change") * pl.col("weight"))
        .sum()
        .alias(f"{category}_stop_changes")
    )
    return changes.join(stops, on="gene_id", how="full", coalesce=True)


def _summarize_events(
    codons: pl.LazyFrame,
    sns_events: pl.LazyFrame,
    snv_events: pl.LazyFrame,
    sns_sites: pl.LazyFrame,
    snv_sites: pl.LazyFrame,
    tie_sites: pl.LazyFrame,
    *,
    prefix: str = "",
) -> pl.LazyFrame:
    denominators = (
        codons.select("gene_id", "codon_id", "reference_codon_code")
        .unique("codon_id")
        .join(codon_site_lookup().lazy(), on="reference_codon_code", how="left")
        .filter(pl.col("valid_reference"))
        .group_by("gene_id")
        .agg(
            pl.len().cast(pl.Int64).alias("callable_codons"),
            pl.col("synonymous_sites").sum(),
            pl.col("nonsynonymous_sites").sum(),
        )
    )
    site_counts = sns_sites.group_by("gene_id").agg(
        pl.len().cast(pl.Int64).alias("sns_count")
    ).join(
        snv_sites.group_by("gene_id").agg(pl.len().cast(pl.Int64).alias("snv_count")),
        on="gene_id",
        how="full",
        coalesce=True,
    )
    tie_counts = tie_sites.group_by("gene_id").agg(
        pl.len().cast(pl.Int64).alias("allele_tie_sites")
    )
    result = (
        denominators
        .join(site_counts, on="gene_id", how="left")
        .join(_event_totals(sns_events, "sns"), on="gene_id", how="left")
        .join(_event_totals(snv_events, "snv"), on="gene_id", how="left")
        .join(tie_counts, on="gene_id", how="left")
        .with_columns(
            pl.col("sns_count", "snv_count", "allele_tie_sites")
            .fill_null(0)
            .cast(pl.Int64),
            pl.col(
                "sns_synonymous_changes", "sns_nonsynonymous_changes",
                "snv_synonymous_changes", "snv_nonsynonymous_changes",
                "sns_stop_changes", "snv_stop_changes",
            ).fill_null(0.0).cast(pl.Float64),
        )
    )

    def proportion(changes: str, sites: str) -> pl.Expr:
        return pl.when(pl.col(sites) > 0).then(pl.col(changes) / pl.col(sites)).otherwise(None)

    def ratio(numerator: str, denominator: str) -> pl.Expr:
        return pl.when(pl.col(denominator) > 0).then(
            pl.col(numerator) / pl.col(denominator)
        ).otherwise(None)

    result = (
        result.with_columns(
            proportion("snv_synonymous_changes", "synonymous_sites").alias("pS"),
            proportion("snv_nonsynonymous_changes", "nonsynonymous_sites").alias("pN"),
            proportion("sns_synonymous_changes", "synonymous_sites").alias("dS"),
            proportion("sns_nonsynonymous_changes", "nonsynonymous_sites").alias("dN"),
        )
        .with_columns(
            ratio("pN", "pS").alias("pN_pS"),
        )
        .with_columns(ratio("dN", "dS").alias("dN_dS"))
    )
    if prefix:
        result = result.rename(
            {column: f"{prefix}{column}" for column in DNDS_RESULT_COLUMNS}
        )
    return result


def reference_dnds(
    codon_profile: Union[str, Path],
    gene_info: Union[str, Path],
    *,
    min_cov: int = 5,
    allele_integration: str = "consensus",
) -> pl.LazyFrame:
    """Calculate reference-relative SNS dN/dS and SNV pN/pS by gene."""
    allele_integration = validate_allele_integration(allele_integration)
    codons = pl.scan_parquet(codon_profile).filter(pl.col("min_coverage") >= min_cov)
    event_frames = _reference_event_frames(codons, allele_integration)
    summary = _summarize_events(codons, *event_frames, prefix="ref_").collect(
        optimizations=DNDS_QUERY_OPTIMIZATIONS
    ).lazy()
    return (
        pl.scan_parquet(gene_info)
        .select("gene_id", "genome", "gene")
        .join(summary, on="gene_id", how="left")
        .drop("gene_id")
    )


def validate_codon_profiles(
    codon_profile_1: Union[str, Path],
    codon_profile_2: Union[str, Path],
) -> None:
    """Require matching codon schema and gene-information contracts."""
    from zipstrain import utils

    metadata_1 = utils._read_custom_parquet_metadata(codon_profile_1)
    metadata_2 = utils._read_custom_parquet_metadata(codon_profile_2)
    for path, metadata in (
        (codon_profile_1, metadata_1),
        (codon_profile_2, metadata_2),
    ):
        if metadata.get(CODON_SCHEMA_METADATA_KEY) != CODON_PROFILE_SCHEMA_VERSION:
            raise ValueError(
                f"Codon profile uses an unsupported schema: {path}. "
                "Re-run profiling with the current ZipStrain version."
            )
    for key in (CODON_SCHEMA_METADATA_KEY, GENE_INFO_HASH_METADATA_KEY):
        value_1 = metadata_1.get(key)
        value_2 = metadata_2.get(key)
        if value_1 is None or value_2 is None:
            raise ValueError(f"Codon profile metadata is missing required field {key!r}.")
        if value_1 != value_2:
            raise ValueError(f"Codon profiles have incompatible {key!r} metadata.")


def pairwise_dnds(
    codon_profile_1: Union[str, Path],
    codon_profile_2: Union[str, Path],
    gene_info: Union[str, Path],
    *,
    min_cov: int = 5,
    allele_integration: str = "consensus",
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> pl.LazyFrame:
    """Calculate pairwise SNS dN/dS and SNV pN/pS from codon profiles."""
    from zipstrain import utils

    validate_codon_profiles(codon_profile_1, codon_profile_2)
    expected_gene_hash = utils.sha256_file(gene_info)
    metadata = utils._read_custom_parquet_metadata(codon_profile_1)
    if metadata.get(GENE_INFO_HASH_METADATA_KEY) != expected_gene_hash:
        raise ValueError(
            "The supplied gene information table does not match the codon profiles."
        )
    allele_integration = validate_allele_integration(allele_integration)
    c1 = pl.scan_parquet(codon_profile_1).filter(
        pl.col("min_coverage") >= min_cov
    ).set_sorted("codon_id")
    c2 = pl.scan_parquet(codon_profile_2).filter(
        pl.col("min_coverage") >= min_cov
    ).select(
        "codon_id",
        *(pl.col(column).alias(f"{column}_right") for column in CODON_PROFILE_COLUMNS[4:]),
    ).set_sorted("codon_id")
    # codon_id is database-wide and unique, so sorted inputs permit Polars to
    # use a merge join rather than building a large hash table.
    codons = c1.join(c2, on="codon_id", how="inner", maintain_order="left")
    event_frames = _pairwise_event_frames(codons, allele_integration)
    summary = _summarize_events(codons, *event_frames).collect(
        optimizations=DNDS_QUERY_OPTIMIZATIONS
    ).lazy()
    genes = pl.scan_parquet(gene_info).select("gene_id", "genome", "gene")
    if genome_scope != "all":
        genes = genes.filter(pl.col("genome") == genome_scope)
    if gene_scope != "all":
        genes = genes.filter(pl.col("gene") == gene_scope)
    return genes.join(summary, on="gene_id", how="inner").drop("gene_id")
