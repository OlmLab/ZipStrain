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
CODON_PROFILE_SCHEMA_VERSION = "1"
CODON_SCHEMA_METADATA_KEY = "zipstrain_codon_profile_schema"
GENE_INFO_HASH_METADATA_KEY = "zipstrain_gene_info_hash"
DEFAULT_MEMORY_LIMIT = "1GB"

CODON_PROFILE_COLUMNS = (
    "gene_id",
    "codon_id",
    "codon_code",
    "reference_codon_code",
    "min_coverage",
    "min_major_freq",
)

DNDS_RESULT_COLUMNS = (
    "callable_codons",
    "synonymous_changes",
    "nonsynonymous_changes",
    "synonymous_sites",
    "nonsynonymous_sites",
    "dS",
    "dN",
    "dN_dS",
    "stop_changes",
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


def _quote(value: Union[str, Path]) -> str:
    return str(value).replace("'", "''")


def write_codon_profile(
    *,
    profile_path: Union[str, Path],
    gene_info_path: Union[str, Path],
    output_path: Union[str, Path],
    temp_directory: Union[str, Path],
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
    metadata: Optional[dict[str, str]] = None,
) -> None:
    """Build a sparse, coordinate-independent codon summary out of core."""
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
                  GREATEST(p.A, p.C, p.G, p.T)::BIGINT AS major_count,
                  CASE
                    WHEN p.A = GREATEST(p.A, p.C, p.G, p.T) THEN 0
                    WHEN p.C = GREATEST(p.A, p.C, p.G, p.T) THEN 1
                    WHEN p.G = GREATEST(p.A, p.C, p.G, p.T) THEN 2
                    ELSE 3
                  END::INTEGER AS genomic_base,
                  CASE p.ref_base_bitmask
                    WHEN 1 THEN 0 WHEN 2 THEN 1 WHEN 4 THEN 2 WHEN 8 THEN 3
                    ELSE NULL
                  END::INTEGER AS genomic_reference_base,
                  ((p.A = GREATEST(p.A, p.C, p.G, p.T))::INTEGER
                   + (p.C = GREATEST(p.A, p.C, p.G, p.T))::INTEGER
                   + (p.G = GREATEST(p.A, p.C, p.G, p.T))::INTEGER
                   + (p.T = GREATEST(p.A, p.C, p.G, p.T))::INTEGER) AS major_ties
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
                  CASE WHEN strand = 1 THEN genomic_base ELSE 3 - genomic_base END AS coding_base,
                  CASE WHEN strand = 1
                    THEN genomic_reference_base
                    ELSE 3 - genomic_reference_base
                  END AS coding_reference_base
                FROM mapped
                WHERE coverage > 0
                  AND genomic_reference_base IS NOT NULL
                  AND major_ties = 1
              ), coded AS (
                SELECT
                  gene_id,
                  first_codon_id + FLOOR(coding_index / 3)::UBIGINT AS codon_id,
                  (coding_index % 3)::INTEGER AS codon_offset,
                  coding_base,
                  coding_reference_base,
                  coverage,
                  major_count::DOUBLE / coverage AS major_freq
                FROM oriented
                WHERE coding_index >= 0
                  AND coding_index < n_codons * 3
              )
              SELECT
                gene_id,
                codon_id,
                SUM(coding_base * CASE codon_offset WHEN 0 THEN 16 WHEN 1 THEN 4 ELSE 1 END)::UTINYINT AS codon_code,
                SUM(coding_reference_base * CASE codon_offset WHEN 0 THEN 16 WHEN 1 THEN 4 ELSE 1 END)::UTINYINT AS reference_codon_code,
                MIN(coverage)::UINTEGER AS min_coverage,
                MIN(major_freq)::FLOAT AS min_major_freq
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


def _summarize_pairs(
    frame: pl.LazyFrame,
    *,
    group_column: str = "gene_id",
    prefix: str = "",
) -> pl.LazyFrame:
    lookup = codon_pair_lookup().lazy()
    annotated = frame.join(lookup, on="pair_code", how="left")
    valid = annotated.filter(pl.col("valid_dnds"))
    totals = valid.group_by(group_column).agg(
        pl.len().cast(pl.Int64).alias(f"{prefix}callable_codons"),
        pl.col("synonymous_changes").sum().alias(f"{prefix}synonymous_changes"),
        pl.col("nonsynonymous_changes").sum().alias(f"{prefix}nonsynonymous_changes"),
        pl.col("synonymous_sites").sum().alias(f"{prefix}synonymous_sites"),
        pl.col("nonsynonymous_sites").sum().alias(f"{prefix}nonsynonymous_sites"),
    )
    stops = annotated.group_by(group_column).agg(
        pl.col("stop_change").sum().cast(pl.Int64).alias(f"{prefix}stop_changes")
    )
    return totals.join(stops, on=group_column, how="full", coalesce=True).with_columns(
        pl.when(pl.col(f"{prefix}synonymous_sites") > 0)
        .then(pl.col(f"{prefix}synonymous_changes") / pl.col(f"{prefix}synonymous_sites"))
        .otherwise(None)
        .cast(pl.Float64)
        .alias(f"{prefix}dS"),
        pl.when(pl.col(f"{prefix}nonsynonymous_sites") > 0)
        .then(pl.col(f"{prefix}nonsynonymous_changes") / pl.col(f"{prefix}nonsynonymous_sites"))
        .otherwise(None)
        .cast(pl.Float64)
        .alias(f"{prefix}dN"),
    ).with_columns(
        pl.when(pl.col(f"{prefix}dS") > 0)
        .then(pl.col(f"{prefix}dN") / pl.col(f"{prefix}dS"))
        .otherwise(None)
        .cast(pl.Float64)
        .alias(f"{prefix}dN_dS")
    )


def reference_dnds(
    codon_profile: Union[str, Path],
    gene_info: Union[str, Path],
    *,
    min_cov: int = 5,
    min_major_freq: float = 0.0,
) -> pl.LazyFrame:
    """Calculate sample-versus-reference dN/dS by gene."""
    codons = pl.scan_parquet(codon_profile).filter(
        (pl.col("min_coverage") >= min_cov)
        & (pl.col("min_major_freq") >= min_major_freq)
    ).with_columns(
        (
            pl.col("reference_codon_code").cast(pl.UInt16) * 64
            + pl.col("codon_code").cast(pl.UInt16)
        ).alias("pair_code")
    )
    summary = _summarize_pairs(codons, prefix="ref_")
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
    min_major_freq: float = 0.0,
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> pl.LazyFrame:
    """Calculate sample-to-sample dN/dS from two sparse codon profiles."""
    from zipstrain import utils

    validate_codon_profiles(codon_profile_1, codon_profile_2)
    expected_gene_hash = utils.sha256_file(gene_info)
    metadata = utils._read_custom_parquet_metadata(codon_profile_1)
    if metadata.get(GENE_INFO_HASH_METADATA_KEY) != expected_gene_hash:
        raise ValueError(
            "The supplied gene information table does not match the codon profiles."
        )
    c1 = (
        pl.scan_parquet(codon_profile_1)
        .filter(
            (pl.col("min_coverage") >= min_cov)
            & (pl.col("min_major_freq") >= min_major_freq)
        )
        .select("gene_id", "codon_id", "codon_code")
        .set_sorted("codon_id")
    )
    c2 = (
        pl.scan_parquet(codon_profile_2)
        .filter(
            (pl.col("min_coverage") >= min_cov)
            & (pl.col("min_major_freq") >= min_major_freq)
        )
        .select("codon_id", pl.col("codon_code").alias("codon_code_2"))
        .set_sorted("codon_id")
    )
    # codon_id is database-wide and unique, so sorted inputs permit Polars to
    # use a merge join rather than building a large hash table.
    pairs = c1.join(c2, on="codon_id", how="inner", maintain_order="left").with_columns(
        (
            pl.col("codon_code").cast(pl.UInt16) * 64
            + pl.col("codon_code_2").cast(pl.UInt16)
        ).alias("pair_code")
    )
    summary = _summarize_pairs(pairs)
    genes = pl.scan_parquet(gene_info).select("gene_id", "genome", "gene")
    if genome_scope != "all":
        genes = genes.filter(pl.col("genome") == genome_scope)
    if gene_scope != "all":
        genes = genes.filter(pl.col("gene") == gene_scope)
    return genes.join(summary, on="gene_id", how="inner").drop("gene_id")
