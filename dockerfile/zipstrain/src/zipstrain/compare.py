"""zipstrain.compare
========================
This module provides all comparison functions for zipstrain.

"""

from pathlib import Path
from typing import Iterable, Literal, Optional, Union

import polars as pl
import duckdb


GENOME_COMPARISON_CALCULATIONS = ("ani", "ibs", "identical_genes")
GENOME_COMPARISON_DEFAULT_CALCULATIONS = GENOME_COMPARISON_CALCULATIONS
GENOME_COMPARISON_CALCULATION_ALIASES = {
    "ani": "ani",
    "popani": "ani",
    "ibs": "ibs",
    "max_block": "ibs",
    "max_consecutive_length": "ibs",
    "identical_genes": "identical_genes",
    "gene_identity": "identical_genes",
    "id_genes": "identical_genes",
}


GenomeCalculate = Literal["all", "ani"]


def parse_genome_calculations(calculate: Optional[Union[str, Iterable[str]]] = None) -> tuple[str, ...]:
    """Parse and normalize genome metric selection tokens.

    Accepted input formats:
      - None -> default ("ani", "ibs", "identical_genes")
      - "ani+ibs+identical_genes"
      - "all"
      - iterable of token strings
    """
    if calculate is None:
        return GENOME_COMPARISON_DEFAULT_CALCULATIONS

    raw_tokens: list[str] = []
    if isinstance(calculate, str):
        for plus_part in calculate.split("+"):
            for comma_part in plus_part.split(","):
                token = comma_part.strip().lower()
                if token:
                    raw_tokens.append(token)
    else:
        raw_tokens = [str(token).strip().lower() for token in calculate if str(token).strip()]

    if not raw_tokens:
        return tuple()
    if "all" in raw_tokens:
        return GENOME_COMPARISON_CALCULATIONS

    normalized: set[str] = set()
    for token in raw_tokens:
        mapped = GENOME_COMPARISON_CALCULATION_ALIASES.get(token)
        if mapped is None:
            supported = "all," + ",".join(GENOME_COMPARISON_CALCULATIONS)
            raise ValueError(f"Unsupported calculation '{token}'. Supported values: {supported}")
        normalized.add(mapped)
    return tuple(metric for metric in GENOME_COMPARISON_CALCULATIONS if metric in normalized)


def _normalize_genome_calculate(calculate: Optional[Union[str, Iterable[str]]] = None) -> GenomeCalculate:
    """Compatibility helper for legacy all-vs-ani code paths."""
    calculations = parse_genome_calculations(calculate)
    if calculations == ("ani",):
        return "ani"
    return "all"


def genome_metric_output_columns(calculate: Optional[Union[str, Iterable[str]]] = None) -> list[str]:
    """Return ordered output columns for selected genome-level calculations."""
    calculations = parse_genome_calculations(calculate)
    cols = ["genome"]
    if "ani" in calculations:
        cols.extend(["total_positions", "share_allele_pos", "genome_pop_ani"])
    if "ibs" in calculations:
        cols.append("max_consecutive_length")
    if "identical_genes" in calculations:
        cols.extend(["shared_genes_count", "identical_gene_count", "perc_id_genes"])
    return cols


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
    MPILE_1_BASES = ["A", "T", "C", "G"]
    MPILE_2_BASES = ["A_2", "T_2", "C_2", "G_2"]

    def popani(self):
        return pl.col("A")*pl.col("A_2") + pl.col("C")*pl.col("C_2") + pl.col("G")*pl.col("G_2") + pl.col("T")*pl.col("T_2")
    
    def conani(self):
        max_base_1=pl.max_horizontal(*[pl.col(base) for base in self.MPILE_1_BASES])
        max_base_2=pl.max_horizontal(*[pl.col(base) for base in self.MPILE_2_BASES])
        return pl.when((pl.col("A")==max_base_1) & (pl.col("A_2")==max_base_2) | 
                       (pl.col("T")==max_base_1) & (pl.col("T_2")==max_base_2) | 
                       (pl.col("C")==max_base_1) & (pl.col("C_2")==max_base_2) | 
                       (pl.col("G")==max_base_1) & (pl.col("G_2")==max_base_2)).then(1).otherwise(0)
    
    def generalized_cos_ani(self,threshold:float=0.4):
        dot_product = pl.col("A")*pl.col("A_2") + pl.col("C")*pl.col("C_2") + pl.col("G")*pl.col("G_2") + pl.col("T")*pl.col("T_2")
        magnitude_1 = (pl.col("A")**2 + pl.col("C")**2 + pl.col("G")**2 + pl.col("T")**2)**0.5
        magnitude_2 = (pl.col("A_2")**2 + pl.col("C_2")**2 + pl.col("G_2")**2 + pl.col("T_2")**2)**0.5
        cos_sim = dot_product / (magnitude_1 * magnitude_2)
        return pl.when(cos_sim >= threshold).then(1).otherwise(0)

    def __getattribute__(self, name):
        if name.startswith("cosani_"):
            try:
                threshold = float(name.split("_")[1])
            except ValueError:
                raise AttributeError(f"Invalid threshold in method name: {name}")
            return lambda: self.generalized_cos_ani(threshold)
        else:
            return super().__getattribute__(name)

def _duckdb_quote_sql_string(value: str) -> str:
    return value.replace("'", "''")


def _duckdb_from_source(conn: duckdb.DuckDBPyConnection, source: Union[str, Path, pl.LazyFrame], view_name: str) -> str:
    if isinstance(source, pl.LazyFrame):
        conn.register(view_name, source.collect().to_arrow())
        return view_name
    if isinstance(source, (str, Path)):
        path = _duckdb_quote_sql_string(str(source))
        return f"read_parquet('{path}')"
    raise TypeError(f"Unsupported source type for DuckDB input: {type(source)}")


def _duckdb_configure_connection(
    conn: duckdb.DuckDBPyConnection,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> None:
    if memory_limit:
        conn.execute(f"SET memory_limit='{_duckdb_quote_sql_string(memory_limit)}'")
    if temp_directory:
        conn.execute(f"SET temp_directory='{_duckdb_quote_sql_string(str(temp_directory))}'")
    if threads is not None:
        if threads < 1:
            raise ValueError("threads must be >= 1")
        conn.execute(f"SET threads={int(threads)}")


def _duckdb_ani_expression(ani_method: str) -> str:
    if ani_method == "popani":
        # Cast at expression time (post-join) to avoid upfront casting in filtered scans.
        return (
            "CAST(p1.A AS INT32)*CAST(p2.A AS INT32) + "
            "CAST(p1.C AS INT32)*CAST(p2.C AS INT32) + "
            "CAST(p1.G AS INT32)*CAST(p2.G AS INT32) + "
            "CAST(p1.T AS INT32)*CAST(p2.T AS INT32)"
        )
    if ani_method == "conani":
        return """
        CASE WHEN (
          (p1.A=GREATEST(p1.A,p1.T,p1.C,p1.G) AND p2.A=GREATEST(p2.A,p2.T,p2.C,p2.G)) OR
          (p1.T=GREATEST(p1.A,p1.T,p1.C,p1.G) AND p2.T=GREATEST(p2.A,p2.T,p2.C,p2.G)) OR
          (p1.C=GREATEST(p1.A,p1.T,p1.C,p1.G) AND p2.C=GREATEST(p2.A,p2.T,p2.C,p2.G)) OR
          (p1.G=GREATEST(p1.A,p1.T,p1.C,p1.G) AND p2.G=GREATEST(p2.A,p2.T,p2.C,p2.G))
        ) THEN 1 ELSE 0 END
        """
    if ani_method.startswith("cosani_"):
        thr = float(ani_method.split("_")[1])
        return f"""
        CASE WHEN (
          (
            CAST(p1.A AS DOUBLE)*CAST(p2.A AS DOUBLE) +
            CAST(p1.C AS DOUBLE)*CAST(p2.C AS DOUBLE) +
            CAST(p1.G AS DOUBLE)*CAST(p2.G AS DOUBLE) +
            CAST(p1.T AS DOUBLE)*CAST(p2.T AS DOUBLE)
          ) /
          NULLIF(
            SQRT(
              POWER(CAST(p1.A AS DOUBLE),2) +
              POWER(CAST(p1.C AS DOUBLE),2) +
              POWER(CAST(p1.G AS DOUBLE),2) +
              POWER(CAST(p1.T AS DOUBLE),2)
            ) *
            SQRT(
              POWER(CAST(p2.A AS DOUBLE),2) +
              POWER(CAST(p2.C AS DOUBLE),2) +
              POWER(CAST(p2.G AS DOUBLE),2) +
              POWER(CAST(p2.T AS DOUBLE),2)
            ),
            0
          )
        ) >= {thr} THEN 1 ELSE 0 END
        """
    raise ValueError(f"Unsupported ani_method: {ani_method}")


def _duckdb_shared_query(
    p1_source: str,
    p2_source: str,
    min_cov: int,
    genome_scope: str,
    ani_method: str,
    gene_scope: str = "all",
) -> str:
    ani_expr = _duckdb_ani_expression(ani_method)
    genome_scope_sql = _duckdb_quote_sql_string(genome_scope)
    gene_scope_sql = _duckdb_quote_sql_string(gene_scope)
    return f"""
    WITH p1f AS (
      SELECT
        chrom,
        pos,
        gene,
        genome,
        A,
        T,
        C,
        G
      FROM {p1_source} p1
      WHERE (A + T + C + G) >= {min_cov}
        AND ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
        AND ('{gene_scope_sql}' = 'all' OR gene = '{gene_scope_sql}')
    ),
    p2f AS (
      SELECT
        chrom,
        pos,
        gene,
        genome,
        A,
        T,
        C,
        G
      FROM {p2_source} p2
      WHERE (A + T + C + G) >= {min_cov}
        AND ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
        AND ('{gene_scope_sql}' = 'all' OR gene = '{gene_scope_sql}')
    )
    SELECT
      {ani_expr} AS surr,
      p1.chrom AS scaffold,
      p1.pos,
      p1.gene,
      p1.genome
    FROM p1f p1
    INNER JOIN p2f p2
      ON p1.genome = p2.genome AND p1.chrom = p2.chrom AND p1.pos = p2.pos
    """


def _duckdb_copy_query_to_parquet(conn: duckdb.DuckDBPyConnection, query: str, output_file: Union[str, Path]) -> None:
    output_sql = _duckdb_quote_sql_string(str(output_file))
    conn.execute(f"COPY ({query}) TO '{output_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)")


def _duckdb_build_query_with_ctes(ctes: list[str], final_select: str) -> str:
    return "WITH\n" + ",\n".join(ctes) + "\n" + final_select


def _finalize_genome_compare_output(genome_comp: pl.LazyFrame, calculate: GenomeCalculate) -> pl.LazyFrame:
    if calculate == "ani":
        return genome_comp.with_columns(
            pl.col("total_positions").fill_null(0).cast(pl.Int64),
            pl.col("share_allele_pos").fill_null(0).cast(pl.Int64),
            pl.col("genome_pop_ani").fill_null(0.0).cast(pl.Float64),
        ).select(
            "genome",
            "total_positions",
            "share_allele_pos",
            "genome_pop_ani",
        )

    return genome_comp.with_columns(
        pl.col("total_positions").fill_null(0).cast(pl.Int64),
        pl.col("share_allele_pos").fill_null(0).cast(pl.Int64),
        pl.col("genome_pop_ani").fill_null(0.0).cast(pl.Float64),
        pl.col("max_consecutive_length").fill_null(0).cast(pl.Int64),
        pl.col("shared_genes_count").fill_null(0).cast(pl.Int64),
        pl.col("identical_gene_count").fill_null(0).cast(pl.Int64),
        pl.col("perc_id_genes").fill_null(0.0).cast(pl.Float64),
    ).select(
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
        "max_consecutive_length",
        "shared_genes_count",
        "identical_gene_count",
        "perc_id_genes",
    )


def _join_all_requested_genomes(
    genome_comp: pl.LazyFrame,
    stb_file: Optional[Union[str, Path]],
    genome_scope: str,
) -> pl.LazyFrame:
    if stb_file is None:
        return genome_comp

    genomes_utf8 = (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(pl.col("column_2").cast(pl.Utf8).alias("genome"))
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


def _duckdb_genomes_scope_cte(stb_sql: str, genome_scope_sql: str) -> str:
    return f"""
genomes AS (
  SELECT DISTINCT CAST(column1 AS VARCHAR) AS genome
  FROM read_csv('{stb_sql}', delim='\\t', header=false)
  WHERE ('{genome_scope_sql}' = 'all' OR CAST(column1 AS VARCHAR) = '{genome_scope_sql}')
)""".strip()


def _duckdb_contig_pop_max_ctes(
    shared_source: str = "shared",
    include_pop: bool = True,
    include_max_blocks: bool = True,
) -> list[str]:
    ctes: list[str] = []
    if not include_pop and not include_max_blocks:
        return ctes
    ctes.extend(
        [
            f"""
contig_base AS (
  SELECT
    s.*,
    LAG(scaffold) OVER (PARTITION BY genome ORDER BY scaffold, pos) AS prev_scaffold
  FROM {shared_source} s
)""".strip(),
            """
contig AS (
  SELECT
    *,
    SUM(
      CASE
        WHEN prev_scaffold IS NULL OR scaffold != prev_scaffold OR surr = 0 THEN 1
        ELSE 0
      END
    ) OVER (PARTITION BY genome ORDER BY scaffold, pos ROWS UNBOUNDED PRECEDING) AS group_id
  FROM contig_base
)""".strip(),
        ]
    )
    if include_pop:
        ctes.append(
            """
pop AS (
  SELECT
    genome,
    COUNT(*)::BIGINT AS total_positions,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END) * 100.0 / NULLIF(COUNT(*), 0) AS genome_pop_ani
  FROM contig
  GROUP BY genome
)""".strip()
        )
    if include_max_blocks:
        ctes.append(
            """
max_blocks AS (
  SELECT
    genome,
    MAX(length)::BIGINT AS max_consecutive_length
  FROM (
    SELECT genome, scaffold, group_id, COUNT(*)::BIGINT AS length
    FROM contig
    GROUP BY genome, scaffold, group_id
  ) blocks
  GROUP BY genome
)""".strip(),
        )
    return ctes


def _duckdb_gene_stats_ctes(min_gene_compare_len: int, contig_source: str = "contig") -> list[str]:
    return [
        f"""
gene_base AS (
  SELECT
    genome,
    gene,
    COUNT(*)::BIGINT AS total_positions,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos
  FROM {contig_source}
  WHERE gene != 'NA'
  GROUP BY genome, gene
  HAVING COUNT(*) >= {min_gene_compare_len}
)""".strip(),
        """
gene_stats AS (
  SELECT
    genome,
    COUNT(*)::BIGINT AS shared_genes_count,
    SUM(CASE WHEN share_allele_pos = total_positions THEN 1 ELSE 0 END)::BIGINT AS identical_gene_count,
    SUM(CASE WHEN share_allele_pos = total_positions THEN 1 ELSE 0 END) * 100.0 / NULLIF(COUNT(*), 0) AS perc_id_genes
  FROM gene_base
  GROUP BY genome
)""".strip(),
    ]


def _duckdb_merged_genome_stats_cte() -> str:
    return """
merged AS (
  SELECT
    g.genome,
    COALESCE(p.total_positions, 0)::BIGINT AS total_positions,
    COALESCE(p.share_allele_pos, 0)::BIGINT AS share_allele_pos,
    COALESCE(p.genome_pop_ani, 0.0)::DOUBLE AS genome_pop_ani,
    COALESCE(m.max_consecutive_length, 0)::BIGINT AS max_consecutive_length,
    COALESCE(gs.shared_genes_count, 0)::BIGINT AS shared_genes_count,
    COALESCE(gs.identical_gene_count, 0)::BIGINT AS identical_gene_count,
    COALESCE(gs.perc_id_genes, 0.0)::DOUBLE AS perc_id_genes
  FROM genomes g
  LEFT JOIN pop p ON g.genome = p.genome
  LEFT JOIN max_blocks m ON g.genome = m.genome
  LEFT JOIN gene_stats gs ON g.genome = gs.genome
)""".strip()


def _duckdb_genome_stats_select() -> str:
    return """
SELECT
  p.genome,
  p.total_positions,
  p.share_allele_pos,
  p.genome_pop_ani,
  COALESCE(m.max_consecutive_length, 0)::BIGINT AS max_consecutive_length
FROM pop p
LEFT JOIN max_blocks m ON p.genome = m.genome
ORDER BY p.genome
""".strip()


def _duckdb_genome_output_select(sample_1_sql: str, sample_2_sql: str) -> str:
    return f"""
SELECT
  genome,
  total_positions,
  share_allele_pos,
  genome_pop_ani,
  max_consecutive_length,
  shared_genes_count,
  identical_gene_count,
  perc_id_genes,
  '{sample_1_sql}' AS sample_1,
  '{sample_2_sql}' AS sample_2
FROM merged
ORDER BY genome
""".strip()


def _duckdb_genomes_from_shared_cte(genome_scope_sql: str) -> str:
    return f"""
genomes AS (
  SELECT DISTINCT CAST(genome AS VARCHAR) AS genome
  FROM shared
  WHERE ('{genome_scope_sql}' = 'all' OR CAST(genome AS VARCHAR) = '{genome_scope_sql}')
)""".strip()


def _duckdb_genome_compare_query(
    *,
    shared_query: str,
    min_gene_compare_len: int,
    genome_scope: str,
    calculate: Optional[Union[str, Iterable[str]]] = None,
    stb_file: Optional[Union[str, Path]] = None,
    sample_1_name: Optional[str] = None,
    sample_2_name: Optional[str] = None,
) -> str:
    calculations = parse_genome_calculations(calculate)
    need_ani = "ani" in calculations
    need_ibs = "ibs" in calculations
    need_identical_genes = "identical_genes" in calculations

    genome_scope_sql = _duckdb_quote_sql_string(genome_scope)
    ctes = [f"shared AS (\n{shared_query}\n)"]
    if stb_file is None:
        ctes.append(_duckdb_genomes_from_shared_cte(genome_scope_sql=genome_scope_sql))
    else:
        stb_sql = _duckdb_quote_sql_string(str(stb_file))
        ctes.append(_duckdb_genomes_scope_cte(stb_sql=stb_sql, genome_scope_sql=genome_scope_sql))

    if need_ani or need_ibs:
        ctes.extend(
            _duckdb_contig_pop_max_ctes(
                shared_source="shared",
                include_pop=need_ani,
                include_max_blocks=need_ibs,
            )
        )
    if need_identical_genes:
        ctes.extend(_duckdb_gene_stats_ctes(min_gene_compare_len=min_gene_compare_len, contig_source="shared"))

    select_parts = ["g.genome AS genome"]
    joins: list[str] = []
    if need_ani:
        select_parts.extend(
            [
                "COALESCE(p.total_positions, 0)::BIGINT AS total_positions",
                "COALESCE(p.share_allele_pos, 0)::BIGINT AS share_allele_pos",
                "COALESCE(p.genome_pop_ani, 0.0)::DOUBLE AS genome_pop_ani",
            ]
        )
        joins.append("LEFT JOIN pop p ON g.genome = p.genome")
    if need_ibs:
        select_parts.append("COALESCE(m.max_consecutive_length, 0)::BIGINT AS max_consecutive_length")
        joins.append("LEFT JOIN max_blocks m ON g.genome = m.genome")
    if need_identical_genes:
        select_parts.extend(
            [
                "COALESCE(gs.shared_genes_count, 0)::BIGINT AS shared_genes_count",
                "COALESCE(gs.identical_gene_count, 0)::BIGINT AS identical_gene_count",
                "COALESCE(gs.perc_id_genes, 0.0)::DOUBLE AS perc_id_genes",
            ]
        )
        joins.append("LEFT JOIN gene_stats gs ON g.genome = gs.genome")
    if sample_1_name is not None and sample_2_name is not None:
        sample_1_sql = _duckdb_quote_sql_string(sample_1_name)
        sample_2_sql = _duckdb_quote_sql_string(sample_2_name)
        select_parts.extend([f"'{sample_1_sql}' AS sample_1", f"'{sample_2_sql}' AS sample_2"])

    final_select = "SELECT\n  " + ",\n  ".join(select_parts) + "\nFROM genomes g\n"
    if joins:
        final_select += "\n".join(joins) + "\n"
    final_select += "ORDER BY g.genome"
    return _duckdb_build_query_with_ctes(ctes, final_select)


def _as_lazy_profile(source: Union[str, Path, pl.LazyFrame]) -> pl.LazyFrame:
    if isinstance(source, pl.LazyFrame):
        return source
    if isinstance(source, (str, Path)):
        return pl.scan_parquet(source)
    raise TypeError(f"Unsupported source type for profile input: {type(source)}")


def _profile_scope_predicate(genome_scope: str = "all", gene_scope: str = "all") -> pl.Expr:
    expr = pl.lit(True)
    if genome_scope != "all":
        expr = expr & (pl.col("genome") == genome_scope)
    if gene_scope != "all":
        expr = expr & (pl.col("gene") == gene_scope)
    return expr


def _filter_profiles_polars(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    cov_expr = (pl.col("A") + pl.col("T") + pl.col("C") + pl.col("G")) >= min_cov
    scope_expr = _profile_scope_predicate(genome_scope=genome_scope, gene_scope=gene_scope)
    p1 = _as_lazy_profile(mpile1).filter(cov_expr & scope_expr)
    p2 = _as_lazy_profile(mpile2).filter(cov_expr & scope_expr)
    return p1, p2


def _shared_loci_polars(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
) -> pl.LazyFrame:
    ani_expr = getattr(PolarsANIExpressions(), ani_method)()
    p1, p2 = _filter_profiles_polars(
        mpile1=mpile1,
        mpile2=mpile2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
    )
    return (
        p1.join(
            p2,
            on=["genome", "chrom", "pos"],
            how="inner",
            suffix="_2",
        )
        .with_columns(ani_expr.alias("surr"))
        .select(
            pl.col("surr"),
            scaffold=pl.col("chrom"),
            pos=pl.col("pos"),
            gene=pl.col("gene"),
            genome=pl.col("genome"),
        )
    )


def duckdb_prefilter_by_scope(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    genome_scope: str = "all",
    gene_scope: str = "all",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Scope-filter both profiles in DuckDB and return in-memory LazyFrames."""
    scope_requested = genome_scope != "all" or gene_scope != "all"
    if not scope_requested:
        return _as_lazy_profile(mpile1), _as_lazy_profile(mpile2)

    con = duckdb.connect()
    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_prefilter")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_prefilter")
        genome_scope_sql = _duckdb_quote_sql_string(genome_scope)
        gene_scope_sql = _duckdb_quote_sql_string(gene_scope)
        where_sql = (
            f"('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}') "
            f"AND ('{gene_scope_sql}' = 'all' OR gene = '{gene_scope_sql}')"
        )
        p1_table = con.execute(
            f"""
            SELECT chrom, pos, gene, genome, A, T, C, G
            FROM {p1_source}
            WHERE {where_sql}
            """
        ).fetch_arrow_table()
        p2_table = con.execute(
            f"""
            SELECT chrom, pos, gene, genome, A, T, C, G
            FROM {p2_source}
            WHERE {where_sql}
            """
        ).fetch_arrow_table()
        return pl.from_arrow(p1_table).lazy(), pl.from_arrow(p2_table).lazy()
    finally:
        con.close()


def duckdb_filter_join(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",   # or specific genome
    ani_method: str = "popani",  # "popani" / "conani" / "cosani_0.4"
    gene_scope: str = "all",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> pl.LazyFrame:
    """Filter two profile sources and inner-join shared loci in DuckDB.

    Inputs can be parquet paths or Polars LazyFrames. Coverage, genome scope,
    and optional gene scope are pushed down into DuckDB. The returned lazy frame
    contains: `surr`, `scaffold`, `pos`, `gene`, and `genome`.
    """
    con = duckdb.connect()

    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_source")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_source")
        query = _duckdb_shared_query(
            p1_source=p1_source,
            p2_source=p2_source,
            min_cov=min_cov,
            genome_scope=genome_scope,
            ani_method=ani_method,
            gene_scope=gene_scope,
        )
        table = con.execute(query).fetch_arrow_table()
        return pl.from_arrow(table).lazy()
    finally:
        con.close()


def duckdb_filter_join_with_genome_stats(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",
    ani_method: str = "popani",
    gene_scope: str = "all",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Build shared loci and genome-level popANI/contiguity stats in DuckDB.

    Returns:
        tuple[pl.LazyFrame, pl.LazyFrame]:
            - shared loci projected for gene aggregation (`surr`, `gene`, `genome`)
            - genome stats (`genome`, `total_positions`, `share_allele_pos`,
              `genome_pop_ani`, `max_consecutive_length`)
    """
    con = duckdb.connect()
    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_source")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_source")
        shared_query = _duckdb_shared_query(
            p1_source=p1_source,
            p2_source=p2_source,
            min_cov=min_cov,
            genome_scope=genome_scope,
            ani_method=ani_method,
            gene_scope=gene_scope,
        )
        con.execute(f"CREATE TEMP TABLE shared AS {shared_query}")

        # Only transfer columns needed for downstream Polars gene aggregation.
        shared_table = con.execute(
            "SELECT surr, gene, genome FROM shared"
        ).fetch_arrow_table()

        genome_stats_query = _duckdb_build_query_with_ctes(
            _duckdb_contig_pop_max_ctes(shared_source="shared"),
            _duckdb_genome_stats_select(),
        )
        genome_stats_table = con.execute(genome_stats_query).fetch_arrow_table()
        return pl.from_arrow(shared_table).lazy(), pl.from_arrow(genome_stats_table).lazy()
    finally:
        con.close()


def duckdb_compare_genomes_to_parquet(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    output_file: Union[str, Path],
    stb_file: Union[str, Path],
    sample_1_name: str,
    sample_2_name: str,
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    calculate: Optional[Union[str, Iterable[str]]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> None:
    """Run genome comparison in DuckDB and write final output directly to parquet.

    This path avoids materializing large intermediate tables in Python memory.
    """
    con = duckdb.connect()
    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_source")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_source")
        shared_query = _duckdb_shared_query(
            p1_source=p1_source,
            p2_source=p2_source,
            min_cov=min_cov,
            genome_scope=genome_scope,
            ani_method=ani_method,
        )
        query = _duckdb_genome_compare_query(
            shared_query=shared_query,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            calculate=calculate,
            stb_file=stb_file,
            sample_1_name=sample_1_name,
            sample_2_name=sample_2_name,
        )
        _duckdb_copy_query_to_parquet(con, query, output_file)
    finally:
        con.close()


def duckdb_compare_genomes(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    calculate: Optional[Union[str, Iterable[str]]] = None,
    stb_file: Optional[Union[str, Path]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> pl.LazyFrame:
    """Run genome comparison in DuckDB and return selected metrics as a LazyFrame."""
    con = duckdb.connect()
    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_source")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_source")
        shared_query = _duckdb_shared_query(
            p1_source=p1_source,
            p2_source=p2_source,
            min_cov=min_cov,
            genome_scope=genome_scope,
            ani_method=ani_method,
        )
        query = _duckdb_genome_compare_query(
            shared_query=shared_query,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            calculate=calculate,
            stb_file=stb_file,
        )
        table = con.execute(query).fetch_arrow_table()
        return pl.from_arrow(table).lazy()
    finally:
        con.close()


def duckdb_compare_genes_to_parquet(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    output_file: Union[str, Path],
    sample_1_name: str,
    sample_2_name: str,
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> None:
    """Run gene comparison in DuckDB and write final output directly to parquet."""
    con = duckdb.connect()
    try:
        _duckdb_configure_connection(
            con,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
        p1_source = _duckdb_from_source(con, mpile1, "mpile1_source")
        p2_source = _duckdb_from_source(con, mpile2, "mpile2_source")
        shared_query = _duckdb_shared_query(
            p1_source=p1_source,
            p2_source=p2_source,
            min_cov=min_cov,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
            ani_method=ani_method,
        )
        sample_1_sql = _duckdb_quote_sql_string(sample_1_name)
        sample_2_sql = _duckdb_quote_sql_string(sample_2_name)
        query = f"""
        WITH shared AS (
          {shared_query}
        )
        SELECT
          genome,
          gene,
          COUNT(*)::BIGINT AS total_positions,
          SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos,
          SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END) * 100.0 / NULLIF(COUNT(*), 0) AS ani,
          '{sample_1_sql}' AS sample_1,
          '{sample_2_sql}' AS sample_2
        FROM shared
        GROUP BY genome, gene
        HAVING COUNT(*) >= {min_gene_compare_len}
        ORDER BY genome, gene
        """
        _duckdb_copy_query_to_parquet(con, query, output_file)
    finally:
        con.close()



def get_shared_locs(mpile_contig_1:pl.LazyFrame, mpile_contig_2:pl.LazyFrame,ani_method:str="popani") -> pl.LazyFrame:
    """
    Returns a lazyframe with ATCG information for shared scaffolds and positions between two mpileup files.

    Args:
        mpile_contig_1 (pl.LazyFrame): The first mpileup LazyFrame.
        mpile_contig_2 (pl.LazyFrame): The second mpileup LazyFrame.
        ani_method (str): The ANI calculation method to use. Default is "popani".
    
    Returns:
        pl.LazyFrame: Merged LazyFrame containing shared scaffolds and positions with ATCG information.
    """
    ani_expr=getattr(PolarsANIExpressions(), ani_method)()

    mpile_contig= mpile_contig_1.join(
        mpile_contig_2,
        on=["chrom", "pos"],
        how="inner",
        suffix="_2"  # To distinguish lf2 columns
    ).with_columns(
        ani_expr.alias("surr")
    ).select(
        pl.col("surr"),
        scaffold=pl.col("chrom"),
        pos=pl.col("pos"),
        gene=pl.col("gene")
    )
    return mpile_contig

def add_contiguity_info(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """ Adds group id information to the lazy frame. If on the same scaffold and not popANI, then they are in the same group.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with group id information added.
    """

    sort_cols = ["scaffold", "pos"]
    break_expr = (pl.col("scaffold") != pl.col("scaffold").shift(1).fill_null(pl.col("scaffold"))) | (pl.col("surr") == 0)
    if "genome" in mpile_contig.collect_schema().names():
        sort_cols = ["genome", "scaffold", "pos"]
        break_expr = (
            (pl.col("genome") != pl.col("genome").shift(1).fill_null(pl.col("genome")))
            | (pl.col("scaffold") != pl.col("scaffold").shift(1).fill_null(pl.col("scaffold")))
            | (pl.col("surr") == 0)
        )
    mpile_contig = mpile_contig.sort(sort_cols)
    mpile_contig = mpile_contig.with_columns(
        break_expr.cast(pl.Int64).cum_sum().alias("group_id")
    )
    return mpile_contig

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
    NOTE: Remember that this function should be applied to the merged mpileup using get_shared_locs.

    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with population ANI information added.
    """
    return mpile_contig.group_by("genome").agg(
            total_positions=pl.len(),
            share_allele_pos=(pl.col("surr") > 0 ).sum()
        ).with_columns(
            genome_pop_ani=pl.col("share_allele_pos")/pl.col("total_positions")*100,
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
    return block_lengths.group_by("genome").agg(pl.col("length").max().alias("max_consecutive_length"))

def get_gene_ani(mpile_contig:pl.LazyFrame, min_gene_compare_len:int) -> pl.LazyFrame:
    """
    Calculates gene ANI (Average Nucleotide Identity) for each gene in each genome.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
        min_gene_compare_len (int): Minimum length of the gene to consider for comparison.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with gene ANI information added.
    """
    return mpile_contig.select(["genome", "gene", "surr"]).filter(
        pl.col("gene") != "NA"
    ).group_by(["genome", "gene"]).agg(
        total_positions=pl.len(),
        share_allele_pos=(pl.col("surr") > 0).sum()
    ).filter(pl.col("total_positions") >= min_gene_compare_len).with_columns(
        identical=(pl.col("share_allele_pos") == pl.col("total_positions")),
    ).group_by("genome").agg(
        shared_genes_count=pl.len(),
        identical_gene_count=pl.col("identical").sum()
    ).with_columns(perc_id_genes=pl.col("identical_gene_count") / pl.col("shared_genes_count") * 100)

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
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
) -> pl.LazyFrame:
    """Compare two profiles fully in Polars and return genome-level statistics."""
    calculations = parse_genome_calculations(calculate)
    shared = _shared_loci_polars(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        ani_method=ani_method,
    )
    genome_comp_parts: list[pl.LazyFrame] = []
    if "ani" in calculations:
        genome_comp_parts.append(calculate_pop_ani(shared))
    if "ibs" in calculations:
        genome_comp_parts.append(get_longest_consecutive_blocks(add_contiguity_info(shared)))
    if "identical_genes" in calculations:
        genome_comp_parts.append(get_gene_ani(shared, min_gene_compare_len))

    if genome_comp_parts:
        genome_comp = genome_comp_parts[0]
        for part in genome_comp_parts[1:]:
            genome_comp = genome_comp.join(part, on="genome", how="left")
    else:
        genome_comp = shared.select("genome").unique()

    if stb_file is not None:
        genomes_utf8 = (
            pl.scan_csv(stb_file, separator="\t", has_header=False)
            .select(pl.col("column_2").cast(pl.Utf8).alias("genome"))
            .unique()
        )

    casts: list[pl.Expr] = []
    if "ani" in calculations:
        casts.extend(
            [
                pl.col("total_positions").fill_null(0).cast(pl.Int64),
                pl.col("share_allele_pos").fill_null(0).cast(pl.Int64),
                pl.col("genome_pop_ani").fill_null(0.0).cast(pl.Float64),
            ]
        )
    if "ibs" in calculations:
        casts.append(pl.col("max_consecutive_length").fill_null(0).cast(pl.Int64))
    if "identical_genes" in calculations:
        casts.extend(
            [
                pl.col("shared_genes_count").fill_null(0).cast(pl.Int64),
                pl.col("identical_gene_count").fill_null(0).cast(pl.Int64),
                pl.col("perc_id_genes").fill_null(0.0).cast(pl.Float64),
            ]
        )
    if casts:
        genome_comp = genome_comp.with_columns(casts)

    return genome_comp.select(genome_metric_output_columns(calculations))


def compare_genes_polars(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
) -> pl.LazyFrame:
    """Compare two profiles fully in Polars and return gene-level ANI statistics."""
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


def _compare_genomes_mixed(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    calculate: GenomeCalculate = "all",
    duckdb_memory_limit: Optional[str] = None,
    duckdb_temp_directory: Optional[Union[str, Path]] = None,
    duckdb_threads: Optional[int] = None,
) -> pl.LazyFrame:
    """DuckDB shared-loci + Polars gene aggregation path."""
    calculate = _normalize_genome_calculate(calculate)
    lf, genome_comp = duckdb_filter_join_with_genome_stats(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        ani_method=ani_method,
        genome_scope=genome_scope,
        memory_limit=duckdb_memory_limit,
        temp_directory=duckdb_temp_directory,
        threads=duckdb_threads,
    )
    if calculate == "ani":
        return _finalize_genome_compare_output(genome_comp, calculate)
    gene = get_gene_ani(lf, min_gene_compare_len)
    return _finalize_genome_compare_output(genome_comp.join(gene, on="genome", how="left"), calculate)


def _compare_genes_mixed(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
    duckdb_memory_limit: Optional[str] = None,
    duckdb_temp_directory: Optional[Union[str, Path]] = None,
    duckdb_threads: Optional[int] = None,
) -> pl.LazyFrame:
    """DuckDB shared-loci + Polars gene aggregation path."""
    lf = duckdb_filter_join(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        ani_method=ani_method,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
        memory_limit=duckdb_memory_limit,
        temp_directory=duckdb_temp_directory,
        threads=duckdb_threads,
    )
    return (
        lf.select(["genome", "gene", "surr"])
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


def compare_genomes(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    duckdb_memory_limit: Optional[str] = None,
    duckdb_temp_directory: Optional[Union[str, Path]] = None,
    duckdb_threads: Optional[int] = None,
    engine: Literal["polars", "duckdb"] = "polars",
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
) -> pl.LazyFrame:
    """Compare two profiles with selectable execution engine."""
    calculations = parse_genome_calculations(calculate)
    if engine == "polars":
        return compare_genomes_polars(
            mpile_contig_1=mpile_contig_1,
            mpile_contig_2=mpile_contig_2,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            ani_method=ani_method,
            stb_file=stb_file,
            calculate=calculations,
        )
    if engine == "duckdb":
        return duckdb_compare_genomes(
            mpile1=mpile_contig_1,
            mpile2=mpile_contig_2,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            ani_method=ani_method,
            calculate=calculations,
            stb_file=stb_file,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )
    raise ValueError(f"Unsupported engine: {engine}")


def compare_genes(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: str = "popani",
    duckdb_memory_limit: Optional[str] = None,
    duckdb_temp_directory: Optional[Union[str, Path]] = None,
    duckdb_threads: Optional[int] = None,
    engine: Literal["polars", "duckdb"] = "polars",
) -> pl.LazyFrame:
    """Compare two profiles at gene level with selectable execution engine."""
    if engine == "polars":
        return compare_genes_polars(
            mpile_contig_1=mpile_contig_1,
            mpile_contig_2=mpile_contig_2,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
            ani_method=ani_method,
        )
    if engine == "duckdb":
        return _compare_genes_mixed(
            mpile_contig_1=mpile_contig_1,
            mpile_contig_2=mpile_contig_2,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
            ani_method=ani_method,
            duckdb_memory_limit=duckdb_memory_limit,
            duckdb_temp_directory=duckdb_temp_directory,
            duckdb_threads=duckdb_threads,
        )
    raise ValueError(f"Unsupported engine: {engine}")
