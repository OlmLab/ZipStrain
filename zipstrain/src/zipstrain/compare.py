"""zipstrain.compare
========================
This module provides all comparison functions for zipstrain.

"""

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable, Literal, Optional, Union
import math
import re

import polars as pl
import duckdb


PROFILE_SORTED_METADATA_KEY = "zipstrain_sorted_by"
PROFILE_SORTED_METADATA_VALUE = "chrom,pos"

# Stride between scaffolds on the global comparison axis. Must exceed the length
# of any single scaffold so that scaffold_id * GENE_AXIS_OFFSET + pos stays
# monotone and no gene range can straddle a scaffold boundary.
GENE_AXIS_OFFSET = 1 << 32


GENOME_COMPARISON_CALCULATIONS = ("genome_ani", "ibs", "gene")
GENOME_COMPARISON_DEFAULT_CALCULATIONS = ("genome_ani", "ibs")
GENOME_COMPARISON_CALCULATION_ALIASES = {
    "genome_ani": "genome_ani",
    "ani": "genome_ani",
    "popani": "genome_ani",
    "ibs": "ibs",
    "max_block": "ibs",
    "max_consecutive_length": "ibs",
    "gene": "gene",
    "genes": "gene",
    "gene_ani": "gene",
}


def parse_ani_methods(
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> tuple[str, ...]:
    """Parse one or more ANI methods while preserving their requested order."""
    if ani_method is None:
        raw_methods = ["popani"]
    elif isinstance(ani_method, str):
        raw_methods = [part.strip().lower() for part in ani_method.split(",") if part.strip()]
    else:
        raw_methods = []
        for value in ani_method:
            raw_methods.extend(
                part.strip().lower() for part in str(value).split(",") if part.strip()
            )
    if not raw_methods:
        raise ValueError("At least one ANI method is required.")

    methods: list[str] = []
    for method in raw_methods:
        if method not in {"popani", "conani"}:
            if not method.startswith("cosani_"):
                raise ValueError(
                    f"Unsupported ANI method {method!r}. "
                    "Supported values: popani, conani, cosani_<threshold>."
                )
            threshold_text = method.split("_", 1)[1]
            try:
                threshold = float(threshold_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid cosANI threshold {threshold_text!r}; "
                    "expected a number between 0 and 1."
                ) from exc
            if not math.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
                raise ValueError("cosANI threshold must be between 0 and 1.")
        if method not in methods:
            methods.append(method)
    return tuple(methods)


def canonical_ani_methods(
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> str:
    """Return the canonical comma-separated ANI method string."""
    return ",".join(parse_ani_methods(ani_method))


def ani_method_suffix(ani_method: str) -> str:
    """Return a stable method suffix suitable for a parquet column name."""
    method = parse_ani_methods(ani_method)
    if len(method) != 1:
        raise ValueError("A column suffix requires exactly one ANI method.")
    return re.sub(r"[^a-z0-9]+", "_", method[0]).strip("_")


def ani_metric_column(base: str, ani_method: str, ani_methods: Iterable[str]) -> str:
    """Name one method-dependent metric, preserving the single-method schema."""
    methods = tuple(ani_methods)
    return base if len(methods) == 1 else f"{base}_{ani_method_suffix(ani_method)}"


def ani_match_column(ani_method: str, ani_methods: Iterable[str]) -> str:
    """Internal per-locus match-signal column for one ANI method."""
    return ani_metric_column("surr", ani_method, ani_methods)


def parse_genome_calculations(
    calculate: Optional[Union[str, Iterable[str]]] = None,
    include_gene_from_all: bool = False,
) -> tuple[str, ...]:
    """Parse and normalize comparison metric selection tokens.

    Accepted input formats:
      - None -> default ("genome_ani", "ibs")
      - "genome_ani+ibs+gene"
      - "all"
      - iterable of token strings

    ``gene`` needs a gene range table, so ``"all"`` only pulls it in when
    ``include_gene_from_all`` says one is available; asking for it explicitly
    without ranges is an error raised by the caller rather than a silent drop.
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

    normalized: set[str] = set()
    if "all" in raw_tokens:
        normalized.update(("genome_ani", "ibs"))
        if include_gene_from_all:
            normalized.add("gene")

    for token in raw_tokens:
        if token == "all":
            continue
        mapped = GENOME_COMPARISON_CALCULATION_ALIASES.get(token)
        if mapped is None:
            supported = "all," + ",".join(GENOME_COMPARISON_CALCULATIONS)
            raise ValueError(f"Unsupported calculation '{token}'. Supported values: {supported}")
        normalized.add(mapped)
    return tuple(metric for metric in GENOME_COMPARISON_CALCULATIONS if metric in normalized)


def comparison_is_gene_grained(calculate: Optional[Union[str, Iterable[str]]] = None) -> bool:
    """True when the output table is one row per gene rather than per genome."""
    return "gene" in parse_genome_calculations(calculate)


def genome_metric_output_columns(
    calculate: Optional[Union[str, Iterable[str]]] = None,
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> list[str]:
    """Return ordered output columns for a genome-grained comparison table."""
    calculations = parse_genome_calculations(calculate)
    ani_methods = parse_ani_methods(ani_method)
    cols = ["genome"]
    if "genome_ani" in calculations:
        cols.append("total_positions")
        for method in ani_methods:
            cols.extend(
                [
                    ani_metric_column("share_allele_pos", method, ani_methods),
                    ani_metric_column("genome_ani", method, ani_methods),
                ]
            )
    if "ibs" in calculations:
        cols.extend(
            ani_metric_column("max_consecutive_length", method, ani_methods)
            for method in ani_methods
        )
    return cols


def gene_metric_output_columns(
    calculate: Optional[Union[str, Iterable[str]]] = None,
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> list[str]:
    """Return ordered output columns for a gene-grained comparison table.

    ``total_positions``/``share_allele_pos`` are the gene's own counts; the
    genome-level metrics that were also requested ride along on every gene row
    so the table can be read without joining back to a genome table.
    """
    calculations = parse_genome_calculations(calculate)
    ani_methods = parse_ani_methods(ani_method)
    cols = ["genome", "gene", "total_positions"]
    for method in ani_methods:
        cols.extend(
            [
                ani_metric_column("share_allele_pos", method, ani_methods),
                ani_metric_column("gene_ani", method, ani_methods),
            ]
        )
    if "genome_ani" in calculations:
        cols.extend(
            ani_metric_column("genome_ani", method, ani_methods)
            for method in ani_methods
        )
    if "ibs" in calculations:
        cols.extend(
            ani_metric_column("max_consecutive_length", method, ani_methods)
            for method in ani_methods
        )
    return cols


def comparison_output_columns(
    calculate: Optional[Union[str, Iterable[str]]] = None,
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> list[str]:
    """Ordered output columns for whichever grain the selection implies."""
    if comparison_is_gene_grained(calculate):
        return gene_metric_output_columns(calculate, ani_method)
    return genome_metric_output_columns(calculate, ani_method)


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
    ani_method: Optional[Union[str, Iterable[str]]],
    gene_scope: str = "all",
) -> str:
    ani_methods = parse_ani_methods(ani_method)
    ani_selects = ",\n      ".join(
        f"{_duckdb_ani_expression(method)} AS {ani_match_column(method, ani_methods)}"
        for method in ani_methods
    )
    genome_scope_sql = _duckdb_quote_sql_string(genome_scope)
    _ = gene_scope  # gene scoping is applied to gene ranges, not to positions
    return f"""
    WITH p1f AS (
      SELECT
        chrom,
        pos,
        genome,
        A,
        T,
        C,
        G
      FROM {p1_source} p1
      WHERE (A + T + C + G) >= {min_cov}
        AND ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
    ),
    p2f AS (
      SELECT
        chrom,
        pos,
        genome,
        A,
        T,
        C,
        G
      FROM {p2_source} p2
      WHERE (A + T + C + G) >= {min_cov}
        AND ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
    )
    SELECT
      {ani_selects},
      p1.chrom AS scaffold,
      p1.pos,
      p1.genome
    FROM p1f p1
    INNER JOIN p2f p2
      ON p1.chrom = p2.chrom AND p1.pos = p2.pos
    """


def _duckdb_copy_query_to_parquet(conn: duckdb.DuckDBPyConnection, query: str, output_file: Union[str, Path]) -> None:
    output_sql = _duckdb_quote_sql_string(str(output_file))
    conn.execute(f"COPY ({query}) TO '{output_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)")


def _duckdb_build_query_with_ctes(ctes: list[str], final_select: str) -> str:
    return "WITH\n" + ",\n".join(ctes) + "\n" + final_select


def _join_stb_requested_genomes(
    genome_comp: pl.LazyFrame,
    stb_file: Union[str, Path],
    genome_scope: str,
) -> pl.LazyFrame:
    genomes_utf8 = (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(pl.col("column_2").cast(pl.Utf8).str.strip_chars().alias("genome"))
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
  SELECT DISTINCT trim(CAST(column1 AS VARCHAR)) AS genome
  FROM read_csv('{stb_sql}', delim='\\t', header=false)
  WHERE ('{genome_scope_sql}' = 'all' OR trim(CAST(column1 AS VARCHAR)) = '{genome_scope_sql}')
)""".strip()


def _duckdb_contig_pop_max_ctes(
    shared_source: str = "shared",
    include_pop: bool = True,
    include_max_blocks: bool = True,
    ani_method: str = "popani",
    ani_methods: Optional[Iterable[str]] = None,
) -> list[str]:
    methods = parse_ani_methods(ani_methods or (ani_method,))
    match_col = ani_match_column(ani_method, methods)
    tag = "" if len(methods) == 1 else f"_{ani_method_suffix(ani_method)}"
    contig_base = f"contig_base{tag}"
    contig = f"contig{tag}"
    pop = f"pop{tag}"
    max_blocks = f"max_blocks{tag}"
    ctes: list[str] = []
    if not include_pop and not include_max_blocks:
        return ctes
    ctes.extend(
        [
            f"""
{contig_base} AS (
  SELECT
    s.*,
    LAG(scaffold) OVER (PARTITION BY genome ORDER BY scaffold, pos) AS prev_scaffold
  FROM {shared_source} s
)""".strip(),
            f"""
{contig} AS (
  SELECT
    *,
    SUM(
      CASE
        WHEN prev_scaffold IS NULL OR scaffold != prev_scaffold OR {match_col} = 0 THEN 1
        ELSE 0
      END
    ) OVER (PARTITION BY genome ORDER BY scaffold, pos ROWS UNBOUNDED PRECEDING) AS group_id
  FROM {contig_base}
)""".strip(),
        ]
    )
    if include_pop:
        ctes.append(
            f"""
{pop} AS (
  SELECT
    genome,
    COUNT(*)::BIGINT AS total_positions,
    SUM(CASE WHEN {match_col} > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos,
    SUM(CASE WHEN {match_col} > 0 THEN 1 ELSE 0 END) * 100.0 / NULLIF(COUNT(*), 0) AS genome_ani
  FROM {contig}
  GROUP BY genome
)""".strip()
        )
    if include_max_blocks:
        ctes.append(
            f"""
{max_blocks} AS (
  SELECT
    genome,
    MAX(length)::BIGINT AS max_consecutive_length
  FROM (
    SELECT genome, scaffold, group_id, COUNT(*)::BIGINT AS length
    FROM {contig}
    GROUP BY genome, scaffold, group_id
  ) blocks
  GROUP BY genome
)""".strip(),
        )
    return ctes


def _attach_genome_metrics(
    gene_frame: pl.LazyFrame,
    genome_frame: pl.LazyFrame,
    calculations: tuple[str, ...],
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> pl.LazyFrame:
    """Repeat the requested genome-level metrics on every gene row.

    Only the genome-level metrics ride along; the count columns on a gene row are
    the gene's own, so the genome-level ones are dropped before the join.
    """
    ani_methods = parse_ani_methods(ani_method)
    carried: list[str] = []
    if "genome_ani" in calculations:
        carried.extend(
            ani_metric_column("genome_ani", method, ani_methods)
            for method in ani_methods
        )
    if "ibs" in calculations:
        carried.extend(
            ani_metric_column("max_consecutive_length", method, ani_methods)
            for method in ani_methods
        )
    gene_casts: list[pl.Expr] = [pl.col("total_positions").cast(pl.Int64)]
    for method in ani_methods:
        gene_casts.extend(
            [
                pl.col(ani_metric_column("share_allele_pos", method, ani_methods)).cast(pl.Int64),
                pl.col(ani_metric_column("gene_ani", method, ani_methods)).cast(pl.Float64),
            ]
        )
    gene_frame = gene_frame.with_columns(gene_casts)
    if carried:
        gene_frame = gene_frame.join(
            genome_frame.select(["genome", *carried]), on="genome", how="left"
        )
    return gene_frame.select(gene_metric_output_columns(calculations, ani_methods))


def _duckdb_gene_frame(
    con: duckdb.DuckDBPyConnection,
    shared_query: str,
    gene_range: Union[str, Path, pl.LazyFrame, pl.DataFrame],
    profile: Union[str, Path, pl.LazyFrame],
    min_gene_compare_len: int,
    tmp_dir: Path,
    stb_file: Optional[Union[str, Path]] = None,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> pl.LazyFrame:
    """Gene-level results for the DuckDB engine.

    DuckDB performs the shared-loci join and the sort, spilling as needed, and
    hands off a sorted parquet. The range sums are then a merge over sorted
    input, which is Polars' ``join_asof``; an ASOF in SQL needs a partition key
    to avoid degrading and is still several times slower even with one.
    """
    axis_map, gene_ranges = _gene_range_axis_inputs(
        gene_range,
        profile=profile,
        stb_file=stb_file,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
    )
    con.register(
        "scaffold_axis",
        pl.DataFrame(
            {"scaffold": list(axis_map.keys()), "scaffold_id": list(axis_map.values())},
            schema={"scaffold": pl.Utf8, "scaffold_id": pl.Int64},
        ).to_arrow(),
    )
    shared_axis_path = tmp_dir / "shared_axis.parquet"
    _duckdb_export_sorted_shared_axis(
        con, shared_query, shared_axis_path, ani_method=ani_method
    )
    return get_gene_ani(
        shared=pl.scan_parquet(shared_axis_path).set_sorted("axis"),
        gene_bounds=gene_boundary_table(gene_ranges, axis_map),
        min_gene_compare_len=min_gene_compare_len,
        ani_method=ani_method,
    )


def _duckdb_export_sorted_shared_axis(
    con: duckdb.DuckDBPyConnection,
    shared_query: str,
    output_path: Path,
    axis_map_source: str = "scaffold_axis",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> None:
    """Write shared positions to parquet on the global axis, sorted.

    DuckDB is used here for what it is genuinely best at: a large join plus an
    out-of-core sort that spills. The gene aggregation that follows is a
    two-pointer merge, which Polars does far better, so it is done there instead
    of with an ASOF join in SQL.
    """
    out_sql = _duckdb_quote_sql_string(str(output_path))
    ani_methods = parse_ani_methods(ani_method)
    match_selects = ",\n            ".join(
        f"s.{ani_match_column(method, ani_methods)}" for method in ani_methods
    )
    con.execute(
        f"""
        COPY (
          WITH shared AS (
{shared_query}
          )
          SELECT
            {match_selects},
            sa.scaffold_id * {GENE_AXIS_OFFSET} + s.pos AS axis
          FROM shared s
          JOIN {axis_map_source} sa ON s.scaffold = sa.scaffold
          ORDER BY axis
        ) TO '{out_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )

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
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> str:
    calculations = parse_genome_calculations(calculate)
    ani_methods = parse_ani_methods(ani_method)
    need_genome_ani = "genome_ani" in calculations
    need_ibs = "ibs" in calculations

    genome_scope_sql = _duckdb_quote_sql_string(genome_scope)
    ctes = [f"shared AS (\n{shared_query}\n)"]
    if stb_file is None:
        ctes.append(_duckdb_genomes_from_shared_cte(genome_scope_sql=genome_scope_sql))
    else:
        stb_sql = _duckdb_quote_sql_string(str(stb_file))
        ctes.append(_duckdb_genomes_scope_cte(stb_sql=stb_sql, genome_scope_sql=genome_scope_sql))

    if need_genome_ani or need_ibs:
        for method in ani_methods:
            ctes.extend(
                _duckdb_contig_pop_max_ctes(
                    shared_source="shared",
                    include_pop=need_genome_ani,
                    include_max_blocks=need_ibs,
                    ani_method=method,
                    ani_methods=ani_methods,
                )
            )
    # Genome-level metrics are selected the same way in both grains; when genes
    # are requested the driving table becomes gene_ani, so each gene row carries
    # its own counts plus the repeated genome-level values.
    genome_selects: list[str] = []
    joins: list[str] = []
    for index, method in enumerate(ani_methods):
        tag = "" if len(ani_methods) == 1 else f"_{ani_method_suffix(method)}"
        alias = f"a{index}"
        if need_genome_ani:
            if index == 0:
                genome_selects.append(
                    f"COALESCE({alias}.total_positions, 0)::BIGINT AS total_positions"
                )
            shared_col = ani_metric_column("share_allele_pos", method, ani_methods)
            genome_ani_col = ani_metric_column("genome_ani", method, ani_methods)
            genome_selects.extend(
                [
                    f"COALESCE({alias}.share_allele_pos, 0)::BIGINT AS {shared_col}",
                    f"COALESCE({alias}.genome_ani, 0.0)::DOUBLE AS {genome_ani_col}",
                ]
            )
            joins.append(f"LEFT JOIN pop{tag} {alias} ON g.genome = {alias}.genome")
    for index, method in enumerate(ani_methods):
        tag = "" if len(ani_methods) == 1 else f"_{ani_method_suffix(method)}"
        if need_ibs:
            ibs_alias = f"b{index}"
            ibs_col = ani_metric_column("max_consecutive_length", method, ani_methods)
            genome_selects.append(
                f"COALESCE({ibs_alias}.max_consecutive_length, 0)::BIGINT AS {ibs_col}"
            )
            joins.append(
                f"LEFT JOIN max_blocks{tag} {ibs_alias} ON g.genome = {ibs_alias}.genome"
            )

    sample_selects: list[str] = []
    if sample_1_name is not None and sample_2_name is not None:
        sample_1_sql = _duckdb_quote_sql_string(sample_1_name)
        sample_2_sql = _duckdb_quote_sql_string(sample_2_name)
        sample_selects = [f"'{sample_1_sql}' AS sample_1", f"'{sample_2_sql}' AS sample_2"]

    select_parts = ["g.genome AS genome", *genome_selects, *sample_selects]
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


GENE_RANGE_TABLE_COLUMNS = ("gene", "scaffold", "start", "end")


def _as_gene_range_frame(
    gene_range: Union[str, Path, pl.LazyFrame, pl.DataFrame],
) -> pl.LazyFrame:
    """Normalize a gene range table into a lazy ``gene, scaffold, start, end`` frame.

    On disk this is the headerless TSV emitted by ``zipstrain gene-range-table``;
    in-memory frames are accepted as-is so callers can build ranges themselves.
    """
    if isinstance(gene_range, pl.DataFrame):
        return gene_range.lazy()
    if isinstance(gene_range, pl.LazyFrame):
        return gene_range
    if isinstance(gene_range, (str, Path)):
        return pl.scan_csv(
            gene_range,
            has_header=False,
            separator="\t",
            new_columns=list(GENE_RANGE_TABLE_COLUMNS),
        ).select(
            pl.col("gene").cast(pl.Utf8),
            pl.col("scaffold").cast(pl.Utf8),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64),
        )
    raise TypeError(f"Unsupported source type for gene range table: {type(gene_range)}")


def _stb_scaffold_genome_map(stb_file: Union[str, Path]) -> pl.LazyFrame:
    """Scaffold-to-genome pairs from the headerless STB TSV (scaffold, genome)."""
    return (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(
            pl.col("column_1").cast(pl.Utf8).str.strip_chars().alias("scaffold"),
            pl.col("column_2").cast(pl.Utf8).str.strip_chars().alias("genome"),
        )
        .unique()
    )


def _scaffold_genome_map(
    source: Union[str, Path, pl.LazyFrame],
    stb_file: Optional[Union[str, Path]] = None,
) -> pl.LazyFrame:
    """Scaffold-to-genome pairs, used to attach a genome to each gene range.

    The gene range table is a headerless ``gene, scaffold, start, end`` TSV with
    no genome column. The STB file already holds that mapping and is tiny, so it
    is preferred; falling back to the profile costs a full pass over its columns.
    """
    if stb_file is not None:
        return _stb_scaffold_genome_map(stb_file)
    return (
        _as_lazy_profile(source)
        .select(pl.col("chrom").alias("scaffold"), pl.col("genome"))
        .unique()
    )


def _gene_range_axis_inputs(
    gene_range: Union[str, Path, pl.LazyFrame, pl.DataFrame],
    profile: Union[str, Path, pl.LazyFrame],
    stb_file: Optional[Union[str, Path]] = None,
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> tuple[dict[str, int], pl.LazyFrame]:
    """Build the axis map and the genome-labelled gene ranges for a comparison.

    The axis map is keyed only on scaffolds that actually carry genes. Positions
    on any other scaffold can be dropped before the cumulative pass without
    changing a single gene's counts: a gene's value is the difference of two
    boundaries on its own scaffold, so anything outside it shifts both ends
    equally and cancels. That keeps this off the profile entirely -- the gene
    range table and STB file are both small.
    """
    gene_ranges = _as_gene_range_frame(gene_range).join(
        _scaffold_genome_map(profile, stb_file=stb_file), on="scaffold", how="inner"
    )
    if genome_scope != "all":
        gene_ranges = gene_ranges.filter(pl.col("genome") == genome_scope)
    if gene_scope != "all":
        gene_ranges = gene_ranges.filter(pl.col("gene") == gene_scope)
    gene_ranges = gene_ranges.collect().lazy()
    axis_map = scaffold_axis_map(
        gene_ranges.select("scaffold").unique().collect()["scaffold"].to_list()
    )
    return axis_map, gene_ranges


def _profile_scope_predicate(genome_scope: str = "all", gene_scope: str = "all") -> pl.Expr:
    expr = pl.lit(True)
    if genome_scope != "all":
        expr = expr & (pl.col("genome") == genome_scope)
    if gene_scope != "all":
        expr = expr & (pl.col("gene") == gene_scope)
    return expr


def _profile_is_coordinate_sorted(source: Union[str, Path, pl.LazyFrame]) -> bool:
    if not isinstance(source, (str, Path)):
        return False
    try:
        metadata = pl.read_parquet_metadata(source)
    except Exception:
        return False
    return metadata.get(PROFILE_SORTED_METADATA_KEY) == PROFILE_SORTED_METADATA_VALUE


def _filter_profiles_polars(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,

) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    cov_expr = (pl.col("A") + pl.col("T") + pl.col("C") + pl.col("G")) >= min_cov
    p1 = _as_lazy_profile(mpile1).filter(cov_expr)
    p2 = _as_lazy_profile(mpile2).filter(cov_expr)
    return p1, p2


def _shared_loci_polars(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> pl.LazyFrame:
    ani_methods = parse_ani_methods(ani_method)
    ani_expressions = [
        getattr(PolarsANIExpressions(), method)().alias(
            ani_match_column(method, ani_methods)
        )
        for method in ani_methods
    ]
    p1, p2 = _filter_profiles_polars(
        mpile1=mpile1,
        mpile2=mpile2,
        min_cov=min_cov,
    )
    if _profile_is_coordinate_sorted(mpile1) and _profile_is_coordinate_sorted(mpile2):
        p1 = p1.set_sorted(["chrom", "pos"])
        p2 = p2.set_sorted(["chrom", "pos"])
    return (
        p1.join(
            p2,
            on=["chrom", "pos"],
            how="inner",
            suffix="_2",
        )
        .with_columns(ani_expressions)
        .select(
            *[pl.col(ani_match_column(method, ani_methods)) for method in ani_methods],
            scaffold=pl.col("chrom"),
            pos=pl.col("pos"),
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


def polars_prefilter_by_scope(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    genome_scope: str = "all",
    gene_scope: str = "all",
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Scope-filter both profiles in Polars and return in-memory LazyFrames."""
    scope_requested = genome_scope != "all" or gene_scope != "all"
    if not scope_requested:
        return _as_lazy_profile(mpile1), _as_lazy_profile(mpile2)

    def _filter_one(source: Union[str, Path, pl.LazyFrame]) -> pl.LazyFrame:
        lf = _as_lazy_profile(source)
        filters: list[pl.Expr] = []
        if genome_scope != "all":
            filters.append(pl.col("genome") == genome_scope)
        if gene_scope != "all":
            filters.append(pl.col("gene") == gene_scope)
        if filters:
            lf = lf.filter(pl.all_horizontal(*filters))
        return lf.collect(engine="streaming").lazy()

    return _filter_one(mpile1), _filter_one(mpile2)


def duckdb_filter_join(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int,
    genome_scope: str = "all",   # or specific genome
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
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


def duckdb_compare_genomes_to_parquet(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    output_file: Union[str, Path],
    stb_file: Optional[Union[str, Path]],
    sample_1_name: str,
    sample_2_name: str,
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
    calculate: Optional[Union[str, Iterable[str]]] = None,
    gene_range: Optional[Union[str, Path, pl.LazyFrame, pl.DataFrame]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> None:
    """Run comparison in DuckDB and write final output directly to parquet.

    This path avoids materializing large intermediate tables in Python memory.
    Emits a gene-grained table when ``gene`` is among the calculations.
    """
    calculations = parse_genome_calculations(
        calculate, include_gene_from_all=gene_range is not None
    )
    ani_methods = parse_ani_methods(ani_method)
    if "gene" in calculations and gene_range is None:
        raise ValueError(
            "Gene comparison was requested but no gene range table was provided. "
            "Pass gene_range=<gene range table> or drop 'gene' from --calculate."
        )
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
            # Pass the already-resolved tuple: re-parsing the raw token here
            # would drop 'gene' from "all", since availability of gene ranges is
            # only known at this level.
            calculate=calculations,
            stb_file=stb_file,
            sample_1_name=sample_1_name,
            sample_2_name=sample_2_name,
            ani_method=ani_methods,
        )
        if "gene" not in calculations:
            _duckdb_copy_query_to_parquet(con, query, output_file)
            return
        genome_frame = pl.from_arrow(con.execute(query).fetch_arrow_table()).lazy()
        with TemporaryDirectory(prefix="zipstrain_gene_compare_") as tmp:
            gene_frame = _duckdb_gene_frame(
                con,
                shared_query=shared_query,
                gene_range=gene_range,
                profile=mpile1,
                min_gene_compare_len=min_gene_compare_len,
                tmp_dir=Path(tmp),
                stb_file=stb_file,
                genome_scope=genome_scope,
                gene_scope=gene_scope,
                ani_method=ani_methods,
            )
            out = _attach_genome_metrics(
                gene_frame, genome_frame, calculations, ani_methods
            )
            if sample_1_name is not None and sample_2_name is not None:
                out = out.with_columns(
                    sample_1=pl.lit(sample_1_name), sample_2=pl.lit(sample_2_name)
                )
            out.sink_parquet(output_file, compression="zstd")
    finally:
        con.close()


def duckdb_compare_genomes(
    mpile1: Union[str, Path, pl.LazyFrame],
    mpile2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    gene_scope: str = "all",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
    calculate: Optional[Union[str, Iterable[str]]] = None,
    gene_range: Optional[Union[str, Path, pl.LazyFrame, pl.DataFrame]] = None,
    stb_file: Optional[Union[str, Path]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> pl.LazyFrame:
    """Run comparison in DuckDB and return selected metrics as a LazyFrame."""
    calculations = parse_genome_calculations(
        calculate, include_gene_from_all=gene_range is not None
    )
    ani_methods = parse_ani_methods(ani_method)
    if "gene" in calculations and gene_range is None:
        raise ValueError(
            "Gene comparison was requested but no gene range table was provided. "
            "Pass gene_range=<gene range table> or drop 'gene' from --calculate."
        )
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
            ani_method=ani_methods,
        )
        query = _duckdb_genome_compare_query(
            shared_query=shared_query,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            calculate=calculations,
            stb_file=stb_file,
            ani_method=ani_methods,
        )
        genome_frame = pl.from_arrow(con.execute(query).fetch_arrow_table()).lazy()
        if "gene" not in calculations:
            return genome_frame
        with TemporaryDirectory(prefix="zipstrain_gene_compare_") as tmp:
            gene_frame = _duckdb_gene_frame(
                con,
                shared_query=shared_query,
                gene_range=gene_range,
                profile=mpile1,
                min_gene_compare_len=min_gene_compare_len,
                tmp_dir=Path(tmp),
                stb_file=stb_file,
                genome_scope=genome_scope,
                gene_scope=gene_scope,
                ani_method=ani_methods,
            )
            return _attach_genome_metrics(
                gene_frame, genome_frame, calculations, ani_methods
            ).collect().lazy()
    finally:
        con.close()



def add_contiguity_info(
    mpile_contig: pl.LazyFrame,
    match_column: str = "surr",
    group_column: str = "group_id",
) -> pl.LazyFrame:
    """ Adds group id information to the lazy frame. If on the same scaffold and not popANI, then they are in the same group.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
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
    # or at a position that is itself an SNV (surr == 0). Marking those break
    # points and taking a cumulative sum numbers the blocks, so `group_id` is
    # constant within each block and increments at every break.
    break_expr = changed("scaffold") | (pl.col(match_column) == 0)
    if has_genome:
        break_expr = (
            changed("genome")
            | changed("scaffold")
            | (pl.col(match_column) == 0)
        )

    return mpile_contig.sort(sort_cols).with_columns(
        break_expr.cast(pl.Int64).cum_sum().alias(group_column)
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

def calculate_ani(
    mpile_contig: pl.LazyFrame,
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> pl.LazyFrame:
    """Aggregate one or more per-position ANI match signals by genome."""
    ani_methods = parse_ani_methods(ani_method)
    aggregations: list[pl.Expr] = [pl.len().alias("total_positions")]
    for method in ani_methods:
        shared_col = ani_metric_column("share_allele_pos", method, ani_methods)
        aggregations.append(
            (pl.col(ani_match_column(method, ani_methods)) > 0).sum().alias(shared_col)
        )
    result = mpile_contig.group_by("genome").agg(aggregations)
    return result.with_columns(
        (
            pl.col(ani_metric_column("share_allele_pos", method, ani_methods))
            * 100.0
            / pl.col("total_positions")
        ).alias(ani_metric_column("genome_ani", method, ani_methods))
        for method in ani_methods
    )


def calculate_pop_ani(mpile_contig:pl.LazyFrame) -> pl.LazyFrame:
    """
    Calculates the population ANI (Average Nucleotide Identity) for the given mpileup LazyFrame.
    NOTE: Remember that this function should be applied to the merged mpileup using `_shared_loci_polars`
    or the equivalent shared-loci helper for the active engine.

    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with population ANI information added.
    """
    return calculate_ani(mpile_contig, "popani")

def get_longest_consecutive_blocks(
    mpile_contig: pl.LazyFrame,
    group_column: str = "group_id",
    output_column: str = "max_consecutive_length",
) -> pl.LazyFrame:
    """
    Calculates the longest consecutive blocks for each genome in the mpileup LazyFrame for any genome.
    
    Args:
        mpile_contig (pl.LazyFrame): The input LazyFrame containing mpileup data.
    
    Returns:
        pl.LazyFrame: Updated LazyFrame with longest consecutive blocks information added.
    """
    block_lengths = (
        mpile_contig.group_by(["genome", "scaffold", group_column])
        .agg(pl.len().alias("length"))
    ) 
    return block_lengths.group_by("genome").agg(
        pl.col("length").max().alias(output_column)
    )

def scaffold_axis_map(scaffolds: Iterable[str]) -> dict[str, int]:
    """Map scaffold names to dense ids in lexicographic order.

    The ids must be assigned in the same order the shared-loci frame is sorted
    by (``chrom`` ascending), so that ``axis`` below is monotonically increasing
    across the whole frame. That single monotone axis is what lets gene ANI use
    a plain running total and a merge-style ``join_asof`` with no ``by=`` groups.

    Ordering is done by Polars rather than Python's ``sorted`` so the ranking
    matches the collation the frame was actually sorted with; the two can differ
    on non-ASCII scaffold names, which would silently break axis monotonicity.
    """
    names = pl.Series("scaffold", list(scaffolds), dtype=pl.Utf8).unique().sort().to_list()
    return {name: idx for idx, name in enumerate(names)}


def with_axis(
    frame: pl.LazyFrame,
    axis_map: dict[str, int],
    *,
    scaffold_col: str = "scaffold",
    pos_col: str = "pos",
    alias: str = "axis",
) -> pl.LazyFrame:
    """Fold (scaffold, pos) into one monotone int64 coordinate.

    ``axis = scaffold_id * GENE_AXIS_OFFSET + pos``. Positions are 1-based and
    always below the offset, so ordering by ``axis`` equals ordering by
    (scaffold, pos), and no genomic range ever straddles two scaffolds.
    """
    return frame.with_columns(
        (
            pl.col(scaffold_col).replace_strict(axis_map, default=None, return_dtype=pl.Int64)
            * GENE_AXIS_OFFSET
            + pl.col(pos_col).cast(pl.Int64)
        ).alias(alias)
    )


def gene_boundary_table(gene_range: pl.LazyFrame, axis_map: dict[str, int]) -> pl.LazyFrame:
    """Turn each gene range into two signed boundary rows on the global axis.

    A gene's counts are a range sum, and a range sum is the difference of two
    cumulative values: ``count(start..end) = count(<= end) - count(< start)``.
    Emitting a ``+1`` row at ``end`` and a ``-1`` row at ``start - 1`` turns that
    subtraction into a single signed sum per gene, so overlapping and nested
    genes simply read the same cumulative array at their own boundaries instead
    of competing for ownership of a position.
    """
    ranges = gene_range.filter(
        pl.col("scaffold").is_in(list(axis_map.keys())) & (pl.col("start") <= pl.col("end"))
    )
    low = ranges.select(
        "gene",
        "genome",
        "scaffold",
        pl.col("start").cast(pl.Int64) - 1,
        sign=pl.lit(-1, pl.Int64),
    ).rename({"start": "pos"})
    high = ranges.select(
        "gene",
        "genome",
        "scaffold",
        pl.col("end").cast(pl.Int64).alias("pos"),
        sign=pl.lit(1, pl.Int64),
    )
    return (
        with_axis(pl.concat([low, high]), axis_map, alias="baxis")
        .drop(["scaffold", "pos"])
        .sort("baxis")
    )


def get_gene_ani(
    shared: pl.LazyFrame,
    gene_bounds: pl.LazyFrame,
    min_gene_compare_len: int,
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
) -> pl.LazyFrame:
    """Per-gene ANI over shared positions, correct for overlapping genes.

    ``shared`` must carry ``axis`` and be sorted by it. One running total over
    the positions gives, at every row, how many shared positions (and how many
    matching ones) lie at or before it; each gene then needs only its two
    boundary values. Cost is one pass over the positions plus two lookups per
    gene, and the merge join keeps it streaming.

    Boundaries need not exist as rows: ``strategy="backward"`` resolves each to
    the last covered position at or before it, which is exactly the cumulative
    count wanted. A ``start - 1`` that precedes the scaffold's first position
    lands on the previous scaffold's last row, and the shared global offset
    means that value still cancels correctly in the subtraction.
    """
    # The running total is only meaningful in axis order, and the inner join that
    # produced `shared` does not guarantee it preserves input ordering, so sort
    # explicitly rather than assuming. Sorting one int64 key is cheap next to the
    # join itself, and getting this wrong is silent: an unsorted prefix yields
    # plausible but wrong counts.
    ani_methods = parse_ani_methods(ani_method)
    cumulative_shared = [
        (pl.col(ani_match_column(method, ani_methods)) > 0)
        .cast(pl.Int64)
        .cum_sum()
        .alias(f"cum_shared_{ani_method_suffix(method)}")
        for method in ani_methods
    ]
    prefix = (
        shared.sort("axis")
        .with_columns(
            cum_total=pl.int_range(1, pl.len() + 1, dtype=pl.Int64),
            *cumulative_shared,
        )
        .select(
            "axis",
            "cum_total",
            *[f"cum_shared_{ani_method_suffix(method)}" for method in ani_methods],
        )
        .set_sorted("axis")
    )
    aggregations: list[pl.Expr] = [
        (pl.col("cum_total").fill_null(0) * pl.col("sign"))
        .sum()
        .alias("total_positions")
    ]
    for method in ani_methods:
        aggregations.append(
            (
                pl.col(f"cum_shared_{ani_method_suffix(method)}").fill_null(0)
                * pl.col("sign")
            )
            .sum()
            .alias(ani_metric_column("share_allele_pos", method, ani_methods))
        )
    result = (
        gene_bounds.join_asof(
            prefix,
            left_on="baxis",
            right_on="axis",
            strategy="backward",
        )
        .group_by(["genome", "gene"])
        .agg(aggregations)
        .filter(pl.col("total_positions") >= min_gene_compare_len)
    )
    return result.with_columns(
        (
            pl.col(ani_metric_column("share_allele_pos", method, ani_methods))
            * 100.0
            / pl.col("total_positions")
        ).alias(ani_metric_column("gene_ani", method, ani_methods))
        for method in ani_methods
    )

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
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
    gene_range: Optional[Union[str, Path, pl.LazyFrame, pl.DataFrame]] = None,
    gene_scope: str = "all",
) -> pl.LazyFrame:
    """Compare two profiles fully in Polars.

    Returns a genome-grained table, or a gene-grained one when ``gene`` is among
    the calculations; in that case the requested genome-level metrics are
    repeated on every gene row.
    """
    calculations = parse_genome_calculations(
        calculate, include_gene_from_all=gene_range is not None
    )
    ani_methods = parse_ani_methods(ani_method)
    if "gene" in calculations and gene_range is None:
        raise ValueError(
            "Gene comparison was requested but no gene range table was provided. "
            "Pass gene_range=<gene range table> or drop 'gene' from --calculate."
        )
    shared = _shared_loci_polars(
        mpile1=mpile_contig_1,
        mpile2=mpile_contig_2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        ani_method=ani_methods,
    )
    genome_comp_parts: list[pl.LazyFrame] = []
    if "genome_ani" in calculations:
        genome_comp_parts.append(calculate_ani(shared, ani_methods))
    if "ibs" in calculations:
        for method in ani_methods:
            match_col = ani_match_column(method, ani_methods)
            group_col = f"group_id_{ani_method_suffix(method)}"
            output_col = ani_metric_column(
                "max_consecutive_length", method, ani_methods
            )
            genome_comp_parts.append(
                get_longest_consecutive_blocks(
                    add_contiguity_info(
                        shared,
                        match_column=match_col,
                        group_column=group_col,
                    ),
                    group_column=group_col,
                    output_column=output_col,
                )
            )

    if genome_comp_parts:
        genome_comp = genome_comp_parts[0]
        for part in genome_comp_parts[1:]:
            genome_comp = genome_comp.join(part, on="genome", how="left")
    else:
        genome_comp = shared.select("genome").unique()

    if stb_file is not None:
        genome_comp = _join_stb_requested_genomes(
            genome_comp=genome_comp,
            stb_file=stb_file,
            genome_scope=genome_scope,
        )

    genome_casts: list[pl.Expr] = []
    if "genome_ani" in calculations:
        genome_casts.append(pl.col("total_positions").fill_null(0).cast(pl.Int64))
        for method in ani_methods:
            genome_casts.extend(
                [
                    pl.col(
                        ani_metric_column("share_allele_pos", method, ani_methods)
                    ).fill_null(0).cast(pl.Int64),
                    pl.col(
                        ani_metric_column("genome_ani", method, ani_methods)
                    ).fill_null(0.0).cast(pl.Float64),
                ]
            )
    if "ibs" in calculations:
        genome_casts.extend(
            pl.col(
                ani_metric_column("max_consecutive_length", method, ani_methods)
            ).fill_null(0).cast(pl.Int64)
            for method in ani_methods
        )
    if genome_casts:
        genome_comp = genome_comp.with_columns(genome_casts)

    if "gene" not in calculations:
        return genome_comp.select(
            genome_metric_output_columns(calculations, ani_methods)
        )

    axis_map, gene_ranges = _gene_range_axis_inputs(
        gene_range,
        profile=mpile_contig_1,
        stb_file=stb_file,
        genome_scope=genome_scope,
        gene_scope=gene_scope,
    )
    gene_comp = get_gene_ani(
        # Positions on scaffolds without genes carry no axis and are dropped;
        # they cancel out of every gene's boundary difference anyway.
        shared=with_axis(shared, axis_map).drop_nulls("axis"),
        gene_bounds=gene_boundary_table(gene_ranges, axis_map),
        min_gene_compare_len=min_gene_compare_len,
        ani_method=ani_methods,
    )

    return _attach_genome_metrics(
        gene_comp, genome_comp, calculations, ani_methods
    )






def compare_genomes(
    mpile_contig_1: Union[str, Path, pl.LazyFrame],
    mpile_contig_2: Union[str, Path, pl.LazyFrame],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: Optional[Union[str, Iterable[str]]] = "popani",
    duckdb_memory_limit: Optional[str] = None,
    duckdb_temp_directory: Optional[Union[str, Path]] = None,
    duckdb_threads: Optional[int] = None,
    engine: Literal["polars", "duckdb"] = "polars",
    stb_file: Optional[Union[str, Path]] = None,
    calculate: Optional[Union[str, Iterable[str]]] = None,
    gene_range: Optional[Union[str, Path, pl.LazyFrame, pl.DataFrame]] = None,
    gene_scope: str = "all",
) -> pl.LazyFrame:
    """Compare two profiles with selectable execution engine.

    Emits a gene-grained table when ``gene`` is among the calculations, which
    requires ``gene_range``; otherwise the table is genome-grained.
    """
    calculations = parse_genome_calculations(
        calculate, include_gene_from_all=gene_range is not None
    )
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
            gene_range=gene_range,
            gene_scope=gene_scope,
        )
    if engine == "duckdb":
        return duckdb_compare_genomes(
            mpile1=mpile_contig_1,
            mpile2=mpile_contig_2,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            gene_scope=gene_scope,
            ani_method=ani_method,
            calculate=calculations,
            gene_range=gene_range,
            stb_file=stb_file,
            memory_limit=duckdb_memory_limit,
            temp_directory=duckdb_temp_directory,
            threads=duckdb_threads,
        )
    raise ValueError(f"Unsupported engine: {engine}")


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
