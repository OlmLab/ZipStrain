"""Lightweight profile/compare pipeline backed by DuckDB tables.

This module is intentionally separate from the existing zipstrain CLI and keeps the
legacy behavior untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, Union

import duckdb
import polars as pl

from zipstrain import compare as cp
from zipstrain import profile as pf


CALCULATION_ORDER = ("ani", "ibs", "genes")
CALCULATION_ALIASES = {
    "ani": "ani",
    "ibs": "ibs",
    "genes": "genes",
    "gene": "genes",
    "identical_genes": "genes",
    "all": "all",
}
CALCULATION_COLUMNS = {
    "ani": ("total_positions", "share_allele_pos", "genome_pop_ani"),
    "ibs": ("max_consecutive_length",),
    "genes": ("shared_genes_count", "identical_gene_count", "perc_id_genes"),
}
PROFILE_DUCKDB_FILENAME = "profile.duckdb"


@dataclass(frozen=True)
class LightProfileSummary:
    """Summary of a generated light profile database."""

    coverage_rows: int
    gene_stats_rows: int = 0
    genome_stats_rows: int = 0


@dataclass(frozen=True)
class _LightProfileInput:
    """Resolved light profile input source."""

    root: Path
    duckdb_path: Path


def _resolve_output_duckdb_path(
    output_file: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Path:
    """Resolve output location to a DuckDB file path.

    `output_dir` is kept for backwards compatibility with the earlier
    directory-based light profile layout.
    """
    if output_file is not None and output_dir is not None:
        raise ValueError("Provide only one of output_file or output_dir")
    if output_file is None and output_dir is None:
        raise ValueError("One of output_file or output_dir is required")

    if output_file is not None:
        path = Path(output_file)
    else:
        path = Path(output_dir) / PROFILE_DUCKDB_FILENAME

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _parse_calculate(calculate: str) -> tuple[str, ...]:
    """Parse modular calculation selection string."""
    raw = (calculate or "").strip().lower()
    if not raw:
        raise ValueError("calculate cannot be empty")
    tokens = [tok.strip() for tok in raw.split("+") if tok.strip()]
    if not tokens:
        raise ValueError("calculate must contain at least one calculation group")

    resolved: set[str] = set()
    for tok in tokens:
        if tok not in CALCULATION_ALIASES:
            valid = ", ".join(sorted(CALCULATION_ALIASES))
            raise ValueError(f"Unknown calculation '{tok}'. Valid values: {valid}")
        mapped = CALCULATION_ALIASES[tok]
        if mapped == "all":
            resolved = set(CALCULATION_ORDER)
            break
        resolved.add(mapped)

    if not resolved:
        raise ValueError("calculate must resolve to at least one group")
    return tuple(calc for calc in CALCULATION_ORDER if calc in resolved)


def _selected_output_columns(calculations: tuple[str, ...]) -> list[str]:
    cols = ["genome"]
    for calc in calculations:
        cols.extend(CALCULATION_COLUMNS[calc])
    return cols


def _configure_duckdb_connection(
    conn: duckdb.DuckDBPyConnection,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> None:
    """Apply optional DuckDB runtime settings."""
    if memory_limit:
        conn.execute(f"SET memory_limit = '{memory_limit}'")
    if temp_directory:
        conn.execute(f"SET temp_directory = '{Path(temp_directory)}'")
    if threads is not None:
        conn.execute(f"SET threads = {int(threads)}")


def _attached_has_table(conn: duckdb.DuckDBPyConnection, schema_name: str, table_name: str) -> bool:
    row = conn.execute(
        """
        SELECT COUNT(*)
        FROM information_schema.tables
        WHERE table_catalog = ? AND table_name = ?
        """,
        [schema_name, table_name],
    ).fetchone()
    return bool(row and row[0] > 0)


def _companion_stats_path(profile_parquet: Union[str, Path], stat_suffix: str) -> Path:
    profile_path = Path(profile_parquet)
    stem = profile_path.stem
    prefix = stem[: -len("_profile")] if stem.endswith("_profile") else stem
    return profile_path.with_name(f"{prefix}_{stat_suffix}.parquet")


def _resolve_requested_stats_paths(
    profile_parquet: Union[str, Path],
    include_gene_stats: bool,
    include_genome_stats: bool,
    gene_stats_parquet: Optional[Union[str, Path]] = None,
    genome_stats_parquet: Optional[Union[str, Path]] = None,
) -> tuple[Optional[Path], Optional[Path]]:
    gene_stats_path: Optional[Path] = None
    genome_stats_path: Optional[Path] = None

    if include_gene_stats:
        gene_stats_path = Path(gene_stats_parquet) if gene_stats_parquet is not None else _companion_stats_path(profile_parquet, "gene_stats")
        if not gene_stats_path.exists():
            if gene_stats_parquet is not None:
                raise FileNotFoundError(f"Requested gene_stats table but file was not found: {gene_stats_path}")
            gene_stats_path = None

    if include_genome_stats:
        genome_stats_path = Path(genome_stats_parquet) if genome_stats_parquet is not None else _companion_stats_path(profile_parquet, "genome_stats")
        if not genome_stats_path.exists():
            if genome_stats_parquet is not None:
                raise FileNotFoundError(f"Requested genome_stats table but file was not found: {genome_stats_path}")
            genome_stats_path = None

    return gene_stats_path, genome_stats_path


def _write_light_duckdb_profile(
    profile_parquet: Union[str, Path],
    min_cov: int,
    output_db: Union[str, Path],
    include_gene_stats: bool = True,
    include_genome_stats: bool = True,
    gene_stats_parquet: Optional[Union[str, Path]] = None,
    genome_stats_parquet: Optional[Union[str, Path]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Write light profile tables to a DuckDB database file from profile parquet."""
    output_path = Path(output_db)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    gene_stats_path, genome_stats_path = _resolve_requested_stats_paths(
        profile_parquet=profile_parquet,
        include_gene_stats=include_gene_stats,
        include_genome_stats=include_genome_stats,
        gene_stats_parquet=gene_stats_parquet,
        genome_stats_parquet=genome_stats_parquet,
    )

    profile_sql = cp._duckdb_quote_sql_string(str(profile_parquet))
    conn = duckdb.connect(str(output_path))
    try:
        _configure_duckdb_connection(conn, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        conn.execute(
            f"""
            CREATE TABLE mpileup AS
            SELECT
              CAST(p.chrom AS VARCHAR) AS chrom,
              CAST(p.pos AS BIGINT) AS pos,
              COALESCE(CAST(p.gene AS VARCHAR), 'NA') AS gene,
              COALESCE(CAST(p.genome AS VARCHAR), 'NA') AS genome,
              (
                CASE WHEN CAST(p.A AS INTEGER) > 0 THEN 1 ELSE 0 END
                + CASE WHEN CAST(p.T AS INTEGER) > 0 THEN 2 ELSE 0 END
                + CASE WHEN CAST(p.C AS INTEGER) > 0 THEN 4 ELSE 0 END
                + CASE WHEN CAST(p.G AS INTEGER) > 0 THEN 8 ELSE 0 END
              )::UTINYINT AS bit_representation
            FROM read_parquet('{profile_sql}') AS p
            WHERE (
              CAST(p.A AS INTEGER)
              + CAST(p.T AS INTEGER)
              + CAST(p.C AS INTEGER)
              + CAST(p.G AS INTEGER)
            ) > {int(min_cov)}
            ORDER BY genome, chrom, pos
            """
        )
        coverage_rows = int(conn.execute("SELECT COUNT(*) FROM mpileup").fetchone()[0])

        gene_stats_rows = 0
        if gene_stats_path is not None:
            gene_stats_sql = cp._duckdb_quote_sql_string(str(gene_stats_path))
            conn.execute(f"CREATE TABLE gene_stats AS SELECT * FROM read_parquet('{gene_stats_sql}')")
            gene_stats_rows = int(conn.execute("SELECT COUNT(*) FROM gene_stats").fetchone()[0])

        genome_stats_rows = 0
        if genome_stats_path is not None:
            genome_stats_sql = cp._duckdb_quote_sql_string(str(genome_stats_path))
            conn.execute(f"CREATE TABLE genome_stats AS SELECT * FROM read_parquet('{genome_stats_sql}')")
            genome_stats_rows = int(conn.execute("SELECT COUNT(*) FROM genome_stats").fetchone()[0])
    finally:
        conn.close()

    return LightProfileSummary(
        coverage_rows=coverage_rows,
        gene_stats_rows=gene_stats_rows,
        genome_stats_rows=genome_stats_rows,
    )


def _resolve_light_profile_input(profile_input: Union[str, Path]) -> _LightProfileInput:
    """Resolve a light profile input path to its DuckDB source."""
    path = Path(profile_input)
    if not path.exists():
        raise FileNotFoundError(f"Light profile path does not exist: {path}")
    if path.is_file():
        return _LightProfileInput(root=path.parent, duckdb_path=path)
    if not path.is_dir():
        raise ValueError(f"Unsupported light profile path: {path}")

    db_file = path / PROFILE_DUCKDB_FILENAME
    duckdb_files = sorted(path.glob("*.duckdb"))
    if not db_file.exists():
        if duckdb_files:
            extras = ", ".join(sorted(p.name for p in duckdb_files))
            raise ValueError(
                f"Invalid light profile directory {path}: expected {PROFILE_DUCKDB_FILENAME}, "
                f"found other db files ({extras})."
            )
        raise ValueError(f"Unrecognized light profile directory {path}: missing {PROFILE_DUCKDB_FILENAME}.")

    extra_db_files = [p for p in duckdb_files if p.name != PROFILE_DUCKDB_FILENAME]
    if extra_db_files:
        extras = ", ".join(sorted(p.name for p in extra_db_files))
        raise ValueError(
            f"Ambiguous light profile directory {path}: "
            f"only {PROFILE_DUCKDB_FILENAME} is allowed, found extra db files ({extras})."
        )

    return _LightProfileInput(root=path, duckdb_path=db_file)


def _register_light_profile_source(
    conn: duckdb.DuckDBPyConnection,
    source: _LightProfileInput,
    schema_name: str,
) -> None:
    """Attach a light profile under `schema_name`."""
    db_sql = cp._duckdb_quote_sql_string(str(source.duckdb_path))
    conn.execute(f"ATTACH '{db_sql}' AS {schema_name} (READ_ONLY)")
    if not _attached_has_table(conn, schema_name, "mpileup"):
        raise ValueError(f"Light profile database is missing required table {schema_name}.mpileup")


def build_light_profile_bundle_from_parquet(
    profile_parquet: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    min_cov: int = 5,
    include_gene_stats: bool = True,
    include_genome_stats: bool = True,
    gene_stats_parquet: Optional[Union[str, Path]] = None,
    genome_stats_parquet: Optional[Union[str, Path]] = None,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Create a light profile DuckDB file from a regular profile parquet."""
    output_db = _resolve_output_duckdb_path(output_file=output_file, output_dir=output_dir)
    if output_db.exists():
        output_db.unlink()

    return _write_light_duckdb_profile(
        profile_parquet=profile_parquet,
        min_cov=min_cov,
        output_db=output_db,
        include_gene_stats=include_gene_stats,
        include_genome_stats=include_genome_stats,
        gene_stats_parquet=gene_stats_parquet,
        genome_stats_parquet=genome_stats_parquet,
        memory_limit=memory_limit,
        temp_directory=temp_directory,
        threads=threads,
    )


def build_light_profile_bundle_from_bam(
    bam_file: Union[str, Path],
    bed_file: Union[str, Path],
    gene_range_table: Union[str, Path],
    stb_file: Union[str, Path],
    null_model: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    min_cov: int = 5,
    include_gene_stats: bool = True,
    include_genome_stats: bool = True,
    num_chunks: int = 24,
    max_concurrency: int = 4,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Profile a BAM and write a light profile DuckDB file."""
    with TemporaryDirectory(prefix="zipstrain_light_profile_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        stb = pl.scan_csv(stb_file, separator="\t", has_header=False).with_columns(
            pl.col("column_1").alias("scaffold"),
            pl.col("column_2").alias("genome"),
        )
        null_model_lf = pl.scan_parquet(null_model)
        pf.profile_bam(
            bed_file=str(bed_file),
            bam_file=str(bam_file),
            gene_range_table=str(gene_range_table),
            stb=stb,
            null_model=null_model_lf,
            output_dir=str(tmp_path),
            num_chunks=num_chunks,
            max_concurrency=max_concurrency,
        )
        sample_name = Path(bam_file).stem
        profile_parquet = tmp_path / f"{sample_name}_profile.parquet"
        if not profile_parquet.exists():
            raise FileNotFoundError(f"Expected profile output not found: {profile_parquet}")

        gene_stats_parquet = tmp_path / f"{sample_name}_gene_stats.parquet"
        genome_stats_parquet = tmp_path / f"{sample_name}_genome_stats.parquet"
        return build_light_profile_bundle_from_parquet(
            profile_parquet=profile_parquet,
            output_file=output_file,
            output_dir=output_dir,
            min_cov=min_cov,
            include_gene_stats=include_gene_stats,
            include_genome_stats=include_genome_stats,
            gene_stats_parquet=gene_stats_parquet if include_gene_stats else None,
            genome_stats_parquet=genome_stats_parquet if include_genome_stats else None,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )


def _light_genomes_cte(genome_scope_sql: str) -> str:
    return f"""
genomes AS (
  SELECT DISTINCT genome
  FROM (
    SELECT CAST(genome AS VARCHAR) AS genome FROM p1.mpileup
    UNION
    SELECT CAST(genome AS VARCHAR) AS genome FROM p2.mpileup
  ) g
  WHERE ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
)""".strip()


def _light_shared_cte(genome_scope_sql: str) -> str:
    return f"""
shared AS (
  SELECT
    CAST(c1.genome AS VARCHAR) AS genome,
    CAST(c1.chrom AS VARCHAR) AS scaffold,
    CAST(c1.pos AS BIGINT) AS pos,
    COALESCE(CAST(c1.gene AS VARCHAR), 'NA') AS gene,
    CASE
      WHEN (CAST(c1.bit_representation AS INTEGER) & CAST(c2.bit_representation AS INTEGER)) > 0 THEN 1
      ELSE 0
    END::TINYINT AS surr
  FROM p1.mpileup c1
  INNER JOIN p2.mpileup c2
    ON CAST(c1.genome AS VARCHAR) = CAST(c2.genome AS VARCHAR)
   AND CAST(c1.chrom AS VARCHAR) = CAST(c2.chrom AS VARCHAR)
   AND CAST(c1.pos AS BIGINT) = CAST(c2.pos AS BIGINT)
  WHERE ('{genome_scope_sql}' = 'all' OR CAST(c1.genome AS VARCHAR) = '{genome_scope_sql}')
)""".strip()


def _light_pop_cte(source: str = "shared") -> str:
    return f"""
pop AS (
  SELECT
    genome,
    COUNT(*)::BIGINT AS total_positions,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END) * 100.0 / NULLIF(COUNT(*), 0) AS genome_pop_ani
  FROM {source}
  GROUP BY genome
)""".strip()


def _light_contig_ctes(shared_source: str = "shared") -> list[str]:
    return [
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


def _light_max_blocks_cte(contig_source: str = "contig") -> str:
    return f"""
max_blocks AS (
  SELECT
    genome,
    MAX(length)::BIGINT AS max_consecutive_length
  FROM (
    SELECT genome, scaffold, group_id, COUNT(*)::BIGINT AS length
    FROM {contig_source}
    GROUP BY genome, scaffold, group_id
  ) blocks
  GROUP BY genome
)""".strip()


def _light_gene_ctes(min_gene_compare_len: int, shared_source: str = "shared") -> list[str]:
    return [
        f"""
gene_base AS (
  SELECT
    genome,
    gene,
    COUNT(*)::BIGINT AS total_positions,
    SUM(CASE WHEN surr > 0 THEN 1 ELSE 0 END)::BIGINT AS share_allele_pos
  FROM {shared_source}
  WHERE gene != 'NA'
  GROUP BY genome, gene
  HAVING COUNT(*) >= {int(min_gene_compare_len)}
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


def _light_final_select_sql(calculations: tuple[str, ...]) -> str:
    select_cols = ["g.genome AS genome"]
    joins = []

    if "ani" in calculations:
        select_cols.extend(
            [
                "COALESCE(p.total_positions, 0)::BIGINT AS total_positions",
                "COALESCE(p.share_allele_pos, 0)::BIGINT AS share_allele_pos",
                "COALESCE(p.genome_pop_ani, 0.0)::DOUBLE AS genome_pop_ani",
            ]
        )
        joins.append("LEFT JOIN pop p ON g.genome = p.genome")

    if "ibs" in calculations:
        select_cols.append("COALESCE(m.max_consecutive_length, 0)::BIGINT AS max_consecutive_length")
        joins.append("LEFT JOIN max_blocks m ON g.genome = m.genome")

    if "genes" in calculations:
        select_cols.extend(
            [
                "COALESCE(gs.shared_genes_count, 0)::BIGINT AS shared_genes_count",
                "COALESCE(gs.identical_gene_count, 0)::BIGINT AS identical_gene_count",
                "COALESCE(gs.perc_id_genes, 0.0)::DOUBLE AS perc_id_genes",
            ]
        )
        joins.append("LEFT JOIN gene_stats gs ON g.genome = gs.genome")

    return (
        "SELECT\n  "
        + ",\n  ".join(select_cols)
        + "\nFROM genomes g\n"
        + ("\n".join(joins) + "\n" if joins else "")
        + "ORDER BY g.genome"
    )


def _light_duckdb_compare_query(
    genome_scope: str,
    min_gene_compare_len: int,
    calculations: tuple[str, ...],
) -> str:
    genome_scope_sql = cp._duckdb_quote_sql_string(genome_scope)
    ctes = [
        _light_genomes_cte(genome_scope_sql=genome_scope_sql),
        _light_shared_cte(genome_scope_sql=genome_scope_sql),
    ]

    if "ani" in calculations:
        ctes.append(_light_pop_cte(source="shared"))

    if "ibs" in calculations:
        ctes.extend(_light_contig_ctes(shared_source="shared"))
        ctes.append(_light_max_blocks_cte(contig_source="contig"))

    if "genes" in calculations:
        ctes.extend(_light_gene_ctes(min_gene_compare_len=min_gene_compare_len, shared_source="shared"))

    return cp._duckdb_build_query_with_ctes(ctes, _light_final_select_sql(calculations=calculations))


def compare_light_genomes(
    profile_db_1: Union[str, Path],
    profile_db_2: Union[str, Path],
    stb_file: Optional[Union[str, Path]] = None,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
    calculate: str = "all",
) -> pl.LazyFrame:
    """Compare two light profiles and return genome-level metrics."""
    if ani_method != "popani":
        raise ValueError("zipstrain-light currently supports only --ani-method popani")

    calculations = _parse_calculate(calculate)
    con = duckdb.connect()
    try:
        _configure_duckdb_connection(con, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        source_1 = _resolve_light_profile_input(profile_db_1)
        source_2 = _resolve_light_profile_input(profile_db_2)
        _register_light_profile_source(con, source_1, "p1")
        _register_light_profile_source(con, source_2, "p2")

        query = _light_duckdb_compare_query(
            genome_scope=genome_scope,
            min_gene_compare_len=min_gene_compare_len,
            calculations=calculations,
        )
        table = con.execute(query).fetch_arrow_table()
        output_columns = _selected_output_columns(calculations)
        return pl.from_arrow(table).lazy().select(output_columns)
    finally:
        con.close()


def compare_light_genomes_to_parquet(
    profile_db_1: Union[str, Path],
    profile_db_2: Union[str, Path],
    output_file: Union[str, Path],
    sample_1_name: str,
    sample_2_name: str,
    stb_file: Optional[Union[str, Path]] = None,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
    calculate: str = "all",
) -> None:
    """Compare two light profiles and write final parquet output."""
    if ani_method != "popani":
        raise ValueError("zipstrain-light currently supports only --ani-method popani")

    calculations = _parse_calculate(calculate)
    con = duckdb.connect()
    try:
        _configure_duckdb_connection(con, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        source_1 = _resolve_light_profile_input(profile_db_1)
        source_2 = _resolve_light_profile_input(profile_db_2)
        _register_light_profile_source(con, source_1, "p1")
        _register_light_profile_source(con, source_2, "p2")

        query = _light_duckdb_compare_query(
            genome_scope=genome_scope,
            min_gene_compare_len=min_gene_compare_len,
            calculations=calculations,
        )
        sample_1_sql = cp._duckdb_quote_sql_string(sample_1_name)
        sample_2_sql = cp._duckdb_quote_sql_string(sample_2_name)
        out_query = (
            "SELECT *, "
            f"'{sample_1_sql}' AS sample_1, "
            f"'{sample_2_sql}' AS sample_2 "
            f"FROM ({query}) q ORDER BY genome"
        )
        cp._duckdb_copy_query_to_parquet(con, out_query, output_file)
    finally:
        con.close()
