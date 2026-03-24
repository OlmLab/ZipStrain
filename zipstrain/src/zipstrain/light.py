"""Lightweight profile/compare pipeline backed by DuckDB tables.

This module is intentionally separate from the existing zipstrain CLI and keeps the
legacy behavior untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal, Optional, Union

import duckdb
import polars as pl

from zipstrain import compare as cp
from zipstrain import profile as pf


BASE_TO_MASK = {"A": 1, "T": 2, "C": 4, "G": 8}
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
NUMERIC_COLUMN_DTYPES = {
    "total_positions": pl.Int64,
    "share_allele_pos": pl.Int64,
    "genome_pop_ani": pl.Float64,
    "max_consecutive_length": pl.Int64,
    "shared_genes_count": pl.Int64,
    "identical_gene_count": pl.Int64,
    "perc_id_genes": pl.Float64,
}
PROFILE_DUCKDB_FILENAME = "profile.duckdb"
PROFILE_COVERAGE_FILENAME = "coverage.parquet"
PROFILE_SNP_FILENAME = "snp.parquet"


@dataclass(frozen=True)
class LightProfileSummary:
    """Summary of a generated light profile database."""

    coverage_rows: int
    snp_rows: int


@dataclass(frozen=True)
class _LightProfileInput:
    """Resolved light profile input source."""

    format: Literal["duckdb", "parquet"]
    root: Path
    duckdb_path: Optional[Path] = None
    coverage_path: Optional[Path] = None
    snp_path: Optional[Path] = None


def _parse_calculate(calculate: str) -> tuple[str, ...]:
    """Parse modular calculation selection string.

    Supported canonical calculation groups:
      - ani
      - ibs
      - genes

    Input format:
      - "all"
      - "<group>"
      - "<group>+<group>+..."
    """
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


def _parse_reference_fasta(reference_fasta: Union[str, Path]) -> dict[str, str]:
    """Read FASTA records into an uppercase scaffold->sequence dictionary."""
    fasta_path = Path(reference_fasta)
    refs: dict[str, list[str]] = {}
    current: Optional[str] = None

    with fasta_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                refs[current] = []
                continue
            if current is None:
                raise ValueError(f"Malformed FASTA {fasta_path}: sequence encountered before header")
            refs[current].append(line.upper())

    return {name: "".join(chunks) for name, chunks in refs.items()}


def _ref_mask_for_position(ref_sequences: dict[str, str], chrom: str, pos: int) -> int:
    """Return one-hot base mask for a 1-based reference position."""
    seq = ref_sequences.get(chrom)
    if seq is None:
        return 0
    idx = int(pos) - 1
    if idx < 0 or idx >= len(seq):
        return 0
    return BASE_TO_MASK.get(seq[idx], 0)


def _profile_with_masks(profile_parquet: Union[str, Path], reference_fasta: Union[str, Path]) -> pl.LazyFrame:
    """Load profile parquet and derive coverage/allele masks plus reference masks."""
    refs = _parse_reference_fasta(reference_fasta)

    lf = (
        pl.scan_parquet(profile_parquet)
        .select("chrom", "pos", "gene", "genome", "A", "T", "C", "G")
        .with_columns(
            pl.col("chrom").cast(pl.Utf8),
            pl.col("gene").cast(pl.Utf8),
            pl.col("genome").cast(pl.Utf8),
            pl.col("pos").cast(pl.Int32),
            cov=(pl.col("A") + pl.col("T") + pl.col("C") + pl.col("G")).cast(pl.Int32),
            allele_mask=(
                (pl.col("A") > 0).cast(pl.UInt8) * pl.lit(BASE_TO_MASK["A"], dtype=pl.UInt8)
                + (pl.col("T") > 0).cast(pl.UInt8) * pl.lit(BASE_TO_MASK["T"], dtype=pl.UInt8)
                + (pl.col("C") > 0).cast(pl.UInt8) * pl.lit(BASE_TO_MASK["C"], dtype=pl.UInt8)
                + (pl.col("G") > 0).cast(pl.UInt8) * pl.lit(BASE_TO_MASK["G"], dtype=pl.UInt8)
            ).cast(pl.UInt8),
        )
        .with_columns(
            ref_mask=pl.struct(["chrom", "pos"]).map_elements(
                lambda row: _ref_mask_for_position(refs, str(row["chrom"]), int(row["pos"])),
                return_dtype=pl.UInt8,
            )
        )
    )
    return lf


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
    """Check whether an attached DuckDB schema contains a table."""
    row = conn.execute(
        """
        SELECT COUNT(*)
        FROM information_schema.tables
        WHERE table_schema = ? AND table_name = ?
        """,
        [schema_name, table_name],
    ).fetchone()
    return bool(row and row[0] > 0)


def _split_light_tables(
    profile_parquet: Union[str, Path],
    reference_fasta: Union[str, Path],
    min_cov: int = 5,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Build coverage and SNP tables from a regular ZipStrain profile parquet.

    Coverage rows are retained only for loci with `cov > min_cov`.
    SNP rows are called only on that filtered coverage set.
    """
    profile_lf = _profile_with_masks(profile_parquet=profile_parquet, reference_fasta=reference_fasta)
    filtered_lf = profile_lf.filter(pl.col("cov") > int(min_cov))
    coverage_lf = filtered_lf.select("genome", "chrom", "pos", "gene", "cov", "ref_mask")
    snp_lf = filtered_lf.filter(pl.col("allele_mask") != pl.col("ref_mask")).select(
        "genome", "chrom", "pos", "gene", "allele_mask"
    )
    return coverage_lf, snp_lf


def _write_light_duckdb_profile(
    profile_parquet: Union[str, Path],
    reference_fasta: Union[str, Path],
    min_cov: int,
    output_db: Union[str, Path],
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Write light profile tables to a DuckDB database file from profile parquet."""
    output_path = Path(output_db)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    refs = _parse_reference_fasta(reference_fasta)
    if not refs:
        raise ValueError(f"No reference sequences found in FASTA: {reference_fasta}")

    profile_sql = cp._duckdb_quote_sql_string(str(profile_parquet))
    conn = duckdb.connect(str(output_path))
    try:
        _configure_duckdb_connection(conn, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        conn.execute("CREATE TEMP TABLE ref_sequences(chrom VARCHAR, seq VARCHAR)")
        conn.executemany(
            "INSERT INTO ref_sequences VALUES (?, ?)",
            [(chrom, seq) for chrom, seq in refs.items()],
        )
        conn.execute(
            f"""
            CREATE TEMP TABLE profile_masks AS
            SELECT
              CAST(p.genome AS VARCHAR) AS genome,
              CAST(p.chrom AS VARCHAR) AS chrom,
              CAST(p.pos AS INTEGER) AS pos,
              CAST(p.gene AS VARCHAR) AS gene,
              (CAST(p.A AS INTEGER) + CAST(p.T AS INTEGER) + CAST(p.C AS INTEGER) + CAST(p.G AS INTEGER)) AS cov,
              (
                CASE WHEN CAST(p.A AS INTEGER) > 0 THEN 1 ELSE 0 END
                + CASE WHEN CAST(p.T AS INTEGER) > 0 THEN 2 ELSE 0 END
                + CASE WHEN CAST(p.C AS INTEGER) > 0 THEN 4 ELSE 0 END
                + CASE WHEN CAST(p.G AS INTEGER) > 0 THEN 8 ELSE 0 END
              )::UTINYINT AS allele_mask,
              (
                CASE UPPER(SUBSTR(COALESCE(r.seq, ''), CAST(p.pos AS INTEGER), 1))
                  WHEN 'A' THEN 1
                  WHEN 'T' THEN 2
                  WHEN 'C' THEN 4
                  WHEN 'G' THEN 8
                  ELSE 0
                END
              )::UTINYINT AS ref_mask
            FROM read_parquet('{profile_sql}') p
            LEFT JOIN ref_sequences r
              ON CAST(p.chrom AS VARCHAR) = r.chrom
            WHERE (CAST(p.A AS INTEGER) + CAST(p.T AS INTEGER) + CAST(p.C AS INTEGER) + CAST(p.G AS INTEGER)) > {int(min_cov)}
            """
        )
        conn.execute("CREATE TABLE coverage AS SELECT genome, chrom, pos, gene, cov, ref_mask FROM profile_masks")
        conn.execute(
            "CREATE TABLE snp AS "
            "SELECT genome, chrom, pos, gene, allele_mask "
            "FROM profile_masks WHERE allele_mask <> ref_mask"
        )
        conn.execute("CREATE TABLE genomes AS SELECT DISTINCT CAST(genome AS VARCHAR) AS genome FROM coverage")
        conn.execute(
            "CREATE TABLE metadata AS SELECT ?::VARCHAR AS reference_fasta, ?::INTEGER AS min_cov",
            [str(reference_fasta), int(min_cov)],
        )
        coverage_rows = int(conn.execute("SELECT COUNT(*) FROM coverage").fetchone()[0])
        snp_rows = int(conn.execute("SELECT COUNT(*) FROM snp").fetchone()[0])
    finally:
        conn.close()
    return LightProfileSummary(coverage_rows=coverage_rows, snp_rows=snp_rows)


def _resolve_light_profile_input(profile_input: Union[str, Path]) -> _LightProfileInput:
    """Resolve a light profile input directory to either DuckDB or parquet source."""
    path = Path(profile_input)
    if not path.exists():
        raise FileNotFoundError(f"Light profile path does not exist: {path}")

    if path.is_file():
        raise ValueError(
            f"Unsupported light profile path: {path}. "
            "Provide a profile directory, not a file."
        )

    if not path.is_dir():
        raise ValueError(f"Unsupported light profile path: {path}")

    db_file = path / PROFILE_DUCKDB_FILENAME
    coverage_file = path / PROFILE_COVERAGE_FILENAME
    snp_file = path / PROFILE_SNP_FILENAME
    duckdb_files = sorted(path.glob("*.duckdb"))

    has_db = db_file.exists()
    has_cov = coverage_file.exists()
    has_snp = snp_file.exists()

    if has_db and not has_cov and not has_snp:
        extra_db_files = [p for p in duckdb_files if p.name != PROFILE_DUCKDB_FILENAME]
        if extra_db_files:
            extras = ", ".join(sorted(p.name for p in extra_db_files))
            raise ValueError(
                f"Ambiguous light profile directory {path}: "
                f"only {PROFILE_DUCKDB_FILENAME} is allowed, found extra db files ({extras})."
            )
        return _LightProfileInput(format="duckdb", root=path, duckdb_path=db_file)

    if has_cov and has_snp and not duckdb_files:
        return _LightProfileInput(
            format="parquet",
            root=path,
            coverage_path=coverage_file,
            snp_path=snp_file,
        )

    if has_cov != has_snp:
        missing = PROFILE_SNP_FILENAME if has_cov else PROFILE_COVERAGE_FILENAME
        raise ValueError(f"Incomplete light profile directory {path}: missing {missing}")

    if duckdb_files and not db_file.exists():
        extras = ", ".join(sorted(p.name for p in duckdb_files))
        raise ValueError(
            f"Invalid light profile directory {path}: expected {PROFILE_DUCKDB_FILENAME}, "
            f"found other db files ({extras})."
        )

    raise ValueError(
        f"Unrecognized light profile directory {path}. "
        f"Expected either {PROFILE_DUCKDB_FILENAME} or "
        f"{PROFILE_COVERAGE_FILENAME}+{PROFILE_SNP_FILENAME}."
    )


def _register_light_profile_source(
    conn: duckdb.DuckDBPyConnection,
    source: _LightProfileInput,
    schema_name: str,
) -> bool:
    """Expose a profile source under `schema_name` as coverage/snp/genomes tables."""
    if source.format == "duckdb":
        db_sql = cp._duckdb_quote_sql_string(str(source.duckdb_path))
        conn.execute(f"ATTACH '{db_sql}' AS {schema_name} (READ_ONLY)")
        return _attached_has_table(conn, schema_name, "genomes")

    coverage_sql = cp._duckdb_quote_sql_string(str(source.coverage_path))
    snp_sql = cp._duckdb_quote_sql_string(str(source.snp_path))
    conn.execute(f"CREATE SCHEMA {schema_name}")
    conn.execute(
        f"CREATE VIEW {schema_name}.coverage AS "
        f"SELECT * FROM read_parquet('{coverage_sql}')"
    )
    conn.execute(
        f"CREATE VIEW {schema_name}.snp AS "
        f"SELECT * FROM read_parquet('{snp_sql}')"
    )
    conn.execute(
        f"CREATE VIEW {schema_name}.genomes AS "
        f"SELECT DISTINCT CAST(genome AS VARCHAR) AS genome FROM {schema_name}.coverage"
    )
    return True


def build_light_profile_bundle_from_parquet(
    profile_parquet: Union[str, Path],
    reference_fasta: Union[str, Path],
    output_dir: Union[str, Path],
    profile_engine: Literal["duckdb", "polars"] = "duckdb",
    min_cov: int = 5,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Create a light profile directory from a regular profile parquet.

    `profile_engine=duckdb` writes one file:
      - `profile.duckdb`

    `profile_engine=polars` writes two files:
      - `coverage.parquet`
      - `snp.parquet`
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keep profile directory deterministic across reruns.
    output_db = out_dir / PROFILE_DUCKDB_FILENAME
    if output_db.exists():
        output_db.unlink()
    for parquet_file in (out_dir / PROFILE_COVERAGE_FILENAME, out_dir / PROFILE_SNP_FILENAME):
        if parquet_file.exists():
            parquet_file.unlink()

    if profile_engine == "duckdb":
        extra_db_files = [p for p in out_dir.glob("*.duckdb") if p.name != PROFILE_DUCKDB_FILENAME]
        if extra_db_files:
            extras = ", ".join(sorted(str(p.name) for p in extra_db_files))
            raise ValueError(
                f"Output directory {out_dir} contains unrelated DuckDB files ({extras}). "
                f"Use an empty directory or remove them."
            )
        return _write_light_duckdb_profile(
            profile_parquet=profile_parquet,
            reference_fasta=reference_fasta,
            min_cov=min_cov,
            output_db=output_db,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
    elif profile_engine == "polars":
        coverage_lf, snp_lf = _split_light_tables(
            profile_parquet=profile_parquet,
            reference_fasta=reference_fasta,
            min_cov=min_cov,
        )
        coverage_path = out_dir / PROFILE_COVERAGE_FILENAME
        snp_path = out_dir / PROFILE_SNP_FILENAME
        coverage_lf.sink_parquet(coverage_path, compression="zstd", engine="streaming")
        snp_lf.sink_parquet(snp_path, compression="zstd", engine="streaming")
        coverage_rows = int(pl.scan_parquet(coverage_path).select(pl.len()).collect().item())
        snp_rows = int(pl.scan_parquet(snp_path).select(pl.len()).collect().item())
    else:
        raise ValueError(f"Unsupported light profile engine: {profile_engine}")

    return LightProfileSummary(
        coverage_rows=coverage_rows,
        snp_rows=snp_rows,
    )


def build_light_profile_bundle_from_bam(
    bam_file: Union[str, Path],
    bed_file: Union[str, Path],
    gene_range_table: Union[str, Path],
    stb_file: Union[str, Path],
    null_model: Union[str, Path],
    reference_fasta: Union[str, Path],
    output_dir: Union[str, Path],
    profile_engine: Literal["duckdb", "polars"] = "duckdb",
    min_cov: int = 5,
    num_chunks: int = 24,
    max_concurrency: int = 4,
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
) -> LightProfileSummary:
    """Profile a BAM and write a light profile directory."""
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
        profile_parquet = tmp_path / f"{Path(bam_file).stem}_profile.parquet"
        if not profile_parquet.exists():
            raise FileNotFoundError(f"Expected profile output not found: {profile_parquet}")
        return build_light_profile_bundle_from_parquet(
            profile_parquet=profile_parquet,
            reference_fasta=reference_fasta,
            output_dir=output_dir,
            profile_engine=profile_engine,
            min_cov=min_cov,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )


def _light_shared_query(min_cov: int, genome_scope: str = "all") -> str:
    """Build shared-loci query for light profile DBs attached as p1/p2."""
    genome_scope_sql = cp._duckdb_quote_sql_string(genome_scope)
    return f"""
    WITH shared_cov AS (
      SELECT
        c1.genome,
        c1.chrom AS scaffold,
        c1.pos,
        c1.gene,
        c1.ref_mask
      FROM p1.coverage c1
      INNER JOIN p2.coverage c2
        ON c1.genome = c2.genome
       AND c1.chrom = c2.chrom
       AND c1.pos = c2.pos
      WHERE c1.cov >= {int(min_cov)}
        AND c2.cov >= {int(min_cov)}
        AND ('{genome_scope_sql}' = 'all' OR c1.genome = '{genome_scope_sql}')
    ),
    snp_union AS (
      SELECT
        sc.genome,
        sc.scaffold,
        sc.pos,
        COALESCE(s1.allele_mask, sc.ref_mask) AS mask_1,
        COALESCE(s2.allele_mask, sc.ref_mask) AS mask_2
      FROM shared_cov sc
      LEFT JOIN p1.snp s1
        ON sc.genome = s1.genome
       AND sc.scaffold = s1.chrom
       AND sc.pos = s1.pos
      LEFT JOIN p2.snp s2
        ON sc.genome = s2.genome
       AND sc.scaffold = s2.chrom
       AND sc.pos = s2.pos
      WHERE s1.pos IS NOT NULL OR s2.pos IS NOT NULL
    ),
    disjoint AS (
      SELECT genome, scaffold, pos
      FROM snp_union
      WHERE (CAST(mask_1 AS INTEGER) & CAST(mask_2 AS INTEGER)) = 0
    )
    SELECT
      sc.genome,
      sc.scaffold,
      sc.pos,
      sc.gene,
      CASE WHEN d.pos IS NULL THEN 1 ELSE 0 END AS surr
    FROM shared_cov sc
    LEFT JOIN disjoint d
      ON sc.genome = d.genome
     AND sc.scaffold = d.scaffold
     AND sc.pos = d.pos
    """


def _light_genomes_scope_cte(
    genome_scope_sql: str,
    stb_file: Union[str, Path],
    use_embedded_genomes: bool = True,
) -> str:
    """Build genomes scope CTE from STB for zipstrain-compatible genome universe."""
    stb_sql = cp._duckdb_quote_sql_string(str(stb_file))
    return cp._duckdb_genomes_scope_cte(stb_sql=stb_sql, genome_scope_sql=genome_scope_sql)


def _light_ani_only_ctes(min_cov: int, genome_scope_sql: str) -> list[str]:
    """CTEs for ANI-only mode without materializing per-position surr rows."""
    return [
        f"""
shared_cov AS (
  SELECT
    c1.genome,
    c1.chrom AS scaffold,
    c1.pos,
    c1.ref_mask
  FROM p1.coverage c1
  INNER JOIN p2.coverage c2
    ON c1.genome = c2.genome
   AND c1.chrom = c2.chrom
   AND c1.pos = c2.pos
  WHERE c1.cov >= {int(min_cov)}
    AND c2.cov >= {int(min_cov)}
    AND ('{genome_scope_sql}' = 'all' OR c1.genome = '{genome_scope_sql}')
)""".strip(),
        """
total_cov AS (
  SELECT
    genome,
    COUNT(*)::BIGINT AS total_positions
  FROM shared_cov
  GROUP BY genome
)""".strip(),
        """
disjoint_counts AS (
  SELECT
    sc.genome,
    COUNT(*)::BIGINT AS disjoint_positions
  FROM shared_cov sc
  LEFT JOIN p1.snp s1
    ON sc.genome = s1.genome
   AND sc.scaffold = s1.chrom
   AND sc.pos = s1.pos
  LEFT JOIN p2.snp s2
    ON sc.genome = s2.genome
   AND sc.scaffold = s2.chrom
   AND sc.pos = s2.pos
  WHERE (s1.pos IS NOT NULL OR s2.pos IS NOT NULL)
    AND (
      CAST(COALESCE(s1.allele_mask, sc.ref_mask) AS INTEGER)
      & CAST(COALESCE(s2.allele_mask, sc.ref_mask) AS INTEGER)
    ) = 0
  GROUP BY sc.genome
)""".strip(),
        """
pop AS (
  SELECT
    t.genome,
    t.total_positions,
    (t.total_positions - COALESCE(d.disjoint_positions, 0))::BIGINT AS share_allele_pos,
    (t.total_positions - COALESCE(d.disjoint_positions, 0)) * 100.0 / NULLIF(t.total_positions, 0) AS genome_pop_ani
  FROM total_cov t
  LEFT JOIN disjoint_counts d ON t.genome = d.genome
)""".strip(),
    ]


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


def _light_pop_cte(source: str) -> str:
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


def _light_gene_ctes(min_gene_compare_len: int, contig_source: str = "contig") -> list[str]:
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
    min_cov: int,
    genome_scope: str,
    min_gene_compare_len: int,
    calculations: tuple[str, ...],
    stb_file: Union[str, Path],
    use_embedded_genomes: bool = True,
) -> str:
    genome_scope_sql = cp._duckdb_quote_sql_string(genome_scope)
    genomes_cte = _light_genomes_scope_cte(
        genome_scope_sql=genome_scope_sql,
        stb_file=stb_file,
        use_embedded_genomes=use_embedded_genomes,
    )
    if calculations == ("ani",):
        ctes = [genomes_cte, *_light_ani_only_ctes(min_cov=min_cov, genome_scope_sql=genome_scope_sql)]
        return cp._duckdb_build_query_with_ctes(ctes, _light_final_select_sql(calculations=calculations))

    shared_query = _light_shared_query(min_cov=min_cov, genome_scope=genome_scope)
    needs_contig = "ibs" in calculations or "genes" in calculations
    ctes = [
        f"shared AS (\n{shared_query}\n)",
        genomes_cte,
    ]

    if needs_contig:
        ctes.extend(_light_contig_ctes(shared_source="shared"))

    if "ani" in calculations:
        pop_source = "contig" if needs_contig else "shared"
        ctes.append(_light_pop_cte(source=pop_source))

    if "ibs" in calculations:
        ctes.append(_light_max_blocks_cte(contig_source="contig"))

    if "genes" in calculations:
        ctes.extend(_light_gene_ctes(min_gene_compare_len=min_gene_compare_len, contig_source="contig"))

    return cp._duckdb_build_query_with_ctes(ctes, _light_final_select_sql(calculations=calculations))


def _light_shared_polars(
    profile_1: Union[str, Path],
    profile_2: Union[str, Path],
    min_cov: int,
    genome_scope: str,
    memory_limit: Optional[str],
    temp_directory: Optional[Union[str, Path]],
    threads: Optional[int],
) -> pl.LazyFrame:
    """Build shared loci (`surr`, scaffold, pos, gene, genome) using Polars joins.

    Data extraction is scope-filtered and projected in DuckDB before transfer
    into Polars.
    """
    con = duckdb.connect()
    try:
        _configure_duckdb_connection(con, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        source_1 = _resolve_light_profile_input(profile_1)
        source_2 = _resolve_light_profile_input(profile_2)
        genome_scope_sql = cp._duckdb_quote_sql_string(genome_scope)
        _register_light_profile_source(con, source_1, "p1")
        _register_light_profile_source(con, source_2, "p2")

        shared_cov = pl.from_arrow(
            con.execute(
                f"""
                SELECT
                  c1.genome,
                  c1.chrom AS scaffold,
                  c1.pos,
                  c1.gene,
                  c1.ref_mask
                FROM p1.coverage c1
                INNER JOIN p2.coverage c2
                  ON c1.genome = c2.genome
                 AND c1.chrom = c2.chrom
                 AND c1.pos = c2.pos
                WHERE c1.cov >= {int(min_cov)}
                  AND c2.cov >= {int(min_cov)}
                  AND ('{genome_scope_sql}' = 'all' OR c1.genome = '{genome_scope_sql}')
                """
            ).fetch_arrow_table()
        ).lazy()

        snp_1 = pl.from_arrow(
            con.execute(
                f"""
                SELECT
                  genome,
                  chrom AS scaffold,
                  pos,
                  allele_mask AS allele_mask_1
                FROM p1.snp
                WHERE ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
                """
            ).fetch_arrow_table()
        ).lazy()

        snp_2 = pl.from_arrow(
            con.execute(
                f"""
                SELECT
                  genome,
                  chrom AS scaffold,
                  pos,
                  allele_mask AS allele_mask_2
                FROM p2.snp
                WHERE ('{genome_scope_sql}' = 'all' OR genome = '{genome_scope_sql}')
                """
            ).fetch_arrow_table()
        ).lazy()
    finally:
        con.close()

    return (
        shared_cov
        .join(snp_1, on=["genome", "scaffold", "pos"], how="left")
        .join(snp_2, on=["genome", "scaffold", "pos"], how="left")
        .with_columns(
            mask_1=pl.coalesce([pl.col("allele_mask_1"), pl.col("ref_mask")]).cast(pl.Int32),
            mask_2=pl.coalesce([pl.col("allele_mask_2"), pl.col("ref_mask")]).cast(pl.Int32),
        )
        .with_columns(
            surr=pl.when(
                pl.col("allele_mask_1").is_null() & pl.col("allele_mask_2").is_null()
            )
            .then(pl.lit(1))
            .otherwise(((pl.col("mask_1") & pl.col("mask_2")) > 0).cast(pl.Int8))
        )
        .select("genome", "scaffold", "pos", "gene", "surr")
    )


def _light_genomes_polars(
    genome_scope: str,
    stb_file: Union[str, Path],
) -> pl.LazyFrame:
    """Load comparison genome list from STB (zipstrain-compatible behavior)."""
    genomes = (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(pl.col("column_2").cast(pl.Utf8).alias("genome"))
        .unique()
    )
    if genome_scope != "all":
        genomes = genomes.filter(pl.col("genome") == genome_scope)
    return genomes


def _compare_light_genomes_polars(
    profile_1: Union[str, Path],
    profile_2: Union[str, Path],
    stb_file: Union[str, Path],
    min_cov: int,
    min_gene_compare_len: int,
    genome_scope: str,
    calculations: tuple[str, ...],
    memory_limit: Optional[str],
    temp_directory: Optional[Union[str, Path]],
    threads: Optional[int],
) -> pl.LazyFrame:
    shared = _light_shared_polars(
        profile_1=profile_1,
        profile_2=profile_2,
        min_cov=min_cov,
        genome_scope=genome_scope,
        memory_limit=memory_limit,
        temp_directory=temp_directory,
        threads=threads,
    )
    needs_contig = "ibs" in calculations or "genes" in calculations
    contig = cp.add_contiguity_info(shared) if needs_contig else None

    metric_frames: list[pl.LazyFrame] = []
    if "ani" in calculations:
        metric_frames.append(cp.calculate_pop_ani(shared if contig is None else contig))
    if "ibs" in calculations:
        metric_frames.append(cp.get_longest_consecutive_blocks(contig))
    if "genes" in calculations:
        metric_frames.append(cp.get_gene_ani(contig, min_gene_compare_len))

    genome_comp = metric_frames[0]
    for frame in metric_frames[1:]:
        genome_comp = genome_comp.join(frame, on="genome", how="left")

    genomes = _light_genomes_polars(
        genome_scope=genome_scope,
        stb_file=stb_file,
    )

    output_columns = _selected_output_columns(calculations)
    fill_exprs = [
        pl.col(col).fill_null(0 if NUMERIC_COLUMN_DTYPES[col] == pl.Int64 else 0.0).cast(NUMERIC_COLUMN_DTYPES[col])
        for col in output_columns
        if col != "genome"
    ]
    return genomes.join(genome_comp, on="genome", how="left").with_columns(fill_exprs).select(output_columns)


def _compare_light_genomes_duckdb(
    profile_1: Union[str, Path],
    profile_2: Union[str, Path],
    stb_file: Union[str, Path],
    min_cov: int,
    min_gene_compare_len: int,
    genome_scope: str,
    calculations: tuple[str, ...],
    memory_limit: Optional[str],
    temp_directory: Optional[Union[str, Path]],
    threads: Optional[int],
) -> pl.LazyFrame:
    """Compare two light profile DBs and return genome-level metrics (DuckDB engine)."""
    con = duckdb.connect()
    try:
        _configure_duckdb_connection(con, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
        source_1 = _resolve_light_profile_input(profile_1)
        source_2 = _resolve_light_profile_input(profile_2)
        has_genomes_1 = _register_light_profile_source(con, source_1, "p1")
        has_genomes_2 = _register_light_profile_source(con, source_2, "p2")
        use_embedded_genomes = has_genomes_1 and has_genomes_2

        query = _light_duckdb_compare_query(
            min_cov=min_cov,
            genome_scope=genome_scope,
            min_gene_compare_len=min_gene_compare_len,
            calculations=calculations,
            stb_file=stb_file,
            use_embedded_genomes=use_embedded_genomes,
        )
        table = con.execute(query).fetch_arrow_table()
        output_columns = _selected_output_columns(calculations)
        return pl.from_arrow(table).lazy().select(output_columns)
    finally:
        con.close()


def compare_light_genomes(
    profile_db_1: Union[str, Path],
    profile_db_2: Union[str, Path],
    stb_file: Union[str, Path],
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
    engine: Literal["duckdb", "polars"] = "duckdb",
    calculate: str = "all",
) -> pl.LazyFrame:
    """Compare two light profiles and return genome-level metrics.

    Each profile input must be a directory containing either:
      - `profile.duckdb`
      - `coverage.parquet` and `snp.parquet`
    """
    if not stb_file:
        raise ValueError("stb_file is required for zipstrain-light compare")
    if ani_method != "popani":
        raise ValueError("zipstrain-light currently supports only --ani-method popani")
    calculations = _parse_calculate(calculate)
    if engine == "duckdb":
        return _compare_light_genomes_duckdb(
            profile_1=profile_db_1,
            profile_2=profile_db_2,
            stb_file=stb_file,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            calculations=calculations,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
    if engine == "polars":
        return _compare_light_genomes_polars(
            profile_1=profile_db_1,
            profile_2=profile_db_2,
            stb_file=stb_file,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            calculations=calculations,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
        )
    raise ValueError(f"Unsupported zipstrain-light engine: {engine}")


def compare_light_genomes_to_parquet(
    profile_db_1: Union[str, Path],
    profile_db_2: Union[str, Path],
    stb_file: Union[str, Path],
    output_file: Union[str, Path],
    sample_1_name: str,
    sample_2_name: str,
    min_cov: int = 5,
    min_gene_compare_len: int = 100,
    genome_scope: str = "all",
    ani_method: str = "popani",
    memory_limit: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
    engine: Literal["duckdb", "polars"] = "duckdb",
    calculate: str = "all",
) -> None:
    """Compare two light profiles and write final parquet output."""
    if not stb_file:
        raise ValueError("stb_file is required for zipstrain-light compare")
    if ani_method != "popani":
        raise ValueError("zipstrain-light currently supports only --ani-method popani")
    calculations = _parse_calculate(calculate)

    if engine == "duckdb":
        con = duckdb.connect()
        try:
            _configure_duckdb_connection(con, memory_limit=memory_limit, temp_directory=temp_directory, threads=threads)
            source_1 = _resolve_light_profile_input(profile_db_1)
            source_2 = _resolve_light_profile_input(profile_db_2)
            has_genomes_1 = _register_light_profile_source(con, source_1, "p1")
            has_genomes_2 = _register_light_profile_source(con, source_2, "p2")
            use_embedded_genomes = has_genomes_1 and has_genomes_2
            query = _light_duckdb_compare_query(
                min_cov=min_cov,
                genome_scope=genome_scope,
                min_gene_compare_len=min_gene_compare_len,
                calculations=calculations,
                stb_file=stb_file,
                use_embedded_genomes=use_embedded_genomes,
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
            return
        finally:
            con.close()

    (
        compare_light_genomes(
            profile_db_1=profile_db_1,
            profile_db_2=profile_db_2,
            stb_file=stb_file,
            min_cov=min_cov,
            min_gene_compare_len=min_gene_compare_len,
            genome_scope=genome_scope,
            ani_method=ani_method,
            memory_limit=memory_limit,
            temp_directory=temp_directory,
            threads=threads,
            engine=engine,
            calculate=calculate,
        )
        .with_columns(
            pl.lit(sample_1_name).alias("sample_1"),
            pl.lit(sample_2_name).alias("sample_2"),
        )
        .sink_parquet(output_file, compression="zstd")
    )
