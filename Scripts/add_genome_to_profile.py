#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import uuid
from pathlib import Path

import duckdb


def _quote_sql_string(value: str) -> str:
    return value.replace("'", "''")


def _read_parquet_columns(conn: duckdb.DuckDBPyConnection, parquet_path: Path) -> list[str]:
    parquet_sql = _quote_sql_string(str(parquet_path))
    rows = conn.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{parquet_sql}')"
    ).fetchall()
    return [row[0] for row in rows]


def _stb_cte(stb_path: Path) -> str:
    stb_sql = _quote_sql_string(str(stb_path))
    return (
        "stb AS ("
        f"SELECT DISTINCT CAST(scaffold AS VARCHAR) AS scaffold, CAST(genome AS VARCHAR) AS genome "
        f"FROM read_csv('{stb_sql}', delim='\\t', header=false, columns={{'scaffold':'VARCHAR','genome':'VARCHAR'}})"
        ")"
    )


def _missing_scaffolds(conn: duckdb.DuckDBPyConnection, parquet_path: Path, stb_path: Path) -> list[str]:
    parquet_sql = _quote_sql_string(str(parquet_path))
    query = f"""
    WITH
      {_stb_cte(stb_path)}
    SELECT DISTINCT CAST(p.chrom AS VARCHAR) AS chrom
    FROM read_parquet('{parquet_sql}') AS p
    LEFT JOIN stb AS s
      ON CAST(p.chrom AS VARCHAR) = s.scaffold
    WHERE s.genome IS NULL
    ORDER BY chrom
    LIMIT 20
    """
    return [row[0] for row in conn.execute(query).fetchall()]


def _select_list(columns: list[str]) -> str:
    select_parts: list[str] = []
    inserted_genome = False
    for column in columns:
        if column == "chrom":
            select_parts.append("p.chrom")
            select_parts.append("CAST(s.genome AS VARCHAR) AS genome")
            inserted_genome = True
        else:
            select_parts.append(f"p.{column}")
    if not inserted_genome:
        raise ValueError("Profile parquet is missing required column 'chrom'.")
    return ",\n      ".join(select_parts)


def add_genome_column(profile_parquet: Path, stb_path: Path) -> int:
    profile_path = profile_parquet.expanduser().resolve()
    stb_file = stb_path.expanduser().resolve()

    if not profile_path.exists():
        raise FileNotFoundError(f"Profile parquet not found: {profile_path}")
    if not stb_file.exists():
        raise FileNotFoundError(f"STB file not found: {stb_file}")

    conn = duckdb.connect()
    try:
        columns = _read_parquet_columns(conn, profile_path)
        if "genome" in columns:
            print(f"{profile_path} already has a genome column; nothing to do.")
            return 0

        missing = _missing_scaffolds(conn, profile_path, stb_file)
        if missing:
            preview = ", ".join(missing)
            raise RuntimeError(
                "Could not map all scaffolds to genomes from the STB file. "
                f"Examples: {preview}"
            )

        tmp_path = profile_path.with_name(
            f"{profile_path.name}.tmp.{uuid.uuid4().hex}.parquet"
        )
        profile_sql = _quote_sql_string(str(profile_path))
        tmp_sql = _quote_sql_string(str(tmp_path))
        select_sql = _select_list(columns)
        query = f"""
        COPY (
          WITH
            {_stb_cte(stb_file)}
          SELECT
            {select_sql}
          FROM read_parquet('{profile_sql}') AS p
          LEFT JOIN stb AS s
            ON CAST(p.chrom AS VARCHAR) = s.scaffold
        ) TO '{tmp_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
        conn.execute(query)
    finally:
        conn.close()

    os.replace(tmp_path, profile_path)
    print(f"Added genome column to {profile_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Add a genome column to an old ZipStrain profile parquet using an STB file. "
            "If the profile already has a genome column, the script exits without changing it."
        )
    )
    parser.add_argument("profile_parquet", help="Path to the profile parquet to rewrite in place.")
    parser.add_argument("stb_file", help="Path to the STB file mapping scaffold to genome.")
    args = parser.parse_args()

    try:
        return add_genome_column(Path(args.profile_parquet), Path(args.stb_file))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
