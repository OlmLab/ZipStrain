"""zipstrain.build_db
=====================
Helpers to build a local genome database from abundance outputs produced by
external tools (for example, Sylph).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import gzip
import pathlib
import random
import re
import shutil
import tempfile
import time
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import urlopen
import zipfile

import polars as pl


# Matches common assembly accessions, with optional version.
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d{6,}(?:\.\d+)?)", flags=re.IGNORECASE)

STRING_DTYPES = {pl.Utf8, pl.Categorical, pl.Enum}
NUMERIC_DTYPES = {
    pl.Int8,
    pl.Int16,
    pl.Int32,
    pl.Int64,
    pl.UInt8,
    pl.UInt16,
    pl.UInt32,
    pl.UInt64,
    pl.Float32,
    pl.Float64,
}

# event, fetched_count, remaining_count, total_count, accession(optional)
ProgressCallback = Callable[[str, int, int, int, str | None], None]


def _normalize_accession(accession: str) -> str:
    return accession.strip().upper()


def _clean_genome_name(raw_value: str) -> str:
    raw = raw_value.strip()
    if raw == "":
        return raw
    parsed = urlparse(raw)
    if parsed.scheme in {"http", "https"}:
        filename = pathlib.PurePosixPath(parsed.path).name
    else:
        filename = pathlib.Path(raw).name
    if filename == "":
        return raw
    return pathlib.Path(filename).stem


def _find_column_case_insensitive(columns: list[str], candidates: list[str]) -> str | None:
    lookup = {col.lower(): col for col in columns}
    for candidate in candidates:
        hit = lookup.get(candidate.lower())
        if hit is not None:
            return hit
    return None


def _abundance_table_base_dir(table_source: str | pathlib.Path | pl.DataFrame | pl.LazyFrame) -> pathlib.Path | None:
    if isinstance(table_source, (str, pathlib.Path)):
        return pathlib.Path(table_source).resolve().parent
    return None


def _sanitize_id(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return text.strip("_") or "unknown"


def _filter_nonzero_abundance_rows(table: pl.LazyFrame) -> pl.LazyFrame:
    """Keep rows where at least one numeric column has a value > 0.

    If no numeric columns exist, the input is returned unchanged.
    """
    schema = table.collect_schema()
    numeric_cols = [name for name, dtype in schema.items() if dtype in NUMERIC_DTYPES]
    if not numeric_cols:
        return table
    preferred = _find_column_case_insensitive(
        numeric_cols,
        [
            "eff_cov",
            "genome_cov",
            "abundance",
            "relative_abundance",
            "rel_abundance",
        ],
    )
    if preferred is not None:
        cols_to_check = [preferred]
    else:
        abundance_like_cols = [
            col
            for col in numeric_cols
            if ("abund" in col.lower()) or ("cov" in col.lower())
        ]
        cols_to_check = abundance_like_cols if abundance_like_cols else numeric_cols
    nonzero_exprs = [
        pl.col(col).cast(pl.Float64, strict=False).fill_null(0.0) > 0.0
        for col in cols_to_check
    ]
    return table.filter(pl.any_horizontal(nonzero_exprs))


def _extract_accessions_generic(
    table: pl.LazyFrame | pl.DataFrame,
    *,
    preferred_columns: list[str] | None = None,
    url_columns: list[str] | None = None,
) -> pl.DataFrame:
    """Extract genome accessions from a table using regex matching.

    Returns columns:
      - accession
      - genome_name
      - download_url
    """
    lf = table.lazy() if isinstance(table, pl.DataFrame) else table
    schema = lf.collect_schema()
    str_columns = [name for name, dtype in schema.items() if dtype in STRING_DTYPES]

    if preferred_columns:
        ordered = [c for c in preferred_columns if c in str_columns]
        ordered.extend([c for c in str_columns if c not in ordered])
        str_columns = ordered

    if not str_columns:
        return pl.DataFrame(
            {
                "accession": [],
                "genome_name": [],
                "download_url": [],
            },
            schema={
                "accession": pl.Utf8,
                "genome_name": pl.Utf8,
                "download_url": pl.Utf8,
            },
        )

    candidates: list[pl.DataFrame] = []
    for col_name in str_columns:
        frame = (
            lf.select(
                pl.col(col_name).cast(pl.Utf8).alias("raw_value"),
                pl.col(col_name)
                .cast(pl.Utf8)
                .str.extract(ACCESSION_PATTERN.pattern, 1)
                .str.to_uppercase()
                .alias("accession"),
            )
            .drop_nulls("accession")
            .collect(engine="streaming")
        )
        if frame.height > 0:
            candidates.append(frame)

    if not candidates:
        return pl.DataFrame(
            {
                "accession": [],
                "genome_name": [],
                "download_url": [],
            },
            schema={
                "accession": pl.Utf8,
                "genome_name": pl.Utf8,
                "download_url": pl.Utf8,
            },
        )

    merged = pl.concat(candidates).unique(subset=["accession"], keep="first")
    genomes = (
        merged.with_columns(
            pl.col("accession").str.strip_chars().str.to_uppercase(),
            pl.col("raw_value").map_elements(_clean_genome_name, return_dtype=pl.Utf8).alias("genome_name"),
        )
        .select("accession", "genome_name")
        .with_columns(pl.lit(None, dtype=pl.Utf8).alias("download_url"))
        .sort("accession")
    )

    if not url_columns:
        return genomes

    existing_cols = [c for c in url_columns if c in str_columns]
    if not existing_cols:
        return genomes

    urls: list[pl.DataFrame] = []
    for col_name in existing_cols:
        frame = (
            lf.select(
                pl.col(col_name).cast(pl.Utf8).alias("raw_url"),
                pl.col(col_name)
                .cast(pl.Utf8)
                .str.extract(ACCESSION_PATTERN.pattern, 1)
                .str.to_uppercase()
                .alias("accession"),
            )
            .filter(pl.col("raw_url").str.starts_with("http"))
            .drop_nulls("accession")
            .collect(engine="streaming")
        )
        if frame.height > 0:
            urls.append(frame)

    if not urls:
        return genomes

    url_df = (
        pl.concat(urls)
        .with_columns(pl.col("accession").map_elements(_normalize_accession, return_dtype=pl.Utf8))
        .rename({"raw_url": "download_url"})
        .select("accession", "download_url")
        .unique(subset=["accession"], keep="first")
    )
    return genomes.join(url_df, on="accession", how="left", suffix="_resolved").with_columns(
        pl.coalesce(["download_url_resolved", "download_url"]).alias("download_url")
    ).select("accession", "genome_name", "download_url")


class AbundanceToolAdapter(ABC):
    """Base adapter for parsing abundance outputs from a specific tool."""

    name: str

    @abstractmethod
    def extract_genomes(self, abundance_table: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
        """Return a dataframe with columns accession/genome_name/download_url."""


class SylphAdapter(AbundanceToolAdapter):
    name = "sylph"

    PREFERRED_COLUMNS = [
        "Genome_file",
        "genome_file",
        "genome",
        "Genome",
        "Reference",
        "reference",
        "target",
    ]

    URL_COLUMNS = [
        "FTP_download",
        "ftp_download",
        "download_url",
        "url",
    ]

    def extract_genomes(self, abundance_table: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
        lf = abundance_table.lazy() if isinstance(abundance_table, pl.DataFrame) else abundance_table
        schema = lf.collect_schema()
        columns = list(schema.keys())

        genome_col = _find_column_case_insensitive(columns, self.PREFERRED_COLUMNS)
        if genome_col is None:
            raise ValueError(
                "Could not find a Sylph genome column. Expected one of: "
                + ", ".join(self.PREFERRED_COLUMNS)
            )

        url_col = _find_column_case_insensitive(columns, self.URL_COLUMNS)

        select_exprs = [
            pl.col(genome_col).cast(pl.Utf8).alias("raw_value"),
            pl.col(genome_col).cast(pl.Utf8).alias("source_path"),
            pl.col(genome_col)
            .cast(pl.Utf8)
            .str.extract(ACCESSION_PATTERN.pattern, 1)
            .str.to_uppercase()
            .alias("accession"),
        ]
        if url_col is not None:
            select_exprs.append(pl.col(url_col).cast(pl.Utf8).alias("download_url"))
        else:
            select_exprs.append(pl.lit(None, dtype=pl.Utf8).alias("download_url"))

        extracted = (
            lf.select(*select_exprs)
            .drop_nulls("accession")
            .with_columns(
                pl.col("accession").str.strip_chars().str.to_uppercase(),
                pl.col("raw_value").map_elements(_clean_genome_name, return_dtype=pl.Utf8).alias("genome_name"),
            )
            .select("accession", "genome_name", "download_url", "source_path")
            .unique(subset=["accession"], keep="first")
            .sort("accession")
            .collect(engine="streaming")
        )
        return extracted


ADAPTERS: dict[str, type[AbundanceToolAdapter]] = {
    SylphAdapter.name: SylphAdapter,
}


def get_tool_adapter(tool_name: str) -> AbundanceToolAdapter:
    key = tool_name.strip().lower()
    if key not in ADAPTERS:
        raise ValueError(f"Unsupported tool '{tool_name}'. Supported tools: {', '.join(sorted(ADAPTERS))}")
    return ADAPTERS[key]()


def get_genome_accessions_from_sylph_output(sylph_output: pl.LazyFrame | pl.DataFrame) -> pl.LazyFrame:
    """Backward-compatible helper for Sylph accession extraction."""
    return SylphAdapter().extract_genomes(sylph_output).select("accession", "genome_name", "download_url").lazy()


class LocalGenomeDB:
    """Parquet-backed local genome database index."""

    DEFAULT_DB_NAME = ".genome_db.parquet"
    SCHEMA = {
        "accession": pl.Utf8,
        "genome_name": pl.Utf8,
        "location": pl.Utf8,
        "download_url": pl.Utf8,
        "source_tool": pl.Utf8,
        "exists": pl.Boolean,
    }

    def __init__(self, db_path: str | pathlib.Path):
        self.db_path = pathlib.Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        if self.db_path.exists():
            self._db = pl.read_parquet(self.db_path)
            self._db = self._db.select(list(self.SCHEMA.keys())).cast(self.SCHEMA)
        else:
            self._db = pl.DataFrame(schema=self.SCHEMA)
            self.save()
        self.sync()

    @property
    def db(self) -> pl.DataFrame:
        return self._db

    def save(self) -> None:
        self._db.write_parquet(self.db_path)

    def reload(self) -> None:
        self._db = pl.read_parquet(self.db_path).select(list(self.SCHEMA.keys())).cast(self.SCHEMA)
        self.sync()

    def add_or_update_genomes(self, genomes: pl.DataFrame, source_tool: str) -> None:
        required = {"accession", "genome_name", "download_url"}
        missing = required - set(genomes.columns)
        if missing:
            raise ValueError(f"Input genomes table is missing columns: {sorted(missing)}")

        incoming = (
            genomes.select("accession", "genome_name", "download_url")
            .with_columns(
                pl.col("accession").str.strip_chars().str.to_uppercase(),
                pl.col("genome_name").fill_null(pl.col("accession")).cast(pl.Utf8),
                pl.col("download_url").cast(pl.Utf8),
                pl.lit(source_tool).alias("source_tool"),
            )
            .unique(subset=["accession"], keep="first")
        )

        if self._db.height == 0:
            self._db = incoming.with_columns(
                pl.lit(None, dtype=pl.Utf8).alias("location"),
                pl.lit(False).alias("exists"),
            ).select(list(self.SCHEMA.keys()))
            self.save()
            self.sync()
            return

        existing = self._db
        new_only = incoming.join(existing.select("accession"), on="accession", how="anti")
        if new_only.height > 0:
            appended = new_only.with_columns(
                pl.lit(None, dtype=pl.Utf8).alias("location"),
                pl.lit(False).alias("exists"),
            ).select(list(self.SCHEMA.keys()))
            existing = pl.concat([existing, appended], how="vertical")

        updates = incoming.rename(
            {
                "genome_name": "genome_name_new",
                "download_url": "download_url_new",
                "source_tool": "source_tool_new",
            }
        )
        self._db = (
            existing.join(updates, on="accession", how="left")
            .with_columns(
                pl.coalesce(["genome_name_new", "genome_name"]).alias("genome_name"),
                pl.coalesce(["download_url_new", "download_url"]).alias("download_url"),
                pl.coalesce(["source_tool_new", "source_tool"]).alias("source_tool"),
            )
            .drop(["genome_name_new", "download_url_new", "source_tool_new"])
            .select(list(self.SCHEMA.keys()))
        )
        self.save()
        self.sync()

    def set_location(self, accession: str, location: str | pathlib.Path, *, download_url: str | None = None) -> None:
        acc = _normalize_accession(accession)
        loc = str(pathlib.Path(location))
        self._db = self._db.with_columns(
            pl.when(pl.col("accession") == acc).then(pl.lit(loc)).otherwise(pl.col("location")).alias("location"),
            pl.when(pl.col("accession") == acc).then(pl.lit(download_url)).otherwise(pl.col("download_url")).alias("download_url"),
        )
        self.sync()
        self.save()

    def sync(self) -> None:
        self._db = self._db.with_columns(
            pl.col("location")
            .map_elements(
                lambda p: pathlib.Path(p).exists() if p is not None and p != "" else False,
                return_dtype=pl.Boolean,
                skip_nulls=False,
            )
            .alias("exists")
        )

    def missing(self) -> pl.DataFrame:
        self.sync()
        return self._db.filter(~pl.col("exists"))


def _ncbi_datasets_genome_download_url(accession: str) -> str:
    acc = _normalize_accession(accession)
    return (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{acc}/download?include_annotation_type=GENOME_FASTA"
    )


def _write_response_to_file(url: str, destination: pathlib.Path, timeout_seconds: int = 120) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=timeout_seconds) as response, open(destination, "wb") as out:
        shutil.copyfileobj(response, out)


def _parse_retry_after_seconds(retry_after_value: str | None) -> float | None:
    if retry_after_value is None:
        return None
    value = str(retry_after_value).strip()
    if value == "":
        return None
    if value.isdigit():
        return float(value)
    try:
        dt = parsedate_to_datetime(value)
    except Exception:
        return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    seconds = (dt - datetime.now(timezone.utc)).total_seconds()
    return max(seconds, 0.0)


def _retry_after_from_rate_limit_error(exc: Exception) -> float | None:
    if not isinstance(exc, HTTPError) or exc.code != 429:
        return None
    headers = getattr(exc, "headers", None) or getattr(exc, "hdrs", None)
    header_value = headers.get("Retry-After") if headers is not None and hasattr(headers, "get") else None
    return _parse_retry_after_seconds(header_value)


def _is_transient_download_error(exc: Exception) -> bool:
    """Return True when a genome download failure is worth retrying.

    We retry only likely transient failures such as rate limiting, server-side
    HTTP errors, and clear network/timeout problems. Permanent failures such as
    missing/retired accessions (for example HTTP 404/410) should fail fast.
    """
    if isinstance(exc, HTTPError):
        return exc.code in {408, 425, 429, 500, 502, 503, 504}
    if isinstance(exc, URLError):
        return True
    if isinstance(exc, (TimeoutError, ConnectionError)):
        return True
    if isinstance(exc, OSError):
        message = str(exc).lower()
        transient_markers = (
            "temporary failure",
            "temporarily unavailable",
            "timed out",
            "timeout",
            "connection reset",
            "connection refused",
            "connection aborted",
            "network is unreachable",
            "name or service not known",
            "nodename nor servname provided",
            "no route to host",
            "remote end closed connection",
        )
        return any(marker in message for marker in transient_markers)
    return False


def _extract_first_fasta_from_zip(zip_path: pathlib.Path, destination_fasta: pathlib.Path) -> pathlib.Path:
    with zipfile.ZipFile(zip_path, "r") as zf:
        candidates = [
            name
            for name in zf.namelist()
            if name.lower().endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"))
        ]
        if not candidates:
            raise ValueError(f"No FASTA files found in archive: {zip_path}")
        first = sorted(candidates)[0]
        with zf.open(first, "r") as src:
            if first.lower().endswith(".gz"):
                with gzip.GzipFile(fileobj=src) as gz_src, open(destination_fasta, "wb") as out:
                    shutil.copyfileobj(gz_src, out)
            else:
                with open(destination_fasta, "wb") as out:
                    shutil.copyfileobj(src, out)
    return destination_fasta


def _download_single_genome(
    url: str,
    accession: str,
    genomes_dir: pathlib.Path,
    *,
    timeout_seconds: int = 120,
) -> pathlib.Path:
    genomes_dir.mkdir(parents=True, exist_ok=True)
    destination_fasta = genomes_dir / f"{_normalize_accession(accession)}.fna"

    with tempfile.TemporaryDirectory(prefix="zipstrain_db_dl_") as tmpdir:
        tmp_dir = pathlib.Path(tmpdir)
        if url.lower().endswith(".gz"):
            tmp_file = tmp_dir / "genome.fna.gz"
            _write_response_to_file(url, tmp_file, timeout_seconds=timeout_seconds)
            with gzip.open(tmp_file, "rb") as src, open(destination_fasta, "wb") as out:
                shutil.copyfileobj(src, out)
            return destination_fasta
        if ".zip" in url.lower() or "download?" in url.lower():
            tmp_file = tmp_dir / "genome.zip"
            _write_response_to_file(url, tmp_file, timeout_seconds=timeout_seconds)
            return _extract_first_fasta_from_zip(tmp_file, destination_fasta)
        tmp_file = tmp_dir / "genome.fna"
        _write_response_to_file(url, tmp_file, timeout_seconds=timeout_seconds)
        shutil.move(tmp_file, destination_fasta)
        return destination_fasta


def _default_url_resolver(accession: str, db_row: dict) -> str:
    if db_row.get("download_url"):
        return str(db_row["download_url"])
    return _ncbi_datasets_genome_download_url(accession)


def _resolve_local_source_path(source_path: str, base_dir: pathlib.Path | None = None) -> pathlib.Path | None:
    raw = source_path.strip()
    if raw == "":
        return None
    parsed = urlparse(raw)
    if parsed.scheme in {"http", "https"}:
        return None
    p = pathlib.Path(raw)
    if p.exists():
        return p
    if (not p.is_absolute()) and base_dir is not None:
        candidate = base_dir / p
        if candidate.exists():
            return candidate
    return None


def _cache_local_source_genome(
    source_path: pathlib.Path,
    accession: str,
    genomes_dir: pathlib.Path,
) -> pathlib.Path:
    genomes_dir.mkdir(parents=True, exist_ok=True)
    destination_fasta = genomes_dir / f"{_normalize_accession(accession)}.fna"
    suffix = source_path.suffix.lower()
    if suffix == ".gz":
        with gzip.open(source_path, "rb") as src, open(destination_fasta, "wb") as out:
            shutil.copyfileobj(src, out)
        return destination_fasta
    if suffix == ".zip":
        return _extract_first_fasta_from_zip(source_path, destination_fasta)
    with open(source_path, "rb") as src, open(destination_fasta, "wb") as out:
        shutil.copyfileobj(src, out)
    return destination_fasta


def stage_genomes_from_source_paths(
    local_db: LocalGenomeDB,
    genomes: pl.DataFrame,
    *,
    genomes_dir: str | pathlib.Path,
    base_dir: pathlib.Path | None = None,
) -> None:
    """Populate cache from local source paths in abundance tables when available."""
    if "source_path" not in genomes.columns:
        return
    missing_set = set(local_db.missing()["accession"].to_list())
    if not missing_set:
        return

    candidates = (
        genomes.select("accession", "source_path")
        .drop_nulls("source_path")
        .unique(subset=["accession"], keep="first")
        .iter_rows(named=True)
    )
    for row in candidates:
        accession = _normalize_accession(str(row["accession"]))
        if accession not in missing_set:
            continue
        source_value = str(row["source_path"])
        resolved_source = _resolve_local_source_path(source_value, base_dir=base_dir)
        if resolved_source is None:
            continue
        cached = _cache_local_source_genome(
            resolved_source,
            accession=accession,
            genomes_dir=pathlib.Path(genomes_dir),
        )
        local_db.set_location(accession, cached)
        missing_set.remove(accession)
        if not missing_set:
            break


def fetch_missing_genomes(
    local_db: LocalGenomeDB,
    genomes_dir: str | pathlib.Path,
    *,
    url_resolver: Callable[[str, dict], str] | None = None,
    downloader: Callable[[str, str, pathlib.Path], pathlib.Path] | None = None,
    overwrite: bool = False,
    only_accessions: list[str] | None = None,
    progress_callback: ProgressCallback | None = None,
    initial_fetched: int = 0,
    total_selected: int | None = None,
    max_download_attempts: int = 3,
    backoff_base_seconds: float = 1.0,
    download_workers: int = 4,
) -> pl.DataFrame:
    """Download missing genomes and update DB locations/existence.

    Returns a report with columns:
      accession, status, location, url, error

    Downloads are executed with a bounded thread pool (`download_workers`).
    """
    resolver = url_resolver or _default_url_resolver
    dl = downloader or _download_single_genome
    genomes_path = pathlib.Path(genomes_dir)
    if max_download_attempts < 1:
        raise ValueError("max_download_attempts must be >= 1")
    if download_workers < 1:
        raise ValueError("download_workers must be >= 1")

    records: list[dict[str, str | None]] = []
    missing_df = local_db.missing()
    if only_accessions is not None:
        normalized = [_normalize_accession(acc) for acc in only_accessions]
        missing_df = missing_df.filter(pl.col("accession").is_in(normalized))

    completed_count = initial_fetched
    selected_total = total_selected if total_selected is not None else initial_fetched + missing_df.height
    if progress_callback is not None:
        progress_callback(
            "start",
            completed_count,
            max(selected_total - completed_count, 0),
            selected_total,
            None,
        )

    records_by_index: dict[int, dict[str, str | None]] = {}
    pending_rows: list[tuple[int, dict[str, str | None]]] = []
    for row_idx, row in enumerate(missing_df.iter_rows(named=True)):
        accession = row["accession"]
        location = row["location"]
        if location and pathlib.Path(location).exists() and not overwrite:
            records_by_index[row_idx] = {
                "accession": accession,
                "status": "already_present",
                "location": location,
                "url": row.get("download_url"),
                "error": None,
            }
            completed_count += 1
            if progress_callback is not None:
                progress_callback(
                    "already_present",
                    completed_count,
                    max(selected_total - completed_count, 0),
                    selected_total,
                    accession,
                )
            continue
        pending_rows.append((row_idx, row))
    pending_rows_by_index = {idx: row for idx, row in pending_rows}

    def _download_with_retries(row: dict[str, str | None]) -> dict[str, str | int | None]:
        accession = _normalize_accession(str(row["accession"]))
        last_exc: Exception | None = None
        url: str | None = None
        default_rate_limit_wait_seconds = 15.0
        for attempt in range(1, max_download_attempts + 1):
            try:
                url = resolver(accession, row)
                out_file = dl(url, accession, genomes_path)
                return {
                    "accession": accession,
                    "status": "downloaded",
                    "location": str(out_file),
                    "url": url,
                    "error": None,
                    "retry_count": attempt - 1,
                }
            except Exception as exc:
                last_exc = exc
                if attempt >= max_download_attempts or not _is_transient_download_error(exc):
                    break
                if backoff_base_seconds > 0:
                    delay_seconds = backoff_base_seconds * (2 ** (attempt - 1))
                    retry_after_seconds = _retry_after_from_rate_limit_error(exc)
                    if retry_after_seconds is not None:
                        delay_seconds = max(delay_seconds, retry_after_seconds, default_rate_limit_wait_seconds)
                    jitter_seconds = random.uniform(0.0, min(1.0, max(delay_seconds * 0.1, 0.1)))
                    time.sleep(delay_seconds + jitter_seconds)
                else:
                    retry_after_seconds = _retry_after_from_rate_limit_error(exc)
                    if retry_after_seconds is not None:
                        jitter_seconds = random.uniform(0.0, 1.0)
                        time.sleep(max(retry_after_seconds, default_rate_limit_wait_seconds) + jitter_seconds)

        failed_url = url if url is not None else row.get("download_url")
        attempts_used = 0 if last_exc is None else min(max_download_attempts, attempt)
        failed_error = (
            f"{last_exc} (after {attempts_used} attempts)"
            if last_exc is not None
            else f"failed after {attempts_used} attempts"
        )
        return {
            "accession": accession,
            "status": "failed",
            "location": None,
            "url": failed_url,
            "error": failed_error,
            "retry_count": max(0, attempts_used - 1),
        }

    if pending_rows:
        max_workers = max(1, min(download_workers, len(pending_rows)))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            future_to_index = {
                pool.submit(_download_with_retries, row): row_idx for row_idx, row in pending_rows
            }
            for future in as_completed(future_to_index):
                row_idx = future_to_index[future]
                try:
                    result = future.result()
                except Exception as exc:
                    row = pending_rows_by_index[row_idx]
                    accession = _normalize_accession(str(row["accession"]))
                    result = {
                        "accession": accession,
                        "status": "failed",
                        "location": None,
                        "url": row.get("download_url"),
                        "error": f"unexpected worker failure: {exc}",
                        "retry_count": 0,
                    }
                accession = str(result["accession"])
                retry_count = int(result.get("retry_count", 0))
                if progress_callback is not None:
                    for _ in range(retry_count):
                        progress_callback(
                            "retry",
                            completed_count,
                            max(selected_total - completed_count, 0),
                            selected_total,
                            accession,
                        )

                status = str(result["status"])
                if status == "downloaded":
                    location = str(result["location"])
                    url = None if result.get("url") is None else str(result["url"])
                    local_db.set_location(accession, location, download_url=url)
                    records_by_index[row_idx] = {
                        "accession": accession,
                        "status": "downloaded",
                        "location": location,
                        "url": url,
                        "error": None,
                    }
                    completed_count += 1
                    if progress_callback is not None:
                        progress_callback(
                            "downloaded",
                            completed_count,
                            max(selected_total - completed_count, 0),
                            selected_total,
                            accession,
                        )
                    continue

                records_by_index[row_idx] = {
                    "accession": accession,
                    "status": "failed",
                    "location": None,
                    "url": None if result.get("url") is None else str(result["url"]),
                    "error": None if result.get("error") is None else str(result["error"]),
                }
                completed_count += 1
                if progress_callback is not None:
                    progress_callback(
                        "failed",
                        completed_count,
                        max(selected_total - completed_count, 0),
                        selected_total,
                        accession,
                    )

    if records_by_index:
        records = [records_by_index[idx] for idx in sorted(records_by_index)]

    local_db.sync()
    local_db.save()

    report_schema = {
        "accession": pl.Utf8,
        "status": pl.Utf8,
        "location": pl.Utf8,
        "url": pl.Utf8,
        "error": pl.Utf8,
    }
    if not records:
        return pl.DataFrame(
            {
                "accession": [],
                "status": [],
                "location": [],
                "url": [],
                "error": [],
            },
            schema=report_schema,
        )

    def _as_str_or_none(value):
        return None if value is None else str(value)

    report_columns = {
        "accession": [_as_str_or_none(row.get("accession")) for row in records],
        "status": [_as_str_or_none(row.get("status")) for row in records],
        "location": [_as_str_or_none(row.get("location")) for row in records],
        "url": [_as_str_or_none(row.get("url")) for row in records],
        "error": [_as_str_or_none(row.get("error")) for row in records],
    }
    out_df = pl.DataFrame(report_columns, schema=report_schema)
    if progress_callback is not None:
        progress_callback(
            "done",
            completed_count,
            max(selected_total - completed_count, 0),
            selected_total,
            None,
        )
    return out_df


def load_abundance_table(table_source: str | pathlib.Path | pl.DataFrame | pl.LazyFrame) -> pl.LazyFrame:
    if isinstance(table_source, pl.LazyFrame):
        return table_source
    if isinstance(table_source, pl.DataFrame):
        return table_source.lazy()
    path = pathlib.Path(table_source)
    if not path.exists():
        raise FileNotFoundError(f"Abundance table does not exist: {path}")
    if path.suffix.lower() == ".parquet":
        return pl.scan_parquet(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pl.scan_csv(path, separator="\t")
    return pl.scan_csv(path)


def build_local_genome_db(
    *,
    tool_name: str,
    abundance_table: str | pathlib.Path | pl.DataFrame | pl.LazyFrame,
    db_path: str | pathlib.Path,
    genomes_dir: str | pathlib.Path,
    download: bool = True,
    overwrite: bool = False,
    nonzero_only: bool = True,
    url_resolver: Callable[[str, dict], str] | None = None,
    downloader: Callable[[str, str, pathlib.Path], pathlib.Path] | None = None,
    progress_callback: ProgressCallback | None = None,
    max_download_attempts: int = 3,
    backoff_base_seconds: float = 1.0,
    download_workers: int = 4,
) -> tuple[LocalGenomeDB, pl.DataFrame, pl.DataFrame]:
    """Build/update a local genome DB from an abundance table.

    Returns:
      - LocalGenomeDB object
      - extracted genomes table (accession/genome_name/download_url)
      - download report table
    """
    adapter = get_tool_adapter(tool_name)
    abundance_lf = load_abundance_table(abundance_table)
    table_base_dir = _abundance_table_base_dir(abundance_table)
    if nonzero_only:
        abundance_lf = _filter_nonzero_abundance_rows(abundance_lf)
    genomes = adapter.extract_genomes(abundance_lf)

    local_db = LocalGenomeDB(db_path)
    local_db.add_or_update_genomes(genomes, source_tool=adapter.name)
    stage_genomes_from_source_paths(
        local_db,
        genomes,
        genomes_dir=genomes_dir,
        base_dir=table_base_dir,
    )

    selected = genomes.select("accession").unique()
    selected_total = selected.height
    if selected_total > 0:
        selected_state = (
            selected.join(
                local_db.db.select("accession", "exists"),
                on="accession",
                how="left",
            )
            .with_columns(pl.col("exists").fill_null(False))
        )
        initial_fetched = selected_state.filter(pl.col("exists")).height
    else:
        initial_fetched = 0

    if progress_callback is not None and selected_total == 0:
        progress_callback("start", 0, 0, 0, None)
        progress_callback("done", 0, 0, 0, None)

    report = pl.DataFrame(
        {
            "accession": [],
            "status": [],
            "location": [],
            "url": [],
            "error": [],
        },
        schema={
            "accession": pl.Utf8,
            "status": pl.Utf8,
            "location": pl.Utf8,
            "url": pl.Utf8,
            "error": pl.Utf8,
        },
    )
    if download:
        report = fetch_missing_genomes(
            local_db,
            genomes_dir=genomes_dir,
            url_resolver=url_resolver,
            downloader=downloader,
            overwrite=overwrite,
            only_accessions=genomes["accession"].to_list(),
            progress_callback=progress_callback,
            initial_fetched=initial_fetched,
            total_selected=selected_total,
            max_download_attempts=max_download_attempts,
            backoff_base_seconds=backoff_base_seconds,
            download_workers=download_workers,
        )
    elif progress_callback is not None:
        progress_callback("start", initial_fetched, max(selected_total - initial_fetched, 0), selected_total, None)
        progress_callback("done", initial_fetched, max(selected_total - initial_fetched, 0), selected_total, None)
    return local_db, genomes, report


def _append_fasta_and_stb(
    *,
    fasta_path: pathlib.Path,
    genome_id: str,
    scaffold_prefix: str,
    fasta_handle,
    stb_handle,
    scaffold_name_counts: dict[str, int],
) -> int:
    """Append one genome FASTA into the concatenated FASTA and emit STB rows.

    Returns the number of scaffolds written.
    """
    scaffold_count = 0
    saw_header = False
    with open(fasta_path, "r") as src:
        for raw_line in src:
            if raw_line.startswith(">"):
                header_token = raw_line[1:].strip().split()[0]
                if header_token == "":
                    continue
                saw_header = True
                base_name = f"{scaffold_prefix}__{_sanitize_id(header_token)}"
                seen = scaffold_name_counts.get(base_name, 0)
                scaffold_name_counts[base_name] = seen + 1
                scaffold_name = base_name if seen == 0 else f"{base_name}__{seen + 1}"
                fasta_handle.write(f">{scaffold_name}\n")
                stb_handle.write(f"{scaffold_name}\t{genome_id}\n")
                scaffold_count += 1
                continue

            seq = raw_line.strip()
            if seq == "":
                continue
            if not saw_header:
                raise ValueError(f"Malformed FASTA (sequence before header): {fasta_path}")
            fasta_handle.write(seq + "\n")

    if not saw_header:
        raise ValueError(f"No FASTA headers found in {fasta_path}")
    return scaffold_count


def write_concatenated_reference(
    genomes: pl.DataFrame,
    *,
    output_fasta: str | pathlib.Path,
    output_stb: str | pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Write a concatenated reference FASTA and its STB mapping.

    Expected columns in `genomes`:
      - accession
      - location
    Optional:
      - genome_name
    """
    required = {"accession", "location"}
    missing = required - set(genomes.columns)
    if missing:
        raise ValueError(f"Input genomes table is missing columns: {sorted(missing)}")

    out_fasta = pathlib.Path(output_fasta)
    out_stb = pathlib.Path(output_stb)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_stb.parent.mkdir(parents=True, exist_ok=True)

    sorted_genomes = genomes.sort("accession")
    scaffold_name_counts: dict[str, int] = {}

    with open(out_fasta, "w") as fasta_handle, open(out_stb, "w") as stb_handle:
        for row in sorted_genomes.iter_rows(named=True):
            accession = _normalize_accession(str(row["accession"]))
            location = pathlib.Path(str(row["location"]))
            if not location.exists():
                raise FileNotFoundError(f"Genome file does not exist for {accession}: {location}")

            genome_id = accession
            scaffold_prefix = _sanitize_id(accession)
            _append_fasta_and_stb(
                fasta_path=location,
                genome_id=genome_id,
                scaffold_prefix=scaffold_prefix,
                fasta_handle=fasta_handle,
                stb_handle=stb_handle,
                scaffold_name_counts=scaffold_name_counts,
            )

    return out_fasta, out_stb


def build_reference_from_abundance(
    *,
    tool_name: str,
    abundance_table: str | pathlib.Path | pl.DataFrame | pl.LazyFrame,
    cache_dir: str | pathlib.Path,
    output_dir: str | pathlib.Path,
    url_resolver: Callable[[str, dict], str] | None = None,
    downloader: Callable[[str, str, pathlib.Path], pathlib.Path] | None = None,
    progress_callback: ProgressCallback | None = None,
    max_download_attempts: int = 3,
    backoff_base_seconds: float = 1.0,
    download_workers: int = 4,
    continue_on_missing: bool = True,
) -> tuple[pathlib.Path, pathlib.Path, pl.DataFrame, pl.DataFrame, dict[str, int]]:
    """Build a reference bundle (concatenated FASTA + STB) from an abundance table.

    Genomes are cached under `cache_dir`, so previously downloaded genomes are reused.
    Only genomes with non-zero abundance in at least one numeric column are included.
    """
    cache_path = pathlib.Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    output_path = pathlib.Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    db_path = cache_path / LocalGenomeDB.DEFAULT_DB_NAME
    genomes_cache_dir = cache_path / "genomes"

    local_db, extracted, report = build_local_genome_db(
        tool_name=tool_name,
        abundance_table=abundance_table,
        db_path=db_path,
        genomes_dir=genomes_cache_dir,
        download=True,
        overwrite=False,
        nonzero_only=True,
        url_resolver=url_resolver,
        downloader=downloader,
        progress_callback=progress_callback,
        max_download_attempts=max_download_attempts,
        backoff_base_seconds=backoff_base_seconds,
        download_workers=download_workers,
    )

    selected = extracted.select("accession").unique()
    resolved = (
        selected.join(
            local_db.db.select("accession", "location", "exists"),
            on="accession",
            how="left",
        )
        .filter(pl.col("exists") & pl.col("location").is_not_null())
        .select("accession", "location")
    )

    missing_count = selected.height - resolved.height
    if resolved.height != selected.height:
        missing = selected.join(resolved.select("accession"), on="accession", how="anti")
        missing_ids = missing["accession"].to_list()
        failed_downloads = {}
        if report.height > 0:
            failed_rows = report.filter(pl.col("status") == "failed")
            failed_downloads = {
                row["accession"]: row["error"]
                for row in failed_rows.iter_rows(named=True)
            }
        details = []
        for acc in missing_ids[:10]:
            if acc in failed_downloads:
                details.append(f"{acc} (download failed: {failed_downloads[acc]})")
            else:
                details.append(f"{acc} (not found in local source path or download)")
        detail_msg = "; ".join(details)
        if not continue_on_missing or resolved.height == 0:
            raise RuntimeError(
                "Failed to resolve all requested genomes in cache/download step. "
                f"Missing count={len(missing_ids)}. "
                f"Examples: {detail_msg}"
            )

    if progress_callback is not None and selected.height > 0:
        progress_callback(
            "assembling_reference",
            selected.height,
            0,
            selected.height,
            None,
        )

    out_fasta, out_stb = write_concatenated_reference(
        resolved,
        output_fasta=output_path / "reference_genomes.fna",
        output_stb=output_path / "reference_genomes.stb",
    )
    selected_genomes = extracted.height
    attempted_downloads = report.height
    downloaded_now = report.filter(pl.col("status") == "downloaded").height if report.height else 0
    failed_downloads = report.filter(pl.col("status") == "failed").height if report.height else 0
    cached_before_download = max(0, selected_genomes - attempted_downloads)
    cached_after_download = resolved.height
    summary = {
        "selected_genomes": selected_genomes,
        "cached_before_download": cached_before_download,
        "attempted_downloads": attempted_downloads,
        "downloaded_now": downloaded_now,
        "failed_downloads": failed_downloads,
        "cached_after_download": cached_after_download,
        "missing_after_retries": missing_count,
        "max_download_attempts": max_download_attempts,
        "download_workers": download_workers,
    }
    if progress_callback is not None and selected.height > 0:
        progress_callback(
            "completed",
            selected.height,
            0,
            selected.height,
            None,
        )
    return out_fasta, out_stb, extracted, report, summary


def write_report_csv(report: pl.DataFrame, output_file: str | pathlib.Path) -> None:
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=report.columns)
        writer.writeheader()
        for row in report.iter_rows(named=True):
            writer.writerow(row)
