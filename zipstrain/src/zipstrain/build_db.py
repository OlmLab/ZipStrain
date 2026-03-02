"""zipstrain.build_db
=====================
Helpers to build a local genome database from abundance outputs produced by
external tools (for example, Sylph).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
import csv
import gzip
import pathlib
import re
import shutil
import tempfile
from typing import Callable
from urllib.parse import urlparse
from urllib.request import urlopen
import zipfile

import polars as pl


# Matches common assembly accessions, with optional version.
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d{6,}(?:\.\d+)?)", flags=re.IGNORECASE)

STRING_DTYPES = {pl.Utf8, pl.Categorical, pl.Enum}


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
        return _extract_accessions_generic(
            abundance_table,
            preferred_columns=self.PREFERRED_COLUMNS,
            url_columns=self.URL_COLUMNS,
        )


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
    return SylphAdapter().extract_genomes(sylph_output).lazy()


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


def fetch_missing_genomes(
    local_db: LocalGenomeDB,
    genomes_dir: str | pathlib.Path,
    *,
    url_resolver: Callable[[str, dict], str] | None = None,
    downloader: Callable[[str, str, pathlib.Path], pathlib.Path] | None = None,
    overwrite: bool = False,
) -> pl.DataFrame:
    """Download missing genomes and update DB locations/existence.

    Returns a report with columns:
      accession, status, location, url, error
    """
    resolver = url_resolver or _default_url_resolver
    dl = downloader or _download_single_genome
    genomes_path = pathlib.Path(genomes_dir)

    records: list[dict[str, str | None]] = []
    for row in local_db.missing().iter_rows(named=True):
        accession = row["accession"]
        location = row["location"]
        if location and pathlib.Path(location).exists() and not overwrite:
            records.append(
                {
                    "accession": accession,
                    "status": "already_present",
                    "location": location,
                    "url": row.get("download_url"),
                    "error": None,
                }
            )
            continue

        try:
            url = resolver(accession, row)
            out_file = dl(url, accession, genomes_path)
            local_db.set_location(accession, out_file, download_url=url)
            records.append(
                {
                    "accession": accession,
                    "status": "downloaded",
                    "location": str(out_file),
                    "url": url,
                    "error": None,
                }
            )
        except Exception as exc:
            records.append(
                {
                    "accession": accession,
                    "status": "failed",
                    "location": None,
                    "url": row.get("download_url"),
                    "error": str(exc),
                }
            )

    local_db.sync()
    local_db.save()
    if not records:
        return pl.DataFrame(
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
    return pl.DataFrame(records)


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
    url_resolver: Callable[[str, dict], str] | None = None,
    downloader: Callable[[str, str, pathlib.Path], pathlib.Path] | None = None,
) -> tuple[LocalGenomeDB, pl.DataFrame, pl.DataFrame]:
    """Build/update a local genome DB from an abundance table.

    Returns:
      - LocalGenomeDB object
      - extracted genomes table (accession/genome_name/download_url)
      - download report table
    """
    adapter = get_tool_adapter(tool_name)
    abundance_lf = load_abundance_table(abundance_table)
    genomes = adapter.extract_genomes(abundance_lf)

    local_db = LocalGenomeDB(db_path)
    local_db.add_or_update_genomes(genomes, source_tool=adapter.name)

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
        )
    return local_db, genomes, report


def write_report_csv(report: pl.DataFrame, output_file: str | pathlib.Path) -> None:
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=report.columns)
        writer.writeheader()
        for row in report.iter_rows(named=True):
            writer.writerow(row)
