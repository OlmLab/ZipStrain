# Build a Genome Database from Sylph Abundances

This guide shows how to build and maintain a local genome database from a Sylph abundance table using:

```bash
zipstrain utilities build-genome-db
```

## What this does

Given a Sylph abundance table, ZipStrain will:

1. Parse genome accessions from the table (for example `GCF_000001405.40`).
2. Store/update them in a local DB parquet file.
3. Optionally download missing reference genomes into a local directory.

## Inputs

- A Sylph abundance table (`.csv`, `.tsv`, or `.parquet`).
- The table should contain genome identifiers that include assembly accessions (for example values like `.../GCF_000001405.40_genomic.fna.gz`).

## Quick start

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table /path/to/sylph_abundance.csv \
  --db-file /path/to/.genome_db.parquet \
  --genomes-dir /path/to/genomes \
  --download \
  --report-file /path/to/genome_download_report.csv
```

## Dry run (no downloads)

Use this to only parse/update the DB index:

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table /path/to/sylph_abundance.csv \
  --db-file /path/to/.genome_db.parquet \
  --genomes-dir /path/to/genomes \
  --no-download
```

## Outputs

### 1) Local genome DB parquet (`--db-file`)

Columns:

- `accession`
- `genome_name`
- `location`
- `download_url`
- `source_tool`
- `exists`

### 2) Genome FASTA files (`--genomes-dir`)

Downloaded genomes are stored as accession-named fasta files (for example `GCF_000001405.40.fna`).

### 3) Optional report CSV (`--report-file`)

Download report columns:

- `accession`
- `status` (`downloaded`, `already_present`, or `failed`)
- `location`
- `url`
- `error`

## Notes

- Current tool support is `sylph`.
- By default, missing genomes are resolved/downloaded using accession-based NCBI Datasets endpoints.
- Re-run the same command on updated Sylph tables to incrementally update the local genome DB.
