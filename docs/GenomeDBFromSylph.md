# Build Reference FASTA/STB from Sylph Abundances

This guide shows how to build a reference bundle directly from a Sylph abundance table using:

```bash
zipstrain utilities build-genome-db
```

## What this does

Given a Sylph abundance table, ZipStrain will:

1. Keep genomes with non-zero abundance in at least one sample.
2. Resolve/download those genomes into a persistent cache directory.
3. Reuse genomes already present in that cache (no redownload).
4. Write a concatenated reference FASTA and STB file in your output directory.

## Command

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table /path/to/sylph_abundance.csv \
  --cache-dir /path/to/genome_cache \
  --output-dir /path/to/reference_bundle \
  --download-retries 3 \
  --retry-backoff-seconds 1.0 \
  --download-workers 4
```

## Inputs

- `--tool`: abundance parser to use (currently `sylph`).
- `--abundance-table`: `.csv`, `.tsv`, or `.parquet` table.
- `--cache-dir`: persistent genome cache directory.
- `--output-dir`: output directory for the final reference bundle.

For Sylph tables, ZipStrain extracts genome accessions from the `Genome_file` column
(case-insensitive), including GTDB-style paths such as:
`gtdb_genomes_reps_r220/database/GCA/.../GCA_949068525.1_genomic.fna.gz`.

If `Genome_file` points to a local file (absolute path or path relative to the abundance-table directory),
ZipStrain loads it directly into cache first, then only downloads what is still missing.

## Outputs

The command writes:

- `/path/to/reference_bundle/reference_genomes.fna`
- `/path/to/reference_bundle/reference_genomes.stb`
- `/path/to/reference_bundle/genome_db_build_report.txt`

### STB format

`reference_genomes.stb` has two columns (tab-separated, no header):

- scaffold ID in the concatenated FASTA
- genome ID

Genome IDs are accessions (for example `GCF_000001405.40`).

## Cache behavior

Inside `--cache-dir`, ZipStrain maintains:

- a local DB index (`.genome_db.parquet`)
- downloaded genomes under `genomes/`

Re-running with the same cache directory avoids redownloading genomes that already exist.

## Retry behavior

For genomes that are not available locally/in-cache, ZipStrain retries each download with exponential backoff (default: up to 3 attempts per genome).  
If a genome still fails after retries, it is skipped, and the reference bundle is built from successfully fetched genomes.
Parallelism for remote fetch is controlled with `--download-workers`.

If you see repeated `Too Many Requests` errors on large runs:

- lower `--download-workers` (for example `1` or `2`)
- increase `--download-retries` (for example `8`)
- increase `--retry-backoff-seconds` (for example `3` to `5`)

## Console summary

`build-genome-db` prints a short run summary:

- selected genomes (non-zero abundance)
- genomes already cached before the run
- new download attempts
- downloaded now / failed
- genomes available in cache after the run

The same summary is saved to `genome_db_build_report.txt` and includes explicit failed accession IDs (with error messages) when downloads fail.
