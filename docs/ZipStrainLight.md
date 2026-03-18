# ZipStrain-Light

`zipstrain-light` is a separate command-line interface for a lighter genome-comparison path that uses compact profile artifacts.

It is designed to keep the original `zipstrain` CLI unchanged while adding a new workflow optimized for comparing profiles that share the same backbone reference.

## Design

Each light profile is a directory with one of two formats:

- `--engine duckdb`:
  - `profile.duckdb` with `coverage`, `snp`, `genomes`, `metadata` tables
- `--engine polars`:
  - `coverage.parquet`
  - `snp.parquet`

Comparison logic:

- Shared-position denominator comes from inner-joined `coverage` rows passing `min_cov` in both samples.
- SNP differences come from the union of `snp` rows at shared covered loci.
- Final metrics are computed with the same downstream aggregation logic as regular ZipStrain genome comparison.

## Commands

## 1) Build Light Profile

From an existing ZipStrain profile parquet:

```bash
zipstrain-light profile \
  --profile-parquet sample_profile.parquet \
  --reference-fasta reference_genomes.fna \
  --engine duckdb \
  --min-cov 5 \
  --output-dir sampleA.light_profile
```

Optional DuckDB runtime controls:

- `--duckdb-memory-limit` (example: `2GB`)
- `--duckdb-temp-directory`
- `--duckdb-threads`
- `--min-cov` (default `5`): keep only loci with `cov > min_cov` in light profile tables.

You can also profile directly from BAM (`--bam-file`) and write light profile artifacts directly into the output directory.

## 2) Compare Two Light Profiles

```bash
zipstrain-light compare \
  --profile-1 sampleA.light_profile \
  --profile-2 sampleB.light_profile \
  --stb-file reference_genomes.stb \
  --output-file sampleA_sampleB_comparison.parquet
```

Common options:

- `--stb-file` (required; keeps compare behavior aligned with regular ZipStrain)
- `--genome all|<genome_id>`
- `--min-cov`
- `--min-gene-compare-len`
- `--engine duckdb|polars` (default: `duckdb`)
- `--calculate ani|ibs|genes|all|<combo>` (default: `all`, for example `ani+ibs`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`
- `--sample-1`, `--sample-2`

## Output

The compare output parquet always includes `genome` plus selected metric columns:

- `genome`
- `total_positions`
- `share_allele_pos`
- `genome_pop_ani`
- `max_consecutive_length`
- `shared_genes_count`
- `identical_gene_count`
- `perc_id_genes`
- `sample_1`
- `sample_2`

## Notes

- Current light compare supports `popani`.
- The regular `zipstrain` command remains unchanged.
