# ZipStrain-Light

`zipstrain-light` is a separate command-line interface for a lighter genome-comparison path that leaves regular `zipstrain` unchanged.

## Design

Each light profile is a single DuckDB file, for example:

- `sample.light_profile.duckdb`

Required table:

- `mpileup`
  - columns: `chrom`, `pos`, `gene`, `genome`, `bit_representation`
  - rows: loci that survive the corrected-coverage filter used at profile time

Optional tables:

- `gene_stats`
- `genome_stats`

`bit_representation` is the allele mask derived from the corrected `A/T/C/G` counts:

- `A = 1`
- `T = 2`
- `C = 4`
- `G = 8`

Comparison logic:

- inner join the two `mpileup` tables on shared `genome,chrom,pos`
- compute allele sharing with a bitwise `AND` on the two masks
- aggregate the shared loci into ANI / IBS / gene-level metrics

Because the coverage filter is applied during profiling, light compare does not apply a second `min_cov` filter.

## Build Light Profile

From an existing ZipStrain profile parquet:

```bash
zipstrain-light profile \
  --profile-parquet sample_profile.parquet \
  --min-cov 5 \
  --output-file sample.light_profile.duckdb
```

Directly from a BAM:

```bash
zipstrain-light profile \
  --bam-file sample.bam \
  --bed-file genomes_bed_file.bed \
  --stb-file reference_genomes.stb \
  --null-model null_model.parquet \
  --gene-range-table gene_range_table.tsv \
  --min-cov 5 \
  --output-file sample.light_profile.duckdb
```

Optional profile flags:

- `--gene-stats/--no-gene-stats`: write `gene_stats` into the DuckDB output file
- `--genome-stats/--no-genome-stats`: write `genome_stats` into the DuckDB output file
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`
- `--num-chunks`
- `--max-concurrency`

`gene_stats` and `genome_stats` are enabled by default.

## Compare Two Light Profiles

```bash
zipstrain-light compare \
  --profile-1 sampleA.light_profile.duckdb \
  --profile-2 sampleB.light_profile.duckdb \
  --calculate ani \
  --output-file sampleA_sampleB_comparison.parquet
```

Common options:

- `--genome all|<genome_id>`
- `--min-gene-compare-len`
- `--calculate ani|ibs|genes|all|<combo>`
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`
- `--sample-1`, `--sample-2`

`--stb-file` is still accepted for backward compatibility, but the current light profile schema already stores `genome` in `mpileup`, so light compare does not require it.

## Output

The compare output parquet always includes `genome` plus the selected metric columns:

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

- `zipstrain-light` currently supports `popani`.
- The regular `zipstrain` command and its profile format are unchanged.
