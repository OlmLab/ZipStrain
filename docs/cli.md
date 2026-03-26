# ZipStrain Command Line Interface

This page is organized by workflow area for easier navigation:

- [Profile](#profile)
- [Comparison](#comparison)
- [Utilities](#utilities)

General usage:

```bash
zipstrain --help
```

For command-specific help:

```bash
zipstrain <command-or-group> --help
zipstrain <group> <command> --help
```

## Profile

### Profile Commands At A Glance

| Command | Purpose |
|---|---|
| `zipstrain profile` | Batch profiling for multiple BAM files |
| `zipstrain utilities prepare_profiling` | Build profiling assets (BED, gene ranges, genome lengths) |
| `zipstrain utilities profile-single` | Profile one BAM file |

### `zipstrain profile`

Run BAM profiling in batch mode.

```bash
zipstrain profile \
  --input-table samples.csv \
  --stb-file mapping.stb \
  --null-model null_model.parquet \
  --gene-range-table gene_range_table.tsv \
  --bed-file genomes_bed_file.bed \
  --genome-length-file genome_lengths.parquet \
  --run-dir profile_run
```

Options:

- `-i, --input-table` (required)
- `-s, --stb-file` (required)
- `-u, --null-model` (required)
- `-g, --gene-range-table` (required)
- `-b, --bed-file` (required)
- `-l, --genome-length-file` (required)
- `-r, --run-dir` (required)
- `-n, --num-procs` (default: `8`)
- `-m, --max-concurrent-batches` (default: `5`)
- `-p, --poll-interval` (default: `1`)
- `-e, --execution-mode` (default: `local`)
- `-c, --slurm-config`
- `-o, --container-engine` (default: `local`)
- `-t, --task-per-batch` (default: `10`)

### `zipstrain utilities prepare_profiling`

Prepare profiling database assets.

```bash
zipstrain utilities prepare_profiling \
  --reference-fasta reference.fasta \
  --gene-fasta genes.fasta \
  --stb-file mapping.stb \
  --output-dir profiling_assets
```

Options:

- `-r, --reference-fasta` (required)
- `-g, --gene-fasta` (required)
- `-s, --stb-file` (required)
- `-o, --output-dir` (required)

### `zipstrain utilities profile-single`

Profile a single BAM.

```bash
zipstrain utilities profile-single \
  --bed-file genomes_bed_file.bed \
  --bam-file sample.bam \
  --stb-file mapping.stb \
  --null-model null_model.parquet \
  --gene-range-table gene_range_table.tsv \
  --genome-length-file genome_lengths.parquet \
  --num-chunks 24 \
  --max-concurrency 4 \
  --output-dir sample_profile
```

Options:

- `-b, --bed-file` (required)
- `-a, --bam-file` (required)
- `-s, --stb-file` (required)
- `-m, --null-model` (required)
- `-g, --gene-range-table` (required)
- `-l, --genome-length-file` (required)
- `-n, --num-chunks` (default: `24`) — number of BED chunks to create
- `-c, --max-concurrency` (default: `4`) — how many chunks run simultaneously
- `-o, --output-dir` (required)

Outputs include:

- `<sample>_profile.parquet`
- `<sample>_genome_stats.parquet`
- `<sample>_gene_stats.parquet`

## Comparison

### Comparison Commands At A Glance

| Command | Purpose |
|---|---|
| `zipstrain compare genomes` | Batch genome-level comparisons |
| `zipstrain compare genes` | Batch gene-level comparisons |
| `zipstrain compare build-comp-database` | Build comparison DB object from profile DB + config |
| `zipstrain utilities single_compare_genome` | Compare one pair at genome level |
| `zipstrain utilities single_compare_gene` | Compare one pair at gene level |
| `zipstrain utilities build-genome-comparison-config` | Build genome comparison config |
| `zipstrain utilities build-gene-comparison-config` | Build gene comparison config |
| `zipstrain utilities to-complete-table` | Emit not-yet-completed pair table |

### `zipstrain compare genomes`

```bash
zipstrain compare genomes \
  --genome-comparison-object genome_comp.json \
  --run-dir compare_run \
  --ani-method popani \
  --engine duckdb \
  --calculate all
```

Options:

- `-g, --genome-comparison-object` (required)
- `-r, --run-dir` (required)
- `-m, --max-concurrent-batches` (default: `5`)
- `-p, --poll-interval` (default: `1`)
- `-e, --execution-mode` (default: `local`)
- `-s, --slurm-config`
- `-c, --container-engine` (default: `local`)
- `-t, --task-per-batch` (default: `10`)
- `-a, --ani-method` (default: `popani`) — ANI method (`popani`, `conani`, `cosani_<threshold>`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--calculate` (`ani`, `ibs`, `identical_genes`, `all`, or `+` combinations like `ani+ibs`, default: `all`)
- `-d, --duckdb-memory-limit`
- `--duckdb-threads`

### `zipstrain compare genes`

```bash
zipstrain compare genes \
  --gene-comparison-object gene_comp.json \
  --run-dir gene_compare_run
```

Options:

- `-g, --gene-comparison-object` (required)
- `-r, --run-dir` (required)
- `-m, --max-concurrent-batches` (default: `5`)
- `-p, --poll-interval` (default: `1`)
- `-e, --execution-mode` (default: `local`)
- `-s, --slurm-config`
- `-c, --container-engine` (default: `local`)
- `-t, --task-per-batch` (default: `10`)
- `-n, --ani-method` (default: `popani`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `-d, --duckdb-memory-limit`
- `--duckdb-threads`

### `zipstrain compare build-comp-database`

```bash
zipstrain compare build-comp-database \
  --profile-db-dir profiles.parquet \
  --config-file comparison_config.json \
  --output-dir comparison_db
```

Options:

- `-p, --profile-db-dir` (required)
- `-c, --config-file` (required)
- `-o, --output-dir` (required)
- `-f, --comp-db-file`

### `zipstrain utilities single_compare_genome`

```bash
zipstrain utilities single_compare_genome \
  --mpileup-contig-1 sample_a.parquet \
  --mpileup-contig-2 sample_b.parquet \
  --stb-file mapping.stb \
  --output-file out.parquet
```

Options:

- `-m1, --mpileup-contig-1` (required)
- `-m2, --mpileup-contig-2` (required)
- `-s, --stb-file` (required)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-o, --output-file` (required)
- `-g, --genome` (default: `all`)
- `-a, --ani-method` (default: `popani`)
- `--calculate` (default: `all`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

### `zipstrain utilities single_compare_gene`

```bash
zipstrain utilities single_compare_gene \
  --mpileup-contig-1 sample_a.parquet \
  --mpileup-contig-2 sample_b.parquet \
  --stb-file mapping.stb \
  --scope all:all \
  --output-file out.parquet
```

Options:

- `-m1, --mpileup-contig-1` (required)
- `-m2, --mpileup-contig-2` (required)
- `-s, --stb-file` (required)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-o, --output-file` (required)
- `-g, --scope` (default: `all:all`)
- `-a, --ani-method` (default: `popani`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

### Comparison Config Helpers

`build-genome-comparison-config` and `build-gene-comparison-config` share the same option pattern:

- `-p, --profile-db` (required)
- `-g, --gene-db-id` (required)
- `-r, --reference-genome-id` (required)
- `-s, --scope` (default: `all` for genome, `all:all` for gene)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `200`)
- `-t, --stb-file-loc` (required)
- `-a, --current-comp-table`
- `-o, --output-file` (required)

### `zipstrain utilities to-complete-table`

```bash
zipstrain utilities to-complete-table \
  --genome-comparison-object genome_comp.json \
  --output-file remaining_pairs.csv
```

## Utilities

### Utility Commands At A Glance

| Command | Purpose |
|---|---|
| `zipstrain utilities build-null-model` | Build sequencing-error null model |
| `zipstrain utilities merge_parquet` | Merge parquet files |
| `zipstrain utilities process_mpileup` | Convert mpileup stream to parquet |
| `zipstrain utilities make_bed` | Build bed chunks from fasta |
| `zipstrain utilities get_genome_lengths` | Genome lengths from STB + BED |
| `zipstrain utilities genome_breadth_matrix` | Per-genome breadth output |
| `zipstrain utilities collect_breadth_tables` | Merge breadth tables |
| `zipstrain utilities strain_heterogeneity` | Strain heterogeneity metrics |
| `zipstrain utilities build-profile-db` | Build profile DB parquet |
| `zipstrain utilities build-genome-db` | Build local genome reference bundle from abundance table |
| `zipstrain utilities presence-profile` | Presence profile from coverage + read locations |
| `zipstrain utilities process-read-locs` | Process read-location stream |
| `zipstrain utilities generate_stb` | Create scaffold-to-genome map from genome files |
| `zipstrain utilities gene-range-table` | Create gene range table |
| `zipstrain utilities gene-loc-table` | Create gene-location table for scaffold list |
| `zipstrain test` | Validate local installation/dependencies |

### `zipstrain utilities build-genome-db`

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table sylph_abundance.tsv \
  --cache-dir genome_cache \
  --output-dir .
```

Important options:

- `--download-retries` (default: `3`)
- `--retry-backoff-seconds` (default: `1.0`)
- `--download-workers` (default: `4`)

### `zipstrain utilities merge_parquet`

```bash
zipstrain utilities merge_parquet \
  --input-dir comps \
  --output-file merged_comparisons.parquet \
  --batch-size 5000
```

Notes:

- `--batch-size -1` keeps the current single-pass merge behavior.
- Positive `--batch-size` values first merge input files batch-by-batch into a temporary directory, then do one final lazy merge over those batch outputs.
- Progress is logged as active line-oriented batch updates and flushed immediately, which is easier to follow in cluster logs.

### Other Utility Commands

Use `--help` on each command for full option details:

```bash
zipstrain utilities build-null-model --help
zipstrain utilities merge_parquet --help
zipstrain utilities process_mpileup --help
zipstrain utilities make_bed --help
zipstrain utilities get_genome_lengths --help
zipstrain utilities genome_breadth_matrix --help
zipstrain utilities collect_breadth_tables --help
zipstrain utilities strain_heterogeneity --help
zipstrain utilities build-profile-db --help
zipstrain utilities presence-profile --help
zipstrain utilities process-read-locs --help
zipstrain utilities generate_stb --help
zipstrain utilities gene-range-table --help
zipstrain utilities gene-loc-table --help
zipstrain test --help
```
