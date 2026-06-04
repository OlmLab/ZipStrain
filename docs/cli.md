# ZipStrain Command Line Interface

This page is the command reference.
If you want end-to-end examples first, use the [Tutorial](./Tutorial.md), which includes:

- a standard workflow using the Python CLI
- a standard workflow using Nextflow
- a matrix workflow for repeated all-vs-all comparison

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
  --bed-file genomes_bed_file.bed \
  --genome-length-file genome_lengths.parquet \
  --run-dir profile_run
```

Options:

- `-i, --input-table` (required)
- `-s, --stb-file` (required)
- `-u, --null-model` (required)
- `-g, --gene-range-table` (optional)
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
  --stb-file mapping.stb \
  --output-dir profiling_assets
```

Options:

- `-r, --reference-fasta` (required)
- `-g, --gene-fasta` (optional)
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
  --num-chunks 24 \
  --max-concurrency 4 \
  --output-dir sample_profile
```

Options:

- `-b, --bed-file` (required)
- `-a, --bam-file` (required)
- `-s, --stb-file` (required)
- `-m, --null-model` (required)
- `-g, --gene-range-table` (optional)
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
| `zipstrain utilities single_compare_genome` | Compare one pair at genome level |
| `zipstrain utilities chunk-genome-compare` | Compare many genome-level pairs in Python-side parallel batches |
| `zipstrain utilities single_compare_gene` | Compare one pair at gene level |
| `zipstrain utilities generate-genome-pairs` | Create all non-redundant standard-profile pairs |
| `zipstrain utilities build-profile-db` | Build a profile DB parquet from `profile_name,profile_location` |
| `zipstrain utilities to-complete-table` | Emit not-yet-completed pair table |

### `zipstrain compare genomes`

```bash
zipstrain compare genomes \
  --profile-db profile_db.parquet \
  --scope all \
  --stb-file mapping.stb \
  --run-dir compare_run \
  --ani-method popani \
  --engine duckdb \
  --calculate all
```

Options:

- `--profile-db` (required)
- `--comp-db-file` (optional current genome comparison parquet)
- `--scope` (default: `all`)
- `--min-cov` (default: `5`)
- `--min-gene-compare-len` (default: `100`)
- `--stb-file` (optional)
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
  --profile-db profile_db.parquet \
  --scope all:all \
  --stb-file mapping.stb \
  --run-dir gene_compare_run
```

Options:

- `--profile-db` (required)
- `--comp-db-file` (optional current gene comparison parquet)
- `--scope` (default: `all:all`)
- `--min-cov` (default: `5`)
- `--min-gene-compare-len` (default: `100`)
- `--stb-file` (optional)
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

### `zipstrain utilities generate-genome-pairs`

```bash
zipstrain utilities generate-genome-pairs \
  --profile-dir profiles \
  --output-file genome_pairs.parquet
```

This writes a parquet table with:

- `sample_name_1`
- `sample_name_2`
- `profile_location_1`
- `profile_location_2`

Options:

- `-p, --profile-dir` (required)
- `-o, --output-file` (required)
- `--write-batch-size` (default: `100000`)

### `zipstrain utilities chunk-genome-compare`

```bash
zipstrain utilities chunk-genome-compare \
  --pair-table genome_pairs.parquet \
  --stb-file mapping.stb \
  --output-file chunk_compare.parquet \
  --workers 8 \
  --engine polars
```

This command runs standard genome comparisons directly inside Python for one pair-table chunk.
It is intended as an experimental utility for benchmarking or ad hoc compare runs, and
does not change the main workflow commands.

Accepted pair-table schemas:

- `sample_name_1`, `sample_name_2`, `profile_location_1`, `profile_location_2`
- `sample_name_1`, `sample_name_2`, `profile_1`, `profile_2`
- `sample_1`, `sample_2`, `profile_1`, `profile_2`
- `profile_location_1`, `profile_location_2`
- `profile_1`, `profile_2`

Options:

- `-p, --pair-table` (required)
- `-s, --stb-file` (required)
- `-o, --output-file` (required)
- `-w, --workers` (defaults to CPU count capped by pair count)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-g, --genome` (default: `all`)
- `-a, --ani-method` (default: `popani`)
- `--calculate` (default: `all`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

The final console summary includes:

- total pairs processed
- total genome-level output rows written
- total elapsed time
- average wall time per pair
- average compute time per pair
- average time per genome-level output row

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

### `zipstrain utilities to-complete-table`

```bash
zipstrain utilities to-complete-table \
  --profile-db profile_db.parquet \
  --comp-db-file current_compare.parquet \
  --output-file remaining_pairs.csv
```

Options:

- `--profile-db` (required)
- `--comp-db-file` (optional)
- `-o, --output-file` (required)

## Utilities

### Utility Commands At A Glance

| Command | Purpose |
|---|---|
| `zipstrain utilities build-null-model` | Build sequencing-error null model |
| `zipstrain utilities merge_parquet` | Merge parquet files |
| `zipstrain utilities merge-stat-tables` | Merge gene/genome stat parquet files with sample labels |
| `zipstrain utilities get-coverage-stats` | Rebuild coverage-only gene/genome stats from a profile parquet |
| `zipstrain utilities process_mpileup` | Convert mpileup stream to parquet |
| `zipstrain utilities make_bed` | Build bed chunks from fasta |
| `zipstrain utilities get_genome_lengths` | Genome lengths from STB + BED |
| `zipstrain utilities generate-genome-pairs` | Create all non-redundant standard-profile pairs |
| `zipstrain utilities chunk-genome-compare` | Compare many genome-level pairs in Python-side parallel batches |
| `zipstrain utilities strain_heterogeneity` | Strain heterogeneity metrics |
| `zipstrain utilities build-profile-db` | Build profile DB parquet |
| `zipstrain utilities build-matrix-db` | Build the current per-sample genome matrix store directly from profile parquets |
| `zipstrain utilities append-matrix-db` | Append new profiles into an existing matrix store |
| `zipstrain utilities matrix-db-to-hdf5` | Convert a legacy matrix DuckDB into the current matrix-store format |
| `zipstrain utilities matrix-compare` | Resumable all-vs-all matrix compare into a DuckDB compare DB |
| `zipstrain utilities matrix-compare-export` | Export a matrix compare DuckDB to parquet |
| `zipstrain utilities build-genome-db` | Build local genome reference bundle from abundance table |
| `zipstrain utilities presence-profile` | Presence profile from coverage + read locations |
| `zipstrain utilities process-read-locs` | Process read-location stream |
| `zipstrain utilities generate_stb` | Create scaffold-to-genome map from genome files |
| `zipstrain utilities gene-range-table` | Create gene range table |
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

- `--download-retries` (default: `8`)
- `--retry-backoff-seconds` (default: `10.0`)
- `--download-workers` (default: `1`)

### `zipstrain utilities build-matrix-db`

```bash
zipstrain utilities build-matrix-db \
  --profile-dir profiles \
  --output-file matrix_db.h5 \
  --bed-file genomes_bed_file.bed \
  --memory-limit-gb 16
```

What it does:

- scans a directory of standard ZipStrain profile parquets
- builds one matrix store directly from those profiles
- stores each genome as one sample-major dense dataset with shape `samples x positions x 4`
- positions with total coverage below `5` are zeroed during matrix build
- can optionally store scaffold-relative gene ranges for later gene ANI
- is intended for repeated cohort-scale comparison runs against the same reference set

Important options:

- `-p, --profile-dir` (required)
- `-o, --output-file` (required)
- `-g, --genome` optional genome scope (default: `all`)
- `-b, --bed-file` optional BED file to define scaffold extents instead of scanning profile min/max positions
- `--count-dtype` stored matrix dtype (`uint16|uint32`, default: `uint16`)
- `--memory-limit-gb` approximate maximum memory budget for the entire build process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)

Notes:

- the output matrix store is intended for `zipstrain utilities matrix-compare`
- new matrix stores are append-friendly on the sample axis
- install matrix support with `pip install "zipstrain[matrix]"`
- the CLI shows a progress bar in an interactive terminal
- in non-interactive runs, the CLI emits throttled structured progress lines to stderr for log files
- if `--bed-file` is provided, scaffold spans come from grouped BED intervals rather than inferred profile min/max positions

### `zipstrain utilities append-matrix-db`

```bash
zipstrain utilities append-matrix-db \
  --profile-dir new_profiles \
  --matrix-db-file matrix_db.h5 \
  --memory-limit-gb 16
```

What it does:

- scans a directory of new standard ZipStrain profile parquets
- validates that they match the existing matrix-store contract
- appends new sample rows and whole-genome matrices into the existing matrix store

Important options:

- `-p, --profile-dir` (required)
- `-m, --matrix-db-file` (required)
- `--memory-limit-gb` approximate maximum memory budget for the append process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB used when rewriting an older fixed-size store (default: `128.0`)

Append requirements:

- sample names must be new
- genomes and scaffolds must already exist in the matrix store
- profile positions must stay within the stored scaffold coordinate ranges

### `zipstrain utilities matrix-db-to-hdf5`

```bash
zipstrain utilities matrix-db-to-hdf5 \
  --matrix-db-file matrix_db.duckdb \
  --output-file matrix_db.h5
```

What it does:

- converts an existing legacy matrix DuckDB into the current matrix-store layout
- preserves sample, genome, and scaffold metadata

Important options:

- `-m, --matrix-db-file` (required)
- `-o, --output-file` optional output HDF5 path; defaults to the same basename with `.h5`
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)

### `zipstrain utilities matrix-compare`

```bash
zipstrain utilities matrix-compare \
  --matrix-db-file matrix_db.h5 \
  --output-file matrix_compare.duckdb \
  --memory-limit-gb 16 \
  --anchor-queue-size 1 \
  --target-queue-size 1 \
  --result-transfer-batch-size 1 \
  --loader-executor thread \
  --writer-executor thread \
  --calculate ani+ibs \
  --backend numpy
```

What it does:

- reads a per-sample genome-matrix store from `build-matrix-db` or `matrix-db-to-hdf5`
- writes results into a DuckDB compare database
- if the compare DB already exists, only pairs not yet marked completed are processed
- loads one anchor sample plus as many target samples as fit the memory budget
- computes genome ANI from dense whole-genome matrices
- computes IBS from the shared-allele boolean mask
- computes gene ANI when gene annotations are present in the matrix store and `gene` is requested
- stores result rows and completion metadata in the compare DB incrementally

Important options:

- `-m, --matrix-db-file` (required)
- `-o, --output-file` (required)
- `-g, --genome` optional genome scope (default: `all`)
- `--memory-limit-gb` approximate compare memory budget
- `--anchor-queue-size` number of torch anchor matrices to keep queued in host RAM while still transferring only one anchor at a time to the GPU (default: `1`)
- `--target-queue-size` number of torch target blocks to keep queued in host RAM; `1` preserves the current synchronous target-load behavior (default: `1`)
- `--result-transfer-batch-size` number of torch compare units to batch before transferring result vectors back to CPU (default: `1`)
- `--loader-executor` executor kind for torch loader prefetch work (`thread|process`, default: `thread`)
- `--writer-executor` executor kind for torch result writing/checkpoint work (`thread|process`, default: `thread`)
- `--calculate` matrix metrics to compute:
  - `ani`
  - `ani+ibs`
  - `+gene` or `gene` for gene ANI
  - `all` means `ani+ibs`, and also `gene` when the matrix store contains gene annotations
- `--backend` compute backend:
  - `numpy`
  - `torch`
  - `torch-cpu`
  - `torch-cuda`
  - `torch-mps`

Notes:

- install Torch support with `pip install "zipstrain[matrix]"`
- on Apple Silicon, the standard `torch` wheel can use MPS
- on Linux with NVIDIA GPUs, replace Torch with the CUDA wheel that matches your system, for example:

  ```
  pip install "zipstrain[matrix]"
  pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
  ```

- MPS requires native macOS; Linux containers cannot expose Apple Metal
- `torch` auto-selects CUDA, then MPS, then CPU
- for torch backends, GPU work stays on the main process while loader and writer stages can run through either thread or process executors
- the compare database is resumable; rerunning the same command on the same output file only processes unfinished sample pairs
- the CLI shows a progress bar in an interactive terminal
- in non-interactive runs, the CLI emits throttled structured progress lines to stderr for log files

### `zipstrain utilities matrix-compare-export`

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare.parquet
```

Export the gene table instead:

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare_gene.parquet \
  --table gene
```

What it does:

- reads `matrix_compare_results` from a matrix compare DuckDB
- exports the standard compare columns to parquet
- uses the stored compare metadata to choose the correct output columns
- can also export `matrix_compare_gene_results` with `--table gene`
- raises an error if `--table gene` is requested but the compare DB has no gene table

Important options:

- `-m, --matrix-compare-db-file` (required)
- `-o, --output-file` (required)
- `--table genome|gene` (optional, default `genome`)

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
- When input files contain ZipStrain compare metadata, those metadata fields must match across inputs unless `--allow-mismatch` is used.
- With `--allow-mismatch`, mismatched compare metadata are rewritten to `NA` in the merged parquet metadata.

### `zipstrain utilities build-profile-db`

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file profile_db.parquet
```

Input CSV columns:

- `profile_name`
- `profile_location`

Notes:

- By default, ZipStrain checks that all listed profiles carry matching embedded contract metadata.
- Use `--allow-mismatch` to skip that validation and build a mixed profile DB intentionally.

### `zipstrain utilities get-coverage-stats`

```bash
zipstrain utilities get-coverage-stats \
  --profile-parquet sample_profile.parquet \
  --gene-bed reference_genomes_gene_ranges.tsv \
  --genome-bed genomes_bed_file.bed \
  --output-dir stats \
  --prefix sample1
```

What it does:

- rebuilds coverage-only genome and gene stats from an existing profile parquet
- writes:
  - `<output-dir>/<prefix>_gene_stats.parquet`
  - `<output-dir>/<prefix>_genome_stats.parquet`
- does not require read-location files
- uses the profile’s existing `gene` and `genome` columns for counts
- uses the supplied gene/genome BED files only to calculate lengths

Output columns:

- gene stats:
  - `genome`
  - `gene`
  - `length`
  - `breadth`
  - `coverage`
  - `5x_cov_sites`
  - `ber`
- genome stats:
  - `genome`
  - `length`
  - `breadth`
  - `coverage`
  - `5x_cov_sites`
  - `ber`

Important options:

- `-p, --profile-parquet` (required)
- `-g, --gene-bed` (required) — supports 4 columns `gene, scaffold, start, end` or 5 columns with genome appended
- `-b, --genome-bed` (required) — supports 3 columns `scaffold, start, end` or 4 columns with genome appended
- `-o, --output-dir` (required)
- `--prefix` (required)

### Other Utility Commands

Use `--help` on each command for full option details:

```bash
zipstrain utilities build-null-model --help
zipstrain utilities merge_parquet --help
zipstrain utilities get-coverage-stats --help
zipstrain utilities process_mpileup --help
zipstrain utilities make_bed --help
zipstrain utilities get_genome_lengths --help
zipstrain utilities generate-genome-pairs --help
zipstrain utilities chunk-genome-compare --help
zipstrain utilities merge-stat-tables --help
zipstrain utilities strain_heterogeneity --help
zipstrain utilities build-profile-db --help
zipstrain utilities build-matrix-db --help
zipstrain utilities matrix-compare --help
zipstrain utilities presence-profile --help
zipstrain utilities process-read-locs --help
zipstrain utilities generate_stb --help
zipstrain utilities gene-range-table --help
zipstrain test --help
```
