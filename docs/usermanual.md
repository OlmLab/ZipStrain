# User Manual

##  An overview of ZipStrain and the data it generates

ZipStrain is a program for microbial metagenomic analysis. The primary use cases for ZipStrain are to accurately determine organism presence / absence in a community, and to perform detailed comparisons between organisms in different samples.

There are two ways of interacting with ZipStrain- NextFlow and the ZipStrain python CLI.

TODO: Make a figure like this inStrain one (or just edit this inStrain one) so people have a vauge idea of what's going on from the bat: https://instrain.readthedocs.io/en/latest/_images/OverviewFigure1_v1.4.png . We can then add this figure lots of places

## ZipStrain Command Line Interface

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

### Profile

**Profile Commands At A Glance**

| Command | Purpose |
|---|---|
| `zipstrain profile` | Batch profiling for multiple BAM files |
| `zipstrain utilities prepare_profiling` | Build profiling assets (BED, gene ranges, genome lengths, null model, profiling contract) |
| `zipstrain utilities profile-single` | Profile one BAM file |
| `zipstrain utilities get-snp-reference` | Generate SNV table |

#### `zipstrain profile`

Run BAM profiling in batch mode.

```bash
zipstrain profile \
  --input-table samples.csv \
  --stb-file mapping.stb \
  --null-model profiling_assets/null_model.parquet \
  --profiling-contract profiling_assets/profiling_contract.json \
  --bed-file profiling_assets/genomes_bed_file.bed \
  --genome-length-file profiling_assets/genome_lengths.parquet \
  --run-dir profile_run
```

Options:

- `-i, --input-table` (required)
- `-s, --stb-file` (required)
- `-u, --null-model` (required)
- `-g, --gene-range-table` (optional)
- `--profiling-contract` (optional)
- `-b, --bed-file` (required)
- `-l, --genome-length-file` (required)
- `-r, --run-dir` (required)
- `-n, --num-procs` (default: `8`)
- `-m, --max-concurrent-batches` (default: `5`)
- `-p, --poll-interval` (default: `1`)
- `-e, --execution-mode` (default: `local`)
- `-c, --slurm-config`
- `-o, --container-engine` (default: `local`)
- `--container-address` (optional) — explicit image/address override for `docker` or `apptainer`; otherwise the current ZipStrain version tag is used
- `-t, --task-per-batch` (default: `10`)
- `--min-mapq` (default: `0`)
- `--min-baseq` (default: `13`)
- `--min-read-ani` (optional) — filters reads before pileup using the BAM `NM` tag and aligned query span
- `--read-inclusion` (`proper-pairs|paired|all-mapped`, default: `all-mapped`)

Read-filter notes:

- `--min-mapq 0` and `--min-baseq 13` match current samtools mpileup defaults
- `--read-inclusion all-mapped` is the least restrictive mode and is the default
- `--min-read-ani` requires BAM alignments with `NM` tags

#### `zipstrain utilities prepare_profiling`

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
- `-e, --error-rate` (default: `0.001`)
- `-m, --max-total-reads` (default: `10000`)
- `-p, --p-threshold` (default: `0.05`)
- `-t, --model-type` (default: `poisson`)
- `-o, --output-dir` (required)

Outputs:

- `genomes_bed_file.bed`
- `gene_range_table.tsv`
- `genome_lengths.parquet`
- `null_model.parquet`
- `profiling_contract.json`

#### `zipstrain utilities profile-single`

Profile a single BAM.

```bash
zipstrain utilities profile-single \
  --bed-file genomes_bed_file.bed \
  --bam-file sample.bam \
  --stb-file mapping.stb \
  --null-model profiling_assets/null_model.parquet \
  --profiling-contract profiling_assets/profiling_contract.json \
  --num-chunks 24 \
  --max-concurrency 4 \
  --output-dir sample_profile
```

Options:

- `-r, --reference-fasta` (optional) — when provided, profiling also records `ref_base_bitmask` and adds `ref_ani` to gene/genome stat tables
- `-b, --bed-file` (required)
- `-a, --bam-file` (required)
- `-s, --stb-file` (required)
- `-m, --null-model` (required)
- `-g, --gene-range-table` (optional)
- `--profiling-contract` (optional)
- `-n, --num-chunks` (default: `24`) — number of BED chunks to create
- `-c, --max-concurrency` (default: `4`) — how many chunks run simultaneously
- `--min-mapq` (default: `0`)
- `--min-baseq` (default: `13`)
- `--min-read-ani` (optional) — filters reads before pileup using the BAM `NM` tag and aligned query span
- `--read-inclusion` (`proper-pairs|paired|all-mapped`, default: `all-mapped`)
- `-o, --output-dir` (required)

Read-inclusion modes:

- `proper-pairs`: keep only mapped read pairs carrying the aligner `PROPER_PAIR` flag
- `paired`: keep mapped paired-end reads even if they are discordant
- `all-mapped`: keep any mapped read, whether paired or single-end

Outputs include:

- `<sample>_profile.parquet`
- `<sample>_genome_stats.parquet`
- `<sample>_gene_stats.parquet`

When `--reference-fasta` is provided during profiling, the profile parquet includes `ref_base_bitmask`.
In the same case, the generated genome and gene stat tables also include a `ref_ani` column.

`ref_ani` is the percentage of covered sites whose observed allele set still contains the reference allele after ZipStrain's sequence-error adjustment.

`ref_base_bitmask` uses this encoding:

- `1` = reference base `A`
- `2` = reference base `C`
- `4` = reference base `G`
- `8` = reference base `T`
- `0` = non-ACGT or unknown reference base

This is a one-hot bitmask, so current profiles are expected to contain only `0`, `1`, `2`, `4`, or `8` in this column.

#### `zipstrain utilities get-snp-reference`

Emit profile-like rows that are SNPs relative to the reference from a classic profile parquet that includes `ref_base_bitmask`.

```bash
zipstrain utilities get-snp-reference \
  --profile-file sample_profile.parquet \
  --min-cov 5 \
  --output-file sample_reference_snps.parquet
```

Options:

- `-p, --profile-file` (required)
- `-c, --min-cov` (default: `5`)
- `-o, --output-file` (required)

The output preserves the input profile-like columns and includes only positions that:

- have coverage `>= min_cov`
- have a known reference base in `ref_base_bitmask`
- do not retain the reference allele after profile sequence-error adjustment

This uses the same reference-sharing logic used to populate `ref_ani` in the gene and genome stat tables.

### Comparison

**Comparison Commands At A Glance**

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

#### `zipstrain compare genomes`

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
- `--container-address` (optional) — explicit image/address override for `docker` or `apptainer`; otherwise the current ZipStrain version tag is used
- `-t, --task-per-batch` (default: `10`)
- `-a, --ani-method` (default: `popani`) — ANI method (`popani`, `conani`, `cosani_<threshold>`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--calculate` (`ani`, `ibs`, `identical_genes`, `all`, or `+` combinations like `ani+ibs`, default: `all`)
- `-d, --duckdb-memory-limit`
- `--duckdb-threads`

#### `zipstrain compare genes`

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
- `--container-address` (optional) — explicit image/address override for `docker` or `apptainer`; otherwise the current ZipStrain version tag is used
- `-t, --task-per-batch` (default: `10`)
- `-n, --ani-method` (default: `popani`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `-d, --duckdb-memory-limit`
- `--duckdb-threads`

#### `zipstrain utilities single_compare_genome`

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

#### `zipstrain utilities generate-genome-pairs`

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

#### `zipstrain utilities chunk-genome-compare`

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

#### `zipstrain utilities single_compare_gene`

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

#### `zipstrain utilities to-complete-table`

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

Output columns:

- `sample_name_1`
- `sample_name_2`
- `profile_location_1`
- `profile_location_2`

Notes:

- this command does not need `--scope`, `--min-cov`, `--min-gene-compare-len`, or `--stb-file`
- it only compares the sample-pair universe implied by the profile DB against the pairs already present in the current genome comparison parquet

### Utilities

**Utility Commands At A Glance**

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
| `zipstrain utilities matrix-db-to-hdf5` | Convert a DuckDB matrix database into the current matrix-store format |
| `zipstrain utilities matrix-compare` | Resumable all-vs-all matrix compare into a DuckDB compare DB |
| `zipstrain utilities matrix-compare-export` | Export a matrix compare DuckDB to parquet |
| `zipstrain utilities build-genome-db` | Build local genome reference bundle from abundance table |
| `zipstrain utilities presence-profile` | Presence profile from coverage + read locations |
| `zipstrain utilities process-read-locs` | Process read-location stream |
| `zipstrain utilities generate_stb` | Create scaffold-to-genome map from genome files |
| `zipstrain utilities gene-range-table` | Create gene range table |
| `zipstrain test` | Validate local installation/dependencies |

#### `zipstrain utilities build-genome-db`

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

#### `zipstrain utilities build-matrix-db`

```bash
zipstrain utilities build-matrix-db \
  --profile-dir profiles \
  --output-file matrix_db.h5 \
  --bed-file genomes_bed_file.bed \
  --stb-file reference.stb \
  --gene-range-table gene_range_table.tsv \
  --memory-limit-gb 16
```

What it does:

- scans a directory of standard ZipStrain profile parquets
- builds one matrix store directly from those profiles
- uses the BED and STB files as the explicit scaffold/genome contract for the store
- stores each genome as one sample-major dense dataset with shape `samples x positions x 4`
- positions with total coverage below `5` are zeroed during matrix build
- can optionally store scaffold-relative gene ranges for later gene ANI
- is intended for repeated cohort-scale comparison runs against the same reference set

Important options:

- `-p, --profile-dir` (required)
- `-o, --output-file` (required)
- `-g, --genome` optional genome scope (default: `all`)
- `-b, --bed-file` (required) BED file defining scaffold coordinate extents for the matrix contract
- `-s, --stb-file` (required) STB file defining scaffold-to-genome membership for the matrix contract
- `--gene-range-table` optional headerless TSV of `gene, scaffold, start, end` for gene ANI support
- `--count-dtype` stored matrix dtype (`uint16|uint32`, default: `uint16`)
- `--memory-limit-gb` approximate maximum memory budget for the entire build process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)
- `--sparse` store genome matrices sparsely in HDF5

Notes:

- the output matrix store is intended for `zipstrain utilities matrix-compare`
- new matrix stores are append-friendly on the sample axis
- every input profile is interpreted against the BED+STB contract you provide here
- install matrix support with `pip install "zipstrain[matrix]"`
- the CLI shows a progress bar in an interactive terminal
- in non-interactive runs, the CLI emits throttled structured progress lines to stderr for log files
- if `--gene-range-table` is omitted, matrix compare can still compute genome ANI and IBS, but not gene ANI
- `--sparse` reduces on-disk HDF5 size, but matrix compare currently materializes sparse storage back into dense arrays when loading for comparison

#### `zipstrain utilities append-matrix-db`

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
- materializes newly encountered genomes when they are still compatible with the stored BED+STB contract
- ignores genomes that fall outside the stored contract and reports the ignored count

Important options:

- `-p, --profile-dir` (required)
- `-m, --matrix-db-file` (required)
- `--memory-limit-gb` approximate maximum memory budget for the append process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB used when rewriting an older fixed-size store (default: `128.0`)

Append behavior:

- sample names must be new
- known scaffolds and coordinate ranges must stay within the stored contract
- compatible genomes can be appended even if no matrix dataset existed for them yet
- genomes outside the stored contract are skipped and counted in the summary output

#### `zipstrain utilities matrix-db-to-hdf5`

```bash
zipstrain utilities matrix-db-to-hdf5 \
  --matrix-db-file matrix_db.duckdb \
  --output-file matrix_db.h5
```

What it does:

- converts an existing DuckDB matrix database into the current matrix-store layout
- preserves sample, genome, and scaffold metadata
- is only needed when you already have a DuckDB-based matrix database from an older workflow

Important options:

- `-m, --matrix-db-file` (required)
- `-o, --output-file` optional output HDF5 path; defaults to the same basename with `.h5`
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)

#### `zipstrain utilities matrix-compare`

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
- `--backend numpy` works without Torch and is the simplest CPU-only path
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

#### `zipstrain utilities matrix-compare-export`

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

#### `zipstrain utilities merge_parquet`

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

#### `zipstrain utilities build-profile-db`

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
- The output parquet stores at least `profile_name` and `profile_location`, plus shared metadata fields derived from the listed profiles.

#### `zipstrain utilities get-coverage-stats`

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

#### Other Utility Commands

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


## Nextflow Implementation

This page reflects the current `zipstrain.nf` workflow in this repository.

TODO: Add a description of why you'd choose the nextflow implementation over the command line implementation

Important scope note:

- the standard profile-based compare workflows are available in Nextflow
- the newer matrix-store workflow is currently CLI-driven and is documented in [Tutorial](./Tutorial.md) and [CLI](./cli.md)
- if you want a worked standard example, the [Tutorial](./Tutorial.md) now includes both a Python/CLI route and a Nextflow route

### General Usage Instructions

TODO: Here's some instructions on how to install Nextflow, what the benefits of nextflow are, and how to set up your config file for this


#### General Running Pattern

```bash
nextflow run zipstrain.nf \
  --mode <mode> \
  --input_table <path/to/input.csv> \
  --output_dir <path/to/output_dir> \
  -c conf.config \
  -profile <docker|alpine|gutbot|blanca|fiji> \
  -resume
```

`conf.config` already defines resources for the current process set and includes example execution profiles.
The `fiji` profile is configured for Slurm plus Singularity, which is useful on clusters that provide Singularity rather than Apptainer.
Review the container paths or tags in your config before running on a new system.

#### ZipStrain Nextflow Pipeline Commands

- Read mapping with Bowtie2 (`map_reads`)
- Profile generation from BAM files (`profile`)
- End-to-end SRA to profile (`from_sra_to_profile`)
- Pairwise genome comparison across profiles (`compare_genomes`)
- Pairwise gene comparison across profiles (`compare_genes`)

#### Key Pipeline Parameters

- `--mode`: `map_reads`, `profile`, `from_sra_to_profile`, `compare_genomes`, `compare_genes`
- `--input_type`: depends on mode (`local`, `sra`, `profile_table`, `pair_table`)
- `--parallel_mode`: `single` or `batched` for comparison workflows
- `--batch_size`: number of pairs per batch when `--parallel_mode batched`
- `--batch_compare_n_parallel`: parallel jobs inside each batched comparison task
- `--compare_genome_scope`: genome scope for genome comparisons (`all` or genome ID)
- `--compare_genome_calculate`: genome metrics to compute (`all` or `ani`)
- `--compare_gene_scope`: gene scope for gene comparisons (`all:all`, `<genome>:all`, `all:<gene>`, `<genome>:<gene>`)
- `--compare_ani_method`: ANI method forwarded to compare tasks (`popani`, `conani`, `cosani_<threshold>`)
- `--compare_engine`: comparison engine for standard compare tasks (`polars` or `duckdb`). Default: `polars`
- `--compare_duckdb_memory_limit`: forwarded to single compare commands
- `--compare_calculate`: genome metrics for genome compare (`ani`, `ibs`, `identical_genes`, `all`, or `+` combinations). Default: `all`

#### Specific Notes

- For `--input_type profile_table`, use `sample_name` and `profile_location`.
- For `--input_type pair_table`, use `sample_name_1`, `sample_name_2`, `profile_location_1`, and `profile_location_2`.
- `--compare_duckdb_memory_limit` is only relevant when `--compare_engine duckdb`.
- The current Nextflow profiling workflow uses the default profiling read filters unless you edit `zipstrain.nf` directly.
- For auto-built references, genome selection comes from the merged Sylph abundance table through `zipstrain utilities build-genome-db`.


### Command 1: Map Reads (`mode=map_reads`)

#### Input Table (`--input_type local`)

Paired-end:

```csv
sample_name,reads1,reads2
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz
```

Single-end:

```csv
sample_name,reads1
S1,/data/S1.fastq.gz
S2,/data/S2.fastq.gz
```

#### Input Table (`--input_type sra`)

```csv
Run
SRR12345678
SRR12345679
```

#### Option 1.A: Use Existing Reference Genome

```bash
nextflow run zipstrain.nf \
  --mode map_reads \
  --input_type local \
  --input_table reads.csv \
  --reference_genome reference_genomes.fna \
  --stb reference_genomes.stb \
  --output_dir out_map \
  -c conf.config \
  -profile docker \
  -resume
```

Optional:

- `--index_files` to reuse existing Bowtie2 index files
- `--bowtie2_non_competitive_mapping true` to pass `-a` to Bowtie2

#### Option 1.B: Build Reference from Sylph Automatically

If `--reference_genome` is not provided, the pipeline does:

1. per-sample `sylph profile`
2. merge all per-sample Sylph abundance tables
3. `zipstrain utilities build-genome-db --tool sylph ...`
4. `prodigal` gene prediction on the generated reference FASTA
5. Bowtie2 indexing
6. mapping

```bash
nextflow run zipstrain.nf \
  --mode map_reads \
  --input_type local \
  --input_table reads.csv \
  --output_dir out_map \
  --genome_db_cache_dir genome_cache \
  --sylph_db /path/to/custom.syldb \
  -c conf.config \
  -profile docker \
  -resume
```

If `--sylph_db` is omitted, `--sylph_db_link` is used for download.

#### Map Outputs

- BAM files: `<output_dir>/*.bam`
- Sylph tables: `<output_dir>/sylph_abundance/`
- Built reference bundle (when auto-built): `<output_dir>/db_from_sylph/`

### Command 2: Generate Profiles from BAM (`mode=profile`)

#### Input Table

```csv
sample_name,bamfile
S1,/data/S1.bam
S2,/data/S2.bam
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode profile \
  --input_table bams.csv \
  --reference_genome reference_genomes.fna \
  --gene_file reference_genomes_gene.fasta \
  --stb reference_genomes.stb \
  --output_dir out_profile \
  -c conf.config \
  -profile docker \
  -resume
```

#### Profile Outputs

- `<output_dir>/*_profile.parquet`
- `<output_dir>/*_genome_stats.parquet`
- `<output_dir>/*_gene_stats.parquet`

Because the current Nextflow profiling modes pass a reference FASTA into profiling, these outputs normally include the reference-aware fields:

- profiles include `ref_base_bitmask`
- gene/genome stat tables include `ref_ani`

### Command 3: End-to-End SRA to Profile (`mode=from_sra_to_profile`)

#### Input Table

```csv
Run
SRR12345678
SRR12345679
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode from_sra_to_profile \
  --input_table sra.csv \
  --reference_genome reference_genomes.fna \
  --gene_file reference_genomes_gene.fasta \
  --stb reference_genomes.stb \
  --output_dir out_sra_profile \
  -c conf.config \
  -profile docker \
  -resume
```

#### Outputs

- `<output_dir>/profiles/*_profile.parquet`
- `<output_dir>/profiles/*_genome_stats.parquet`
- `<output_dir>/profiles/*_gene_stats.parquet`

### Command 4: Compare Genomes (`mode=compare_genomes`)

#### Input Option 4.A: All-vs-All from Profile List (`--input_type profile_table`)

```csv
sample_name,profile_location
S1,/profiles/S1_profile.parquet
S2,/profiles/S2_profile.parquet
S3,/profiles/S3_profile.parquet
```

#### Input Option 4.B: Explicit Pairs (`--input_type pair_table`)

```csv
sample_name_1,sample_name_2,profile_location_1,profile_location_2
S1,S2,/profiles/S1_profile.parquet,/profiles/S2_profile.parquet
S1,S3,/profiles/S1_profile.parquet,/profiles/S3_profile.parquet
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genomes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_engine polars \
  --compare_genome_scope all \
  --compare_calculate ani+ibs+identical_genes \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir out_compare_genomes \
  -c conf.config \
  -profile docker \
  -resume
```

### Command 5: Compare Genes (`mode=compare_genes`)

Input-table formats are the same as genome compare (`profile_table` or `pair_table`).

#### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_engine polars \
  --compare_gene_scope all:all \
  --compare_ani_method popani \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir out_compare_genes \
  -c conf.config \
  -profile docker \
  -resume
```

#### Comparison Outputs

- Final merged table: `<output_dir>/merged_comparisons.parquet`
- Intermediate batched outputs (when `parallel_mode=batched`): `<output_dir>/batch_comparisons/`

## Automatically Build Reference Database

This guide shows how to build a reference bundle directly from a Sylph abundance table using:

```bash
zipstrain utilities build-genome-db
```

Note- This functionality can also be performed at-scale with the helper program [MetaTrawl](https://github.com/OlmLab/MetaTrawl)

### What this does

Given a Sylph abundance table, ZipStrain will:

1. Keep genomes with non-zero abundance in at least one sample.
2. Resolve/download those genomes into a persistent cache directory.
3. Reuse genomes already present in that cache (no redownload).
4. Write a concatenated reference FASTA and STB file in your output directory.

### Command

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table /path/to/sylph_abundance.csv \
  --cache-dir /path/to/genome_cache \
  --output-dir /path/to/reference_bundle \
  --download-retries 8 \
  --retry-backoff-seconds 10.0 \
  --download-workers 1
```

### Inputs

- `--tool`: abundance parser to use (currently `sylph`).
- `--abundance-table`: `.csv`, `.tsv`, or `.parquet` table.
- `--cache-dir`: persistent genome cache directory.
- `--output-dir`: output directory for the final reference bundle.

For Sylph tables, ZipStrain extracts genome accessions from the `Genome_file` column
(case-insensitive), including GTDB-style paths such as:
`gtdb_genomes_reps_r220/database/GCA/.../GCA_949068525.1_genomic.fna.gz`.

If `Genome_file` points to a local file (absolute path or path relative to the abundance-table directory),
ZipStrain loads it directly into cache first, then only downloads what is still missing.

### Outputs

The command writes:

- `/path/to/reference_bundle/reference_genomes.fna`
- `/path/to/reference_bundle/reference_genomes.stb`
- `/path/to/reference_bundle/genome_db_build_report.txt`

#### STB format

`reference_genomes.stb` has two columns (tab-separated, no header):

- scaffold ID in the concatenated FASTA
- genome ID

Genome IDs are accessions (for example `GCF_000001405.40`).

### Cache behavior

Inside `--cache-dir`, ZipStrain maintains:

- a local DB index (`.genome_db.parquet`)
- downloaded genomes under `genomes/`

Re-running with the same cache directory avoids redownloading genomes that already exist.

### Retry behavior

For genomes that are not available locally/in-cache, ZipStrain retries each download with exponential backoff (default: up to 8 attempts per genome).  
If a genome still fails after retries, it is skipped, and the reference bundle is built from successfully fetched genomes.
Parallelism for remote fetch is controlled with `--download-workers`.

If you see repeated `Too Many Requests` errors on large runs:

- lower `--download-workers` (for example `1` or `2`)
- increase `--download-retries` (for example `8`)
- increase `--retry-backoff-seconds` (for example `10` to `20`)

### Console summary

`build-genome-db` prints a short run summary:

- selected genomes (non-zero abundance)
- genomes already cached before the run
- new download attempts
- downloaded now / failed
- genomes available in cache after the run

The same summary is saved to `genome_db_build_report.txt` and includes explicit failed accession IDs (with error messages) when downloads fail.
