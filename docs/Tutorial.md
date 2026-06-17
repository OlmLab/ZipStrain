# Tutorial

This tutorial walks through the current ZipStrain workflow from mapped reads to comparison outputs.

It covers two current comparison methods:

1. The **standard** method, which compares profile tables directly with table operations
2. The **matrix** method, which first builds a reusable whole-genome matrix store and then compares from that store

The most important difference is not “old versus new.” It is when each method is a better fit.

- The **standard** method is easier to start with, more general, and usually the better choice when you want to compare many genomes at the same time. Its parallel work is naturally spread across genomes.
- The **matrix** method is worth the extra setup when you expect repeated comparison runs across many samples but only need one genome or a small set of genomes at a time. Its parallel work is naturally spread across sample pairs once the target genome matrices are loaded.

## Before You Start

You should already have one of these two starting points:

1. mapped BAM files for your samples, plus the reference bundle they were mapped against
2. raw reads for your samples, plus enough information to build the reference bundle first

Important:

- BAM files used for ZipStrain profiling should be coordinate-sorted.
- In practice, profile generation assumes BAM records are ordered by reference coordinate.
- If your BAMs are not sorted, sort them before profiling.

The reference bundle means:

- a reference genome FASTA
- an STB file mapping scaffolds to genomes

You can provide that bundle directly, or generate it from Sylph abundance output.

Optionally, you may also have:

- a Prodigal-style gene FASTA
- or a scaffold-relative gene range table

If you need installation details first, see [installation](./installation.md).
If you need pipeline-specific details, see [NextflowPipeline](./NextflowPipeline.md).
If you want command-by-command help, see [cli](./cli.md).

## Workflow Overview

A typical ZipStrain project now looks like this:

1. Prepare a reference bundle, or build one from Sylph.
2. Map reads to that reference genome database.
3. Build profiling assets.
4. Build a null model.
5. Generate one profile parquet per sample.
6. Choose a comparison route:
   - standard compare directly from profile parquets
   - matrix compare from a prebuilt matrix store
7. Export or analyze the results.

If you are starting from BAM files that are already mapped against the right reference bundle, you can skip the bundle-building and mapping steps.

## Route Selection

If you are deciding which route to use before doing anything else, use this rule of thumb:

- Choose **standard** if you want the simplest path from profiles to results, or if you want to compare many genomes in one run.
- Choose **matrix** if you expect to reuse the same profile cohort for repeated comparison jobs centered on one genome or a small set of genomes.

Use the standard workflow when:

- you want to compare many genomes at the same time
- you only need a few comparisons
- you want the fewest moving parts
- you want direct profile-table operations without building another artifact first
- you want the most general route across genomes and workflows

Use the matrix workflow when:

- you have many samples against the same reference set
- you mostly care about one genome or a small set of genomes per run
- you want resumable all-vs-all comparison runs
- you expect to append new samples over time
- you want one prebuilt store that can be reused for many comparison jobs
- you want parallel work to concentrate on sample pairs rather than genome fan-out

In practice:

- standard compare is simpler for small ad hoc work
- standard compare usually scales better when one run needs many genomes at once
- matrix compare is usually not worth it for only a few samples, because you first have to build the matrix store
- matrix compare becomes attractive once repeated one-genome or few-genome comparison cost becomes the bottleneck

One important clarification:

- the **standard comparison route** can be run in two practical ways:
  - a Python/CLI workflow
  - a Nextflow workflow
- the **matrix comparison route** is currently a CLI workflow built on top of standard profile outputs

## Worked Example Map

If you want a concrete starting point instead of reading the whole reference flow first, use one of these:

1. [Worked Example A: Standard workflow with the Python CLI](#worked-example-a-standard-workflow-with-the-python-cli)
2. [Worked Example B: Standard workflow with Nextflow](#worked-example-b-standard-workflow-with-nextflow)
3. [Worked Example C: Matrix workflow for repeated all-vs-all comparison](#worked-example-c-matrix-workflow-for-repeated-all-vs-all-comparison)

## Step 1: Prepare the Reference Bundle and Map Reads

ZipStrain expects reads to be mapped to a concatenated reference genome set.
The output of this step should be BAM files.

There are two normal ways to get there:

1. start with an existing reference FASTA and STB
2. build the reference bundle from Sylph, then map reads

### Option A: Use an Existing Reference Bundle

If you already have:

- `reference_genomes.fna`
- `reference_genomes.stb`

you can map reads directly with either the CLI/your own mapper or the Nextflow pipeline.

Example with the Nextflow pipeline:

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

### Option B: Build the Reference Bundle from Sylph

If you do not already have a prepared reference bundle, ZipStrain can build one from Sylph abundance output.

You can do that in either workflow style:

- with the Python CLI using `zipstrain utilities build-genome-db`
- or inside Nextflow by omitting `--reference_genome` and letting the pipeline build the bundle from Sylph

CLI example:

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table sylph_abundance.csv \
  --cache-dir genome_cache \
  --output-dir reference_bundle
```

This creates:

- `reference_bundle/reference_genomes.fna`
- `reference_bundle/reference_genomes.stb`

Nextflow example:

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

In that mode, Nextflow runs Sylph, builds the reference bundle, and then maps reads.

### Reads Table Example

```csv
sample_name,reads1,reads2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz
```

This step produces BAM files for each sample.

## Step 2: Prepare Profiling Assets

If you are using the CLI profile workflow, generate the supporting assets once:

```bash
zipstrain utilities prepare_profiling \
  --reference-fasta reference_genomes.fna \
  --gene-fasta reference_genes.fna \
  --stb-file reference_genomes.stb \
  --output-dir profiling_assets
```

This creates:

- `genomes_bed_file.bed`
- `genome_lengths.parquet`
- `gene_range_table.tsv`
- `null_model.parquet`
- `profiling_contract.json`

Notes:

- `--gene-fasta` is optional, but recommended if you want gene-aware outputs later.
- The generated `gene_range_table.tsv` is scaffold-relative and can also be reused by the matrix workflow.
- The generated `null_model.parquet` is ready to use for profiling immediately.
- The generated `profiling_contract.json` lets later profiling runs stamp reference/gene/null-model/STB hashes into profile metadata.

If you want non-default null-model parameters, set them directly on `prepare_profiling` with `--error-rate`, `--max-total-reads`, and `--p-threshold`.

## Step 3: Generate Profiles

Profiles are parquet tables with nucleotide counts at covered positions.
Each sample profile always contains these columns:

```text
chrom | genome | gene | pos | A | C | G | T
```

Where:

- `chrom` is the scaffold
- `genome` is the genome name
- `gene` is the gene label if gene context is available
- `pos` is the coordinate on the scaffold
- `A/C/G/T` are nucleotide counts

If you provide `--reference-fasta` during profiling, the profile also includes:

```text
ref_base_bitmask
```

`ref_base_bitmask` is the reference base encoded as a one-hot bitmask:
  - `1` = `A`
  - `2` = `C`
  - `4` = `G`
  - `8` = `T`
  - `0` = non-ACGT or unknown reference base

This choice also affects the stat tables:

- with `--reference-fasta`: `*_gene_stats.parquet` and `*_genome_stats.parquet` include `ref_ani`
- without `--reference-fasta`: those stat tables do not include `ref_ani`

Profiling input requirement:

- the BAM files passed into profiling should be coordinate-sorted BAMs
- this applies to both the CLI profile workflow and the Nextflow profile workflow
- if BAMs are unsorted, sort and index them before running profiling

### CLI Profile Workflow

Prepare a CSV of BAM files:

```csv
sample_name,bamfile
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
```

Then run:

```bash
zipstrain profile \
  --input-table bams.csv \
  --reference-fasta reference_genomes.fna \
  --stb-file reference_genomes.stb \
  --null-model profiling_assets/null_model.parquet \
  --profiling-contract profiling_assets/profiling_contract.json \
  --gene-range-table profiling_assets/gene_range_table.tsv \
  --bed-file profiling_assets/genomes_bed_file.bed \
  --genome-length-file profiling_assets/genome_lengths.parquet \
  --run-dir out_profile
```

Notes:

- `--reference-fasta` is optional on `zipstrain profile`
- if you provide it, profiles contain `ref_base_bitmask` and stat tables contain `ref_ani`
- if you omit it, profiling still works, but those reference-aware fields are omitted
- `prepare_profiling` still requires a reference FASTA because it builds the BED, genome lengths, and profiling contract assets

### Nextflow Profile Workflow

```bash
nextflow run zipstrain.nf \
  --mode profile \
  --input_table bams.csv \
  --reference_genome reference_genomes.fna \
  --gene_file reference_genes.fna \
  --stb reference_genomes.stb \
  --output_dir out_profile \
  -c conf.config \
  -profile docker \
  -resume
```

Outputs include:

- `*_profile.parquet`
- `*_genome_stats.parquet`
- `*_gene_stats.parquet`

Because the Nextflow profile workflow takes `--reference_genome`, these outputs normally include the reference-aware fields:

- profiles include `ref_base_bitmask`
- gene/genome stat tables include `ref_ani`

## Step 5A: Standard Comparison Route

This route compares profile parquets directly.
It is the better fit when one run needs to evaluate many genomes at once, because the work parallelizes efficiently across genomes.

### Standard Route Inputs

The standard route expects standard ZipStrain profile parquets with required columns:

```text
chrom | genome | gene | pos | A | C | G | T
```

Optional profile column:

```text
ref_base_bitmask
```

`ref_base_bitmask` is not required for standard pairwise compare. It is only present when profiling was run with `--reference-fasta`.

Depending on which command you use, the immediate inputs are:

- one pair of profile parquets for `single_compare_genome`
- a profile DB plus direct compare arguments for `zipstrain compare genomes`
- optionally an STB when you want the requested genome universe carried through explicitly

### Standard Route Outputs

Genome-level standard comparison tables contain one row per `sample_1`, `sample_2`, `genome` combination.

When ANI is requested, the output columns include:

- `sample_1`
- `sample_2`
- `genome`
- `total_positions`
- `share_allele_pos`
- `genome_pop_ani`

When IBS is also requested, add:

- `max_consecutive_length`

When `identical_genes` is also requested, add:

- `shared_genes_count`
- `identical_gene_count`
- `perc_id_genes`

### Single Pair Compare

For one pair of samples:

```bash
zipstrain utilities single_compare_genome \
  --mpileup-contig-1 sample1_profile.parquet \
  --mpileup-contig-2 sample2_profile.parquet \
  --stb-file reference_genomes.stb \
  --calculate ani+ibs \
  --engine duckdb \
  --output-file sample1_vs_sample2.parquet
```

### Large Standard Compare Runs

For many standard comparisons, the direct profile-table workflow is still:

1. build a profile DB
2. optionally point at an existing comparison parquet
3. run `zipstrain compare genomes`

Example:

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file profile_db.parquet
```

```bash
zipstrain compare genomes \
  --profile-db profile_db.parquet \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file reference_genomes.stb \
  --run-dir compare_run \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

This route is still valid, but it is not the best option once you repeatedly compare many samples against the same reference set.
It remains especially strong when you need broad genome scope in the same run.

## Step 5B: Matrix Comparison Route

This route converts standard sample profiles into a reusable matrix store, then compares all non-redundant sample pairs from that store.
It is the better fit when you have many samples and want repeated comparison work on one genome or a small set of genomes, because the hot parallel work is across pairs rather than across genome fan-out.

### Why Use the Matrix Route

The matrix route is useful when:

- you want to compare many samples repeatedly
- you usually care about one genome or a small set of genomes per compare run
- you want appendable sample growth over time
- you want resumable compare state in a DuckDB result database
- you want ANI, IBS, and optional gene ANI from the same store

### Matrix Route Overview

The matrix route has four steps:

1. build the matrix store
2. optionally append new samples later
3. run matrix compare
4. export the compare DB to parquet

### Build the Matrix Store

Expected input:

- a directory of standard ZipStrain profile parquets
- each profile parquet should contain:

```text
chrom | genome | gene | pos | A | C | G | T
```

```bash
zipstrain utilities build-matrix-db \
  --profile-dir out_profile \
  --output-file matrix_db.h5 \
  --gene-range-table profiling_assets/gene_range_table.tsv \
  --memory-limit-gb 16
```

What this does:

- scans all standard profile parquets in `out_profile`
- builds one HDF5-backed matrix store
- keeps one whole-genome matrix per sample per genome
- stores gene coordinate metadata if `--gene-range-table` is provided

Important notes:

- `matrix_db.h5` is the current matrix-store format
- dense storage is the default; add `--sparse` if you want sparse HDF5 storage
- this store is appendable on the sample axis
- if gene ranges are included here, matrix compare can also compute gene ANI later

This step does not yet produce comparison result tables.
Its output is the reusable matrix store itself: `matrix_db.h5`.

### Append New Samples Later

If you later generate more profiles against the same reference set:

```bash
zipstrain utilities append-matrix-db \
  --profile-dir out_profile_new \
  --matrix-db-file matrix_db.h5 \
  --memory-limit-gb 16
```

This validates that the new profiles match the stored genome/scaffold contract and appends only new sample rows.

### Convert a Legacy Matrix DuckDB

If you already have an old matrix DuckDB from a previous workflow:

```bash
zipstrain utilities matrix-db-to-hdf5 \
  --matrix-db-file legacy_matrix.duckdb \
  --output-file matrix_db.h5
```

This is only needed for older matrix builds. New matrix workflows should start directly from `build-matrix-db`.

### Run Matrix Compare

Run pairwise sample comparison from the matrix store.
If you want to focus on one genome, add `--genome <genome_id>`.
If you want a small set of genomes, run a few scoped matrix-compare jobs, one genome at a time:

```bash
zipstrain utilities matrix-compare \
  --matrix-db-file matrix_db.h5 \
  --output-file matrix_compare.duckdb \
  --memory-limit-gb 16 \
  --loader-executor thread \
  --writer-executor thread \
  --backend torch-cpu \
  --calculate all
```

What `--calculate` means here:

- `ani`: genome ANI only
- `ani+ibs`: genome ANI plus IBS
- `+gene` or `gene`: gene ANI, which also implies ANI
- `all`: `ani+ibs`, and also `gene` when the matrix store contains gene annotations

Backend choices:

- `numpy`: simple CPU path
- `torch-cpu`: torch on CPU
- `torch-cuda`: NVIDIA GPU
- `torch-mps`: Apple Silicon GPU
- `torch`: auto-select device

Recommended starting points:

- Apple Silicon: `--backend torch-mps`
- NVIDIA GPU: `--backend torch-cuda`
- debugging or CPU-only: `--backend torch-cpu`

Important operational notes:

- the compare output is a DuckDB database, not parquet
- the compare DB is resumable
- rerunning the same command on the same output file only processes unfinished sample pairs
- `--memory-limit-gb` is the main throughput control for target block size

Expected input:

- one matrix store such as `matrix_db.h5`
- optionally one genome scope via `--genome`
- optionally gene metadata already embedded in that store if you want gene ANI

Expected output:

- a DuckDB database such as `matrix_compare.duckdb`
- a genome-level result table named `matrix_compare_results`
- optionally a gene-level result table named `matrix_compare_gene_results`

`matrix_compare_results` contains one row per `sample_1`, `sample_2`, `genome` combination.

When ANI is requested, the table includes:

- `sample_1`
- `sample_2`
- `genome`
- `total_positions`
- `share_allele_pos`
- `genome_pop_ani`

When IBS is also requested, add:

- `max_consecutive_length`

`matrix_compare_gene_results` contains one row per `sample_1`, `sample_2`, `genome`, `gene` combination with:

- `sample_1`
- `sample_2`
- `genome`
- `gene`
- `gene_pop_ani`

### Export Matrix Compare Results

Export genome-level rows:

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare.parquet
```

Export gene-level ANI rows:

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare_gene.parquet \
  --table gene
```

If the compare DB does not contain gene results, the gene export command raises a clear error.

These export commands are the normal way to extract the relevant result tables from the compare DB:

- export `matrix_compare_results` with the default `matrix-compare-export` command
- export `matrix_compare_gene_results` with `--table gene`

## Understanding Standard and Matrix Inputs and Outputs

### Standard Compare Input

The standard route reads classic ZipStrain profile parquets directly.
Those parquets contain per-position nucleotide counts and annotation columns:

- `chrom`
- `genome`
- `gene`
- `pos`
- `A`
- `C`
- `G`
- `T`

### Standard Compare Output

Standard genome comparison outputs are parquet tables with one row per sample pair and genome.

Core genome ANI columns:

- `sample_1`
- `sample_2`
- `genome`
- `total_positions`
- `share_allele_pos`
- `genome_pop_ani`

Optional additional columns:

- `max_consecutive_length` for IBS
- `shared_genes_count`
- `identical_gene_count`
- `perc_id_genes` for identical-gene summaries

### Matrix Store Input

The matrix store is:

- built from standard profile parquets
- stored as HDF5
- organized by genome
- reused across many compare runs
- stored densely by default, or sparsely if `--sparse` is used during build or conversion

### Matrix Compare Output

The compare DB contains:

- genome-level results in `matrix_compare_results`
- optional gene-level results in `matrix_compare_gene_results`
- resume metadata for completed sample-pair/genome work

Genome parquet exports contain name-based columns such as:

- `sample_1`
- `sample_2`
- `genome`
- `total_positions`
- `share_allele_pos`
- `genome_pop_ani`
- `max_consecutive_length` when IBS is requested

Gene parquet exports contain:

- `sample_1`
- `sample_2`
- `genome`
- `gene`
- `gene_pop_ani`

## Choosing Between Standard and Matrix Compare

Use standard compare when:

- you want to compare many genomes in the same run
- you only need a limited number of pairwise comparisons
- you want direct profile-to-profile comparison without building another store
- your existing workflow is already based on direct profile comparisons or Nextflow standard compare jobs

Use matrix compare when:

- you want to compare all samples against all samples repeatedly
- you are usually working on one genome or a small set of genomes at a time
- you want resumability at the compare-DB level
- you expect to append new profiles over time
- you want one reusable comparison substrate for ANI, IBS, and gene ANI

## Current Limitations and Practical Notes

- The matrix route is currently a CLI workflow, not a Nextflow mode.
- HDF5 matrix input currently requires a torch backend.
- Gene ANI from matrix compare requires that the matrix store was built with `--gene-range-table`.
- The matrix route does not replace the standard `zipstrain compare` workflow everywhere yet; both routes remain supported.

## Recommended Starting Point

For a new project:

1. generate standard sample profiles
2. if your main analyses will span many genomes at once, start with the standard compare route
3. if your main analyses will revisit one genome or a small set of genomes across many samples, build one matrix store and use the matrix route
4. append new samples into the matrix store later only if you chose the matrix route

This gives you the right optimization for your actual workload:

- standard route: simpler direct profile-based comparison across many genomes
- matrix route: reusable comparison input plus resumable output state for repeated pair-heavy genome-specific work

## Worked Example A: Standard Workflow with the Python CLI

This is the most direct standard ZipStrain workflow when you already have mapped BAM files and want to stay in the CLI.

If you do not already have `reference_genomes.fna` and `reference_genomes.stb`, build them first with `zipstrain utilities build-genome-db --tool sylph ...` as shown earlier in Step 1.

### When to use this route

Use this when:

- you already have BAM files
- you want profile-parquet outputs directly
- you only need a modest number of comparisons
- you want the standard profile-table compare behavior without building a matrix store

### Example layout

```text
project/
├── mapped_bams/
│   ├── sample1.bam
│   ├── sample2.bam
│   └── sample3.bam
├── reference/
│   ├── reference_genomes.fna
│   ├── reference_genes.fna
│   └── reference_genomes.stb
└── outputs/
```

### 1. Prepare profiling assets once

```bash
zipstrain utilities prepare_profiling \
  --reference-fasta reference/reference_genomes.fna \
  --gene-fasta reference/reference_genes.fna \
  --stb-file reference/reference_genomes.stb \
  --output-dir outputs/profiling_assets
```

`prepare_profiling` also writes:

- `outputs/profiling_assets/null_model.parquet`
- `outputs/profiling_assets/profiling_contract.json`

### 2. Make the BAM input table

`bams.csv`

```csv
sample_name,bamfile
sample1,/abs/path/project/mapped_bams/sample1.bam
sample2,/abs/path/project/mapped_bams/sample2.bam
sample3,/abs/path/project/mapped_bams/sample3.bam
```

### 3. Generate profiles

```bash
zipstrain profile \
  --input-table bams.csv \
  --stb-file reference/reference_genomes.stb \
  --null-model outputs/profiling_assets/null_model.parquet \
  --profiling-contract outputs/profiling_assets/profiling_contract.json \
  --gene-range-table outputs/profiling_assets/gene_range_table.tsv \
  --bed-file outputs/profiling_assets/genomes_bed_file.bed \
  --genome-length-file outputs/profiling_assets/genome_lengths.parquet \
  --run-dir outputs/profile_run
```

This gives you one profile parquet per sample plus gene and genome stats tables.

### 4. Build the profile DB table

Create `profiles.csv`:

```csv
profile_name,profile_location
sample1,/abs/path/project/outputs/profile_run/sample1_profile.parquet
sample2,/abs/path/project/outputs/profile_run/sample2_profile.parquet
sample3,/abs/path/project/outputs/profile_run/sample3_profile.parquet
```

Then build the DB:

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file outputs/profile_db.parquet
```

### 5. Run the standard compare

```bash
zipstrain compare genomes \
  --profile-db outputs/profile_db.parquet \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file reference/reference_genomes.stb \
  --run-dir outputs/standard_compare \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

### Result

This route gives you standard profile-based comparison outputs without building a matrix store. It is the best tutorial example when you want the simplest direct route from profiles to compare results, especially when the same run needs many genomes at once.

## Worked Example B: Standard Workflow with Nextflow

This is the best route when you want orchestration, batch execution, resumability at the pipeline level, and easier deployment on clusters.

### When to use this route

Use this when:

- you want mapping, profiling, and comparison wrapped in one pipeline framework
- you are running on a cluster or container-based environment
- you want the standard workflow but do not want to script every CLI step yourself

### 1. Prepare the read input table

`reads.csv`

```csv
sample_name,reads1,reads2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz
sample3,/data/sample3_R1.fastq.gz,/data/sample3_R2.fastq.gz
```

### 2. Map reads

If you already have a reference bundle, use it directly:

```bash
nextflow run zipstrain.nf \
  --mode map_reads \
  --input_type local \
  --input_table reads.csv \
  --reference_genome reference/reference_genomes.fna \
  --stb reference/reference_genomes.stb \
  --output_dir outputs/nf_map \
  -c conf.config \
  -profile docker \
  -resume
```

If you do not already have a prepared reference bundle, use the Sylph-backed route from Step 1 Option B instead and let Nextflow build it first.

### 3. Generate profiles

Build a BAM table that points at the BAM outputs from the map step, then run:

```bash
nextflow run zipstrain.nf \
  --mode profile \
  --input_table bams.csv \
  --reference_genome reference/reference_genomes.fna \
  --gene_file reference/reference_genes.fna \
  --stb reference/reference_genomes.stb \
  --output_dir outputs/nf_profile \
  -c conf.config \
  -profile docker \
  -resume
```

This writes profile parquets under `outputs/nf_profile/profiles/`.

### 4. Build the profile list for compare

`profiles.csv`

```csv
sample_name,profile_location
sample1,/abs/path/project/outputs/nf_profile/profiles/sample1_profile.parquet
sample2,/abs/path/project/outputs/nf_profile/profiles/sample2_profile.parquet
sample3,/abs/path/project/outputs/nf_profile/profiles/sample3_profile.parquet
```

### 5. Run genome comparison in Nextflow

```bash
nextflow run zipstrain.nf \
  --mode compare_genomes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference/reference_genomes.stb \
  --compare_engine polars \
  --compare_genome_scope all \
  --compare_calculate ani+ibs+identical_genes \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir outputs/nf_compare_genomes \
  -c conf.config \
  -profile docker \
  -resume
```

### 6. Optional gene compare

```bash
nextflow run zipstrain.nf \
  --mode compare_genes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference/reference_genomes.stb \
  --compare_engine polars \
  --compare_gene_scope all:all \
  --compare_ani_method popani \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir outputs/nf_compare_genes \
  -c conf.config \
  -profile docker \
  -resume
```

### Result

This route keeps you in the standard profile-based world, but uses Nextflow to orchestrate the work instead of calling the underlying CLI tools manually.

## Worked Example C: Matrix Workflow for Repeated All-vs-All Comparison

This is the best route when you expect to compare many samples repeatedly against the same reference set, but each run is usually focused on one genome or a small set of genomes.

### When to use this route

Use this when:

- you already have standard profile parquets
- you want resumable all-vs-all compare runs
- you expect to add new samples over time
- you usually care about one genome or a small set of genomes per run
- you want genome ANI, IBS, and optional gene ANI from one reusable matrix store

### 1. Start from standard profiles

This route assumes you already produced profile parquets using either:

- the Python CLI profile workflow
- the Nextflow profile workflow

For example:

```text
outputs/profile_run/sample1_profile.parquet
outputs/profile_run/sample2_profile.parquet
outputs/profile_run/sample3_profile.parquet
```

### 2. Build the matrix store

```bash
zipstrain utilities build-matrix-db \
  --profile-dir outputs/profile_run \
  --output-file outputs/matrix_db.h5 \
  --gene-range-table outputs/profiling_assets/gene_range_table.tsv \
  --memory-limit-gb 16
```

### 3. Run matrix compare

```bash
zipstrain utilities matrix-compare \
  --matrix-db-file outputs/matrix_db.h5 \
  --output-file outputs/matrix_compare.duckdb \
  --genome target_genome_name \
  --memory-limit-gb 16 \
  --loader-executor thread \
  --writer-executor thread \
  --backend torch-cpu \
  --calculate all
```

Typical backend substitutions:

- Apple Silicon: `--backend torch-mps`
- NVIDIA GPU: `--backend torch-cuda`
- CPU only: `--backend torch-cpu`

If you want all genomes instead, omit `--genome`.

### 4. Export genome-level results

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file outputs/matrix_compare.duckdb \
  --output-file outputs/matrix_compare.parquet
```

### 5. Export gene-level results

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file outputs/matrix_compare.duckdb \
  --output-file outputs/matrix_compare_gene.parquet \
  --table gene
```

### 6. Append new samples later

When more profiles arrive against the same reference set:

```bash
zipstrain utilities append-matrix-db \
  --profile-dir outputs/new_profile_run \
  --matrix-db-file outputs/matrix_db.h5 \
  --memory-limit-gb 16
```

Then rerun `matrix-compare` against the same compare DB. Already completed sample pairs remain marked complete, and only unfinished work is processed.

### Result

This route turns standard per-sample profile parquets into a reusable comparison substrate. It is the best example for large or growing cohorts where repeated pair-heavy work on one genome or a small set of genomes matters more than one-time matrix construction.

The relevant result tables can then be extracted with:

- `zipstrain utilities matrix-compare-export --matrix-compare-db-file outputs/matrix_compare.duckdb --output-file outputs/matrix_compare.parquet`
- `zipstrain utilities matrix-compare-export --matrix-compare-db-file outputs/matrix_compare.duckdb --output-file outputs/matrix_compare_gene.parquet --table gene`

## Optional: Add a New Sample Later and Run Only the Remaining Pairs

This is the practical “my cohort grew by one sample” scenario.
The two methods handle it differently.

### Standard Method

The standard method does not maintain a dedicated resumable compare database in the same way the matrix route does. Instead, you:

1. add the new sample profile to the profile DB
2. point the next run at the previous comparison table
3. regenerate only the remaining pairs
4. rerun compare with the same direct arguments

#### 1. Update the profile DB input table

Extend your `profiles.csv` to include the new sample:

```csv
profile_name,profile_location
sample1,/abs/path/project/outputs/profile_run/sample1_profile.parquet
sample2,/abs/path/project/outputs/profile_run/sample2_profile.parquet
sample3,/abs/path/project/outputs/profile_run/sample3_profile.parquet
sample4,/abs/path/project/outputs/profile_run/sample4_profile.parquet
```

Rebuild the profile DB:

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file outputs/profile_db.parquet
```

#### 2. Rebuild the remaining-pairs table with the current comparison parquet

Point `--comp-db-file` at the existing genome comparison parquet from your previous standard run. The path below is an example placeholder; use the actual parquet produced by your earlier standard compare workflow:

```bash
zipstrain utilities to-complete-table \
  --profile-db outputs/profile_db.parquet \
  --comp-db-file outputs/standard_compare/genome_compare_results.parquet \
  --output-file outputs/remaining_pairs.csv
```

This gives you the not-yet-completed sample pairs implied by the updated profile DB plus existing comparison parquet.

#### 3. Rerun standard compare

```bash
zipstrain compare genomes \
  --profile-db outputs/profile_db.parquet \
  --comp-db-file outputs/standard_compare/genome_compare_results.parquet \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file reference/reference_genomes.stb \
  --run-dir outputs/standard_compare_update \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

In other words, the standard route handles incremental updates by carrying forward the previous comparison parquet directly, without a separate comparison-config JSON layer.

### Matrix Method

The matrix method is simpler here because resumability is built into the compare DB itself.

#### 1. Append the new sample into the matrix store

```bash
zipstrain utilities append-matrix-db \
  --profile-dir outputs/new_profile_run \
  --matrix-db-file outputs/matrix_db.h5 \
  --memory-limit-gb 16
```

#### 2. Rerun matrix compare against the same compare DB

```bash
zipstrain utilities matrix-compare \
  --matrix-db-file outputs/matrix_db.h5 \
  --output-file outputs/matrix_compare.duckdb \
  --memory-limit-gb 16 \
  --loader-executor thread \
  --writer-executor thread \
  --backend torch-cpu \
  --calculate all
```

Because the compare DB tracks completed sample-pair/genome work, the rerun only processes unfinished work. That is one of the main reasons the matrix route is attractive for growing cohorts.

## Appendix: Input and Output Formats

This section is a practical reference for what the core ZipStrain inputs and outputs should look like on disk.

The main idea is:

- plain-text helper tables are usually flat tab-separated or comma-separated tables
- parquet outputs are typed columnar tables
- some inputs require headers, and some do not

If a file is malformed, the most common causes are:

- wrong delimiter
- a header row where ZipStrain expects no header
- missing required columns
- columns in the wrong order

### Reference FASTA

What it is:

- a standard multi-record FASTA of scaffolds/contigs used as the mapping reference

Expected form:

- plain-text FASTA
- each record header starts with `>`
- the scaffold name used by ZipStrain is the first token in the FASTA header

Example:

```text
>contigA
ACGTACGTACGT
>contigB
TTGCAATGCAAA
```

Quick ways to make it:

- provide your existing reference FASTA directly
- or generate a reference bundle with:

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table sylph_abundance.csv \
  --cache-dir genome_cache \
  --output-dir reference_bundle
```

### STB File

What it is:

- a scaffold-to-genome mapping table

Expected form:

- plain-text TSV
- no header
- exactly 2 columns
- column order:
  - `scaffold`
  - `genome`

Example:

```text
contigA	genome_1
contigB	genome_1
contigC	genome_2
```

Important:

- do not include a header row for the main profiling and compare workflows
- scaffold names must match the reference FASTA scaffold names

Quick ways to make it:

- produced by `zipstrain utilities build-genome-db`
- or generated from a directory of genome FASTAs with:

```bash
zipstrain utilities generate_stb \
  --genomes-dir-file genomes_dir \
  --output-file reference.stb
```

### Gene FASTA

What it is:

- a nucleotide gene FASTA, typically from Prodigal

Expected form:

- plain-text FASTA
- for the main `gene_range_table` builder, the headers should follow the Prodigal-style format that encodes scaffold, start, and end

Example:

```text
>contigA_1 # 1 # 300 # 1 # ID=1_1;partial=00
ATG...
>contigA_2 # 450 # 900 # -1 # ID=1_2;partial=00
ATG...
```

Quick ways to make it:

- provide your existing Prodigal nucleotide FASTA
- or run Prodigal yourself, for example:

```bash
prodigal -i reference_genomes.fna -d reference_genes.fna -p meta
```

### `gene_range_table.tsv`

What it is:

- a scaffold-relative gene interval table used during profiling and by the matrix workflow

Expected form:

- plain-text TSV
- no header
- exactly 4 columns
- column order:
  - `gene`
  - `scaffold`
  - `start`
  - `end`

Example:

```text
contigA_1	contigA	1	300
contigA_2	contigA	450	900
contigB_1	contigB	10	250
```

Important:

- this is a flat table, not JSON and not nested
- the coordinates are scaffold-relative integer coordinates
- for the main profiling workflow, ZipStrain expects no header row
- if you were thinking of a 3-column file here, that is a different object; the current `gene_range_table.tsv` used by profiling is 4 columns

Quick ways to make it:

- created automatically by `zipstrain utilities prepare_profiling`
- or created directly from a Prodigal-style gene FASTA with:

```bash
zipstrain utilities gene-range-table \
  --gene-file reference_genes.fna \
  --output-file gene_range_table.tsv
```

### `genomes_bed_file.bed`

What it is:

- the list of scaffold intervals to profile

Expected form:

- plain-text TSV/BED-style file
- no header
- exactly 3 columns
- column order:
  - `scaffold`
  - `start`
  - `end`

Example:

```text
contigA	0	500000
contigA	500000	1000000
contigB	0	240123
```

Important:

- this follows BED-style interval logic: scaffold plus integer start/end coordinates
- long scaffolds may appear in multiple rows
- scaffold names must match the reference FASTA and STB

Quick ways to make it:

- created automatically by `zipstrain utilities prepare_profiling`
- or created directly from a reference FASTA with:

```bash
zipstrain utilities make_bed \
  --db-fasta-dir reference_genomes.fna \
  --output-file genomes_bed_file.bed
```

### `genome_lengths.parquet`

What it is:

- per-genome total reference length derived from the BED plus STB

Expected columns:

```text
genome | genome_length
```

Example rows:

```text
genome_1 | 2384921
genome_2 | 1912044
```

Quick ways to make it:

- created automatically by `zipstrain utilities prepare_profiling`
- or created from an STB plus BED with:

```bash
zipstrain utilities get_genome_lengths \
  --stb-file reference.stb \
  --bed-file genomes_bed_file.bed \
  --output-file genome_lengths.parquet
```

### `null_model.parquet`

What it is:

- the sequencing-error threshold table used during profile adjustment

Expected columns:

```text
cov | max_error_count
```

Example rows:

```text
0 | 0
1 | 0
2 | 0
...
```

Quick ways to make it:

- created automatically by `zipstrain utilities prepare_profiling`
- or created directly with:

```bash
zipstrain utilities build-null-model \
  --error-rate 0.001 \
  --max-total-reads 10000 \
  --p-threshold 0.05 \
  --output-file null_model.parquet
```

### `profiling_contract.json`

What it is:

- a small JSON metadata manifest describing the reference-side assets used for profiling

Expected form:

- JSON object
- current keys:
  - `reference_hash`
  - `gene_hash`
  - `null_model_hash`
  - `stb_hash`

Example:

```json
{
  "gene_hash": "a1b2c3...",
  "null_model_hash": "d4e5f6...",
  "reference_hash": "112233...",
  "stb_hash": "445566..."
}
```

Quick ways to make it:

- created automatically by `zipstrain utilities prepare_profiling`

### BAM Input Table for `zipstrain profile`

What it is:

- the sample table for batch profiling from already mapped BAM files

Expected form:

- CSV
- header required
- required columns:
  - `sample_name`
  - `bamfile`

Example:

```csv
sample_name,bamfile
sample1,/abs/path/sample1.bam
sample2,/abs/path/sample2.bam
sample3,/abs/path/sample3.bam
```

Important:

- this table does require a header
- BAMs should already be coordinate-sorted and indexed

Quick way to make it with shell:

```bash
printf "sample_name,bamfile\n" > bams.csv
find /path/to/bams -name '*.bam' | sort | while read -r bam; do
  sample=$(basename "$bam" .bam)
  printf "%s,%s\n" "$sample" "$bam" >> bams.csv
done
```

### Reads Table for Nextflow Mapping

What it is:

- the sample table for starting from raw reads

Expected form:

- CSV
- header required

Paired-end example:

```csv
sample_name,reads1,reads2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz
```

Single-end example:

```csv
sample_name,reads1
sample1,/data/sample1.fastq.gz
sample2,/data/sample2.fastq.gz
```

SRA example:

```csv
Run
SRR12345678
SRR12345679
```

### Classic Profile Parquet

What it is:

- the main per-sample profile output

Expected columns:

```text
chrom | genome | gene | pos | A | C | G | T
```

Notes:

- `ref_base_bitmask` is an optional extra column when profiling includes `--reference-fasta`
- rows are expected to be sorted by `(chrom, pos)`
- `A/C/G/T` are post-adjustment counts

Example:

```text
contigA | genome_1 | contigA_1 | 1 | 10 | 0 | 0 | 0 | 1
contigA | genome_1 | contigA_1 | 2 | 0  | 7 | 0 | 0 | 2
```

This example assumes profiling was run with `--reference-fasta`, so the final value shown is `ref_base_bitmask`.

Quick way to make it:

- `zipstrain utilities profile-single`
- or batch mode with `zipstrain profile`

### `*_genome_stats.parquet`

What it is:

- per-sample genome summary statistics

Expected columns:

```text
genome | coverage | breadth | genome_length | gap_mean | gap_std | 5x_cov_sites | heterogeneity | ber | fug | reads_mapped
```

Notes:

- if profiling used `--reference-fasta`, this file also includes:

```text
ref_ani
```

- `ref_ani` is the percent of covered positions that still contain the reference allele after sequence-error adjustment

### `*_gene_stats.parquet`

What it is:

- per-sample gene summary statistics

Expected columns:

```text
genome | gene | length | breadth | coverage
```

Notes:

- if profiling used `--reference-fasta`, this file also includes:

```text
ref_ani
```

- `ref_ani` is the percent of covered positions in the gene that still contain the reference allele after sequence-error adjustment

### `profiles.csv` for Building a Profile DB

What it is:

- the normal CSV input used by `zipstrain utilities build-profile-db`

Expected form:

- CSV
- header required
- required columns:
  - `profile_name`
  - `profile_location`

Example:

```csv
profile_name,profile_location
sample1,/abs/path/outputs/profile_run/sample1_profile.parquet
sample2,/abs/path/outputs/profile_run/sample2_profile.parquet
sample3,/abs/path/outputs/profile_run/sample3_profile.parquet
```

Quick way to make it with shell:

```bash
printf "profile_name,profile_location\n" > profiles.csv
find /abs/path/outputs/profile_run -name '*_profile.parquet' | sort | while read -r pf; do
  sample=$(basename "$pf" _profile.parquet)
  printf "%s,%s\n" "$sample" "$pf" >> profiles.csv
done
```

Then build the DB:

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file profile_db.parquet
```

### `profile_db.parquet`

What it is:

- the parquet form of the profile database

Expected columns:

```text
profile_name | profile_location
```

Quick way to make it:

- built from `profiles.csv` using `zipstrain utilities build-profile-db`

### Remaining-Pairs Table

What it is:

- the table of sample pairs still needing comparison in the standard workflow

Expected columns:

```text
sample_name_1 | sample_name_2 | profile_location_1 | profile_location_2
```

Quick way to make it:

```bash
zipstrain utilities to-complete-table \
  --profile-db profile_db.parquet \
  --comp-db-file current_compare.parquet \
  --output-file remaining_pairs.csv
```

### Genome Comparison Parquet

What it is:

- the standard genome-level comparison output

Core columns:

```text
genome | total_positions | share_allele_pos | sample_1 | sample_2
```

Common additional columns:

```text
genome_pop_ani | max_consecutive_length | shared_genes_count | identical_gene_count | perc_id_genes
```

Important:

- the exact set of metric columns depends on `--calculate`
- sample pairs are represented by `sample_1` and `sample_2`
- one row usually corresponds to one sample-pair/genome result

### Gene Comparison Parquet

What it is:

- the standard gene-level comparison output

Expected columns:

```text
genome | gene | total_positions | share_allele_pos | ani | sample_1 | sample_2
```

### Matrix Compare Export Outputs

What they are:

- the exported parquet tables from the matrix compare DuckDB database

Typical forms:

- genome export: same practical row meaning as the standard genome comparison parquet
- gene export: same practical row meaning as the standard gene comparison parquet

Quick ways to make them:

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file outputs/matrix_compare.duckdb \
  --output-file outputs/matrix_compare.parquet
```

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file outputs/matrix_compare.duckdb \
  --output-file outputs/matrix_compare_gene.parquet \
  --table gene
```

## Further Reading

- [CLI Reference](./cli.md)
- [Installation](./installation.md)
- [Nextflow Pipeline](./NextflowPipeline.md)
- [API Reference](./api.md)
