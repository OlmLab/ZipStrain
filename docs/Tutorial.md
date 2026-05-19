# Tutorial

This tutorial walks through the current ZipStrain workflow from mapped reads to comparison outputs.

It covers two current comparison methods:

1. The **standard** method, which compares profile tables directly with table operations
2. The **matrix** method, which first builds dense whole-genome matrices and then compares from that store

The most important difference is not “old versus new.” It is when each method is a better fit.

- The **standard** method is easier to start with, more general, and usually the better choice when you only have a few samples or only need a modest number of comparisons.
- The **matrix** method is worth the extra setup when you expect repeated all-vs-all comparison runs against the same reference set, especially when your analysis is centered on one genome or a small set of target genomes across many samples.

## Before You Start

You should already have one of these two starting points:

1. mapped BAM files for your samples, plus the reference bundle they were mapped against
2. raw reads for your samples, plus enough information to build the reference bundle first

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

- Choose **standard** if you want the simplest path from profiles to results.
- Choose **matrix** if you expect to reuse the same profile cohort for repeated comparison jobs.

Use the standard workflow when:

- you only need a few comparisons
- you want the fewest moving parts
- you want direct profile-table operations without building another artifact first
- you want the most general route across genomes and workflows

Use the matrix workflow when:

- you have many samples against the same reference set
- you want resumable all-vs-all comparison runs
- you expect to append new samples over time
- you want one prebuilt store that can be reused for many comparison jobs
- you are mainly focused on one genome or a small set of genomes and want repeated comparison throughput

In practice:

- standard compare is simpler for small ad hoc work
- matrix compare is usually not worth it for only a few samples, because you first have to build the matrix store
- matrix compare becomes attractive once repeated comparison cost becomes the bottleneck

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

Notes:

- `--gene-fasta` is optional, but recommended if you want gene-aware outputs later.
- The generated `gene_range_table.tsv` is scaffold-relative and can also be reused by the matrix workflow.

## Step 3: Build the Null Model

ZipStrain uses a null model for sequencing-error adjustment during profiling.
Build it once from any representative BAM file:

```bash
zipstrain utilities build-null-model \
  --bam-file /path/to/example.bam \
  --output-file null_model.parquet
```

## Step 4: Generate Profiles

Profiles are parquet tables with nucleotide counts at covered positions.
Each sample profile contains columns like:

```text
chrom | genome | gene | pos | A | C | G | T
```

Where:

- `chrom` is the scaffold
- `genome` is the genome name
- `gene` is the gene label if gene context is available
- `pos` is the coordinate on the scaffold
- `A/C/G/T` are nucleotide counts

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
  --stb-file reference_genomes.stb \
  --null-model null_model.parquet \
  --gene-range-table profiling_assets/gene_range_table.tsv \
  --bed-file profiling_assets/genomes_bed_file.bed \
  --genome-length-file profiling_assets/genome_lengths.parquet \
  --run-dir out_profile
```

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

## Step 5A: Standard Comparison Route

This route compares profile parquets directly.

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
2. build a comparison config
3. run `zipstrain compare genomes`

Example:

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file profile_db.parquet
```

```bash
zipstrain utilities build-genome-comparison-config \
  --profile-db profile_db.parquet \
  --gene-db-id ref_genes_v1 \
  --reference-genome-id ref_genomes_v1 \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file-loc reference_genomes.stb \
  --output-file genome_compare_config.json
```

```bash
zipstrain compare genomes \
  --genome-comparison-object genome_compare_config.json \
  --run-dir compare_run \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

This route is still valid, but it is not the best option once you repeatedly compare many samples against the same reference set.

## Step 5B: Matrix Comparison Route

This route converts standard sample profiles into a reusable matrix store, then compares all non-redundant sample pairs from that store.

### Why Use the Matrix Route

The matrix route is useful when:

- you want to compare many samples repeatedly
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
- keeps one dense whole-genome matrix per sample per genome
- stores gene coordinate metadata if `--gene-range-table` is provided

Important notes:

- `matrix_db.h5` is the current matrix-store format
- this store is appendable on the sample axis
- if gene ranges are included here, matrix compare can also compute gene ANI later

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

Run all-vs-all sample comparison from the matrix store:

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
- `--position-tile-size` is currently a legacy compatibility option and is not used by the current full-genome compare path

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

## Understanding Matrix Inputs and Outputs

### Matrix Store Input

The matrix store is:

- built from standard profile parquets
- stored as HDF5
- organized by genome
- reused across many compare runs

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

- you only need a limited number of pairwise comparisons
- you want direct profile-to-profile comparison without building another store
- your existing workflow is already based on comparison config objects or Nextflow standard compare jobs

Use matrix compare when:

- you want to compare all samples against all samples repeatedly
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
2. build one matrix store from those profiles
3. run matrix compare into a resumable DuckDB compare DB
4. append new samples into the matrix store as the cohort grows

This gives you:

- reusable comparison input
- resumable output state
- fast repeated cohort comparison runs
- one place to export genome and gene comparison tables

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

### 2. Build the null model once

```bash
zipstrain utilities build-null-model \
  --bam-file mapped_bams/sample1.bam \
  --output-file outputs/null_model.parquet
```

### 3. Make the BAM input table

`bams.csv`

```csv
sample_name,bamfile
sample1,/abs/path/project/mapped_bams/sample1.bam
sample2,/abs/path/project/mapped_bams/sample2.bam
sample3,/abs/path/project/mapped_bams/sample3.bam
```

### 4. Generate profiles

```bash
zipstrain profile \
  --input-table bams.csv \
  --stb-file reference/reference_genomes.stb \
  --null-model outputs/null_model.parquet \
  --gene-range-table outputs/profiling_assets/gene_range_table.tsv \
  --bed-file outputs/profiling_assets/genomes_bed_file.bed \
  --genome-length-file outputs/profiling_assets/genome_lengths.parquet \
  --run-dir outputs/profile_run
```

This gives you one profile parquet per sample plus gene and genome stats tables.

### 5. Build the profile DB table

Create `profiles.csv`:

```csv
sample_name,profile_loc
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

### 6. Build the standard compare config

```bash
zipstrain utilities build-genome-comparison-config \
  --profile-db outputs/profile_db.parquet \
  --gene-db-id tutorial_genes_v1 \
  --reference-genome-id tutorial_reference_v1 \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file-loc reference/reference_genomes.stb \
  --output-file outputs/genome_compare_config.json
```

### 7. Run the standard compare

```bash
zipstrain compare genomes \
  --genome-comparison-object outputs/genome_compare_config.json \
  --run-dir outputs/standard_compare \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

### Result

This route gives you standard profile-based comparison outputs without building a matrix store. It is the best tutorial example when you want the simplest direct route from profiles to compare results.

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
sample_names,mpileup_files
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
  --compare_genome_scope all \
  --compare_calculate ani+ibs+identical_genes \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --compare_duckdb_memory_limit 4GB \
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
  --compare_gene_scope all:all \
  --compare_ani_method popani \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --compare_duckdb_memory_limit 4GB \
  --output_dir outputs/nf_compare_genes \
  -c conf.config \
  -profile docker \
  -resume
```

### Result

This route keeps you in the standard profile-based world, but uses Nextflow to orchestrate the work instead of calling the underlying CLI tools manually.

## Worked Example C: Matrix Workflow for Repeated All-vs-All Comparison

This is the best route when you expect to compare many samples repeatedly against the same reference set.

### When to use this route

Use this when:

- you already have standard profile parquets
- you want resumable all-vs-all compare runs
- you expect to add new samples over time
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

This route turns standard per-sample profile parquets into a reusable comparison substrate. It is the best example for large or growing cohorts where repeated comparison cost matters more than one-time matrix construction.

## Optional: Add a New Sample Later and Run Only the Remaining Pairs

This is the practical “my cohort grew by one sample” scenario.
The two methods handle it differently.

### Standard Method

The standard method does not maintain a dedicated resumable compare database in the same way the matrix route does. Instead, you:

1. add the new sample profile to the profile DB
2. point the new comparison config at the previous comparison table
3. regenerate only the remaining pairs
4. rerun compare on the updated config

#### 1. Update the profile DB input table

Extend your `profiles.csv` to include the new sample:

```csv
sample_name,profile_loc
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

#### 2. Rebuild the standard compare config with the current comparison table

Point `--current-comp-table` at the existing genome comparison parquet from your previous standard run. The path below is an example placeholder; use the actual parquet produced by your earlier standard compare workflow:

```bash
zipstrain utilities build-genome-comparison-config \
  --profile-db outputs/profile_db.parquet \
  --gene-db-id tutorial_genes_v1 \
  --reference-genome-id tutorial_reference_v1 \
  --scope all \
  --min-cov 5 \
  --min-gene-compare-len 200 \
  --stb-file-loc reference/reference_genomes.stb \
  --current-comp-table outputs/standard_compare/genome_compare_results.parquet \
  --output-file outputs/genome_compare_config_updated.json
```

#### 3. Optionally inspect the remaining pairs

```bash
zipstrain utilities to-complete-table \
  --genome-comparison-object outputs/genome_compare_config_updated.json \
  --output-file outputs/remaining_pairs.csv
```

This gives you the not-yet-completed sample pairs implied by the updated config.

#### 4. Rerun standard compare

```bash
zipstrain compare genomes \
  --genome-comparison-object outputs/genome_compare_config_updated.json \
  --run-dir outputs/standard_compare_update \
  --calculate ani+ibs+identical_genes \
  --ani-method popani \
  --engine duckdb \
  --duckdb-threads 8
```

In other words, the standard route handles incremental updates by carrying forward the previous comparison table into a new config.

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

## Further Reading

- [CLI Reference](./cli.md)
- [Installation](./installation.md)
- [Nextflow Pipeline](./NextflowPipeline.md)
- [API Reference](./api.md)
