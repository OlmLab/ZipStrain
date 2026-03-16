# Nextflow Pipeline for ZipStrain

This page reflects the current `zipstrain.nf` workflow in this repository.

## What the Pipeline Supports

- Read mapping with Bowtie2 (`map_reads`)
- Profile generation from BAM files (`profile`)
- End-to-end SRA to profile (`from_sra_to_profile`)
- Pairwise genome comparison across profiles (`compare_genomes`)
- Pairwise gene comparison across profiles (`compare_genes`)

## Running Pattern

```bash
nextflow run zipstrain.nf \
  --mode <mode> \
  --input_table <path/to/input.csv> \
  --output_dir <path/to/output_dir> \
  -c conf.config \
  -profile <docker|alpine|gutbot|blanca> \
  -resume
```

`conf.config` already defines resources for the current process set and includes example execution profiles.

## Key Pipeline Parameters

- `--mode`: `map_reads`, `profile`, `from_sra_to_profile`, `compare_genomes`, `compare_genes`
- `--input_type`: depends on mode (`local`, `sra`, `profile_table`, `pair_table`)
- `--parallel_mode`: `single` or `batched` for comparison workflows
- `--batch_size`: number of pairs per batch when `--parallel_mode batched`
- `--batch_compare_n_parallel`: parallel jobs inside each batched comparison task
- `--compare_genome_scope`: genome scope for genome comparisons (`all` or genome ID)
- `--compare_gene_scope`: gene scope for gene comparisons (`all:all`, `<genome>:all`, `all:<gene>`, `<genome>:<gene>`)
- `--compare_duckdb_memory_limit`: forwarded to single compare commands

## 1) Map Reads (`mode=map_reads`)

### Input Table (`--input_type local`)

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

### Input Table (`--input_type sra`)

```csv
Run
SRR12345678
SRR12345679
```

### A) Use Existing Reference Genome

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

### B) Build Reference from Sylph Automatically

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

### Map Outputs

- BAM files: `<output_dir>/*.bam`
- Sylph tables: `<output_dir>/sylph_abundance/`
- Built reference bundle (when auto-built): `<output_dir>/db_from_sylph/`

## 2) Generate Profiles from BAM (`mode=profile`)

### Input Table

```csv
sample_name,bamfile
S1,/data/S1.bam
S2,/data/S2.bam
```

### Command

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

### Profile Outputs

- `<output_dir>/*_profile.parquet`
- `<output_dir>/*_genome_stats.parquet`
- `<output_dir>/*_gene_stats.parquet`

## 3) End-to-End SRA to Profile (`mode=from_sra_to_profile`)

### Input Table

```csv
Run
SRR12345678
SRR12345679
```

### Command

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

### Outputs

- `<output_dir>/profiles/*_profile.parquet`
- `<output_dir>/profiles/*_genome_stats.parquet`
- `<output_dir>/profiles/*_gene_stats.parquet`

## 4) Compare Genomes (`mode=compare_genomes`)

### Input Option A: All-vs-All from Profile List (`--input_type profile_table`)

```csv
sample_names,mpileup_files
S1,/profiles/S1_profile.parquet
S2,/profiles/S2_profile.parquet
S3,/profiles/S3_profile.parquet
```

### Input Option B: Explicit Pairs (`--input_type pair_table`)

```csv
sample_name_1,sample_name_2,profile_location_1,profile_location_2
S1,S2,/profiles/S1_profile.parquet,/profiles/S2_profile.parquet
S1,S3,/profiles/S1_profile.parquet,/profiles/S3_profile.parquet
```

### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genomes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_genome_scope all \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --compare_duckdb_memory_limit 4GB \
  --output_dir out_compare_genomes \
  -c conf.config \
  -profile docker \
  -resume
```

## 5) Compare Genes (`mode=compare_genes`)

Input-table formats are the same as genome compare (`profile_table` or `pair_table`).

### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_gene_scope all:all \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --compare_duckdb_memory_limit 4GB \
  --output_dir out_compare_genes \
  -c conf.config \
  -profile docker \
  -resume
```

## Comparison Outputs

- Final merged table: `<output_dir>/merged_comparisons.parquet`
- Intermediate batched outputs (when `parallel_mode=batched`): `<output_dir>/batch_comparisons/`

## Important Notes

- The old `--compare_memory_mode` and `--compare_chrom_batch_size` parameters are not part of the current `zipstrain.nf`.
- The pipeline currently forwards DuckDB memory limit via `--compare_duckdb_memory_limit` but does not expose compare engine/threads as Nextflow params in this script.
- For auto-built references, genome selection comes from the merged Sylph abundance table through `zipstrain utilities build-genome-db`.
