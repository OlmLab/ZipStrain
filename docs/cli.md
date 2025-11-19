# ZipStrain Command Line Interface

ZipStrain provides a comprehensive command-line interface for strain-level analysis of metagenomic data. The CLI is organized into several command groups for different functionalities.

## Installation

First, ensure ZipStrain is installed (see [Installation.md](installation.md) for details).

## Usage

The main CLI command is `zipstrain` with various subcommands organized into functional groups:

```
zipstrain [GROUP] [COMMAND] [OPTIONS]
```

## Command Groups

### 1. Utilities (`utilities`)

General utility commands for data processing and preparation.

#### Build Null Model

Create a null model for sequencing errors based on statistical distributions:

```
zipstrain utilities build-null-model [OPTIONS]
```

**Options:**

- `-e, --error-rate FLOAT`: Error rate for the sequencing technology (default: 0.001)
- `-m, --max-total-reads INTEGER`: Maximum coverage to consider for a base (default: 10000)
- `-p, --p-threshold FLOAT`: Significance threshold for the Poisson distribution (default: 0.05)
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]
- `-t, --model-type [poisson]`: Type of null model to build (default: poisson)

**Example:**

```
zipstrain utilities build-null-model -e 0.001 -m 10000 -p 0.05 -o null_model.parquet
```

#### Merge Parquet Files

Merge multiple Parquet files into a single file:

```
zipstrain utilities merge_parquet -i INPUT_DIR -o OUTPUT_FILE
```

**Options:**

- `-i, --input-dir TEXT`: Directory containing Parquet files to merge [required]
- `-o, --output-file TEXT`: Path to save the merged Parquet file [required]

#### Process Mpileup

Process mpileup files and save results in Parquet format:

```
zipstrain utilities process_mpileup [OPTIONS]
```

**Options:**

- `-g, --gene-range-table-loc TEXT`: Location of the gene range table in TSV format [required]
- `-b, --batch-bed TEXT`: Location of the batch BED file [required]
- `-s, --batch-size INTEGER`: Buffer size for processing stdin from samtools (default: 10000)
- `-o, --output-file TEXT`: Location to save the output Parquet file [required]

#### Make BED File

Create a BED file from a FASTA database:

```
zipstrain utilities make_bed -d DB_FASTA_DIR -o OUTPUT_FILE
```

**Options:**

- `-d, --db-fasta-dir TEXT`: Path to the database in FASTA format [required]
- `-m, --max-scaffold-length INTEGER`: Maximum scaffold length to split into multiple entries (default: 500000)
- `-o, --output-file TEXT`: Path to save the output BED file [required]

#### Get Genome Lengths

Extract genome length information from scaffold-to-genome mapping:

```
zipstrain utilities get_genome_lengths -s STB_FILE -b BED_FILE -o OUTPUT_FILE
```

**Options:**

- `-s, --stb-file TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-b, --bed-file TEXT`: Path to the BED file [required]
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]

#### Genome Breadth Matrix

Generate a genome breadth matrix from profiles:

```
zipstrain utilities genome_breadth_matrix [OPTIONS]
```

**Options:**

- `-p, --profile TEXT`: Path to the profile Parquet file [required]
- `-g, --genome-length TEXT`: Path to the genome length Parquet file [required]
- `-s, --stb TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-c, --min-cov INTEGER`: Minimum coverage to consider a position (default: 1)
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]

#### Collect Breadth Tables

Collect multiple genome breadth tables into a single file:

```
zipstrain utilities collect_breadth_tables -d BREADTH_TABLES_DIR -o OUTPUT_FILE
```

**Options:**

- `-d, --breadth-tables-dir TEXT`: Directory containing breadth tables in Parquet format [required]
- `-e, --extension TEXT`: File extension of the breadth tables (default: parquet)
- `-o, --output-file TEXT`: Path to save the collected breadth tables [required]

#### Strain Heterogeneity

Calculate strain heterogeneity for each genome based on nucleotide frequencies:

```
zipstrain utilities strain_heterogeneity [OPTIONS]
```

**Options:**

- `-p, --profile-file TEXT`: Path to the profile Parquet file [required]
- `-s, --stb-file TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-c, --min-cov INTEGER`: Minimum coverage to consider a position (default: 5)
- `-f, --freq-threshold FLOAT`: Frequency threshold to define dominant nucleotide (default: 0.8)
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]

#### Build Profile Database

Build a profile database from a CSV file:

```
zipstrain utilities build-profile-db -p PROFILE_DB_CSV -o OUTPUT_FILE
```

**Options:**

- `-p, --profile-db-csv TEXT`: Path to the profile database CSV file [required]
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]

#### Build Comparison Configuration

Build a comparison configuration JSON file:

```
zipstrain utilities build-comparison-config [OPTIONS]
```

**Options:**

- `-p, --profile-db TEXT`: Path to the profile database Parquet file [required]
- `-g, --gene-db-id TEXT`: Gene database ID [required]
- `-r, --reference-db-id TEXT`: Reference fasta ID [required]
- `-s, --scope TEXT`: Genome scope for comparison (default: "all")
- `-c, --min-cov INTEGER`: Minimum coverage to consider a position (default: 5)
- `-l, --min-gene-compare-len INTEGER`: Minimum gene length to consider for comparison (default: 200)
- `-n, --null-model-p-value FLOAT`: P-value threshold for the null model (default: 0.05)
- `-t, --stb-file-loc TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-m, --null-model-loc TEXT`: Path to the null model Parquet file [required]
- `-a, --current-comp-table TEXT`: Path to existing comparison table in Parquet format (optional)
- `-o, --output-file TEXT`: Path to save the output configuration JSON file [required]

#### To Complete Table

Generate a table of remaining genome comparison pairs:

```
zipstrain utilities to-complete-table -g GENOME_COMPARISON_OBJECT -o OUTPUT_FILE
```

**Options:**

- `-g, --genome-comparison-object TEXT`: Path to the genome comparison object in JSON format [required]
- `-o, --output-file TEXT`: Path to save the completed pairs CSV file [required]

### 2. Gene Tools (`gene_tools`)

Commands for working with gene annotations and locations.

#### Gene Range Table

Build a gene range table from gene files:

```
zipstrain gene_tools gene-range-table -g GENE_FILE -o OUTPUT_FILE
```

**Options:**

- `-g, --gene-file TEXT`: Location of gene file (Prodigal's nucleotide FASTA output) [required]
- `-o, --output-file TEXT`: Location to save output TSV file [required]

#### Gene Location Table

Build a gene location table for specified scaffolds:

```
zipstrain gene_tools gene-loc-table -g GENE_FILE -s SCAFFOLD_LIST -o OUTPUT_FILE
```

**Options:**

- `-g, --gene-file TEXT`: Location of gene file (Prodigal's nucleotide FASTA output) [required]
- `-s, --scaffold-list TEXT`: Location of scaffold list (text file with scaffold names) [required]
- `-o, --output-file TEXT`: Location to save output Parquet file [required]

### 3. Compare (`compare`)

Commands for comparing genomes and samples.

#### Single Genome Compare

Compare two mpileup files for genome analysis:

```
zipstrain compare single_compare_genome [OPTIONS]
```

**Options:**

- `-m1, --mpileup-contig-1 TEXT`: Path to the first mpileup file [required]
- `-m2, --mpileup-contig-2 TEXT`: Path to the second mpileup file [required]
- `-s1, --scaffolds-1 TEXT`: Path to the list of scaffolds for the first mpileup file [required]
- `-s2, --scaffolds-2 TEXT`: Path to the list of scaffolds for the second mpileup file [required]
- `-n, --null-model TEXT`: Path to the null model Parquet file [required]
- `-s, --stb-file TEXT`: Path to the scaffold to genome mapping file [required]
- `-c, --min-cov INTEGER`: Minimum coverage to consider a position (default: 5)
- `-l, --min-gene-compare-len INTEGER`: Minimum gene length to consider for comparison (default: 100)
- `-m, --memory-mode [heavy|light]`: Memory mode for processing (default: heavy)
- `-b, --chrom-batch-size INTEGER`: Batch size for processing chromosomes (default: 10000)
- `-o, --output-file TEXT`: Path to save the parquet file [required]
- `-g, --genome TEXT`: If provided, do comparison only for the specified genome (default: all)
- `-e, --engine [streaming|gpu|auto]`: Engine to use for processing (default: streaming)

### 4. Profile (`profile`)

Commands for profiling BAM files.

#### Prepare Profiling

Prepare files needed for profiling BAM files:

```
zipstrain profile prepare_profiling [OPTIONS]
```

**Options:**

- `-r, --reference-fasta TEXT`: Path to the reference genome in FASTA format [required]
- `-g, --gene-fasta TEXT`: Path to the gene annotations in FASTA format [required]
- `-s, --stb-file TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-o, --output-dir TEXT`: Directory to save the profiling database [required]

#### Profile Single

Profile a single BAM file:

```
zipstrain profile profile-single [OPTIONS]
```

**Options:**

- `-b, --bed-file TEXT`: Path to the BED file describing regions to be profiled [required]
- `-a, --bam-file TEXT`: Path to the BAM file to be profiled [required]
- `-g, --gene-range-table TEXT`: Path to the gene range table [required]
- `-n, --num-workers INTEGER`: Number of workers to use for profiling (default: 1)
- `-o, --output-dir TEXT`: Directory to save the profiling output [required]

### 5. Run (`run`)

Commands for running large-scale analyses and managing workflows.

#### Profile

Run BAM file profiling in batches:

```
zipstrain run profile [OPTIONS]
```

**Options:**

- `-i, --input-table TEXT`: Path to the input table in TSV format containing sample names and BAM file paths [required]
- `-s, --stb-file TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-g, --gene-range-table TEXT`: Path to the gene range table file [required]
- `-b, --bed-file TEXT`: Path to the BED file for profiling regions [required]
- `-l, --genome-length-file TEXT`: Path to the genome length file [required]
- `-r, --run-dir TEXT`: Directory to save the run data [required]
- `-n, --num-procs INTEGER`: Number of processors to use for each profiling task (default: 8)
- `-m, --max-concurrent-batches INTEGER`: Maximum number of concurrent batches to run (default: 5)
- `-p, --poll-interval INTEGER`: Polling interval in seconds to check batch status (default: 1)
- `-e, --execution-mode TEXT`: Execution mode: 'local' or 'slurm' (default: local)
- `-c, --slurm-config TEXT`: Path to the SLURM configuration file in JSON format (required if execution mode is 'slurm')
- `-o, --container-engine TEXT`: Container engine to use: 'local', 'docker' or 'apptainer' (default: local)
- `-t, --task-per-batch INTEGER`: Number of tasks to include in each batch (default: 10)

#### Compare Genomes

Run genome comparisons in batches:

```
zipstrain run compare_genomes [OPTIONS]
```

**Options:**

- `-g, --genome-comparison-object TEXT`: Path to the genome comparison object in JSON format [required]
- `-r, --run-dir TEXT`: Directory to save the run data [required]
- `-m, --max-concurrent-batches INTEGER`: Maximum number of concurrent batches to run (default: 5)
- `-p, --poll-interval INTEGER`: Polling interval in seconds to check batch status (default: 1)
- `-e, --execution-mode TEXT`: Execution mode: 'local' or 'slurm' (default: local)
- `-s, --slurm-config TEXT`: Path to the SLURM configuration file in JSON format (required if execution mode is 'slurm')
- `-c, --container-engine TEXT`: Container engine to use: 'local', 'docker' or 'apptainer' (default: local)
- `-t, --task-per-batch INTEGER`: Number of tasks to include in each batch (default: 10)
- `-a, --polars-engine [streaming|gpu|auto]`: Polars engine to use (default: streaming)
- `-b, --chrom-batch-size INTEGER`: Batch size for processing chromosomes (default: 10000)
- `-h, --memory-mode [heavy|light]`: Memory mode for processing (default: heavy)

#### Build Comparison Database

Build a genome comparison database from profiles:

```
zipstrain run build-comp-database [OPTIONS]
```

**Options:**

- `-p, --profile-db-dir TEXT`: Directory containing profile in Parquet format [required]
- `-c, --config-file TEXT`: Path to the genome comparison database config file in JSON format [required]
- `-o, --output-dir TEXT`: Directory to save genome comparison database object [required]
- `-f, --comp-db-file TEXT`: Path to initial database file (optional)

## Test Command

### Test Installation

Verify ZipStrain setup and dependencies:

```
zipstrain test
```

This command checks if all required dependencies (like samtools) are properly installed and accessible.

## Example Workflows

### Basic Profiling and Analysis Workflow

1. **Prepare profiling files:**

```
zipstrain profile prepare_profiling \
  -r reference.fasta \
  -g genes.fasta \
  -s scaffold_to_genome.tsv \
  -o profiling_db
```

2. **Build null model:**

```
zipstrain utilities build-null-model \
  -e 0.001 -m 10000 -p 0.05 \
  -o null_model.parquet
```

3. **Profile single BAM file:**

```
zipstrain profile profile-single \
  -b profiling_db/genomes_bed_file.bed \
  -a sample1.bam \
  -g profiling_db/gene_range_table.tsv \
  -n 4 \
  -o sample1_profile
```

4. **Compare two profiles:**

```
zipstrain compare single_compare_genome \
  -m1 sample1_profile/sample1.parquet \
  -m2 sample2_profile/sample2.parquet \
  -s1 sample1_profile/sample1.parquet.scaffolds \
  -s2 sample2_profile/sample2.parquet.scaffolds \
  -n null_model.parquet \
  -s scaffold_to_genome.tsv \
  -o comparison_results.parquet
```

### Large-scale Analysis Workflow

1. **Build profile database:**

```
zipstrain utilities build-profile-db \
  -p profiles.csv \
  -o profiles.parquet
```

2. **Build comparison configuration:**

```
zipstrain utilities build-comparison-config \
  -p profiles.parquet \
  -g genes_v1 \
  -r gtdb_r214 \
  -s all \
  -c 5 \
  -l 200 \
  -t scaffold_to_genome.tsv \
  -m null_model.parquet \
  -o comparison_config.json
```

3. **Build comparison database:**

```
zipstrain run build-comp-database \
  -p profiles.parquet \
  -c comparison_config.json \
  -o comparison_db.json
```

4. **Run batch comparisons:**

```
zipstrain run compare_genomes \
  -g comparison_db.json \
  -r run_results \
  -m 10 \
  -e slurm \
  -s slurm_config.json \
  -c apptainer
```

### Strain Analysis Workflow

1. **Calculate strain heterogeneity:**

```
zipstrain utilities strain_heterogeneity \
  -p sample_profile.parquet \
  -s scaffold_to_genome.tsv \
  -c 5 \
  -f 0.8 \
  -o strain_het.parquet
```

2. **Generate genome breadth matrix:**

```
zipstrain utilities genome_breadth_matrix \
  -p sample_profile.parquet \
  -g genome_lengths.parquet \
  -s scaffold_to_genome.tsv \
  -c 1 \
  -o breadth_matrix.parquet
```

## File Formats

- **Parquet**: Primary format for structured data storage and processing
- **TSV**: Tab-separated values for simple tabular data
- **BED**: Standard genomic interval format
- **JSON**: Configuration files and object serialization
- **FASTA**: Sequence data input
- **CSV**: Simple data tables for input/output

## Performance Considerations

- Use `--memory-mode light` for large datasets with limited memory
- Choose appropriate `--engine` based on available hardware (GPU vs CPU)
- Adjust `--chrom-batch-size` based on available memory
- Use container engines for consistent environments across platforms
- Consider SLURM execution mode for large-scale HPC analyses

## Getting Help

For detailed help on any command:

```
zipstrain [GROUP] [COMMAND] --help
```

For general help:

```
zipstrain --help
```

To test your installation:

```
zipstrain test
```
