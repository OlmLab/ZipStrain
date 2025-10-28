# ZipStrain Command Line Interface

ZipStrain provides a comprehensive command-line interface for strain-level analysis of metagenomic data. The CLI is organized into several command groups for different functionalities.

## Installation

First, ensure ZipStrain is installed (see [Installation.md](Installation.md) for details).

## Usage

The main CLI command is `zipstrain` with various subcommands organized into functional groups:

```bash
zipstrain [GROUP] [COMMAND] [OPTIONS]
```

## Command Groups

### 1. Utilities (`utilities`)

General utility commands for data processing and preparation.

#### Build Null Model

Create a null model for sequencing errors based on statistical distributions:

```bash
zipstrain utilities build-null-model [OPTIONS]
```

**Options:**
- `-e, --error-rate FLOAT`: Error rate for the sequencing technology (default: 0.001)
- `-m, --max-total-reads INTEGER`: Maximum coverage to consider for a base (default: 10000)
- `-p, --p-threshold FLOAT`: Significance threshold for the Poisson distribution (default: 0.05)
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]
- `-t, --model-type [poisson]`: Type of null model to build (default: poisson)

**Example:**
```bash
zipstrain utilities build-null-model -e 0.001 -m 10000 -p 0.05 -o null_model.parquet
```

#### Merge Parquet Files

Merge multiple Parquet files into a single file:

```bash
zipstrain utilities merge_parquet -i INPUT_DIR -o OUTPUT_FILE
```

**Options:**
- `-i, --input-dir TEXT`: Directory containing Parquet files to merge [required]
- `-o, --output-file TEXT`: Path to save the merged Parquet file [required]

#### Process Mpileup

Process mpileup files and save results in Parquet format:

```bash
zipstrain utilities process_mpileup [OPTIONS]
```

**Options:**
- `-g, --gene-range-table-loc TEXT`: Location of the gene range table in TSV format [required]
- `-b, --batch-bed TEXT`: Location of the batch BED file [required]
- `-s, --batch-size INTEGER`: Buffer size for processing stdin from samtools (default: 10000)
- `-o, --output-file TEXT`: Location to save the output Parquet file [required]

#### Make BED File

Create a BED file from a FASTA database:

```bash
zipstrain utilities make_bed -d DB_FASTA_DIR -o OUTPUT_FILE
```

**Options:**
- `-d, --db-fasta-dir TEXT`: Path to the database in FASTA format [required]
- `-m, --max-scaffold-length INTEGER`: Maximum scaffold length to split into multiple entries (default: 500000)
- `-o, --output-file TEXT`: Path to save the output BED file [required]

#### Get Genome Lengths

Extract genome length information from scaffold-to-genome mapping:

```bash
zipstrain utilities get_genome_lengths -s STB_FILE -b BED_FILE -o OUTPUT_FILE
```

**Options:**
- `-s, --stb-file TEXT`: Path to the scaffold-to-genome mapping file [required]
- `-b, --bed-file TEXT`: Path to the BED file [required]
- `-o, --output-file TEXT`: Path to save the output Parquet file [required]

#### Genome Breadth Matrix

Generate a genome breadth matrix from profiles:

```bash
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

```bash
zipstrain utilities collect_breadth_tables -d BREADTH_TABLES_DIR -o OUTPUT_FILE
```

**Options:**
- `-d, --breadth-tables-dir TEXT`: Directory containing breadth tables in Parquet format [required]
- `-e, --extension TEXT`: File extension of the breadth tables (default: parquet)
- `-o, --output-file TEXT`: Path to save the collected breadth tables [required]

### 2. Gene Tools (`gene_tools`)

Commands for working with gene annotations and locations.

#### Gene Range Table

Build a gene range table from gene files:

```bash
zipstrain gene_tools gene-range-table -g GENE_FILE -l GENE_LIST -o OUTPUT_FILE
```

**Options:**
- `-g, --gene-file TEXT`: Location of gene file (Prodigal's nucleotide FASTA output) [required]
- `-l, --gene-list TEXT`: Location of gene list (text file with gene names) [required]
- `-o, --output-file TEXT`: Location to save output TSV file [required]

#### Gene Location Table

Build a gene location table for specified scaffolds:

```bash
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

```bash
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

### 4. Run (`run`)

Commands for running large-scale analyses and managing workflows.

#### Compare Genomes

Run genome comparisons in batches:

```bash
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

```bash
zipstrain run build-comp-database [OPTIONS]
```

**Options:**
- `-p, --profile-db-dir TEXT`: Directory containing profile in Parquet format [required]
- `-c, --config-file TEXT`: Path to the genome comparison database config file in JSON format [required]
- `-o, --output-dir TEXT`: Directory to genome comparison database object [required]
- `-f, --comp-db-file TEXT`: The initial database file (optional)

#### To Complete Table

Generate a table of completed genome comparison pairs:

```bash
zipstrain run to-complete-table -g GENOME_COMPARISON_OBJECT -o OUTPUT_FILE
```

**Options:**
- `-g, --genome-comparison-object TEXT`: Path to the genome comparison object in JSON format [required]
- `-o, --output-file TEXT`: Path to save the completed pairs Parquet file [required]

## Example Workflows

### Basic Strain Analysis Workflow

1. **Prepare null model:**
```bash
zipstrain utilities build-null-model -e 0.001 -m 10000 -o null_model.parquet
```

2. **Create BED file from database:**
```bash
zipstrain utilities make_bed -d database.fasta -o database.bed
```

3. **Process mpileup data:**
```bash
zipstrain utilities process_mpileup -g gene_ranges.tsv -b batch.bed -o processed.parquet
```

4. **Compare genomes:**
```bash
zipstrain compare single_compare_genome -m1 sample1.parquet -m2 sample2.parquet \
  -s1 scaffolds1.txt -s2 scaffolds2.txt -n null_model.parquet \
  -s scaffold_to_genome.tsv -o comparison_results.parquet
```

### Large-scale Analysis Workflow

1. **Build comparison database:**
```bash
zipstrain run build-comp-database -p profiles_dir -c config.json -o comp_db
```

2. **Run batch comparisons:**
```bash
zipstrain run compare_genomes -g comp_db/comparison_object.json -r run_results \
  -m 10 -e slurm -s slurm_config.json -c apptainer
```

3. **Generate results table:**
```bash
zipstrain run to-complete-table -g comp_db/comparison_object.json -o final_results.parquet
```

## File Formats

- **Parquet**: Primary format for structured data storage and processing
- **TSV**: Tab-separated values for simple tabular data
- **BED**: Standard genomic interval format
- **JSON**: Configuration files and object serialization
- **FASTA**: Sequence data input

## Performance Considerations

- Use `--memory-mode light` for large datasets with limited memory
- Choose appropriate `--engine` based on available hardware (GPU vs CPU)
- Adjust `--chrom-batch-size` based on available memory
- Use container engines for consistent environments across platforms

## Getting Help

For detailed help on any command:
```bash
zipstrain [GROUP] [COMMAND] --help
```

For general help:
```bash
zipstrain --help
```