# Tutorial
This tutorial will guide you through the basics of using ZipStrain.

!!! warning "Work in progress"
    The following sections are still being completed:

    - [ ] Download links for example input/output files throughout the tutorial
    - [ ] Step 6: Comparing profiled samples at the genome level
    - [ ] Step 7: Comparing profiled samples at the gene level
    - [ ] Code examples for strain sharing and downstream analyses
    - [ ] Completing the "Recommended Workflow" section

## Introduction to ZipStrain
ZipStrain is a tool for profiling a metagenomics sample against a database of reference genomes and performing comparisons between the profiles as well as downstream analyses. 

The typical workflow consists of the following steps:

1- Mapping reads to a reference database using any read mapper of your choice (e.g., BWA, Bowtie2, Minimap2).

2- Generating a profile from the mapped reads using ZipStrain. This step will generate a parquet file containing nucleotide frequencies for all positions in the reference genomes.

3- Comparing the generated profiles against each other and/or against a reference profile.

4- Performing downstream analyses, such as Identity-By-State (IBS), StrainSharing, and clustering.

### Mapping Reads

In this step, you will map your metagenomics reads to a reference database using a read mapper of your choice. The reference database is a concatenation of all reference genomes you want to profile against. The output of this step should be a BAM file. It is recommended to discard unmapped reads to reduce file size. Currently, you can use the nextflow pipeline accompanying ZipStrain to perform this step with Bowtie2. Here is an example command using Bowtie2, but for more information about the nextflow pipeline, please refer to the [Nextflow Pipeline Documentation](./NextflowPipeline.md).

```

nextflow run zipstrain.nf --mode 'map_reads' --input_type 'sra|local' --input_table 'path/to/your/input_table.tsv' --reference_genome <path/to/your/reference_db.fasta> --output_dir <path/to/your/output_dir> -c conf.config -profile <your_profile> -resume

```

NOTE: map_reads mode can work with either SRA accessions or local FASTQ files. Please refer to the Nextflow Pipeline Documentation for more details.

This step will generate BAM files for each of your samples in the specified output directory.

### Generating Profiles

In this step, you will generate profiles from the mapped reads using ZipStrain. A profile is simply a table in parquet format that contains the 
following columns:

```

|chrom|genome|gene|pos|A|C|G|T|

```

Where `chrom` is the scaffold, `genome` is the genome name (derived from the STB file), `gene` is the gene name (if a gene file is provided), `pos` is the position in the reference genome, and `A`, `C`, `G`, and `T` are the counts of each nucleotide at that position based on the mapped reads.

The information will be generated for every position in the reference genomes that has at least one read mapped to it. There are two ways to generate profiles using ZipStrain: 

- ZipStrain as a command-line tool

- ZipStrain Nextflow pipeline.

You can profile multiple BAM files using either the ZipStrain command-line interface (CLI) or Nextflow Pipeline. Below are examples of both methods.
For both you need to prepare a csv file containing the paths to the bam files to be profiled and the sample names. Example of such a csv file:


```csv
sample_name,bamfile
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
```

#### ZipStrain CLI

To profile multiple BAM files using the ZipStrain CLI, you should first prepare some files:

```bash
zipstrain utilities prepare_profiling  --reference-fasta <path/to/reference/fasta> --gene-fasta <path/to/reference/fasta/genes> --stb-file  <path/to/stb/file> --output-dir <directory/to/save/outputs>
```

Your output directory should contain the following files:

-   genomes_bed_file.bed
-   genome_lengths.parquet
-   gene_range_table.tsv

You also need to build a null model for sequencing error adjustment:

```bash
zipstrain utilities build-null-model --bam-file <path/to/any/sample.bam> --output-file <path/to/null_model.parquet>
```

Now you can profile your bam files:

```bash
zipstrain profile --input-table <path/to/bam/csv> --stb-file <path/to/stb/file> --null-model <path/to/null_model.parquet> --gene-range-table <path/to/gene/range> --bed-file <path/to/bed/file> --genome-length-file <path/to/genome_lengths.parquet> --run-dir <path/to/save/generated/files>
```


#### Nextflow Pipeline

To profile multiple BAM files using Nextflow, you can create a Nextflow script as follows:

```

nextflow run zipstrain.nf --mode "profile" --input_table <path/to/bam/csv>  --gene_file <path/to/reference/fasta/genes> --stb <path/to/stb/file>  --output_dir <path/to/save/generated/files> --reference_genome <path/to/reference/fasta> -c conf.config -profile <your/system/specific/profile> -resume

```

**Note**  With the nextflow pipeline, you don't need the preparation step and those will be made along the way.

**Note** `profile` mode requires a gene file and a STB file. The gene file MUST BE following Prodigal NUCLEOTIDE format. For generating STB file, use the following command:

```zipstrain utilities generate_stb --genomes-dir-file <path/to/genomes_dir_file>```

If you have trouble generating the STB file, you can simply make one using a custom script. The STB file is a tab-separated file with two columns: scaffold name and genome file name. 

### Compare genomes in multiple profiled samples 

You can compare multiple profiled samples using either the ZipStrain command-line interface (CLI) or Nextflow Pipeline. Below are examples of both methods. The two approaches are slightly different in terms of input requirements. First let's see how this is done using the ZipStrain CLI.

#### ZipStrain CLI
To compare multiple profiled samples using the ZipStrain CLI, first you need to build a profile database that contains all the profiled samples you want to compare. You can do this by running the following command:

```bash
zipstrain utilities build-profile-db --profile-db-csv <path/to/profiles/csv> --output-file <path/to/save/profile/db>
```
The input CSV file should have the following columns:

    
    - profile_name: An arbitrary name given to the profile (Usually sample name or name of the parquet file)
    
    - profile_location: The location of the profile
    
    - reference_db_id: The ID of the reference database. This could be the name or any other identifier for the database that the reads are mapped to.
    
    - gene_db_id: The ID of the gene database in fasta format. This could be the name or any other identifier for the database that the reads are mapped to.

Running this command will perform the necessary checks and if successful, it will create a profile database in parquet format at the specified output location.

Next, You use this profile database build a configuration json file that will be used to calculate the pairs that need to be compared. In this step you need to define the required parameters for comparison such as min_coverage, etc. You can make the configuration file by running the following command:

```bash
zipstrain utilities build-genome-comparison-config \
--profile-db <path/to/profile/db> \
--gene-db-id <gene_db_id_used_in_profile_db> \
--reference-genome-id <reference_db_id_used_in_profile_db> \
--scope "all" \
--min-cov 5 \
--min-gene-compare-len 200 \
--stb-file-loc <path/to/stb/file> \
--current-comp-table <path/to/current/comparison/table.parquet> \
--output-file <path/to/save/comparison/config.json>
```
Note that providing current-comp-table is optional. If provided, the comparison config will only include pairs that are not already compared in the current comparison table.

Finally, you can run the comparison using the generated configuration file and the profile database:

```bash
zipstrain compare genomes \
--genome-comparison-object <path/to/comparison/config.json> \
--run-dir <path/to/save/comparison/outputs> \
--calculate ani \
--ani-method popani \
--engine duckdb \
--calculate ani+ibs+identical_genes \
--max-concurrent-batches 1 \
--duckdb-threads 8
```

`single_compare_genome` supports `--engine polars|duckdb` (default: `polars`) and metric selection via `--calculate` (default: `all`). For lower-memory machines, set DuckDB's memory limit, and for CPU control set DuckDB threads:

```bash
zipstrain utilities single_compare_genome \
--mpileup-contig-1 <profile_1.parquet> \
--mpileup-contig-2 <profile_2.parquet> \
--stb-file <path/to/stb.tsv> \
--calculate ani+ibs+identical_genes \
--engine duckdb \
--output-file <out.parquet> \
--duckdb-memory-limit 2GB \
--duckdb-threads 8 \
--duckdb-temp-directory /tmp
```

#### Nextflow Workflow

There are two ways you can run the comparison workflow in ZipStrain Nextflow pipeline. 

1- I have a comparison config file already made using the ZipStrain CLI as explained above.

In this case you obtain a table of remaining pairs to be compared from the comparison config file and run the comparison as follows:

```
zipstrain utilities to-complete-table --genome-comparison-object <path/to/comparison/config.json> --output-file <path/to/output/remaining_pairs.csv>
```

Then you can run the comparison using the following command:

```
nextflow run zipstrain.nf --mode compare_genomes \
 --input_table <path/to/output/remaining_pairs.csv> \
 --input_type "pair_table" --gene_file <path/to/gene/fasta/file> \
 --reference_genome <path/to/reference/genome.fasta> \  
 --stb <path/to/stb/file.stb> -c conf.config \
 --output_dir "<path/to/output/directory>" \
 --compare_genome_scope "all" \
 --compare_memory_mode "heavy" \
 --parallel_mode "batched" \
 --batch_size 2000 -profile <profile_name> \
 --batch_compare_n_parallel 3 -qs 200 -resume
```

2- I don't have a comparison config file and I want to run the comparison directly from the profile lists.

In this case, the nextflow pipeline will run all non-redundant pairwise comparisons between the provided profiles. Here is an example command:

```bash
nextflow run zipstrain.nf --mode compare_genomes \
 --input_table <path/to/profiles/csv> \
 --input_type "profile_table" --gene_file <path/to/gene/fasta/file> \
 --reference_genome <path/to/reference/genome.fasta> \
 --stb <path/to/stb/file.stb> -c conf.config \
 --output_dir "<path/to/output/directory>" \
 --compare_genome_scope "all" \
 --compare_memory_mode "heavy" \
 --parallel_mode "batched" \
 --batch_size 2000 -profile <profile_name> \
 --batch_compare_n_parallel 3 -qs 200 -resume
```

For more information, refer to the [Nextflow Pipeline Documentation](./NextflowPipeline.md).

#### Output files

Comparing profiles yields a results table containing all compared profile pairs. When using the ZipStrain CLI, the output is at `<run_dir>/Outputs/all_comparisons.parquet`; with Nextflow, it is `merged_comparisons.parquet`. The table has the following columns:

|genome|total_positions|share_allele_pos|genome_pop_ani|max_consecutive_length|shared_genes_count|identical_gene_count|perc_id_genes|sample_1|sample_2|
|-----|---------------|----------------|--------------|---------------------|------------------|--------------------|-------------|--------|--------|



### Downstream Analyses

In this step, you can use the perform statiscal analyses on the comparison results generated in the previous step and visualize them using ZipStrain's Python API. All of the functionalities in the visualization module are explained below:

#### Strain Sharing Analysis

Definition of strain sharing is somewhat losely defined in the literature. In ZipStrain, we use the following steps to define strainsharing:

-   For each genome in the reference fasta for which the profiling is done, popANI is calculated between each pair of samples in the comparison step.

-   If the popANI between two samples for a given genome is above a certain threshold (default: 99.9%), and the breadth of coverage for that genome in both samples is above a certain threshold (default: 0.5), then the two samples are considered to share a strain for that genome.

-  The strain sharing is checked for all genomes in the reference fasta.

-  Each sample is mapped to a group using a sample to population mapping file provided by the user (See the example below). 

-  Finally, the strain sharing between group A,B is calculated as the number of genome strains shared between samples in group A and samples in group B divided by the total number of genomes in group A. As a result the order of the groups matters here. The main justification for this definition is that in many cases this ordering is biologically relevant. For example, when looking at mother-infant pairs, if only 3 genomes are present in infants, but 100 genomes are present in mothers, and they share all 3 genomes in the infants, we would like to see that the infants have 100% of their strains shared with mothers, while mothers only have 3% of their strains shared with infants.

### Recommended Workflow

ZipStrain's Python API provides functionalities needed to create and organize your profiles and comparisons tables. It is highly recommended to use this workflow to manage your profiles and comparisons if you plan to use ZipStrain for multiple studies to utilize the full potential of ZipStrain. Here is a recommended workflow:

#### 1- Create a centralized database for your profiles

Database module in ZipStrain provides functionalities to create and manage a centralized database for your profiles. The easiest way to create a profile database is to provide a CSV file that has necessary information about each profile. Your CSV file should have the following columns (and only these columns):

- profile_name: An arbitrary name given to the profile (Usually sample name or name of the parquet file)
    
- profile_location: The location of the profile in parquet format. This is the main parquet file generated during the profiling step.

- reference_db_id: The ID of the reference database. This could be the name or any other identifier for the database that the reads are mapped to. Could be used for filtering profiles based on reference database later on.
    
- gene_db_id: The ID of the gene database in fasta format. This could be the name or any other identifier of your choice for the database that the reads are mapped to. Could be used for filtering profiles based on gene database later on.



## An Example Workflow

This example guides you through a complete workflow of using ZipStrain from building a genome database all the way to perform pairwise genome and gene comparisons and between multiple samples. Here we use MGnify Genomes mouse gut catalogue v1.0 as our reference genome database and use some example reads provided below to demonstrate the workflow. For each step, you can download the necessary input files from the provided links and follow the instructions to produce the outputs or if you want to skip a step, you can directly download the output files from the provided links.

### Step 1- Prepare Reference Database

|Inputs|Link|
|------|-----|
| MGnify Genomes mouse gut catalogue v1.0 metadata | [Link](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/mouse-gut/v1.0/genomes-all_metadata.tsv) |

|Outputs|Link|
|-------|-----|
| Concatenated reference genome fasta file | [Link](TOBEIMPLEMENTED) |
| STB file | [Link](TOBEIMPLEMENTED) |

We first need to extract the accessions of the species representatives genomes from the metadata file provided above. You can use any tool of your choice to extract the accessions. Here is an example command using polars library (installed when you install ZipStrain) in python:

```python
import polars as pl
metadata_df = pl.read_csv("genomes-all_metadata.tsv", separator="\t") #Load metadata table downloaded from the provided link
metadata_df.join(pl.read_csv("genomes-all_metadata.tsv",separator="\t").select('Species_rep',"FTP_download").unique("Species_rep"),left_on="Genome",right_on="Species_rep",how="inner").select("Genome").with_columns("https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/mouse-gut/v1.0/species_catalogue/"+
pl.col("Genome").str.slice(0,11)+"/"+pl.col("Genome")+"/genome/"+pl.col("Genome")+".fna").select("literal"
).write_csv("sp_rep_genomes.csv",include_header=False) # A one(ish) liner that generates the link of each species representative genome
```

This will a text file containing the link to download each genome. Using the tool of your choice, download all the genomes. Here is an example command using `wget`:

```bash
mkdir genomes
for link in $(cat genome_accessions.txt); do
    wget $link -P genomes/
done
```

Now you can use ZipStrain to build the STB file:

```bash
zipstrain utilities generate_stb -g genomes/ -o  mgnify_mouse_gut_genomes.stb --extension ".fna"
```

Finally, concatenate all the genomes into a single fasta file:

```bash
cat genomes/*.fna > mgnify_mouse_gut_genomes.fa 
```

This will be your reference genome database for profiling.

### Step 2- Find genes in the reference genomes by running Prodigal

|Inputs|Link|
|------|-----|
| MGnify Genomes mouse gut catalogue v1.0 concatenated reference genome fasta file | [Link](TOBEIMPLEMENTED) |

|Outputs|Link|
|-------|-----|
| Prodigal gene fasta file | [Link](TOBEIMPLEMENTED) |  

For this example, we can directly use the concatenated genome fasta file from Step 1. You can run Prodigal to find genes in the reference genomes using the following command:

```bash
prodigal -i mgnify_mouse_gut_genomes.fa -d mgnify_mouse_gut_genes.fasta  -p meta
```

### Step 3- Map example reads to the reference database

|Inputs|Link|
|------|-----|
| Example metagenomics reads | [Link](TOBEIMPLEMENTED) |
| MGnify Genomes mouse gut catalogue v1.0 concatenated reference genome fasta file | [Link](TOBEIMPLEMENTED) |

|Outputs|Link|
|-------|-----|
| Mapped BAM files | [Link](TOBEIMPLEMENTED) |

You can use any read mapper of your choice to map the reads to the reference database. Here we use the bowtie2 option provided in the ZipStrain Nextflow pipeline to perform the mapping. First, prepare an input CSV file containing the sample names and paths to the FASTQ files:

```csv
sample_name,reads1,reads2
sample1,/path/to/sample1_R1.fastq,/path/to/sample1_R2.fastq
sample2,/path/to/sample2_R1.fastq,/path/to/sample2_R2.fastq
```

Then run the following Nextflow command to perform the mapping:

```bash
nextflow run zipstrain.nf --mode 'map_reads' --input_type 'local' --input_table 'path/to/your/input_table.csv' --reference_genome mgnify_mouse_gut_genomes.fa --output_dir mapping_output/ -c conf.config -profile <your_profile> -resume
```

This will generate BAM files for each sample in the `mapping_output/` directory.

### Step 4- Prepare necessary files for profiling

|Inputs|Link|
|------|-----|
| MGnify Genomes mouse gut catalogue v1.0 concatenated reference genome fasta file | [Link](TOBEIMPLEMENTED) |
| Prodigal gene fasta file | [Link](TOBEIMPLEMENTED) |
| STB file | [Link](TOBEIMPLEMENTED)|

|Outputs|Link|
|-------|-----|
| Bed file | [Link](TOBEIMPLEMENTED) |
| Genome lengths parquet file | [Link](TOBEIMPLEMENTED) |
| Gene range table TSV file | [Link](TOBEIMPLEMENTED) |

You can use ZipStrain to prepare the necessary files for profiling using the following command:

```bash
zipstrain utilities prepare_profiling -r mgnify_mouse_gut_genomes.fa -g mgnify_mouse_gut_genes.fasta -s mgnify_mouse_gut_genomes.stb  -o preprofiles
``` 

This will generate the following files in the `preprofiles/` directory:

- genomes_bed_file.bed
- genome_lengths.parquet
- gene_range_table.tsv

### Step 5- Profile the mapped BAM files

|Inputs                       |Link                                 |
|-----------------------------|-------------------------------------|
| Mapped BAM files            | [mapped_bam](TOBEIMPLEMENTED)       |
| Bed file                    | [bed_file](TOBEIMPLEMENTED)         |
| Genome lengths parquet file | [genome_lengths](TOBEIMPLEMENTED)   |
| Gene range table TSV file   | [gene_range_table](TOBEIMPLEMENTED) |
| STB file                    | [stb_file](TOBEIMPLEMENTED)         |
| Null model parquet file     | [null_model](TOBEIMPLEMENTED)       |

|Outputs                          |Link                                 |
|---------------------------------|-------------------------------------|
| Profile parquet files           |[profile_link](TOBEIMPLEMENTED)      |
| Genome Statistics parquet files |[genome_stats_link](TOBEIMPLEMENTED) |
| Gene Statistics parquet files   |[gene_stats_link](TOBEIMPLEMENTED)   |

Now we have to make a CSV file containing the sample names and paths to the BAM files:

```csv
sample_name,bamfile
sample1,/path/to/mapping_output/sample1.bam
sample2,/path/to/mapping_output/sample2.bam
```

First, build the null model for sequencing error adjustment (you only need to do this once per reference database, using any one of your BAM files):

```bash
zipstrain utilities build-null-model --bam-file mapping_output/sample1.bam --output-file null_model.parquet
```

You can profile the mapped BAM files using ZipStrain with the following command:

```bash
zipstrain profile --input-table <path/to/bam/csv> --stb-file mgnify_mouse_gut_genomes.stb --null-model null_model.parquet --gene-range-table preprofiles/gene_range_table.tsv --bed-file preprofiles/genomes_bed_file.bed --genome-length-file preprofiles/genome_lengths.parquet --run-dir profiling_output/
```

This will generate profile parquet files, genome statistics parquet files, and gene statistics parquet files for each sample in the `profiling_output/` directory.

As an alternative, you can use the Nextflow pipeline to perform the profiling:

```bash
nextflow run zipstrain.nf --mode "profile" --input_table <path/to/bam/csv>  --gene_file mgnify_mouse_gut_genes.fasta --stb mgnify_mouse_gut_genomes.stb  --output_dir profiling_output/ --reference_genome mgnify_mouse_gut_genomes.fa -c conf.config -profile <your/system/specific/profile> -resume
```

### Step 6- Compare the profiled samples at the genome level

### Step 7- Compare the profiled samples at the gene level
