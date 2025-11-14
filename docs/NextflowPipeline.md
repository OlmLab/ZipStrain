# Nextflow Pipeline for ZipStrain

The nextflow pipeline for ZipStrain supports:

- Mapping metagenomics reads to a reference fasta file using Bowtie2.

- Generating ZipStrain profiles from the mapped reads (BAM files).

- Comparing ZipStrain profiles with each other and generating a merged comparison file.

In the following sections, we will describe how to use the nextflow pipeline to perform each of these steps.

## Installation of Nextflow

You can install Nextflow by following the instructions on the [Nextflow website](https://www.nextflow.io/docs/latest/getstarted.html#installation).

An easy way to install Nextflow creating a conda environment and installing Nextflow using conda:

```
conda create -n nextflow -c bioconda nextflow
```
Activate the environment using:

```
conda activate nextflow
```

Now create a directory for your Nextflow pipeline and navigate to it:

```
mkdir zipstrain_nextflow
cd zipstrain_nextflow
```
You can then download the `zipstrain.nf` file and the `conf.config` file from the ZipStrain GitHub repository into this directory.

## Deciding on the Execution Environment

ZipStrain is provieded as a containerized application using Docker and can be run using Apptainer as well. For this to work, you need to have either Docker or Apptainer installed on your system. You can run all the nextflow commands by using the config file provided in the ZipStrain GitHub repository (conf.config).

```
Nextflow run zipstrain.nf <command> -c conf.config -profile <docker|apptainer> -resume
```

Nextflow also supports running on various cloud platforms. Please refer to the Nextflow documentation for more details on setting up your execution environment.
Once you build a profile for your execution environment, you can put that in another config file and  use the commands above with adding that config file to ZipStrain's config file.

```
Nextflow run zipstrain.nf <command> -c conf.config,<your_config_file> -profile <your_profile_defined_in_<your_config_file>> -resume
```

## Mapping Reads

ZipStrain Nextflow pipeline supports mapping reads using Bowtie2. Here is an example command to map reads:

```
nextflow run zipstrain.nf --mode "map_reads" --input_table <input_table.csv> --input_type <"local"|"sra"> --reference_genome <reference_genomes.fasta> -c conf.config,<your_config_file> -profile <your_profile_defined_in_<your_config_file>> -resume
```
If you have your reads, in this command, replace `<input_table.csv>` with the path to your input table file. This will be a table in CSV format with three columns: `sample_id`, `reads1`, and `reads2`. The `reads1` and `reads2` columns should contain the paths to the FASTQ files for each sample. In this case `<input_type>` should be set to `local`. Alternatively, if you have SRA accessions, you can provide those in a table that must have a "Run" column containing the SRA accessions (In complience with SRA metadata). In this case, set `<input_type>` to `sra`.

Reference genomes should be the result of concatenating all the reference genome FASTA files into a single FASTA file. 

Running this command will generate BAM files for each of your samples in the specified output directory.

## Profile Generation

Once you have your BAM files, you can generate profiles using the nextflow pipeline. For this, you need to prepare a simple input table in CSV format as follows:

```csv
sample_name,bamfile
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
```

Now you can run the following command to generate profiles:
nextflow run zipstrain.nf --mode "fast_profile" --input_table <path/to/bam/csv>  --gene_file <path/to/reference/fasta/genes> --stb <path/to/stb/file>  --output_dir <path/to/save/generated/files> --reference_genome <path/to/reference/fasta> -c conf.config,<your_config_file> -profile <your_profile_defined_in_<your_config_file>> -resume

```

- `--gene_file`: Path to the gene file in Prodigal NUCLEOTIDE format.
- `--stb`: Path to the STB file. This is a tab-separated file with two columns: scaffold name and genome file name. This table should not have a header line.
```
This command will generate profile, breadth, and scaffols files for each of your samples in the specified output directory.

## From SRA to profiles

To save some time and space, you can directly go from SRA accessions to profiles using the following command:

```

nextflow run instrain2.nf --mode "from_sra_to_profile" --input_table <path/to/sra/csv> --gene_file <path/to/reference/fasta/genes>  --reference_genome <path/to/reference/fasta>   --stb <path/to/stb/file> -c conf.config,<your_config_file>   -profile <your_profile_defined_in_<your_config_file>> --output_dir <path/to/save/generated/files>   -resume

```
This command will download the SRA files, map the reads to the reference genome, and generate profiles in one step and then DELETEs all the intermediate files to save space.

## Comparing Profiles

Once you have generated profiles, you have two options to compare them using the nextflow pipeline.

### Comparing a list of profile pairs

First, you have a CSV file that includes the pairs of profiles to compare. The table must have the following columns:

- sample_name_1

- sample_name_2

- profile_location_1

- scaffold_location_1

- profile_location_2

- scaffold_location_2


```
nextflow run zipstrain.nf --mode fast_compare \
 --input_table <path/to/output/remaining_pairs.csv> \
 --input_type "pair_table" --gene_file <path/to/gene/fasta/file> \
 --reference_genome <path/to/reference/genome.fasta> \
 --stb <path/to/stb/file.stb> -c conf.config,<your_config_file> \
 --output_dir "<path/to/output/directory>" \
 --compare_genome_scope "all" \
 --compare_memory_mode "heavy" \
 --parallel_mode "batched" \
 --batch_size <batch_size> -profile <profile_name> \
 -resume

```

If you want to perform the comparison for a specific genome only, you can use the ` --compare_genome_scope` and give the genome name as follows:

```
--compare_genome_scope "<genome_name>"
```
<genome_name> should be the exact name of the genome as it appears in the STB file.

If your reference database is very big, you might hit memory issues during comparison. In that case, you can use the `--compare_memory_mode` option to reduce memory usage. There are two modes available: "light" and "heavy". The "light" mode uses less memory but is slower. If this still happens you can provide also --compare_chrom_batch_size parameter to limit the number of scaffolds loaded into memory at once.

Finally, you can also control the parallelization of the comparison step using the `--parallel_mode` option. There are two modes available: "single" and "batched". In batched mode, multiple number of pairs will be processed in a single task controlled by `--batch_size` parameter. "single" mode submits one comparison per task. When performin huge number of comparisons, batched mode is recommended to reduce the overhead of task submission as well as number of small output files generated.

### Comparing all profiles against each other

In this case, you provide a CSV file that has all the profiles you want to compare and nextflow will do every possible non-redundant pairwise comparison.

```
nextflow run zipstrain.nf --mode fast_compare \
 --input_table <path/to/profiles/csv> \
 --input_type "profile_table" --gene_file <path/to/gene/fasta/file> \
 --reference_genome <path/to/reference/genome.fasta> \
 --stb <path/to/stb/file.stb> -c conf.config \
 --output_dir "<path/to/output/directory>" \
 --compare_genome_scope "all" \
 --compare_memory_mode "heavy" \
 --parallel_mode "batched" \
 --batch_size <batch_size> -profile <profile_name> \
 -resume

```