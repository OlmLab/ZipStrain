# Tutorial
This tutorial will guide you through the basics of using ZipStrain.

## Introduction to ZipStrain
ZipStrain is a tool for profiling a metagenomics sample against a database of reference genomes and performing comparisons between the profiles as well as downstream analyses. The typical workflow consists of the following steps:

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

|chrom|pos|gene|A|C|G|T|

```

Where `chrom` is the scaffold, `pos` is the position in the reference genome, `gene` is the gene name (if a gene file is provided), and `A`, `C`, `G`, and `T` are the counts of each nucleotide at that position based on the mapped reads.

The information will be generated for every position in the reference genomes that has at least one read mapped to it. There are two ways to generate profiles using ZipStrain: 

- ZipStrain as a command-line tool
- ZipStrain Nextflow pipeline.

#### Using ZipStrain as a Command-Line Tool

UNDER CONSTRUCTION

#### Using the Nextflow Pipeline

You can also use the nextflow pipeline to generate profiles. Here is an example command:

```
nextflow run zipstrain.nf --mode "fast_profile" --input_table all_bams.csv  --gene_file <path/to/genefile> --stb <path/to/stbfile>   --output_dir <path/to/output/dir>  --reference_genome <path/to/your/reference_db.fasta> -c conf.config -profile <your_profile> -resume
```
This step will generate parquet files for each of your samples in the specified output directory.

NOTE: fast_profile mode requires a gene file and a STB file. The gene file MUST BE following Prodigal NUCLEOTIDE format. For generating STB file, use the following command:

```zipstrain utilities generate_stb --genomes-dir-file <path/to/genomes_dir_file>```

If you have trouble generating the STB file, you can simply make one using a custom script. The STB file is a tab-separated file with two columns: scaffold name and genome file name. 

