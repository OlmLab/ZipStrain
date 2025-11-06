# Tutorial
This tutorial will guide you through the basics of using ZipStrain.

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

### Comparing Profiles 
In this step, you compare each profile at nucleotide level against another profile and aggregate the comparison results at genome level. There are two ways to perform comparisons using ZipStrain:

- ZipStrain as a command-line tool

- ZipStrain Nextflow pipeline.

#### Using ZipStrain as a Command-Line Tool

UNDER CONSTRUCTION

#### Using the Nextflow Pipeline

You can also use the nextflow pipeline to perform comparisons. Here is an example command:

```
nextflow run zipstrain.nf --mode "fast_compare" --input_table all_profiles.csv  --null_model <path/to/null_model.parquet> --stb <path/to/stbfile>   --output_dir <path/to/output/dir>  --reference_genome <path/to/your/reference_db.fasta> -c conf.config -profile <your_profile> -resume
``` 
This step will generate comparison parquet files for a every batch of profile pairs as well as a merged comparison parquet file in the specified output directory. For more information about the Nextflow pipeline, please refer to the [Nextflow Pipeline Documentation](./NextflowPipeline.md).

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

- scaffold_location: The location of the scaffold file in TSV format that is generated during the profiling step. This file basically contains any scaffold names present in the reference database which have at least one read mapped to them in the sample.

- reference_db_id: The ID of the reference database. This could be the name or any other identifier for the database that the reads are mapped to. Could be used for filtering profiles based on reference database later on.
    
- gene_db_id: The ID of the gene database in fasta format. This could be the name or any other identifier of your choice for the database that the reads are mapped to. Could be used for filtering profiles based on gene database later on.