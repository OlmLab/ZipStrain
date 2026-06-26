# Tutorial

## Tutorial #1) Running ZipStrain on provided test data using the CLI

The following tutorial goes through an example run of ZipStrain using the python CLI. You can follow along with your own data, or use a small set of reads that are included in the inStrain source code for testing. To download them directly, run:

```bash
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1.R1.fastq.gz
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1.R2.fastq.gz
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1_scaffold_min1000.fa
```

Or browse and download them manually from [GitHub](https://github.com/MrOlm/inStrain/tree/master/test/test_data). The files we'll use are the forward and reverse metagenomic reads (`N5_271_010G1.R1.fastq.gz` and 1N5_271_010G1.R2.fastq.gz1) and a .fasta file to map to (`N5_271_010G1_scaffold_min1000.fa`). In case you're curious, these reads come from a premature infant fecal sample.

### Preparing input files

After downloading the genome file that you would like to profile (the fasta file) and at least one set of paired reads, the first thing to do is to map the reads to the .fasta file in order to generate a bam file.

In this tutorial we will do that using the mapping program Bowtie 2 as follows:

```bash
mkdir bt2

bowtie2-build N5_271_010G1_scaffold_min1000.fa bt2/N5_271_010G1_scaffold_min1000.fa

bowtie2 -p 6 -x bt2/N5_271_010G1_scaffold_min1000.fa -1 N5_271_010G1.R1.fastq.gz -2 N5_271_010G1.R2.fastq.gz > N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam
```