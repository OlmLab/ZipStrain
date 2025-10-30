# Tutorial
This tutorial will guide you through the basics of using ZipStrain.

## Introduction to ZipStrain
ZipStrain is a tool for profiling a metagenomics sample against a database of reference genomes and performing comparisons between the profiles as well as downstream analyses. The typical workflow consists of the following steps:

1- Mapping reads to a reference database using any read mapper of your choice (e.g., BWA, Bowtie2, Minimap2).

2- Generating a profile from the mapped reads using ZipStrain. This step will generate a parquet file containing nucleotide frequencies for all positions in the reference genomes.

3- Comparing the generated profiles against each other and/or against a reference profile.

4- Performing downstream analyses, such as Identity-By-State (IBS), StrainSharing, and clustering.

### Mapping Reads

