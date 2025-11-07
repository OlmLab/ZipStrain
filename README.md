# ZipStrain

Official Repository for ZipStrain python package. ZipStrain is a bioinformatics tool designed for rapid profiling of metagenomic samples as well as comparative analysis of strain-level variations within microbial communities. 

![ZipStrain Logo](docs/Zipstrain.svg)

## Quick Start

### Installation

For installation and usage instructions, please refer to the [Documentation](docs/installation.md) or the documentation website [here](https://zipstrain.readthedocs.io/en/latest/).

### Profile multiple bam files

You can profile multiple BAM files using either the ZipStrain command-line interface (CLI) or Nextflow Pipeline. Below are examples of both methods.
For both you need to prepare a csv file containing the paths to the bam files to be profiled and the sample names. Example of such a csv file:

```csv
sample_name,bamfile
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
```
#### ZipStrain CLI

To profile multiple BAM files using the ZipStrain CLI, you can use the following command:

```bash
zipstrain run prepare-profiling  blah
zipstrain run profile blah
```



#### Nextflow Pipeline

To profile multiple BAM files using Nextflow, you can create a Nextflow script as follows:

### Compare multiple profiled samples

To compare multiple profiled samples, you can use either the ZipStrain CLI or the Nextflow Pipeline. Below are examples of both methods.

#### ZipStrain CLI
To compare multiple profiled samples using the ZipStrain CLI, you can use the following command:

```bash
zipstrain run compare blah

```
#### Nextflow Workflow
To compare multiple profiled samples using Nextflow, you can create a Nextflow script as follows
```nextflow
// Nextflow script to compare multiple profiled samples

```

### Downstream Analysis and Visualization

ZipStrain provides a Python API that allows for basic downstream analysis and visualization of the profiling and comparison results. Here's an example of how to use the API for visualization:

```python
import zipstrain as zs  
blah
```

For more detailed examples and documentation, please refer to Tutorials section in the [Documentation](docs/tutorials.md) or visit the documentation website [here](https://zipstrain.readthedocs.io/en/latest/).