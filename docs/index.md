# Welcome to ZipStrain

ZipStrain is a toolkit for strain-level metagenomic profiling and comparison.
It supports both the standard profile-table comparison workflow and the matrix-store workflow for repeated large all-vs-all comparisons.

![ZipStrain Logo](Zipstrain.svg)

## What ZipStrain Does

ZipStrain helps you:

- profile mapped metagenomic samples at nucleotide resolution
- compare samples at the genome level with ANI and IBS
- compare samples at the gene level
- scale to larger cohorts with resumable matrix-store comparisons
- run either from the CLI, the Python API, or the bundled Nextflow pipeline

## Comparison Routes

ZipStrain now supports two comparison styles:

1. **Standard profile compare**
   - compares profile parquets directly through table operations
   - best when one run needs to evaluate many genomes at the same time

2. **Matrix compare**
   - builds a reusable matrix store from profile parquets
   - best for repeated many-sample comparisons focused on one genome or a small set of genomes
   - supports resumable compare databases and optional gene ANI

If you are new to the project, start with the [Tutorial](./Tutorial.md).
It now includes worked examples for:

- the standard Python/CLI workflow
- the standard Nextflow workflow
- the matrix-store workflow

## Documentation Map

- [Tutorial](./Tutorial.md): end-to-end workflow plus worked examples for standard CLI, standard Nextflow, and matrix compare
- [CLI Reference](./cli.md): command-by-command usage
- [Installation](./installation.md): environment setup and optional matrix dependencies
- [Nextflow Pipeline](./NextflowPipeline.md): pipeline execution details
- [API Reference](./api.md): Python modules and programmatic entry points

## Highlights

- **Strain-level resolution**: profile and compare samples at nucleotide resolution
- **Multiple execution routes**: standard profile compare, matrix compare, and Nextflow workflows
- **Resumability**: matrix compare writes resumable DuckDB result databases
- **Scalability**: supports CPU, CUDA, and Apple Silicon torch backends for matrix compare
- **Container support**: Docker and Apptainer/Singularity friendly
