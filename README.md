# ZipStrain

**Strain-resolution metagenomics at scale.**

ZipStrain profiles mapped reads into per-position nucleotide counts and compares profiles at both genome and gene levels — across hundreds of samples, on a laptop or an HPC cluster.

![ZipStrain Logo](docs/Zipstrain.svg)

## Why ZipStrain?

- **Nucleotide-resolution profiling** — per-position A/C/G/T counts for every genome and gene in your reference database
- **Flexible comparisons** — popANI, conANI, cosANI, identity-by-state (IBS), and identical-gene metrics, all in one tool
- **Built for scale** — Polars and DuckDB engines, automatic batching, Slurm support, and resumable runs
- **Multiple execution modes** — local multi-process, Slurm batch, Docker, or Apptainer
- **Nextflow integration** — end-to-end pipelines from raw reads (or SRA accessions) to comparisons
- **Automatic reference building** — build reference databases directly from Sylph abundance tables
- **Python API** — programmatic access for custom pipelines, downstream analyses, and visualizations

Developed by Parsa Ghadermazi and team at the [Olm Lab](https://www.olmlab.org/), University of Colorado Boulder.

## Quick install

```bash
pip install zipstrain
zipstrain test
```

> Requires Python 3.12+. Samtools must be installed separately ([htslib.org](http://www.htslib.org/download/)).
> Also available via Conda (`conda install -c conda-forge -c bioconda zipstrain`) and Docker (`parsaghadermazi/zipstrain`).

## Documentation

Full documentation, tutorials, and API reference:

**[https://OlmLab.github.io/ZipStrain/](https://OlmLab.github.io/ZipStrain/)**

| Page | What you'll find |
|------|-----------------|
| [Installation](https://OlmLab.github.io/ZipStrain/installation/) | pip, Conda, Docker, and Apptainer setup |
| [Tutorial](https://OlmLab.github.io/ZipStrain/Tutorial/) | End-to-end walkthrough from reads to comparisons |
| [CLI Reference](https://OlmLab.github.io/ZipStrain/cli/) | Every command and option |
| [Nextflow Pipelines](https://OlmLab.github.io/ZipStrain/NextflowPipeline/) | HPC-ready workflows |
| [Build Genome DB](https://OlmLab.github.io/ZipStrain/GenomeDBFromSylph/) | Automatic reference building from Sylph |
| [Python API](https://OlmLab.github.io/ZipStrain/api/) | Programmatic access and downstream analysis |

## Minimal example

```bash
# 1. Prepare profiling assets
zipstrain utilities prepare_profiling \
  --reference-fasta reference.fasta \
  --gene-fasta genes.fna \
  --stb-file mapping.stb \
  --output-dir assets/

# 2. Profile BAM files
zipstrain profile \
  --input-table samples.csv \
  --stb-file mapping.stb \
  --null-model null_model.parquet \
  --gene-range-table assets/gene_range_table.tsv \
  --bed-file assets/genomes_bed_file.bed \
  --genome-length-file assets/genome_lengths.parquet \
  --run-dir profiles/

# 3. Compare genomes
zipstrain compare genomes \
  --genome-comparison-object genome_compare.json \
  --run-dir comparisons/ \
  --calculate all \
  --engine duckdb
```

## Citation

If you use ZipStrain in your research, please cite:

> *Citation coming soon.*

## License

See [LICENSE](LICENSE) for details.
