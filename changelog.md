# Changelog

Entries are brief by design and describe changes relative to the previous released version.

## 1.0.2

Compared with `1.0.1`:

- Release metadata: bumped the Python package, citation metadata, Nextflow container defaults, and installation documentation to `1.0.2`.

## 1.0.1

Compared with `1.0.0`:

- Readability: clarified several dense methods without changing behavior — the ANI expressions in `PolarsANIExpressions` (`popani`, `conani`, `generalized_cos_ani`) now cast base counts to `Int64` once up front instead of repeating the cast inline on every column, and the contiguous-block grouping in `add_contiguity_info` spells out its break conditions. Outputs are unchanged.
- Profiling: fixed the null-model boundary so counts equal to the maximum plausible error count are removed; set the default assumed error rate to 0.1% and raised the default coverage ceiling to 50,000; coverage above that ceiling now fails explicitly instead of being silently zeroed.
- Profiling: added `--min-freq` (default `0.01`) to the CLI, task-manager, and Nextflow profiling paths, requiring a 1% within-position allele frequency by default.
- Visualization: fixed the scikit-learn silhouette path for thresholds that produce one singleton cluster per sample; these undefined scores are now recorded as `NaN` instead of raising an error.
- Matrix comparison: added compact packed-allele bitmask storage for popANI and exact `A,T,C,G` count storage for popANI, conANI, and thresholded cosANI on NumPy or Torch backends. Matrix coverage filtering is now fixed at build time, count dtypes are overflow-checked (`auto|uint16|uint32`), and sparse stores preserve nonzero values.

## 1.0.0

First stable release. This version consolidates ZipStrain around three top-level commands — `map`, `profile`, and `compare` — plus a `test` environment check and a `utilities` group of lower-level helpers. Compared with `0.11.3`:

### Command-line interface

- New streamlined top level: `zipstrain map | profile | compare | test`, with the previous lower-level helpers gathered under `zipstrain utilities`.
- `compare` is now a single command: it defaults to genome-level comparison, takes `--compare-genes` for gene-level comparison, and selects the engine with `--method standard|matrix`. The old `zipstrain compare genomes` / `zipstrain compare genes` subcommands are gone.
- `compare --profile-db` accepts a CSV of `profile_name,profile_location` directly, so a separate `build-profile-db` step is no longer required.
- `compare` runs are resumable and extendable: re-running with the same `--run-dir` and a profiles table that adds samples computes only the new pairs. Both methods write a single `all_comparisons.parquet` (or `all_gene_comparisons.parquet`) at the top of the run directory.
- `compare --calculate` now accepts fine-grained metric combos (for example `ani+ibs+identical_genes`) through the batched standard path, not just `all`/`ani`.
- Utility cleanup: `generate-genome-pairs` is now `generate-sample-pair`; redundant coverage stats, BED, genome-length, and strain-heterogeneity helper commands were removed.
- Grouped, sectioned `--help` output for the main commands.

### Mapping (`zipstrain map`)

- New command that turns a reads table into sorted, indexed BAMs ready for `profile`, and emits the reference FASTA + STB and a `samples.txt`.
- Omitting `--reference-fasta` runs Sylph to auto-pick, download, and cache reference genomes from the reads, then maps against the built reference; a `reference_genomes_taxonomy.tsv` is written that `profile` auto-discovers.
- Resumable: completed Sylph tables, the built reference, the Bowtie2 index, and each sorted+indexed BAM are reused on re-run (`--force` to rebuild). Optional `--predict-genes` (Prodigal) and `--non-competitive` mapping.

### Profiling (`zipstrain profile`)

- Profiling assets (null model, BED, genome lengths, gene range table, profiling contract) are auto-generated into `<run-dir>/profiling_assets` and cached, so a minimal run needs only `--input-table`, `--reference-fasta`, and `--stb-file`.
- Reference-aware profiling: `--reference-fasta` adds `ref_base_bitmask` to profiles and `ref_ani` to gene/genome stats.
- Per-sample SNV/SNP calling (`<sample>_SNVs.parquet`), with inStrain-style stats in genome/gene stats: `SNS_count`, `SNV_count`, `conANI_reference`, `ber`, `fug`, coverage summaries, and an automated `present`/`absent` `presence` call (BER/FUG thresholds).
- Read/base filters aligned with inStrain and samtools: `--min-mapq`, `--min-baseq`, `--min-read-ani` (NM-tag based), and `--read-inclusion proper-pairs|paired|all-mapped` (default `paired`).
- Optional `--genome-taxonomy` (auto-discovered from the Sylph route) adds a `genome_taxonomy` column; profile outputs are parquet-only.

### Comparison outputs and metrics

- Genome comparison tables keep the `0.11.2` metric schema, but the ANI column is now named `genome_ani`; the parquet header stores `zipstrain_compare_ani_method` so the metric can be interpreted as `popani`, `conani`, or `cosani_<threshold>`.
- Matrix comparison is integrated as `compare --method matrix` (HDF5 matrix store, numpy/torch CPU/CUDA/MPS backends, optional sparse storage), while the low-level `utilities build-matrix-db` / `append-matrix-db` / `matrix-compare` / `matrix-compare-export` / `matrix-db-to-hdf5` remain available.
- `utilities get-snp-reference` can emit either profile-like parquet or a site-only VCF.
- `utilities parquet-to-csv` converts any parquet table to CSV explicitly.

### Reliability and logging

- Universal run logging: `map`, `profile`, and `compare` write a `zipstrain_run.log` and `zipstrain_status.json` so a run's state (running / finished / crashed) is always inspectable.
- Profiling subprocess hardening: each task runs in its own process group with a timeout and group-kill on failure; fixed a profiling concurrency deadlock and now require same-scaffold read pairs by default.

### Nextflow

- Added a root `nextflow.config` so `nextflow run OlmLab/ZipStrain` works out of the box with Docker on a laptop; HPC Singularity/SLURM profiles live in `conf.config`.
- Nextflow modes (`map_reads`, `from_sra_to_profile`, `profile`, `compare_genomes`, `compare_genes`) call the updated CLI and inherit run logging.

### Documentation and packaging

- Reorganized docs into a Quick Start, a comprehensive User Manual (CLI + Nextflow reference), an Expected Output reference, a Downstream Analyses guide, and updated Installation/Tutorial pages; added quickstart GIFs.
- Fleshed-out `pyproject.toml` metadata with `matrix` (torch, h5py) and `test` extras.

## 0.11.3

Compared with `0.11.2`:

- Profiling: each profiling run now uses a private run-local reference FASTA path and waits for `samtools faidx` before launching parallel mpileup chunks, avoiding shared `.fai` races.
- Docker: stabilized the Conda install step by using explicit strict `conda-forge`/`bioconda` channel order in a single solve.

## 0.11.2

Compared with `0.11.1`:

- Profiling: replaced asyncio subprocess orchestration with thread-pooled synchronous subprocess pipelines for raw mpileup/read-location chunk generation before CPU-heavy postprocessing begins.
- Different linkage method support for vizualization module

## 0.11.1

Compared with `0.11.0`:

- Visualization: split silhouette analysis into a compute step and a plotting step.
- Visualization: added peak summarization for silhouette curves, including candidate peaks and best-peak selection.
- Visualization: kept `get_silhouette_plot` as a convenience wrapper around the new compute/plot pair.
- Matrix comparison: reverted the IBS torch path to the 0.11.0 behavior while keeping the visualization split/peak work.
- Reference SNP export: `zipstrain utilities get-snp-reference` can now emit either profile-like parquet or site-only VCF output.
- Visualization: clustering helpers now accept a configurable hierarchical linkage method while keeping `average` as the default.

## 0.11.0

Compared with `0.10.2`:

- Visualization: replaced the slow Polars pivot path in similarity-matrix preparation with direct NumPy matrix construction.
- Visualization: silhouette analysis now uses `scikit-learn` when available and warns when falling back to the manual implementation.
- Profiling: added explicit read filters for MAPQ, base quality, read ANI, and read inclusion mode.
- Profiling/reference: profiles can include `ref_base_bitmask`, and gene/genome stats can include `ref_ani` when a reference FASTA is provided.
- Matrix workflow: matrix builds now require an explicit BED+STB contract; sparse HDF5 storage is supported.
- Documentation: refreshed CLI, tutorial, installation, and Nextflow docs to match current behavior.
