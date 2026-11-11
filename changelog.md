# Changelog

## 2.0.0b1

First public beta of ZipStrain 2.0, compared with 1.3.0. Install it in a Python 3.12 or newer environment with `python -m pip install "zipstrain==2.0.0b1"`. For the container, use `parsaghadermazi/zipstrain:2.0.0b1`. To run the matching Nextflow pipeline, use `nextflow run OlmLab/ZipStrain -r v2.0.0b1` with your usual inputs.

### Profile and gene data

- Preprofile now produces `gene_info.parquet` and `reference_codons.parquet`. The gene catalog records coordinates, strand, length, and a content-derived gene ID. Reference and annotation fingerprints let ZipStrain reject incompatible assets before comparing samples.
- Position profiles omit the gene column. The separate `<sample>_gene_stats.parquet` contains gene and genome identifiers, coordinates, length, coverage, breadth, reference ANI, and the six selection fields `dN`, `dS`, `dN_dS`, `pN`, `pS`, and `pN_pS`. Genes with overlapping intervals are counted independently.
- `--collect-codons` saves `<sample>_codon_counts.parquet` and enables selection metrics against the reference. Codon counts come from complete triplets observed on reads. Overlapping mates count once, and retained codons use alleles present in the position profile. Consensus and weighted integration preserve read phasing for synonymous and nonsynonymous counts.
- The V2 Rust profiler source supports gene-free position output and native codon extraction. It shares the gene statistics finalizer with Python. The backends now use the same BAQ mode and account for trailing scaffold gaps. For the V2 Rust backend, see the source installation instructions in the [V2 guide](docs/v2.md).

### Comparison and migration from 1.x

- `zipstrain compare` offers `--method table` (the default, using Polars) and `--method matrix`. Both support genome ANI, gene ANI, maximum identical stretch, and gene selection. The matrix `anchor` kernel supports gene metrics. The `block` kernel is limited to genome ANI.
- Select metrics with `--calculate genome_ani`, `max_identical_stretch`, `gene_ani`, or `gene_selection`. Combine them with `+`. `all` runs the metrics supported by the inputs and warns about omissions. When selection is requested, codon evidence must cover every sample. Replace `--dnds` with `--calculate gene_selection` and `--compare-genes` with `--calculate gene_ani`. Old metric aliases and the selectable DuckDB table engine have been removed.
- Each run writes one comparison Parquet, with CSV available via `--force-csv`. Genome-only results have one row per sample pair and genome. Gene requests add `gene_id` to the key, and requested genome metrics repeat across the corresponding gene rows. `max_consecutive_length` is now `max_identical_stretch`.
- Comparison tables store `dN`, `dS`, `pN`, and `pS` independently. Derive the ratios with `gene_info.with_selection_ratios`. The downstream IBS and identical-gene analyses now derive gene summaries from gene rows and count each sample pair once. `--min-gene-compare-len` defaults to one shared position, and those analyses apply their own default of 100.
- Use `zipstrain utilities compare-pair` for a single pair or `--pair-table` for a batch of genome-only comparisons. It replaces `single_compare_genome`, `single_compare_gene`, and `chunk-genome-compare`. Use `--output-file` for the result path and `--scope` for a genome or, on a single pair, `GENOME:GENE`.
- A comparison can use compatible 1.x inputs together or compatible V2 inputs together. Mixed versions are rejected. The 1.x comparison path retains its stretch counting convention, and existing matrix checkpoints and legacy comparison tables are upgraded when extending a run. For V2 gene selection, prepare the shared V2 assets and profile each sample with `--collect-codons`.

### Python API and execution

- `PairComparison` validates two profiles once and exposes their aligned positions, matches, and analysis methods. Use `pair.collect(...)` to share the alignment across compatible lazy queries. `pair.compare(...)` and `pair.identical_stretches(...)` run their own streaming passes at collection time.
- `identical_stretches` returns one-based inclusive coordinates and significance estimates against each genome's background divergence, accounting for scaffold boundaries and coverage gaps. `mosaic_regions` returns unusually low- or high-divergence intervals from a three-state model with a refitted background. Both are available from Python.
- The legacy Python helpers `add_contiguity_info`, `get_longest_consecutive_blocks`, and `calculate_pop_ani` were removed. Use `PairComparison` for pair analysis and `identical_stretches` for stretch coordinates.
- `--threads` on profiling and comparison commands sets the CPU budget for local workers and native pools. Slurm jobs use it when the job configuration leaves CPUs per task unset. Python callers can use `zipstrain.set_threads(...)`. The Python profiling defaults now match the CLI at `min_read_ani=0.95` and `read_inclusion="paired"`.
- Nextflow uses `--mode compare` and `--compare_scope`, replacing the separate genome and gene comparison modes and scope options. Single and batched jobs pass codon evidence through to the pair command. Nextflow honors `--compare_ani_method` for genome comparisons and preserves V2 comparison metadata when merging results.

### Performance, reliability, and requirements

- The table method reads the Parquet sort marker and uses a Polars streaming merge join for sorted profiles. Each pair computes position metrics in one pass. Collection comparisons validate samples once and reuse completed pairs when resuming or appending to compatible runs.
- Gene statistics use streaming interval aggregation. Gene selection processes whole-gene ranges with a compact reference-codon lookup shared across comparisons. Profiling writes the position Parquet with its final metadata in one pass. The CLI loads command dependencies on demand, and Python profiling parses mpileup chunks in parallel.
- Atomic output renames retry transient busy-file errors on shared filesystems for up to 30 seconds. If the file remains busy, the error identifies the affected path.
- The package requires Polars `>=1.36.1,<3`. Tests cover the V2 Parquet schemas, codon calculations, table and matrix parity, and the optional Rust backend. Release benchmarks measure elapsed time and peak memory for whole collections and individual pairs, across single and multiple genomes and both allele-integration modes. See the [V2 guide](docs/v2.md) and [benchmark charts](docs/benchmarks.md) for examples and results.

## 1.3.0

Compared with `1.2.3`:

- Comparison metric options now use the descriptive names `max_identical_stretch`, `gene_ani`, and `genome_ani` across the CLI, Python API, task manager, and Nextflow. The previous names remain accepted; existing output columns and resume contracts are preserved.

- The main `zipstrain compare --method matrix` command now honors `--calculate` and accepts `--compare-kernel`, making ANI-only block comparisons available without the lower-level utilities command. Gene comparisons automatically include gene ANI, and resumed runs validate the selected calculations before appending samples.
- Renamed comparison kernels to `anchor` (default, all metrics) and `block` (genome-level popANI/conANI only). The old `pairwise` and `gemm` names remain accepted with deprecation warnings; output values and schemas are unchanged.
- Matrix comparison progress now refreshes `zipstrain_status.json` during the all-vs-all step without adding a log line for every progress event. Eligible ANI-only anchor runs display a tip about the block kernel.

## 1.2.3

Compared with `1.2.2`:

- Matrix comparison: added an opt-in GEMM kernel (`matrix_compare(..., compare_kernel="gemm")`, CLI `--compare-kernel gemm`) for genome-level popANI and conANI. It computes whole blocks of sample pairs with exact `float32` matrix products and writes the same rows and completion keys as the pairwise kernel, so runs stay resumable and interchangeable. It requires `--calculate ani`; IBS, gene ANI and cosANI still use the pairwise kernel.
- GEMM count-store reads use bounded conversion batches when filling the allele-mask cache, avoiding a full raw-count allocation before packing. Dense count rows can be tiled across positions; sparse rows that exceed the conversion workspace fail before allocation.
- Genome statistics: include trailing scaffold gaps in gap-uniformity estimates without inflating mapped-read counts.

## 1.2.2

Compared with `1.2.1`:

- Support Polars 2.0 without breaking Polars 1.x by naming columns explicitly when reading headerless STB, BED, and gene-range tables. Remove an incorrect sortedness assertion in genome annotation, and test both Polars major versions in CI.

## 1.2.1

Compared with `1.2.0`:

- Matrix comparison: preserve sample and genome IDs when resuming against rebuilt, reordered, or expanded HDF5 matrices. Completed comparisons, including zero-overlap pairs and gene results, are reused without rewriting prior rows.

## 1.2.0

Compared with `1.1.1`:

- Profiling: added an optional plugin backend interface and the separately installable `zipstrain_rust_profiler` package; Python remains the default backend.
- Documentation: explained the Poisson null model and its assumptions in the parameter-selection guide; clarified the experiments, fit all six figures to the tutorial width without internal scrolling, reduced chart sizes, added missing figure numbers, and simplified labels and explanations.

## 1.1.1

Compared with `1.1.0`:

- Keep DuckDB unloaded on Polars-only CLI paths and bound both DuckDB's default pool and explicit connections at creation, avoiding node-sized thread pools on shared HPC nodes.

## 1.1.0

Compared with `1.0.4`:

- Resource limits: bound ZipStrain's Polars, DuckDB and other native worker pools to each job's CPU allocation; divide Nextflow batched-comparison CPUs among concurrent pairs (default one), and cap in-process profiling/chunk-compare workers. Nextflow also assigns a per-pair DuckDB memory limit from task memory unless overridden.
- Profiling: tightened the default null-model `p_threshold` to `1e-6` across Python, CLI and Nextflow; retained `error_rate=0.001` and `min_freq=0.01`. Rebuild existing null models to apply the new threshold. Stricter filtering can reduce sensitivity to rare alleles.
- Documentation: added simulated-read accuracy and minor-allele sensitivity results to the parameter-selection guide, with interactive figures, source data and a regeneration script; distinguished depth stability from accuracy and updated the recommendation.

## 1.0.4

Compared with `1.0.3`:

- Matrix comparison: reduced gene-prefix memory, reused target preprocessing, overlapped CPU IBS with Torch computation, and bounded result queues with a persistent batched Arrow/DuckDB writer. Memory planning now accounts for metrics, storage, and queues; fixed noncontiguous anchor loading when resuming. Calculation semantics and output columns are unchanged.

## 1.0.3

Compared with `1.0.2`:

- Standard comparison now breaks `max_consecutive_length` runs at gaps in genomic coordinates, matching matrix comparison semantics.
- Tutorial #1 comparison output was regenerated with current defaults and now documents coverage-intersection and cross-method output semantics accurately.
- Container builds pin each Torch variant and can resolve CUDA dependencies missing from PyTorch's version-specific indexes through PyPI.

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
