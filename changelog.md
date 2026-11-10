# Changelog

Entries are brief by design and describe changes relative to the previous released version.

## Unreleased

V2 development (`2.0.0.dev0`):

- Preprofile writes an authenticated `gene_info.parquet` catalog and reference-codon definitions. Content-derived gene IDs and reference/annotation fingerprints replace unscoped row numbering.
- New position profiles omit gene labels; gene summaries use overlap-safe interval aggregation and the agreed 15-column schema.
- The V2 Rust profiler writes gene-free worker chunks and final position tables directly, removing native gene-range loading and redundant gene summaries. Shared V2 code produces gene statistics and optional selection metrics.
- Optional observed codon counts preserve read phasing, deduplicate overlapping mates, and feed SNS dN/dS and SNV pN/pS using the `gene-overlap` consensus/weighted integration definitions. Full-codon context fixes the shared-nonreference-background classification bug; formulas and exclusions are documented in the V2 guide.
- Explicit-input Polars pair and collection comparisons use shared preprofile assets, validate source evidence, and retain content-addressed completed-pair results for resuming/appending samples. CLI, Nextflow, matrix readers, and the native-backend adapter understand the new files.
- V2 compares into one Parquet. Genome-only requests use genome rows. Requesting gene metrics uses gene rows and includes any requested genome metrics on those rows.
- Breaking: `--calculate` is the one way to choose metrics, on every entry point. It accepts `genome_ani`, `max_identical_stretch`, `gene_ani`, `gene_selection`, and `all`. The `--dnds` and `--compare-genes` flags and the old metric names (`ani`, `popani`, `ibs`, `max_block`, `identical_genes`, and others) are removed, and an old name points to its replacement.
- Breaking: `all` means every metric the inputs support and warns about what it skips, including `gene_selection` unless every sample has a codon table. Codon tables for only some samples are an error.
- Breaking: comparison tables no longer contain `shared_genes_count`, `identical_gene_count`, `perc_id_genes`, `dN_dS`, or `pN_pS`. Use `gene_info.with_selection_ratios` for the ratios. `calculate_ibs` and `calculate_identical_frac_vs_popani` now read gene rows, take `min_gene_positions`, and count each sample pair once.
- Breaking: `max_consecutive_length` is now `max_identical_stretch` in every output. Existing matrix compare databases and 1.x comparison tables are upgraded when a run extends them.
- Fixed: pair comparisons of V2 profiles without assets no longer produce zero gene summaries, and mixed 1.x and V2 inputs are rejected on every route.
- Behaviour change: `--min-gene-compare-len` now defaults to 1 everywhere (previously 1, 100, or 200 by entry point), with analyses applying 100. It no longer invalidates cached genome-only pairs. Python API read filters now match the CLI (`min_read_ani=0.95`, `read_inclusion=paired`).
- Matrix comparison supports `gene_selection` with weighted and consensus integration. It processes observed codons in NumPy batches, caches completed pairs, and applies the same gene position counts and comparable-length filter as table comparison.
- Comparison offers `--method table` (the default) and `--method matrix`. The table method uses the profile's Parquet sort marker for a Polars streaming merge join. The former DuckDB table engine and its options have been removed.
- The PyPI release workflow benchmarks whole comparisons and isolated pairs on synthetic single and multiple genome workloads, and updates the [performance charts](docs/benchmarks.md) with measured runtime and peak resident memory.
- Release benchmarks compare matrix and table selection workloads in both weighted and consensus modes, and verify parity before recording the results.
- Faster comparison with the same results: a pair now computes all its metrics in one pass over the shared positions, and one selection routine serves reference, table, and matrix comparison. On a test collection of 15 pairs with `gene_selection`, the run went from 505 s to 24 s. Selection values can differ from earlier development builds only by floating-point rounding.
- `zipstrain compare` runs pairs in worker processes within `--memory-limit-gb`, and validates each sample once. The Python API stays serial unless given `workers`.
- Fixed: Slurm and container runs of `zipstrain compare` now accept V2 profiles, using the same pair command as Nextflow. A run directory holds one comparison: rerunning it with the same inputs and settings reuses the result, and other settings need a new `--run-dir`.
- Nextflow compare jobs now receive codon tables, so `gene_selection` works in single and batched modes with values equal to `zipstrain compare`.
- Breaking: `utilities compare-pair` is the one pair-level comparison command, and `--calculate` decides genome or gene rows. It replaces `single_compare_genome`, `single_compare_gene`, and `chunk-genome-compare` (now `--pair-table`), takes `--output-file` instead of `--output-dir`, and `--scope` (`all`, a genome, or `GENOME:GENE`) instead of `--genome`. With 1.x profiles, `gene_ani` on its own gives the legacy gene table. It finds codon tables beside the profiles and accepts `--codon-counts-1/2`, `--allele-integration`, and `--sample-name-1/2`.
- New Python function `compare.identical_stretches` returns every identical stretch of at least `min_length` positions between two profiles, with its genome, scaffold, and one-based inclusive start and end. Its longest stretch equals `max_identical_stretch`. It has no CLI command.
- `identical_stretches` rows now carry `background_divergence`, `expected_runs` (an E-value), and `p_value`, which test each stretch against popANI mismatches falling at random at the genome's rate.
- New Python function and `PairComparison` method `mosaic_regions` finds blocks of unusually low divergence (such as recent recombination) or high divergence (such as imports from a more distant lineage). A three-state hidden Markov model runs over the runs of matches between mismatches, so it costs one step per mismatch, and refits each genome's background rate without the blocks. It returns a lazy frame that shares the pair's join with the other analyses. It has no CLI command.
- New Python class `compare.PairComparison` validates a pair once and exposes its lazy join (`aligned`, `matches`), `compare`, and `identical_stretches`. `collect` runs several analyses together so Polars joins the profiles once. `compare_profiles`, `identical_stretches`, and `align_profiles` now wrap it, with unchanged results.
- Faster start-up: the CLI loads each ZipStrain module only when a command uses it, and `scipy` only when profiling builds its error model. `zipstrain utilities compare-pair` starts in 0.2 s instead of 1.0 s, so every scheduled comparison job finishes about 0.6 s sooner. Shared option defaults and choices now live in `zipstrain.options`; the old module-level names remain as aliases.
- Fixed: `gene_selection` no longer runs out of memory on large catalogs. Codon evidence is checked in streaming passes, then packed and matched in ranges of whole genes, and a one-byte-per-codon reference lookup replaces the join against `reference_codons.parquet`. A pair of 1 GB profiles (128 M codons each) now finishes in about 50 s, where it was previously killed. Sorted codon rows are grouped without hashing, codons that are one identical codon in both samples skip the change kernel, and table pairs and reference selection no longer stage packed files. Opportunity totals are now summed exactly, so results can differ from before in the last digit (within 1e-12).
- Breaking: Nextflow `--mode compare` replaces `compare_genomes` and `compare_genes`, and `--compare_scope` replaces `--compare_genome_scope` and `--compare_gene_scope`. The old modes stop with a pointer to the new one. Nextflow genome comparisons now honour `--compare_ani_method`, which they used to ignore, and Slurm gene runs honour `--scope`.
- Fixed: merging Nextflow comparison tables keeps the V2 identity metadata. Batch file names no longer collide on systems without `sha1sum`. Table comparisons record `zipstrain_compare_engine`. Gene and genome statistics no longer claim a sort order they do not have.
- Fixed: the Rust profiler disagreed with Python on real data because it used a different BAQ mode and left out the trailing scaffold gap. Both backends now give identical positions, codon tables, and gene statistics, and genome statistics that differ only by rounding. Native-backend genome statistics now carry the same metadata as Python's.
- One `--threads N` option on `profile`, `compare`, `utilities profile-single` and `utilities compare-pair` sets the most CPUs the command uses in total; ZipStrain shares them out between worker processes, native pools and concurrent local batches, and in Slurm mode uses them as each job's `cpus-per-task` when the config sets none. `zipstrain.set_threads(N)` does the same from Python. `--num-procs`, `--max-concurrency` and `--workers` are hidden and capped by it.
  - Changed: without `--threads`, `profile` gives each concurrent local batch an equal share of the machine's CPUs instead of 8 CPUs each, and `profile-single` uses every CPU instead of 4 chunks.
  - Fixed: an explicit pair-worker count could exceed the CPU budget.
- Review fixes for the V2 beta:
  - Fixed: `identical_stretches` p-values were too small for genomes split over many scaffolds or with coverage gaps. `expected_runs` now counts one run per mismatch, scaffold start and coordinate gap (`M + S`) instead of `M + 1`. The stretches are found in one streaming NumPy pass that keeps only run boundaries, so they no longer hash every matching position. Like `compare`, the method now runs its own pass when collected.
  - Fixed: `mosaic_regions` no longer returns NaN scores for a genome where every shared position differs.
  - Profile gene statistics (`gene_info.range_sums`) stream the profile once with prefix sums and binary search instead of sorting it, so profiling a 1 GB sample no longer holds the whole profile in memory. Profiles are written once, with their final metadata, instead of three times. Coordinate-order checks now stream too.
  - `PairComparison.compare`, `pairwise_selection` and `reference_selection` are lazy: no codon is read and no profile is hashed until the frame is collected.
  - A gene catalog is authenticated once per process and file version, profile digests are cached the same way, and the matrix route builds each reference codon lookup once per process instead of once per sample.
  - Codon tables whose row groups are not in gene order are sorted once into a temporary copy, with a warning, instead of being read in full once per gene range.
  - Profiling parses mpileup output about 4× faster: `utils.clean_bases` copies the text between read markers and indels in slices instead of matching every character, and `utils.count_mpileup_bases` counts with `str.count`. Counts are unchanged. Each chunk's parse also runs in a spawned worker process, so `--max-concurrency` chunks parse in parallel instead of taking turns on Python's interpreter lock (about 3.4× more with 6 chunks on an 8-core Mac).
  - 1.x genome comparisons use the same streaming scan as V2, with the 1.x stretch definition (a run opened by a mismatch counts it), so `max_identical_stretch` values are unchanged. Breaking: the unused helpers `compare.add_contiguity_info`, `compare.get_longest_consecutive_blocks` and `compare.calculate_pop_ani` are removed.
- The V2 output formats are frozen for 2.0 and covered by contract tests; see [V2 development](docs/v2.md).
- Requires Polars 1.36.1 or newer; earlier versions cannot run V2 comparison. The lock file moves to Polars 2.1.0.
- Every new function has tests and documentation; see [V2 development](docs/v2.md).

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
