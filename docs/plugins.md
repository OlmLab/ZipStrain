# Plugins

ZipStrain offers a plugin system for extending its core workflows with optional, separately installed components. Plugins can replace or add functionality in a specific part of the pipeline (for example, a faster profiling engine) without changing how the rest of ZipStrain works. The plugin catalog will grow over time.

## Available plugins

| Plugin | Extends | Select with | What it does | Install |
|---|---|---|---|---|
| [`zipstrain_rust_profiler`](#rust-profiler) | Profiling | `--backend rust_profiler` | Profiles BAMs with a native, multithreaded Rust implementation; up to ~10x faster with lower memory use on large BAMs | `pip install zipstrain_rust_profiler` |

ZipStrain works fully without any plugins. Each plugin is installed separately and only takes effect when you select it by name at run time; otherwise ZipStrain uses its built-in behavior.

### Naming conventions

Plugin distributions follow the `zipstrain_<capability>` naming pattern. The
distribution and Python import use the same underscore-separated name; the
run-time selector is the shorter capability name. Python package indexes treat
underscores and hyphens as equivalent in distribution names.

### What is not a plugin

- `python` is the built-in profiling backend, not a plugin.
- The `zipstrain[matrix]` extra installs dependencies for matrix comparison; it is not a plugin.
- `--backend` on `zipstrain compare` selects a **comparison** compute backend, not a profiling plugin.

## Rust profiler

**Extends:** profiling · **Selector:** `--backend rust_profiler`

### Why use it

Profiling time and memory grow with BAM size. The Rust profiler is a drop-in replacement for the built-in Python profiler that:

- **Runs up to ~10x faster**, using a native, multithreaded implementation.
- **Manages memory better on large BAMs**, keeping peak memory lower and more predictable so deeply sequenced samples are less likely to exhaust available RAM on a workstation or HPC node.
- **Produces the same outputs**, so downstream comparison and analysis steps work unchanged.

It is most useful for large or deeply sequenced BAMs, large sample sets, and memory-constrained cluster jobs.

### How it fits into ZipStrain

The Rust profiler reads a coordinate-sorted, indexed BAM and writes the same three core Parquet outputs as Python profiling: `<sample>_profile.parquet`, `<sample>_gene_stats.parquet`, and `<sample>_genome_stats.parquet`. When you use the high-level `zipstrain profile` command, ZipStrain still prepares profiling assets and performs its usual output finalization after the plugin finishes. The comparison commands do not change.

The plugin is optional: installing ZipStrain alone does not install it, and leaving out `--backend` keeps the Python behavior.

### Install from PyPI

In a Python 3.12+ environment, install the plugin with pip:

```bash
python -m pip install zipstrain_rust_profiler
```

This also installs ZipStrain 1.2.0 or newer if needed. On a supported platform, pip downloads a prebuilt wheel, so Rust is not needed at installation time. Building from source requires a Rust toolchain and native build prerequisites, including libclang on Linux.

If using the high-level `profile` workflow, keep `samtools` available as described in [Installation](installation.md). The standalone Rust profiler reads the BAM directly, but the ZipStrain workflow still invokes `samtools index`.

### Select it in the CLI

Use the normal profiling command and add one flag:

```bash
zipstrain profile \
  --input-table mapped/samples.txt \
  --reference-fasta ref.fna \
  --stb-file ref.stb \
  --run-dir profiled \
  --backend rust_profiler
```

`mapped/samples.txt` is the same CSV used by Python profiling, with `sample_name,bamfile` columns. Keep the BAMs coordinate-sorted. `--num-procs` controls the maximum Rust workers per sample; the plugin also respects ZipStrain's detected CPU allocation. Existing profiling filters such as `--min-baseq`, `--min-read-ani`, and `--read-inclusion` are passed through.

For one BAM with already prepared assets, use `zipstrain utilities profile-single`:

```bash
zipstrain utilities profile-single \
  --bam-file sample.bam \
  --bed-file genomes.bed \
  --stb-file reference.stb \
  --null-model null_model.parquet \
  --output-dir sample_profile \
  --backend rust_profiler
```

This command requires an indexed BAM, BED file, STB file, and null-model Parquet. `--reference-fasta`, `--gene-range-table`, and `--profiling-contract` remain optional; provide the reference FASTA if you want `ref_base_bitmask` in the profile. Its worker limit is `--max-concurrency`.

### Select it in Nextflow

For Nextflow modes that run profiling, set:

```bash
--profile_backend rust_profiler
```

For example, add that parameter to the [`--mode profile` command](usermanual.md#command-2-profile-bams-mode-profile). **Install the plugin inside the environment that runs each profiling task.** The standard ZipStrain Docker/Apptainer images do not include this optional package, so a containerized run needs a custom image containing both ZipStrain and the Rust plugin. Installing the plugin only on the machine that launches Nextflow is not sufficient. With the default `python` backend, Nextflow does not pass a plugin flag to the standard images.

If ZipStrain says `Profiling backend 'rust_profiler' is not installed`, check the Python environment or task container where `profile-single` actually runs. ZipStrain does not silently fall back to Python when you explicitly request a plugin.
