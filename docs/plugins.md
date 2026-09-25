# Plugins

ZipStrain works without plugins. Profiling uses its built-in Python backend unless you select another one. Plugins are installed separately from ZipStrain and selected by name at run time.

## Available plugins

| Plugin | Select with | What it changes | Availability |
|---|---|---|---|
| `zipstrain_rust_profiler` | `--backend rust_profiler` | Profiles BAMs with a native, multithreaded Rust implementation | Separate PyPI package |

Plugin distributions follow the `zipstrain_<capability>` naming pattern. The
distribution and Python import use the same underscore-separated name; the
backend selector is the shorter capability name. Python package indexes treat
underscores and hyphens as equivalent in distribution names.

`python` is the built-in profiling backend, not a plugin. The `zipstrain[matrix]` extra is also not a plugin: it installs dependencies for matrix comparison. Likewise, `--backend` on `zipstrain compare` selects a **comparison** compute backend, not a profiling plugin.

## Rust profiler

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
