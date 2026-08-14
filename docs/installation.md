# Installation

ZipStrain runs on **Linux and macOS** (including Apple Silicon). There are two ways to use it, and they install differently — pick the one you need:

- the **Python CLI** (`zipstrain map | profile | compare`) — install ZipStrain plus a few external tools. See [Python CLI Installation](#python-cli-installation).
- the **Nextflow pipeline** (`nextflow run OlmLab/ZipStrain`) — install Nextflow and a container engine, which then pull ZipStrain and its tools for you. See [Nextflow Installation](#nextflow-installation).

For the difference between the two, see the [User Manual](./usermanual.md#zipstrain-command-line-interface). Whichever you choose, finish by confirming it works.

!!! tip "Recommended quick path (CLI)"
    `conda create -n zipstrain -c conda-forge -c bioconda python=3.12 zipstrain bowtie2 samtools sylph` installs the CLI and every external tool `map`/`profile`/`compare` need in one command. Then `conda activate zipstrain && zipstrain test`.

---

## Python CLI Installation

### Requirements

| Component | Needed for | Notes |
|---|---|---|
| **Python ≥ 3.12** | everything | 3.12 is the confirmed baseline |
| **samtools** | `profile` (and `map`) | not bundled with the pip wheel — install separately |
| **bowtie2** | `map` | read alignment |
| **sylph** | `map` (auto-reference route) | picks reference genomes from reads |
| **prodigal** | `map --predict-genes` | gene prediction; optional |
| **torch**, **h5py** | `compare --method matrix` | the `[matrix]` extra |

Profiling and comparison need only ZipStrain + samtools. The aligners (bowtie2, sylph, prodigal) are used exclusively by `zipstrain map`.

### Install ZipStrain

#### Option A — Conda (recommended)

Conda installs ZipStrain **and** samtools together, and can pull the aligners in the same command:

```bash
# ZipStrain + samtools + the tools `zipstrain map` needs
conda create -n zipstrain -c conda-forge -c bioconda \
  python=3.12 zipstrain bowtie2 samtools sylph
conda activate zipstrain
zipstrain test
```

Add `prodigal` to that list if you plan to use `zipstrain map --predict-genes`.

Always list the channels in this order — `conda-forge` before `bioconda` — as [Bioconda requires](https://bioconda.github.io/). Do not add the `defaults` channel. If solving is slow, [Miniforge](https://github.com/conda-forge/miniforge) (which ships the fast `mamba` solver and defaults to conda-forge) is a good base; swap `conda` for `mamba` in the commands above.

!!! note "Apple Silicon (M1/M2/M3…)"
    ZipStrain is a `noarch` (pure-Python) package, and its dependencies — `samtools`, `bowtie2`, `sylph`, `prodigal` — all have native `osx-arm64` builds on Bioconda, so a **native Apple Silicon environment is fully supported and preferred** (faster, no emulation).

    The catch is that a Conda installed as **Intel (osx-64)** builds every environment as Intel-under-Rosetta, even on an Apple-Silicon Mac. Check yours:

    ```bash
    conda info | grep platform
    ```

    - `platform : osx-arm64` → you are already native; the command above installs the Apple Silicon builds. ✅
    - `platform : osx-64` on an Apple-Silicon Mac → either install a **native arm64 Conda** ([Miniforge](https://github.com/conda-forge/miniforge) arm64 is the easiest), or force this one environment to arm64:

      ```bash
      CONDA_SUBDIR=osx-arm64 conda create -n zipstrain -c conda-forge -c bioconda \
        python=3.12 zipstrain bowtie2 samtools sylph
      conda activate zipstrain
      conda config --env --set subdir osx-arm64   # keep future installs in this env native
      zipstrain test
      ```

#### Option B — pip

The pip wheel installs ZipStrain and all of its **Python** dependencies, but **not samtools** (or the other external tools). Install into a fresh virtual environment so ZipStrain's pinned dependencies don't clash with other projects:

```bash
python3.12 -m venv zipstrain-env
source zipstrain-env/bin/activate
pip install zipstrain
```

Then install samtools separately — from your package manager (`brew install samtools`, `apt install samtools`), from [htslib.org](http://www.htslib.org/download/), or from bioconda. If you will run `zipstrain map`, also install `bowtie2` and `sylph` (and `prodigal` for `--predict-genes`); see [External tools for `zipstrain map`](#external-tools-for-zipstrain-map).

Verify:

```bash
zipstrain test
```

#### Matrix workflow extra

`zipstrain compare --method matrix` needs `torch` and `h5py`. They are an optional extra so the base install stays light:

```bash
pip install "zipstrain[matrix]"
```

- **Apple Silicon:** the standard `torch` wheel uses Metal through the MPS backend (`--backend torch-mps`). MPS requires a **native** macOS Python environment — Linux containers cannot expose Apple Metal.
- **Linux + NVIDIA GPU:** install the extra, then replace Torch with the CUDA wheel matching your driver:

  ```bash
  pip install "zipstrain[matrix]"
  pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
  ```

- **CPU-only** is always available with `--backend numpy` and needs no Torch at all.

After installing the extra, `zipstrain test` confirms the matrix dependencies are importable.

### External tools for `zipstrain map`

`zipstrain map` turns reads into BAMs by shelling out to external aligners; it is the only part of the CLI that needs them. Install from bioconda:

```bash
conda install -c conda-forge -c bioconda bowtie2 samtools sylph
# only if you use --predict-genes:
conda install -c conda-forge -c bioconda prodigal
```

When `zipstrain map` runs **without** `--reference-fasta`, it uses Sylph to pick reference genomes automatically:

- the Sylph database is downloaded to the path you pass as `--sylph-db` if it is missing (the default GTDB r220 database is **~14 GB**);
- genomes Sylph selects are downloaded and cached under `--genome-cache-dir` for reuse across runs;
- GTDB taxonomy tables are downloaded once (also cached under `--genome-cache-dir`) so `zipstrain profile` can add a `genome_taxonomy` column.

`zipstrain map` is resumable, so an interrupted download/mapping run picks up where it left off on the next invocation (see the [User Manual](./usermanual.md#map)).

### Containerized CLI (Docker / Apptainer)

You can also run the CLI from the published image, without a local Python install.

**Docker:**

```bash
docker run -it parsaghadermazi/zipstrain:<version> zipstrain test
# mount local data to work with your files
docker run -it -v /path/to/your/data:/data parsaghadermazi/zipstrain:<version> bash
```

Optional GPU images for Linux/NVIDIA hosts (for the matrix workflow):

```bash
docker run -it --gpus all parsaghadermazi/zipstrain:<version>-cuda12.4 zipstrain test
```

Available GPU tags: `:<version>-cuda11.8` (CUDA 11 era) and `:<version>-cuda12.4` (CUDA 12.4). Images are published for `linux/amd64`; the base image includes the CPU implementation of every command.

**Apptainer / Singularity** can pull directly from the Docker image — useful on HPC clusters where Docker is unavailable:

```bash
apptainer run docker://parsaghadermazi/zipstrain:<version> zipstrain test
```

### Confirmed working versions

A known-good combination (the current ZipStrain development environment). Newer point releases within the ranges in `pyproject.toml` are expected to work; if you hit trouble, fall back to these.

| Package / tool | Version | Used by |
|---|---|---|
| zipstrain | 1.0.3 | — |
| Python | 3.12 | everything |
| samtools | 1.23 | `profile`, `map` |
| bowtie2 | 2.5.5 | `map` |
| sylph | 0.9.0 | `map` (auto-reference) |
| prodigal | 2.6.3 | `map --predict-genes` |
| torch | 2.x (2.12 tested) | `[matrix]` extra |
| h5py | 3.16 | `[matrix]` extra |

Key Python dependencies (installed automatically): polars 1.42, numpy 2.5, pyarrow 22.0, duckdb 1.5, pandas 2.3, scipy 1.18, pydantic 2.13, click 8.4, rich 14.3, psutil 7.2.

### Troubleshooting

- **`samtools` not found** — required for profiling. Make sure it is on your `PATH`; check with `samtools --version`. With a pip install, samtools is not bundled — install it separately.
- **`zipstrain map` reports a missing tool** — `bowtie2`, `sylph`, or `prodigal` is not on `PATH`. Install the relevant tool from bioconda (see [External tools](#external-tools-for-zipstrain-map)).
- **Python version errors** — ZipStrain requires Python 3.12+. Check with `python --version`.
- **`zipstrain test` fails** — a dependency is missing. Re-check samtools (pip installs) and, for the matrix workflow, that the `[matrix]` extra is installed.

---

## Nextflow Installation

The [Nextflow pipeline](./usermanual.md#nextflow-implementation) is an alternative to the CLI, best for large cohorts and clusters. It needs three things: Nextflow, a Java runtime, and a container engine — which then pull ZipStrain and its tools for you (no separate Python install required).

### Requirements

| Component | Notes |
|---|---|
| **Java ≥ 17** | required by Nextflow (up to Java 26) |
| **Nextflow** | tested against 24.10.0 |
| **Docker** or **Apptainer/Singularity** | runs each pipeline step in the ZipStrain container |

### 1. Java

Nextflow requires **Java 17 or newer**. Check what you have:

```bash
java -version
```

If it is older than 17, install a modern JDK — for example via conda (`conda install -c conda-forge openjdk=17`), [SDKMAN!](https://sdkman.io/), or your OS package manager.

### 2. Nextflow

Install with the official installer or via conda:

```bash
# Official installer (honors NXF_VER to pin a version)
curl -s https://get.nextflow.io | bash
# or
conda install -c bioconda nextflow
```

Verify:

```bash
nextflow -version
```

### 3. Container engine

The pipeline runs each step inside a ZipStrain container, so you need one of:

- **Docker** — on a laptop/workstation;
- **Singularity / Apptainer** — on most HPC systems (no root needed).

The repo ships a `nextflow.config` that enables **Docker by default**, so a laptop run needs no profile — you can run the pipeline straight from GitHub:

```bash
nextflow run OlmLab/ZipStrain --mode ... --input_table ... --output_dir ... -resume
```

On a cluster, keep site-specific SLURM profiles in a local ignored config such as `conf.local.config`, then run with `-c conf.local.config -profile <name>`. The public config only ships generic Docker/Apptainer profiles and uses the published image tag from `nextflow.config`.

### Confirmed working versions

| Component | Version |
|---|---|
| Nextflow | 24.10.0 |
| Java (JDK) | 17+ |
| ZipStrain container | `parsaghadermazi/zipstrain:1.0.3` |

ZipStrain's continuous integration runs the pipeline against **Nextflow 24.10.0**, so that is the confirmed-working version. Pin it with `export NXF_VER=24.10.0` before running to match CI exactly. Newer Nextflow releases (26.x) tightened DSL syntax — the bundled `zipstrain.nf` is kept compatible, but pin to a tested version if you hit parse errors.

### Troubleshooting

- **"Cannot find Java or it's the wrong version"** — install Java 17+ and make sure `java -version` reports it (Nextflow reads `JAVA_HOME`/`PATH`).
- **Nextflow parse/DSL errors** — pin to the tested release with `export NXF_VER=24.10.0` before running.
- **Container not found / pull failures** — confirm your engine is installed and running, and that the image tag in `nextflow.config` or your local site config is reachable from your host or cluster.
