# Installation

Follow the steps below to set up ZipStrain on your system.

## Requirements

- **Python 3.12 or higher**
- **samtools** (not bundled with the pip install — see below)

## Install using pip

1. Make sure you have Python 3.12+ installed and your environment is activated.

2. Install ZipStrain:

    ```bash
    pip install zipstrain
    ```

    This installs all Python dependencies **except samtools**. Install samtools separately from [htslib.org](http://www.htslib.org/download/) or via your package manager (e.g., `brew install samtools`, `apt install samtools`).

3. Verify the installation:

    ```bash
    zipstrain test
    ```

## Install using Conda

Conda installs ZipStrain **and** samtools together:

```bash
conda install -c conda-forge -c bioconda -c defaults zipstrain
```

Verify with:

```bash
zipstrain test
```

## Docker

1. Make sure Docker is installed on your system.

2. Run ZipStrain:

    ```bash
    docker run -it parsaghadermazi/zipstrain:latest zipstrain test
    ```

To use ZipStrain interactively or with local data, mount your data directory:

```bash
docker run -it -v /path/to/your/data:/data parsaghadermazi/zipstrain:latest bash
```

## Apptainer (Singularity)

Apptainer can pull directly from the Docker image:

```bash
apptainer run docker://parsaghadermazi/zipstrain:latest zipstrain test
```

This is especially useful on HPC clusters where Docker is not available.

## Troubleshooting

- **`samtools` not found**: ZipStrain requires samtools for profiling. Make sure it is on your `PATH`. Run `samtools --version` to check.
- **Python version errors**: ZipStrain requires Python 3.12+. Check with `python --version`.
- **`zipstrain test` fails**: Ensure all dependencies are installed. If using pip, double-check that samtools is available.
