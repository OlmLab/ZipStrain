# Installation

ZipStrain is a python program that works on mac and linux environments. Follow the steps below to set up ZipStrain on your system.

## Requirements

- **Python 3.12 or higher**
- **samtools** (not bundled with the pip install — see below)

TODO: Add Nextflow versioning

## Install using pip

1. Make sure you have Python 3.12+ installed and your environment is activated.

2. Install ZipStrain:

    ```
    pip install zipstrain
    ```

    This installs all Python dependencies **except samtools**. Install samtools separately from [htslib.org](http://www.htslib.org/download/) or via your package manager (e.g., `brew install samtools`, `apt install samtools`).

3. If you want the matrix-store workflow from a pip install, add the optional matrix dependencies:

    ```
    pip install "zipstrain[matrix]"
    ```

    Notes:

    - On Apple Silicon, the standard `torch` wheel can use Metal through the MPS backend.
    - On Linux with NVIDIA GPUs, install the ZipStrain extra and then replace Torch with the CUDA wheel that matches your system. For example:

      ```
      pip install "zipstrain[matrix]"
      pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
      ```

    - MPS requires a native macOS Python environment. Docker, Singularity, and Apptainer containers on Linux will not expose Apple Metal.

4. Verify the installation:

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

If you plan to use the matrix-store workflow from Conda, `zipstrain test` is the quickest way to confirm that the environment also has the expected matrix dependencies available.

## Docker

1. Make sure Docker is installed on your system.

2. Run ZipStrain:

    ```bash
    docker run -it parsaghadermazi/zipstrain:<version> zipstrain test
    ```

To use ZipStrain interactively or with local data, mount your data directory:

```bash
docker run -it -v /path/to/your/data:/data parsaghadermazi/zipstrain:<version> bash
```

Optional CUDA image for Linux/NVIDIA hosts:

```bash
docker run -it --gpus all parsaghadermazi/zipstrain:<version>-cuda11.8 zipstrain test
```

Available GPU image tags:

- `:<version>-cuda11.8` for older CUDA 11 era drivers
- `:<version>-cuda12.4` for newer CUDA 12.4 capable drivers

Notes:

- Docker images are published for `linux/amd64`.
- The base image includes the CPU implementation of all supported ZipStrain commands.

## Apptainer (Singularity)

Apptainer can pull directly from the Docker image:

```bash
apptainer run docker://parsaghadermazi/zipstrain:<version> zipstrain test
```

This is especially useful on HPC clusters where Docker is not available.

## Troubleshooting

- **`samtools` not found**: ZipStrain requires samtools for profiling. Make sure it is on your `PATH`. Run `samtools --version` to check.
- **Python version errors**: ZipStrain requires Python 3.12+. Check with `python --version`.
- **`zipstrain test` fails**: Ensure all dependencies are installed. If using pip, double-check that samtools is available.
