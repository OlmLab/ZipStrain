# ZipStrain

ZipStrain is a strain-resolution metagenomics toolkit for profiling mapped reads
into nucleotide-count tables and comparing samples at genome and gene
resolution.

Documentation:

- Docs: https://OlmLab.github.io/ZipStrain/
- Repository: https://github.com/OlmLab/ZipStrain

## Install

Base install:

```bash
pip install zipstrain
```

Experimental matrix utilities with Torch-backed execution:

```bash
pip install "zipstrain[matrix]"
```

Notes:

- Apple Silicon can use the standard `torch` wheel with the MPS backend.
- Linux CUDA installs should replace Torch with the matching CUDA wheel from
  PyTorch:

  ```bash
  pip install "zipstrain[matrix]"
  pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
  ```

- Profiling requires `samtools` to be installed separately when using pip.

## Verify

```bash
zipstrain --version
zipstrain test
```
