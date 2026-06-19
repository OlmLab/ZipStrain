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

Matrix-store workflow dependencies:

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

## Cite

If you use ZipStrain in your research, please cite the preprint:

Ghadermazi P, Emerson JB, Olm MR. 2026. *ZipStrain Enables Rapid and Precise
Strain-Resolved Metagenomics*. bioRxiv. DOI:
`10.64898/2026.05.20.726564`

GitHub citation metadata is provided in `CITATION.cff`:

- https://github.com/OlmLab/ZipStrain/blob/main/CITATION.cff

## License

ZipStrain is distributed under the MIT License.
