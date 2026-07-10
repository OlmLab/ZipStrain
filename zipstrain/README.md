# ZipStrain

Fast strain-level metagenomics in three commands: **`map`** reads to BAMs,
**`profile`** them into per-position nucleotide-count tables, and **`compare`**
samples by ANI to tell whether they share a strain.

Documentation:

- Docs: https://OlmLab.github.io/ZipStrain/
- Repository: https://github.com/OlmLab/ZipStrain

## Install

Conda is the easiest path — it brings in `samtools` and the `map` aligners too:

```bash
conda create -n zipstrain -c conda-forge -c bioconda \
  python=3.12 zipstrain bowtie2 samtools sylph
```

Or with pip (into a fresh virtual environment):

```bash
pip install zipstrain
```

Matrix-store comparison dependencies (`compare --method matrix`):

```bash
pip install "zipstrain[matrix]"
```

Notes:

- With pip, install `samtools` separately (profiling needs it); `zipstrain map`
  additionally needs `bowtie2` and `sylph` (and `prodigal` for `--predict-genes`).
- On Apple Silicon, use a native `osx-arm64` Conda so dependencies install
  natively; the standard `torch` wheel uses the MPS backend.
- Linux CUDA installs should replace Torch with the matching CUDA wheel from
  PyTorch:

  ```bash
  pip install "zipstrain[matrix]"
  pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
  ```

See the [installation guide](https://OlmLab.github.io/ZipStrain/installation/)
for full details.

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
