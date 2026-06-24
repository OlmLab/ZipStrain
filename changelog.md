# Changelog

Entries are brief by design and describe changes relative to the previous released version.

## 0.11.1

Compared with `0.11.0`:

- Visualization: split silhouette analysis into a compute step and a plotting step.
- Visualization: added peak summarization for silhouette curves, including candidate peaks and best-peak selection.
- Visualization: kept `get_silhouette_plot` as a convenience wrapper around the new compute/plot pair.
- Matrix comparison: IBS torch runs now stay on-device until the final reduced per-pair vector is transferred back to the CPU.

## 0.11.0

Compared with `0.10.2`:

- Visualization: replaced the slow Polars pivot path in similarity-matrix preparation with direct NumPy matrix construction.
- Visualization: silhouette analysis now uses `scikit-learn` when available and warns when falling back to the manual implementation.
- Profiling: added explicit read filters for MAPQ, base quality, read ANI, and read inclusion mode.
- Profiling/reference: profiles can include `ref_base_bitmask`, and gene/genome stats can include `ref_ani` when a reference FASTA is provided.
- Matrix workflow: matrix builds now require an explicit BED+STB contract; sparse HDF5 storage is supported.
- Documentation: refreshed CLI, tutorial, installation, and Nextflow docs to match current behavior.
