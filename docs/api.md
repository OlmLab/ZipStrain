# API

ZipStrain provides a Python API for building workflows programmatically.
Most users will interact with the matrix workflow through the CLI, but the underlying matrix-store helpers also live in Python.

## Core Modules

### Database

:::zipstrain.database

### Profile

:::zipstrain.profile

### Compare

:::zipstrain.compare

### Utils

:::zipstrain.utils

### Visualize

:::zipstrain.visualize

## Matrix Store Workflow

The current matrix workflow is implemented in `zipstrain.matrix_pairs`.
The main public entry points are:

- `parse_matrix_calculations(...)`
- `build_matrix_hdf5(...)`
- `append_matrix_hdf5(...)`
- `export_matrix_db_hdf5(...)`
- `matrix_compare(...)`
- `export_matrix_compare_parquet(...)`

These power the CLI commands:

- `zipstrain utilities build-matrix-db`
- `zipstrain utilities append-matrix-db`
- `zipstrain utilities matrix-db-to-hdf5`
- `zipstrain utilities matrix-compare`
- `zipstrain utilities matrix-compare-export`

`build_matrix_hdf5(...)` defaults to the compact popANI bitmask layout
(`A=1`, `T=2`, `C=4`, `G=8`). Use `storage_mode="counts"` with
`count_dtype="auto"` to retain `A,T,C,G` counts for conANI or cosANI.
`min_cov` is applied at build time and persisted in the store contract.
`matrix_compare(..., ani_method=...)` accepts `popani`, `conani`, and
`cosani_<threshold>`; the latter two require count storage.
`matrix_compare(..., compare_kernel="block", calculate="genome_ani")` selects the
blocked matrix-multiplication kernel in `zipstrain.matrix_gemm` for
genome-level popANI/conANI; results are identical to the default
`compare_kernel="anchor"`. The old names `gemm` and `pairwise` remain
accepted as deprecated aliases for `block` and `anchor`, respectively.

`calculate` accepts `max_identical_stretch`, `gene_ani`, `genome_ani`,
`all`, or `+`/comma combinations. Legacy spellings remain accepted; stored
comparison contracts and output columns are unchanged.

`matrix_workflow.run_matrix_compare(...)` also accepts `calculate`,
`compare_kernel`, and `compare_progress_callback`, exposing the same metric
selection and kernels through the build/append/compare/export workflow.

For usage patterns and examples, see the [Tutorial](./Tutorial.md) and the [CLI reference](./usermanual.md#zipstrain-command-line-interface).

## Pair analyses

`zipstrain.compare.PairComparison` validates two V2 profiles once and joins them on their
shared positions. Its analyses return lazy Polars frames, and `PairComparison.collect` runs
several of them over that one join:

```python
from zipstrain.compare import PairComparison

pair = PairComparison("A_profile.parquet", "B_profile.parquet", ani_method="popani")
table = pair.compare("genome_ani+max_identical_stretch").collect()
blocks, stretches = pair.collect(pair.mosaic_regions(), pair.identical_stretches(min_length=1_000))
```

- `compare(calculate)` gives the same rows and columns as `zipstrain compare` for one pair.
  Pass `preprofile_dir=` to the constructor for gene metrics.
- `identical_stretches(min_length)` lists runs of matching shared positions, each with a
  p-value for its length.
- `mosaic_regions(factor=20, switch_probability=1e-3)` finds blocks of unusually low or high
  divergence: candidate recombined and imported segments.

The [recombination mosaics guide](parameter_selection/recombination/index.md) shows
`identical_stretches` and `mosaic_regions` on simulated data and how to choose their
parameters. The [V2 guide](./v2.md) documents the annotation, observed-codon,
gene-statistics, and explicit-input comparison functions, their schemas, selection formulas,
and corresponding tests.

## Task Manager

For more advanced users, custom workflows can also be assembled through the task manager layer:

:::zipstrain.task_manager
