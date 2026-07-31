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
The high-level `zipstrain compare --method matrix --ani-method ...` command also
accepts a comma-separated list. It reuses one count store, keeps one resumable
compare database per method, and joins the method result columns for export.

For usage patterns and examples, see the [Tutorial](./Tutorial.md) and the [CLI reference](./usermanual.md#zipstrain-command-line-interface).

## Task Manager

For more advanced users, custom workflows can also be assembled through the task manager layer:

:::zipstrain.task_manager
