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

For usage patterns and examples, see the [Tutorial](./Tutorial.md) and [CLI Reference](./cli.md).

## Task Manager

For more advanced users, custom workflows can also be assembled through the task manager layer:

:::zipstrain.task_manager
