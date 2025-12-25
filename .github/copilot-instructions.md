This file gives concise, repo-specific guidance to AI coding assistants working on ZipStrain.

Please keep responses actionable, reference files shown below, and avoid making assumptions about external systems (Slurm, Docker) unless the user says they exist.

- Project summary
  - ZipStrain is a Python CLI/library for profiling and comparing genomes. The Python package sources are under `zipstrain/src/zipstrain` (there is a duplicate `zipstrain/` copy in the repo; prefer editing `zipstrain/src/zipstrain/`).

- Big-picture architecture (what to read first)
  - `zipstrain/cli.py` — the CLI entrypoint (installed via `pyproject.toml` as `zipstrain`).
  - `task_manager.py` — core execution engine: Task / Batch / Runner abstractions, Engines (Docker/Apptainer/Local), and batching logic. Key classes: `Task`, `Batch`, `Runner`, `ProfileRunner`, `CompareRunner`, `GeneCompareRunner`.
    - Pattern: Task objects define a `TEMPLATE_CMD` that uses angle-bracket placeholders (e.g. `<bam-file>`). `map_io()` maps `Input`/`Output` objects into the command before execution.
    - `Engine.wrap(command, file_inputs)` should return a wrapper command (Docker/Apptainer will mount file inputs; LocalEngine returns raw command).
  - `database.py` — data model and helper classes for Profile and Comparison databases (e.g. `ProfileDatabase`, `GenomeComparisonDatabase` and pydantic models used for validation).
  - `compare.py`, `profile.py`, `utils.py`, `visualize.py` — domain logic for profile & compare computations.

- Important patterns & conventions (concrete examples)
  - Inputs/Outputs: `FileInput`, `StringInput`, `IntInput`, `FileOutput`, `BatchFileOutput` are used to validate and represent I/O. `FileInput` validates file existence and stores absolute paths.
  - Tasks: implement `TEMPLATE_CMD` and (optionally) `pre_run` to set status or change directories. Example: `CollectComps` uses `pre_run` to write `.status`.
  - Batches: `LocalBatch` builds a shell script and runs `bash`; `SlurmBatch` uses `sbatch`/`sacct` and therefore requires Slurm available on the host.
  - Runners: `Runner.run()` launches asyncio coroutines and a Rich Live UI; it registers signal handlers (SIGINT/SIGTERM). Unit tests create instances of Task/Engine directly — prefer headless runs in tests.

- Developer workflows (how to run & test locally)
  - Python requirement: project requires Python >= 3.12 (see `zipstrain/pyproject.toml`).
  - Local development (recommended):
    - python -m venv .venv
    - source .venv/bin/activate
    - cd zipstrain
    - pip install -e .
    - pip install pytest
    - pytest -q
  - Run the CLI after install: `zipstrain --help` (entrypoint in `pyproject.toml`).
  - If using Slurm features: ensure `sbatch` and `sacct` are available; SlurmBatch._check_slurm_works will call them.
  - Docker/Apptainer engines: `DockerEngine.wrap` and `ApptainerEngine.wrap` mount file paths. Do not assume the host has images unless user says so.

- Tests & CI notes
  - Tests live under `zipstrain/tests/` (run `pytest` from the `zipstrain` directory). Tests directly import `zipstrain.task_manager` and exercise Input/Engine behaviors.
  - Tests are standard pytest; no special tox/CI config detected in repository root — adapt if the project adds GitHub Actions later.

- Integration touchpoints
  - CLI -> `zipstrain.cli` -> calls into `task_manager` and domain modules.
  - `task_manager` integrates with system tools: `samtools` (used in TEMPLATE_CMDs), Slurm (`sbatch`/`sacct`), Docker and Apptainer runtimes.
  - `database.py` exposes factories used by runner helpers (e.g., `comps_db.to_complete_input_table()` used when building Compare tasks).

- Guidance for code edits
  - When editing tasks/batches, preserve the `TEMPLATE_CMD` placeholder convention and ensure `map_io()` still replaces placeholders with mapped inputs/outputs.
  - Avoid changing public method signatures in `Runner`/`Task`/`Batch` without updating corresponding tests and the CLI wiring.
  - If adding Slurm/Docker behavior, document the runtime expectation (host binary, container image) and fail fast with clear error messages when missing.

If anything here is incomplete or you'd like more details (CI commands, preferred dev folder, or examples for running a small profile/compare), tell me which area to expand and I'll iterate.
