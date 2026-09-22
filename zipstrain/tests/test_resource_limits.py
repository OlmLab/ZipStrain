import os
import subprocess
import sys

import duckdb
import pytest

from zipstrain import compare
from zipstrain import matrix_pairs
from zipstrain import resource_limits


def test_cpu_budget_uses_smaller_explicit_and_slurm_allocation(monkeypatch):
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "3")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    assert resource_limits.cpu_budget() == 3
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    assert resource_limits.cpu_budget() == 2


def test_slurm_job_without_cpu_declaration_defaults_to_one(monkeypatch):
    monkeypatch.delenv("ZIPSTRAIN_CPU_BUDGET", raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    assert resource_limits.cpu_budget() == 1


@pytest.mark.parametrize("command,pool", [
    (["utilities", "single_compare_genome"], 3),
    (["utilities", "profile-single", "--max-concurrency", "8"], 2),
    (["profile", "--num-procs", "3"], 3),
])
def test_cli_initializes_polars_with_job_budget(command, pool):
    env = os.environ.copy()
    if command[0] == "profile":
        env.pop("ZIPSTRAIN_CPU_BUDGET", None)
    else:
        env["ZIPSTRAIN_CPU_BUDGET"] = "8" if "profile-single" in command else "3"
    env.pop("SLURM_CPUS_PER_TASK", None)
    env.pop("SLURM_JOB_ID", None)
    env.pop("ZIPSTRAIN_POLARS_THREADS", None)
    env.pop("POLARS_MAX_THREADS", None)
    code = (
        "import sys; "
        f"sys.argv = {['zipstrain', *command]!r}; "
        "import zipstrain.cli; import polars as pl; import pyarrow as pa; "
        "print(pl.thread_pool_size(), pa.cpu_count(), pa.io_thread_count())"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    assert tuple(map(int, result.stdout.split())) == (pool, pool, 1)


def test_duckdb_connection_uses_job_budget(monkeypatch):
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "2")
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    conn = duckdb.connect()
    try:
        compare._duckdb_configure_connection(conn, threads=32)
        assert conn.execute("SELECT current_setting('threads')").fetchone()[0] == 2
    finally:
        conn.close()


def test_matrix_build_duckdb_connection_uses_job_budget(monkeypatch, tmp_path):
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "2")
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    conn = matrix_pairs._open_matrix_build_connection(
        tmp_path / "matrix.duckdb", duckdb_memory_limit_bytes=100_000_000,
    )
    try:
        assert conn.execute("SELECT current_setting('threads')").fetchone()[0] == 2
    finally:
        conn.close()
