import os
import subprocess
import sys

import duckdb
import polars as pl
import pytest

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
    (["utilities", "compare-pair"], 3),
    (["utilities", "profile-single", "--max-concurrency", "8"], 2),
    (["profile", "--threads", "3"], 3),
    (["compare", "--threads", "2"], 2),
    (["utilities", "compare-pair", "--threads=5"], 5),
    (["utilities", "profile-single", "--threads", "8"], 2),
])
def test_cli_initializes_polars_with_job_budget(command, pool):
    env = os.environ.copy()
    if "--threads" in " ".join(command):
        env["ZIPSTRAIN_CPU_BUDGET"] = "64"  # --threads overrides an inherited budget
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


def test_polars_cli_import_does_not_initialize_duckdb():
    env = os.environ.copy()
    env["ZIPSTRAIN_CPU_BUDGET"] = "2"
    code = (
        "import sys, psutil; "
        "before = psutil.Process().num_threads(); "
        "import zipstrain.cli; "
        "after = psutil.Process().num_threads(); "
        "print(int('duckdb' in sys.modules), before, after)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], env=env, text=True, capture_output=True, check=True,
    )
    loaded, before, after = map(int, result.stdout.split())
    assert loaded == 0
    assert after < before + 8


def test_polars_pair_comparison_does_not_initialize_duckdb(tmp_path):
    frame = pl.DataFrame({
        "chrom": ["chr1"], "pos": [1], "gene": ["NA"], "genome": ["genome1"],
        "A": [10], "C": [0], "G": [0], "T": [0],
    })
    profile_1 = tmp_path / "a_profile.parquet"
    profile_2 = tmp_path / "b_profile.parquet"
    output = tmp_path / "comparison.parquet"
    frame.write_parquet(profile_1)
    frame.write_parquet(profile_2)
    args = [
        "utilities", "compare-pair", "--profile-location-1", str(profile_1),
        "--profile-location-2", str(profile_2),
        "--calculate", "genome_ani", "--output-file", str(output),
    ]
    code = (
        "import sys; from zipstrain.cli import cli; "
        f"cli.main(args={args!r}, standalone_mode=False); "
        "print('duckdb_loaded=' + str(int('duckdb' in sys.modules)))"
    )
    env = os.environ.copy()
    env["ZIPSTRAIN_CPU_BUDGET"] = "2"
    result = subprocess.run(
        [sys.executable, "-c", code], env=env, text=True, capture_output=True, check=True,
    )
    assert "duckdb_loaded=0" in result.stdout
    assert pl.read_parquet(output).get_column("genome_ani").to_list() == [100.0]


def test_duckdb_connection_and_default_pool_are_bounded_at_creation():
    env = os.environ.copy()
    env["ZIPSTRAIN_CPU_BUDGET"] = "2"
    code = (
        "import psutil; "
        "before = psutil.Process().num_threads(); "
        "from zipstrain.resource_limits import connect_duckdb; "
        "conn = connect_duckdb(threads=32); "
        "import duckdb; "
        "active = psutil.Process().num_threads(); "
        "print(before, active, "
        "conn.execute(\"SELECT current_setting('threads')\").fetchone()[0], "
        "duckdb.sql(\"SELECT current_setting('threads')\").fetchone()[0])"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], env=env, text=True, capture_output=True, check=True,
    )
    before, active, explicit_threads, default_threads = map(int, result.stdout.split())
    assert (explicit_threads, default_threads) == (2, 2)
    assert active < before + 8


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


def test_batch_share_splits_budget_between_running_batches(monkeypatch):
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "16")
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    assert resource_limits.batch_share(5) == 3
    assert resource_limits.batch_share(5, batch_count=2) == 8
    assert resource_limits.batch_share(64) == 1


def test_set_threads_refuses_a_larger_running_polars_pool(monkeypatch):
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    if pl.thread_pool_size() > 1:
        with pytest.raises(RuntimeError, match="set_threads"):
            resource_limits.set_threads(1)
    with pytest.raises(ValueError):
        resource_limits.set_threads(0)


def test_set_threads_before_polars_bounds_its_pool():
    env = {key: value for key, value in os.environ.items()
           if key not in ("ZIPSTRAIN_CPU_BUDGET", "POLARS_MAX_THREADS", "SLURM_CPUS_PER_TASK", "SLURM_JOB_ID")}
    code = "import zipstrain; zipstrain.set_threads(2); import polars as pl; print(pl.thread_pool_size())"
    result = subprocess.run([sys.executable, "-c", code], env=env, text=True, capture_output=True, check=True)
    assert result.stdout.strip() == "2"


def test_requested_pair_workers_never_exceed_the_budget(monkeypatch):
    from zipstrain import compare
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "4")
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    assert compare._pair_workers(16, 100, 64.0, 1000) == (4, 1)
