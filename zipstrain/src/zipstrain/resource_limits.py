"""Set native worker limits before importing Polars or other compute libraries."""

import os
import sys


def _positive_int(value: str | None) -> int | None:
    try:
        parsed = int(value) if value is not None else None
    except ValueError:
        return None
    return parsed if parsed is not None and parsed > 0 else None


def cpu_budget() -> int | None:
    """Return the requested allocation, never the host's reported CPU count."""
    explicit = _positive_int(os.environ.get("ZIPSTRAIN_CPU_BUDGET"))
    slurm = _positive_int(os.environ.get("SLURM_CPUS_PER_TASK"))
    if explicit is not None and slurm is not None:
        return min(explicit, slurm)
    if explicit is not None or slurm is not None:
        return explicit or slurm
    # A Slurm job without --cpus-per-task has one CPU per task by default.
    return 1 if os.environ.get("SLURM_JOB_ID") else None


def _option_int(name: str) -> int | None:
    for index, arg in enumerate(sys.argv[1:], 1):
        if arg == name and index + 1 < len(sys.argv):
            return _positive_int(sys.argv[index + 1])
        if arg.startswith(name + "="):
            return _positive_int(arg.partition("=")[2])
    return None


def configure_cli_threads() -> None:
    """Bound native pools for one CLI process before the CLI imports Polars."""
    budget = cpu_budget()
    profiling = "profile-single" in sys.argv[1:]
    profile_runner = sys.argv[1:2] == ["profile"]
    if budget is None and profiling:
        budget = _option_int("--max-concurrency") or 4
        os.environ["ZIPSTRAIN_CPU_BUDGET"] = str(budget)
    if budget is None and profile_runner:
        budget = _option_int("--num-procs") or _option_int("-n") or 8
        os.environ["ZIPSTRAIN_CPU_BUDGET"] = str(budget)
    if budget is None:
        return

    requested_pool = _positive_int(os.environ.get("ZIPSTRAIN_POLARS_THREADS"))
    if requested_pool is None:
        requested_pool = max(1, budget // 4) if profiling else budget
    pool = min(budget, requested_pool)
    for name in (
        "POLARS_MAX_THREADS", "RAYON_NUM_THREADS", "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        existing = _positive_int(os.environ.get(name))
        os.environ[name] = str(min(pool, existing) if existing is not None else pool)

    import pyarrow as pa

    pa.set_cpu_count(pool)
    pa.set_io_thread_count(1)
