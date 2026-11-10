"""Set native worker limits before importing Polars or other compute libraries."""

from contextlib import contextmanager
import os
import sys
import threading


_duckdb_default_lock = threading.Lock()
_duckdb_default_threads: int | None = None


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


def available_cpus() -> int:
    """Return the CPU budget, or the CPUs this process may run on when there is none."""
    budget = cpu_budget()
    if budget is None:
        budget = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    return max(1, budget or 1)


def batch_share(concurrent_batches: int, batch_count: int | None = None) -> int:
    """Split :func:`available_cpus` evenly between the batches that run at the same time.

    ``batch_count`` caps the divisor, so a run with fewer batches than ``concurrent_batches``
    gives each batch more CPUs.
    """
    running = max(1, int(concurrent_batches))
    if batch_count is not None:
        running = max(1, min(running, int(batch_count)))
    return max(1, available_cpus() // running)


def duckdb_thread_budget(requested: int | None = None) -> int:
    """Bound DuckDB by the job budget, or by affinity outside a job."""
    if requested is not None and requested < 1:
        raise ValueError("DuckDB threads must be >= 1")
    budget = available_cpus()
    return min(requested, budget) if requested is not None else budget


def connect_duckdb(*args, threads: int | None = None, **kwargs):
    """Open a bounded connection without leaving DuckDB's default pool unbounded."""
    global _duckdb_default_threads

    limit = duckdb_thread_budget(threads)
    import duckdb

    # DuckDB initializes a module-level connection on import. Keep its worker
    # pool bounded even though ZipStrain uses explicit connections for queries.
    with _duckdb_default_lock:
        if _duckdb_default_threads is None or limit < _duckdb_default_threads:
            duckdb.sql(f"SET threads={limit}")
            _duckdb_default_threads = limit

    config = dict(kwargs.pop("config", {}))
    config["threads"] = limit
    return duckdb.connect(*args, config=config, **kwargs)


def _option_int(name: str) -> int | None:
    for index, arg in enumerate(sys.argv[1:], 1):
        if arg == name and index + 1 < len(sys.argv):
            return _positive_int(sys.argv[index + 1])
        if arg.startswith(name + "="):
            return _positive_int(arg.partition("=")[2])
    return None


def configure_cli_threads() -> None:
    """Bound native pools for one CLI process before the CLI imports Polars.

    ``--threads N`` on any command becomes the CPU budget (still capped by a Slurm
    allocation). Without it, ``profile-single`` falls back to ``--max-concurrency`` or every CPU.
    """
    threads = _option_int("--threads")
    if threads is not None:
        os.environ["ZIPSTRAIN_CPU_BUDGET"] = str(threads)
    budget = cpu_budget()
    profiling = "profile-single" in sys.argv[1:]
    if budget is None and profiling:
        budget = _option_int("--max-concurrency") or available_cpus()
        os.environ["ZIPSTRAIN_CPU_BUDGET"] = str(budget)
    if budget is None:
        return
    _limit_native_pools(budget, profiling)


def _limit_native_pools(budget: int, profiling: bool = False) -> None:
    """Size the Polars, Rayon, BLAS and PyArrow pools for ``budget`` CPUs."""
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


def set_threads(threads: int) -> None:
    """Limit ZipStrain to ``threads`` CPUs in this Python process; call it before any work.

    This is the library form of ``--threads``. Polars sizes its pool once, when it is first
    imported, so a larger pool that is already running cannot be shrunk: that raises
    ``RuntimeError`` instead of being silently ignored.
    """
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    polars = sys.modules.get("polars")
    if polars is not None and polars.thread_pool_size() > threads:
        raise RuntimeError(
            f"Polars is already running {polars.thread_pool_size()} threads; call "
            "zipstrain.set_threads() before importing polars or any zipstrain module that uses it.")
    os.environ["ZIPSTRAIN_CPU_BUDGET"] = str(threads)
    _limit_native_pools(cpu_budget() or threads)


_NATIVE_THREAD_VARIABLES = (
    "POLARS_MAX_THREADS", "RAYON_NUM_THREADS", "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
)


@contextmanager
def job_pool(workers: int, threads: int):
    """Yield a spawned process pool whose workers each use ``threads`` native threads.

    One worker yields ``None``: the caller runs its jobs in this process, on its own pools.
    Native pools read their size once, when a library loads, so the limits are set in the
    environment that workers inherit and restored afterwards.
    """
    if workers <= 1:
        yield None
        return
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    names = (*_NATIVE_THREAD_VARIABLES, "ZIPSTRAIN_CPU_BUDGET")
    saved = {name: os.environ.get(name) for name in names}
    os.environ.update({name: str(threads) for name in names})
    try:
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            yield pool
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def map_jobs(pool, function, jobs) -> list:
    """Run jobs in ``pool``, or here when there is no pool, and return results in order."""
    return list(pool.map(function, jobs) if pool is not None else map(function, jobs))
