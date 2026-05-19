#!/usr/bin/env python3
"""Profile zipstrain matrix compare on the local machine.

This script mirrors the main matrix-compare CLI arguments, wraps the compare
call in cProfile, and prints profiling summaries even if the run is cancelled
with Ctrl+C. It also collects lightweight wall-clock timings for key
zipstrain.matrix_pairs functions inside the current process.

Notes:
- cProfile only sees Python execution in the current process.
- The wall-clock wrappers also only see the current process.
- If you use `--loader-executor process` or `--writer-executor process`,
  subprocess work will not appear as per-function timings in this report.
  For the most complete host-side attribution, prefer thread executors while
  profiling.
"""

from __future__ import annotations

import argparse
import cProfile
import io
import pstats
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_SRC = REPO_ROOT / "zipstrain" / "src"
if str(PACKAGE_SRC) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SRC))

from zipstrain import matrix_pairs as mp  # noqa: E402


@dataclass
class FunctionWallStat:
    calls: int = 0
    wall_seconds: float = 0.0


class WallTimerRegistry:
    def __init__(self) -> None:
        self._stats: dict[str, FunctionWallStat] = {}
        self._lock = threading.Lock()

    def record(self, name: str, seconds: float) -> None:
        with self._lock:
            stat = self._stats.setdefault(name, FunctionWallStat())
            stat.calls += 1
            stat.wall_seconds += seconds

    def snapshot(self) -> list[tuple[str, FunctionWallStat]]:
        with self._lock:
            return sorted(
                ((name, FunctionWallStat(stat.calls, stat.wall_seconds)) for name, stat in self._stats.items()),
                key=lambda item: item[1].wall_seconds,
                reverse=True,
            )


class MonkeyPatchWallTimers:
    def __init__(self, module, function_names: list[str], registry: WallTimerRegistry) -> None:
        self.module = module
        self.function_names = function_names
        self.registry = registry
        self.originals: dict[str, Callable] = {}

    def __enter__(self):
        for name in self.function_names:
            if not hasattr(self.module, name):
                continue
            original = getattr(self.module, name)
            self.originals[name] = original

            def make_wrapper(func_name: str, func: Callable):
                def wrapper(*args, **kwargs):
                    started = time.perf_counter()
                    try:
                        return func(*args, **kwargs)
                    finally:
                        self.registry.record(func_name, time.perf_counter() - started)

                wrapper.__name__ = getattr(func, "__name__", func_name)
                wrapper.__doc__ = getattr(func, "__doc__", None)
                wrapper.__qualname__ = getattr(func, "__qualname__", func_name)
                return wrapper

            setattr(self.module, name, make_wrapper(name, original))
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        for name, original in self.originals.items():
            setattr(self.module, name, original)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Profile zipstrain matrix compare and print useful host-side timings on completion or Ctrl+C.",
    )
    parser.add_argument("--matrix-db-file", "-m", required=True, type=Path, help="Input matrix store (HDF5-backed current format).")
    parser.add_argument("--output-file", "-o", required=True, type=Path, help="Output DuckDB compare database.")
    parser.add_argument("--genome", "-g", default="all", help="Optional genome scope.")
    parser.add_argument("--memory-limit-gb", type=float, default=16.0, help="Approximate memory budget for compare.")
    parser.add_argument(
        "--anchor-queue-size",
        type=int,
        default=1,
        help="Host-side torch anchor queue size. Only one anchor is transferred to the GPU at a time.",
    )
    parser.add_argument(
        "--target-queue-size",
        type=int,
        default=1,
        help="Host-side torch target queue size. 1 keeps the current synchronous target-load behavior.",
    )
    parser.add_argument(
        "--result-transfer-batch-size",
        type=int,
        default=mp.MATRIX_COMPARE_RESULT_TRANSFER_BATCH_SIZE_DEFAULT,
        help="Number of torch compare units to batch before transferring result vectors back to CPU.",
    )
    parser.add_argument(
        "--loader-executor",
        choices=mp.MATRIX_IO_EXECUTOR_KINDS,
        default="thread",
        help="Executor used for torch loader prefetch work.",
    )
    parser.add_argument(
        "--writer-executor",
        choices=mp.MATRIX_IO_EXECUTOR_KINDS,
        default="thread",
        help="Executor used for torch result writing/checkpoint work.",
    )
    parser.add_argument(
        "--calculate",
        default="all",
        help="Matrix metrics to compute: ani or ani+ibs. 'all' currently means ani+ibs.",
    )
    parser.add_argument(
        "--backend",
        choices=mp.MATRIX_PAIR_BACKENDS,
        default="numpy",
        help="Compute backend. Torch backends use CPU/CUDA/MPS depending on selection.",
    )
    parser.add_argument(
        "--sort",
        choices=("cumulative", "time", "calls"),
        default="cumulative",
        help="Sort order for the cProfile report.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=120,
        help="How many rows to print from each profiling table.",
    )
    parser.add_argument(
        "--profile-output",
        type=Path,
        default=None,
        help="Optional path to save the raw cProfile stats file.",
    )
    parser.add_argument(
        "--zipstrain-only",
        action="store_true",
        help="Only print cProfile rows whose file path contains 'zipstrain'.",
    )
    return parser


def format_seconds(seconds: float) -> str:
    return f"{seconds:.3f}s"


def print_header(title: str) -> None:
    print(f"\n{title}")
    print("-" * len(title))


def render_cprofile_report(
    profiler: cProfile.Profile,
    *,
    limit: int,
    sort_key: str,
    zipstrain_only: bool,
) -> str:
    stats = pstats.Stats(profiler)
    stats.sort_stats(sort_key)
    rows: list[tuple[tuple[str, int, str], tuple[int, int, float, float, dict]]] = list(stats.stats.items())
    if zipstrain_only:
        rows = [item for item in rows if "zipstrain" in item[0][0]]
    if sort_key == "cumulative":
        rows.sort(key=lambda item: item[1][3], reverse=True)
    elif sort_key == "time":
        rows.sort(key=lambda item: item[1][2], reverse=True)
    else:
        rows.sort(key=lambda item: item[1][1], reverse=True)
    stream = io.StringIO()
    stream.write(f"{'ncalls':>12} {'tottime':>10} {'cumtime':>10} function\n")
    stream.write(f"{'-'*12} {'-'*10} {'-'*10} {'-'*50}\n")
    for (filename, lineno, funcname), (cc, nc, tt, ct, _callers) in rows[:limit]:
        function_label = f"{Path(filename).name}:{lineno}:{funcname}"
        stream.write(f"{nc:>12} {tt:>10.4f} {ct:>10.4f} {function_label}\n")
    return stream.getvalue().rstrip()


def print_wall_timer_report(registry: WallTimerRegistry, overall_wall_seconds: float, limit: int) -> None:
    snapshot = registry.snapshot()
    print_header("Wall-Time Summary (Current Process)")
    print(f"overall_wall_time={format_seconds(overall_wall_seconds)}")
    if not snapshot:
        print("No function wall timings were captured.")
        return
    print(f"{'function':<48} {'calls':>8} {'wall_s':>12} {'share':>9}")
    print(f"{'-'*48} {'-'*8} {'-'*12} {'-'*9}")
    for name, stat in snapshot[:limit]:
        share = (stat.wall_seconds / overall_wall_seconds * 100.0) if overall_wall_seconds > 0 else 0.0
        print(f"{name:<48} {stat.calls:>8} {stat.wall_seconds:>12.4f} {share:>8.1f}%")


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.loader_executor == "process" or args.writer_executor == "process":
        print(
            "note: process executors are enabled; subprocess work will not appear in per-function timings in this report.",
            file=sys.stderr,
        )

    wall_registry = WallTimerRegistry()
    wrapped_functions = [
        "_load_matrix_compare_resume_state",
        "_plan_chunk_sizes",
        "_load_target_queue_block_for_hdf5_torch",
        "_load_target_prefetch_unit_for_hdf5_torch",
        "_load_anchor_queue_batch_for_hdf5_torch",
        "_prepare_torch_matrix",
        "_compare_tile_presence_torch_tensors",
        "_compare_tile_presence_torch_tensors_with_mask",
        "_compare_anchor_against_target_chunk_torch_device",
        "_accumulate_gene_counts_from_full_torch_masks",
        "_max_ibs_from_torch_mask",
        "_update_ibs_numpy",
        "_download_torch_result_tensor_batch",
        "_make_arrow_table_from_compare_payload",
        "_write_matrix_compare_payload_batch",
    ]

    profiler = cProfile.Profile()
    started = time.perf_counter()
    interrupted = False
    summary = None
    error: Optional[BaseException] = None

    with MonkeyPatchWallTimers(mp, wrapped_functions, wall_registry):
        profiler.enable()
        try:
            summary = mp.matrix_compare(
                matrix_db_file=args.matrix_db_file.resolve(),
                output_file=args.output_file.resolve(),
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome=args.genome,
                memory_limit_gb=args.memory_limit_gb,
                anchor_queue_size=args.anchor_queue_size,
                target_queue_size=args.target_queue_size,
                result_transfer_batch_size=args.result_transfer_batch_size,
                loader_executor_kind=args.loader_executor,
                writer_executor_kind=args.writer_executor,
                backend=args.backend,
                calculate=args.calculate,
                emit_writer_logs=False,
                progress_callback=None,
            )
        except KeyboardInterrupt as exc:
            interrupted = True
            error = exc
        except BaseException as exc:  # noqa: BLE001
            error = exc
        finally:
            profiler.disable()

    ended = time.perf_counter()
    total_wall = ended - started

    if args.profile_output is not None:
        args.profile_output.parent.mkdir(parents=True, exist_ok=True)
        profiler.dump_stats(str(args.profile_output))

    print_header("Run Summary")
    print(f"interrupted={interrupted}")
    print(f"overall_wall_time={format_seconds(total_wall)}")
    print(f"overall_profiled_python_time={format_seconds(pstats.Stats(profiler).total_tt)}")
    if args.profile_output is not None:
        print(f"profile_output={args.profile_output}")
    if summary is not None:
        print(
            " ".join(
                [
                    f"requested_pairs={summary.requested_pairs}",
                    f"rows={summary.written_rows}",
                    f"anchor_groups={summary.anchor_groups}",
                    f"target_chunks={summary.target_chunks}",
                    f"backend={summary.backend}",
                    f"device={summary.device}",
                ]
            )
        )

    print_wall_timer_report(wall_registry, overall_wall_seconds=total_wall, limit=args.limit)

    print_header(f"cProfile ({'zipstrain only' if args.zipstrain_only else 'all functions'})")
    print(render_cprofile_report(profiler, limit=args.limit, sort_key=args.sort, zipstrain_only=args.zipstrain_only))

    if interrupted:
        print("\nInterrupted by user. Profiling results above include work completed before Ctrl+C.", file=sys.stderr)
        return 130

    if error is not None:
        print(f"\nRun failed: {error!r}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
