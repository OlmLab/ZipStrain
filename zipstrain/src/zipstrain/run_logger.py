"""zipstrain.run_logger
========================
A small, universal run logger shared by the top-level ZipStrain commands
(``map``, ``profile``, ``compare``).

Each run writes two files into its output directory so a user (or a script) can
tell, at any moment, whether the program is still running, finished cleanly, or
crashed — and, if it crashed, what it was doing at the time:

* ``zipstrain_run.log``    — a human-readable, append-only timeline of the run.
* ``zipstrain_status.json`` — a machine-readable snapshot of the current status,
  the last step that started, the process id, timestamps, and (on failure) the
  error and traceback.

The logger is deliberately defensive: logging must never be the reason a command
fails, so all file I/O is wrapped and degrades to a warning if something goes
wrong.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import pathlib
import traceback
from typing import Any, Callable

from zipstrain import __version__

RUN_LOG_NAME = "zipstrain_run.log"
STATUS_JSON_NAME = "zipstrain_status.json"

# Status values written to zipstrain_status.json["status"].
STATUS_RUNNING = "running"
STATUS_COMPLETED = "completed"
STATUS_CRASHED = "crashed"
STATUS_ABORTED = "aborted"  # user interrupt (Ctrl-C) or click.Abort


def _now() -> _dt.datetime:
    return _dt.datetime.now().astimezone()


def _ts(moment: _dt.datetime) -> str:
    return moment.isoformat(timespec="seconds")


class RunLogger:
    """Context manager that records a command's lifecycle to the output dir.

    Usage::

        with RunLogger(out_dir, command="map", argv=sys.argv) as log:
            log.step("Building Bowtie2 index")
            ...

    Entering writes a header and marks the run ``running``. Leaving normally
    marks it ``completed``; leaving via an exception marks it ``crashed`` (or
    ``aborted`` for ``KeyboardInterrupt``/``click.Abort``), records the last step
    and the traceback, and then re-raises so the CLI behaves exactly as before.
    """

    def __init__(
        self,
        output_dir: str | pathlib.Path,
        *,
        command: str,
        argv: list[str] | None = None,
    ) -> None:
        self.output_dir = pathlib.Path(output_dir)
        self.command = command
        self.argv = list(argv) if argv is not None else None
        self.pid = os.getpid()
        self.started_at = _now()
        self.updated_at = self.started_at
        self.last_step: str | None = None
        self.steps_completed = 0
        self.status = STATUS_RUNNING
        self._disabled = False

    # -- paths -------------------------------------------------------------
    @property
    def log_path(self) -> pathlib.Path:
        return self.output_dir / RUN_LOG_NAME

    @property
    def status_path(self) -> pathlib.Path:
        return self.output_dir / STATUS_JSON_NAME

    # -- low-level, never-raising writers ----------------------------------
    def _append_log(self, line: str) -> None:
        if self._disabled:
            return
        try:
            with open(self.log_path, "a", encoding="utf-8") as handle:
                handle.write(line + "\n")
        except OSError as exc:  # pragma: no cover - defensive
            self._disabled = True
            print(f"Warning: could not write run log ({exc}); continuing without it.")

    def _write_status(self, extra: dict[str, Any] | None = None) -> None:
        if self._disabled:
            return
        payload: dict[str, Any] = {
            "command": self.command,
            "status": self.status,
            "pid": self.pid,
            "zipstrain_version": __version__,
            "started_at": _ts(self.started_at),
            "updated_at": _ts(self.updated_at),
            "elapsed_seconds": round((self.updated_at - self.started_at).total_seconds(), 1),
            "steps_completed": self.steps_completed,
            "last_step": self.last_step,
            "output_dir": str(self.output_dir.absolute()),
            "cwd": os.getcwd(),
        }
        if self.argv is not None:
            payload["argv"] = self.argv
        if extra:
            payload.update(extra)
        try:
            tmp = self.status_path.with_suffix(self.status_path.suffix + ".part")
            with open(tmp, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
                handle.write("\n")
            os.replace(tmp, self.status_path)
        except OSError as exc:  # pragma: no cover - defensive
            self._disabled = True
            print(f"Warning: could not write run status ({exc}); continuing without it.")

    # -- public API --------------------------------------------------------
    def start(self) -> "RunLogger":
        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:  # pragma: no cover - defensive
            self._disabled = True
            print(f"Warning: could not create output dir for run log ({exc}).")
            return self
        self._append_log("=" * 72)
        self._append_log(f"ZipStrain {__version__} — `zipstrain {self.command}`")
        self._append_log(f"Started: {_ts(self.started_at)}  (pid {self.pid})")
        if self.argv is not None:
            self._append_log("Command: " + " ".join(self.argv))
        self._append_log(f"Working dir: {os.getcwd()}")
        self._append_log("-" * 72)
        self._append_log(f"[{_ts(self.started_at)}] RUNNING")
        self._write_status()
        return self

    def step(self, message: str) -> None:
        """Record that a new step has started; updates the 'last step' status."""
        self.updated_at = _now()
        self.last_step = message
        self.steps_completed += 1
        self._append_log(f"[{_ts(self.updated_at)}] STEP {self.steps_completed}: {message}")
        self._write_status()

    def note(self, message: str) -> None:
        """Record an informational line without advancing the step counter."""
        self.updated_at = _now()
        self._append_log(f"[{_ts(self.updated_at)}] NOTE: {message}")
        self._write_status()

    def as_callback(self, wrapped: Callable[[str], None] | None = None) -> Callable[[str], None]:
        """Return a progress callback that logs each message via :meth:`step`.

        If ``wrapped`` is given (e.g. the console printer), it is called first so
        the terminal output is unchanged; the message is then logged.
        """

        def _callback(message: str) -> None:
            if wrapped is not None:
                wrapped(message)
            self.step(message)

        return _callback

    def complete(self) -> None:
        self.status = STATUS_COMPLETED
        self.updated_at = _now()
        elapsed = (self.updated_at - self.started_at).total_seconds()
        self._append_log(f"[{_ts(self.updated_at)}] COMPLETED (elapsed {elapsed:.1f}s)")
        self._write_status(extra={"finished_at": _ts(self.updated_at)})

    def fail(self, exc: BaseException, *, aborted: bool = False) -> None:
        self.status = STATUS_ABORTED if aborted else STATUS_CRASHED
        self.updated_at = _now()
        label = "ABORTED" if aborted else "CRASHED"
        during = f" during: {self.last_step}" if self.last_step else " before the first step"
        self._append_log(f"[{_ts(self.updated_at)}] {label}{during}")
        self._append_log(f"    {type(exc).__name__}: {exc}")
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        for tb_line in tb.rstrip().splitlines():
            self._append_log("    " + tb_line)
        self._write_status(
            extra={
                "finished_at": _ts(self.updated_at),
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "failed_during_step": self.last_step,
            }
        )

    # -- context manager ---------------------------------------------------
    def __enter__(self) -> "RunLogger":
        return self.start()

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc is None:
            self.complete()
            return False
        # KeyboardInterrupt is not a subclass of Exception; click raises Abort.
        aborted = isinstance(exc, KeyboardInterrupt) or type(exc).__name__ == "Abort"
        self.fail(exc, aborted=aborted)
        return False  # propagate the exception
