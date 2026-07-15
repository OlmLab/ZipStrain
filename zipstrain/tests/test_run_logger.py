import json

import pytest

from zipstrain.run_logger import (
    RUN_LOG_NAME,
    STATUS_JSON_NAME,
    STATUS_ABORTED,
    STATUS_COMPLETED,
    STATUS_CRASHED,
    RunLogger,
)


def _status(out_dir):
    return json.loads((out_dir / STATUS_JSON_NAME).read_text())


def test_completed_run_records_steps_and_status(tmp_path):
    out = tmp_path / "out"
    with RunLogger(out, command="demo", argv=["zipstrain", "demo"]) as log:
        log.step("first")
        log.step("second")

    status = _status(out)
    assert status["status"] == STATUS_COMPLETED
    assert status["steps_completed"] == 2
    assert status["last_step"] == "second"
    assert status["command"] == "demo"
    assert "finished_at" in status

    log_text = (out / RUN_LOG_NAME).read_text()
    assert "RUNNING" in log_text
    assert "STEP 1: first" in log_text
    assert "STEP 2: second" in log_text
    assert "COMPLETED" in log_text


def test_crashed_run_records_last_step_and_traceback(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(RuntimeError, match="boom"):
        with RunLogger(out, command="demo") as log:
            log.step("only step")
            raise RuntimeError("boom")

    status = _status(out)
    assert status["status"] == STATUS_CRASHED
    assert status["failed_during_step"] == "only step"
    assert status["error_type"] == "RuntimeError"
    assert status["error_message"] == "boom"

    log_text = (out / RUN_LOG_NAME).read_text()
    assert "CRASHED during: only step" in log_text
    assert "RuntimeError: boom" in log_text
    assert "Traceback (most recent call last)" in log_text


def test_keyboard_interrupt_marks_aborted(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(KeyboardInterrupt):
        with RunLogger(out, command="demo") as log:
            log.step("interrupted here")
            raise KeyboardInterrupt()

    status = _status(out)
    assert status["status"] == STATUS_ABORTED
    assert status["failed_during_step"] == "interrupted here"
    assert "ABORTED during: interrupted here" in (out / RUN_LOG_NAME).read_text()


def test_as_callback_calls_wrapped_then_logs(tmp_path):
    out = tmp_path / "out"
    seen = []
    with RunLogger(out, command="demo") as log:
        cb = log.as_callback(seen.append)
        cb("step via callback")

    assert seen == ["step via callback"]
    status = _status(out)
    assert status["last_step"] == "step via callback"
    assert status["steps_completed"] == 1


def test_reruns_append_rather_than_truncate(tmp_path):
    out = tmp_path / "out"
    with RunLogger(out, command="demo") as log:
        log.step("run one")
    with RunLogger(out, command="demo") as log:
        log.step("run two")

    log_text = (out / RUN_LOG_NAME).read_text()
    assert "run one" in log_text and "run two" in log_text
    # The second run's status overwrites the first (latest state wins).
    assert _status(out)["last_step"] == "run two"
