from pathlib import Path
from types import SimpleNamespace

import pytest

from zipstrain_rust_profiler import RustProfileBackend
import zipstrain_rust_profiler


def _request(tmp_path, **overrides):
    values = dict(
        bam_file=tmp_path / "sample.bam",
        bed_file=tmp_path / "regions.bed",
        stb_file=tmp_path / "mapping.stb",
        null_model=tmp_path / "null.parquet",
        output_dir=tmp_path / "out",
        max_concurrency=3,
        num_chunks=8,
        min_mapq=10,
        min_baseq=20,
        min_freq=0.02,
        read_inclusion="paired",
        reference_fasta=None,
        gene_range_table=None,
        profiling_contract=None,
        min_read_ani=None,
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def test_adapter_passes_settings_without_optional_inputs(tmp_path, monkeypatch):
    executable = tmp_path / "zipstrain-rust-profiler"
    executable.touch()
    monkeypatch.setattr(zipstrain_rust_profiler.sysconfig, "get_path", lambda key: str(tmp_path))
    calls = []
    monkeypatch.setattr(zipstrain_rust_profiler.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs)))
    RustProfileBackend().run(_request(tmp_path))
    command = calls[0][0][0]
    assert command[0] == str(executable)
    assert command[command.index("--threads") + 1] == "3"
    assert command[command.index("--min-freq") + 1] == "0.02"
    assert "--reference" not in command
    assert "--min-read-ani" not in command
    assert calls[0][1]["check"] is True


def test_adapter_passes_optional_inputs(tmp_path, monkeypatch):
    (tmp_path / "zipstrain-rust-profiler").touch()
    monkeypatch.setattr(zipstrain_rust_profiler.sysconfig, "get_path", lambda key: str(tmp_path))
    commands = []
    monkeypatch.setattr(zipstrain_rust_profiler.subprocess, "run", lambda command, **kwargs: commands.append(command))
    request = _request(
        tmp_path,
        reference_fasta=tmp_path / "reference.fasta",
        gene_range_table=tmp_path / "genes.tsv",
        profiling_contract=tmp_path / "contract.json",
        min_read_ani=0.95,
    )
    RustProfileBackend().run(request)
    command = commands[0]
    for flag, expected in (
        ("--reference", request.reference_fasta),
        ("--gene-range", request.gene_range_table),
        ("--profiling-contract", request.profiling_contract),
        ("--min-read-ani", 0.95),
    ):
        assert command[command.index(flag) + 1] == str(expected)


def test_adapter_respects_cpu_budget(tmp_path, monkeypatch):
    (tmp_path / "zipstrain-rust-profiler").touch()
    monkeypatch.setattr(zipstrain_rust_profiler.sysconfig, "get_path", lambda key: str(tmp_path))
    monkeypatch.setattr(zipstrain_rust_profiler, "cpu_budget", lambda: 2)
    commands = []
    monkeypatch.setattr(zipstrain_rust_profiler.subprocess, "run", lambda command, **kwargs: commands.append(command))
    RustProfileBackend().run(_request(tmp_path, max_concurrency=8))
    command = commands[0]
    assert command[command.index("--threads") + 1] == "2"


def test_adapter_reports_missing_binary(tmp_path, monkeypatch):
    monkeypatch.setattr(zipstrain_rust_profiler.sysconfig, "get_path", lambda key: str(tmp_path))
    monkeypatch.setattr(zipstrain_rust_profiler.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="binary is missing"):
        RustProfileBackend().run(_request(tmp_path))
