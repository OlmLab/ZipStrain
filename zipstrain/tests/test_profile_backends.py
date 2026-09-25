from pathlib import Path
from importlib.metadata import EntryPoints

from click.testing import CliRunner
import pytest

from zipstrain import cli, profile_backends


class _EntryPoint:
    def __init__(self, backend):
        self.backend = backend

    def load(self):
        return lambda: self.backend


class _Backend:
    API_VERSION = 1

    def __init__(self):
        self.requests = []

    def run(self, request):
        self.requests.append(request)
        for suffix in ("profile", "gene_stats", "genome_stats"):
            (request.output_dir / f"{request.bam_file.stem}_{suffix}.parquet").touch()


def _request(tmp_path):
    return profile_backends.ProfileRequest(
        bam_file=tmp_path / "sample.bam",
        bed_file=tmp_path / "regions.bed",
        stb_file=tmp_path / "mapping.stb",
        null_model=tmp_path / "null.parquet",
        output_dir=tmp_path / "out",
    )


def test_python_backend_does_not_discover_plugins(monkeypatch):
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: pytest.fail("discovered plugins"))
    assert profile_backends.load_backend("python") is None


def test_missing_backend_has_actionable_error(monkeypatch):
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: EntryPoints(()))
    with pytest.raises(ValueError, match="not installed"):
        profile_backends.load_backend("rust_profiler")


def test_plugin_contract_rejects_duplicates_and_api_mismatch(monkeypatch):
    backend = _Backend()
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: [_EntryPoint(backend)] * 2)
    with pytest.raises(ValueError, match="Multiple"):
        profile_backends.load_backend("rust_profiler")
    backend.API_VERSION = 99
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: [_EntryPoint(backend)])
    with pytest.raises(ValueError, match="incompatible API"):
        profile_backends.load_backend("rust_profiler")


def test_plugin_receives_request_and_required_outputs_are_checked(tmp_path, monkeypatch):
    backend = _Backend()
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: [_EntryPoint(backend)])
    request = _request(tmp_path)
    profile_backends.run_profile(request, "rust_profiler")
    assert backend.requests == [request]

    request = profile_backends.ProfileRequest(**{
        **request.__dict__, "output_dir": tmp_path / "missing",
    })
    backend.run = lambda request: None
    with pytest.raises(profile_backends.ProfileBackendError, match="did not write"):
        profile_backends.run_profile(request, "rust_profiler")


def test_profile_single_reports_missing_plugin_without_traceback(tmp_path, monkeypatch):
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: EntryPoints(()))
    result = CliRunner().invoke(cli.cli, [
        "utilities", "profile-single", "--bam-file", str(tmp_path / "sample.bam"),
        "--bed-file", str(tmp_path / "regions.bed"), "--stb-file", str(tmp_path / "mapping.stb"),
        "--null-model", str(tmp_path / "null.parquet"), "--output-dir", str(tmp_path / "out"),
        "--backend", "rust_profiler",
    ])
    assert result.exit_code != 0
    assert "is not installed" in result.output


def test_profile_single_passes_backend_to_plugin(tmp_path, monkeypatch):
    backend = _Backend()
    monkeypatch.setattr(profile_backends, "entry_points", lambda **kwargs: [_EntryPoint(backend)])
    result = CliRunner().invoke(cli.cli, [
        "utilities", "profile-single", "--bam-file", str(tmp_path / "sample.bam"),
        "--bed-file", str(tmp_path / "regions.bed"), "--stb-file", str(tmp_path / "mapping.stb"),
        "--null-model", str(tmp_path / "null.parquet"), "--output-dir", str(tmp_path / "out"),
        "--backend", "rust_profiler", "--max-concurrency", "2",
    ])
    assert result.exit_code == 0, result.output
    assert backend.requests[0].max_concurrency == 2
    assert backend.requests[0].bam_file == Path(tmp_path / "sample.bam")
