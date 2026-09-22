import inspect
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from zipstrain import cli, profile, utils


@pytest.mark.parametrize("command", [cli.profile, cli.prepare_profiling, cli.build_null_model])
def test_cli_null_model_defaults(command):
    defaults = {param.name: param.default for param in command.params}
    assert defaults["error_rate"] == 0.001
    assert defaults["p_threshold"] == 1e-6
    if "min_freq" in defaults:
        assert defaults["min_freq"] == 0.01


def test_python_and_nextflow_defaults_agree():
    for function in (utils.build_null_poisson, profile.prepare_profiling_assets):
        params = inspect.signature(function).parameters
        assert params["error_rate"].default == 0.001
        assert params["p_threshold"].default == 1e-6
    assert profile.PROFILE_MIN_FREQ_DEFAULT == 0.01
    text = (Path(__file__).resolve().parents[2] / "zipstrain.nf").read_text()
    assert "params.error_rate=0.001\n" in text
    assert "params.p_threshold=0.000001\n" in text
    assert "params.min_freq=0.01\n" in text


def test_default_model_removes_spurious_shared_allele_and_allows_override(tmp_path):
    raw = pl.DataFrame({"A": [0], "C": [78], "G": [0], "T": [1]}).lazy()
    for threshold, expected_t in ((None, 0), ("0.05", 1)):
        output = tmp_path / f"null_{threshold}.parquet"
        args = ["--max-total-reads", "100", "--output-file", str(output)]
        if threshold is not None:
            args += ["--p-threshold", threshold]
        result = CliRunner().invoke(cli.build_null_model, args)
        assert result.exit_code == 0, result.output
        adjusted = profile.adjust_for_sequence_errors(raw, pl.scan_parquet(output)).collect()
        assert adjusted["C"].item() == 78
        assert adjusted["T"].item() == expected_t
