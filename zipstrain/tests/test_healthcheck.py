from click.testing import CliRunner

from zipstrain import cli
from zipstrain import healthcheck as hc


def _check_map(report: hc.HealthReport) -> dict[tuple[str, str], hc.HealthCheck]:
    return {(check.category, check.name): check for check in report.checks}


def test_collect_health_report_all_ready():
    def executable_probe(command: str) -> tuple[bool, str]:
        return True, f"/usr/local/bin/{command}"

    def package_probe(import_name: str, package_name: str | None) -> tuple[bool, str]:
        return True, f"installed ({import_name}-1.0)"

    report = hc.collect_health_report(
        executable_probe=executable_probe,
        package_probe=package_probe,
    )
    checks = _check_map(report)

    assert checks[("External tools", "samtools")].status == "ok"
    assert checks[("External tools", "prodigal")].status == "ok"
    assert checks[("External tools", "sylph")].status == "ok"
    assert checks[("Python packages", "torch")].status == "ok"
    assert checks[("Python packages", "h5py")].status == "ok"
    assert checks[("Workflow readiness", "BAM profiling")].status == "ok"
    assert checks[("Workflow readiness", "Standard comparison")].status == "ok"
    assert checks[("Workflow readiness", "Matrix comparison")].status == "ok"


def test_collect_health_report_marks_optional_and_workflow_gaps():
    def executable_probe(command: str) -> tuple[bool, str]:
        if command == "samtools":
            return False, "not found in PATH"
        return False, "not found in PATH"

    def package_probe(import_name: str, package_name: str | None) -> tuple[bool, str]:
        return False, "not installed"

    report = hc.collect_health_report(
        executable_probe=executable_probe,
        package_probe=package_probe,
    )
    checks = _check_map(report)

    assert checks[("External tools", "samtools")].status == "missing"
    assert checks[("External tools", "prodigal")].status == "missing"
    assert checks[("External tools", "sylph")].status == "missing"
    assert checks[("Python packages", "torch")].status == "missing"
    assert checks[("Python packages", "h5py")].status == "missing"
    assert checks[("Workflow readiness", "BAM profiling")].status == "warn"
    assert checks[("Workflow readiness", "Standard comparison")].status == "ok"
    assert checks[("Workflow readiness", "Matrix comparison")].status == "warn"


def test_cli_test_renders_friendly_health_report(monkeypatch):
    fake_report = hc.HealthReport(
        checks=(
            hc.HealthCheck(
                category="External tools",
                name="samtools",
                status="ok",
                detail="/usr/local/bin/samtools",
                note="Required for BAM -> profile workflows.",
            ),
            hc.HealthCheck(
                category="Workflow readiness",
                name="Matrix comparison",
                status="warn",
                detail="not ready",
                note='Needs both torch and h5py. Install with: pip install "zipstrain[matrix]"',
            ),
        )
    )

    monkeypatch.setattr(hc, "collect_health_report", lambda: fake_report)

    runner = CliRunner()
    result = runner.invoke(cli.cli, ["test"])

    assert result.exit_code == 0
    assert "ZipStrain health check" in result.output
    assert "samtools" in result.output
    assert "Matrix comparison" in result.output
    assert "MISSING" not in result.output
