"""User-facing ZipStrain environment and workflow health checks."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
from importlib import metadata
import shutil
from typing import Callable

from rich.console import Console
from rich.panel import Panel
from rich.table import Table


@dataclass(frozen=True)
class HealthCheck:
    """A single dependency or workflow readiness check."""

    category: str
    name: str
    status: str
    detail: str
    note: str | None = None


@dataclass(frozen=True)
class HealthReport:
    """Structured health-check report for CLI rendering and tests."""

    checks: tuple[HealthCheck, ...]


ExecutableProbe = Callable[[str], tuple[bool, str]]
PackageProbe = Callable[[str, str | None], tuple[bool, str]]


def _probe_executable(command: str) -> tuple[bool, str]:
    """Return whether an executable is available in PATH and its resolved path."""
    resolved = shutil.which(command)
    if resolved is None:
        return False, "not found in PATH"
    return True, resolved


def _probe_python_package(import_name: str, package_name: str | None = None) -> tuple[bool, str]:
    """Return whether a Python package can be imported and its installed version."""
    try:
        module = importlib.import_module(import_name)
    except ImportError:
        return False, "not installed"

    try:
        version = metadata.version(package_name or import_name)
    except metadata.PackageNotFoundError:
        version = getattr(module, "__version__", None)

    if version:
        return True, f"installed ({version})"
    return True, "installed"


def collect_health_report(
    *,
    executable_probe: ExecutableProbe = _probe_executable,
    package_probe: PackageProbe = _probe_python_package,
) -> HealthReport:
    """Collect a lightweight environment and workflow readiness report."""
    checks: list[HealthCheck] = []

    samtools_ok, samtools_detail = executable_probe("samtools")
    prodigal_ok, prodigal_detail = executable_probe("prodigal")
    sylph_ok, sylph_detail = executable_probe("sylph")
    torch_ok, torch_detail = package_probe("torch", "torch")
    h5py_ok, h5py_detail = package_probe("h5py", "h5py")

    checks.extend(
        [
            HealthCheck(
                category="External tools",
                name="samtools",
                status="ok" if samtools_ok else "missing",
                detail=samtools_detail,
                note="Required for BAM -> profile workflows.",
            ),
            HealthCheck(
                category="External tools",
                name="prodigal",
                status="ok" if prodigal_ok else "missing",
                detail=prodigal_detail,
                note="Optional. Helpful when you need to generate gene annotations from reference FASTA files.",
            ),
            HealthCheck(
                category="External tools",
                name="sylph",
                status="ok" if sylph_ok else "missing",
                detail=sylph_detail,
                note="Optional. Needed only for the Sylph-based reference-build route.",
            ),
            HealthCheck(
                category="Python packages",
                name="torch",
                status="ok" if torch_ok else "missing",
                detail=torch_detail,
                note='Optional. Needed for matrix comparison. Install with: pip install "zipstrain[matrix]"',
            ),
            HealthCheck(
                category="Python packages",
                name="h5py",
                status="ok" if h5py_ok else "missing",
                detail=h5py_detail,
                note='Optional. Needed for matrix-store IO. Install with: pip install "zipstrain[matrix]"',
            ),
        ]
    )

    checks.extend(
        [
            HealthCheck(
                category="Workflow readiness",
                name="BAM profiling",
                status="ok" if samtools_ok else "warn",
                detail="ready" if samtools_ok else "not ready",
                note=(
                    "You can profile BAM inputs."
                    if samtools_ok
                    else "Needs samtools in PATH before `zipstrain profile ...` can process BAM files."
                ),
            ),
            HealthCheck(
                category="Workflow readiness",
                name="Standard comparison",
                status="ok",
                detail="ready",
                note="Core profile-table comparison dependencies are installed with ZipStrain.",
            ),
            HealthCheck(
                category="Workflow readiness",
                name="Matrix comparison",
                status="ok" if torch_ok and h5py_ok else "warn",
                detail="ready" if torch_ok and h5py_ok else "not ready",
                note=(
                    "Torch and HDF5 support are available."
                    if torch_ok and h5py_ok
                    else 'Needs both torch and h5py. Install with: pip install "zipstrain[matrix]"'
                ),
            ),
        ]
    )

    return HealthReport(checks=tuple(checks))


def render_health_report(report: HealthReport, *, console: Console | None = None) -> None:
    """Render the health report in a compact, user-friendly CLI format."""
    console = console or Console()

    checks_by_category: dict[str, list[HealthCheck]] = {}
    for check in report.checks:
        checks_by_category.setdefault(check.category, []).append(check)

    ok_count = sum(1 for check in report.checks if check.status == "ok")
    warn_count = sum(1 for check in report.checks if check.status == "warn")
    missing_count = sum(1 for check in report.checks if check.status == "missing")

    console.print(
        Panel.fit(
            (
                "[bold]ZipStrain health check[/bold]\n"
                "This report highlights what is ready now and what is optional but missing."
            ),
            border_style="blue",
        )
    )

    for category, checks in checks_by_category.items():
        table = Table(title=category, show_header=True, header_style="bold")
        table.add_column("Check", style="cyan")
        table.add_column("Status", style="bold")
        table.add_column("Details")
        table.add_column("Notes")
        for check in checks:
            if check.status == "ok":
                status_label = "[green]OK[/green]"
            elif check.status == "warn":
                status_label = "[yellow]WARN[/yellow]"
            else:
                status_label = "[red]MISSING[/red]"
            table.add_row(check.name, status_label, check.detail, check.note or "")
        console.print(table)

    summary_bits = [f"{ok_count} ok"]
    if warn_count:
        summary_bits.append(f"{warn_count} warnings")
    if missing_count:
        summary_bits.append(f"{missing_count} missing")
    console.print(
        Panel.fit(
            "[bold]Summary:[/bold] " + ", ".join(summary_bits),
            border_style="green" if missing_count == 0 else "yellow",
        )
    )
