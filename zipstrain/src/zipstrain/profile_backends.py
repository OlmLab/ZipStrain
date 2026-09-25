"""Profiling backend boundary shared by the CLI and optional distributions."""

from dataclasses import dataclass
from importlib.metadata import entry_points
from pathlib import Path


PROFILE_BACKEND_API_VERSION = 1
PROFILE_BACKEND_GROUP = "zipstrain.profile_backends"


class ProfileBackendError(ValueError):
    """Invalid or unavailable optional profiling backend."""


@dataclass(frozen=True)
class ProfileRequest:
    bam_file: Path
    bed_file: Path
    stb_file: Path
    null_model: Path
    output_dir: Path
    reference_fasta: Path | None = None
    gene_range_table: Path | None = None
    profiling_contract: Path | None = None
    num_chunks: int = 24
    max_concurrency: int = 4
    min_mapq: int = 0
    min_baseq: int = 13
    min_freq: float = 0.01
    min_read_ani: float | None = None
    read_inclusion: str = "all-mapped"


def load_backend(name: str):
    """Load only the selected backend; an absent plugin never affects Python."""
    if name == "python":
        return None
    if not name:
        raise ProfileBackendError(f"Invalid profiling backend name: {name!r}")
    matches = list(entry_points(group=PROFILE_BACKEND_GROUP, name=name))
    if not matches:
        raise ProfileBackendError(
            f"Profiling backend {name!r} is not installed. Install its plugin "
            f"in the environment running profile-single (entry-point group {PROFILE_BACKEND_GROUP})."
        )
    if len(matches) != 1:
        raise ProfileBackendError(f"Multiple profiling backends are registered as {name!r}")
    backend = matches[0].load()()
    if getattr(backend, "API_VERSION", None) != PROFILE_BACKEND_API_VERSION:
        raise ProfileBackendError(
            f"Profiling backend {name!r} has an incompatible API version; "
            f"ZipStrain requires {PROFILE_BACKEND_API_VERSION}."
        )
    if not callable(getattr(backend, "run", None)):
        raise ProfileBackendError(f"Profiling backend {name!r} must provide run(request)")
    return backend


def run_profile(request: ProfileRequest, backend_name: str = "python") -> None:
    backend = load_backend(backend_name)
    request.output_dir.mkdir(parents=True, exist_ok=True)
    if backend is None:
        import polars as pl

        from zipstrain import profile, utils

        contract = (
            utils.read_profile_contract_file(request.profiling_contract)
            if request.profiling_contract is not None
            else None
        )
        profile.profile_bam(
            bed_file=request.bed_file,
            bam_file=request.bam_file,
            reference_fasta=request.reference_fasta,
            gene_range_table=request.gene_range_table,
            stb=profile.read_stb(request.stb_file),
            null_model=pl.scan_parquet(request.null_model),
            output_dir=request.output_dir,
            num_chunks=request.num_chunks,
            max_concurrency=request.max_concurrency,
            profile_contract=contract,
            min_mapq=request.min_mapq,
            min_baseq=request.min_baseq,
            min_freq=request.min_freq,
            min_read_ani=request.min_read_ani,
            read_inclusion=request.read_inclusion,
        )
        return

    backend.run(request)
    stem = request.bam_file.stem
    missing = [
        str(request.output_dir / f"{stem}_{suffix}.parquet")
        for suffix in ("profile", "gene_stats", "genome_stats")
        if not (request.output_dir / f"{stem}_{suffix}.parquet").is_file()
    ]
    if missing:
        raise ProfileBackendError(f"Profiling backend {backend_name!r} did not write: {', '.join(missing)}")
