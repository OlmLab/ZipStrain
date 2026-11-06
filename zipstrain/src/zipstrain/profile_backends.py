"""Profiling backend boundary shared by the CLI and optional distributions."""

from dataclasses import dataclass, replace
from importlib.metadata import entry_points
from pathlib import Path

from zipstrain import options


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
    min_read_ani: float | None = options.PROFILE_MIN_READ_ANI
    read_inclusion: str = options.PROFILE_READ_INCLUSION
    collect_codons: bool = False
    allele_integration: str = "consensus"


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
    """Run a selected profiler and normalize all outputs to the V2 file contract.

    API-v1 native backends produce position counts without gene annotation;
    ``finish_profile_outputs`` applies the authenticated catalog and optional
    observed codon preparation consistently. The Python backend finalizes inside
    its profiling function. Read and allele-integration settings are preserved.
    """
    backend = load_backend(backend_name)
    runtime_reference_id = None
    if backend is not None and request.reference_fasta is not None:
        from zipstrain import gene_info, profile
        runtime_reference_id = gene_info.reference_identity(request.reference_fasta, profile.read_stb(request.stb_file))
    if backend is not None and request.gene_range_table is not None:
        from zipstrain import gene_info, profile
        gene_info.read_gene_info(request.gene_range_table)
        if request.reference_fasta is not None:
            import polars as pl
            if runtime_reference_id != pl.read_parquet_metadata(request.gene_range_table)[gene_info.REFERENCE_ID]:
                raise ProfileBackendError("Profiling reference/STB do not match gene_info.parquet")
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
            collect_codons=request.collect_codons,
            allele_integration=request.allele_integration,
        )
        return

    # V2 resolves all gene annotations in the shared finalizer. Passing no gene
    # ranges also avoids redundant annotation work in older API-v1 plugins.
    backend.run(replace(request, gene_range_table=None))
    stem = request.bam_file.stem
    missing = [
        str(request.output_dir / f"{stem}_{suffix}.parquet")
        for suffix in ("profile", "genome_stats")
        if not (request.output_dir / f"{stem}_{suffix}.parquet").is_file()
    ]
    if missing:
        raise ProfileBackendError(f"Profiling backend {backend_name!r} did not write: {', '.join(missing)}")
    from zipstrain.profile import finish_profile_outputs
    finish_profile_outputs(request.output_dir / f"{stem}_profile.parquet", request.gene_range_table,
        bam_file=request.bam_file, collect_codons=request.collect_codons,
        min_mapq=request.min_mapq, min_baseq=request.min_baseq,
        min_read_ani=request.min_read_ani, read_inclusion=request.read_inclusion,
        min_freq=request.min_freq, allele_integration=request.allele_integration,
        reference_id=runtime_reference_id)
