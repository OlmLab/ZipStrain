"""Adapter for the optional ZipStrain Rust profiling executable."""

from pathlib import Path
import shutil
import subprocess
import sysconfig
import tempfile

from zipstrain.resource_limits import cpu_budget


def _binary() -> Path:
    executable_name = "zipstrain-rust-profiler"
    if sysconfig.get_platform().startswith("win"):
        executable_name += ".exe"
    binary = Path(sysconfig.get_path("scripts")) / executable_name
    if binary.is_file():
        return binary
    found = shutil.which(executable_name)
    if found is None:
        raise RuntimeError(
            f"The {executable_name} binary is missing. Reinstall the "
            "zipstrain_rust_profiler wheel in this Python environment."
        )
    return Path(found)


class RustProfileBackend:
    API_VERSION = 1

    def extract_codons(self, bam_file, genes, output_file, *, min_baseq, min_mapq, min_read_ani,
                       read_inclusion, threads) -> None:
        """Write raw per-fragment codon counts. ZipStrain gates and finalizes them."""
        budget = cpu_budget()
        threads = min(threads, budget) if budget is not None else threads
        with tempfile.TemporaryDirectory(prefix="zipstrain_rust_codons_") as staging:
            table = Path(staging) / "genes.tsv"
            genes.select("gene_id", "chrom", "start", "end", "length", "strand", "phase").write_csv(
                table, separator="\t", include_header=False)
            command = [
                str(_binary()), "codons",
                "--bam", str(bam_file), "--genes", str(table), "--output", str(output_file),
                "--threads", str(threads), "--min-baseq", str(min_baseq), "--min-mapq", str(min_mapq),
                "--read-inclusion", str(read_inclusion),
            ]
            if min_read_ani is not None:
                command.extend(("--min-read-ani", str(min_read_ani)))
            subprocess.run(command, check=True)

    def run(self, request) -> None:
        budget = cpu_budget()
        threads = min(request.max_concurrency, budget) if budget is not None else request.max_concurrency
        command = [
            str(_binary()),
            "--bam", str(request.bam_file),
            "--bed", str(request.bed_file),
            "--stb", str(request.stb_file),
            "--null-model", str(request.null_model),
            "--output-dir", str(request.output_dir),
            "--threads", str(threads),
            "--num-chunks", str(request.num_chunks),
            "--min-mapq", str(request.min_mapq),
            "--min-baseq", str(request.min_baseq),
            "--min-freq", str(request.min_freq),
            "--read-inclusion", str(request.read_inclusion),
        ]
        for option, value in (
            ("--reference", request.reference_fasta),
            ("--profiling-contract", request.profiling_contract),
            ("--min-read-ani", request.min_read_ani),
        ):
            if value is not None:
                command.extend((option, str(value)))
        subprocess.run(command, check=True)
