"""Adapter for the optional ZipStrain Rust profiling executable."""

from pathlib import Path
import shutil
import subprocess
import sysconfig

from zipstrain.resource_limits import cpu_budget


class RustProfileBackend:
    API_VERSION = 1

    def run(self, request) -> None:
        executable_name = "zipstrain-rust-profiler"
        if sysconfig.get_platform().startswith("win"):
            executable_name += ".exe"
        binary = Path(sysconfig.get_path("scripts")) / executable_name
        if not binary.is_file():
            found = shutil.which(executable_name)
            if found is None:
                raise RuntimeError(
                    f"The {executable_name} binary is missing. Reinstall the "
                    "zipstrain_rust_profiler wheel in this Python environment."
                )
            binary = Path(found)

        budget = cpu_budget()
        threads = min(request.max_concurrency, budget) if budget is not None else request.max_concurrency
        command = [
            str(binary),
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
            ("--gene-range", request.gene_range_table),
            ("--profiling-contract", request.profiling_contract),
            ("--min-read-ani", request.min_read_ani),
        ):
            if value is not None:
                command.extend((option, str(value)))
        subprocess.run(command, check=True)
