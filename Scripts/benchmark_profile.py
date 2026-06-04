#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import random
import re
import shlex
import shutil
import stat
import subprocess
import sys
import time

import polars as pl


PRESETS: dict[str, dict[str, int]] = {
    "tiny": {
        "genome_count": 2,
        "contigs_per_genome": 1,
        "contig_length": 2_000,
        "read_length": 150,
        "tile_step": 75,
        "read_copies": 1,
        "gene_length": 300,
        "gene_step": 450,
    },
    "small": {
        "genome_count": 2,
        "contigs_per_genome": 2,
        "contig_length": 20_000,
        "read_length": 150,
        "tile_step": 40,
        "read_copies": 2,
        "gene_length": 900,
        "gene_step": 1_200,
    },
    "medium": {
        "genome_count": 4,
        "contigs_per_genome": 3,
        "contig_length": 100_000,
        "read_length": 150,
        "tile_step": 30,
        "read_copies": 2,
        "gene_length": 900,
        "gene_step": 1_200,
    },
    "large": {
        "genome_count": 6,
        "contigs_per_genome": 4,
        "contig_length": 150_000,
        "read_length": 150,
        "tile_step": 25,
        "read_copies": 2,
        "gene_length": 900,
        "gene_step": 1_200,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark `zipstrain utilities profile-single` on a reusable synthetic BAM fixture. "
            "Create the dataset once, then rerun the exact same inputs across commits."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("work/profile_benchmark"),
        help="Directory used for the reusable benchmark fixture and benchmark outputs.",
    )
    parser.add_argument(
        "--zipstrain-command",
        default=f"{shlex.quote(sys.executable)} -m zipstrain.cli",
        help=(
            "Command prefix used to invoke ZipStrain, for example "
            "`zipstrain` or `/path/to/python -m zipstrain.cli`."
        ),
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Label for the benchmark result file. Defaults to the current git short SHA if available.",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=3,
        help="Number of timed profile-single runs to record.",
    )
    parser.add_argument(
        "--num-chunks",
        type=int,
        default=24,
        help="Number of BED chunks to pass to profile-single.",
    )
    parser.add_argument(
        "--max-concurrency",
        type=int,
        default=4,
        help="Chunk concurrency to pass to profile-single.",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default="medium",
        help="Synthetic dataset size preset. `medium` is intended to be large enough to show real differences.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=13,
        help="Random seed used to generate the synthetic reference sequence.",
    )
    parser.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Regenerate the fixture even if it already exists.",
    )
    parser.add_argument(
        "--skip-setup",
        action="store_true",
        help="Reuse an existing fixture and do not regenerate inputs.",
    )
    parser.add_argument(
        "--setup-only",
        action="store_true",
        help="Generate the fixture and exit without running the benchmark.",
    )
    parser.add_argument(
        "--keep-run-dirs",
        action="store_true",
        help="Keep per-run output directories instead of deleting them after each timed run.",
    )
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("--repetitions must be >= 1")
    if args.num_chunks < 1:
        parser.error("--num-chunks must be >= 1")
    if args.max_concurrency < 1:
        parser.error("--max-concurrency must be >= 1")
    return args


def _git_short_sha(cwd: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def _ensure_zipstrain_wrapper(bin_dir: Path, zipstrain_command: str) -> Path:
    bin_dir.mkdir(parents=True, exist_ok=True)
    wrapper_path = bin_dir / "zipstrain"
    wrapper_path.write_text(
        "#!/bin/sh\n"
        f"exec {zipstrain_command} \"$@\"\n"
    )
    wrapper_path.chmod(wrapper_path.stat().st_mode | stat.S_IEXEC)
    return wrapper_path


def _choose_time_wrapper() -> tuple[list[str], str | None]:
    which_time = shutil.which("time")
    candidates = [Path("/usr/bin/time")]
    if which_time:
        candidates.append(Path(which_time))
    for candidate in candidates:
        if not str(candidate):
            continue
        if candidate.exists():
            result = subprocess.run(
                [str(candidate), "--version"],
                capture_output=True,
                text=True,
            )
            version_text = f"{result.stdout}\n{result.stderr}".lower()
            if "gnu time" in version_text:
                return [str(candidate), "-v"], "gnu"
            return [str(candidate), "-l"], "bsd"
    return [], None


def _parse_max_rss_kb(stderr: str, flavor: str | None) -> int | None:
    if flavor == "bsd":
        match = re.search(r"(?m)^(\d+)\s+maximum resident set size$", stderr)
        return int(match.group(1)) if match else None
    if flavor == "gnu":
        match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", stderr)
        return int(match.group(1)) if match else None
    return None


def _random_sequence(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _iter_contig_specs(genome_count: int, contigs_per_genome: int) -> list[tuple[str, str]]:
    specs: list[tuple[str, str]] = []
    for genome_idx in range(genome_count):
        genome_name = f"genome_{genome_idx + 1:02d}"
        for contig_idx in range(contigs_per_genome):
            contig_name = f"{genome_name}_contig_{contig_idx + 1:02d}"
            specs.append((genome_name, contig_name))
    return specs


def _write_reference_files(
    fixture_dir: Path,
    *,
    genome_count: int,
    contigs_per_genome: int,
    contig_length: int,
    gene_length: int,
    gene_step: int,
    seed: int,
) -> dict[str, object]:
    rng = random.Random(seed)
    reference_fasta = fixture_dir / "reference.fna"
    bed_file = fixture_dir / "genomes.bed"
    gene_range_table = fixture_dir / "gene_range_table.tsv"
    stb_file = fixture_dir / "reference.stb"
    metadata_path = fixture_dir / "fixture_metadata.json"

    contig_specs = _iter_contig_specs(genome_count, contigs_per_genome)
    reference_lines: list[str] = []
    bed_lines: list[str] = []
    gene_lines: list[str] = []
    stb_lines: list[str] = []
    sequences: dict[str, str] = {}
    gene_counter = 0

    for genome_name, contig_name in contig_specs:
        sequence = _random_sequence(contig_length, rng)
        sequences[contig_name] = sequence
        reference_lines.extend([f">{contig_name}", sequence])
        bed_lines.append(f"{contig_name}\t0\t{contig_length}")
        stb_lines.append(f"{contig_name}\t{genome_name}")
        gene_start = 1
        while gene_start <= contig_length:
            gene_end = min(gene_start + gene_length - 1, contig_length)
            gene_counter += 1
            gene_name = f"{contig_name}_gene_{gene_counter:05d}"
            gene_lines.append(f"{gene_name}\t{contig_name}\t{gene_start}\t{gene_end}")
            gene_start += gene_step

    reference_fasta.write_text("\n".join(reference_lines) + "\n")
    bed_file.write_text("\n".join(bed_lines) + "\n")
    gene_range_table.write_text("\n".join(gene_lines) + "\n")
    stb_file.write_text("\n".join(stb_lines) + "\n")

    metadata = {
        "genome_count": genome_count,
        "contigs_per_genome": contigs_per_genome,
        "contig_length": contig_length,
        "gene_length": gene_length,
        "gene_step": gene_step,
        "seed": seed,
        "contig_count": len(contig_specs),
        "position_count": len(contig_specs) * contig_length,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return {
        "reference_fasta": reference_fasta,
        "bed_file": bed_file,
        "gene_range_table": gene_range_table,
        "stb_file": stb_file,
        "sequences": sequences,
        "metadata": metadata,
    }


def _write_null_model(null_model_path: Path, max_cov: int = 512) -> None:
    pl.DataFrame(
        {
            "cov": list(range(max_cov + 1)),
            "max_error_count": [0] * (max_cov + 1),
        }
    ).write_parquet(null_model_path)


def _stream_sam_to_sorted_bam(
    *,
    bam_path: Path,
    contig_sequences: dict[str, str],
    read_length: int,
    tile_step: int,
    read_copies: int,
) -> dict[str, int]:
    samtools = shutil.which("samtools")
    if samtools is None:
        raise RuntimeError("samtools is required to generate the benchmark BAM.")
    unsorted_bam = bam_path.with_suffix(".unsorted.bam")
    header_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for contig_name, sequence in contig_sequences.items():
        header_lines.append(f"@SQ\tSN:{contig_name}\tLN:{len(sequence)}")

    total_reads = 0
    total_bases = 0
    read_starts_per_contig = 0
    quality = "I" * read_length
    with subprocess.Popen(
        [samtools, "view", "-bS", "-o", str(unsorted_bam), "-"],
        stdin=subprocess.PIPE,
        text=True,
    ) as proc:
        assert proc.stdin is not None
        proc.stdin.write("\n".join(header_lines) + "\n")
        for contig_name, sequence in contig_sequences.items():
            contig_read_starts = 0
            last_start = len(sequence) - read_length
            if last_start < 0:
                raise ValueError(f"read_length={read_length} exceeds contig length for {contig_name}")
            for start0 in range(0, last_start + 1, tile_step):
                contig_read_starts += 1
                read_seq = sequence[start0:start0 + read_length]
                for copy_idx in range(read_copies):
                    total_reads += 1
                    total_bases += read_length
                    qname = f"benchmark_{contig_name}_{start0 + 1:07d}_{copy_idx + 1:02d}"
                    proc.stdin.write(
                        f"{qname}\t0\t{contig_name}\t{start0 + 1}\t60\t{read_length}M\t*\t0\t0\t{read_seq}\t{quality}\n"
                    )
            read_starts_per_contig = contig_read_starts
        proc.stdin.close()
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, proc.args)

    subprocess.run([samtools, "sort", "-o", str(bam_path), str(unsorted_bam)], check=True)
    subprocess.run([samtools, "index", str(bam_path)], check=True)
    unsorted_bam.unlink(missing_ok=True)
    return {
        "total_reads": total_reads,
        "total_read_bases": total_bases,
        "read_starts_per_contig": read_starts_per_contig,
    }


def ensure_fixture(
    fixture_dir: Path,
    *,
    preset_name: str,
    seed: int,
    force_rebuild: bool,
) -> dict[str, object]:
    if force_rebuild and fixture_dir.exists():
        shutil.rmtree(fixture_dir)
    fixture_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = fixture_dir / "fixture_metadata.json"
    bam_path = fixture_dir / "benchmark_sample.bam"
    null_model_path = fixture_dir / "null_model.parquet"
    preset = PRESETS[preset_name]

    if metadata_path.exists() and bam_path.exists() and null_model_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("preset") != preset_name or metadata.get("seed") != seed:
            raise ValueError(
                "Existing fixture metadata does not match the requested preset/seed. "
                "Use --force-rebuild or change --work-dir."
            )
        return {
            "fixture_dir": fixture_dir,
            "metadata": metadata,
            "bam_file": bam_path,
            "bed_file": fixture_dir / "genomes.bed",
            "stb_file": fixture_dir / "reference.stb",
            "gene_range_table": fixture_dir / "gene_range_table.tsv",
            "null_model": null_model_path,
        }

    reference_data = _write_reference_files(
        fixture_dir,
        genome_count=preset["genome_count"],
        contigs_per_genome=preset["contigs_per_genome"],
        contig_length=preset["contig_length"],
        gene_length=preset["gene_length"],
        gene_step=preset["gene_step"],
        seed=seed,
    )
    _write_null_model(null_model_path)
    bam_stats = _stream_sam_to_sorted_bam(
        bam_path=bam_path,
        contig_sequences=reference_data["sequences"],
        read_length=preset["read_length"],
        tile_step=preset["tile_step"],
        read_copies=preset["read_copies"],
    )
    metadata = {
        **reference_data["metadata"],
        "preset": preset_name,
        "read_length": preset["read_length"],
        "tile_step": preset["tile_step"],
        "read_copies": preset["read_copies"],
        **bam_stats,
        "approx_mean_depth": round((preset["read_length"] / preset["tile_step"]) * preset["read_copies"], 2),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return {
        "fixture_dir": fixture_dir,
        "metadata": metadata,
        "bam_file": bam_path,
        "bed_file": reference_data["bed_file"],
        "stb_file": reference_data["stb_file"],
        "gene_range_table": reference_data["gene_range_table"],
        "null_model": null_model_path,
    }


def _build_profile_command(
    zipstrain_command: str,
    *,
    fixture: dict[str, object],
    output_dir: Path,
    num_chunks: int,
    max_concurrency: int,
) -> list[str]:
    return [
        *shlex.split(zipstrain_command),
        "utilities",
        "profile-single",
        "--bed-file",
        str(fixture["bed_file"]),
        "--bam-file",
        str(fixture["bam_file"]),
        "--stb-file",
        str(fixture["stb_file"]),
        "--null-model",
        str(fixture["null_model"]),
        "--gene-range-table",
        str(fixture["gene_range_table"]),
        "--num-chunks",
        str(num_chunks),
        "--max-concurrency",
        str(max_concurrency),
        "--output-dir",
        str(output_dir),
    ]


def _run_timed_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> dict[str, object]:
    time_prefix, flavor = _choose_time_wrapper()
    attempts: list[tuple[list[str], str | None]] = []
    if time_prefix:
        attempts.append(([*time_prefix, *cmd], flavor))
    attempts.append((cmd, None))
    last_error: RuntimeError | None = None

    for full_cmd, attempt_flavor in attempts:
        wall_start = time.perf_counter()
        proc = subprocess.run(full_cmd, cwd=cwd, env=env, capture_output=True, text=True)
        wall_seconds = time.perf_counter() - wall_start
        stderr = proc.stderr or ""
        max_rss_kb = _parse_max_rss_kb(stderr, attempt_flavor)
        if proc.returncode == 0:
            return {
                "wall_seconds": wall_seconds,
                "max_rss_kb": max_rss_kb,
                "stdout": proc.stdout,
                "stderr": stderr,
                "timing_flavor": attempt_flavor,
            }
        last_error = RuntimeError(
            f"Benchmark command failed with exit code {proc.returncode}\n"
            f"CMD: {shlex.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{stderr}"
        )

    assert last_error is not None
    raise last_error


def _summarize_profile_outputs(output_dir: Path) -> dict[str, object]:
    profile_path = next(output_dir.glob("*_profile.parquet"))
    gene_stats_path = next(output_dir.glob("*_gene_stats.parquet"))
    genome_stats_path = next(output_dir.glob("*_genome_stats.parquet"))
    profile_df = pl.scan_parquet(profile_path).select(pl.len().alias("rows")).collect(engine="streaming")
    gene_df = pl.scan_parquet(gene_stats_path).select(pl.len().alias("rows")).collect(engine="streaming")
    genome_df = pl.scan_parquet(genome_stats_path).select(pl.len().alias("rows")).collect(engine="streaming")
    return {
        "profile_file": str(profile_path),
        "gene_stats_file": str(gene_stats_path),
        "genome_stats_file": str(genome_stats_path),
        "profile_rows": int(profile_df.item()),
        "gene_stats_rows": int(gene_df.item()),
        "genome_stats_rows": int(genome_df.item()),
        "profile_size_bytes": profile_path.stat().st_size,
        "gene_stats_size_bytes": gene_stats_path.stat().st_size,
        "genome_stats_size_bytes": genome_stats_path.stat().st_size,
    }


def benchmark_profile(args: argparse.Namespace) -> dict[str, object]:
    repo_root = Path(__file__).resolve().parents[1]
    work_dir = args.work_dir.resolve()
    fixture_dir = work_dir / "fixture"
    runs_dir = work_dir / "runs"
    results_dir = work_dir / "results"
    runs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    if args.skip_setup and not fixture_dir.exists():
        raise FileNotFoundError(f"--skip-setup was provided but fixture directory does not exist: {fixture_dir}")

    fixture = ensure_fixture(
        fixture_dir=fixture_dir,
        preset_name=args.preset,
        seed=args.seed,
        force_rebuild=args.force_rebuild and not args.skip_setup,
    )
    if args.setup_only:
        return {
            "fixture_only": True,
            "fixture_dir": str(fixture_dir),
            "fixture_metadata": fixture["metadata"],
        }

    label = args.label or _git_short_sha(repo_root) or "unlabeled"
    wrapper_dir = work_dir / "bin"
    _ensure_zipstrain_wrapper(wrapper_dir, args.zipstrain_command)
    run_env = os.environ.copy()
    run_env["PATH"] = f"{wrapper_dir}{os.pathsep}{run_env.get('PATH', '')}"
    benchmark_runs: list[dict[str, object]] = []
    for rep_idx in range(args.repetitions):
        run_dir = runs_dir / f"{label}_rep_{rep_idx + 1:02d}"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        cmd = _build_profile_command(
            args.zipstrain_command,
            fixture=fixture,
            output_dir=run_dir,
            num_chunks=args.num_chunks,
            max_concurrency=args.max_concurrency,
        )
        timing = _run_timed_command(cmd, cwd=repo_root, env=run_env)
        outputs = _summarize_profile_outputs(run_dir)
        benchmark_runs.append(
            {
                "rep": rep_idx + 1,
                "command": shlex.join(cmd),
                **timing,
                **outputs,
            }
        )
        if not args.keep_run_dirs:
            shutil.rmtree(run_dir)

    wall_values = [float(run["wall_seconds"]) for run in benchmark_runs]
    rss_values = [run["max_rss_kb"] for run in benchmark_runs if run["max_rss_kb"] is not None]
    result = {
        "label": label,
        "zipstrain_command": args.zipstrain_command,
        "work_dir": str(work_dir),
        "fixture_dir": str(fixture_dir),
        "fixture_metadata": fixture["metadata"],
        "num_chunks": args.num_chunks,
        "max_concurrency": args.max_concurrency,
        "repetitions": args.repetitions,
        "runs": benchmark_runs,
        "summary": {
            "min_wall_seconds": min(wall_values),
            "mean_wall_seconds": sum(wall_values) / len(wall_values),
            "max_wall_seconds": max(wall_values),
            "min_max_rss_kb": min(rss_values) if rss_values else None,
            "mean_max_rss_kb": (sum(rss_values) / len(rss_values)) if rss_values else None,
            "max_max_rss_kb": max(rss_values) if rss_values else None,
        },
    }
    output_json = results_dir / f"{label}.json"
    output_json.write_text(json.dumps(result, indent=2) + "\n")
    result["result_json"] = str(output_json)
    return result


def main() -> int:
    args = parse_args()
    result = benchmark_profile(args)
    if result.get("fixture_only"):
        print(f"fixture_dir={result['fixture_dir']}")
        metadata = result["fixture_metadata"]
        print(
            "fixture_positions="
            f"{metadata['position_count']} fixture_reads={metadata['total_reads']} "
            f"approx_depth={metadata['approx_mean_depth']}"
        )
        return 0

    summary = result["summary"]
    metadata = result["fixture_metadata"]
    print(f"label={result['label']}")
    print(f"result_json={result['result_json']}")
    print(
        "fixture_positions="
        f"{metadata['position_count']} fixture_reads={metadata['total_reads']} "
        f"contigs={metadata['contig_count']} approx_depth={metadata['approx_mean_depth']}"
    )
    print(
        "wall_seconds "
        f"min={summary['min_wall_seconds']:.3f} "
        f"mean={summary['mean_wall_seconds']:.3f} "
        f"max={summary['max_wall_seconds']:.3f}"
    )
    if summary["mean_max_rss_kb"] is not None:
        print(
            "max_rss_kb "
            f"min={int(summary['min_max_rss_kb'])} "
            f"mean={summary['mean_max_rss_kb']:.1f} "
            f"max={int(summary['max_max_rss_kb'])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
