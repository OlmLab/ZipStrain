#!/usr/bin/env python3
"""Measure ZipStrain V2 comparison on repeatable synthetic profiles."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import version
from itertools import combinations
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
import polars as pl
from polars.testing import assert_frame_equal

from zipstrain import __version__, codons, gene_info, utils
from zipstrain.profile import PROFILE_SORTED_METADATA_KEY, PROFILE_SORTED_METADATA_VALUE


SUITE_VERSION = 2
THREADS = 2
CALCULATIONS = {
    "genome": "genome_ani+max_identical_stretch",
    "mixed": "genome_ani+max_identical_stretch+gene_ani",
}
LAYOUTS = {
    "single": (1, 120_000),
    "multiple": (8, 30_000),
}
_CODON_TEXT = np.array([a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"])

# Run the normal Click entry point in a fresh process. The table method and
# NumPy matrix method do their work in this process and its threads.
TIMED_CLI = """
import json
from pathlib import Path
import resource
import sys
from zipstrain.cli import cli

cli.main(args=sys.argv[2:], prog_name="zipstrain", standalone_mode=False)
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
Path(sys.argv[1]).write_text(json.dumps(rss if sys.platform == "darwin" else rss * 1024))
"""


def simulated_codon_counts(genes: pl.DataFrame, calls: np.ndarray, alternate: np.ndarray,
                           depth: np.ndarray, keep: np.ndarray, rng: np.random.Generator) -> pl.DataFrame:
    """Write joint counts from two simulated haplotypes at a subset of covered codons."""
    indices = (genes["start"].to_numpy()[:, None] - 1 + np.arange(0, 300, 3)).reshape(-1)
    callable_codons = keep[indices] & keep[indices + 1] & keep[indices + 2]
    callable_codons &= rng.random(len(indices)) < 0.25
    gene_id = np.repeat(genes["gene_id"].to_numpy(), 100)[callable_codons]
    codon_index = np.tile(np.arange(100, dtype=np.uint32), genes.height)[callable_codons]
    bases = indices[callable_codons, None] + np.arange(3)
    primary = calls[bases] @ np.array([16, 4, 1])
    secondary = alternate[bases] @ np.array([16, 4, 1])
    minor_depth = depth[bases] // 3
    major_reads = (depth[bases] - minor_depth).min(axis=1)
    minor_reads = minor_depth.min(axis=1)
    polymorphic = primary != secondary
    return pl.DataFrame({
        "gene_id": np.concatenate((gene_id, gene_id[polymorphic])),
        "codon_index": np.concatenate((codon_index, codon_index[polymorphic])),
        "codon": _CODON_TEXT[np.concatenate((primary, secondary[polymorphic]))],
        "count": np.concatenate((major_reads + np.where(polymorphic, 0, minor_reads), minor_reads[polymorphic])),
    }, schema=codons.COUNT_SCHEMA)


def make_fixture(root: Path, genomes: int, length: int, samples: int, *, selection=False) -> list[tuple[str, Path]]:
    """Write deterministic position profiles and optional phased evidence outside the timer."""
    root.mkdir()
    rng = np.random.default_rng(8291)
    alphabet = np.asarray(list("ACGT"))
    positions = np.arange(1, length + 1, dtype=np.int32)
    scaffolds = [f"scaffold_{index:02d}" for index in range(genomes)]
    genomes_by_scaffold = {chrom: f"genome_{index:02d}" for index, chrom in enumerate(scaffolds)}
    references = rng.integers(0, 4, size=(genomes, length), dtype=np.int32)

    reference = root / "reference.fna"
    genes = root / "genes.fna"
    stb = root / "reference.stb"
    with reference.open("w") as fasta, genes.open("w") as cds, stb.open("w") as mapping:
        for chrom, bases in zip(scaffolds, references, strict=True):
            sequence = "".join(alphabet[bases])
            fasta.write(f">{chrom}\n{sequence}\n")
            mapping.write(f"{chrom}\t{genomes_by_scaffold[chrom]}\n")
            for index, start in enumerate(range(0, length - 299, 270), 1):
                cds.write(f">{chrom}_{index} # {start + 1} # {start + 300} # 1 # partial=00\n")
                cds.write(sequence[start:start + 300] + "\n")
    catalog = gene_info.write_gene_info(reference, genes, stb, root / "gene_info.parquet")
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata[PROFILE_SORTED_METADATA_KEY] = PROFILE_SORTED_METADATA_VALUE
    if selection:
        genes_by_scaffold = gene_info.read_gene_info(catalog).partition_by("chrom", as_dict=True)
        evidence_rng = np.random.default_rng(2375)
        (root / "regions.bed").write_text("".join(f"{chrom}\t0\t{length}\n" for chrom in scaffolds))

    profiles = []
    for sample in range(samples):
        name = f"sample_{sample:02d}"
        frames = []
        evidence = []
        for chrom, bases in zip(scaffolds, references, strict=True):
            substitutions = rng.random(length) < 0.015
            calls = (bases + substitutions * rng.integers(1, 4, size=length)) % 4
            depth = rng.integers(5, 36, size=length, dtype=np.int32)
            keep = rng.random(length) >= 0.04
            if selection:
                mutations = evidence_rng.random(length) < 0.035
                alternate = (calls + mutations * evidence_rng.integers(1, 4, size=length)) % 4
                minor_reads = depth // 3
                position_counts = {
                    base: ((depth - minor_reads) * (calls == index) + minor_reads * (alternate == index))[keep]
                    for index, base in enumerate(alphabet)
                }
                evidence.append(simulated_codon_counts(
                    genes_by_scaffold[(chrom,)], calls, alternate, depth, keep, evidence_rng,
                ))
            else:
                position_counts = {base: (depth * (calls == index))[keep] for index, base in enumerate(alphabet)}
            frames.append(pl.DataFrame({
                "pos": positions[keep],
                **position_counts,
            }).with_columns(
                pl.lit(chrom).alias("chrom"),
                pl.lit(genomes_by_scaffold[chrom]).alias("genome"),
            ).select("chrom", "genome", "pos", "A", "C", "G", "T"))
        path = root / f"{name}_profile.parquet"
        pl.concat(frames).write_parquet(path, compression="zstd", metadata=metadata)
        if selection:
            pl.concat(evidence).sort(["gene_id", "codon_index", "codon"]).write_parquet(
                root / f"{name}_codon_counts.parquet", compression="zstd", metadata={
                    **utils._read_custom_parquet_metadata(catalog),
                    "zipstrain_codon_evidence": "read_phased_v1",
                    "zipstrain_source_profile_sha256": utils.sha256_file(path),
                },
            )
        profiles.append((name, path))
    with (root / "profiles.csv").open("w", newline="") as handle:
        csv.writer(handle).writerows([("profile_name", "profile_location"), *profiles])
    return profiles


def timed_compare(args: list[str], output: Path, expected_pairs: set[tuple[str, str]],
                  expected_rows_per_pair: int, mixed: bool, *, selection=False) -> dict:
    """Measure a clean CLI run and verify its result outside the timed interval."""
    rss_path = output.parent.with_suffix(".rss.json")
    env = dict(os.environ, ZIPSTRAIN_CPU_BUDGET=str(THREADS), POLARS_MAX_THREADS=str(THREADS),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    start = perf_counter()
    process = subprocess.run(
        [sys.executable, "-c", TIMED_CLI, str(rss_path), *args],
        env=env, capture_output=True, text=True,
    )
    seconds = perf_counter() - start
    if process.returncode:
        raise RuntimeError(f"Comparison failed ({process.returncode})\n{process.stdout}\n{process.stderr}")
    peak_bytes = json.loads(rss_path.read_text())
    if peak_bytes <= 0:
        raise ValueError("The OS did not report a positive peak RSS")
    table = pl.read_parquet(output)
    keys = ["sample_1", "sample_2", "genome"] + (["gene_id"] if mixed else [])
    required = {*keys, "genome_ani", "max_consecutive_length"}
    if mixed:
        required.add("gene_ani")
    if selection:
        required.update(gene_info.SELECTION_COLUMNS)
    if not required.issubset(table.columns) or ("gene_id" in table.columns) != mixed:
        raise ValueError(f"Unexpected columns in {output}: {table.columns}")
    if selection and not table.select(pl.any_horizontal(pl.col("pN", "pS") > 0).any()).item():
        raise ValueError(f"No polymorphic selection events measured in {output}")
    pairs = table.group_by("sample_1", "sample_2").len()
    if (
        set(pairs.select("sample_1", "sample_2").iter_rows()) != expected_pairs
        or pairs["len"].to_list() != [expected_rows_per_pair] * len(expected_pairs)
        or table.select(keys).n_unique() != table.height
        or list(output.parent.glob("*.parquet")) != [output]
    ):
        raise ValueError(f"Unexpected comparison rows or output files in {output.parent}")
    return {"seconds": seconds, "peak_rss_bytes": peak_bytes}


def run_suite(root: Path, repetitions: int, quick: bool) -> list[dict]:
    """Benchmark every layout and metric set using fresh outputs for each run."""
    results = []
    samples = 3 if quick else 4
    for layout, (genomes, length) in LAYOUTS.items():
        if quick:
            length = 3_000
        fixture = root / layout
        profiles = make_fixture(fixture, genomes, length, samples)
        selection_fixture = root / f"{layout}-selection"
        phased_profiles = make_fixture(selection_fixture, genomes, length, samples, selection=True)
        workloads = [("table", metric_set, None) for metric_set in CALCULATIONS]
        workloads += [
            (method, "selection", integration)
            for integration in codons.ALLELE_INTEGRATION_MODES
            for method in ("table", "matrix")
        ]
        table_selection_outputs = {}
        for method, metric_set, integration in workloads:
            selection = metric_set == "selection"
            inputs = selection_fixture if selection else fixture
            pairs = list(combinations(phased_profiles if selection else profiles, 2))
            names = {(a, b) for (a, _), (b, _) in pairs}
            calculate = CALCULATIONS["mixed" if selection else metric_set]
            mixed = metric_set != "genome"
            rows = genomes * ((length - 300) // 270 + 1) if mixed else genomes
            common = ["--preprofile-dir", str(inputs), "--calculate", calculate,
                      "--min-cov", "5", "--min-gene-compare-len", "100"]
            if selection:
                common += ["--dnds", "--allele-integration", integration]
            if method == "matrix":
                common += ["--stb-file", str(inputs / "reference.stb"),
                           "--bed-file", str(inputs / "regions.bed"), "--backend", "numpy",
                           "--compare-kernel", "anchor", "--memory-limit-gb", "1"]
            compare_args = ["compare", "--method", method, "--no-csv", *common]
            case = {
                "layout": layout, "metric_set": metric_set, "method": method,
                "allele_integration": integration, "genomes": genomes,
                "bases_per_genome": length, "samples": samples, "pairs": len(pairs), "runs": [],
            }
            label = "-".join(value for value in (layout, method, metric_set, integration) if value)
            for repeat in range(repetitions):
                run_dir = root / f"{label}-{repeat}"
                run_dir.mkdir()
                collection = run_dir / "collection"
                output = collection / "all_comparisons.parquet"
                total = timed_compare(
                    [*compare_args, "--profile-db", str(inputs / "profiles.csv"), "--run-dir", str(collection)],
                    output, names, rows, mixed, selection=selection,
                )
                if selection and method == "table":
                    table_selection_outputs[integration, repeat] = output
                elif selection:
                    keys = ["sample_1", "sample_2", "genome", "gene_id"]
                    assert_frame_equal(
                        pl.read_parquet(output).sort(keys),
                        pl.read_parquet(table_selection_outputs[integration, repeat]).sort(keys),
                    )
                isolated = []
                for (a, profile_a), (b, profile_b) in pairs:
                    pair_dir = run_dir / f"{a}-{b}"
                    if selection:
                        pair_csv = pair_dir.with_suffix(".csv")
                        with pair_csv.open("w", newline="") as handle:
                            csv.writer(handle).writerows([
                                ("profile_name", "profile_location"), (a, profile_a), (b, profile_b),
                            ])
                        pair_args = [*compare_args, "--profile-db", str(pair_csv), "--run-dir", str(pair_dir)]
                        pair_output = pair_dir / "all_comparisons.parquet"
                    else:
                        pair_args = ["utilities", "compare-pair", "--profile-1", str(profile_a),
                                     "--profile-2", str(profile_b), "--output-dir", str(pair_dir), *common]
                        pair_output = pair_dir / "comparisons.parquet"
                    observation = timed_compare(pair_args, pair_output, {(a, b)}, rows, mixed, selection=selection)
                    isolated.append({"sample_1": a, "sample_2": b, **observation})
                case["runs"].append({"collection": total, "isolated_pairs": isolated})
            results.append(case)
            print(f"Measured {label} ({len(pairs)} pairs × {repetitions} repeats)", flush=True)
    return results


def main() -> None:
    """Run the comparison suite and write its raw timings and environment as JSON."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--quick", action="store_true", help="Use small fixtures to smoke-test the suite.")
    parser.add_argument("--release", action="store_true", help="Label measurements from the release runner.")
    args = parser.parse_args()
    if args.repetitions < 1 or (args.release and args.quick):
        parser.error("Use positive repetitions and full-sized fixtures for releases.")
    if sys.platform not in {"linux", "darwin"}:
        parser.error("Peak RSS measurement supports Linux and macOS.")
    cpu = platform.processor() or platform.machine()
    if sys.platform == "linux":
        cpu = next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
                    if line.startswith("model name")), cpu)

    with TemporaryDirectory(prefix="zipstrain-benchmarks-") as temporary:
        cases = run_suite(Path(temporary), args.repetitions, args.quick)
    result = {
        "suite_version": SUITE_VERSION, "quick": args.quick,
        "series": "release" if args.release else "local",
        "version": __version__,
        "revision": os.environ.get("GITHUB_SHA"),
        "measured_at": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "platform": platform.platform(), "cpu": cpu, "threads": THREADS,
            "python": platform.python_version(), "polars": pl.__version__,
            "numpy": np.__version__, "h5py": version("h5py"),
            "runner_image": os.environ.get("ImageVersion"),
        },
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
