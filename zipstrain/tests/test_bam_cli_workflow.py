from __future__ import annotations

import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys

from click.testing import CliRunner
import polars as pl
import pytest

from zipstrain import cli as cli_module


READ_DEPTH = 3
REFERENCE_ID = "toy_reference"
GENE_DB_ID = "toy_genes"
SCAFFOLD_TO_GENOME = {
    "contigA": "genome_a",
    "contigB": "genome_b",
}
SAMPLE_SEQUENCES = {
    "sample_alpha": {
        "contigA": "ACGTACGTAA",
        "contigB": "TTGCAACGTA",
    },
    "sample_beta": {
        "contigA": "ACGTACGTAA",
        "contigB": "TTGCATCGTA",
    },
}


def _require_samtools() -> str:
    samtools = shutil.which("samtools")
    assert samtools is not None, "samtools is required for BAM-backed CLI workflow tests."
    return samtools


def _write_zipstrain_wrapper(bin_dir: Path) -> None:
    src_root = Path(cli_module.__file__).resolve().parents[1]
    wrapper = bin_dir / "zipstrain"
    wrapper.write_text(
        "\n".join(
            [
                f"#!{sys.executable}",
                "import sys",
                f"sys.path.insert(0, {str(src_root)!r})",
                "from zipstrain.cli import cli",
                "cli()",
                "",
            ]
        )
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)


def _write_reference_fasta(reference_fasta: Path) -> None:
    reference_fasta.write_text(
        ">contigA\nNNNNNNNNNN\n>contigB\nNNNNNNNNNN\n"
    )


def _write_gene_fasta(gene_fasta: Path) -> None:
    gene_fasta.write_text(
        ">contigA_1 # 1 # 10 # 1 # ID=1_1;partial=00\nACGTACGTAA\n"
        ">contigB_1 # 1 # 10 # 1 # ID=2_1;partial=00\nTTGCAACGTA\n"
    )


def _write_stb(stb_file: Path) -> None:
    stb_file.write_text(
        "contigA\tgenome_a\n"
        "contigB\tgenome_b\n"
    )


def _write_bam(sample_name: str, sequences: dict[str, str], bam_file: Path) -> None:
    sam_file = bam_file.with_suffix(".sam")
    unsorted_bam = bam_file.with_suffix(".unsorted.bam")
    header_lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:contigA\tLN:10",
        "@SQ\tSN:contigB\tLN:10",
    ]
    read_lines: list[str] = []
    for contig in ("contigA", "contigB"):
        read_seq = sequences[contig]
        qual = "I" * len(read_seq)
        for read_idx in range(READ_DEPTH):
            read_lines.append(
                f"{sample_name}_{contig}_{read_idx}\t0\t{contig}\t1\t60\t10M\t*\t0\t0\t{read_seq}\t{qual}"
            )
    sam_file.write_text("\n".join(header_lines + read_lines) + "\n")
    subprocess.run(
        ["samtools", "view", "-bS", str(sam_file), "-o", str(unsorted_bam)],
        check=True,
    )
    subprocess.run(
        ["samtools", "sort", "-o", str(bam_file), str(unsorted_bam)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(bam_file)], check=True)


def _expected_profile_frame(sample_name: str) -> pl.DataFrame:
    rows: list[dict[str, object]] = []
    for chrom in ("contigA", "contigB"):
        genome = SCAFFOLD_TO_GENOME[chrom]
        gene = f"{chrom}_1"
        for pos, base in enumerate(SAMPLE_SEQUENCES[sample_name][chrom], start=1):
            rows.append(
                {
                    "chrom": chrom,
                    "genome": genome,
                    "gene": gene,
                    "pos": pos,
                    "A": READ_DEPTH if base == "A" else 0,
                    "C": READ_DEPTH if base == "C" else 0,
                    "G": READ_DEPTH if base == "G" else 0,
                    "T": READ_DEPTH if base == "T" else 0,
                    "ref_base_bitmask": 0,
                }
            )
    return pl.DataFrame(rows).sort(["chrom", "pos"])


def _workflow_input_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    _require_samtools()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    _write_zipstrain_wrapper(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    reference_fasta = tmp_path / "reference.fna"
    gene_fasta = tmp_path / "genes.fna"
    stb_file = tmp_path / "reference.stb"
    prep_dir = tmp_path / "prep"
    profiles_dir = tmp_path / "profiles"
    prep_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir.mkdir(parents=True, exist_ok=True)

    _write_reference_fasta(reference_fasta)
    _write_gene_fasta(gene_fasta)
    _write_stb(stb_file)

    for sample_name, sequences in SAMPLE_SEQUENCES.items():
        _write_bam(sample_name, sequences, tmp_path / f"{sample_name}.bam")

    runner = CliRunner()
    prep_result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "prepare_profiling",
            "--reference-fasta",
            str(reference_fasta),
            "--gene-fasta",
            str(gene_fasta),
            "--stb-file",
            str(stb_file),
            "--output-dir",
            str(prep_dir),
        ],
    )
    assert prep_result.exit_code == 0, prep_result.output

    return {
        "reference_fasta": reference_fasta,
        "gene_fasta": gene_fasta,
        "stb_file": stb_file,
        "null_model": prep_dir / "null_model.parquet",
        "profiling_contract": prep_dir / "profiling_contract.json",
        "prep_dir": prep_dir,
        "profiles_dir": profiles_dir,
        "sample_alpha_bam": tmp_path / "sample_alpha.bam",
        "sample_beta_bam": tmp_path / "sample_beta.bam",
    }


def _profile_dir_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    paths = _workflow_input_fixture(tmp_path, monkeypatch)
    runner = CliRunner()

    for sample_name in SAMPLE_SEQUENCES:
        result = runner.invoke(
            cli_module.cli,
            [
                "utilities",
                "profile-single",
                "--reference-fasta",
                str(paths["reference_fasta"]),
                "--bed-file",
                str(paths["prep_dir"] / "genomes_bed_file.bed"),
                "--bam-file",
                str(paths[f"{sample_name}_bam"]),
                "--stb-file",
                str(paths["stb_file"]),
                "--null-model",
                str(paths["null_model"]),
                "--profiling-contract",
                str(paths["profiling_contract"]),
                "--gene-range-table",
                str(paths["prep_dir"] / "gene_range_table.tsv"),
                "--num-chunks",
                "2",
                "--max-concurrency",
                "2",
                "--output-dir",
                str(paths["profiles_dir"]),
            ],
        )
        assert result.exit_code == 0, result.output

    paths["sample_alpha_profile"] = paths["profiles_dir"] / "sample_alpha_profile.parquet"
    paths["sample_beta_profile"] = paths["profiles_dir"] / "sample_beta_profile.parquet"
    return paths


def test_cli_profile_from_real_bams_matches_expected_counts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = _profile_dir_fixture(tmp_path, monkeypatch)

    gene_ranges = pl.read_csv(
        paths["prep_dir"] / "gene_range_table.tsv",
        separator="\t",
        has_header=False,
        new_columns=["gene", "scaffold", "start", "end"],
    ).sort("gene")
    assert gene_ranges.to_dicts() == [
        {"gene": "contigA_1", "scaffold": "contigA", "start": 1, "end": 10},
        {"gene": "contigB_1", "scaffold": "contigB", "start": 1, "end": 10},
    ]

    bed_rows = pl.read_csv(
        paths["prep_dir"] / "genomes_bed_file.bed",
        separator="\t",
        has_header=False,
        new_columns=["scaffold", "start", "end"],
    ).sort("scaffold")
    assert bed_rows.to_dicts() == [
        {"scaffold": "contigA", "start": 0, "end": 10},
        {"scaffold": "contigB", "start": 0, "end": 10},
    ]

    alpha_profile = pl.read_parquet(paths["sample_alpha_profile"]).sort(["chrom", "pos"])
    beta_profile = pl.read_parquet(paths["sample_beta_profile"]).sort(["chrom", "pos"])

    assert alpha_profile.equals(_expected_profile_frame("sample_alpha"))
    assert beta_profile.equals(_expected_profile_frame("sample_beta"))

    contract = cli_module.ut.read_profile_contract_file(paths["profiling_contract"])
    alpha_metadata = pl.read_parquet_metadata(paths["sample_alpha_profile"])
    for logical_name, metadata_key in cli_module.ut.PROFILE_CONTRACT_METADATA_KEYS.items():
        assert alpha_metadata[metadata_key] == contract[logical_name]

    alpha_gene_stats = pl.read_parquet(paths["profiles_dir"] / "sample_alpha_gene_stats.parquet").sort("gene")
    beta_gene_stats = pl.read_parquet(paths["profiles_dir"] / "sample_beta_gene_stats.parquet").sort("gene")
    for gene_stats in (alpha_gene_stats, beta_gene_stats):
        assert gene_stats.select("gene", "length").to_dicts() == [
            {"gene": "contigA_1", "length": 10},
            {"gene": "contigB_1", "length": 10},
        ]
        assert gene_stats["breadth"].to_list() == pytest.approx([1.0, 1.0])
        assert gene_stats["coverage"].to_list() == pytest.approx([float(READ_DEPTH), float(READ_DEPTH)])
        assert "ref_ani" in gene_stats.columns

    alpha_genome_stats = pl.read_parquet(paths["profiles_dir"] / "sample_alpha_genome_stats.parquet").sort("genome")
    alpha_by_genome = alpha_genome_stats.rows_by_key("genome", unique=True, named=True)
    assert alpha_by_genome["genome_a"]["coverage"] == pytest.approx(float(READ_DEPTH))
    assert alpha_by_genome["genome_b"]["coverage"] == pytest.approx(float(READ_DEPTH))
    assert alpha_by_genome["genome_a"]["breadth"] == pytest.approx(1.0)
    assert alpha_by_genome["genome_b"]["breadth"] == pytest.approx(1.0)
    assert "ref_ani" in alpha_genome_stats.columns


def test_cli_profile_workflow_from_real_bams(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = _workflow_input_fixture(tmp_path, monkeypatch)
    runner = CliRunner()
    input_table = tmp_path / "bam_inputs.csv"
    input_table.write_text(
        "\n".join(
            [
                "sample_name,bamfile",
                f"sample_alpha,{paths['sample_alpha_bam']}",
                f"sample_beta,{paths['sample_beta_bam']}",
            ]
        )
        + "\n"
    )
    run_dir = tmp_path / "profile_run"
    result = runner.invoke(
        cli_module.cli,
            [
                "profile",
                "--input-table",
                str(input_table),
                "--reference-fasta",
                str(paths["reference_fasta"]),
                "--stb-file",
                str(paths["stb_file"]),
            "--null-model",
            str(paths["null_model"]),
            "--gene-range-table",
            str(paths["prep_dir"] / "gene_range_table.tsv"),
            "--bed-file",
            str(paths["prep_dir"] / "genomes_bed_file.bed"),
            "--genome-length-file",
            str(paths["prep_dir"] / "genome_lengths.parquet"),
            "--run-dir",
            str(run_dir),
            "--num-procs",
            "2",
            "--task-per-batch",
            "2",
            "--max-concurrent-batches",
            "1",
            "--poll-interval",
            "1",
            "--execution-mode",
            "local",
            "--container-engine",
            "local",
        ],
    )
    assert result.exit_code == 0, result.output

    # After a successful run the batch scaffolding is flattened away: sample
    # folders sit directly under run_dir, and batch logs land in profiling_assets/log.
    assert not (run_dir / "batch_0").exists()
    log_dir = run_dir / "profiling_assets" / "log"
    assert (log_dir / "batch_events.log").exists()

    for sample_name in SAMPLE_SEQUENCES:
        sample_dir = run_dir / sample_name
        assert sample_dir.exists()
        # The real output parquets live at the top of the sample dir; everything
        # else is in intermediate_files/.
        top_level = {p.name for p in sample_dir.iterdir()}
        assert top_level == {
            f"{sample_name}_profile.parquet",
            f"{sample_name}_genome_stats.parquet",
            f"{sample_name}_gene_stats.parquet",
            f"{sample_name}_SNVs.parquet",
            "intermediate_files",
        }
        assert (sample_dir / "intermediate_files").is_dir()
        assert (sample_dir / "intermediate_files" / ".status").read_text().strip() == "done"
        prof = pl.read_parquet(sample_dir / f"{sample_name}_profile.parquet").sort(["chrom", "pos"])
        assert prof.equals(_expected_profile_frame(sample_name))
        genome_stats = pl.read_parquet(sample_dir / f"{sample_name}_genome_stats.parquet").sort("genome")
        by_genome = genome_stats.rows_by_key("genome", unique=True, named=True)
        assert by_genome["genome_a"]["coverage"] == pytest.approx(float(READ_DEPTH))
        assert by_genome["genome_b"]["coverage"] == pytest.approx(float(READ_DEPTH))
        assert "ref_ani" in genome_stats.columns
        gene_stats = pl.read_parquet(sample_dir / f"{sample_name}_gene_stats.parquet").sort("gene")
        assert gene_stats["breadth"].to_list() == pytest.approx([1.0, 1.0])
        assert "ref_ani" in gene_stats.columns


def test_cli_standard_compare_workflow_from_real_bams(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = _profile_dir_fixture(tmp_path, monkeypatch)
    runner = CliRunner()

    profile_db_csv = tmp_path / "profiles.csv"
    profile_db_csv.write_text(
        "\n".join(
            [
                "profile_name,profile_location",
                f"sample_alpha,{paths['sample_alpha_profile']}",
                f"sample_beta,{paths['sample_beta_profile']}",
            ]
        )
        + "\n"
    )
    profile_db_parquet = tmp_path / "profile_db.parquet"
    build_db_result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "build-profile-db",
            "--profile-db-csv",
            str(profile_db_csv),
            "--output-file",
            str(profile_db_parquet),
        ],
    )
    assert build_db_result.exit_code == 0, build_db_result.output
    profile_db = pl.read_parquet(profile_db_parquet).sort("profile_name")
    assert profile_db.to_dicts() == [
        {
            "profile_name": "sample_alpha",
            "profile_location": str(paths["sample_alpha_profile"]),
        },
        {
            "profile_name": "sample_beta",
            "profile_location": str(paths["sample_beta_profile"]),
        },
    ]

    pair_table = tmp_path / "pairs.parquet"
    pair_result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "generate-sample-pair",
            "--profile-dir",
            str(paths["profiles_dir"]),
            "--output-file",
            str(pair_table),
        ],
    )
    assert pair_result.exit_code == 0, pair_result.output
    pairs = pl.read_parquet(pair_table)
    assert pairs.height == 1

    genome_compare_out = tmp_path / "genome_compare.parquet"
    compare_result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "chunk-genome-compare",
            "--pair-table",
            str(pair_table),
            "--stb-file",
            str(paths["stb_file"]),
            "--output-file",
            str(genome_compare_out),
            "--workers",
            "1",
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--calculate",
            "ani+ibs",
            "--engine",
            "polars",
        ],
    )
    assert compare_result.exit_code == 0, compare_result.output

    genome_compare = pl.read_parquet(genome_compare_out).sort("genome")
    by_genome = genome_compare.rows_by_key("genome", unique=True, named=True)
    assert by_genome["genome_a"]["total_positions"] == 10
    assert by_genome["genome_a"]["share_allele_pos"] == 10
    assert by_genome["genome_a"]["genome_ani"] == pytest.approx(100.0)
    assert by_genome["genome_a"]["max_consecutive_length"] == 10
    assert by_genome["genome_b"]["total_positions"] == 10
    assert by_genome["genome_b"]["share_allele_pos"] == 9
    assert by_genome["genome_b"]["genome_ani"] == pytest.approx(90.0)
    assert by_genome["genome_b"]["max_consecutive_length"] == 5
    assert {by_genome["genome_a"]["sample_1"], by_genome["genome_a"]["sample_2"]} == {
        "sample_alpha",
        "sample_beta",
    }

    gene_compare_out = tmp_path / "gene_compare.parquet"
    gene_result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "single_compare_gene",
            "--profile-location-1",
            str(paths["sample_alpha_profile"]),
            "--profile-location-2",
            str(paths["sample_beta_profile"]),
            "--stb-file",
            str(paths["stb_file"]),
            "--min-cov",
            "1",
            "--min-gene-compare-len",
            "1",
            "--scope",
            "all:all",
            "--output-file",
            str(gene_compare_out),
            "--engine",
            "polars",
        ],
    )
    assert gene_result.exit_code == 0, gene_result.output

    gene_compare = pl.read_parquet(gene_compare_out).sort(["genome", "gene"])
    gene_by_name = gene_compare.rows_by_key("gene", unique=True, named=True)
    assert gene_by_name["contigA_1"]["genome"] == "genome_a"
    assert gene_by_name["contigA_1"]["total_positions"] == 10
    assert gene_by_name["contigA_1"]["share_allele_pos"] == 10
    assert gene_by_name["contigA_1"]["ani"] == pytest.approx(100.0)
    assert gene_by_name["contigB_1"]["genome"] == "genome_b"
    assert gene_by_name["contigB_1"]["total_positions"] == 10
    assert gene_by_name["contigB_1"]["share_allele_pos"] == 9
    assert gene_by_name["contigB_1"]["ani"] == pytest.approx(90.0)
