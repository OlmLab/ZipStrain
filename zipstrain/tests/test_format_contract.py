"""Pin the V2 output formats from a real run. See "Output formats frozen for 2.0" in docs/v2.md.

A failure here means a released format changed. If that is intended, change the
``zipstrain_profile_schema`` value and the guide along with these literals.
"""

import shutil

import polars as pl
import pytest
from click.testing import CliRunner

from zipstrain import cli, gene_info as gi

from .test_v2 import write_bam

STRINGS, INT64, FLOAT64 = pl.String, pl.Int64, pl.Float64
CATALOG_IDENTITY = {"zipstrain_profile_schema", "zipstrain_reference_id", "zipstrain_annotation_id"}
PROFILE_KEYS = {"zipstrain_profile_filters", "zipstrain_selection_method", "zipstrain_reference_hash",
                "zipstrain_gene_hash", "zipstrain_null_model_hash", "zipstrain_stb_hash"} | CATALOG_IDENTITY
COMPARISON_KEYS = CATALOG_IDENTITY | {
    "zipstrain_compare_calculate", "zipstrain_compare_kind", "zipstrain_compare_engine",
    "zipstrain_compare_ani_method", "zipstrain_compare_min_cov", "zipstrain_compare_min_gene_compare_len",
    "zipstrain_compare_allele_integration", "zipstrain_compare_selection_method",
    "zipstrain_compare_reference_id", "zipstrain_compare_annotation_id",
    "zipstrain_allele_integration", "zipstrain_selection_method",
}
COMPARISON_COLUMNS = {
    "sample_1": STRINGS, "sample_2": STRINGS, "genome": STRINGS, "gene_id": pl.UInt64,
    "gene_total_positions": INT64, "gene_share_allele_pos": INT64, "gene_ani": FLOAT64,
    "dN": FLOAT64, "dS": FLOAT64, "pN": FLOAT64, "pS": FLOAT64,
    "total_positions": INT64, "share_allele_pos": INT64, "genome_ani": FLOAT64,
    "max_identical_stretch": INT64,
}
FORMATS = {
    "assets/gene_info.parquet": {
        "gene_id": pl.UInt64, "gene": STRINGS, "genome": STRINGS, "chrom": STRINGS, "start": INT64, "end": INT64,
        "length": INT64, "strand": pl.Int8, "phase": pl.UInt8, "genetic_code": pl.UInt8,
        "partial_5p": pl.Boolean, "partial_3p": pl.Boolean,
    },
    "assets/reference_codons.parquet": {"gene_id": pl.UInt64, "codon_index": pl.UInt32, "reference_codon": STRINGS},
    "sample/reads_profile.parquet": {
        "chrom": STRINGS, "genome": STRINGS, "pos": pl.Int32, "A": pl.Int32, "C": pl.Int32, "G": pl.Int32,
        "T": pl.Int32, "ref_base_bitmask": pl.UInt8,
    },
    "sample/reads_codon_counts.parquet": {
        "gene_id": pl.UInt64, "codon_index": pl.UInt32, "codon": STRINGS, "count": pl.UInt64,
    },
    "sample/reads_gene_stats.parquet": {
        "gene_id": pl.UInt64, "genome": STRINGS, "chrom": STRINGS, "start": INT64, "end": INT64, "length": INT64,
        "coverage": FLOAT64, "breadth": FLOAT64, "ref_ani": FLOAT64, "dN": FLOAT64, "dS": FLOAT64, "dN_dS": FLOAT64,
        "pN": FLOAT64, "pS": FLOAT64, "pN_pS": FLOAT64,
    },
    "table/all_comparisons.parquet": COMPARISON_COLUMNS,
    "matrix/all_comparisons.parquet": COMPARISON_COLUMNS,
}
REQUIRED_METADATA = {
    "assets/gene_info.parquet": CATALOG_IDENTITY | {"zipstrain_reference_codons_sha256"},
    "assets/reference_codons.parquet": CATALOG_IDENTITY,
    "sample/reads_profile.parquet": PROFILE_KEYS | {"zipstrain_sorted_by"},
    "sample/reads_codon_counts.parquet": CATALOG_IDENTITY | {
        "zipstrain_codon_evidence", "zipstrain_codon_filters", "zipstrain_codon_quality",
        "zipstrain_selection_method", "zipstrain_source_profile_sha256",
    },
    "sample/reads_gene_stats.parquet": PROFILE_KEYS | {"zipstrain_allele_integration"},
    "sample/reads_genome_stats.parquet": PROFILE_KEYS,
    "table/all_comparisons.parquet": COMPARISON_KEYS,
    "matrix/all_comparisons.parquet": COMPARISON_KEYS,
}


def _invoke(*args):
    result = CliRunner().invoke(cli.cli, [str(arg) for arg in args])
    assert result.exit_code == 0, result.output + repr(result.exception)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """Profile a BAM with collected codons, then compare two samples by both methods."""
    pytest.importorskip("h5py")
    root = tmp_path_factory.mktemp("formats")
    (root/"ref.fa").write_text(">chr1\nAAAGGATTTCCCAAAGGTTAA\n>chr2\nGGAGGAGGA\n")
    (root/"genes.fa").write_text(">chr1_1 # 1 # 9 # 1 # partial=00\nAAAGGATTT\n>chr2_1 # 1 # 9 # 1 # partial=00\nGGAGGAGGA\n")
    (root/"ref.stb").write_text("chr1\tgenomeA\nchr2\tgenomeB\n")
    assets, sample = root/"assets", root/"sample"
    _invoke("utilities", "prepare_profiling", "--reference-fasta", root/"ref.fa", "--gene-fasta", root/"genes.fa",
            "--stb-file", root/"ref.stb", "--output-dir", assets)
    records = [f"read{i}\t0\tchr1\t1\t60\t21M\t*\t0\t0\t{'AAAGGATTTCCCAAAGGTTAA' if i < 10 else 'AAAGGTTTTCCCAAAGGTTAA'}"
               f"\t{'I'*21}\tNM:i:{int(i >= 10)}" for i in range(50)]
    _invoke("utilities", "profile-single", "--reference-fasta", root/"ref.fa", "--bed-file", assets/"genomes_bed_file.bed",
            "--bam-file", write_bam(root, records), "--stb-file", root/"ref.stb", "--null-model", assets/"null_model.parquet",
            "--gene-info-table", assets/"gene_info.parquet", "--output-dir", sample, "--num-chunks", 1,
            "--max-concurrency", 1, "--read-inclusion", "all-mapped", "--min-read-ani", 0, "--collect-codons")
    for name in ("a", "b"):
        shutil.copyfile(sample/"reads_profile.parquet", sample/f"{name}_profile.parquet")
        shutil.copyfile(sample/"reads_codon_counts.parquet", sample/f"{name}_codon_counts.parquet")
    (root/"profiles.csv").write_text(
        f"profile_name,profile_location\na,{sample/'a_profile.parquet'}\nb,{sample/'b_profile.parquet'}\n")
    common = ("compare", "--profile-db", root/"profiles.csv", "--preprofile-dir", assets, "--no-csv")
    _invoke(*common, "--method", "table", "--run-dir", root/"table")
    _invoke(*common, "--method", "matrix", "--stb-file", root/"ref.stb", "--bed-file", assets/"genomes_bed_file.bed",
            "--run-dir", root/"matrix")
    return root


@pytest.mark.parametrize("name", FORMATS)
def test_columns_types_and_order_are_frozen(run, name):
    assert list(pl.read_parquet_schema(run/name).items()) == list(FORMATS[name].items())


@pytest.mark.parametrize("name", REQUIRED_METADATA)
def test_metadata_keys_are_frozen(run, name):
    metadata = pl.read_parquet_metadata(run/name)
    assert not REQUIRED_METADATA[name] - set(metadata)
    assert metadata["zipstrain_profile_schema"] == "2"


def test_only_position_profiles_claim_coordinate_order(run):
    assert pl.read_parquet_metadata(run/"sample/reads_profile.parquet")["zipstrain_sorted_by"] == "chrom,pos"
    for name in ("reads_gene_stats.parquet", "reads_genome_stats.parquet"):
        assert "zipstrain_sorted_by" not in pl.read_parquet_metadata(run/"sample"/name)


@pytest.mark.parametrize(("method", "engine"), [("table", "table"), ("matrix", "matrix")])
def test_comparison_records_how_it_was_made(run, method, engine):
    metadata = pl.read_parquet_metadata(run/method/"all_comparisons.parquet")
    assert metadata["zipstrain_compare_engine"] == engine
    assert metadata["zipstrain_compare_kind"] == "gene"
    assert metadata["zipstrain_compare_calculate"] == "genome_ani+max_identical_stretch+gene_ani+gene_selection"
    assert metadata["zipstrain_compare_selection_method"] == "sns_snv_reference_opportunities_v2"


def test_table_and_matrix_write_the_same_rows(run):
    from polars.testing import assert_frame_equal

    keys = ["sample_1", "sample_2", "gene_id"]
    assert_frame_equal(pl.read_parquet(run/"table/all_comparisons.parquet").sort(keys),
                       pl.read_parquet(run/"matrix/all_comparisons.parquet").sort(keys))


def test_readers_reject_an_unknown_schema_version(run, tmp_path):
    source = run/"assets/gene_info.parquet"
    metadata = dict(pl.read_parquet_metadata(source)) | {gi.SCHEMA_KEY: "3"}
    future = tmp_path/"gene_info.parquet"
    pl.read_parquet(source).write_parquet(future, metadata=metadata)
    with pytest.raises(ValueError, match="schema"):
        gi.read_gene_info(future)


def _stand_in_slurm(directory):
    """Provide an sbatch that runs the batch script at once and an sacct that reports success."""
    directory.mkdir()
    (directory/"sbatch").write_text('#!/bin/sh\n[ "$1" = "--version" ] && exit 0\nsh "$2" > job.out 2> job.err\necho 4242\n')
    (directory/"sacct").write_text('#!/bin/sh\n[ "$1" = "--version" ] && exit 0\necho COMPLETED\n')
    for name in ("sbatch", "sacct"):
        (directory/name).chmod(0o755)


def test_scheduled_v2_run_matches_local_and_keeps_the_format(run, tmp_path):
    """Slurm tasks run the real pair command. Only the two Slurm programs are stand-ins."""
    import json
    import os
    import subprocess
    import sys

    from polars.testing import assert_frame_equal

    _stand_in_slurm(tmp_path/"bin")
    config = tmp_path/"slurm.json"
    config.write_text(json.dumps({"partition": "x", "time": "1:00:00"}))
    environment = {**os.environ, "PATH": f"{tmp_path/'bin'}{os.pathsep}{os.path.dirname(sys.executable)}{os.pathsep}{os.environ['PATH']}"}
    output = tmp_path/"scheduled"

    def compare(*extra):
        return subprocess.run(
            [sys.executable, "-m", "zipstrain.cli", "compare", "--profile-db", run/"profiles.csv",
             "--preprofile-dir", run/"assets", "--no-csv", "--run-dir", output, "--execution-mode", "slurm",
             "--slurm-config", config, "--poll-interval", "1", *extra],
            cwd=tmp_path, env=environment, capture_output=True, text=True)

    result = compare()
    assert result.returncode == 0, result.stdout + result.stderr
    table = output/"all_comparisons.parquet"
    keys = ["sample_1", "sample_2", "gene_id"]
    assert_frame_equal(pl.read_parquet(table).sort(keys), pl.read_parquet(run/"table/all_comparisons.parquet").sort(keys))
    assert list(pl.read_parquet_schema(table).items()) == list(COMPARISON_COLUMNS.items())
    assert not COMPARISON_KEYS - set(pl.read_parquet_metadata(table))

    before = table.stat().st_mtime_ns
    assert compare().returncode == 0
    assert table.stat().st_mtime_ns == before
    refused = compare("--allele-integration", "weighted")
    assert refused.returncode != 0 and "new --run-dir" in refused.stdout + refused.stderr


def test_scheduled_v2_tasks_mount_every_input(run):
    """Containers mount only file inputs, so assets and codon tables must be among them."""
    from zipstrain import task_manager as tm

    generator = tm.CompareTaskGenerator.__new__(tm.CompareTaskGenerator)
    generator.preprofile_dir, generator.allele_integration = run/"assets", "weighted"
    generator.codon_tables = {"a": run/"sample/a_codon_counts.parquet", "b": run/"sample/b_codon_counts.parquet"}
    inputs = generator._v2_inputs({"sample_name_1": "a", "sample_name_2": "b"})
    mounted = {value.get_value() for value in inputs.values() if isinstance(value, tm.FileInput)}
    assert mounted == {str(run/"assets"), str(run/"sample/a_codon_counts.parquet"), str(run/"sample/b_codon_counts.parquet")}
    arguments = inputs["v2-args"].get_value()
    assert "--allele-integration weighted" in arguments and "--sample-name-2 b" in arguments
