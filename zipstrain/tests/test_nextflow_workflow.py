import os
import pathlib
import re
import shutil
import subprocess
import sys

import polars as pl
import pytest


ROOT_DIR = pathlib.Path(__file__).resolve().parents[2]
NEXTFLOW_FILE = ROOT_DIR / "zipstrain.nf"
CONF_FILE = ROOT_DIR / "conf.config"

# Resolve nextflow + a JVM for it to use. The conda env running these tests
# (`<env>/bin/nextflow` + `<env>/lib/jvm`) is checked *first* and preferred
# when present, so a stray/older `nextflow` earlier on a developer's PATH
# never gets picked up instead of the pinned env copy. Falls back to a plain
# PATH lookup for environments without that conda layout (e.g. CI, where
# actions/setup-java + nextflow-io/setup-nextflow put java/nextflow on PATH
# and export JAVA_HOME directly).
_ENV_PREFIX = pathlib.Path(sys.prefix)
_CONDA_NEXTFLOW_BIN = _ENV_PREFIX / "bin" / "nextflow"
_CONDA_JAVA_HOME = _ENV_PREFIX / "lib" / "jvm"
_CONDA_LAYOUT_AVAILABLE = _CONDA_NEXTFLOW_BIN.exists() and (_CONDA_JAVA_HOME / "bin" / "java").exists()

if _CONDA_LAYOUT_AVAILABLE:
    _NEXTFLOW_BIN = _CONDA_NEXTFLOW_BIN
else:
    _nextflow_on_path = shutil.which("nextflow")
    _NEXTFLOW_BIN = pathlib.Path(_nextflow_on_path) if _nextflow_on_path else _CONDA_NEXTFLOW_BIN

_JAVA_ON_PATH = shutil.which("java") is not None
_JAVA_AVAILABLE = _CONDA_LAYOUT_AVAILABLE or _JAVA_ON_PATH

_NEXTFLOW_AVAILABLE = _NEXTFLOW_BIN.exists() and _JAVA_AVAILABLE


def _nextflow_subprocess_env() -> dict[str, str]:
    """Builds the env for running nextflow: sets JAVA_HOME to the conda env's
    JVM when using the conda layout (its `java` isn't on PATH by itself),
    and makes sure the resolved nextflow binary's directory is on PATH."""
    env = dict(os.environ)
    if _CONDA_LAYOUT_AVAILABLE:
        env["JAVA_HOME"] = str(_CONDA_JAVA_HOME)
    nextflow_dir = str(_NEXTFLOW_BIN.parent)
    env["PATH"] = f"{nextflow_dir}{os.pathsep}{env.get('PATH', '')}"
    return env


def _docker_available() -> bool:
    docker_bin = shutil.which("docker")
    if docker_bin is None:
        return False
    try:
        result = subprocess.run(
            [docker_bin, "info"], capture_output=True, timeout=10
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


_DOCKER_AVAILABLE = _docker_available()


def _extract_output_section(process_name: str) -> str:
    text = NEXTFLOW_FILE.read_text()
    match = re.search(
        rf"process\s+{re.escape(process_name)}\s*\{{.*?output:\s*(.*?)\s*script:",
        text,
        flags=re.DOTALL,
    )
    assert match is not None, f"process {process_name} not found in zipstrain.nf"
    return match.group(1)


def _nextflow_process_names() -> set[str]:
    return set(re.findall(r"^process\s+([A-Za-z0-9_]+)\s*\{", NEXTFLOW_FILE.read_text(), flags=re.MULTILINE))


def _conf_process_names() -> set[str]:
    return set(re.findall(r"withName:\s*'([^']+)'", CONF_FILE.read_text()))


def test_profile_bam_process_emits_gene_stats():
    output_section = _extract_output_section("profile_bam")
    assert 'path "${bamfile.baseName}_gene_stats.parquet", emit: gene_stats' in output_section


def test_from_sra_to_profile_process_emits_gene_stats():
    output_section = _extract_output_section("fromSRAtoProfile")
    assert 'path "${sra_id}_gene_stats.parquet", emit: gene_stats' in output_section


def test_from_sra_to_profile_build_db_process_emits_gene_stats_and_sylph_output():
    output_section = _extract_output_section("fromSRAtoProfileBuildDb")
    assert 'path "${sra_id}_gene_stats.parquet", emit: gene_stats' in output_section
    assert 'path "${sra_id}_sylph_abundance.tsv", emit: sylph_abundance' in output_section


def test_nextflow_calls_utilities_single_compare_commands():
    text = NEXTFLOW_FILE.read_text()
    assert "zipstrain utilities single_compare_genome" in text
    assert "zipstrain utilities single_compare_gene" in text


def test_nextflow_calls_updated_profile_commands():
    text = NEXTFLOW_FILE.read_text()
    assert "zipstrain utilities prepare_profiling" in text
    assert "zipstrain utilities profile-single" in text


def test_nextflow_compare_profile_tables_use_profile_location_names_and_engine_param():
    text = NEXTFLOW_FILE.read_text()
    assert 'params.compare_engine="polars"' in text
    assert "getProfileLocationsTableColumn" in text
    assert "getProfileSampleNamesTableColumn" in text
    assert "profile_location" in text
    assert "--engine ${params.compare_engine}" in text


def test_nextflow_from_sra_to_profile_auto_builds_reference_without_genes():
    text = NEXTFLOW_FILE.read_text()
    assert 'if (!params.reference_genome)' in text
    assert "fromSRAtoProfileBuildDb(sra_ids, sylph_db)" in text
    assert "zipstrain utilities build-genome-db" in text
    assert "--reference-fasta reference_genomes.fna" in text
    assert "--stb-file reference_genomes.stb" in text


def test_nextflow_has_no_gene_prepare_profile_path():
    text = NEXTFLOW_FILE.read_text()
    assert "process prepare_profile_no_genes{" in text
    assert "prepare_profile_no_genes(reference_genome, file(params.stb))" in text


def test_profile_mode_uses_profile_workflow_and_supports_no_gene_path():
    text = NEXTFLOW_FILE.read_text()
    assert "workflow profile{" in text
    assert "fast_profile(" not in text
    assert "profile(bamfiles, sample_names, gene_file, reference_genome)" in text
    assert "gene_file = params.gene_file ? file(params.gene_file) : null" in text


def test_prepare_profile_process_emits_contract_and_null_model():
    output_section = _extract_output_section("prepare_profile")
    assert 'path "null_model.parquet", emit: null_model' in output_section
    assert 'path "profiling_contract.json", emit: profiling_contract' in output_section


def test_prepare_profile_no_genes_process_emits_contract_and_null_model():
    output_section = _extract_output_section("prepare_profile_no_genes")
    assert 'path "null_model.parquet", emit: null_model' in output_section
    assert 'path "profiling_contract.json", emit: profiling_contract' in output_section


def test_nextflow_profile_processes_pass_profiling_contract():
    text = NEXTFLOW_FILE.read_text()
    assert "--reference-fasta ${reference_fasta}" in text
    assert "--reference-fasta ${reference_genome}" in text
    assert "--reference-fasta reference_genomes.fna" in text
    assert "--profiling-contract ${profiling_contract}" in text
    assert "--profiling-contract profiling_contract.json" in text
    assert "prepare_profile.out.profiling_contract" in text
    assert "prepare_profile_no_genes.out.profiling_contract" in text


def test_nextflow_processes_match_conf_config():
    assert _nextflow_process_names() == _conf_process_names()


def _write_compare_genomes_fixtures(tmp_path: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    """Writes two tiny, identical profile parquets + an stb file and returns
    (profiles_csv, stb_path) for a `mode=compare_genomes input_type=profile_table` run."""
    profile_schema = {
        "chrom": ["chr1", "chr1", "chr1"],
        "genome": ["genome1", "genome1", "genome1"],
        "pos": [0, 1, 2],
        "gene": ["NA", "NA", "NA"],
        "A": [5, 0, 5],
        "T": [0, 5, 0],
        "C": [0, 0, 0],
        "G": [0, 0, 0],
    }
    profile_1 = tmp_path / "sample1_profile.parquet"
    profile_2 = tmp_path / "sample2_profile.parquet"
    pl.DataFrame(profile_schema).write_parquet(profile_1)
    pl.DataFrame(profile_schema).write_parquet(profile_2)

    stb_path = tmp_path / "reference.stb"
    pl.DataFrame({"scaffold": ["chr1"], "genome": ["genome1"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )

    profiles_csv = tmp_path / "profiles.csv"
    profiles_csv.write_text(
        "sample_name,profile_location\n"
        f"sample1,{profile_1}\n"
        f"sample2,{profile_2}\n"
    )
    return profiles_csv, stb_path


def _assert_merged_comparison_is_correct(output_dir: pathlib.Path, nextflow_log: str) -> None:
    merged_path = output_dir / "merged_comparisons.parquet"
    assert merged_path.exists(), f"missing output; nextflow log:\n{nextflow_log}"

    merged = pl.read_parquet(merged_path)
    assert merged.height == 1
    row = merged.row(0, named=True)
    assert row["genome"] == "genome1"
    assert row["genome_pop_ani"] == 100.0
    assert {row["sample_1"], row["sample_2"]} == {"sample1", "sample2"}


@pytest.mark.skipif(
    not _NEXTFLOW_AVAILABLE,
    reason="nextflow + a JVM must be installed in the current Python environment "
    "(e.g. `conda install -c bioconda nextflow openjdk`) to run the pipeline for real.",
)
def test_compare_genomes_mode_runs_end_to_end(tmp_path):
    """
    Actually compiles and executes zipstrain.nf, rather than grepping its text.

    This uses `mode=compare_genomes` with `input_type=profile_table` because it is
    the only mode that needs no containers, no samtools/bowtie2/sylph, and no
    network access: it just calls the `zipstrain` CLI (already on PATH in this
    environment) on two tiny profile parquets and merges the result.
    """
    profiles_csv, stb_path = _write_compare_genomes_fixtures(tmp_path)
    output_dir = tmp_path / "out"

    # The repo's nextflow.config enables Docker by default and includes conf.config,
    # which sizes these processes at 8 cpus for cluster nodes. This test is meant to
    # run container-free with the local `zipstrain` CLI on PATH, so disable Docker and
    # shrink the cpu/memory requests to something a 4-core CI runner can satisfy.
    override = tmp_path / "test_local_resources.config"
    override.write_text(
        "docker { enabled = false }\n"
        "process {\n"
        "    withName: 'compare_genome_fast_profiles_single' { cpus = 1; memory = '512 MB' }\n"
        "    withName: 'merge_comparison_tables' { cpus = 1; memory = '512 MB' }\n"
        "}\n"
    )

    result = subprocess.run(
        [
            str(_NEXTFLOW_BIN),
            "run",
            str(NEXTFLOW_FILE),
            "-c", str(override),
            "--mode", "compare_genomes",
            "--input_type", "profile_table",
            "--input_table", str(profiles_csv),
            "--stb", str(stb_path),
            "--parallel_mode", "single",
            "--output_dir", str(output_dir),
        ],
        cwd=tmp_path,
        env=_nextflow_subprocess_env(),
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, (
        f"nextflow run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    _assert_merged_comparison_is_correct(output_dir, result.stdout)


@pytest.mark.skipif(
    not _NEXTFLOW_AVAILABLE,
    reason="nextflow + a JVM must be installed in the current Python environment "
    "(e.g. `conda install -c bioconda nextflow openjdk`) to run the pipeline for real.",
)
@pytest.mark.skipif(
    not _DOCKER_AVAILABLE,
    reason="a running Docker daemon is required to exercise the `-profile docker` path.",
)
def test_compare_genomes_mode_runs_end_to_end_via_docker(tmp_path):
    """
    Same scenario as test_compare_genomes_mode_runs_end_to_end, but run under
    `-profile docker` so every process (including compare_genome_fast_profiles_single)
    executes inside the published parsaghadermazi/zipstrain image. This is the
    cheapest way to catch a broken/stale published image (missing tools, wrong
    zipstrain version, wrong CPU architecture) without the slow/flaky SRA + Sylph
    download path.
    """
    profiles_csv, stb_path = _write_compare_genomes_fixtures(tmp_path)
    output_dir = tmp_path / "out"

    # conf.config requests 8 cpus / 40 GB for compare_genome_fast_profiles_single,
    # sized for cluster nodes. Override to something any laptop's Docker Desktop
    # allocation can satisfy, since this test only cares whether the container
    # itself works, not real-world throughput.
    resource_override = tmp_path / "test_docker_resources.config"
    resource_override.write_text(
        "process {\n"
        "    withName: 'compare_genome_fast_profiles_single' { cpus = 1; memory = '512 MB' }\n"
        "    withName: 'merge_comparison_tables' { cpus = 1; memory = '512 MB' }\n"
        "}\n"
    )

    result = subprocess.run(
        [
            str(_NEXTFLOW_BIN),
            "run",
            str(NEXTFLOW_FILE),
            "-c", str(CONF_FILE),
            "-c", str(resource_override),
            "-profile", "docker",
            "--mode", "compare_genomes",
            "--input_type", "profile_table",
            "--input_table", str(profiles_csv),
            "--stb", str(stb_path),
            "--parallel_mode", "single",
            "--output_dir", str(output_dir),
        ],
        cwd=tmp_path,
        env=_nextflow_subprocess_env(),
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"nextflow run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    _assert_merged_comparison_is_correct(output_dir, result.stdout)
