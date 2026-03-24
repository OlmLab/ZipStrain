import pathlib
import re


ROOT_DIR = pathlib.Path(__file__).resolve().parents[2]
NEXTFLOW_FILE = ROOT_DIR / "zipstrain.nf"


def _extract_output_section(process_name: str) -> str:
    text = NEXTFLOW_FILE.read_text()
    match = re.search(
        rf"process\s+{re.escape(process_name)}\s*\{{.*?output:\s*(.*?)\s*script:",
        text,
        flags=re.DOTALL,
    )
    assert match is not None, f"process {process_name} not found in zipstrain.nf"
    return match.group(1)


def test_profile_bam_process_emits_gene_stats():
    output_section = _extract_output_section("profile_bam")
    assert 'path "${bamfile.baseName}_gene_stats.parquet", emit: gene_stats' in output_section


def test_from_sra_to_profile_process_emits_gene_stats():
    output_section = _extract_output_section("fromSRAtoProfile")
    assert 'path "${sra_id}_gene_stats.parquet", emit: gene_stats' in output_section


def test_nextflow_calls_utilities_single_compare_commands():
    text = NEXTFLOW_FILE.read_text()
    assert "zipstrain utilities single_compare_genome" in text
    assert "zipstrain utilities single_compare_gene" in text


def test_nextflow_calls_updated_profile_commands():
    text = NEXTFLOW_FILE.read_text()
    assert "zipstrain utilities prepare_profiling" in text
    assert "zipstrain utilities profile-single" in text
