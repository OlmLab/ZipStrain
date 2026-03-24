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


def test_compare_genome_process_supports_calculate_and_ani_method():
    text = NEXTFLOW_FILE.read_text()
    assert 'params.compare_genome_calculate="all"' in text
    assert 'params.compare_ani_method="popani"' in text
    assert "--ani-method ${params.compare_ani_method}" in text
    assert "--calculate ${params.compare_genome_calculate}" in text
