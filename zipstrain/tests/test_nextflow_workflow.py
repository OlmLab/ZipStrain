import pathlib
import re


ROOT_DIR = pathlib.Path(__file__).resolve().parents[2]
NEXTFLOW_FILE = ROOT_DIR / "zipstrain.nf"
CONF_FILE = ROOT_DIR / "conf.config"


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


def test_profile_bam_no_reference_process_emits_gene_stats():
    output_section = _extract_output_section("profile_bam_no_reference")
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


def test_profile_mode_supports_precomputed_assets_or_reference_prep():
    text = NEXTFLOW_FILE.read_text()
    assert "workflow profile{" in text
    assert "fast_profile(" not in text
    assert "profile_uses_prepared_assets = params.bed_file || params.gene_range_table || params.null_model || params.profiling_contract" in text
    assert 'required_profile_asset_params = ["stb", "bed_file", "gene_range_table", "null_model", "profiling_contract"]' in text
    assert "mode=profile requires --reference_genome unless you provide precomputed profiling assets." in text
    assert "profile_bam_no_reference(sample_names, bamfiles, bed_file, stb_file, gene_range_table, null_model, profiling_contract)" in text


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
    assert "--reference-fasta ${reference_genome}" in text
    assert "--reference-fasta reference_genomes.fna" in text
    assert "--profiling-contract ${profiling_contract}" in text
    assert "--profiling-contract profiling_contract.json" in text
    assert "--min-mapq ${params.min_mapq}" in text
    assert "--min-baseq ${params.min_baseq}" in text
    assert "--read-inclusion ${params.read_inclusion}" in text
    assert "getOptionalReadAniArg" in text
    assert "prepare_profile.out.profiling_contract" in text
    assert "prepare_profile_no_genes.out.profiling_contract" in text


def test_nextflow_processes_match_conf_config():
    assert _nextflow_process_names() == _conf_process_names()
