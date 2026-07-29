from zipstrain import profile
from zipstrain import cli as cli_module
from zipstrain import utils
import polars as pl
import pytest
from click.testing import CliRunner
from pathlib import Path
import re


a_chr1=     [ 1, 0, 1, 0, 0, 0, 4, 4, 2, 1]
t_chr1=     [ 0, 2, 3, 1, 0, 0, 1, 0, 0, 1]
c_chr1=     [10, 1, 2, 0, 0, 0, 1, 0, 0, 0]
g_chr1=     [ 0, 2, 0, 2, 3, 1, 0, 2, 7, 0]
a_chr2=     [ 0, 0, 4, 3, 2, 4, 0, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0]
a_chr2_1mut=[ 0, 0, 4, 3, 2, 4, 8, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0] # A single mutation at index (0->1)
t_chr2=     [ 1, 0, 0, 1, 0, 0, 7, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2]
t_chr2_1mut=[ 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2] # A single mutation at index (19->20)
c_chr2=     [ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
c_chr2_1mut=[ 0, 0, 1, 2, 1, 1, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0] # A single mutation at index (18->19)
g_chr2=     [ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
g_chr2_1mut=[ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
a_chr3=     [ 0, 0, 1, 0, 0, 0,15, 0, 0, 0, 0, 0, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_chr3=     [ 1, 0, 0, 1, 0,17, 1, 0,10,11, 0, 2, 3, 1, 3, 2,12,10, 1, 6, 1, 0, 0, 1, 0, 0, 8, 5, 0, 1]
c_chr3=     [ 0, 0,19,21,18, 1, 2,10, 0, 0,12, 1, 0, 0, 0,13, 0, 0,11, 0, 1, 0, 1, 2, 1,20, 0, 1, 0, 0]
g_chr3=     [20,20, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 2, 3, 1, 1, 2, 7, 0]

gene_locs = (["NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA","NA"])+ \
            (["NA","NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA",
              "NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA","NA"])+ \
            (["NA","NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA",
              "NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA","NA",
              "NA","NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA"])


@pytest.fixture(scope="module")
def profile_1()->pl.LazyFrame:
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()


def test_null_model_defaults_and_poisson_boundary():
    assert utils.NULL_MODEL_ERROR_RATE_DEFAULT == 0.001
    assert utils.NULL_MODEL_MAX_COVERAGE_DEFAULT == 50_000
    records = utils.build_null_poisson(max_total_reads=1_000)
    assert len(records) == 1_000
    # At 0.1% total error and 1,000x depth, one observation is still
    # compatible with error; an allele must have at least two to survive.
    assert records[-1] == (1_000, 1)


def test_adjust_for_sequence_errors_uses_strict_ceiling_and_min_freq():
    raw = pl.DataFrame(
        {"A": [90], "C": [5], "G": [4], "T": [1]}
    ).lazy()
    null_model = pl.DataFrame(
        {"cov": [100], "max_error_count": [4]}
    ).lazy()

    adjusted = profile.adjust_for_sequence_errors(
        raw,
        null_model,
        min_freq=0.05,
    ).collect()

    # C is exactly at the frequency cutoff and survives; G is exactly at the
    # null-model ceiling and is removed.
    assert adjusted.select("A", "C", "G", "T").row(0) == (90, 5, 0, 0)


def test_adjust_for_sequence_errors_rejects_coverage_above_model_limit():
    raw = pl.DataFrame(
        {"A": [101], "C": [0], "G": [0], "T": [0]}
    ).lazy()
    null_model = pl.DataFrame(
        {"cov": [100], "max_error_count": [4]}
    ).lazy()

    with pytest.raises(ValueError, match="Observed coverage 101 exceeds"):
        profile.adjust_for_sequence_errors(raw, null_model)


@pytest.mark.parametrize("min_freq", [-0.01, 1.01])
def test_adjust_for_sequence_errors_rejects_invalid_min_freq(min_freq):
    raw = pl.DataFrame(
        {"A": [1], "C": [0], "G": [0], "T": [0]}
    ).lazy()
    null_model = pl.DataFrame(
        {"cov": [1], "max_error_count": [0]}
    ).lazy()

    with pytest.raises(ValueError, match="min_freq must be between 0 and 1"):
        profile.adjust_for_sequence_errors(raw, null_model, min_freq=min_freq)
    

def test_duckdb_chunk_annotation_matches_polars_for_unsorted_chunk(tmp_path: Path):
    adjusted = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "pos": [1, 2, 5, 2, 4],
            "A": [10, 0, 1, 0, 2],
            "C": [0, 8, 0, 0, 0],
            "G": [0, 0, 0, 9, 0],
            "T": [0, 0, 0, 0, 1],
        }
    )
    adjusted_pf = tmp_path / "adj.parquet"
    out_pf = tmp_path / "ann.parquet"
    adjusted.sort(["chrom", "pos"], descending=[True, True]).write_parquet(adjusted_pf)

    stb_lf = pl.DataFrame(
        {
            "scaffold": ["chr1", "chr2"],
            "genome": ["genome1", "genome2"],
        }
    ).lazy()
    gene_lf = pl.DataFrame(
        {
            "gene": ["geneA", "geneB"],
            "scaffold": ["chr1", "chr2"],
            "start": [1, 2],
            "end": [3, 2],
        }
    ).lazy()

    profile._annotate_mpileup_chunk_with_duckdb(
        adjusted_mpileup_parquet=adjusted_pf,
        output_parquet=out_pf,
        scaffold_to_genome=stb_lf,
        gene_range=gene_lf,
    )

    got = pl.read_parquet(out_pf).sort(["genome", "chrom", "pos"])
    expected = (
        profile.add_gene_info_to_mpileup(
            mpileup_df=profile.add_genome_info_to_mpileup(
                mpileup_df=pl.scan_parquet(adjusted_pf),
                scaffold_to_genome=stb_lf,
            ),
            gene_range=gene_lf,
        )
        .drop(["start", "end"])
        .select(["chrom", "genome", "gene", "pos", "A", "C", "G", "T"])
        .collect()
        .sort(["genome", "chrom", "pos"])
    )
    assert got.to_dicts() == expected.to_dicts()


def test_write_sorted_profile_with_metadata_skips_sort_for_sorted_input(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    output_file = tmp_path / "sorted_profile.parquet"
    candidate = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "genome": ["genome1", "genome1", "genome2"],
            "gene": ["geneA", "geneA", "geneB"],
            "pos": [1, 2, 1],
            "A": [1, 0, 0],
            "C": [0, 2, 0],
            "G": [0, 0, 3],
            "T": [0, 0, 0],
        }
    ).lazy()

    def _unexpected_sort(*args, **kwargs):
        raise AssertionError("sorted input should not need fallback sorting")

    monkeypatch.setattr(profile, "_sort_existing_profile_parquet", _unexpected_sort)

    profile._write_sorted_profile_with_metadata(
        profile_lf=candidate,
        output_file=output_file,
        tmp_dir=tmp_path,
    )

    got = pl.read_parquet(output_file)
    assert got.to_dicts() == got.sort(["chrom", "pos"]).to_dicts()
    assert (
        pl.read_parquet_metadata(output_file)[profile.PROFILE_SORTED_METADATA_KEY]
        == profile.PROFILE_SORTED_METADATA_VALUE
    )


def test_write_sorted_profile_with_metadata_falls_back_to_sort_for_unsorted_input(tmp_path: Path):
    output_file = tmp_path / "unsorted_profile.parquet"
    candidate = pl.DataFrame(
        {
            "chrom": ["chr2", "chr1", "chr1"],
            "genome": ["genome2", "genome1", "genome1"],
            "gene": ["geneB", "geneA", "geneA"],
            "pos": [1, 2, 1],
            "A": [0, 0, 1],
            "C": [0, 2, 0],
            "G": [3, 0, 0],
            "T": [0, 0, 0],
        }
    ).lazy()

    profile._write_sorted_profile_with_metadata(
        profile_lf=candidate,
        output_file=output_file,
        tmp_dir=tmp_path,
    )

    got = pl.read_parquet(output_file)
    assert got.to_dicts() == got.sort(["chrom", "pos"]).to_dicts()
    assert (
        pl.read_parquet_metadata(output_file)[profile.PROFILE_SORTED_METADATA_KEY]
        == profile.PROFILE_SORTED_METADATA_VALUE
    )


def test_profile_filter_helpers_handle_read_inclusion_and_read_ani():
    assert profile._aligned_query_bases_from_cigar("10S90M5I2D") == 95

    assert profile._sam_alignment_passes_profile_filters(
        flag=3,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PROPER_PAIRS,
    )
    assert not profile._sam_alignment_passes_profile_filters(
        flag=3,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:4"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PROPER_PAIRS,
    )
    assert profile._sam_alignment_passes_profile_filters(
        flag=1,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
    )
    # 'paired' is single-end-safe: a genuinely single-end read (not flagged
    # paired) is kept, since it cannot be part of a discordant pair.
    assert profile._sam_alignment_passes_profile_filters(
        flag=0,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
    )
    # ...but a half-mapped orphan (paired flag + mate unmapped) is dropped.
    assert not profile._sam_alignment_passes_profile_filters(
        flag=profile.SAM_FLAG_PAIRED | profile.SAM_FLAG_MUNMAP,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
    )
    # 'paired' matches inStrain 'paired_only': a paired read whose mate maps to the
    # SAME scaffold (RNEXT == "=") is kept...
    assert profile._sam_alignment_passes_profile_filters(
        flag=profile.SAM_FLAG_PAIRED,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
        rname="chr1",
        rnext="=",
    )
    # ...while a paired read whose mate maps to a DIFFERENT scaffold is dropped
    # (these cross-scaffold pairs otherwise inflate coverage relative to inStrain).
    assert not profile._sam_alignment_passes_profile_filters(
        flag=profile.SAM_FLAG_PAIRED,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
        rname="chr1",
        rnext="chr2",
    )
    assert profile._sam_alignment_passes_profile_filters(
        flag=0,
        mapq=42,
        cigar="100M",
        optional_fields=["NM:i:2"],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_ALL_MAPPED,
    )
    # A read without an NM tag is kept rather than failing the run when a
    # read-ANI floor is set (mappers that omit NM just don't get ANI filtering).
    assert profile._sam_alignment_passes_profile_filters(
        flag=0,
        mapq=42,
        cigar="100M",
        optional_fields=[],
        min_mapq=10,
        min_read_ani=0.97,
        read_inclusion=profile.READ_INCLUSION_ALL_MAPPED,
    )


def test_read_stb_strips_whitespace(tmp_path):
    # An STB with stray spaces around the tab (which inStrain tolerates) must
    # still yield clean scaffold/genome names so scaffolds are not dropped.
    stb = tmp_path / "padded.stb"
    stb.write_text("scaffold_0 \t genomeA\nscaffold_1\tgenomeA \n")
    df = profile.read_stb(stb).collect()
    assert df.columns == ["scaffold", "genome"]
    assert df["scaffold"].to_list() == ["scaffold_0", "scaffold_1"]
    assert df["genome"].to_list() == ["genomeA", "genomeA"]


def test_profile_defaults_match_instrain_leaning_filters():
    assert profile.PROFILE_MIN_READ_ANI_DEFAULT == 0.95
    assert profile.PROFILE_READ_INCLUSION_DEFAULT == profile.READ_INCLUSION_PAIRED


def _write_profile_test_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, Path, pl.LazyFrame]:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACG\n>chr2\nTAC\n")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t3\nchr2\t0\t3\n")
    gene_range_table = tmp_path / "gene_ranges.tsv"
    gene_range_table.write_text("geneA\tchr1\t1\t3\ngeneB\tchr2\t1\t3\n")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\nchr2\tgenome2\n")
    null_model_file = tmp_path / "null_model.parquet"
    pl.DataFrame(
        {
            "cov": [0, 1, 5, 6, 7, 8, 9, 10],
            "max_error_count": [0, 0, 0, 0, 0, 0, 0, 0],
        }
    ).write_parquet(null_model_file)
    stb_lf = pl.scan_csv(stb_file, separator="\t", has_header=False).with_columns(
        pl.col("column_1").alias("scaffold"),
        pl.col("column_2").alias("genome"),
    )
    return bam_file, reference_fasta, bed_file, gene_range_table, stb_file, null_model_file, stb_lf


def _install_fake_profile_subprocess(monkeypatch: pytest.MonkeyPatch, observed_commands: list[str] | None = None):
    observed_mpileups: list[str] = []
    observed_reference_flags: list[bool] = []

    class _FakeStdout:
        """Iterable/closeable stand-in for Popen.stdout (yields byte lines)."""
        def __init__(self, lines: list[bytes] | None = None):
            self._lines = list(lines or [])

        def __iter__(self):
            return iter(self._lines)

        def close(self):
            pass

    class _FakePopen:
        def __init__(self, *, stdout_lines: list[bytes] | None = None, returncode: int = 0):
            self.returncode = returncode
            self.stdout = _FakeStdout(stdout_lines)
            self.pid = -1

        def wait(self, timeout=None):
            return self.returncode

        def poll(self):
            return self.returncode

    class _FakeCompleted:
        def __init__(self, *, returncode: int = 0, stdout: bytes = b"", stderr: bytes = b""):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def _read_scaffolds_from_bed(bed_path: Path) -> list[str]:
        return [line.strip().split("\t")[0] for line in bed_path.read_text().splitlines() if line.strip()]

    def _fake_mpileup_rows(scaffold: str) -> list[dict]:
        if scaffold == "chr1":
            return [
                {"chrom": "chr1", "pos": 1, "ref": "A", "A": 10, "C": 0, "G": 0, "T": 0, "bases": "." * 10},
                {"chrom": "chr1", "pos": 2, "ref": "C", "A": 0, "C": 6, "G": 0, "T": 0, "bases": "." * 6},
                {"chrom": "chr1", "pos": 3, "ref": "G", "A": 0, "C": 0, "G": 8, "T": 0, "bases": "." * 8},
            ]
        if scaffold == "chr2":
            return [
                {"chrom": "chr2", "pos": 1, "ref": "T", "A": 0, "C": 0, "G": 0, "T": 7, "bases": "." * 7},
                {"chrom": "chr2", "pos": 2, "ref": "A", "A": 5, "C": 0, "G": 0, "T": 0, "bases": "." * 5},
                {"chrom": "chr2", "pos": 3, "ref": "C", "A": 0, "C": 9, "G": 0, "T": 0, "bases": "." * 9},
            ]
        raise AssertionError(f"Unexpected scaffold in fake mpileup data: {scaffold}")

    def _fake_read_loc_rows(scaffold: str) -> list[dict]:
        return [
            {"chrom": scaffold, "pos": 1},
            {"chrom": scaffold, "pos": 2},
            {"chrom": scaffold, "pos": 3},
        ]

    def _fake_mpileup_line(row: dict, *, use_reference_matches: bool) -> bytes:
        bases = row["bases"] if use_reference_matches else row["ref"] * len(row["bases"])
        depth = len(bases)
        return f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{depth}\t{bases}\t*\n".encode()

    def _fake_popen(cmd, stdout=None, stderr=None, **kwargs):
        # Profiling launches `samtools mpileup` via subprocess.Popen with a list argv.
        cmd_list = list(cmd) if isinstance(cmd, (list, tuple)) else re.split(r"\s+", cmd)
        command = " ".join(str(c) for c in cmd_list)
        if observed_commands is not None:
            observed_commands.append(command)
        if cmd_list[:2] == ["samtools", "mpileup"]:
            observed_mpileups.append(command)
            has_reference_fasta = "-f" in cmd_list
            observed_reference_flags.append(has_reference_fasta)
            bed_path = cmd_list[cmd_list.index("-l") + 1]
            scaffolds = _read_scaffolds_from_bed(Path(bed_path))
            lines = []
            for scaffold in scaffolds:
                lines.extend(
                    _fake_mpileup_line(row, use_reference_matches=has_reference_fasta)
                    for row in _fake_mpileup_rows(scaffold)
                )
            return _FakePopen(stdout_lines=lines)
        raise AssertionError(f"Unexpected Popen command in fake subprocess: {command}")

    def _fake_run(cmd, shell=False, stdout=None, stderr=None, **kwargs):
        command = cmd if isinstance(cmd, str) else " ".join(str(c) for c in cmd)
        if observed_commands is not None:
            observed_commands.append(command)
        if "process-read-locs" in command:
            # samtools view | zipstrain utilities process-read-locs --output-file <abs>
            out_match = re.search(r"--output-file\s+([^\s]+)", command)
            assert out_match is not None, command
            output_file = Path(out_match.group(1))
            output_file.parent.mkdir(parents=True, exist_ok=True)
            bed_match = re.search(r"-L\s+([^\s]+)", command)
            assert bed_match is not None, command
            scaffolds = _read_scaffolds_from_bed(Path(bed_match.group(1)))
            rows = []
            for scaffold in scaffolds:
                rows.extend(_fake_read_loc_rows(scaffold))
            pl.DataFrame(rows).write_parquet(output_file)
            return _FakeCompleted()
        if isinstance(cmd, (list, tuple)) and list(cmd[:2]) == ["samtools", "faidx"]:
            return _FakeCompleted()
        raise AssertionError(f"Unexpected run command in fake subprocess: {command}")

    monkeypatch.setattr(profile.subprocess, "Popen", _fake_popen)
    monkeypatch.setattr(profile.subprocess, "run", _fake_run)

    # The default profiling filters now trigger a samtools prefilter pass; these
    # fake-subprocess tests don't exercise real samtools, so no-op it. Tests that
    # want to observe the prefilter override this setattr afterwards.
    def _noop_filter_bam_for_profiling(**kwargs):
        kwargs["output_bam"].write_bytes(kwargs["input_bam"].read_bytes())
        bai = kwargs["output_bam"].with_suffix(kwargs["output_bam"].suffix + ".bai")
        bai.write_text("")

    monkeypatch.setattr(profile, "_filter_bam_for_profiling", _noop_filter_bam_for_profiling)
    return observed_mpileups, observed_reference_flags


def test_profile_bam_end_to_end_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    observed_limits, observed_reference_flags = _install_fake_profile_subprocess(monkeypatch)
    bam_file, reference_fasta, bed_file, gene_range_table, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)
    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
        reference_fasta=str(reference_fasta),
        gene_range_table=str(gene_range_table),
        stb=stb_lf,
        null_model=pl.scan_parquet(null_model_file),
        output_dir=str(tmp_path),
        num_chunks=2,
        max_concurrency=2,
    )

    profile_file = tmp_path / "sample_profile.parquet"
    gene_stats_file = tmp_path / "sample_gene_stats.parquet"
    genome_stats_file = tmp_path / "sample_genome_stats.parquet"
    assert profile_file.exists()
    assert gene_stats_file.exists()
    assert genome_stats_file.exists()
    assert not (tmp_path / "tmp").exists()

    prof = pl.read_parquet(profile_file)
    assert prof.columns == ["chrom", "genome", "gene", "pos", "A", "C", "G", "T", "ref_base_bitmask"]
    assert prof.height == 6
    assert prof.schema["chrom"] == pl.Utf8
    assert prof.schema["genome"] == pl.Utf8
    assert observed_limits
    assert observed_reference_flags
    assert all(observed_reference_flags)
    assert prof.schema["gene"] == pl.Utf8
    assert prof.schema["ref_base_bitmask"] == pl.UInt8
    assert prof.to_dicts() == prof.sort(["chrom", "pos"]).to_dicts()
    assert prof["ref_base_bitmask"].to_list() == [1, 2, 4, 8, 1, 2]
    assert (
        pl.read_parquet_metadata(profile_file)[profile.PROFILE_SORTED_METADATA_KEY]
        == profile.PROFILE_SORTED_METADATA_VALUE
    )

    gene_stats = pl.read_parquet(gene_stats_file)
    g1 = gene_stats.filter((pl.col("genome") == "genome1") & (pl.col("gene") == "geneA")).to_dicts()[0]
    g2 = gene_stats.filter((pl.col("genome") == "genome2") & (pl.col("gene") == "geneB")).to_dicts()[0]
    assert g1["length"] == 3
    assert g1["breadth"] == pytest.approx(1.0)
    assert g1["coverage"] == pytest.approx(8.0)
    assert g1["ref_ani"] == pytest.approx(100.0)
    assert g2["length"] == 3
    assert g2["breadth"] == pytest.approx(1.0)
    assert g2["coverage"] == pytest.approx(7.0)
    assert g2["ref_ani"] == pytest.approx(100.0)

    genome_stats = pl.read_parquet(genome_stats_file)
    assert set(genome_stats["genome"].to_list()) == {"genome1", "genome2"}
    by_genome = genome_stats.rows_by_key("genome", unique=True, named=True)
    assert by_genome["genome1"]["coverage"] == pytest.approx(8.0)
    assert by_genome["genome2"]["coverage"] == pytest.approx(7.0)
    assert by_genome["genome1"]["breadth"] == pytest.approx(1.0)
    assert by_genome["genome2"]["breadth"] == pytest.approx(1.0)
    assert by_genome["genome1"]["ref_ani"] == pytest.approx(100.0)
    assert by_genome["genome2"]["ref_ani"] == pytest.approx(100.0)


def test_profile_bam_prefilters_alignments_and_applies_mpileup_thresholds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    observed_commands: list[str] = []
    observed_limits, _ = _install_fake_profile_subprocess(monkeypatch, observed_commands=observed_commands)
    bam_file, reference_fasta, bed_file, gene_range_table, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)

    captured_prefilter = {}

    def _fake_filter_bam_for_profiling(**kwargs):
        captured_prefilter.update(kwargs)
        kwargs["output_bam"].write_text("")
        kwargs["output_bam"].with_suffix(".bam.bai").write_text("")

    monkeypatch.setattr(profile, "_filter_bam_for_profiling", _fake_filter_bam_for_profiling)

    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
        reference_fasta=str(reference_fasta),
        gene_range_table=str(gene_range_table),
        stb=stb_lf,
        null_model=pl.scan_parquet(null_model_file),
        output_dir=str(tmp_path),
        num_chunks=2,
        max_concurrency=2,
        min_mapq=7,
        min_baseq=23,
        min_read_ani=0.96,
        read_inclusion=profile.READ_INCLUSION_PAIRED,
    )

    assert captured_prefilter["min_mapq"] == 7
    assert captured_prefilter["min_read_ani"] == 0.96
    assert captured_prefilter["read_inclusion"] == profile.READ_INCLUSION_PAIRED
    assert captured_prefilter["output_bam"].name.endswith(".filtered.bam")
    mpileup_commands = [cmd for cmd in observed_commands if cmd.startswith("samtools mpileup")]
    assert mpileup_commands
    assert all("-q 7" in cmd for cmd in mpileup_commands)
    assert all("-Q 23" in cmd for cmd in mpileup_commands)
    assert all(".filtered.bam" in cmd for cmd in mpileup_commands)
    assert observed_limits


def test_profile_bam_end_to_end_outputs_without_gene_ranges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, reference_fasta, bed_file, _, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)
    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
        reference_fasta=str(reference_fasta),
        gene_range_table=None,
        stb=stb_lf,
        null_model=pl.scan_parquet(null_model_file),
        output_dir=str(tmp_path),
        num_chunks=2,
        max_concurrency=2,
    )

    profile_file = tmp_path / "sample_profile.parquet"
    gene_stats_file = tmp_path / "sample_gene_stats.parquet"
    genome_stats_file = tmp_path / "sample_genome_stats.parquet"
    assert profile_file.exists()
    assert gene_stats_file.exists()
    assert genome_stats_file.exists()

    prof = pl.read_parquet(profile_file)
    assert set(prof["gene"].unique().to_list()) == {"NA"}
    assert "ref_base_bitmask" in prof.columns
    assert (
        pl.read_parquet_metadata(profile_file)[profile.PROFILE_SORTED_METADATA_KEY]
        == profile.PROFILE_SORTED_METADATA_VALUE
    )

    gene_stats = pl.read_parquet(gene_stats_file)
    assert gene_stats.columns == ["genome", "gene", "length", "breadth", "coverage", "ref_ani"]
    assert gene_stats.height == 0

    genome_stats = pl.read_parquet(genome_stats_file)
    assert set(genome_stats["genome"].to_list()) == {"genome1", "genome2"}
    assert "ref_ani" in genome_stats.columns


def test_profile_bam_end_to_end_outputs_without_reference_base_column(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _, observed_reference_flags = _install_fake_profile_subprocess(monkeypatch)
    bam_file, _, bed_file, gene_range_table, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)
    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
        reference_fasta=None,
        gene_range_table=str(gene_range_table),
        stb=stb_lf,
        null_model=pl.scan_parquet(null_model_file),
        output_dir=str(tmp_path),
        num_chunks=2,
        max_concurrency=2,
    )

    prof = pl.read_parquet(tmp_path / "sample_profile.parquet")
    assert "ref_base_bitmask" not in prof.columns
    gene_stats = pl.read_parquet(tmp_path / "sample_gene_stats.parquet")
    genome_stats = pl.read_parquet(tmp_path / "sample_genome_stats.parquet")
    assert "ref_ani" not in gene_stats.columns
    assert "ref_ani" not in genome_stats.columns
    assert observed_reference_flags
    assert not any(observed_reference_flags)


def test_cli_profile_single_bam_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, reference_fasta, bed_file, gene_range_table, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
    out_dir = tmp_path / "cli_out"

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "profile-single",
            "--reference-fasta",
            str(reference_fasta),
            "--bed-file",
            str(bed_file),
            "--bam-file",
            str(bam_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model_file),
            "--gene-range-table",
            str(gene_range_table),
            "--num-chunks",
            "2",
            "--max-concurrency",
            "2",
            "--output-dir",
            str(out_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (out_dir / "sample_profile.parquet").exists()
    assert (out_dir / "sample_gene_stats.parquet").exists()
    assert (out_dir / "sample_genome_stats.parquet").exists()


def test_cli_profile_single_bam_outputs_without_gene_ranges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, reference_fasta, bed_file, _, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
    out_dir = tmp_path / "cli_out_no_gene"

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "profile-single",
            "--reference-fasta",
            str(reference_fasta),
            "--bed-file",
            str(bed_file),
            "--bam-file",
            str(bam_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model_file),
            "--num-chunks",
            "2",
            "--max-concurrency",
            "2",
            "--output-dir",
            str(out_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    prof = pl.read_parquet(out_dir / "sample_profile.parquet")
    assert set(prof["gene"].unique().to_list()) == {"NA"}
    assert "ref_base_bitmask" in prof.columns


def test_cli_profile_single_bam_without_reference_base_column(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, _, bed_file, gene_range_table, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
    out_dir = tmp_path / "cli_out_no_ref"

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "utilities",
            "profile-single",
            "--bed-file",
            str(bed_file),
            "--bam-file",
            str(bam_file),
            "--stb-file",
            str(stb_file),
            "--null-model",
            str(null_model_file),
            "--gene-range-table",
            str(gene_range_table),
            "--num-chunks",
            "2",
            "--max-concurrency",
            "2",
            "--output-dir",
            str(out_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    prof = pl.read_parquet(out_dir / "sample_profile.parquet")
    assert "ref_base_bitmask" not in prof.columns


# ---------------------------------------------------------------------------
# Profiling-asset preparation and caching
# ---------------------------------------------------------------------------

import json


# Building the null model dominates asset-prep runtime; a small max_total_reads
# keeps these caching-focused tests fast without changing what they verify.
_FAST_NULL_MAX_READS = 100


def _write_reference_and_stb(tmp_path):
    """Write a tiny 2-scaffold reference FASTA + matching stb, return their paths."""
    reference = tmp_path / "reference.fna"
    reference.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTTGGGGCC\n")
    stb = tmp_path / "reference.stb"
    stb.write_text("chr1\tgenome1\nchr2\tgenome1\n")
    return reference, stb


def _asset_names(assets_dir):
    return {p.name for p in assets_dir.iterdir()}


def test_prepare_profiling_assets_writes_all_files(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    out_dir = tmp_path / "assets"

    assets = profile.prepare_profiling_assets(
        reference_fasta=reference,
        stb_file=stb,
        output_dir=out_dir,
        max_total_reads=_FAST_NULL_MAX_READS,
    )

    for path in (
        assets.bed_file,
        assets.gene_range_table,
        assets.genome_length_file,
        assets.null_model_file,
        assets.profiling_contract_file,
    ):
        assert path.exists(), path

    # No gene fasta -> empty gene range table, contract records gene hash as missing.
    contract = json.loads(assets.profiling_contract_file.read_text())
    assert contract["gene_hash"] == "NA"
    assert contract["reference_hash"] != "NA"


def test_resolve_profiling_assets_auto_generates_into_run_dir(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    assets = profile.resolve_profiling_assets(
        run_dir=run_dir,
        reference_fasta=reference,
        stb_file=stb,
        max_total_reads=_FAST_NULL_MAX_READS,
    )

    assets_dir = run_dir / profile.DEFAULT_PROFILING_ASSETS_DIRNAME
    assert assets_dir.is_dir()
    assert assets.null_model_file == assets_dir / profile.ASSET_NULL_MODEL_FILENAME
    assert assets.bed_file == assets_dir / profile.ASSET_BED_FILENAME
    assert (assets_dir / profile.ASSET_CACHE_MANIFEST_FILENAME).exists()


def test_resolve_profiling_assets_reuses_cache_when_inputs_unchanged(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    first = profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    null_model_mtime = first.null_model_file.stat().st_mtime_ns

    # Second call with identical inputs must NOT regenerate (mtime unchanged).
    profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    assert first.null_model_file.stat().st_mtime_ns == null_model_mtime


def test_resolve_profiling_assets_regenerates_when_null_model_param_changes(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    first = profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb, error_rate=0.001,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    original_mtime = first.null_model_file.stat().st_mtime_ns

    # A different null-model parameter invalidates the cache -> regeneration.
    profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb, error_rate=0.05,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    assert first.null_model_file.stat().st_mtime_ns != original_mtime


def test_resolve_profiling_assets_force_prepare_regenerates(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    first = profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    original_mtime = first.bed_file.stat().st_mtime_ns

    profile.resolve_profiling_assets(
        run_dir=run_dir, reference_fasta=reference, stb_file=stb, force_prepare=True,
        max_total_reads=_FAST_NULL_MAX_READS,
    )
    assert first.bed_file.stat().st_mtime_ns != original_mtime


def test_resolve_profiling_assets_explicit_paths_win_and_skip_generation(tmp_path):
    reference, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    explicit_null = tmp_path / "my_null.parquet"
    explicit_bed = tmp_path / "my.bed"
    explicit_len = tmp_path / "my_len.parquet"
    for p in (explicit_null, explicit_bed, explicit_len):
        p.write_text("placeholder")

    assets = profile.resolve_profiling_assets(
        run_dir=run_dir,
        reference_fasta=reference,
        stb_file=stb,
        null_model_file=explicit_null,
        bed_file=explicit_bed,
        genome_length_file=explicit_len,
        max_total_reads=_FAST_NULL_MAX_READS,
    )

    # All three required assets supplied -> nothing generated, no assets dir.
    assert assets.null_model_file == explicit_null
    assert assets.bed_file == explicit_bed
    assert assets.genome_length_file == explicit_len
    assert not (run_dir / profile.DEFAULT_PROFILING_ASSETS_DIRNAME).exists()


def test_resolve_profiling_assets_requires_reference_for_autogen(tmp_path):
    _, stb = _write_reference_and_stb(tmp_path)
    run_dir = tmp_path / "run"

    with pytest.raises(ValueError, match="reference-fasta is required"):
        profile.resolve_profiling_assets(
            run_dir=run_dir,
            reference_fasta=None,
            stb_file=stb,
        )
