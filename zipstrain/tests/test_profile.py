from zipstrain import profile
from zipstrain import cli as cli_module
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
    

@pytest.fixture(scope="module")
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()
    
@pytest.mark.parametrize("min_cov,freq_threshold", [
    (5, 0.8),
    (10, 0.5)
])
def test_get_strain_hetrogeneity(profile_1, stb, min_cov, freq_threshold):
    het_profile=profile.get_strain_hetrogeneity(profile_1, stb, min_cov=min_cov, freq_threshold=freq_threshold).collect().rows_by_key("genome",unique=True,named=True)
    genome1_freq=[
        a_chr1 + a_chr2 ,
        t_chr1 + t_chr2 ,
        c_chr1 + c_chr2 ,
        g_chr1 + g_chr2 ,
    ]
    genome2_freq=[
        a_chr3,
        t_chr3,
        c_chr3,
        g_chr3,
    ]
    assert het_profile["genome1"][f"total_sites_at_{min_cov}_coverage"]==sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov)
    assert het_profile["genome2"][f"total_sites_at_{min_cov}_coverage"]==sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov)
    assert het_profile["genome1"]["heterogeneous_sites"]== sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome1_freq))[pos])/sum(list(zip(*genome1_freq))[pos]) < freq_threshold)
    assert het_profile["genome2"]["heterogeneous_sites"]== sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome2_freq))[pos])/sum(list(zip(*genome2_freq))[pos]) < freq_threshold)
    
    assert het_profile["genome1"]["strain_heterogeneity"]== sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome1_freq))[pos])/sum(list(zip(*genome1_freq))[pos]) < freq_threshold) / \
                                                        sum(1 for pos in range(len(genome1_freq[0])) 
                                                        if sum(list(zip(*genome1_freq))[pos]) >= min_cov)
    
    assert het_profile["genome2"]["strain_heterogeneity"]== sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov and 
                                                           max(list(zip(*genome2_freq))[pos])/sum(list(zip(*genome2_freq))[pos]) < freq_threshold) / \
                                                        sum(1 for pos in range(len(genome2_freq[0])) 
                                                        if sum(list(zip(*genome2_freq))[pos]) >= min_cov)


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


def _write_profile_test_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, pl.LazyFrame]:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("")
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
    return bam_file, bed_file, gene_range_table, stb_file, null_model_file, stb_lf


def _install_fake_profile_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    observed_limits: list[int | None] = []

    class _FakeStream:
        def __init__(self, *, lines: list[bytes] | None = None, blob: bytes = b""):
            self._lines = list(lines or [])
            self._line_index = 0
            self._blob = blob

        async def readline(self):
            if self._line_index >= len(self._lines):
                return b""
            line = self._lines[self._line_index]
            self._line_index += 1
            return line

        async def read(self):
            if self._lines:
                if self._line_index >= len(self._lines):
                    return b""
                remaining = b"".join(self._lines[self._line_index :])
                self._line_index = len(self._lines)
                return remaining
            blob = self._blob
            self._blob = b""
            return blob

    class _FakeProc:
        def __init__(self, *, stdout_lines: list[bytes] | None = None, stderr: bytes = b"", returncode: int = 0):
            self.returncode = returncode
            self.stdout = _FakeStream(lines=stdout_lines)
            self.stderr = _FakeStream(blob=stderr)

        async def communicate(self):
            return await self.stdout.read(), await self.stderr.read()

        async def wait(self):
            return self.returncode

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

    def _fake_mpileup_line(row: dict) -> bytes:
        bases = row["bases"]
        depth = len(bases)
        return f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{depth}\t{bases}\t*\n".encode()

    async def _fake_create_subprocess_shell(command: str, stdout=None, stderr=None, cwd=None, limit=None):
        if command.startswith("samtools mpileup"):
            observed_limits.append(limit)
            bed_match = re.search(r"-l\s+([^\s]+)", command)
            assert bed_match is not None, command
            scaffolds = _read_scaffolds_from_bed(Path(bed_match.group(1)))
            lines = []
            for scaffold in scaffolds:
                lines.extend(_fake_mpileup_line(row) for row in _fake_mpileup_rows(scaffold))
            return _FakeProc(stdout_lines=lines)
        if "process-read-locs" in command:
            out_match = re.search(r"--output-file\s+([^\s]+)", command)
            assert out_match is not None, command
            output_file = Path(cwd) / out_match.group(1)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            bed_match = re.search(r"-L\s+([^\s]+)", command)
            assert bed_match is not None, command
            scaffolds = _read_scaffolds_from_bed(Path(bed_match.group(1)))
            rows = []
            for scaffold in scaffolds:
                rows.extend(_fake_read_loc_rows(scaffold))
            pl.DataFrame(rows).write_parquet(output_file)
            return _FakeProc()
        raise AssertionError(f"Unexpected command in fake subprocess: {command}")

    monkeypatch.setattr(profile.asyncio, "create_subprocess_shell", _fake_create_subprocess_shell)
    return observed_limits


def test_profile_bam_end_to_end_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    observed_limits = _install_fake_profile_subprocess(monkeypatch)
    bam_file, bed_file, gene_range_table, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)
    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
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
    assert all(limit == profile.MPILEUP_ASYNCIO_STREAM_LIMIT_BYTES for limit in observed_limits)
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
    assert g2["length"] == 3
    assert g2["breadth"] == pytest.approx(1.0)
    assert g2["coverage"] == pytest.approx(7.0)

    genome_stats = pl.read_parquet(genome_stats_file)
    assert set(genome_stats["genome"].to_list()) == {"genome1", "genome2"}
    by_genome = genome_stats.rows_by_key("genome", unique=True, named=True)
    assert by_genome["genome1"]["coverage"] == pytest.approx(8.0)
    assert by_genome["genome2"]["coverage"] == pytest.approx(7.0)
    assert by_genome["genome1"]["breadth"] == pytest.approx(1.0)
    assert by_genome["genome2"]["breadth"] == pytest.approx(1.0)


def test_profile_bam_end_to_end_outputs_without_gene_ranges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, bed_file, _, _, null_model_file, stb_lf = _write_profile_test_inputs(tmp_path)
    profile.profile_bam(
        bed_file=str(bed_file),
        bam_file=str(bam_file),
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
    assert gene_stats.columns == ["genome", "gene", "length", "breadth", "coverage"]
    assert gene_stats.height == 0

    genome_stats = pl.read_parquet(genome_stats_file)
    assert set(genome_stats["genome"].to_list()) == {"genome1", "genome2"}


def test_cli_profile_single_bam_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, bed_file, gene_range_table, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
    out_dir = tmp_path / "cli_out"

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
    assert (out_dir / "sample_profile.parquet").exists()
    assert (out_dir / "sample_gene_stats.parquet").exists()
    assert (out_dir / "sample_genome_stats.parquet").exists()


def test_cli_profile_single_bam_outputs_without_gene_ranges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _install_fake_profile_subprocess(monkeypatch)
    bam_file, bed_file, _, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
    out_dir = tmp_path / "cli_out_no_gene"

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
    bam_file, bed_file, gene_range_table, stb_file, null_model_file, _ = _write_profile_test_inputs(tmp_path)
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
            "--no-include-reference-base",
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
