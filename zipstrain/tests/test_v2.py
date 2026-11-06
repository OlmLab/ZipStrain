"""V2 contract, overlap, read-phasing, and scientific regression tests."""

from pathlib import Path
from collections import Counter
import json
import shutil
import subprocess

import polars as pl
import pytest
from click.testing import CliRunner
from polars.testing import assert_frame_equal

from zipstrain import cli, codons, compare, gene_info as gi, profile, utils


@pytest.fixture
def assets(tmp_path):
    root = tmp_path / "assets"
    root.mkdir()
    reference = root / "ref.fa"
    reference.write_text(">chr1\nAAAGGATTTCCCAAAGGTTAA\n>chr2\nGGAGGAGGA\n")
    genes = root / "genes.fa"
    genes.write_text(
        ">chr1_1 # 1 # 9 # 1 # partial=00\nAAAGGATTT\n"
        ">chr1_2 # 4 # 12 # 1 # partial=00\nGGATTTCCC\n"
        ">chr1_3 # 1 # 9 # -1 # partial=00\nAAATCCTTT\n"
        ">chr2_1 # 1 # 9 # 1 # partial=00\nGGAGGAGGA\n")
    stb = root / "ref.stb"
    stb.write_text("chr1\tgenomeA\nchr2\tgenomeB\n")
    catalog = gi.write_gene_info(reference, genes, stb, root / "gene_info.parquet")
    return root, reference, genes, stb, catalog


def write_sample(root, catalog, name="sample", sequence="AAAGGATTTCCCAAAGGTTAA", depth=10):
    rows = [("chr1", "genomeA", i+1, *[depth if base == b else 0 for b in "ACGT"], {"A":1,"C":2,"G":4,"T":8}[base])
            for i, base in enumerate(sequence)]
    schema = {"chrom":pl.String,"genome":pl.String,"pos":pl.Int64,"A":pl.Int64,"C":pl.Int64,"G":pl.Int64,"T":pl.Int64,"ref_base_bitmask":pl.UInt8}
    frame = pl.DataFrame(rows, schema=schema, orient="row")
    path = root / f"{name}_profile.parquet"
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata[profile.PROFILE_SORTED_METADATA_KEY] = profile.PROFILE_SORTED_METADATA_VALUE
    frame.write_parquet(path, metadata=metadata)
    return path


def write_counts(root, catalog, name, rows, *, source_profile=None):
    path = root / f"{name}_codon_counts.parquet"
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata["zipstrain_codon_evidence"] = "read_phased_v1"
    if source_profile is not None:
        metadata["zipstrain_source_profile_sha256"] = utils.sha256_file(source_profile)
    pl.DataFrame(rows, schema=codons.COUNT_SCHEMA, orient="row").write_parquet(path, metadata=metadata)
    return path


def write_bam(root, records):
    sam = root / "reads.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:21\n@SQ\tSN:chr2\tLN:9\n" + "\n".join(records) + "\n")
    bam = root / "reads.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
    return bam


def test_fasta_and_content_identity(tmp_path):
    path = tmp_path / "a.fa"
    path.write_text(">a description\naa\ncc\n\n>b\ngg\n")
    assert list(gi.fasta_records(path)) == [("a description", "AACC"), ("b", "GG")]
    assert gi.content_digest({"a":1,"b":2}) == gi.content_digest({"b":2,"a":1})
    for content in ["AAA\n", ">a\n>b\nAAA\n", ">a\n"]:
        path.write_text(content)
        with pytest.raises(ValueError):
            list(gi.fasta_records(path))


def test_catalog_schema_ids_reference_codons_and_tampering(assets):
    root, reference, genes, stb, catalog = assets
    frame = gi.read_gene_info(catalog)
    assert frame.schema == gi.GENE_SCHEMA
    assert frame.height == 4 and frame["gene_id"].n_unique() == 4
    codon_frame = pl.read_parquet(root / "reference_codons.parquet")
    reverse_id = frame.filter(pl.col("gene") == "chr1_3")["gene_id"][0]
    assert codon_frame.filter(pl.col("gene_id") == reverse_id)["reference_codon"].to_list() == ["AAA", "TCC", "TTT"]
    old = frame.sort("gene")
    records = list(gi.fasta_records(genes))[::-1]
    genes.write_text("".join(f">{h}\n{s}\n" for h,s in records))
    gi.write_gene_info(reference, genes, stb, catalog)
    assert gi.read_gene_info(catalog).sort("gene").equals(old)
    assert gi.read_gene_ranges(catalog).collect().columns == ["gene", "scaffold", "start", "end"]
    legacy = root / "ranges.tsv"
    legacy.write_text("g\tchr1\t1\t9\n")
    assert gi.read_gene_ranges(legacy).collect()["gene"].to_list() == ["g"]
    metadata = utils._read_custom_parquet_metadata(catalog)
    frame.with_columns((pl.col("start")+1).alias("start")).write_parquet(catalog, metadata=metadata)
    with pytest.raises(ValueError, match="identity"):
        gi.read_gene_info(catalog)


@pytest.mark.parametrize("header", ["bad", "chr1_1 # 0 # 9 # 1", "chr1_1 # 1 # 90 # 1",
    "chr1_1 # 1 # 9 # 0", "chr1_1 # 1 # 9 # 1 # phase=3", "chr1_1 # 1 # 9 # 1 # transl_table=4",
    "chr1_1 # 1 # 9 # 1 # partial=xy", "missing_1 # 1 # 9 # 1"])
def test_catalog_rejects_bad_annotations(assets, header):
    root, reference, genes, stb, catalog = assets
    genes.write_text(f">{header}\nAAAGGATTT\n")
    with pytest.raises(ValueError):
        gi.write_gene_info(reference, genes, stb, catalog)


def test_empty_catalog_and_phase_partial(assets):
    root, reference, genes, stb, catalog = assets
    gi.write_gene_info(reference, None, stb, catalog)
    assert gi.read_gene_info(catalog).height == 0
    assert pl.read_parquet(root / "reference_codons.parquet").schema == gi.CODON_SCHEMA
    genes.write_text(">chr1_1 # 1 # 9 # -1 # phase=1;partial=10\nAATCCTTT\n")
    gi.write_gene_info(reference, genes, stb, catalog)
    gene = gi.read_gene_info(catalog).row(0, named=True)
    assert gene["phase"] == 1 and gene["partial_3p"] and not gene["partial_5p"]
    assert pl.read_parquet(root / "reference_codons.parquet")["reference_codon"].to_list() == ["AAT", "CCT"]


def test_overlap_statistics_and_missing_reference(assets):
    root, _, _, _, catalog = assets
    path = write_sample(root, catalog)
    stats = gi.gene_statistics(pl.scan_parquet(path), catalog).collect()
    assert stats.columns == list(gi.GENE_STATS_COLUMNS)
    assert stats.filter(pl.col("genome") == "genomeA")["coverage"].to_list() == [10.0]*3
    assert stats.filter(pl.col("genome") == "genomeA")["breadth"].to_list() == [1.0]*3
    assert stats.filter(pl.col("genome") == "genomeA")["ref_ani"].to_list() == [100.0]*3
    assert stats.filter(pl.col("genome") == "genomeB")["breadth"].to_list() == [0.0]
    assert stats["dN"].null_count() == 4
    absent = gi.gene_statistics(pl.scan_parquet(path).drop("ref_base_bitmask"), catalog).collect()
    assert absent["ref_ani"].null_count() == 4
    genes = gi.read_gene_info(catalog)
    sums = gi.range_sums(pl.scan_parquet(path).select("chrom", "pos", pl.lit(1).alias("sites")), genes, ["sites"]).collect()
    assert sorted(sums["sites"].to_list()) == [0.0, 9.0, 9.0, 9.0]
    with pytest.raises(ValueError):
        gi.gene_statistics(pl.scan_parquet(path), catalog, min_cov=0)


def test_profile_identity_and_finalizer(assets):
    root, _, _, _, catalog = assets
    path = write_sample(root, catalog)
    metadata = gi.validate_profile_assets([path], catalog)
    assert metadata[gi.SCHEMA_KEY] == "2"
    legacy = pl.read_parquet(path).with_columns(pl.lit("wrong").alias("gene"))
    legacy.write_parquet(path, metadata={"keep": "yes"})
    profile.finish_profile_outputs(path, catalog)
    assert "gene" not in pl.read_parquet_schema(path)
    assert pl.read_parquet_metadata(path)["keep"] == "yes"
    gi.validate_profile_assets([path], catalog)
    assert pl.read_parquet(root / "sample_gene_stats.parquet").columns == list(gi.GENE_STATS_COLUMNS)
    with pytest.raises(ValueError, match="require"):
        profile.finish_profile_outputs(path, collect_codons=True)
    pl.read_parquet(path).write_parquet(path)
    with pytest.raises(ValueError, match="Incompatible|schema"):
        gi.validate_profile_assets([path], catalog)


def test_reference_identity_and_resumable_collection(assets):
    root, reference, _, stb, catalog = assets
    assert gi.reference_identity(reference, profile.read_stb(stb)) == pl.read_parquet_metadata(catalog)[gi.REFERENCE_ID]
    a, b = write_sample(root, catalog, "a"), write_sample(root, catalog, "b")
    run = root / "compare"
    output = compare.compare_profile_collection([("a", a), ("b", b)], run, preprofile_dir=root)
    assert output == run / "all_comparisons.parquet"
    assert pl.read_parquet(output).height == 3
    assert list(run.glob("*.parquet")) == [output]
    cached = list((run / "intermediate_files" / "v2_pairs").glob("*.parquet"))
    timestamps = {path: path.stat().st_mtime_ns for path in cached}
    c = write_sample(root, catalog, "c")
    compare.compare_profile_collection([("a", a), ("b", b), ("c", c)], run, preprofile_dir=root)
    assert pl.read_parquet(output).height == 9
    assert all(path.stat().st_mtime_ns == timestamp for path, timestamp in timestamps.items())
    with pytest.raises(ValueError, match="two"):
        compare.compare_profile_collection([("a", a)], run)
    result = CliRunner().invoke(cli.cli, ["utilities", "compare-pair", "--profile-1", str(a), "--profile-2", str(b),
        "--preprofile-dir", str(root), "--output-dir", str(root / "pair")])
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert list((root / "pair").glob("*.parquet")) == [root / "pair" / "comparisons.parquet"]
    assert pl.read_parquet(root / "pair" / "comparisons.parquet").height == 3


def test_read_alignment_and_observed_codons(assets):
    root, _, _, _, catalog = assets
    bases = codons.aligned_read_bases(1, "2S3M1I2M1D1M", "NNAAACGGT", "I"*9)
    assert bases[1] == ("A", 40, 2) and bases[7][0] == "T" and 6 not in bases
    with pytest.raises(ValueError):
        codons.aligned_read_bases(1, "9Q", "AAA", "III")
    with pytest.raises(ValueError):
        codons.aligned_read_bases(1, "3M", "AAA", "*")
    with pytest.raises(ValueError):
        codons.aligned_read_bases(1, "4M", "AAA", "III")
    gene = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1").row(0, named=True)
    assert list(codons.read_codon_observations(bases, gene)) == [(0, "AAA", 40)]
    reverse = {**gene, "strand": -1}
    full = codons.aligned_read_bases(1, "9M", "AAAGGATTT", "I"*9)
    assert [row[1] for row in codons.read_codon_observations(full, reverse)] == ["AAA", "TCC", "TTT"]
    assert list(codons.read_codon_observations({}, gene)) == []
    assert list(codons.read_codon_observations(full, gene, min_baseq=41)) == []


def test_observed_codon_counts_and_mate_conflicts(assets):
    root, _, _, _, catalog = assets
    path = write_sample(root, catalog)
    records = [
        "pair\t65\tchr1\t1\t60\t9M\t=\t1\t0\tAAAGGATTT\tIIIIIIIII",
        "pair\t129\tchr1\t1\t60\t9M\t=\t1\t0\tAAAGGATTT\tIIIIIIIII",
        "low\t0\tchr1\t1\t1\t9M\t*\t0\t0\tAAAGGATTT\tIIIIIIIII",
        "dup\t1024\tchr1\t1\t60\t9M\t*\t0\t0\tAAAGGATTT\tIIIIIIIII",
        "other\t0\tchr1\t1\t60\t9M\t*\t0\t0\tAAAGGATTT\tIIIIIIIII",
    ]
    bam = write_bam(root, records)
    output = codons.write_codon_counts(bam, path, catalog, root / "counts.parquet", min_mapq=10)
    counts = codons.read_codon_counts(output, catalog).collect()
    assert counts["count"].to_list() == [2]*8  # gene 1/3 have 3 codons, gene 2 has 2 covered codons
    assert counts.schema == codons.COUNT_SCHEMA
    assert not output.with_suffix(".raw.parquet").exists()


def test_sense_codon_counting_and_shared_nonreference_context():
    assert codons.codon_changes("GGA", "GGT") == (1.0, 0.0)
    assert codons.codon_changes("AAA", "AGA") == (0.0, 1.0)
    assert codons.codon_changes("AAA", "AGG") == (1.0, 1.0)
    assert codons.codon_changes("AAA", "AAA") == (0.0, 0.0)
    assert codons.synonymous_sites("GGA") == pytest.approx(1.0)
    assert codons.synonymous_sites("AAA") == pytest.approx(1/3)
    with pytest.raises(ValueError):
        codons.codon_changes("TAA", "AAA")
    with pytest.raises(ValueError):
        codons.synonymous_sites("NNN")


def test_pair_alignment_genome_gene_comparison_and_annotation_guards(assets):
    root, _, _, _, catalog = assets
    a = write_sample(root, catalog, "a")
    b = write_sample(root, catalog, "b", "AAAGGTTTTCCCAAAGGTTAA")
    aligned = compare.align_profiles(a, b).collect()
    assert aligned.height == 21 and "A_2" in aligned.columns and "gene" not in aligned.columns
    genome = compare.compare_profiles(a, b, preprofile_dir=root, calculate="genome_ani").collect()
    assert "gene_id" not in genome.columns and "max_consecutive_length" not in genome.columns
    assert set(genome["genome"]) == {"genomeA", "genomeB"}
    assert genome.filter(pl.col("genome") == "genomeA")["total_positions"][0] == 21
    gene_only = compare.compare_profiles(a, b, preprofile_dir=root, calculate="gene_ani").collect()
    assert "genome_ani" not in gene_only.columns
    gene = compare.compare_profiles(a, b, preprofile_dir=root,
        calculate="genome_ani+max_identical_stretch+gene_ani").collect()
    assert gene.height == 3 and gene["gene_total_positions"].to_list() == [9]*3
    assert gene["total_positions"].to_list() == [21]*3
    assert gene["share_allele_pos"].to_list() == [20]*3
    assert gene["max_consecutive_length"].to_list() == [15]*3
    assert gene.select("sample_1", "sample_2", "genome", "gene_id").n_unique() == 3
    assert gene["gene_ani"].to_list() == pytest.approx([100*8/9]*3)
    with pytest.raises(ValueError, match="requires"):
        compare.compare_profiles(a, b, calculate="gene_ani")
    with pytest.raises(ValueError, match="requires"):
        compare.compare_profiles(a, b, preprofile_dir=root, dnds=True)
    frame = pl.read_parquet(b)
    pl.concat([frame, frame.head(1)]).write_parquet(b, metadata=utils._read_custom_parquet_metadata(catalog))
    with pytest.raises(ValueError, match="Duplicate"):
        compare.align_profiles(a, b)


def test_v2_matrix_compare_exports_one_gene_table(assets):
    pytest.importorskip("h5py")
    root, _, _, stb, catalog = assets
    a, b = write_sample(root, catalog, "a"), write_sample(root, catalog, "b")
    bed = root / "regions.bed"
    bed.write_text("chr1\t0\t21\nchr2\t0\t9\n")
    profiles = root / "profiles.csv"
    profiles.write_text(f"profile_name,profile_location\na,{a}\nb,{b}\n")
    run = root / "matrix_comparison"
    result = CliRunner().invoke(cli.cli, [
        "compare", "--method", "matrix", "--profile-db", str(profiles),
        "--stb-file", str(stb), "--bed-file", str(bed), "--preprofile-dir", str(root),
        "--run-dir", str(run), "--calculate", "gene_ani+genome_ani+max_identical_stretch",
        "--min-gene-compare-len", "1",
    ])
    assert result.exit_code == 0, result.output + repr(result.exception)
    output = run / "all_comparisons.parquet"
    assert list(run.glob("*.parquet")) == [output]
    assert not output.with_suffix(".csv").exists()
    frame = pl.read_parquet(output)
    assert frame.height == 3
    assert frame.select("sample_1", "sample_2", "genome", "gene_id").n_unique() == 3
    assert frame["gene_ani"].to_list() == frame["genome_ani"].to_list() == [100.0] * 3
    assert frame["max_consecutive_length"].to_list() == [21] * 3
    assert pl.read_parquet_metadata(output)["zipstrain_compare_kind"] == "gene"


@pytest.mark.parametrize(("integration", "calculate"), [("weighted", "all"), ("consensus", "genome_ani")])
def test_matrix_selection_matches_table_and_reuses_pairs(assets, integration, calculate):
    pytest.importorskip("h5py")
    root, _, _, stb, catalog = assets
    ids = dict(gi.read_gene_info(catalog).select("gene", "gene_id").iter_rows())
    g1, g2 = ids["chr1_1"], ids["chr1_2"]
    a = write_sample(root, catalog, "a")
    b = write_sample(root, catalog, "b", "AAAGGTTTTCCCAAAGGTTAA")
    left = [(g1, 0, "AAA", 5), (g1, 0, "AAG", 5), (g1, 0, "TAG", 5),
            (g1, 1, "GGA", 10), (g1, 2, "TTT", 10), (g2, 0, "GGA", 10)]
    right = [(g1, 0, "AAC", 5), (g1, 0, "AAG", 5), (g1, 0, "AAA", 5),
             (g1, 1, "GGT", 10), (g1, 2, "CTC", 10), (g2, 0, "GGT", 2)]
    write_counts(root, catalog, "a", left, source_profile=a)
    write_counts(root, catalog, "b", right, source_profile=b)
    bed = root / "regions.bed"
    bed.write_text("chr1\t0\t21\nchr2\t0\t9\n")
    profiles = root / "profiles.csv"
    profiles.write_text(f"profile_name,profile_location\na,{a}\nb,{b}\n")
    common = ["compare", "--profile-db", str(profiles), "--preprofile-dir", str(root),
              "--dnds", "--allele-integration", integration, "--calculate", calculate,
              "--min-gene-compare-len", "1"]
    matrix = root / "matrix"
    matrix_args = [*common, "--method", "matrix", "--stb-file", str(stb), "--bed-file", str(bed),
                   "--run-dir", str(matrix)]
    table = root / "table"
    table_args = [*common, "--method", "table", "--run-dir", str(table)]
    runner = CliRunner()

    def check_outputs():
        for args in (matrix_args, table_args):
            result = runner.invoke(cli.cli, args)
            assert result.exit_code == 0, result.output + repr(result.exception)
        keys = ["sample_1", "sample_2", "gene_id"]
        actual = pl.read_parquet(matrix / "all_comparisons.parquet").sort(keys)
        expected = pl.read_parquet(table / "all_comparisons.parquet").sort(keys)
        assert_frame_equal(actual, expected)
        assert actual.filter(pl.col("gene_id") == g1)["pS"].min() > 0
        assert actual.filter(pl.col("gene_id") == g2)["pS"].null_count() > 0
        return actual

    assert check_outputs().height == 3
    cache = matrix / "intermediate_files" / "matrix_selection"
    before = {path: path.stat().st_mtime_ns for path in cache.glob("pair-*.parquet")}
    assert len(before) == 1
    if integration == "weighted":
        from zipstrain.resource_limits import connect_duckdb

        # Simulate a result database made before per-gene counts were saved.
        with connect_duckdb(str(cache.parent / "matrix_compare.duckdb")) as conn:
            conn.execute("ALTER TABLE matrix_compare_gene_results DROP COLUMN gene_total_positions")
            conn.execute("ALTER TABLE matrix_compare_gene_results DROP COLUMN gene_share_allele_pos")
        assert check_outputs().height == 3
    c = write_sample(root, catalog, "0")
    write_counts(root, catalog, "0", left, source_profile=c)
    with profiles.open("a") as handle:
        handle.write(f"0,{c}\n")
    assert check_outputs().height == 9
    assert all(path.stat().st_mtime_ns == stamp for path, stamp in before.items())
    assert len(list(cache.glob("pair-*.parquet"))) == 3

    result = runner.invoke(cli.cli, [*matrix_args, "--min-gene-compare-len", "10"])
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert pl.read_parquet(matrix / "all_comparisons.parquet").is_empty()
    write_counts(root, catalog, "0", left, source_profile=b)
    result = runner.invoke(cli.cli, matrix_args)
    assert result.exit_code != 0 and "does not match its sample position profile" in result.output


def test_alignment_uses_parquet_order_marker_without_scanning_rows(tmp_path, monkeypatch):
    def make_profile(chromosomes, positions, counts):
        return pl.DataFrame({
            "chrom": chromosomes, "genome": ["g" + chrom for chrom in chromosomes],
            "pos": positions, "A": counts, "C": 0, "G": 0, "T": 0,
        })

    left = make_profile(["a", "c", "c", "z"], [1, 1, 2**32, 1], [10, 11, 12, 13])
    right = make_profile(["b", "c", "c", "z"], [1, 1, 2**32, 1], [20, 21, 22, 23])
    sorted_left = tmp_path / "left.parquet"
    sorted_right = tmp_path / "right.parquet"
    unsorted_right = tmp_path / "unmarked.parquet"
    metadata = {profile.PROFILE_SORTED_METADATA_KEY: profile.PROFILE_SORTED_METADATA_VALUE}
    left.write_parquet(sorted_left, metadata=metadata)
    right.write_parquet(sorted_right, metadata=metadata)
    right.reverse().write_parquet(unsorted_right)

    def unexpected_collect(*args, **kwargs):
        pytest.fail("Planning a marked pair must not scan profile rows")

    with monkeypatch.context() as patch:
        patch.setattr(pl.LazyFrame, "collect", unexpected_collect)
        plan = compare.align_profiles(sorted_left, sorted_right)
    ordered = plan.collect(engine="streaming")
    fallback = compare.align_profiles(sorted_left, unsorted_right).collect(engine="streaming")
    assert ordered.equals(fallback)
    assert ordered.select("chrom", "pos", "A", "A_2").rows() == [
        ("c", 1, 11, 21), ("c", 2**32, 12, 22), ("z", 1, 13, 23),
    ]


def test_selection_evidence_and_pairwise_regression(assets):
    root, _, _, _, catalog = assets
    gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    a = write_counts(root, catalog, "a", [(gene_id, 1, "GGA", 10)])
    b = write_counts(root, catalog, "b", [(gene_id, 1, "GGT", 10)])
    pair = codons.pairwise_selection(a, b, catalog).collect().row(0, named=True)
    assert pair["dN"] == 0 and pair["dS"] == 1 and pair["dN_dS"] == 0
    ref = codons.reference_selection(b, catalog).collect().row(0, named=True)
    assert ref["dN"] == 0 and ref["dS"] == 1 and ref["pS"] == 0
    poly = write_counts(root, catalog, "poly", [(gene_id, 1, "GGA", 10), (gene_id, 1, "GGT", 2)])
    stats = codons.reference_selection(poly, catalog).collect().row(0, named=True)
    assert stats["dS"] == 0 and stats["pS"] == 1 and stats["pN"] == 0
    assert stats["dN_dS"] is None
    with pytest.raises(ValueError):
        codons.reference_selection(a, catalog, min_freq=2)
    pl.read_parquet(a).write_parquet(a, metadata=utils._read_custom_parquet_metadata(catalog))
    with pytest.raises(ValueError, match="read-phased"):
        codons.read_codon_counts(a, catalog)


def test_branch_allele_integration_defaults_ties_and_frequencies():
    source = Counter({"A": 7, "G": 2, "C": 1})
    category, events = codons.allele_events(source, None, "A")
    assert category == "snv" and events == [("A", "G", 1.0)]
    assert codons.allele_events(source, None, "A", "weighted")[1] == [("A", "G", .2), ("A", "C", .1)]
    assert codons.allele_events(Counter({"G":5,"C":5}), None, "A")[1] == [("A", "G", .5), ("A", "C", .5)]
    assert codons.allele_events(Counter({"G":10}), None, "A")[0] == "sns"
    assert codons.allele_events(Counter({"A":10}), None, "A")[1] == []
    a, b = Counter({"A":8,"G":2}), Counter({"A":5,"C":5})
    assert codons.allele_events(a,b,"A")[1] == [("A","C",1.0)]
    assert codons.allele_events(a,b,"A","weighted")[1] == [("A","C",.4),("G","A",.1),("G","C",.1)]
    with pytest.raises(ValueError):
        codons.allele_events(a,b,"A","unknown")


def test_site_path_effects_and_selection_frame():
    effects = [codons.site_changes("AAA", "AGG", i) for i in range(3)]
    assert tuple(map(sum, zip(*effects))) == codons.codon_changes("AAA", "AGG")
    assert codons.site_changes("GGA", "GGT", 2) == (1.0, 0.0)
    with pytest.raises(ValueError):
        codons.site_changes("AAA", "AAG", 3)
    with pytest.raises(ValueError):
        codons.site_changes("TAA", "AAA", 0)
    rows = codons.selection_frame({1:[1,2,3,4,2,1],2:[0,0,0,0,0,0]}).collect().sort("gene_id")
    assert rows.row(0, named=True) == dict(gene_id=1,dN=1.,dS=1.,dN_dS=1.,pN=3.,pS=2.,pN_pS=1.5)
    assert all(value is None for value in rows.row(1)[1:])


def test_reference_and_pairwise_branch_integration_modes(assets):
    root, _, _, _, catalog = assets
    gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    variants = write_counts(root,catalog,"mixed",[(gene_id,0,"AAG",8),(gene_id,0,"AAC",2)])
    consensus = codons.reference_selection(variants,catalog).collect().row(0,named=True)
    weighted = codons.reference_selection(variants,catalog,allele_integration="weighted").collect().row(0,named=True)
    assert consensus["pS"] == pytest.approx(3.) and consensus["pN"] == 0
    assert weighted["pS"] == pytest.approx(2.4) and weighted["pN"] == pytest.approx(.2/(8/3))
    reference = write_counts(root,catalog,"ref",[(gene_id,0,"AAA",10)])
    forward = codons.pairwise_selection(reference,variants,catalog,allele_integration="weighted").collect()
    reverse = codons.pairwise_selection(variants,reference,catalog,allele_integration="weighted").collect()
    assert forward.equals(reverse)
    assert forward["pS"][0] == pytest.approx(2.4)
    assert forward["pN"][0] == pytest.approx(.2/(8/3))
    # Both samples share two nonreference bases. Only the last differs, synonymously.
    a = write_counts(root,catalog,"gga",[(gene_id,0,"GGA",10)])
    b = write_counts(root,catalog,"ggt",[(gene_id,0,"GGT",10)])
    result = codons.pairwise_selection(a,b,catalog).collect()
    assert result["dN"][0] == 0 and result["dS"][0] == pytest.approx(3.)
    groups = list(codons.codon_groups(variants,catalog,min_cov=5))
    assert groups == [((gene_id,0),"AAA",Counter({"AAG":8,"AAC":2}))]
    assert list(codons.codon_groups(variants,catalog,min_cov=100)) == []


def test_observed_evidence_rejects_bad_rows_and_reference_mixups(assets):
    root, _, _, _, catalog = assets
    gene_id = gi.read_gene_info(catalog)["gene_id"][0]
    for rows in [[(gene_id,0,"AAA",0)],[(gene_id,0,"NNN",10)],[(gene_id,999,"AAA",10)],
                 [(gene_id,0,"AAA",1),(gene_id,0,"AAA",2)],[(0,0,"AAA",10)]]:
        path = write_counts(root,catalog,"invalid",rows)
        with pytest.raises(ValueError):
            codons.read_codon_counts(path,catalog).collect()
    reference = root / "reference_codons.parquet"
    frame = pl.read_parquet(reference)
    frame.with_columns(pl.lit("AAA").alias("reference_codon")).write_parquet(reference,
        metadata=utils._read_custom_parquet_metadata(reference))
    path = write_counts(root,catalog,"invalid",[(gene_id,0,"AAA",10)])
    with pytest.raises(ValueError, match="Reference codon contents"):
        codons.read_codon_counts(path,catalog)


@pytest.mark.parametrize("backend", ["python", "rust_profiler"])
def test_real_bam_selection_cli_and_collection(assets, backend):
    if backend == "rust_profiler":
        from importlib.metadata import entry_points
        if not entry_points(group="zipstrain.profile_backends", name=backend):
            pytest.skip("The compiled Rust plugin wheel is not installed")
    root, reference, _, stb, catalog = assets
    records = []
    for index in range(15):
        sequence = "AAAGGATTTCCCAAAGGTTAA" if index < 10 else "AAAGGTTTTCCCAAAGGTTAA"
        records.append(f"read{index}\t0\tchr1\t1\t60\t21M\t*\t0\t0\t{sequence}\t{'I'*21}\tNM:i:{int(index>=10)}")
    bam = write_bam(root, records)
    bed = root / "regions.bed"
    bed.write_text("chr1\t0\t21\nchr2\t0\t9\n")
    model = root / "null.parquet"
    pl.DataFrame({"cov":list(range(101)),"max_error_count":[0]*101}).write_parquet(model)
    out = root / "sample_out"
    args = ["utilities","profile-single","--reference-fasta",str(reference),"--bed-file",str(bed),
        "--bam-file",str(bam),"--stb-file",str(stb),"--null-model",str(model),
        "--gene-info-table",str(catalog),"--output-dir",str(out),"--num-chunks","1",
        "--max-concurrency","1","--read-inclusion","all-mapped","--min-read-ani","0",
        "--collect-codons","--allele-integration","weighted","--backend",backend]
    result = CliRunner().invoke(cli.cli,args)
    assert result.exit_code == 0, result.output+repr(result.exception)
    positions = out / "reads_profile.parquet"
    assert not {"gene", "gene_id"} & set(pl.read_parquet_schema(positions))
    counts = out / "reads_codon_counts.parquet"
    stats = pl.read_parquet(out / "reads_gene_stats.parquet")
    assert stats.columns == list(gi.GENE_STATS_COLUMNS)
    gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    row = stats.filter(pl.col("gene_id") == gene_id).row(0,named=True)
    assert row["pS"] == pytest.approx(.2)
    assert row["dN"] == row["dS"] == row["pN"] == 0
    assert pl.read_parquet_metadata(out / "reads_gene_stats.parquet")["zipstrain_allele_integration"] == "weighted"
    assert pl.read_parquet(counts)["count"].min() >= 5
    b = out / "b_profile.parquet"
    b_counts = out / "b_codon_counts.parquet"
    shutil.copyfile(positions,b)
    shutil.copyfile(counts,b_counts)
    pair = root / "pair_dnds"
    result = CliRunner().invoke(cli.cli,["utilities","compare-pair","--profile-1",str(positions),
        "--profile-2",str(b),"--preprofile-dir",str(root),"--output-dir",str(pair),
        "--codon-counts-1",str(counts),"--codon-counts-2",str(b_counts),"--dnds",
        "--allele-integration","weighted"])
    assert result.exit_code == 0,result.output+repr(result.exception)
    genes = pl.read_parquet(pair / "comparisons.parquet")
    assert all(name in genes.columns for name in gi.SELECTION_COLUMNS)
    assert genes.filter(pl.col("gene_id") == gene_id)["pS"][0] > 0
    inputs = root / "profiles.csv"
    inputs.write_text(f"profile_name,profile_location\na,{positions}\nb,{b}\n")
    run = root / "collection"
    result = CliRunner().invoke(cli.cli,["compare","--profile-db",str(inputs),"--run-dir",str(run),
        "--preprofile-dir",str(root),"--dnds","--min-gene-compare-len","1","--no-csv"])
    assert result.exit_code == 0,result.output+repr(result.exception)
    assert pl.read_parquet(run / "all_comparisons.parquet").height == 3


def test_new_function_docstrings_and_guide_entries():
    import ast
    root = Path(__file__).resolve().parents[2]
    guide = (root / "docs" / "v2.md").read_text()
    modules = {"gene_info":None, "codons":None, "matrix_selection":None,
               "profile":{"finish_profile_outputs"},
               "compare":{"align_profiles","compare_profiles","compare_profile_collection"},
               "cli":{"compare_pair"}}
    for module, names in modules.items():
        tree = ast.parse((root / "zipstrain" / "src" / "zipstrain" / f"{module}.py").read_text())
        for node in tree.body:
            if isinstance(node,ast.FunctionDef) and (names is None or node.name in names):
                assert ast.get_docstring(node),f"Missing function documentation: {module}.{node.name}"
                assert f"{module}.{node.name}" in guide,f"Missing guide entry: {module}.{node.name}"


def test_nextflow_profiles_with_shared_catalog_and_observed_codons(assets):
    from .test_nextflow_workflow import _NEXTFLOW_AVAILABLE, _NEXTFLOW_BIN, _nextflow_subprocess_env
    if not _NEXTFLOW_AVAILABLE:
        pytest.skip("Nextflow and Java are required")
    root, reference, genes, stb, _catalog = assets
    records = [f"r{i}\t0\tchr1\t1\t60\t21M\t*\t0\t0\tAAAGGATTTCCCAAAGGTTAA\t{'I'*21}\tNM:i:0" for i in range(10)]
    bam = write_bam(root,records)
    inputs = root / "bams.csv"
    inputs.write_text(f"sample_name,bamfile\nreads,{bam}\n")
    config = root / "local.config"
    config.write_text("docker.enabled = false\nprocess {\n"
        "withName: 'prepare_profile' { cpus = 1; memory = '1 GB' }\n"
        "withName: 'profile_bam' { cpus = 1; memory = '1 GB' }\n}\n")
    output = root / "nextflow_out"
    repository = Path(__file__).resolve().parents[2]
    result = subprocess.run([str(_NEXTFLOW_BIN),"run",str(repository/"zipstrain.nf"),"-c",str(config),
        "--mode","profile","--input_table",str(inputs),"--reference_genome",str(reference),
        "--gene_file",str(genes),"--stb",str(stb),"--output_dir",str(output),
        "--max_total_reads","100","--collect_codons","true","--allele_integration","weighted"],
        cwd=root,env=_nextflow_subprocess_env(),capture_output=True,text=True,timeout=180)
    assert result.returncode == 0,result.stdout+result.stderr
    assert "gene" not in pl.read_parquet_schema(output/"reads_profile.parquet")
    assert (output/"reads_codon_counts.parquet").exists()
    assert (output/"gene_info.parquet").exists() and (output/"reference_codons.parquet").exists()
    assert pl.read_parquet(output/"reads_gene_stats.parquet").columns == list(gi.GENE_STATS_COLUMNS)


def test_mate_disagreement_and_quality_resolution(assets):
    root, _, _, _, catalog = assets
    positions = write_sample(root,catalog)
    frame = pl.read_parquet(positions).with_columns(
        pl.when(pl.col("pos") == 6).then(5).otherwise(pl.col("T")).alias("T"))
    frame.write_parquet(positions,metadata=utils._read_custom_parquet_metadata(catalog))
    records = []
    for name, quality in [("tie","I"),("quality","?")]:
        records.extend([
            f"{name}\t65\tchr1\t1\t60\t9M\t=\t1\t0\tAAAGGATTT\t{quality*9}",
            f"{name}\t129\tchr1\t1\t60\t9M\t=\t1\t0\tAAAGGTTTT\t{'I'*9}"])
    bam = write_bam(root,records)
    counts = codons.write_codon_counts(bam,positions,catalog,root/"resolved.parquet")
    gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    rows = pl.read_parquet(counts).filter(pl.col("gene_id") == gene_id).sort("codon_index")
    assert rows.select("codon_index","codon","count").rows() == [(0,"AAA",2),(1,"GGT",1),(2,"TTT",2)]


def test_phasing_preserves_different_populations_with_same_position_counts(assets):
    root, _, _, _, catalog = assets
    positions = write_sample(root,catalog)
    frame = pl.read_parquet(positions).with_columns(
        *[pl.when(pl.col("pos") <= 3).then(5).otherwise(pl.col(b)).alias(b) for b in "AG"])
    frame.write_parquet(positions,metadata=utils._read_custom_parquet_metadata(catalog))
    gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    observed = []
    for label, haplotypes in [("a",["AAA","GGG"]),("b",["AGA","GAG"])]:
        folder = root/label
        folder.mkdir()
        records = [f"r{i}\t0\tchr1\t1\t60\t3M\t*\t0\t0\t{haplotypes[i%2]}\tIII" for i in range(10)]
        bam = write_bam(folder,records)
        counts = codons.write_codon_counts(bam,positions,catalog,folder/"counts.parquet")
        rows = pl.read_parquet(counts).filter((pl.col("gene_id") == gene_id) & (pl.col("codon_index") == 0))
        observed.append(dict(rows.select("codon","count").rows()))
    assert observed == [{"AAA":5,"GGG":5},{"AGA":5,"GAG":5}]


def test_v2_stretches_exclude_mismatches_and_empty_intersections():
    source = pl.DataFrame({"chrom":["chr"]*5,"genome":["g"]*5,"pos":[1,2,3,4,5],
                          "A":[10]*5,"C":[0]*5,"G":[0]*5,"T":[0]*5}).lazy()
    target = source.with_columns(pl.when(pl.col("pos") == 1).then(0).otherwise(10).alias("A"),
                                 pl.when(pl.col("pos") == 1).then(10).otherwise(0).alias("T"))
    row = compare.compare_profiles(source,target).collect().row(0,named=True)
    assert row["max_consecutive_length"] == 4 and row["genome_ani"] == 80
    mismatch = source.with_columns(pl.lit(0).alias("A"),pl.lit(10).alias("T"))
    row = compare.compare_profiles(source,mismatch).collect().row(0,named=True)
    assert row["max_consecutive_length"] == 0 and row["genome_ani"] == 0
    no_overlap = source.with_columns((pl.col("pos")+100).alias("pos"))
    row = compare.compare_profiles(source,no_overlap).collect().row(0,named=True)
    assert row["total_positions"] == row["max_consecutive_length"] == 0 and row["genome_ani"] is None


def test_profile_validation_rejects_filter_and_reference_mismatches(assets):
    root, _, _, _, catalog = assets
    a,b = write_sample(root,catalog,"a"),write_sample(root,catalog,"b")
    frame = pl.read_parquet(b)
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata[gi.REFERENCE_ID] = "different"
    frame.write_parquet(b,metadata=metadata)
    with pytest.raises(ValueError,match="reference"):
        compare.align_profiles(a,b)
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata["zipstrain_profile_filters"] = '{"min_freq":0.1}'
    frame.write_parquet(b,metadata=metadata)
    with pytest.raises(ValueError,match="filters"):
        compare.align_profiles(a,b)
    with pytest.raises(ValueError,match="integer"):
        compare.align_profiles(pl.scan_parquet(a).with_columns(pl.col("A").cast(pl.Float64)),pl.scan_parquet(a))
    with pytest.raises(ValueError,match="positive"):
        compare.align_profiles(a,b,min_cov=0)
