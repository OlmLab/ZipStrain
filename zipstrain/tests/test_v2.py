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

from .selection_oracle import allele_events, codon_groups


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
        "--preprofile-dir", str(root), "-o", str(root / "pair" / "comparisons.parquet")])
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
    assert "gene_id" not in genome.columns and "max_identical_stretch" not in genome.columns
    assert set(genome["genome"]) == {"genomeA", "genomeB"}
    assert genome.filter(pl.col("genome") == "genomeA")["total_positions"][0] == 21
    gene_only = compare.compare_profiles(a, b, preprofile_dir=root, calculate="gene_ani").collect()
    assert "genome_ani" not in gene_only.columns
    gene = compare.compare_profiles(a, b, preprofile_dir=root,
        calculate="genome_ani+max_identical_stretch+gene_ani").collect()
    assert gene.height == 3 and gene["gene_total_positions"].to_list() == [9]*3
    assert gene["total_positions"].to_list() == [21]*3
    assert gene["share_allele_pos"].to_list() == [20]*3
    assert gene["max_identical_stretch"].to_list() == [15]*3
    assert gene.select("sample_1", "sample_2", "genome", "gene_id").n_unique() == 3
    assert gene["gene_ani"].to_list() == pytest.approx([100*8/9]*3)
    with pytest.raises(ValueError, match="requires"):
        compare.compare_profiles(a, b, calculate="gene_ani")
    with pytest.raises(ValueError, match="requires"):
        compare.compare_profiles(a, b, preprofile_dir=root, calculate="gene_selection")
    frame = pl.read_parquet(b)
    pl.concat([frame, frame.head(1)]).write_parquet(b, metadata=utils._read_custom_parquet_metadata(catalog))
    with pytest.raises(ValueError, match="Duplicate"):
        compare.align_profiles(a, b)


def write_two_scaffold_sample(root, catalog, name, chr1, *, missing=()):
    """Write chr1 (genomeA) and chr2 (genomeB) rows, leaving out the ``missing`` chr1 positions."""
    rows = [(chrom, genome, i+1, *[10 if base == b else 0 for b in "ACGT"], {"A":1,"C":2,"G":4,"T":8}[base])
            for chrom, genome, sequence in (("chr1", "genomeA", chr1), ("chr2", "genomeB", "GGAGGAGGA"))
            for i, base in enumerate(sequence) if (chrom, i+1) not in {("chr1", p) for p in missing}]
    schema = {"chrom":pl.String,"genome":pl.String,"pos":pl.Int64,"A":pl.Int64,"C":pl.Int64,"G":pl.Int64,"T":pl.Int64,"ref_base_bitmask":pl.UInt8}
    path = root / f"{name}_profile.parquet"
    metadata = utils._read_custom_parquet_metadata(catalog)
    metadata[profile.PROFILE_SORTED_METADATA_KEY] = profile.PROFILE_SORTED_METADATA_VALUE
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(path, metadata=metadata)
    return path


def test_identical_stretches_report_coordinates_and_match_max_stretch(assets):
    root, _, _, _, catalog = assets
    a = write_two_scaffold_sample(root, catalog, "a", "AAAGGATTTCCCAAAGGTTAA")
    # A mismatch at chr1:6 and no coverage at chr1:12 split chr1 into 1-5, 7-11 and 13-21.
    b = write_two_scaffold_sample(root, catalog, "b", "AAAGGCTTTCCCAAAGGTTAA", missing=(12,))
    stretches = compare.identical_stretches(a, b, min_length=1).collect()
    assert stretches.columns == ["sample_1", "sample_2", "genome", "chrom", "stretch_id", "start", "end", "length",
                                 "background_divergence", "expected_runs", "p_value"]
    # genomeA: 1 mismatch in 20 shared positions, so 2 runs and q = 0.95 ** length.
    nine = stretches.filter(pl.col("length") == 9, pl.col("genome") == "genomeA").row(0, named=True)
    assert nine["background_divergence"] == pytest.approx(0.05)
    assert nine["expected_runs"] == pytest.approx(2 * 0.95 ** 9)
    assert nine["p_value"] == pytest.approx(1 - (1 - 0.95 ** 9) ** 2)
    # genomeB matches everywhere, so no run length is surprising.
    assert stretches.filter(pl.col("genome") == "genomeB")["p_value"].to_list() == [1.0]
    assert stretches.select("genome", "chrom", "start", "end", "length").rows() == [
        ("genomeA", "chr1", 1, 5, 5), ("genomeA", "chr1", 7, 11, 5),
        ("genomeA", "chr1", 13, 21, 9), ("genomeB", "chr2", 1, 9, 9)]
    assert stretches["stretch_id"].is_unique().all() and stretches["stretch_id"].is_sorted()
    assert set(stretches["sample_1"]) == {"a"} and set(stretches["sample_2"]) == {"b"}
    longest = compare.compare_profiles(a, b, calculate="max_identical_stretch").collect()
    assert dict(longest.select("genome", "max_identical_stretch").rows()) == dict(
        stretches.group_by("genome").agg(pl.col("length").max()).rows())
    long_only = compare.identical_stretches(a, b, min_length=6, sample_1="x", sample_2="y").collect()
    assert long_only.select("sample_1", "sample_2", "chrom", "start", "end").rows() == [
        ("x", "y", "chr1", 13, 21), ("x", "y", "chr2", 1, 9)]
    assert long_only["stretch_id"].to_list() == stretches["stretch_id"].to_list()[2:]
    assert compare.identical_stretches(a, b, min_length=50).collect().is_empty()
    with pytest.raises(ValueError, match="positive"):
        compare.identical_stretches(a, b, min_length=0)


def _planted_matches(rng, d=1e-3):
    """A genome with random background mismatches, a planted identical and a planted divergent block."""
    import numpy as np

    frames = []
    for chrom, length in (("c1", 300_000), ("c2", 100_000)):
        pos = np.arange(1, length + 1)
        mismatch = rng.random(length) < d
        if chrom == "c1":
            mismatch[(pos >= 100_001) & (pos <= 160_000)] = False                 # recombined: identical
            divergent = (pos >= 220_001) & (pos <= 230_000)
            mismatch[divergent] = rng.random(int(divergent.sum())) < 0.03           # imported: 30x denser
            keep = (pos < 50_001) | (pos > 51_000)                                  # a coverage gap
            pos, mismatch = pos[keep], mismatch[keep]
        frames.append(pl.DataFrame({"chrom": chrom, "genome": "g1", "pos": pos, "match": ~mismatch}))
    return pl.concat(frames)


def test_mosaic_regions_find_planted_blocks_and_nothing_else():
    import numpy as np

    matches = _planted_matches(np.random.default_rng(5))
    events = compare._mosaic_events(matches.lazy()).collect()
    assert events.filter(pl.col("kind") == 1).height == (~matches["match"]).sum()
    blocks = compare._decode_mosaics(events, factor=20, switch_probability=1e-3)
    assert blocks["state"].to_list() == ["low_divergence", "high_divergence"]
    low, high = blocks.rows(named=True)
    assert low["chrom"] == high["chrom"] == "c1"
    assert abs(low["start"] - 100_001) < 3_000 and abs(low["end"] - 160_000) < 3_000 and low["mismatches"] <= 2
    assert abs(high["start"] - 220_001) < 500 and abs(high["end"] - 230_000) < 500
    assert high["divergence"] == pytest.approx(0.03, rel=0.25)
    assert blocks["background_divergence"][0] == pytest.approx(1e-3, rel=0.15)
    assert (blocks["log_likelihood_ratio"] > 0).all()
    # Coverage gaps are skipped, not counted: lengths are in shared positions.
    assert low["callable_positions"] == matches.filter(
        pl.col("chrom") == "c1", pl.col("pos").is_between(low["start"], low["end"])).height
    # Without mismatches nothing is reported, and the background alone yields no blocks.
    assert compare._decode_mosaics(events.filter(pl.col("kind") != 1), factor=20, switch_probability=1e-3).is_empty()
    plain = matches.filter(pl.col("chrom") == "c2")
    assert compare._decode_mosaics(compare._mosaic_events(plain.lazy()).collect(),
                                   factor=20, switch_probability=1e-3).is_empty()


def test_mosaic_regions_through_the_pair_api(assets):
    root, _, _, _, catalog = assets
    a = write_two_scaffold_sample(root, catalog, "a", "AAAGGATTTCCCAAAGGTTAA")
    b = write_two_scaffold_sample(root, catalog, "b", "AAAGGCTTTCCCAAAGGTTAA", missing=(12,))
    pair = compare.PairComparison(a, b)
    blocks, stretches = pair.collect(pair.mosaic_regions(), pair.identical_stretches(min_length=1))
    assert blocks.columns == ["sample_1", "sample_2", *compare._MOSAIC_SCHEMA]
    assert_frame_equal(blocks, compare.mosaic_regions(a, b).collect())
    assert stretches.height == 4
    with pytest.raises(ValueError, match="factor"):
        pair.mosaic_regions(factor=1)
    with pytest.raises(ValueError, match="switch_probability"):
        pair.mosaic_regions(switch_probability=0)


def test_pair_comparison_reuses_one_join_for_every_analysis(assets):
    root, _, _, _, catalog = assets
    a = write_two_scaffold_sample(root, catalog, "a", "AAAGGATTTCCCAAAGGTTAA")
    b = write_two_scaffold_sample(root, catalog, "b", "AAAGGCTTTCCCAAAGGTTAA", missing=(12,))
    pair = compare.PairComparison(a, b, preprofile_dir=root)
    assert pair.names == ("a", "b")
    assert_frame_equal(pair.aligned.collect(), compare.align_profiles(a, b).collect())
    counts, stretches = pair.collect(
        pair.matches.group_by("genome").agg(pl.len(), pl.col("match").sum()).sort("genome"),
        pair.identical_stretches(min_length=6),
    )
    assert counts.rows() == [("genomeA", 20, 19), ("genomeB", 9, 9)]
    assert_frame_equal(stretches, compare.identical_stretches(a, b, min_length=6).collect())
    keys = ["sample_1", "sample_2", "genome", "gene_id"]
    for calculate in ("genome_ani+max_identical_stretch", "all"):
        assert_frame_equal(
            pair.compare(calculate).collect().sort(keys[:3] if calculate != "all" else keys),
            compare.compare_profiles(a, b, preprofile_dir=root, calculate=calculate).collect()
            .sort(keys[:3] if calculate != "all" else keys))
    with pytest.raises(ValueError, match="min_cov"):
        compare.PairComparison(a, b, min_cov=0)
    with pytest.raises(ValueError, match="Unsupported"):
        compare.PairComparison(a, b, ani_method="bogus")


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
    assert frame["max_identical_stretch"].to_list() == [21] * 3
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
              "--allele-integration", integration, "--calculate", f"{calculate}+gene_selection",
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
    category, events = allele_events(source, None, "A")
    assert category == "snv" and events == [("A", "G", 1.0)]
    assert allele_events(source, None, "A", "weighted")[1] == [("A", "G", .2), ("A", "C", .1)]
    assert allele_events(Counter({"G":5,"C":5}), None, "A")[1] == [("A", "G", .5), ("A", "C", .5)]
    assert allele_events(Counter({"G":10}), None, "A")[0] == "sns"
    assert allele_events(Counter({"A":10}), None, "A")[1] == []
    a, b = Counter({"A":8,"G":2}), Counter({"A":5,"C":5})
    assert allele_events(a,b,"A")[1] == [("A","C",1.0)]
    assert allele_events(a,b,"A","weighted")[1] == [("A","C",.4),("G","A",.1),("G","C",.1)]
    with pytest.raises(ValueError):
        allele_events(a,b,"A","unknown")


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
    groups = list(codon_groups(variants,catalog,min_cov=5))
    assert groups == [((gene_id,0),"AAA",Counter({"AAG":8,"AAC":2}))]
    assert list(codon_groups(variants,catalog,min_cov=100)) == []


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
        "--profile-2",str(b),"--preprofile-dir",str(root),"-o",str(pair / "comparisons.parquet"),
        "--codon-counts-1",str(counts),"--codon-counts-2",str(b_counts),"--calculate","gene_selection",
        "--allele-integration","weighted"])
    assert result.exit_code == 0,result.output+repr(result.exception)
    genes = pl.read_parquet(pair / "comparisons.parquet")
    assert all(name in genes.columns for name in gi.PAIR_SELECTION_COLUMNS)
    assert "dN_dS" not in genes.columns and "pN_pS" not in genes.columns
    assert genes.filter(pl.col("gene_id") == gene_id)["pS"][0] > 0
    inputs = root / "profiles.csv"
    inputs.write_text(f"profile_name,profile_location\na,{positions}\nb,{b}\n")
    run = root / "collection"
    result = CliRunner().invoke(cli.cli,["compare","--profile-db",str(inputs),"--run-dir",str(run),
        "--preprofile-dir",str(root),"--calculate","gene_selection","--min-gene-compare-len","1","--no-csv"])
    assert result.exit_code == 0,result.output+repr(result.exception)
    assert pl.read_parquet(run / "all_comparisons.parquet").height == 3


def test_new_function_docstrings_and_guide_entries():
    import ast
    root = Path(__file__).resolve().parents[2]
    guide = (root / "docs" / "v2.md").read_text()
    modules = {"gene_info":None, "codons":None, "pair_selection":None,
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
    assert row["max_identical_stretch"] == 4 and row["genome_ani"] == 80
    mismatch = source.with_columns(pl.lit(0).alias("A"),pl.lit(10).alias("T"))
    row = compare.compare_profiles(source,mismatch).collect().row(0,named=True)
    assert row["max_identical_stretch"] == 0 and row["genome_ani"] == 0
    no_overlap = source.with_columns((pl.col("pos")+100).alias("pos"))
    row = compare.compare_profiles(source,no_overlap).collect().row(0,named=True)
    assert row["total_positions"] == row["max_identical_stretch"] == 0 and row["genome_ani"] is None


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


def _random_codon_table(root, catalog, name, rng, *, stop_rate):
    """Write random phased evidence with reference neighbours, stops, and multi-codon sites."""
    references = pl.read_parquet(catalog.with_name("reference_codons.parquet"))
    stops = [codon for codon, acid in codons.GENETIC_CODE.items() if acid == "*"]
    everything = list(codons.GENETIC_CODE)
    rows = {}
    for gene_id, index, reference in references.iter_rows():
        for _ in range(int(rng.integers(1, 7))):
            draw = rng.random()
            if draw < stop_rate:
                codon = stops[int(rng.integers(len(stops)))]
            elif draw < .7:
                offset = int(rng.integers(3))
                codon = reference[:offset]+"ACGT"[int(rng.integers(4))]+reference[offset+1:]
            else:
                codon = everything[int(rng.integers(64))]
            rows[gene_id, index, codon] = rows.get((gene_id, index, codon), 0) + int(rng.integers(1, 40))
    return write_counts(root, catalog, name, [(*key, count) for key, count in rows.items()])


@pytest.mark.parametrize("integration", ["consensus", "weighted"])
@pytest.mark.parametrize(("stop_rate", "min_freq"), [(0.0, 0.0), (0.25, 0.0), (0.05, 0.3)])
def test_selection_kernel_matches_python_oracle(assets, integration, stop_rate, min_freq):
    import numpy as np
    from . import selection_oracle as oracle

    root, _, _, _, catalog = assets
    rng = np.random.default_rng(7)
    options = dict(allele_integration=integration, min_freq=min_freq, min_cov=2)
    a = _random_codon_table(root, catalog, "a", rng, stop_rate=stop_rate)
    b = _random_codon_table(root, catalog, "b", rng, stop_rate=stop_rate)
    for kernel, expected in (
        (codons.pairwise_selection(a, b, catalog, **options), oracle.pairwise_selection(a, b, catalog, **options)),
        (codons.reference_selection(a, catalog, **options), oracle.reference_selection(a, catalog, **options)),
    ):
        assert_frame_equal(kernel.collect().sort("gene_id"), expected.collect().sort("gene_id"), abs_tol=1e-9)


def test_selection_results_do_not_depend_on_batch_size(assets):
    import numpy as np

    root, _, _, _, catalog = assets
    rng = np.random.default_rng(11)
    a = _random_codon_table(root, catalog, "a", rng, stop_rate=.05)
    b = _random_codon_table(root, catalog, "b", rng, stop_rate=.05)
    for integration in codons.ALLELE_INTEGRATION_MODES:
        results = [codons.pairwise_selection(a, b, catalog, allele_integration=integration, batch_size=size).collect()
                   for size in (1, 2, 3, 1024)]
        assert all(result.equals(results[0]) for result in results[1:])
    with pytest.raises(ValueError, match="batch size"):
        codons.pairwise_selection(a, b, catalog, batch_size=0)


def test_selection_results_do_not_depend_on_gene_range_size(assets, monkeypatch):
    """Packing and matching in gene-ID ranges must give the single-range answer."""
    import numpy as np
    import pyarrow.parquet as pq

    root, _, _, _, catalog = assets
    rng = np.random.default_rng(13)
    a = _random_codon_table(root, catalog, "a", rng, stop_rate=.05)
    b = _random_codon_table(root, catalog, "b", rng, stop_rate=.05)
    calls = {"pair": lambda: codons.pairwise_selection(a, b, catalog, min_cov=2).collect().sort("gene_id"),
             "reference": lambda: codons.reference_selection(a, catalog, min_cov=2).collect().sort("gene_id")}
    expected = {name: call() for name, call in calls.items()}
    assert expected["pair"].height > 1
    real_writer = pq.ParquetWriter

    class TinyRowGroups(real_writer):  # split genes across row groups, so cuts must not split them
        def write_table(self, table, row_group_size=None):
            super().write_table(table, row_group_size=3)

    monkeypatch.setattr(pq, "ParquetWriter", TinyRowGroups)
    for chunk in (1, 2, 5, 64):
        monkeypatch.setattr(codons, "SELECTION_CHUNK_CODONS", chunk)
        for name, call in calls.items():
            assert_frame_equal(call(), expected[name]), (name, chunk)


def test_shared_position_scan_is_independent_of_chunks_and_matches_brute_force(assets):
    """Stretches and gene counts must survive batch edges, scaffold changes, and coordinate gaps."""
    import numpy as np

    rng = np.random.default_rng(3)
    rows = []
    for chrom, genome, length in (("c1", "g1", 700), ("c2", "g1", 300), ("c3", "g2", 500)):
        for pos in np.flatnonzero(rng.random(length) < .96) + 1:
            rows.append((chrom, genome, int(pos), bool(rng.random() < .93)))
    shared = pl.DataFrame(rows, schema={"chrom": pl.String, "genome": pl.String, "pos": pl.Int64, "match": pl.Boolean},
                          orient="row")
    genes = pl.DataFrame(
        {"gene_id": [1, 2, 3, 4], "genome": ["g1", "g1", "g1", "g2"], "chrom": ["c1", "c1", "c2", "c3"],
         "start": [10, 50, 1, 100], "end": [120, 400, 300, 90_000]},
        schema_overrides={"gene_id": pl.UInt64},
    )

    expected_genome = {}
    for genome in ("g1", "g2"):
        total = share = best = 0
        for chrom in shared.filter(pl.col("genome") == genome)["chrom"].unique():
            run, previous = 0, None
            for pos, match in shared.filter(pl.col("chrom") == chrom).select("pos", "match").iter_rows():
                total += 1
                share += match
                run = (run + 1 if previous == pos - 1 else 1) if match else 0
                best, previous = max(best, run), pos
        expected_genome[genome] = (total, share, best)
    expected_gene = {}
    for gene_id, _, chrom, start, end in genes.iter_rows():
        inside = shared.filter((pl.col("chrom") == chrom) & pl.col("pos").is_between(start, end))
        expected_gene[gene_id] = (inside.height, int(inside["match"].sum()))

    for chunk_size in (1, 2, 3, 17, 256, 1 << 20):
        genome, gene = compare._scan_shared_positions(shared.lazy(), genes, stretches=True, chunk_size=chunk_size)
        assert {row[0]: row[1:] for row in genome.iter_rows()} == expected_genome, chunk_size
        assert {row[0]: (row[2], row[3]) for row in gene.iter_rows()} == expected_gene, chunk_size
    longest = compare._identical_stretch_rows(shared.lazy()).group_by("genome").agg(pl.col("length").max())
    assert dict(longest.collect().rows()) == {genome: best for genome, (_, _, best) in expected_genome.items()}


def test_collection_validates_each_sample_once_and_parallel_pairs_match(assets, monkeypatch):
    root, _, _, _, catalog = assets
    profiles = []
    for name, sequence in (("a", "AAAGGATTTCCCAAAGGTTAA"), ("b", "AAAGGTTTTCCCAAAGGTTAA"), ("c", "AAAGGATTTCCTAAAGGTTAA")):
        path = write_sample(root, catalog, name, sequence)
        gene_id = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
        write_counts(root, catalog, name, [(gene_id, 0, sequence[:3], 10), (gene_id, 1, sequence[3:6], 10)], source_profile=path)
        profiles.append((name, path))

    real, calls = codons.read_codon_counts, []
    monkeypatch.setattr(codons, "read_codon_counts", lambda path, info: calls.append(path) or real(path, info))
    keys = ["sample_1", "sample_2", "gene_id"]
    serial = pl.read_parquet(compare.compare_profile_collection(profiles, root/"serial", preprofile_dir=root, workers=1)).sort(keys)
    assert len(calls) == 3
    # Spawned workers re-import the test module, so this exercises the real process path.
    parallel = pl.read_parquet(compare.compare_profile_collection(profiles, root/"parallel", preprofile_dir=root, workers=2)).sort(keys)
    assert serial.equals(parallel)
    assert not list((root/"serial").glob("intermediate_files/pair_selection/sample-*"))
    assert len(list((root/"serial").glob("intermediate_files/pair_selection/pair-*"))) == 3
    assert serial.height > 0 and serial["dN"].null_count() < serial.height
    for pair in ("a", "b"), ("a", "c"), ("b", "c"):
        single = compare.compare_profiles(root/f"{pair[0]}_profile.parquet", root/f"{pair[1]}_profile.parquet", preprofile_dir=root,
            codon_counts_1=root/f"{pair[0]}_codon_counts.parquet", codon_counts_2=root/f"{pair[1]}_codon_counts.parquet",
            sample_1=pair[0], sample_2=pair[1]).collect().sort(keys)
        assert_frame_equal(serial.filter((pl.col("sample_1") == pair[0]) & (pl.col("sample_2") == pair[1])), single)
    with pytest.raises(ValueError, match="workers"):
        compare.compare_profile_collection(profiles, root/"bad", preprofile_dir=root, workers=0, calculate="gene_ani")


def test_pair_workers_respect_cpu_and_memory(monkeypatch):
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "8")
    assert compare._pair_workers("auto", 100, 16.0, 1_000_000) == (4, 2)
    assert compare._pair_workers("auto", 2, 16.0, 30_000_000) == (2, 4)
    assert compare._pair_workers(16, 100, 16.0, 1000) == (16, 1)
    # Too little work to repay starting processes.
    assert compare._pair_workers("auto", 3, 16.0, 1000) == (1, 8)
    # A pair of 40 M rows holds about 5 GiB, so 16 GiB fits three.
    assert compare._pair_workers("auto", 100, 16.0, 40_000_000) == (3, 2)
    assert compare._pair_workers("auto", 100, 1.0, 40_000_000) == (1, 8)
    monkeypatch.setenv("ZIPSTRAIN_CPU_BUDGET", "1")
    assert compare._pair_workers("auto", 100, 16.0, 1_000_000) == (1, 1)
