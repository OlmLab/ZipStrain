import pathlib

import pytest

from zipstrain import mapping


def _write_reads_table(tmp_path, rows, header="sample_name,reads1,reads2"):
    path = tmp_path / "reads.csv"
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return path


def test_read_reads_table_paired_and_single(tmp_path):
    table = _write_reads_table(
        tmp_path,
        ["s1,/data/s1_R1.fq.gz,/data/s1_R2.fq.gz", "s2,/data/s2.fq.gz,"],
    )
    samples = mapping.read_reads_table(table)
    assert [s.sample_name for s in samples] == ["s1", "s2"]
    assert samples[0].is_paired
    assert samples[0].reads2 == pathlib.Path("/data/s1_R2.fq.gz")
    assert not samples[1].is_paired
    assert samples[1].reads2 is None


def test_read_reads_table_single_column(tmp_path):
    table = _write_reads_table(tmp_path, ["s1,/data/s1.fq.gz"], header="sample_name,reads1")
    samples = mapping.read_reads_table(table)
    assert len(samples) == 1
    assert not samples[0].is_paired


def test_read_reads_table_requires_columns(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("name,file\ns1,/x\n")
    with pytest.raises(ValueError, match="sample_name"):
        mapping.read_reads_table(bad)


def test_read_reads_table_rejects_incomplete_rows(tmp_path):
    table = _write_reads_table(tmp_path, ["s1,"], header="sample_name,reads1")
    with pytest.raises(ValueError, match="reads1"):
        mapping.read_reads_table(table)


def test_merge_sylph_abundances_keeps_single_header(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    a.write_text("header\ta\nrow1\t1\n")
    b.write_text("header\ta\nrow2\t2\n")
    merged = mapping.merge_sylph_abundances([a, b], tmp_path / "merged.tsv")
    lines = merged.read_text().splitlines()
    assert lines == ["header\ta", "row1\t1", "row2\t2"]


def test_require_tool_missing_raises():
    with pytest.raises(RuntimeError, match="not found on PATH"):
        mapping.require_tool("definitely-not-a-real-tool-xyz")


def test_run_map_requires_stb_with_reference(tmp_path):
    table = _write_reads_table(tmp_path, ["s1,/data/s1.fq.gz"], header="sample_name,reads1")
    with pytest.raises(ValueError, match="--stb-file is required"):
        mapping.run_map(
            reads_table=table,
            output_dir=tmp_path / "out",
            reference_fasta=tmp_path / "ref.fna",
            stb_file=None,
        )


def test_run_map_requires_cache_and_db_without_reference(tmp_path):
    table = _write_reads_table(tmp_path, ["s1,/data/s1.fq.gz"], header="sample_name,reads1")
    with pytest.raises(ValueError, match="--genome-cache-dir is required"):
        mapping.run_map(reads_table=table, output_dir=tmp_path / "out")

    with pytest.raises(ValueError, match="--sylph-db is required"):
        mapping.run_map(
            reads_table=table,
            output_dir=tmp_path / "out",
            genome_cache_dir=tmp_path / "cache",
        )


def test_run_map_provided_reference_orchestration(tmp_path, monkeypatch):
    """Exercise the full orchestration for the provided-reference route with the
    external-tool calls mocked, so we validate wiring and the samples table."""
    table = _write_reads_table(
        tmp_path,
        ["s1,/data/s1_R1.fq.gz,/data/s1_R2.fq.gz", "s2,/data/s2.fq.gz,"],
    )
    ref = tmp_path / "ref.fna"
    ref.write_text(">chr1\nACGT\n")
    stb = tmp_path / "ref.stb"
    stb.write_text("chr1\tgenome1\n")

    calls = {"index": 0, "mapped": [], "index_prefix": None}

    def _fake_index(*, reference_fasta, index_prefix, threads):
        calls["index"] += 1
        calls["index_prefix"] = index_prefix
        return index_prefix

    def _fake_map(*, sample, index_prefix, output_bam, threads, non_competitive=False):
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        output_bam.write_text("bam")
        calls["mapped"].append(sample.sample_name)
        return output_bam

    monkeypatch.setattr(mapping, "build_bowtie2_index", _fake_index)
    monkeypatch.setattr(mapping, "map_sample", _fake_map)

    steps: list[str] = []
    results = mapping.run_map(
        reads_table=table,
        output_dir=tmp_path / "out",
        reference_fasta=ref,
        stb_file=stb,
        threads=2,
        progress_callback=steps.append,
    )

    assert calls["index"] == 1
    assert calls["mapped"] == ["s1", "s2"]
    assert results["reference_fasta"] == ref
    assert results["stb_file"] == stb

    # Progress callback was told about each stage in order.
    joined = " | ".join(steps)
    assert "Building Bowtie2 index" in joined
    assert "Mapping s1 (1/2)" in joined
    assert "Mapping s2 (2/2)" in joined

    # Bowtie2 index lands in intermediate_files/bt2, not the working directory.
    assert calls["index_prefix"].parent == tmp_path / "out" / "intermediate_files" / "bt2"

    samples_txt = (tmp_path / "out" / "samples.txt").read_text().splitlines()
    assert samples_txt[0] == "sample_name,bamfile"
    assert samples_txt[1].startswith("s1,")
    assert samples_txt[2].startswith("s2,")


def test_ensure_sylph_db_reuses_existing(tmp_path):
    db = tmp_path / "existing.syldb"
    db.write_text("already here")
    # Should return immediately without attempting any download.
    assert mapping.ensure_sylph_db(db, url="http://invalid.invalid/x") == db


def test_run_map_sylph_empty_reference_raises_clear_error(tmp_path, monkeypatch):
    """When Sylph detects nothing, the built reference is empty; map should say so."""
    table = _write_reads_table(tmp_path, ["s1,/data/s1.fq.gz"], header="sample_name,reads1")
    db = tmp_path / "db.syldb"
    db.write_text("db")

    monkeypatch.setattr(mapping, "ensure_sylph_db", lambda sylph_db, url=None: pathlib.Path(sylph_db))

    def _fake_sylph_profile(*, sylph_db, sample, output_tsv, threads):
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        output_tsv.write_text("header\n")  # header only -> nothing detected
        return output_tsv

    monkeypatch.setattr(mapping, "run_sylph_profile", _fake_sylph_profile)

    def _fake_build_reference(**kwargs):
        # Simulate build-genome-db producing an empty reference (no genomes).
        out = pathlib.Path(kwargs["output_dir"]) / mapping.REFERENCE_FASTA_NAME
        out.write_text("")
        stb = pathlib.Path(kwargs["output_dir"]) / mapping.REFERENCE_STB_NAME
        stb.write_text("")
        return out, stb, None, None, {}

    monkeypatch.setattr(mapping.bdb, "build_reference_from_abundance", _fake_build_reference)

    with pytest.raises(RuntimeError, match="did not detect any reference genomes"):
        mapping.run_map(
            reads_table=table,
            output_dir=tmp_path / "out",
            sylph_db=db,
            genome_cache_dir=tmp_path / "cache",
        )


def test_write_reference_taxonomy_maps_accessions(tmp_path, monkeypatch):
    stb = tmp_path / "ref.stb"
    stb.write_text("scaffold_a\tGCF_000001.1\nscaffold_b\tGCA_000002.1\n")

    gtdb = tmp_path / "bac120.tsv"
    gtdb.write_text(
        "RS_GCF_000001.1\td__Bacteria;s__Foo bar\n"
        "GB_GCA_000002.1\td__Bacteria;s__Baz qux\n"
        "RS_GCF_999999.9\td__Bacteria;s__Not wanted\n"
    )
    monkeypatch.setattr(mapping, "ensure_gtdb_taxonomy", lambda cache_dir, urls=None: [gtdb])

    out = mapping.write_reference_taxonomy(
        stb_file=stb, cache_dir=tmp_path / "cache", output_file=tmp_path / "tax.tsv"
    )
    import polars as pl

    d = pl.read_csv(out, separator="\t")
    got = {r["genome"]: r["genome_taxonomy"] for r in d.iter_rows(named=True)}
    assert got == {"GCF_000001.1": "d__Bacteria;s__Foo bar", "GCA_000002.1": "d__Bacteria;s__Baz qux"}


def test_write_reference_taxonomy_returns_none_when_no_match(tmp_path, monkeypatch):
    stb = tmp_path / "ref.stb"
    stb.write_text("scaffold_a\tlocal_genome_1\n")
    gtdb = tmp_path / "bac120.tsv"
    gtdb.write_text("RS_GCF_000001.1\td__Bacteria;s__Foo bar\n")
    monkeypatch.setattr(mapping, "ensure_gtdb_taxonomy", lambda cache_dir, urls=None: [gtdb])
    assert mapping.write_reference_taxonomy(stb_file=stb, cache_dir=tmp_path, output_file=tmp_path / "t.tsv") is None
