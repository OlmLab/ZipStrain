import pathlib

import polars as pl

from zipstrain import build_db


def test_extract_sylph_accessions():
    df = pl.DataFrame(
        {
            "Genome_file": [
                "/ref/GTDB/GCF_000001405.40_genomic.fna.gz",
                "/ref/GTDB/GCA_123456.1_genomic.fna.gz",
                "no_accession_here",
            ],
            "other": ["foo", "bar", "GCF_000009999.2 random_label"],
        }
    )
    out = build_db.get_genome_accessions_from_sylph_output(df).collect()
    assert set(out["accession"].to_list()) == {
        "GCF_000001405.40",
        "GCA_123456.1",
        "GCF_000009999.2",
    }
    assert set(out.columns) == {"accession", "genome_name", "download_url"}


def test_local_genome_db_add_update_sync(tmp_path):
    db_path = tmp_path / ".genome_db.parquet"
    local_db = build_db.LocalGenomeDB(db_path)
    genomes = pl.DataFrame(
        {
            "accession": ["GCF_000001405.40", "GCA_123456.1"],
            "genome_name": ["human_ref", "test_ref"],
            "download_url": [None, None],
        }
    )
    local_db.add_or_update_genomes(genomes, source_tool="sylph")
    missing = local_db.missing()
    assert missing.height == 2

    genome_file = tmp_path / "genomes" / "GCF_000001405.40.fna"
    genome_file.parent.mkdir(parents=True, exist_ok=True)
    genome_file.write_text(">chr1\nATCG\n")
    local_db.set_location("GCF_000001405.40", genome_file)
    local_db.sync()

    rows = local_db.db.rows_by_key("accession", unique=True, named=True)
    assert rows["GCF_000001405.40"]["exists"] is True
    assert rows["GCA_123456.1"]["exists"] is False


def test_build_local_genome_db_with_mock_download(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.6\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.4\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def fake_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    local_db, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=fake_downloader,
    )

    assert extracted.height == 2
    assert set(report["status"].to_list()) == {"downloaded"}
    assert local_db.db.filter(pl.col("exists")).height == 2
    assert (genomes_dir / "GCF_000001405.40.fna").exists()
    assert (genomes_dir / "GCA_123456.1.fna").exists()

