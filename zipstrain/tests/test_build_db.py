import pathlib
import gzip
import threading
import time
from urllib.error import HTTPError

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
    }
    assert set(out.columns) == {"accession", "genome_name", "download_url"}


def test_extract_sylph_accessions_from_gtdb_path_layout():
    df = pl.DataFrame(
        {
            "Sample_file": ["sample1", "sample2"],
            "Genome_file": [
                "gtdb_genomes_reps_r220/database/GCA/949/068/525/GCA_949068525.1_genomic.fna.gz",
                "gtdb_genomes_reps_r226/database/GCF/000/001/405/GCF_000001405.40_genomic.fna.gz",
            ],
            "Eff_cov": [0.5, 1.1],
        }
    )
    out = build_db.get_genome_accessions_from_sylph_output(df).collect()
    assert set(out["accession"].to_list()) == {
        "GCA_949068525.1",
        "GCF_000001405.40",
    }


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


def test_build_local_genome_db_report_handles_mixed_success_failure(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.6\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.4\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.zip"

    def mixed_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        if accession == "GCF_000001405.40":
            out = destination_dir / f"{accession}.fna"
            out.write_text(f">{accession}\nATCG\n")
            return out
        raise ValueError("No FASTA files found in archive: fake.zip")

    local_db, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=mixed_downloader,
    )

    assert extracted.height == 2
    assert report.height == 2
    assert set(report["status"].to_list()) == {"downloaded", "failed"}
    error_map = {r["accession"]: r["error"] for r in report.to_dicts()}
    assert error_map["GCF_000001405.40"] is None
    assert "No FASTA files found in archive" in error_map["GCA_123456.1"]
    assert local_db.db.filter(pl.col("exists")).height == 1


def test_build_local_genome_db_retries_download_until_success(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    attempts = {"count": 0}

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def flaky_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        attempts["count"] += 1
        if attempts["count"] < 3:
            raise ConnectionError("temporary network failure")
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=flaky_downloader,
        max_download_attempts=3,
        backoff_base_seconds=0,
    )

    assert extracted.height == 1
    assert attempts["count"] == 3
    assert report.height == 1
    assert report["status"].to_list() == ["downloaded"]
    assert (genomes_dir / "GCF_000001405.40.fna").exists()


def test_build_local_genome_db_uses_parallel_download_workers(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
        "/ref/GCA_123456.1_genomic.fna.gz,1.0\n"
        "/ref/GCF_000009999.1_genomic.fna.gz,1.0\n"
        "/ref/GCA_999999.1_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    thread_ids: set[int] = set()
    lock = threading.Lock()

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def slow_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        with lock:
            thread_ids.add(threading.get_ident())
        time.sleep(0.05)
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=slow_downloader,
        download_workers=3,
        backoff_base_seconds=0,
    )

    assert extracted.height == 4
    assert report.height == 4
    assert set(report["status"].to_list()) == {"downloaded"}
    assert len(thread_ids) >= 2


def test_build_local_genome_db_progress_reaches_total_with_failures(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
        "/ref/GCA_123456.1_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    progress_events: list[tuple[str, int, int, int, str | None]] = []

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def mixed_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        if accession == "GCA_123456.1":
            raise OSError("network issue")
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    def progress_callback(event: str, completed: int, remaining: int, total: int, accession: str | None) -> None:
        progress_events.append((event, completed, remaining, total, accession))

    _, _, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=mixed_downloader,
        progress_callback=progress_callback,
        max_download_attempts=1,
        backoff_base_seconds=0,
        download_workers=2,
    )

    assert report.height == 2
    assert set(report["status"].to_list()) == {"downloaded", "failed"}
    assert progress_events[0][0] == "start"
    assert progress_events[0][1] == 0
    assert progress_events[-1][0] == "done"
    assert progress_events[-1][1] == 2
    assert progress_events[-1][2] == 0
    assert progress_events[-1][3] == 2


def test_build_local_genome_db_honors_retry_after_on_rate_limit(tmp_path, monkeypatch):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    attempts = {"count": 0}
    sleep_calls: list[float] = []

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.zip"

    def rate_limited_then_ok(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise HTTPError(
                url=f"https://example.test/{accession}.zip",
                code=429,
                msg="Too Many Requests",
                hdrs={"Retry-After": "2"},
                fp=None,
            )
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    monkeypatch.setattr(build_db.random, "uniform", lambda _a, _b: 0.0)
    monkeypatch.setattr(build_db.time, "sleep", lambda seconds: sleep_calls.append(seconds))

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=rate_limited_then_ok,
        max_download_attempts=3,
        backoff_base_seconds=0,
        download_workers=1,
    )

    assert extracted.height == 1
    assert report.height == 1
    assert report["status"].to_list() == ["downloaded"]
    assert attempts["count"] == 2
    assert len(sleep_calls) == 1
    assert sleep_calls[0] >= 15.0


def test_build_local_genome_db_does_not_retry_missing_accession_http_404(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    attempts = {"count": 0}

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.zip"

    def missing_accession(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = (url, accession, destination_dir)
        attempts["count"] += 1
        raise HTTPError(
            url=f"https://example.test/{accession}.zip",
            code=404,
            msg="Not Found",
            hdrs={},
            fp=None,
        )

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=missing_accession,
        max_download_attempts=5,
        backoff_base_seconds=0,
        download_workers=1,
    )

    assert extracted.height == 1
    assert attempts["count"] == 1
    assert report.height == 1
    assert report["status"].to_list() == ["failed"]
    assert "404" in report["error"].to_list()[0]
    assert "after 1 attempts" in report["error"].to_list()[0]


def test_build_local_genome_db_retries_server_http_503(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,abundance\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"
    attempts = {"count": 0}

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.zip"

    def busy_server_then_ok(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        attempts["count"] += 1
        if attempts["count"] < 3:
            raise HTTPError(
                url=f"https://example.test/{accession}.zip",
                code=503,
                msg="Service Unavailable",
                hdrs={},
                fp=None,
            )
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}\nATCG\n")
        return out

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=busy_server_then_ok,
        max_download_attempts=3,
        backoff_base_seconds=0,
        download_workers=1,
    )

    assert extracted.height == 1
    assert attempts["count"] == 3
    assert report.height == 1
    assert report["status"].to_list() == ["downloaded"]


def test_build_local_genome_db_filters_zero_abundance_rows(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,sample_a,sample_b\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.0,0.0\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.4,0.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def fake_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}_contig1\nATCG\n")
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

    assert extracted.height == 1
    assert extracted["accession"].to_list() == ["GCA_123456.1"]
    assert report.height == 1
    assert report["accession"].to_list() == ["GCA_123456.1"]
    assert local_db.db.filter(pl.col("exists")).height == 1


def test_build_local_genome_db_prefers_eff_cov_for_nonzero_filter(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,Eff_cov,Mean_cov\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.0,2.5\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.4,1.0\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def fake_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}_contig1\nATCG\n")
        return out

    _, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        url_resolver=fake_resolver,
        downloader=fake_downloader,
    )

    assert extracted.height == 1
    assert extracted["accession"].to_list() == ["GCA_123456.1"]
    assert report.height == 1


def test_build_local_genome_db_uses_local_source_paths_before_download(tmp_path):
    genomes_root = tmp_path / "gtdb_genomes_reps_r220" / "database" / "GCA" / "949" / "068" / "525"
    genomes_root.mkdir(parents=True, exist_ok=True)
    local_gz = genomes_root / "GCA_949068525.1_genomic.fna.gz"
    with gzip.open(local_gz, "wt") as handle:
        handle.write(">contig1\nATCG\n")

    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,Eff_cov\n"
        "gtdb_genomes_reps_r220/database/GCA/949/068/525/GCA_949068525.1_genomic.fna.gz,0.8\n"
    )
    db_path = tmp_path / ".genome_db.parquet"
    genomes_dir = tmp_path / "genomes"

    def failing_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = (url, accession, destination_dir)
        raise RuntimeError("downloader should not be called when local source exists")

    local_db, extracted, report = build_db.build_local_genome_db(
        tool_name="sylph",
        abundance_table=input_csv,
        db_path=db_path,
        genomes_dir=genomes_dir,
        download=True,
        downloader=failing_downloader,
    )

    assert extracted.height == 1
    assert extracted["accession"].to_list() == ["GCA_949068525.1"]
    assert report.height == 0
    cached_fasta = genomes_dir / "GCA_949068525.1.fna"
    assert cached_fasta.exists()
    assert "ATCG" in cached_fasta.read_text()
    db_row = local_db.db.filter(pl.col("accession") == "GCA_949068525.1").to_dicts()[0]
    assert db_row["exists"] is True


def test_build_reference_from_abundance_writes_concat_and_stb(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,sample_1,sample_2\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,0.0,1.2\n"
        "/ref/GCA_123456.1_genomic.fna.gz,0.2,0.0\n"
    )
    cache_dir = tmp_path / "cache"
    output_dir = tmp_path / "output"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def fake_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">contigA\nATCG\n>contigB\nGGTT\n")
        return out

    out_fasta, out_stb, extracted, report, summary = build_db.build_reference_from_abundance(
        tool_name="sylph",
        abundance_table=input_csv,
        cache_dir=cache_dir,
        output_dir=output_dir,
        url_resolver=fake_resolver,
        downloader=fake_downloader,
    )

    assert out_fasta.exists()
    assert out_stb.exists()
    assert extracted.height == 2
    assert report.height == 2
    assert summary["selected_genomes"] == 2
    assert summary["cached_before_download"] == 0
    assert summary["attempted_downloads"] == 2
    assert summary["downloaded_now"] == 2
    assert summary["failed_downloads"] == 0
    assert summary["cached_after_download"] == 2

    fasta_text = out_fasta.read_text()
    assert ">GCA_123456.1__contigA" in fasta_text
    assert ">GCF_000001405.40__contigB" in fasta_text

    stb = pl.read_csv(out_stb, separator="\t", has_header=False).rename(
        {"column_1": "scaffold", "column_2": "genome"}
    )
    assert stb.height == 4
    assert set(stb["genome"].to_list()) == {"GCA_123456.1", "GCF_000001405.40"}


def test_build_reference_from_abundance_skips_failed_after_retries(tmp_path):
    input_csv = tmp_path / "sylph_abundance.csv"
    input_csv.write_text(
        "Genome_file,sample_1\n"
        "/ref/GCF_000001405.40_genomic.fna.gz,1.0\n"
        "/ref/GCA_123456.1_genomic.fna.gz,1.0\n"
    )
    cache_dir = tmp_path / "cache"
    output_dir = tmp_path / "output"

    def fake_resolver(accession: str, _: dict) -> str:
        return f"https://example.test/{accession}.fna"

    def partly_failing_downloader(url: str, accession: str, destination_dir: pathlib.Path) -> pathlib.Path:
        _ = url
        destination_dir.mkdir(parents=True, exist_ok=True)
        if accession == "GCA_123456.1":
            raise ConnectionError("persistent download failure")
        out = destination_dir / f"{accession}.fna"
        out.write_text(f">{accession}_contig1\nATCG\n")
        return out

    out_fasta, out_stb, extracted, report, summary = build_db.build_reference_from_abundance(
        tool_name="sylph",
        abundance_table=input_csv,
        cache_dir=cache_dir,
        output_dir=output_dir,
        url_resolver=fake_resolver,
        downloader=partly_failing_downloader,
        max_download_attempts=3,
        backoff_base_seconds=0,
        continue_on_missing=True,
    )

    assert extracted.height == 2
    assert report.height == 2
    assert set(report["status"].to_list()) == {"downloaded", "failed"}
    assert summary["missing_after_retries"] == 1
    assert summary["cached_after_download"] == 1
    assert out_fasta.exists()
    assert out_stb.exists()
    fasta_text = out_fasta.read_text()
    assert "GCF_000001405.40" in fasta_text
    assert "GCA_123456.1" not in fasta_text
