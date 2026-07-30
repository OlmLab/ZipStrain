from pathlib import Path

import polars as pl
import pytest

from zipstrain import database
from zipstrain import utils


def _write_profile(path: Path, *, reference_hash: str = "ref", gene_hash: str = "gene") -> None:
    pl.DataFrame(
        {
            "chrom": ["chr1"],
            "genome": ["genome1"],
            "pos": [1],
            "gene": ["NA"],
            "A": [1],
            "T": [0],
            "C": [0],
            "G": [0],
        }
    ).write_parquet(
        path,
        metadata={
            utils.PROFILE_REFERENCE_HASH_METADATA_KEY: reference_hash,
            utils.PROFILE_GENE_HASH_METADATA_KEY: gene_hash,
        },
    )


def test_profile_database_empty_schema():
    prof_db = database.ProfileDatabase()
    collected = prof_db.db.collect()
    assert collected.columns == ["profile_name", "profile_location"]
    assert collected.height == 0


def test_profile_database_add_profile(tmp_path):
    profile_path = tmp_path / "sample.parquet"
    _write_profile(profile_path)

    prof_db = database.ProfileDatabase()
    prof_db.add_profile(
        {
            "profile_name": "sample_a",
            "profile_location": str(profile_path),
        }
    )
    collected = prof_db.db.collect()
    assert collected.rows() == [("sample_a", str(profile_path))]


def test_profile_database_from_csv_checks_contracts(tmp_path):
    profile_a = tmp_path / "a.parquet"
    profile_b = tmp_path / "b.parquet"
    _write_profile(profile_a, reference_hash="ref_a")
    _write_profile(profile_b, reference_hash="ref_b")

    csv_path = tmp_path / "profiles.csv"
    pl.DataFrame(
        {
            "profile_name": ["a", "b"],
            "profile_location": [str(profile_a), str(profile_b)],
        }
    ).write_csv(csv_path)

    with pytest.raises(ValueError, match="Profile contract metadata do not match"):
        database.ProfileDatabase.from_csv(csv_path)

    prof_db = database.ProfileDatabase.from_csv(csv_path, allow_mismatch=True)
    assert prof_db.db.collect().height == 2


def test_genome_compare_config_scope_compatibility():
    config_all = database.GenomeComparisonConfig(scope="all", min_cov=5, min_gene_compare_len=100)
    config_subset = database.GenomeComparisonConfig(scope="genome1,genome2", min_cov=5, min_gene_compare_len=100)
    config_other = database.GenomeComparisonConfig(scope="genome3", min_cov=5, min_gene_compare_len=100)

    assert config_all.is_compatible(config_subset)
    assert not config_subset.is_compatible(config_other)


def test_gene_compare_config_scope_compatibility():
    config_all = database.GeneComparisonConfig(scope="all:all", min_cov=5, min_gene_compare_len=100)
    config_subset = database.GeneComparisonConfig(scope="genome1:gene1", min_cov=5, min_gene_compare_len=100)
    config_other = database.GeneComparisonConfig(scope="genome2:gene1", min_cov=5, min_gene_compare_len=100)

    assert config_all.is_compatible(config_subset)
    assert not config_subset.is_compatible(config_other)


def test_genome_comparison_database_remaining_pairs(tmp_path):
    profile_a = tmp_path / "a.parquet"
    profile_b = tmp_path / "b.parquet"
    profile_c = tmp_path / "c.parquet"
    _write_profile(profile_a)
    _write_profile(profile_b)
    _write_profile(profile_c)

    profile_db_path = tmp_path / "profile_db.parquet"
    pl.DataFrame(
        {
            "profile_name": ["a", "b", "c"],
            "profile_location": [str(profile_a), str(profile_b), str(profile_c)],
        }
    ).write_parquet(profile_db_path)

    compare_path = tmp_path / "compare.parquet"
    pl.DataFrame(
        {
            "genome": ["genome1"],
            "total_positions": [10],
            "share_allele_pos": [9],
            "genome_ani": [90.0],
            "max_consecutive_length": [5],
            "sample_1": ["a"],
            "sample_2": ["b"],
        }
    ).write_parquet(compare_path)

    comp_db = database.GenomeComparisonDatabase(
        profile_db=database.ProfileDatabase(str(profile_db_path)),
        config=database.GenomeComparisonConfig(scope="all", min_cov=5, min_gene_compare_len=100),
        comp_db_loc=str(compare_path),
    )
    remaining = comp_db.to_complete_input_table().collect().sort(["sample_name_1", "sample_name_2"])
    assert remaining.select(["sample_name_1", "sample_name_2"]).rows() == [("a", "c"), ("b", "c")]


def test_gene_comparison_database_remaining_pairs(tmp_path):
    profile_a = tmp_path / "a.parquet"
    profile_b = tmp_path / "b.parquet"
    profile_c = tmp_path / "c.parquet"
    _write_profile(profile_a)
    _write_profile(profile_b)
    _write_profile(profile_c)

    profile_db_path = tmp_path / "profile_db.parquet"
    pl.DataFrame(
        {
            "profile_name": ["a", "b", "c"],
            "profile_location": [str(profile_a), str(profile_b), str(profile_c)],
        }
    ).write_parquet(profile_db_path)

    compare_path = tmp_path / "gene_compare.parquet"
    pl.DataFrame(
        {
            "genome": ["genome1"],
            "gene": ["gene1"],
            "total_positions": [10],
            "share_allele_pos": [9],
            "ani": [90.0],
            "sample_1": ["a"],
            "sample_2": ["b"],
        }
    ).write_parquet(compare_path)

    comp_db = database.GeneComparisonDatabase(
        profile_db=database.ProfileDatabase(str(profile_db_path)),
        config=database.GeneComparisonConfig(scope="all:all", min_cov=5, min_gene_compare_len=100),
        comp_db_loc=str(compare_path),
    )
    remaining = comp_db.to_complete_input_table().collect().sort(["sample_name_1", "sample_name_2"])
    assert remaining.select(["sample_name_1", "sample_name_2"]).rows() == [("a", "c"), ("b", "c")]
