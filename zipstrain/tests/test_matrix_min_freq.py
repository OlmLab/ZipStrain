"""Build-time allele-frequency filter (min_freq) for matrix stores."""

from pathlib import Path

from click.testing import CliRunner
import duckdb
import numpy as np
import polars as pl
import pytest

from zipstrain import cli
from zipstrain import matrix_pairs as mp
from zipstrain import matrix_workflow as mw

# position -> (A, T, C, G); axis index is pos - 1 for a single 0-based BED interval.
SAMPLE_1 = {
    1: (95, 5, 0, 0),   # minor allele at exactly 5%: kept at min_freq=0.05
    2: (96, 4, 0, 0),   # 4%: dropped at 0.05
    3: (50, 0, 50, 0),  # two major alleles
    4: (3, 1, 0, 0),    # depth 4 < min_cov 5: position stored as zero
    5: (99, 1, 0, 0),   # 1%: kept at 0.01, dropped at 0.05
    6: (0, 0, 0, 20),
}
SAMPLE_2 = {
    1: (0, 100, 0, 0),
    2: (0, 100, 0, 0),
    3: (0, 0, 100, 0),
    4: (0, 0, 0, 0),
    5: (0, 100, 0, 0),
    6: (0, 0, 0, 20),
}
A, T, C, G = mp.BITMASK_BASE_BITS


def _write_profile(path: Path, rows: dict[int, tuple[int, int, int, int]]) -> None:
    positions = [pos for pos, counts in rows.items() if sum(counts) > 0]
    pl.DataFrame(
        {
            "chrom": ["chr1"] * len(positions),
            "pos": positions,
            "gene": ["gene1"] * len(positions),
            "genome": ["genome1"] * len(positions),
            **{base: [rows[pos][i] for pos in positions] for i, base in enumerate("ATCG")},
        }
    ).write_parquet(path)


def _inputs(tmp_path: Path, samples: dict[str, dict]) -> tuple[Path, Path, Path]:
    profiles = tmp_path / "profiles"
    profiles.mkdir(parents=True, exist_ok=True)
    for name, rows in samples.items():
        _write_profile(profiles / f"{name}.parquet", rows)
    bed, stb = tmp_path / "ref.bed", tmp_path / "ref.stb"
    bed.write_text("chr1\t0\t6\n")
    stb.write_text("chr1\tgenome1\n")
    return profiles, bed, stb


def _build(tmp_path: Path, name: str, samples: dict[str, dict], **kwargs) -> Path:
    profiles, bed, stb = _inputs(tmp_path, samples)
    out = tmp_path / f"{name}.h5"
    mp.build_matrix_hdf5(profile_dir=profiles, output_file=out, bed_file=bed, stb_file=stb,
                         memory_limit_gb=1, **kwargs)
    return out


def _rows(store: Path) -> np.ndarray:
    metadata = mp._load_matrix_hdf5_metadata(store)
    spec = mp._load_matrix_hdf5_genomes(store, genome=None)[0]
    dataset = mp._Hdf5GenomeMatrixNumpyDataset(
        store, spec.genome_idx, spec.matrix_length,
        mp._matrix_dtype_name_from_metadata(metadata), metadata["matrix_value_semantics"],
    )
    try:
        return dataset._load_batch_numpy(np.arange(len(mp._load_matrix_hdf5_samples(store))))
    finally:
        dataset.close()


def test_default_is_five_percent():
    assert mp.MATRIX_BUILD_MIN_FREQ_DEFAULT == 0.05


@pytest.mark.parametrize(
    "min_freq,expected",
    [
        (None, [A | T, A, A | C, 0, A, G]),  # default 0.05
        (0.05, [A | T, A, A | C, 0, A, G]),
        (0.01, [A | T, A | T, A | C, 0, A | T, G]),
        (0.0, [A | T, A | T, A | C, 0, A | T, G]),
        (0.5, [A, A, A | C, 0, A, G]),
    ],
)
def test_bitmask_store_applies_min_freq(tmp_path, min_freq, expected):
    kwargs = {} if min_freq is None else {"min_freq": min_freq}
    store = _build(tmp_path, "m", {"s1": SAMPLE_1}, storage_mode="bitmask", **kwargs)
    assert _rows(store)[0].tolist() == expected
    recorded = mp._load_matrix_hdf5_metadata(store)[mp.MATRIX_MIN_FREQ_METADATA_KEY]
    assert float(recorded) == (0.05 if min_freq is None else min_freq)


def test_count_store_zeroes_filtered_alleles(tmp_path):
    store = _build(tmp_path, "c", {"s1": SAMPLE_1}, storage_mode="counts", count_dtype="auto")
    counts = _rows(store)[0]
    assert counts[0].tolist() == [95, 5, 0, 0]
    assert counts[1].tolist() == [96, 0, 0, 0]
    assert counts[3].tolist() == [0, 0, 0, 0]
    assert counts[4].tolist() == [99, 0, 0, 0]


@pytest.mark.parametrize("bad", [-0.01, 1.01, float("nan")])
def test_invalid_min_freq_rejected(tmp_path, bad):
    with pytest.raises(ValueError, match="min_freq"):
        _build(tmp_path, "bad", {"s1": SAMPLE_1}, min_freq=bad)


@pytest.mark.parametrize("store_min_freq,expected_pos2", [(0.0, A | T), (0.05, A)])
def test_append_uses_store_min_freq(tmp_path, store_min_freq, expected_pos2):
    store = _build(tmp_path, "a", {"s2": SAMPLE_2}, min_freq=store_min_freq)
    extra = tmp_path / "extra"
    extra.mkdir()
    _write_profile(extra / "s1.parquet", SAMPLE_1)
    mp.append_matrix_hdf5(profile_dir=extra, matrix_hdf5_file=store, memory_limit_gb=1)
    names = [name for _idx, name in mp._load_matrix_hdf5_samples(store)]
    assert names == ["s2", "s1"]
    assert _rows(store)[1][1] == expected_pos2
    assert float(mp._load_matrix_hdf5_metadata(store)[mp.MATRIX_MIN_FREQ_METADATA_KEY]) == store_min_freq


def test_metadata_without_key_reads_as_unfiltered():
    assert mp._matrix_min_freq_from_metadata({}) == 0.0
    assert mp._matrix_min_freq_from_metadata({mp.MATRIX_MIN_FREQ_METADATA_KEY: "0.05"}) == 0.05
    with pytest.raises(ValueError, match="allele-frequency"):
        mp._matrix_min_freq_from_metadata({mp.MATRIX_MIN_FREQ_METADATA_KEY: "2"})


def test_min_freq_changes_popani_and_is_recorded(tmp_path):
    samples = {"s1": SAMPLE_1, "s2": SAMPLE_2}
    results = {}
    for min_freq in (0.0, 0.05):
        store = _build(tmp_path / str(min_freq), "m", samples, min_freq=min_freq)
        out = tmp_path / f"cmp_{min_freq}.duckdb"
        mp.matrix_compare(store, out, calculate="ani", memory_limit_gb=1)
        with duckdb.connect(str(out), read_only=True) as conn:
            results[min_freq] = conn.execute(
                "SELECT total_positions, share_allele_pos FROM matrix_compare_results"
            ).fetchone()
            recorded = dict(conn.execute("SELECT key, value FROM matrix_compare_metadata").fetchall())
        assert float(recorded["matrix_min_freq"]) == min_freq
        parquet = tmp_path / f"cmp_{min_freq}.parquet"
        mp.export_matrix_compare_parquet(out, parquet)
        meta = pl.read_parquet_metadata(parquet)
        assert float(meta[cli.ut.COMPARE_MIN_FREQ_METADATA_KEY]) == min_freq
    # Covered in both at positions 1, 2, 3, 5, 6. Unfiltered, the minor T alleles
    # at 1, 2 and 5 match sample 2; at 5% only position 1 (exactly 5%) survives.
    assert results[0.0] == (5, 5)
    assert results[0.05] == (5, 3)


def test_compare_db_rejects_store_with_different_min_freq(tmp_path):
    samples = {"s1": SAMPLE_1, "s2": SAMPLE_2}
    out = tmp_path / "cmp.duckdb"
    mp.matrix_compare(_build(tmp_path / "a", "m", samples, min_freq=0.05), out, calculate="ani", memory_limit_gb=1)
    rebuilt = _build(tmp_path / "b", "m", samples, min_freq=0.0)
    with pytest.raises(ValueError, match="matrix_min_freq"):
        mp.matrix_compare(rebuilt, out, calculate="ani", memory_limit_gb=1)


def test_compare_db_without_key_resumes_against_unfiltered_store(tmp_path):
    """Compare DBs written before the key existed used unfiltered stores."""
    samples = {"s1": SAMPLE_1, "s2": SAMPLE_2}
    store = _build(tmp_path, "m", samples, min_freq=0.0)
    out = tmp_path / "cmp.duckdb"
    mp.matrix_compare(store, out, calculate="ani", memory_limit_gb=1)
    with duckdb.connect(str(out)) as conn:
        conn.execute("DELETE FROM matrix_compare_metadata WHERE key = 'matrix_min_freq'")
    assert mp.matrix_compare(store, out, calculate="ani", memory_limit_gb=1).requested_pairs == 0
    filtered = _build(tmp_path / "f", "m", samples, min_freq=0.05)
    with pytest.raises(ValueError, match="matrix_min_freq"):
        mp.matrix_compare(filtered, out, calculate="ani", memory_limit_gb=1)


def test_build_matrix_db_cli_min_freq(tmp_path):
    profiles, bed, stb = _inputs(tmp_path, {"s1": SAMPLE_1})
    out = tmp_path / "cli.h5"
    args = ["utilities", "build-matrix-db", "--profile-dir", str(profiles), "--output-file", str(out),
            "--bed-file", str(bed), "--stb-file", str(stb), "--memory-limit-gb", "1"]
    result = CliRunner().invoke(cli.cli, args + ["--min-freq", "0"])
    assert result.exit_code == 0, result.output
    assert "min_freq=0.0" in result.output
    assert _rows(out)[0].tolist() == [A | T, A | T, A | C, 0, A | T, G]
    bad = CliRunner().invoke(cli.cli, args + ["--min-freq", "1.5"])
    assert bad.exit_code != 0


def test_matrix_workflow_rejects_conflicting_min_freq(tmp_path):
    store = _build(tmp_path, "m", {"s1": SAMPLE_1}, min_freq=0.05)
    kwargs = dict(store=store, compare_db=tmp_path / "missing.duckdb", scope="all", min_cov=5, ani_method="popani")
    mw._validate_existing_matrix_run(**kwargs)
    mw._validate_existing_matrix_run(**kwargs, min_freq=0.05)
    with pytest.raises(ValueError, match="min_freq=0.05"):
        mw._validate_existing_matrix_run(**kwargs, min_freq=0.01)
