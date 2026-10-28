import pathlib

import polars as pl
import pytest

from zipstrain import matrix_workflow


def _write_profile(path, genomes):
    pl.DataFrame(
        {
            "chrom": ["chr1"] * len(genomes),
            "genome": genomes,
            "gene": ["NA"] * len(genomes),
            "pos": list(range(len(genomes))),
            "A": [1] * len(genomes),
            "T": [0] * len(genomes),
            "C": [0] * len(genomes),
            "G": [0] * len(genomes),
        }
    ).write_parquet(path)


def test_discover_bed_file_walks_up_to_profiling_assets(tmp_path):
    run = tmp_path / "out_profile"
    (run / "profiling_assets").mkdir(parents=True)
    bed = run / "profiling_assets" / "genomes_bed_file.bed"
    bed.write_text("chr1\t0\t10\n")
    (run / "sampleA").mkdir()
    profile = run / "sampleA" / "sampleA_profile.parquet"
    profile.write_text("x")

    assert matrix_workflow.discover_bed_file([profile]) == bed
    assert matrix_workflow.discover_bed_file([tmp_path / "nowhere.parquet"]) is None


def test_clean_stb_strips_whitespace(tmp_path):
    stb = tmp_path / "in.stb"
    stb.write_text("scaffold_0 \t genomeA\nscaffold_1\tgenomeA \n\n")
    cleaned = matrix_workflow._clean_stb(stb, tmp_path / "clean.stb")
    assert cleaned.read_text().splitlines() == ["scaffold_0\tgenomeA", "scaffold_1\tgenomeA"]


def test_stage_profiles_drops_unbinned_and_symlinks_clean(tmp_path):
    with_na = tmp_path / "with_na.parquet"
    clean = tmp_path / "clean.parquet"
    _write_profile(with_na, ["NA", "genomeA", "NA"])
    _write_profile(clean, ["genomeA", "genomeA"])

    dest = tmp_path / "stage"
    matrix_workflow._stage_profiles([("s1", with_na), ("s2", clean)], dest)

    # s1 had NA rows -> filtered copy (real file, no NA rows).
    s1 = pl.read_parquet(dest / "s1.parquet")
    assert "NA" not in s1["genome"].to_list()
    assert s1.height == 1
    # s2 had no NA -> symlinked to the original.
    assert (dest / "s2.parquet").is_symlink()


def test_run_matrix_compare_requires_stb(tmp_path):
    with pytest.raises(ValueError, match="--stb-file is required"):
        matrix_workflow.run_matrix_compare(
            profiles=[("s1", tmp_path / "p1.parquet")],
            run_dir=tmp_path / "run",
            stb_file=None,
        )


def test_run_matrix_compare_errors_when_bed_not_found(tmp_path):
    stb = tmp_path / "ref.stb"
    stb.write_text("chr1\tgenomeA\n")
    profile = tmp_path / "p1.parquet"
    _write_profile(profile, ["genomeA"])
    with pytest.raises(ValueError, match="Could not auto-discover a bed file"):
        matrix_workflow.run_matrix_compare(
            profiles=[("s1", profile)],
            run_dir=tmp_path / "run",
            stb_file=stb,
        )


def test_run_matrix_compare_builds_count_store_for_conani(tmp_path, monkeypatch):
    stb = tmp_path / "ref.stb"
    bed = tmp_path / "ref.bed"
    profile = tmp_path / "p1.parquet"
    stb.write_text("chr1\tgenomeA\n")
    bed.write_text("chr1\t0\t1\n")
    _write_profile(profile, ["genomeA"])
    captured = {}

    def fake_build(**kwargs):
        captured["build"] = kwargs
        pathlib.Path(kwargs["output_file"]).touch()

    def fake_compare(**kwargs):
        captured["compare"] = kwargs

    def fake_export(**kwargs):
        output = pathlib.Path(kwargs["output_file"])
        pl.DataFrame({"genome": ["genomeA"], "genome_ani": [100.0]}).write_parquet(output)
        return output

    monkeypatch.setattr(matrix_workflow.mp, "build_matrix_hdf5", fake_build)
    monkeypatch.setattr(matrix_workflow.mp, "matrix_compare", fake_compare)
    monkeypatch.setattr(matrix_workflow.mp, "export_matrix_compare_parquet", fake_export)

    matrix_workflow.run_matrix_compare(
        profiles=[("s1", profile)],
        run_dir=tmp_path / "run",
        stb_file=stb,
        bed_file=bed,
        ani_method="conani",
        min_cov=7,
    )

    assert captured["build"]["storage_mode"] == matrix_workflow.mp.MATRIX_STORAGE_COUNTS
    assert captured["build"]["count_dtype"] is None
    assert captured["build"]["min_cov"] == 7
    assert captured["compare"]["ani_method"] == "conani"
    assert captured["compare"]["min_cov"] == 7


def test_run_matrix_compare_validates_existing_contract_before_append(tmp_path):
    stb = tmp_path / "ref.stb"
    bed = tmp_path / "ref.bed"
    run_dir = tmp_path / "run"
    stb.write_text("chr1\tgenomeA\n")
    bed.write_text("chr1\t0\t1\n")
    profiles = []
    for sample_name in ("s1", "s2", "s3"):
        profile = tmp_path / f"{sample_name}.parquet"
        _write_profile(profile, ["genomeA"])
        pl.read_parquet(profile).with_columns(pl.lit(5).alias("A")).write_parquet(profile)
        profiles.append((sample_name, profile))

    matrix_workflow.run_matrix_compare(
        profiles=profiles[:2],
        run_dir=run_dir,
        stb_file=stb,
        bed_file=bed,
        ani_method="conani",
        min_cov=5,
    )

    with pytest.raises(ValueError, match="ani_method expected='popani' actual='conani'"):
        matrix_workflow.run_matrix_compare(
            profiles=profiles,
            run_dir=run_dir,
            stb_file=stb,
            bed_file=bed,
            ani_method="popani",
            min_cov=5,
        )

    store = run_dir / matrix_workflow.INTERMEDIATE_DIRNAME / matrix_workflow.MATRIX_STORE_NAME
    assert matrix_workflow._existing_store_sample_names(store) == {"s1", "s2"}


@pytest.mark.parametrize(
    ("kernel", "calculate", "ani_method", "expect_hint"),
    [
        ("anchor", "ani", "popani", True),
        ("anchor", "ani", "conani", True),
        ("anchor", "all", "popani", False),
        ("anchor", "ani+ibs", "popani", False),
        ("anchor", "ani", "cosani_0.4", False),
        ("block", "ani", "popani", False),
    ],
)
def test_block_kernel_hint(kernel, calculate, ani_method, expect_hint):
    hint = matrix_workflow.block_kernel_hint(kernel, calculate=calculate, ani_method=ani_method)
    assert (hint is not None) == expect_hint
    if expect_hint:
        assert "--compare-kernel block" in hint


def test_resolve_calculate_adds_gene_for_gene_runs():
    assert matrix_workflow.resolve_calculate("ani", compare_genes=True) == "ani+gene"
    assert matrix_workflow.resolve_calculate("ani", compare_genes=False) == "ani"
    with pytest.raises(ValueError, match="Unsupported matrix calculation"):
        matrix_workflow.resolve_calculate("identical_genes", compare_genes=False)


def test_run_matrix_compare_forwards_calculate_kernel_and_progress(tmp_path, monkeypatch):
    stb = tmp_path / "ref.stb"
    bed = tmp_path / "ref.bed"
    profile = tmp_path / "p1.parquet"
    stb.write_text("chr1\tgenomeA\n")
    bed.write_text("chr1\t0\t1\n")
    _write_profile(profile, ["genomeA"])
    captured = {}

    monkeypatch.setattr(
        matrix_workflow.mp,
        "build_matrix_hdf5",
        lambda **kwargs: pathlib.Path(kwargs["output_file"]).touch(),
    )
    monkeypatch.setattr(matrix_workflow.mp, "matrix_compare", lambda **kwargs: captured.update(kwargs))

    def fake_export(**kwargs):
        output = pathlib.Path(kwargs["output_file"])
        pl.DataFrame({"genome": ["genomeA"], "genome_ani": [100.0]}).write_parquet(output)
        return output

    monkeypatch.setattr(matrix_workflow.mp, "export_matrix_compare_parquet", fake_export)

    def on_progress(event):
        return None

    matrix_workflow.run_matrix_compare(
        profiles=[("s1", profile)],
        run_dir=tmp_path / "run",
        stb_file=stb,
        bed_file=bed,
        calculate="ani",
        compare_kernel="block",
        compare_progress_callback=on_progress,
    )

    assert captured["calculate"] == "ani"
    assert captured["compare_kernel"] == "block"
    assert captured["progress_callback"] is on_progress


@pytest.mark.parametrize("ani_method", ["popani", "conani"])
def test_ani_only_workflow_resumes_with_block_kernel_and_new_samples(tmp_path, ani_method):
    stb = tmp_path / "ref.stb"
    bed = tmp_path / "ref.bed"
    stb.write_text("chr1\tgenomeA\n")
    bed.write_text("chr1\t0\t4\n")
    profiles = []
    for name, a, t in [
        ("s1", [10, 10, 0, 10], [0, 0, 10, 0]),
        ("s2", [10, 0, 0, 10], [0, 10, 10, 0]),
        ("s3", [10, 10, 0, 0], [0, 0, 10, 10]),
    ]:
        path = tmp_path / f"{name}.parquet"
        pl.DataFrame({
            "chrom": ["chr1"] * 4, "genome": ["genomeA"] * 4,
            "gene": ["NA"] * 4, "pos": [1, 2, 3, 4],
            "A": a, "T": t, "C": [0] * 4, "G": [0] * 4,
        }).write_parquet(path)
        profiles.append((name, path))
    options = dict(stb_file=stb, bed_file=bed, calculate="ani", ani_method=ani_method,
                   memory_limit_gb=1)
    run = tmp_path / "incremental"
    initial = matrix_workflow.run_matrix_compare(profiles=profiles[:2], run_dir=run, **options)
    before = pl.read_parquet(initial)
    events = []
    resumed = matrix_workflow.run_matrix_compare(
        profiles=profiles, run_dir=run, compare_kernel="block",
        compare_progress_callback=events.append, **options,
    )
    actual = pl.read_parquet(resumed).sort(["sample_1", "sample_2", "genome"])
    fresh = matrix_workflow.run_matrix_compare(profiles=profiles, run_dir=tmp_path / "fresh", **options)
    expected = pl.read_parquet(fresh).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)
    assert actual.height == 3
    assert before.height == 1
    assert actual.filter((pl.col("sample_1") == "s1") & (pl.col("sample_2") == "s2")).equals(before)
    assert events and events[-1]["completed"] == events[-1]["total"]

    # Reject a changed calculation contract before adding a fourth sample.
    with pytest.raises(ValueError, match="calculate expected='ani\\+ibs' actual='ani'"):
        matrix_workflow.run_matrix_compare(
            profiles=[*profiles, ("s4", profiles[0][1])], run_dir=run,
            **{**options, "calculate": "all"},
        )
    store = run / matrix_workflow.INTERMEDIATE_DIRNAME / matrix_workflow.MATRIX_STORE_NAME
    assert matrix_workflow._existing_store_sample_names(store) == {"s1", "s2", "s3"}
