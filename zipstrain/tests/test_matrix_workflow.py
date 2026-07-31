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


def test_run_matrix_compare_combines_multiple_ani_methods(tmp_path):
    stb = tmp_path / "ref.stb"
    bed = tmp_path / "ref.bed"
    run_dir = tmp_path / "run"
    stb.write_text("chr1\tgenomeA\n")
    bed.write_text("chr1\t0\t4\n")
    common = {
        "chrom": ["chr1"] * 4,
        "genome": ["genomeA"] * 4,
        "gene": ["NA"] * 4,
        "pos": [0, 1, 2, 3],
    }
    p1 = tmp_path / "s1.parquet"
    p2 = tmp_path / "s2.parquet"
    pl.DataFrame(
        {
            **common,
            "A": [10, 6, 5, 0],
            "T": [0, 0, 0, 10],
            "C": [0, 4, 5, 0],
            "G": [0, 0, 0, 0],
        }
    ).write_parquet(p1)
    pl.DataFrame(
        {
            **common,
            "A": [10, 4, 0, 0],
            "T": [0, 0, 0, 10],
            "C": [0, 6, 0, 0],
            "G": [0, 0, 10, 0],
        }
    ).write_parquet(p2)

    output = matrix_workflow.run_matrix_compare(
        profiles=[("s1", p1), ("s2", p2)],
        run_dir=run_dir,
        stb_file=stb,
        bed_file=bed,
        ani_method="popani,conani,cosani_0.95",
        calculate="ani+ibs",
        min_cov=1,
        backend="numpy",
    )

    out = pl.read_parquet(output)
    assert out.columns == [
        "sample_1",
        "sample_2",
        "genome",
        "total_positions",
        "share_allele_pos_popani",
        "genome_ani_popani",
        "share_allele_pos_conani",
        "genome_ani_conani",
        "share_allele_pos_cosani_0_95",
        "genome_ani_cosani_0_95",
        "max_consecutive_length_popani",
        "max_consecutive_length_conani",
        "max_consecutive_length_cosani_0_95",
    ]
    assert out["share_allele_pos_popani"].to_list() == [3]
    assert out["share_allele_pos_conani"].to_list() == [2]
    assert out["share_allele_pos_cosani_0_95"].to_list() == [2]
    metadata = pl.read_parquet_metadata(output)
    assert metadata["zipstrain_compare_ani_method"] == "popani,conani,cosani_0.95"
    intermediate = run_dir / matrix_workflow.INTERMEDIATE_DIRNAME
    assert (intermediate / "matrix_compare_popani.duckdb").exists()
    assert (intermediate / "matrix_compare_conani.duckdb").exists()
    assert (intermediate / "matrix_compare_cosani_0_95.duckdb").exists()


def test_combine_matrix_method_exports_suffixes_gene_ani(tmp_path):
    exports = []
    for method, value in (("popani", 100.0), ("conani", 75.0)):
        path = tmp_path / f"{method}.parquet"
        pl.DataFrame(
            {
                "sample_1": ["s1"],
                "sample_2": ["s2"],
                "genome": ["g1"],
                "gene": ["gene1"],
                "gene_pop_ani": [value],
            }
        ).write_parquet(
            path,
            metadata={"zipstrain_compare_ani_method": method},
        )
        exports.append((method, path))

    output = matrix_workflow._combine_matrix_method_exports(
        method_exports=exports,
        output=tmp_path / "combined.parquet",
        table="gene",
        calculate="ani+gene",
    )

    out = pl.read_parquet(output)
    assert out.columns == [
        "sample_1",
        "sample_2",
        "genome",
        "gene",
        "gene_ani_popani",
        "gene_ani_conani",
    ]
    assert out["gene_ani_popani"].to_list() == [100.0]
    assert out["gene_ani_conani"].to_list() == [75.0]


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
