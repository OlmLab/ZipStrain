"""Metric selection, shared defaults, comparison schema, and route consistency."""

import inspect
import re
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from zipstrain import cli, codons, compare, gene_info as gi, matrix_workflow, options, profile, profile_backends, utils
from zipstrain import visualize as vz

SUMMARY_COLUMNS = {"shared_genes_count", "identical_gene_count", "perc_id_genes"}
FULL_SCHEMA = {
    "sample_1": pl.String, "sample_2": pl.String, "genome": pl.String, "gene_id": pl.UInt64,
    "gene_total_positions": pl.Int64, "gene_share_allele_pos": pl.Int64, "gene_ani": pl.Float64,
    "dN": pl.Float64, "dS": pl.Float64, "pN": pl.Float64, "pS": pl.Float64,
    "total_positions": pl.Int64, "share_allele_pos": pl.Int64, "genome_ani": pl.Float64,
    "max_identical_stretch": pl.Int64,
}


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


def invoke(*args):
    return CliRunner().invoke(cli.cli, [str(arg) for arg in args])


def write_v1_copy(path, name):
    """Write a 1.x shaped profile, which labels each row with its gene."""
    target = path.with_name(f"{name}_profile.parquet")
    pl.read_parquet(path).with_columns(pl.lit("chr1_1").alias("gene")).write_parquet(target)
    return target


@pytest.fixture
def pair(assets):
    root, _, _, stb, catalog = assets
    return (root, stb, catalog, write_sample(root, catalog, "a"),
            write_sample(root, catalog, "b", "AAAGGTTTTCCCAAAGGTTAA"))


def write_pair_codons(root, catalog, a, b):
    gene = gi.read_gene_info(catalog).filter(pl.col("gene") == "chr1_1")["gene_id"][0]
    write_counts(root, catalog, "a", [(gene, 0, "AAA", 10), (gene, 1, "GGA", 10)], source_profile=a)
    write_counts(root, catalog, "b", [(gene, 0, "AAG", 10), (gene, 1, "GGT", 10)], source_profile=b)
    return root / "a_codon_counts.parquet", root / "b_codon_counts.parquet"


# Metric selection.
def test_only_public_metric_names_are_accepted():
    assert options.parse_calculate("all") == options.METRICS
    assert options.parse_calculate(None) == options.METRICS
    assert options.parse_calculate("gene_ani+genome_ani") == ("genome_ani", "gene_ani")
    assert options.parse_calculate("gene_selection,max_identical_stretch") == ("max_identical_stretch", "gene_selection")


@pytest.mark.parametrize(("removed", "current"), [
    ("dnds", "gene_selection"), ("ani", "genome_ani"), ("popani", "genome_ani"), ("ibs", "max_identical_stretch"),
    ("max_block", "max_identical_stretch"), ("max_consecutive_length", "max_identical_stretch"),
    ("identical_genes", "gene_ani"), ("gene_identity", "gene_ani"), ("id_genes", "gene_ani"),
])
def test_removed_metric_names_point_to_the_current_name(removed, current):
    with pytest.raises(ValueError, match=rf"Unknown metric '{removed}'\. Did you mean '{current}'\?"):
        options.parse_calculate(removed)


def test_unknown_metric_lists_choices():
    with pytest.raises(ValueError, match="Choose from: all, genome_ani"):
        options.parse_calculate("bogus")


def test_all_keeps_supported_metrics_and_warns_about_skipped_ones():
    assert options.resolve_metrics("all", genes_available=True, codon_tables_present=["a", "b"]) == options.METRICS
    with pytest.warns(UserWarning, match="Skipping gene_selection"):
        assert options.resolve_metrics(
            "all", genes_available=True, codon_tables_missing=["a", "b"],
        ) == ("genome_ani", "max_identical_stretch", "gene_ani")
    with pytest.warns(UserWarning, match="Skipping gene_ani and gene_selection"):
        assert options.resolve_metrics("all", genes_available=False) == options.GENOME_METRICS


def test_partial_codon_tables_are_an_error_even_for_all():
    with pytest.raises(ValueError, match="Missing for: b"):
        options.resolve_metrics("all", genes_available=True, codon_tables_present=["a"], codon_tables_missing=["b"])


@pytest.mark.parametrize("metric", options.GENE_METRICS)
def test_explicit_gene_metric_without_catalog_is_an_error(metric):
    with pytest.raises(ValueError, match="--preprofile-dir"):
        options.resolve_metrics(metric, genes_available=False)


def test_explicit_gene_selection_names_samples_without_codon_tables():
    with pytest.raises(ValueError, match="Missing for: a, b"):
        options.resolve_metrics("gene_selection", genes_available=True, codon_tables_missing=["b", "a"])


# Shared defaults.
@pytest.mark.parametrize("function", [
    compare.compare_profiles, compare.compare_profile_collection, compare.compare_genomes,
    compare.compare_genes, compare.compare_genomes_polars, compare.compare_genes_polars,
    matrix_workflow.run_matrix_compare, utils.chunk_genome_compare,
])
def test_gene_threshold_default_is_shared(function):
    assert inspect.signature(function).parameters["min_gene_compare_len"].default == options.MIN_GENE_COMPARE_LEN == 1


@pytest.mark.parametrize("command", ["compare", "utilities single_compare_genome", "utilities single_compare_gene",
                                     "utilities compare-pair", "utilities chunk-genome-compare"])
def test_gene_threshold_cli_default_and_help_are_shared(command):
    group = cli.cli
    for name in command.split():
        group = group.commands[name]
    option = next(param for param in group.params if param.name == "min_gene_compare_len")
    assert option.default == options.MIN_GENE_COMPARE_LEN
    assert "shared positions" in option.help


def test_nextflow_defaults_match_shared_defaults():
    workflow = (Path(__file__).resolve().parents[2] / "zipstrain.nf").read_text()
    value = lambda name: re.search(rf"^params\.{name}\s*=\s*\"?([^\"\n]+)\"?$", workflow, re.M).group(1)
    assert int(value("min_gene_compare_len")) == options.MIN_GENE_COMPARE_LEN
    assert int(value("min_cov")) == options.COMPARE_MIN_COV
    assert value("compare_calculate") == options.CALCULATE
    assert value("compare_ani_method") == options.ANI_METHOD
    assert value("allele_integration") == options.ALLELE_INTEGRATION
    assert float(value("min_read_ani")) == options.PROFILE_MIN_READ_ANI
    assert value("read_inclusion") == options.PROFILE_READ_INCLUSION


def test_python_api_read_filters_match_cli_defaults():
    for function in (profile.finish_profile_outputs, codons.write_codon_counts):
        parameters = inspect.signature(function).parameters
        assert parameters["min_read_ani"].default == options.PROFILE_MIN_READ_ANI
        assert parameters["read_inclusion"].default == options.PROFILE_READ_INCLUSION
    fields = profile_backends.ProfileRequest.__dataclass_fields__
    assert fields["min_read_ani"].default == options.PROFILE_MIN_READ_ANI
    assert fields["read_inclusion"].default == options.PROFILE_READ_INCLUSION


# Routes that used to write zeros for V2 profiles.
def test_single_compare_genome_all_writes_genome_metrics_and_warns(pair, tmp_path):
    _, _, _, a, b = pair
    output = tmp_path / "out.parquet"
    with pytest.warns(UserWarning, match="genome metrics only"):
        result = invoke("utilities", "single_compare_genome", "--profile-1", a, "--profile-2", b, "-o", output)
    assert result.exit_code == 0, result.output
    frame = pl.read_parquet(output)
    assert not SUMMARY_COLUMNS & set(frame.columns)
    row = frame.row(0, named=True)
    assert row["total_positions"] == 21 and row["share_allele_pos"] == 20 and row["max_identical_stretch"] == 15


@pytest.mark.parametrize("extra", [(), ("-s", "STB")])
def test_single_compare_genome_explicit_gene_ani_without_assets_fails(pair, tmp_path, extra):
    _, stb, _, a, b = pair
    output = tmp_path / "out.parquet"
    extra = [stb if value == "STB" else value for value in extra]
    result = invoke("utilities", "single_compare_genome", "--profile-1", a, "--profile-2", b,
                    "--calculate", "gene_ani", *extra, "-o", output)
    assert result.exit_code != 0 and "--preprofile-dir" in result.output
    assert not output.exists()


def test_single_compare_genome_scope_and_stb_on_v2(pair, tmp_path):
    _, stb, _, a, b = pair
    output = tmp_path / "out.parquet"
    result = invoke("utilities", "single_compare_genome", "--profile-1", a, "--profile-2", b, "-s", stb,
                    "--calculate", "genome_ani", "-g", "genomeB", "-o", output)
    assert result.exit_code == 0, result.output
    frame = pl.read_parquet(output)
    assert frame["genome"].to_list() == ["genomeB"]
    assert frame["total_positions"].to_list() == [0] and frame["genome_ani"].to_list() == [None]


def test_single_compare_genome_v2_matches_compare_pair(pair, tmp_path):
    _, _, _, a, b = pair
    selection = "genome_ani+max_identical_stretch"
    single, double = tmp_path / "single.parquet", tmp_path / "pair"
    assert invoke("utilities", "single_compare_genome", "--profile-1", a, "--profile-2", b,
                  "--calculate", selection, "-o", single).exit_code == 0
    assert invoke("utilities", "compare-pair", "--profile-1", a, "--profile-2", b,
                  "--calculate", selection, "--output-dir", double).exit_code == 0
    columns = ["genome", "total_positions", "share_allele_pos", "genome_ani", "max_identical_stretch"]
    assert pl.read_parquet(single).select(columns).equals(pl.read_parquet(double / "comparisons.parquet").select(columns))


def chunk_table(tmp_path, a, b):
    table = tmp_path / "pairs.parquet"
    pl.DataFrame({"sample_name_1": ["a"], "sample_name_2": ["b"],
                  "profile_location_1": [str(a)], "profile_location_2": [str(b)]}).write_parquet(table)
    return table


def test_chunk_genome_compare_all_on_v2_writes_genome_metrics(pair, tmp_path):
    _, _, _, a, b = pair
    output = tmp_path / "chunk.parquet"
    with pytest.warns(UserWarning, match="genome metrics only"):
        result = invoke("utilities", "chunk-genome-compare", "--pair-table", chunk_table(tmp_path, a, b),
                        "--output-file", output)
    assert result.exit_code == 0, result.output
    frame = pl.read_parquet(output)
    assert not SUMMARY_COLUMNS & set(frame.columns)
    assert frame["share_allele_pos"].to_list() == [20] and frame["max_identical_stretch"].to_list() == [15]


def test_chunk_genome_compare_rejects_gene_metrics(pair, tmp_path):
    _, _, _, a, b = pair
    result = invoke("utilities", "chunk-genome-compare", "--pair-table", chunk_table(tmp_path, a, b),
                    "--calculate", "gene_ani", "--output-file", tmp_path / "chunk.parquet")
    assert result.exit_code != 0 and "per gene" in result.output


def test_single_compare_gene_v2_without_assets_fails(pair, tmp_path):
    _, _, _, a, b = pair
    result = invoke("utilities", "single_compare_gene", "--profile-1", a, "--profile-2", b, "-o", tmp_path / "out.parquet")
    assert result.exit_code != 0 and "--preprofile-dir" in result.output


def test_no_route_adds_a_fake_gene_column(pair):
    _, _, _, a, b = pair
    first, second = compare.polars_prefilter_by_scope(a, b, genome_scope="genomeA")
    assert "gene" not in first.collect_schema().names() and "gene" not in second.collect_schema().names()


@pytest.mark.parametrize("v1_first", [True, False])
@pytest.mark.parametrize("route", ["single_genome", "single_gene", "pair", "chunk", "collection", "matrix"])
def test_mixed_generations_are_rejected_in_both_orders(pair, tmp_path, v1_first, route):
    _, stb, _, a, b = pair
    v1 = write_v1_copy(b, "v1")
    first, second = (v1, a) if v1_first else (a, v1)
    if route == "matrix":
        with pytest.raises(ValueError, match="same version"):
            matrix_workflow.run_matrix_compare(profiles=[("x", first), ("y", second)], run_dir=tmp_path / "run",
                                               stb_file=stb, bed_file=stb)
        return
    if route == "collection":
        with pytest.raises(ValueError, match="same version"):
            compare.compare_profile_collection([("x", first), ("y", second)], tmp_path / "run")
        return
    if route == "chunk":
        args = ["chunk-genome-compare", "--pair-table", chunk_table(tmp_path, first, second), "--calculate", "genome_ani",
                "--output-file", tmp_path / "out.parquet"]
    elif route == "pair":
        args = ["compare-pair", "--profile-1", first, "--profile-2", second, "--calculate", "genome_ani", "--output-dir", tmp_path / "out"]
    else:
        command = "single_compare_genome" if route == "single_genome" else "single_compare_gene"
        args = [command, "--profile-1", first, "--profile-2", second, "-o", tmp_path / "out.parquet"]
    result = invoke("utilities", *args)
    assert result.exit_code != 0
    assert "same version" in result.output


def test_v1_all_means_genome_metrics_and_gene_ani_needs_the_gene_route(pair):
    _, _, _, a, b = pair
    v1_a, v1_b = write_v1_copy(a, "v1a"), write_v1_copy(b, "v1b")
    assert compare.compare_genomes(v1_a, v1_b, calculate="all").collect_schema().names() == [
        "genome", "total_positions", "share_allele_pos", "genome_ani", "max_identical_stretch"]
    with pytest.raises(ValueError, match="single_compare_gene"):
        compare.compare_genomes(v1_a, v1_b, calculate="gene_ani")


def test_compare_v1_gene_ani_cannot_combine_with_genome_metrics(pair, tmp_path):
    _, _, _, a, b = pair
    table = tmp_path / "profiles.csv"
    table.write_text(f"profile_name,profile_location\nx,{write_v1_copy(a, 'v1a')}\ny,{write_v1_copy(b, 'v1b')}\n")
    for selection, message in (("gene_ani+genome_ani", "on its own"), ("gene_selection", "requires V2"), ("dnds", "Did you mean")):
        result = invoke("compare", "--profile-db", table, "--allow-mismatch", "--run-dir", tmp_path / "run", "--calculate", selection)
        assert result.exit_code != 0 and message in result.output


@pytest.mark.parametrize("flag", ["--dnds", "--compare-genes"])
def test_old_flags_are_gone(flag):
    assert "No such option" in CliRunner().invoke(cli.cli, ["compare", flag]).output


# Comparison table.
def test_full_table_schema_and_order(pair):
    root, _, catalog, a, b = pair
    counts_a, counts_b = write_pair_codons(root, catalog, a, b)
    frame = compare.compare_profiles(a, b, preprofile_dir=root, codon_counts_1=counts_a, codon_counts_2=counts_b).collect()
    assert dict(frame.schema) == FULL_SCHEMA
    assert list(frame.columns) == list(FULL_SCHEMA)


@pytest.mark.parametrize(("selection", "expected"), [
    ("gene_ani", list(FULL_SCHEMA)[:7]),
    ("gene_selection", [*list(FULL_SCHEMA)[:6], "dN", "dS", "pN", "pS"]),
    ("genome_ani", ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_ani"]),
    ("max_identical_stretch", ["sample_1", "sample_2", "genome", "max_identical_stretch"]),
])
def test_each_metric_adds_its_columns(pair, selection, expected):
    root, _, catalog, a, b = pair
    counts_a, counts_b = write_pair_codons(root, catalog, a, b)
    frame = compare.compare_profiles(a, b, preprofile_dir=root, calculate=selection,
                                     codon_counts_1=counts_a, codon_counts_2=counts_b)
    assert frame.collect_schema().names() == expected


def test_genome_ani_is_null_without_shared_positions(pair):
    root, _, _, a, b = pair
    frame = compare.compare_profiles(a, b, preprofile_dir=root, calculate="genome_ani").collect()
    row = frame.filter(pl.col("genome") == "genomeB").row(0, named=True)
    assert row["total_positions"] == 0 and row["genome_ani"] is None


def test_selection_ratios_are_derived_and_null_on_zero():
    frame = gi.with_selection_ratios(pl.DataFrame({"dN": [1.0, 1.0, None], "dS": [2.0, 0.0, 1.0],
                                                   "pN": [3.0, 0.0, 1.0], "pS": [2.0, 0.0, None]}))
    assert frame["dN_dS"].to_list() == [0.5, None, None]
    assert frame["pN_pS"].to_list() == [1.5, None, None]


def collection_table(root, a, b, name):
    table = root / f"{name}.csv"
    table.write_text(f"profile_name,profile_location\na,{a}\nb,{b}\n")
    return table


def test_compare_all_uses_codon_tables_when_every_sample_has_one(pair, tmp_path):
    root, _, catalog, a, b = pair
    write_pair_codons(root, catalog, a, b)
    result = invoke("compare", "--profile-db", collection_table(root, a, b, "full"), "--preprofile-dir", root,
                    "--run-dir", tmp_path / "run", "--no-csv")
    assert result.exit_code == 0, result.output
    output = tmp_path / "run" / "all_comparisons.parquet"
    assert list(pl.read_parquet(output).columns) == list(FULL_SCHEMA)
    assert pl.read_parquet_metadata(output)["zipstrain_compare_calculate"] == "+".join(options.METRICS)


def test_compare_all_without_codon_tables_warns_and_records_what_ran(pair, tmp_path):
    root, _, _, a, b = pair
    with pytest.warns(UserWarning, match="Skipping gene_selection"):
        result = invoke("compare", "--profile-db", collection_table(root, a, b, "none"), "--preprofile-dir", root,
                        "--run-dir", tmp_path / "run", "--no-csv")
    assert result.exit_code == 0, result.output
    output = tmp_path / "run" / "all_comparisons.parquet"
    assert pl.read_parquet_metadata(output)["zipstrain_compare_calculate"] == "genome_ani+max_identical_stretch+gene_ani"
    assert "dN" not in pl.read_parquet(output).columns


def test_compare_with_codon_tables_for_some_samples_fails(pair, tmp_path):
    root, _, catalog, a, b = pair
    counts_a, counts_b = write_pair_codons(root, catalog, a, b)
    counts_b.unlink()
    result = invoke("compare", "--profile-db", collection_table(root, a, b, "some"), "--preprofile-dir", root,
                    "--run-dir", tmp_path / "run", "--no-csv")
    assert result.exit_code != 0 and "Missing for: b" in result.output


def test_threshold_is_not_part_of_the_cache_key_for_genome_rows(pair, tmp_path):
    root, _, _, a, b = pair
    cache = tmp_path / "run" / "intermediate_files" / "v2_pairs"
    run = lambda selection, threshold: compare.compare_profile_collection(
        [("a", a), ("b", b)], tmp_path / "run", preprofile_dir=root, calculate=selection, min_gene_compare_len=threshold)
    first = run("genome_ani", 1)
    assert pl.read_parquet_metadata(first)["zipstrain_compare_min_gene_compare_len"] == "NA"
    files = set(cache.glob("*.parquet"))
    run("genome_ani", 50)
    assert set(cache.glob("*.parquet")) == files
    run("gene_ani", 1)
    gene_files = set(cache.glob("*.parquet"))
    run("gene_ani", 2)
    assert len(set(cache.glob("*.parquet"))) == len(gene_files) + 1


# Plotting functions read gene rows.
def gene_rows(sample_2="s2", genome="g", identical=2, others=2, positions=200, genome_ani=99.0, stretch=5000, total=20000):
    size = identical + others
    return pl.DataFrame({
        "sample_1": ["s1"] * size, "sample_2": [sample_2] * size, "genome": [genome] * size,
        "gene_id": list(range(size)), "gene_total_positions": [positions] * size,
        "gene_share_allele_pos": [positions] * identical + [positions - 1] * others,
        "total_positions": [total] * size, "genome_ani": [genome_ani] * size,
        "max_identical_stretch": [stretch] * size,
    })


POPULATIONS = pl.DataFrame({"sample": ["s1", "s2", "s3"], "population": ["p", "p", "p"]}).lazy()


def test_summary_counts_each_pair_once_and_applies_the_threshold():
    rows = pl.concat([gene_rows(identical=1, others=3), gene_rows(identical=1, others=0, positions=99).with_columns(pl.lit(50, dtype=pl.Int64).alias("gene_id"))])
    summary = vz._summarize_gene_rows(rows.lazy(), 100, ("genome_ani",)).collect()
    assert summary.height == 1
    assert summary.row(0, named=True) == {"sample_1": "s1", "sample_2": "s2", "genome": "g", "shared_genes_count": 4,
                                          "identical_gene_count": 1, "genome_ani": 99.0, "perc_id_genes": 25.0}


def test_calculate_ibs_keeps_one_value_per_pair():
    rows = pl.concat([gene_rows(identical=0, others=500, stretch=11), gene_rows("s3", identical=500, others=0, stretch=99)])
    result = vz.calculate_ibs(POPULATIONS, rows.lazy())
    assert result["within_population_p|p"].to_list() == [[11]]


def test_identical_fraction_is_one_point_per_pair():
    rows = pl.concat([gene_rows(identical=75, others=75, genome_ani=98.0), gene_rows("s3", identical=30, others=120, genome_ani=97.0)])
    result = vz.calculate_identical_frac_vs_popani("g", "p", "p", POPULATIONS, rows.lazy()).row(0, named=True)
    assert sorted(zip(result["genome_ani"], result["perc_id_genes"])) == [(97.0, 20.0), (98.0, 50.0)]


@pytest.mark.parametrize(("drop", "metric"), [
    ("gene_total_positions", "gene_ani"), ("max_identical_stretch", "max_identical_stretch"), ("total_positions", "genome_ani"),
])
def test_calculate_ibs_names_the_missing_metric(drop, metric):
    with pytest.raises(ValueError, match=rf"no {metric} results.*--calculate {metric}"):
        vz.calculate_ibs(POPULATIONS, gene_rows().drop(drop).lazy())


def test_identical_fraction_names_the_missing_metric():
    with pytest.raises(ValueError, match="no genome_ani results"):
        vz.calculate_identical_frac_vs_popani("g", "p", "p", POPULATIONS, gene_rows().drop("genome_ani").lazy())


def test_threshold_below_the_table_threshold_is_an_error(tmp_path):
    path = tmp_path / "table.parquet"
    gene_rows().write_parquet(path, metadata={"zipstrain_compare_min_gene_compare_len": "200"})
    with pytest.raises(ValueError, match="below the threshold"):
        vz._summarize_gene_rows(path, 100).collect()
    assert vz._summarize_gene_rows(path, 200).collect().height == 1


# Tables written before the rename.
def test_1x_comparison_table_upgrades_when_extended(tmp_path):
    from zipstrain import database as db
    location = tmp_path / "old.parquet"
    pl.DataFrame({"genome": ["g"], "total_positions": [10], "share_allele_pos": [10], "genome_ani": [100.0],
                  "max_consecutive_length": [7], "shared_genes_count": [1], "identical_gene_count": [1],
                  "perc_id_genes": [100.0], "sample_1": ["a"], "sample_2": ["b"]}).write_parquet(location)
    frame = db.GenomeComparisonDatabase._normalize_comp_db_columns(
        db.GenomeComparisonDatabase, pl.scan_parquet(location)).collect()
    assert frame.columns == db.GenomeComparisonDatabase.COLUMN_NAMES
    assert frame["max_identical_stretch"].to_list() == [7]


def test_matrix_store_with_the_previous_column_name_resumes(pair, tmp_path):
    pytest.importorskip("h5py")
    duckdb = pytest.importorskip("duckdb")
    root, stb, _, a, b = pair
    bed = root / "regions.bed"
    bed.write_text("chr1\t0\t21\nchr2\t0\t9\n")
    run = tmp_path / "matrix"
    args = ["compare", "--method", "matrix", "--profile-db", collection_table(root, a, b, "matrix"), "--stb-file", stb,
            "--bed-file", bed, "--preprofile-dir", root, "--run-dir", run, "--calculate", "genome_ani+max_identical_stretch"]
    assert invoke(*args).exit_code == 0
    output = run / "all_comparisons.parquet"
    expected = pl.read_parquet(output).sort("genome")
    database = next((run / "intermediate_files").glob("*.duckdb"))
    connection = duckdb.connect(str(database))
    connection.execute("ALTER TABLE matrix_compare_results RENAME COLUMN max_identical_stretch TO max_consecutive_length")
    connection.close()
    output.unlink()
    result = invoke(*args)
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert pl.read_parquet(output).sort("genome").equals(expected)
