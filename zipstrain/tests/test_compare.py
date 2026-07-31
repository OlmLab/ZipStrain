import random

import pytest
import polars as pl
from zipstrain import compare


a_chr1=[ 1, 0, 1, 0, 0, 0, 4, 4, 2, 1]
t_chr1=[ 0, 2, 3, 1, 0, 0, 1, 0, 0, 1]
c_chr1=[10, 1, 2, 0, 0, 0, 1, 0, 0, 0]
g_chr1=[ 0, 2, 0, 2, 3, 1, 0, 2, 7, 0]
a_chr2=[ 0, 0, 4, 3, 2, 4, 6, 1, 2, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0]
t_chr2=[ 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 3, 1, 3, 2, 2, 1, 1, 2]
c_chr2=[ 0, 0, 1, 2, 1, 1, 2, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
g_chr2=[ 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1]
a_chr3=[ 0, 0, 1, 0, 0, 0,15, 0, 0, 0, 0, 0, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_chr3=[ 1, 0, 0, 1, 0,17, 1, 0,10,11, 0, 2, 3, 1, 3, 2,12,10, 1, 6, 1, 0, 0, 1, 0, 0, 8, 5, 0, 1]
c_chr3=[ 0, 0,19,21,18, 1, 2, 10, 0, 0,12, 1, 0, 0, 0,13, 0, 0,11, 0, 1, 0, 1, 2, 1,20, 0, 1, 0, 0]
g_chr3=[20,20, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 2, 3, 1, 1, 2, 7, 0]   
gene_locs = (["NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA","NA"])+ \
            (["NA","NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA",
              "NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA","NA"])+ \
            (["NA","NA","NA","gene1","gene1","gene1","gene1","NA","NA","NA",
              "NA","gene2","gene2","gene2","gene2","gene2","NA","NA","NA","NA",
              "NA","NA","gene3","gene3","gene3","gene3","gene3","NA","NA","NA"])

@pytest.fixture
def profile_1()->pl.LazyFrame:
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "genome": ["genome1"]*30 + ["genome2"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()

@pytest.fixture
def scaffold_1()->str:
    return "chr1\nchr2\nchr3"

@pytest.fixture
def profile_2()->pl.LazyFrame:
    """Exactly the same as profile_1"""
    return pl.DataFrame({
        "chrom": ["chr1"]*10+["chr2"]*20+["chr3"]*30,
        "genome": ["genome1"]*30 + ["genome2"]*30,
        "pos":list(range(0,len(a_chr1)))+list(range(0,len(a_chr2)))+list(range(0,len(a_chr3))),
        "gene": gene_locs,
        "A": a_chr1 + a_chr2 + a_chr3,
        "T": t_chr1 + t_chr2 + t_chr3,
        "C": c_chr1 + c_chr2 + c_chr3,
        "G": g_chr1 + g_chr2 + g_chr3,
    }).lazy()

@pytest.fixture
def scaffold_2():
    return "chr1\nchr3"

@pytest.fixture
def stb()->pl.LazyFrame:
    return pl.DataFrame({
        "scaffold":["chr1","chr2","chr3"],
        "genome":["genome1","genome1","genome2"],
    }).lazy()

@pytest.fixture
def large_profile_1()->pl.LazyFrame:

    chr1=pl.DataFrame({
        "chrom": ["chr1"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr2=pl.DataFrame({
        "chrom": ["chr2"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr3=pl.DataFrame({
        "chrom": ["chr3"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    return pl.concat([chr1,chr2,chr3])
@pytest.fixture
def large_profile_2()->pl.LazyFrame:
    """A slightly different version of large_profile_1"""
    chr1=pl.DataFrame({
        "chrom": ["chr1"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr2=pl.DataFrame({
        "chrom": ["chr2"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    chr3=pl.DataFrame({
        "chrom": ["chr3"]*1000000,
        "pos":list(range(1000000)),
        "gene": ["NA"]*1000000,
        "A": [random.randint(0,10) for _ in range(1000000)],
        "T": [random.randint(0,10) for _ in range(1000000)],
        "C": [random.randint(0,10) for _ in range(1000000)],
        "G": [random.randint(0,10) for _ in range(1000000)],
    }).lazy()
    return pl.concat([chr1,chr2,chr3])
    
@pytest.mark.parametrize("min_cov,min_gene_compare_len", [(1, 1), (5, 1), (1, 3), (5, 3),(5,5)])
def test_compare_profiles_profile_1_2_mc_mgcl(profile_1,profile_2,min_cov,min_gene_compare_len):
    res_dict=compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=min_cov ,
        min_gene_compare_len=min_gene_compare_len,
        calculate="all",
    ).collect().fill_null(-1).rows_by_key(key="genome",unique=True,named=True)

    a_genome_1=a_chr1 + a_chr2 
    t_genome_1=t_chr1 + t_chr2 
    c_genome_1=c_chr1 + c_chr2 
    g_genome_1=g_chr1 + g_chr2 

    assert res_dict["genome1"]["total_positions"]==len([i for i in zip(a_genome_1,t_genome_1,c_genome_1,g_genome_1) if sum(i)>=min_cov])

    a_genome_2=a_chr3
    t_genome_2=t_chr3
    c_genome_2=c_chr3
    g_genome_2=g_chr3

    assert res_dict["genome2"]["total_positions"]==len([i for i in zip(a_genome_2,t_genome_2,c_genome_2,g_genome_2) if sum(i)>=min_cov])

    assert res_dict["genome1"]["share_allele_pos"]==len([i for i in zip(a_genome_1,t_genome_1,c_genome_1,g_genome_1) if sum(i)>=min_cov])

    assert res_dict["genome2"]["share_allele_pos"]==len([i for i in zip(a_genome_2,t_genome_2,c_genome_2,g_genome_2) if sum(i)>=min_cov])

    assert res_dict["genome1"]["genome_ani"]==100.0
    
    assert res_dict["genome2"]["genome_ani"]==100.0

    assert res_dict["genome1"]["max_consecutive_length"]==max([len([i for i in zip(a_chr1,t_chr1,c_chr1,g_chr1) if sum(i)>=min_cov]),len([i for i in zip(a_chr2,t_chr2,c_chr2,g_chr2) if sum(i)>=min_cov])])

    assert res_dict["genome2"]["max_consecutive_length"]==max([len([i for i in zip(a_chr3,t_chr3,c_chr3,g_chr3) if sum(i)>=min_cov])])


    # Without a gene range table "all" resolves to genome-level metrics only, so
    # the table stays genome-grained and carries no gene columns.
    assert "gene" not in res_dict["genome1"]
    assert "gene_ani" not in res_dict["genome1"]

def test_con_ani_expression():
    profile_1=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [10,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [0,0,5],
    }).lazy()
    profile_2=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [5,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [5,0,5],
    }).lazy()
    profile_3=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [0,0,0],
        "T": [10,10,0],
        "C": [0,0,0],
        "G": [0,0,10],
    }).lazy()

    assert compare._shared_loci_polars(profile_1, profile_2, min_cov=1, ani_method="conani").select(pl.col("surr")).sum().collect()[0,0]==3
    assert compare._shared_loci_polars(profile_3, profile_2, min_cov=1, ani_method="conani").select(pl.col("surr")).sum().collect()[0,0]==2

@pytest.mark.parametrize("threshold1,threshold2", [(0.1,0.8),(0.2,0.9),(0.3,0.95),(0.4,0.99),(0.1,0.2),(0.5,0.6),(0.7,0.8),(0,1)])
def test_cos_ani_expression(threshold1,threshold2):
    profile_1=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [10,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [0,0,5],
    }).lazy()
    profile_2=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [5,0,5],
        "T": [0,10,0],
        "C": [0,0,0],
        "G": [5,0,5],
    }).lazy()
    profile_3=pl.DataFrame({
        "chrom": ["chr1"]*3,
        "genome": ["genome1"]*3,
        "pos":[0,1,2],
        "gene": ["gene1","gene1","gene1"],
        "A": [0,0,0],
        "T": [10,10,0],
        "C": [0,0,0],
        "G": [0,0,10],
    }).lazy()

    assert compare._shared_loci_polars(profile_1, profile_2, min_cov=1, ani_method=f"cosani_{threshold1}").select(pl.col("surr")).sum().collect()[0,0]>=compare._shared_loci_polars(profile_3, profile_2, min_cov=1, ani_method=f"cosani_{threshold2}").select(pl.col("surr")).sum().collect()[0,0]


def test_parse_multiple_ani_methods_normalizes_and_deduplicates():
    assert compare.parse_ani_methods(" popani,conani,popani,cosani_0.4 ") == (
        "popani",
        "conani",
        "cosani_0.4",
    )
    assert compare.canonical_ani_methods("conani, popani") == "conani,popani"
    with pytest.raises(ValueError, match="between 0 and 1"):
        compare.parse_ani_methods("popani,cosani_1.1")


@pytest.fixture
def multi_ani_profiles(tmp_path):
    common = {
        "chrom": ["s1"] * 4,
        "genome": ["genome1"] * 4,
        "pos": [1, 2, 3, 4],
    }
    p1 = pl.DataFrame(
        {
            **common,
            "A": [10, 6, 5, 0],
            "C": [0, 4, 5, 0],
            "G": [0, 0, 0, 0],
            "T": [0, 0, 0, 10],
        }
    )
    p2 = pl.DataFrame(
        {
            **common,
            "A": [10, 4, 0, 0],
            "C": [0, 6, 0, 0],
            "G": [0, 0, 10, 0],
            "T": [0, 0, 0, 10],
        }
    )
    p1_path = tmp_path / "multi_1.parquet"
    p2_path = tmp_path / "multi_2.parquet"
    p1.write_parquet(p1_path)
    p2.write_parquet(p2_path)
    return p1_path, p2_path


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_multiple_ani_methods_match_independent_runs(multi_ani_profiles, engine):
    p1, p2 = multi_ani_profiles
    methods = ("popani", "conani", "cosani_0.95")
    combined = compare.compare_genomes(
        p1,
        p2,
        min_cov=1,
        calculate="all",
        ani_method=",".join(methods),
        engine=engine,
    ).collect()

    assert combined.columns == compare.genome_metric_output_columns("all", methods)
    for method in methods:
        single = compare.compare_genomes(
            p1,
            p2,
            min_cov=1,
            calculate="all",
            ani_method=method,
            engine=engine,
        ).collect()
        suffix = compare.ani_method_suffix(method)
        assert combined[f"share_allele_pos_{suffix}"].to_list() == single[
            "share_allele_pos"
        ].to_list()
        assert combined[f"genome_ani_{suffix}"].to_list() == single[
            "genome_ani"
        ].to_list()
        assert combined[f"max_consecutive_length_{suffix}"].to_list() == single[
            "max_consecutive_length"
        ].to_list()


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_multiple_ani_methods_support_gene_output(
    multi_ani_profiles, engine, tmp_path
):
    p1, p2 = multi_ani_profiles
    ranges = tmp_path / "multi_genes.tsv"
    pl.DataFrame(
        {
            "gene": ["gene1"],
            "scaffold": ["s1"],
            "start": [1],
            "end": [4],
        }
    ).write_csv(ranges, separator="\t", include_header=False)

    out = compare.compare_genomes(
        p1,
        p2,
        min_cov=1,
        min_gene_compare_len=1,
        calculate="all",
        gene_range=ranges,
        ani_method="popani,conani",
        engine=engine,
    ).collect()

    assert out.columns == compare.gene_metric_output_columns(
        "all", "popani,conani"
    )
    assert out["share_allele_pos_popani"].to_list() == [3]
    assert out["share_allele_pos_conani"].to_list() == [2]
    assert out["gene_ani_popani"].to_list() == [75.0]
    assert out["gene_ani_conani"].to_list() == [50.0]


def test_duckdb_compare_genomes_to_parquet(profile_1, profile_2, stb, tmp_path):
    p1 = tmp_path / "p1.parquet"
    p2 = tmp_path / "p2.parquet"
    stb_path = tmp_path / "stb.tsv"
    out_path = tmp_path / "genome_comp.parquet"
    profile_1.sink_parquet(p1)
    profile_2.sink_parquet(p2)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    compare.duckdb_compare_genomes_to_parquet(
        mpile1=p1,
        mpile2=p2,
        output_file=out_path,
        stb_file=stb_path,
        sample_1_name="s1",
        sample_2_name="s2",
        min_cov=5,
        min_gene_compare_len=3,
        genome_scope="all",
        ani_method="popani",
        calculate="all",
        memory_limit="512MB",
        temp_directory=tmp_path,
    )

    out = pl.read_parquet(out_path)
    assert out.shape[0] == 2
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_ani",
        "max_consecutive_length",
        "sample_1",
        "sample_2",
    }
    assert set(out["sample_1"].to_list()) == {"s1"}
    assert set(out["sample_2"].to_list()) == {"s2"}


def test_polars_prefilter_by_scope_filters_genome_and_gene(profile_1):
    p1, p2 = compare.polars_prefilter_by_scope(
        profile_1,
        profile_1,
        genome_scope="genome1",
        gene_scope="gene2",
    )
    out1 = p1.collect()
    out2 = p2.collect()

    assert out1.height > 0
    assert out2.height > 0
    assert out1.get_column("genome").unique().to_list() == ["genome1"]
    assert out2.get_column("genome").unique().to_list() == ["genome1"]
    assert out1.get_column("gene").unique().to_list() == ["gene2"]
    assert out2.get_column("gene").unique().to_list() == ["gene2"]


def test_profile_is_coordinate_sorted_uses_parquet_metadata(tmp_path):
    profile_df = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "genome": ["genome1", "genome1"],
            "pos": [1, 2],
            "gene": ["gene1", "gene1"],
            "A": [5, 0],
            "T": [0, 5],
            "C": [0, 0],
            "G": [0, 0],
        }
    )
    sorted_path = tmp_path / "sorted_profile.parquet"
    legacy_path = tmp_path / "legacy_profile.parquet"

    profile_df.lazy().sink_parquet(
        sorted_path,
        metadata={
            compare.PROFILE_SORTED_METADATA_KEY: compare.PROFILE_SORTED_METADATA_VALUE,
        },
    )
    profile_df.write_parquet(legacy_path)

    assert compare._profile_is_coordinate_sorted(sorted_path) is True
    assert compare._profile_is_coordinate_sorted(legacy_path) is False
    assert compare._profile_is_coordinate_sorted(profile_df.lazy()) is False


def test_shared_loci_polars_sets_sorted_only_when_both_inputs_are_marked(tmp_path, monkeypatch):
    profile_df = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "genome": ["genome1", "genome1"],
            "pos": [1, 2],
            "gene": ["gene1", "gene1"],
            "A": [5, 0],
            "T": [0, 5],
            "C": [0, 0],
            "G": [0, 0],
        }
    )
    sorted_1 = tmp_path / "sorted_1.parquet"
    sorted_2 = tmp_path / "sorted_2.parquet"
    unmarked = tmp_path / "unmarked.parquet"
    metadata = {
        compare.PROFILE_SORTED_METADATA_KEY: compare.PROFILE_SORTED_METADATA_VALUE,
    }
    profile_df.lazy().sink_parquet(sorted_1, metadata=metadata)
    profile_df.lazy().sink_parquet(sorted_2, metadata=metadata)
    profile_df.write_parquet(unmarked)

    original_set_sorted = pl.LazyFrame.set_sorted
    calls: list[tuple[str, ...]] = []

    def _wrapped_set_sorted(self, by, *args, **kwargs):
        if isinstance(by, str):
            calls.append((by,))
        else:
            calls.append(tuple(by))
        return original_set_sorted(self, by, *args, **kwargs)

    monkeypatch.setattr(pl.LazyFrame, "set_sorted", _wrapped_set_sorted)

    compare._shared_loci_polars(sorted_1, sorted_2, min_cov=1).collect()
    assert calls == [("chrom", "pos"), ("chrom", "pos")]

    calls.clear()
    compare._shared_loci_polars(sorted_1, unmarked, min_cov=1).collect()
    assert calls == []


def test_compare_genomes_polars_matches_between_sorted_marked_and_unsorted_unmarked(profile_1, profile_2, tmp_path):
    sorted_left = tmp_path / "sorted_left.parquet"
    sorted_right = tmp_path / "sorted_right.parquet"
    unsorted_left = tmp_path / "unsorted_left.parquet"
    unsorted_right = tmp_path / "unsorted_right.parquet"
    metadata = {
        compare.PROFILE_SORTED_METADATA_KEY: compare.PROFILE_SORTED_METADATA_VALUE,
    }

    profile_1_df = profile_1.collect().sort(["chrom", "pos"])
    profile_2_df = profile_2.collect().sort(["chrom", "pos"])
    profile_1_df.lazy().sink_parquet(sorted_left, metadata=metadata)
    profile_2_df.lazy().sink_parquet(sorted_right, metadata=metadata)

    profile_1_df.reverse().write_parquet(unsorted_left)
    profile_2_df.reverse().write_parquet(unsorted_right)

    sorted_result = (
        compare.compare_genomes(
            mpile_contig_1=sorted_left,
            mpile_contig_2=sorted_right,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )
    unsorted_result = (
        compare.compare_genomes(
            mpile_contig_1=unsorted_left,
            mpile_contig_2=unsorted_right,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            calculate="all",
        )
        .collect()
        .sort("genome")
    )

    assert unsorted_result.equals(sorted_result)


def test_duckdb_compare_genomes_to_parquet_ani_only(profile_1, profile_2, stb, tmp_path):
    p1 = tmp_path / "p1.parquet"
    p2 = tmp_path / "p2.parquet"
    stb_path = tmp_path / "stb.tsv"
    out_path = tmp_path / "genome_comp_ani.parquet"
    profile_1.sink_parquet(p1)
    profile_2.sink_parquet(p2)
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    compare.duckdb_compare_genomes_to_parquet(
        mpile1=p1,
        mpile2=p2,
        output_file=out_path,
        stb_file=stb_path,
        sample_1_name="s1",
        sample_2_name="s2",
        min_cov=5,
        min_gene_compare_len=3,
        genome_scope="all",
        ani_method="popani",
        calculate="ani",
        memory_limit="512MB",
        temp_directory=tmp_path,
    )

    out = pl.read_parquet(out_path)
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_ani",
        "sample_1",
        "sample_2",
    }


def test_compare_genomes_ani_only_matches_full_subset(profile_1, profile_2, stb, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    stb.sink_csv(stb_path, separator="\t", include_header=False)

    full = (
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            stb_file=stb_path,
            calculate="all",
        )
        .collect()
        .select("genome", "total_positions", "share_allele_pos", "genome_ani")
        .sort("genome")
    )
    ani_only = (
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=5,
            min_gene_compare_len=3,
            engine="polars",
            stb_file=stb_path,
            calculate="ani",
        )
        .collect()
        .sort("genome")
    )
    assert ani_only.equals(full)





def test_compare_genomes_default_calculate_is_all_metrics(profile_1, profile_2, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    pl.DataFrame({"scaffold": ["chr1", "chr2", "chr3"], "genome": ["genome1", "genome1", "genome2"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )
    out = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="polars",
        stb_file=stb_path,
    ).collect()
    assert set(out.columns) == {
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_ani",
        "max_consecutive_length",
    }


def test_compare_genomes_calculate_all_matches_between_engines(profile_1, profile_2, tmp_path):
    stb_path = tmp_path / "stb.tsv"
    pl.DataFrame({"scaffold": ["chr1", "chr2", "chr3"], "genome": ["genome1", "genome1", "genome2"]}).write_csv(
        stb_path, separator="\t", include_header=False
    )
    out_polars = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="polars",
        calculate="all",
        stb_file=stb_path,
    ).collect().sort("genome")
    out_duckdb = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=5,
        min_gene_compare_len=3,
        engine="duckdb",
        calculate="all",
        stb_file=stb_path,
    ).collect().sort("genome")
    assert out_polars.equals(out_duckdb)


# --- gene-grained comparison over gene ranges -------------------------------
#
# Gene ANI is computed from a gene range table rather than a per-position gene
# label, so a position covered by several genes counts toward every one of them.
# The label-based approach could only ever attribute a position to one gene.


@pytest.fixture
def overlap_profiles(tmp_path):
    """Two profiles over one scaffold; every 7th position is an SNV."""
    n = 300
    positions = list(range(1, n + 1))
    alt = [10 if (p - 1) % 7 == 0 else 0 for p in positions]
    base = {"chrom": ["s1"] * n, "genome": ["genome1"] * n, "pos": positions}
    p1 = pl.DataFrame({**base, "A": [10] * n, "C": [0] * n, "G": [0] * n, "T": [0] * n})
    p2 = pl.DataFrame(
        {**base, "A": [10 - a for a in alt], "C": alt, "G": [0] * n, "T": [0] * n}
    )
    path_1 = tmp_path / "ov1.parquet"
    path_2 = tmp_path / "ov2.parquet"
    p1.write_parquet(path_1)
    p2.write_parquet(path_2)
    return path_1, path_2


@pytest.fixture
def overlap_gene_range(tmp_path):
    """gene_b and gene_c sit strictly inside gene_a."""
    path = tmp_path / "genes.tsv"
    pl.DataFrame(
        {
            "gene": ["gene_a", "gene_b", "gene_c"],
            "scaffold": ["s1", "s1", "s1"],
            "start": [10, 50, 200],
            "end": [280, 90, 260],
        }
    ).write_csv(path, separator="\t", include_header=False)
    return path


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_nested_genes_each_get_their_full_range(
    overlap_profiles, overlap_gene_range, engine
):
    profile_1, profile_2 = overlap_profiles
    out = (
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            min_cov=1,
            min_gene_compare_len=1,
            engine=engine,
            calculate="all",
            gene_range=overlap_gene_range,
        )
        .collect()
        .sort("gene")
        .rows_by_key(key="gene", unique=True, named=True)
    )
    # Every position in [start, end] belongs to the gene, including those the
    # enclosing gene shares with a nested one.
    assert out["gene_a"]["total_positions"] == 280 - 10 + 1
    assert out["gene_b"]["total_positions"] == 90 - 50 + 1
    assert out["gene_c"]["total_positions"] == 260 - 200 + 1


@pytest.mark.parametrize("engine", ["polars", "duckdb"])
def test_gene_grained_output_carries_genome_metrics(
    overlap_profiles, overlap_gene_range, engine
):
    profile_1, profile_2 = overlap_profiles
    out = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=1,
        min_gene_compare_len=1,
        engine=engine,
        calculate="all",
        gene_range=overlap_gene_range,
    ).collect()
    assert set(out.columns) == {
        "genome",
        "gene",
        "total_positions",
        "share_allele_pos",
        "gene_ani",
        "genome_ani",
        "max_consecutive_length",
    }
    # Genome-level values are repeated on each gene row, so they must be constant.
    assert out["genome_ani"].n_unique() == 1
    assert out["max_consecutive_length"].n_unique() == 1


def test_gene_grained_output_matches_between_engines(
    overlap_profiles, overlap_gene_range
):
    profile_1, profile_2 = overlap_profiles
    kwargs = dict(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=1,
        min_gene_compare_len=1,
        calculate="all",
        gene_range=overlap_gene_range,
    )
    out_polars = compare.compare_genomes(engine="polars", **kwargs).collect().sort("gene")
    out_duckdb = compare.compare_genomes(engine="duckdb", **kwargs).collect().sort("gene")
    assert out_polars.equals(out_duckdb)


def test_gene_calculation_without_ranges_is_rejected(overlap_profiles):
    profile_1, profile_2 = overlap_profiles
    with pytest.raises(ValueError, match="gene range table"):
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            engine="polars",
            calculate="gene",
        ).collect()


def test_all_without_gene_ranges_stays_genome_grained(overlap_profiles):
    profile_1, profile_2 = overlap_profiles
    out = compare.compare_genomes(
        mpile_contig_1=profile_1,
        mpile_contig_2=profile_2,
        min_cov=1,
        engine="polars",
        calculate="all",
    ).collect()
    assert "gene" not in out.columns
    assert "genome_ani" in out.columns


def test_identical_genes_calculation_is_removed(overlap_profiles):
    profile_1, profile_2 = overlap_profiles
    with pytest.raises(ValueError, match="Unsupported calculation"):
        compare.compare_genomes(
            mpile_contig_1=profile_1,
            mpile_contig_2=profile_2,
            engine="polars",
            calculate="identical_genes",
        ).collect()
