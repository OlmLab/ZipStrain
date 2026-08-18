from pathlib import Path

import polars as pl
import math

import pytest
from click.testing import CliRunner

from zipstrain import cli, dnds, profile


def _gene_info(path: Path) -> Path:
    pl.DataFrame(
        {
            "gene_id": pl.Series([0, 1], dtype=pl.UInt32),
            "gene": ["plus", "minus"],
            "genome": ["genome", "genome"],
            "scaffold": ["plus_scaffold", "minus_scaffold"],
            "start": [1, 1],
            "end": [6, 3],
            "strand": pl.Series([1, -1], dtype=pl.Int8),
            "phase": pl.Series([0, 0], dtype=pl.UInt8),
            "genetic_code": pl.Series([11, 11], dtype=pl.UInt8),
            "partial_5p": [False, False],
            "partial_3p": [False, False],
            "first_codon_id": pl.Series([0, 2], dtype=pl.UInt64),
            "n_codons": pl.Series([2, 1], dtype=pl.UInt32),
        }
    ).write_parquet(path)
    return path


def _profile(path: Path, *, sample_2: bool = False) -> Path:
    # plus reference: ATG GAA. minus reference genomic CAT -> coding ATG.
    plus = "ATAGAG" if sample_2 else "ATGGAA"
    minus = "CAA" if sample_2 else "CAT"
    reference = "ATGGAACAT"
    rows = []
    for scaffold, sample, ref in (
        ("plus_scaffold", plus, reference[:6]),
        ("minus_scaffold", minus, reference[6:]),
    ):
        for pos, (base, ref_base) in enumerate(zip(sample, ref), start=1):
            counts = {value: 0 for value in "ACGT"}
            counts[base] = 10
            rows.append(
                {
                    "chrom": scaffold,
                    "genome": "genome",
                    "pos": pos,
                    **counts,
                    "ref_base_bitmask": {"A": 1, "C": 2, "G": 4, "T": 8}[ref_base],
                }
            )
    pl.DataFrame(rows).write_parquet(path)
    return path


def _set_counts(
    path: Path,
    *,
    scaffold: str,
    pos: int,
    A: int = 0,
    C: int = 0,
    G: int = 0,
    T: int = 0,
) -> Path:
    at_site = (pl.col("chrom") == scaffold) & (pl.col("pos") == pos)
    pl.read_parquet(path).with_columns(
        *(
            pl.when(at_site).then(value).otherwise(pl.col(base)).alias(base)
            for base, value in {"A": A, "C": C, "G": G, "T": T}.items()
        )
    ).write_parquet(path)
    return path


def _sidecar_consensus_codon(row: dict) -> str:
    return "".join(
        max(BASES, key=lambda base: row[f"{base}_{offset}"])
        for offset in range(3)
    )


BASES = "ACGT"


def test_codon_lookup_classifies_synonymous_and_nonsynonymous_changes():
    lookup = dnds.codon_pair_lookup()

    atg_ata = dnds._encode_codon("ATG") * 64 + dnds._encode_codon("ATA")
    gaa_gag = dnds._encode_codon("GAA") * 64 + dnds._encode_codon("GAG")

    nonsynonymous = lookup.filter(pl.col("pair_code") == atg_ata).row(0, named=True)
    synonymous = lookup.filter(pl.col("pair_code") == gaa_gag).row(0, named=True)
    assert nonsynonymous["nonsynonymous_changes"] == 1.0
    assert nonsynonymous["synonymous_changes"] == 0.0
    assert synonymous["synonymous_changes"] == 1.0
    assert synonymous["nonsynonymous_changes"] == 0.0


def test_build_gene_info_table_records_coding_orientation(tmp_path):
    genes = tmp_path / "genes.fna"
    genes.write_text(
        ">contig_1 # 1 # 6 # 1 # partial=00;transl_table=11\nATGGAA\n"
        ">contig_2 # 2 # 7 # -1 # partial=10;transl_table=11\nTTTCAT\n"
    )
    stb = pl.DataFrame(
        {"scaffold": ["contig", "contig"], "genome": ["g", "g"]}
    ).lazy()
    result = profile.build_gene_info_table(genes, stb=stb)

    assert result.columns == [
        "gene_id", "gene", "genome", "scaffold", "start", "end", "strand",
        "phase", "genetic_code", "partial_5p", "partial_3p",
        "first_codon_id", "n_codons",
    ]
    assert result["strand"].to_list() == [1, -1]
    assert result["partial_5p"].to_list() == [False, True]
    assert result["first_codon_id"].to_list() == [0, 2]


def test_streaming_gene_info_writer_matches_in_memory_builder(tmp_path):
    genes = tmp_path / "genes.fna"
    genes.write_text(
        ">contig_2 # 7 # 12 # -1 # partial=10;transl_table=11\nTTTCAT\n"
        ">contig_1 # 1 # 6 # 1 # partial=00;transl_table=11\nATGGAA\n"
    )
    stb_path = tmp_path / "reference.stb"
    stb_path.write_text("contig\tgenome\n")
    output = tmp_path / "gene_info.parquet"

    profile.write_gene_info_table(genes, output, stb_file=stb_path)
    expected = profile.build_gene_info_table(
        genes,
        stb=profile.read_stb(stb_path),
    )

    assert pl.read_parquet(output).equals(expected)


def test_sparse_codon_profile_and_reference_dnds_handle_both_strands(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample = _profile(tmp_path / "sample_profile.parquet", sample_2=True)
    codons = tmp_path / "sample_codon_profile.parquet"

    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=genes,
        output_path=codons,
        temp_directory=tmp_path / "tmp",
    )
    codon_frame = pl.read_parquet(codons).sort("codon_id")
    assert [_sidecar_consensus_codon(row) for row in codon_frame.to_dicts()] == [
        "ATA", "GAG", "TTG",
    ]
    assert codon_frame["reference_codon_code"].to_list() == [
        dnds._encode_codon("ATG"),
        dnds._encode_codon("GAA"),
        dnds._encode_codon("ATG"),
    ]

    stats = dnds.reference_dnds(codons, genes, min_cov=5).collect().sort("gene")
    plus = stats.filter(pl.col("gene") == "plus").row(0, named=True)
    minus = stats.filter(pl.col("gene") == "minus").row(0, named=True)
    assert plus["ref_callable_codons"] == 2
    assert plus["ref_sns_count"] == 2
    assert plus["ref_snv_count"] == 0
    assert plus["ref_sns_synonymous_changes"] == 1.0
    assert plus["ref_sns_nonsynonymous_changes"] == 1.0
    assert minus["ref_callable_codons"] == 1
    assert minus["ref_sns_nonsynonymous_changes"] == 1.0


def test_pairwise_dnds_uses_sample_codons_and_is_symmetric(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample_1 = _profile(tmp_path / "one_profile.parquet")
    sample_2 = _profile(tmp_path / "two_profile.parquet", sample_2=True)
    codons_1 = tmp_path / "one_codon_profile.parquet"
    codons_2 = tmp_path / "two_codon_profile.parquet"
    for source, output in ((sample_1, codons_1), (sample_2, codons_2)):
        dnds.write_codon_profile(
            profile_path=source,
            gene_info_path=genes,
            output_path=output,
            temp_directory=tmp_path / output.stem,
        )

    forward = dnds.pairwise_dnds(codons_1, codons_2, genes).collect().sort("gene")
    reverse = dnds.pairwise_dnds(codons_2, codons_1, genes).collect().sort("gene")
    assert forward.select(dnds.DNDS_RESULT_COLUMNS).to_dicts() == reverse.select(
        dnds.DNDS_RESULT_COLUMNS
    ).to_dicts()
    plus = forward.filter(pl.col("gene") == "plus").row(0, named=True)
    assert plus["sns_count"] == 2
    assert plus["snv_count"] == 0
    assert plus["sns_synonymous_changes"] == 1.0
    assert plus["sns_nonsynonymous_changes"] == 1.0
    assert plus["pS"] == 0.0
    assert plus["pN"] == 0.0
    assert plus["pN_pS"] is None

    # dN/dS applies only to fixed SNS substitutions. Tiny genes can saturate
    # under the Jukes-Cantor correction at an observed proportion >= 3/4.
    def jukes_cantor(proportion: float) -> float | None:
        argument = 1.0 - (4.0 / 3.0) * proportion
        return -0.75 * math.log(argument) if argument > 0 else None

    sns_proportions = {
        "dS": plus["sns_synonymous_changes"] / plus["synonymous_sites"],
        "dN": plus["sns_nonsynonymous_changes"] / plus["nonsynonymous_sites"],
    }
    for corrected, proportion in sns_proportions.items():
        expected = jukes_cantor(proportion)
        if expected is None:
            assert proportion >= 0.75
            assert plus[corrected] is None
        else:
            assert plus[corrected] == pytest.approx(expected)

    if plus["dS"] is None or plus["dN"] is None or plus["dS"] == 0:
        assert plus["dN_dS"] is None
    else:
        assert plus["dN_dS"] == pytest.approx(plus["dN"] / plus["dS"])


def test_reference_consensus_averages_only_exact_top_ties(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample = _set_counts(
        _profile(tmp_path / "sample_profile.parquet"),
        scaffold="plus_scaffold",
        pos=6,
        G=5,
        T=5,
    )
    codons = dnds.codon_profile_path(sample)
    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=genes,
        output_path=codons,
        temp_directory=tmp_path / "tmp",
    )

    result = dnds.reference_dnds(
        codons, genes, min_cov=5, allele_integration="consensus"
    ).collect().filter(pl.col("gene") == "plus").row(0, named=True)

    assert result["ref_sns_count"] == 0
    assert result["ref_snv_count"] == 1
    assert result["ref_snv_synonymous_changes"] == pytest.approx(0.5)
    assert result["ref_snv_nonsynonymous_changes"] == pytest.approx(0.5)
    assert result["ref_allele_tie_sites"] == 1


def test_reference_weighted_uses_observed_allele_frequencies(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample = _set_counts(
        _profile(tmp_path / "sample_profile.parquet"),
        scaffold="plus_scaffold",
        pos=6,
        G=8,
        T=2,
    )
    codons = dnds.codon_profile_path(sample)
    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=genes,
        output_path=codons,
        temp_directory=tmp_path / "tmp",
    )

    consensus = dnds.reference_dnds(
        codons, genes, min_cov=5, allele_integration="consensus"
    ).collect().filter(pl.col("gene") == "plus").row(0, named=True)
    weighted = dnds.reference_dnds(
        codons, genes, min_cov=5, allele_integration="weighted"
    ).collect().filter(pl.col("gene") == "plus").row(0, named=True)

    assert consensus["ref_snv_synonymous_changes"] == 1.0
    assert consensus["ref_snv_nonsynonymous_changes"] == 0.0
    assert weighted["ref_snv_synonymous_changes"] == pytest.approx(0.8)
    assert weighted["ref_snv_nonsynonymous_changes"] == pytest.approx(0.2)


def test_pairwise_consensus_tie_is_symmetric(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample_1 = _profile(tmp_path / "one_profile.parquet")
    sample_2 = _set_counts(
        _profile(tmp_path / "two_profile.parquet"),
        scaffold="plus_scaffold",
        pos=6,
        G=5,
        T=5,
    )
    codons = []
    for sample in (sample_1, sample_2):
        output = dnds.codon_profile_path(sample)
        dnds.write_codon_profile(
            profile_path=sample,
            gene_info_path=genes,
            output_path=output,
            temp_directory=tmp_path / output.stem,
        )
        codons.append(output)

    forward = dnds.pairwise_dnds(
        codons[0], codons[1], genes, allele_integration="consensus"
    ).collect().sort("gene")
    reverse = dnds.pairwise_dnds(
        codons[1], codons[0], genes, allele_integration="consensus"
    ).collect().sort("gene")
    assert forward.select(dnds.DNDS_RESULT_COLUMNS).to_dicts() == reverse.select(
        dnds.DNDS_RESULT_COLUMNS
    ).to_dicts()
    plus = forward.filter(pl.col("gene") == "plus").row(0, named=True)
    assert plus["snv_count"] == 1
    assert plus["snv_synonymous_changes"] == pytest.approx(0.5)
    assert plus["snv_nonsynonymous_changes"] == pytest.approx(0.5)
    assert plus["allele_tie_sites"] == 1


def test_zero_divergence_reports_zero_not_negative_zero(tmp_path):
    """log(1) yields -0.0; the emitted distance should be a plain zero."""
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample = _profile(tmp_path / "sample_profile.parquet")
    codons = dnds.codon_profile_path(sample)
    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=genes,
        output_path=codons,
        temp_directory=tmp_path / "tmp",
    )
    result = (
        dnds.pairwise_dnds(codons, codons, genes)
        .collect()
        .filter(pl.col("gene") == "plus")
        .row(0, named=True)
    )
    assert result["pS"] == 0.0
    assert result["dS"] == 0.0
    assert math.copysign(1.0, result["dS"]) > 0


def test_pairwise_dnds_rejects_different_gene_contracts(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample = _profile(tmp_path / "sample_profile.parquet")
    codons_1 = tmp_path / "one_codon_profile.parquet"
    codons_2 = tmp_path / "two_codon_profile.parquet"
    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=genes,
        output_path=codons_1,
        temp_directory=tmp_path / "one",
    )
    changed_genes = pl.read_parquet(genes).with_columns(pl.lit("changed").alias("gene"))
    changed_path = tmp_path / "changed_gene_info.parquet"
    changed_genes.write_parquet(changed_path)
    dnds.write_codon_profile(
        profile_path=sample,
        gene_info_path=changed_path,
        output_path=codons_2,
        temp_directory=tmp_path / "two",
    )

    with pytest.raises(ValueError, match="incompatible"):
        dnds.pairwise_dnds(codons_1, codons_2, genes).collect()


def test_codon_profile_rejects_unsupported_genetic_code(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    pl.read_parquet(genes).with_columns(
        pl.lit(4, dtype=pl.UInt8).alias("genetic_code")
    ).write_parquet(genes)
    sample = _profile(tmp_path / "sample_profile.parquet")

    with pytest.raises(ValueError, match="supports genetic codes 1 and 11"):
        dnds.write_codon_profile(
            profile_path=sample,
            gene_info_path=genes,
            output_path=tmp_path / "sample_codon_profile.parquet",
            temp_directory=tmp_path / "tmp",
        )


def test_single_compare_cli_adds_pairwise_dnds_columns(tmp_path):
    genes = _gene_info(tmp_path / "gene_info.parquet")
    sample_1 = _profile(tmp_path / "one_profile.parquet")
    sample_2 = _profile(tmp_path / "two_profile.parquet", sample_2=True)
    for source in (sample_1, sample_2):
        dnds.write_codon_profile(
            profile_path=source,
            gene_info_path=genes,
            output_path=dnds.codon_profile_path(source),
            temp_directory=tmp_path / f"tmp_{source.stem}",
        )
    stb = tmp_path / "reference.stb"
    stb.write_text("plus_scaffold\tgenome\nminus_scaffold\tgenome\n")
    output = tmp_path / "comparison.parquet"

    result = CliRunner().invoke(
        cli.cli,
        [
            "utilities",
            "single-compare",
            "--profile-location-1", str(sample_1),
            "--profile-location-2", str(sample_2),
            "--stb-file", str(stb),
            "--gene-info-table", str(genes),
            "--calculate", "gene+dnds",
            "--allele-integration", "weighted",
            "--min-cov", "5",
            "--min-gene-compare-len", "1",
            "--output-file", str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    result_frame = pl.read_parquet(output)
    assert set(dnds.DNDS_RESULT_COLUMNS).issubset(result_frame.columns)
    plus = result_frame.filter(pl.col("gene") == "plus").row(0, named=True)
    assert plus["sns_synonymous_changes"] == 1.0
    assert plus["sns_nonsynonymous_changes"] == 1.0
    assert pl.read_parquet_metadata(output)[
        "zipstrain_compare_allele_integration"
    ] == "weighted"
