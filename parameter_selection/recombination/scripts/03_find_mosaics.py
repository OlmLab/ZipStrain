"""Find mosaic regions and identical stretches between pairs of ZipStrain profiles.

For each comparison this writes, under results/<comparison>/:
  mosaic_regions.tsv      blocks of unusually low or high divergence (PairComparison.mosaic_regions)
  identical_stretches.tsv identical stretches of at least 1 kb, with p-values
  windows.tsv             shared positions and mismatches in 10 kb windows (for plotting)
  mismatches.tsv          every shared position where the two samples do not match
and results/evaluation.tsv compares the A-vs-B blocks with the planted truth.
"""

import argparse
from pathlib import Path

import polars as pl
from zipstrain.compare import PairComparison

COMPARISONS = [
    # name,               sample 1,   sample 2,   ANI definition of a "match"
    ("A_vs_B_popani", "sample_A", "sample_B", "popani"),
    ("A_vs_M_popani", "sample_A", "sample_M", "popani"),
    ("A_vs_M_conani", "sample_A", "sample_M", "conani"),
]
WINDOW = 10_000


def profile(profiles: Path, sample: str) -> Path:
    return profiles / sample / f"{sample}_profile.parquet"


def analyse(pair: PairComparison):
    """Run every analysis of one pair; Polars shares the profile join between them."""
    matches = pair.matches
    windows = (matches.with_columns(window_start=(pl.col("pos") - 1) // WINDOW * WINDOW + 1)
               .group_by("genome", "chrom", "window_start")
               .agg(pl.len().alias("shared_positions"), (~pl.col("match")).sum().alias("mismatches")))
    mismatches = matches.filter(~pl.col("match")).select("genome", "chrom", "pos")
    mosaics, stretches, windows, mismatches = pair.collect(
        pair.mosaic_regions(), pair.identical_stretches(min_length=1_000), windows, mismatches)
    # Sort after collecting: the analyses run their own NumPy passes when collected.
    return (mosaics.sort("genome", "chrom", "start"),
            stretches.sort("length", descending=True),
            windows.sort("genome", "chrom", "window_start"),
            mismatches.sort("genome", "chrom", "pos"))


def evaluate(truth: pl.DataFrame, mosaics: pl.DataFrame) -> pl.DataFrame:
    """Match each planted block to the detected block of the same kind that overlaps it most."""
    rows = []
    for block in truth.iter_rows(named=True):
        hits = mosaics.filter((pl.col("state") == block["kind"]) & (pl.col("chrom") == block["chrom"]) &
                              (pl.col("start") <= block["end"]) & (pl.col("end") >= block["start"]))
        row = {key: block[key] for key in ("block", "kind", "start", "end", "length")}
        if hits.is_empty():
            rows.append({**row, "detected": False, "detected_start": None, "detected_end": None,
                         "start_offset": None, "end_offset": None, "log_likelihood_ratio": None})
            continue
        overlap = (pl.min_horizontal(pl.col("end"), block["end"]) -
                   pl.max_horizontal(pl.col("start"), block["start"]))
        best = hits.sort(overlap, descending=True).row(0, named=True)
        rows.append({**row, "detected": True, "detected_start": best["start"], "detected_end": best["end"],
                     "start_offset": best["start"] - block["start"], "end_offset": best["end"] - block["end"],
                     "log_likelihood_ratio": round(best["log_likelihood_ratio"], 1)})
    return pl.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profiles", default="profiles", type=Path)
    parser.add_argument("--truth", default="data/truth/planted_blocks.tsv", type=Path)
    parser.add_argument("--outdir", default="results", type=Path)
    args = parser.parse_args()

    for name, sample_1, sample_2, ani_method in COMPARISONS:
        pair = PairComparison(profile(args.profiles, sample_1), profile(args.profiles, sample_2),
                              ani_method=ani_method)
        mosaics, stretches, windows, mismatches = analyse(pair)
        outdir = args.outdir / name
        outdir.mkdir(parents=True, exist_ok=True)
        mosaics.write_csv(outdir / "mosaic_regions.tsv", separator="\t")
        stretches.write_csv(outdir / "identical_stretches.tsv", separator="\t")
        windows.write_csv(outdir / "windows.tsv", separator="\t")
        mismatches.write_csv(outdir / "mismatches.tsv", separator="\t")
        shared = int(windows["shared_positions"].sum())
        print(f"{name}: {shared:,} shared positions, {mismatches.height:,} mismatches, "
              f"{mosaics.height} mosaic blocks, {stretches.height} stretches >= 1 kb")
        if name == "A_vs_B_popani":
            evaluation = evaluate(pl.read_csv(args.truth, separator="\t"), mosaics)

    evaluation.write_csv(args.outdir / "evaluation.tsv", separator="\t")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200, tbl_hide_dataframe_shape=True):
        print(evaluation)


if __name__ == "__main__":
    main()
