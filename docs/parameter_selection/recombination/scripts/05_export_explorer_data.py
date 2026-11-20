"""Export the A-vs-B comparison for the interactive explorer (interactive/mosaic_explorer.html).

Writes interactive/tutorial_data.js with every mismatch between samples A and B, its rank
along the shared positions (the number of shared positions before it), 10 kb window counts and
the planted blocks, so the page can rerun mosaic_regions on exactly the tutorial's data.
"""

import argparse
import json
from pathlib import Path

import polars as pl
from zipstrain.compare import PairComparison

WINDOW = 10_000


def fasta_length(path: Path) -> int:
    with open(path) as handle:
        return sum(len(line.strip()) for line in handle if not line.startswith(">"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profiles", default="profiles", type=Path)
    parser.add_argument("--truth", default="data/truth/planted_blocks.tsv", type=Path)
    parser.add_argument("--reference", default="data/strain_A.fna", type=Path)
    parser.add_argument("--out", default="interactive/tutorial_data.js", type=Path)
    args = parser.parse_args()

    pair = PairComparison(args.profiles / "sample_A" / "sample_A_profile.parquet",
                          args.profiles / "sample_B" / "sample_B_profile.parquet", ani_method="popani")
    shared = (pair.collect(pair.matches)[0].sort("pos")
              .with_row_index("rank").with_columns(pl.col("rank").cast(pl.Int64)))
    if shared["chrom"].n_unique() != 1:
        raise SystemExit("The explorer expects a single chromosome.")
    mismatches = shared.filter(~pl.col("match"))
    windows = (shared.with_columns(window=(pl.col("pos") - 1) // WINDOW)
               .group_by("window").agg(pl.len().alias("shared"), (~pl.col("match")).sum().alias("mismatches"))
               .sort("window"))
    n_windows = (fasta_length(args.reference) + WINDOW - 1) // WINDOW
    counts = {row["window"]: row for row in windows.iter_rows(named=True)}
    truth = pl.read_csv(args.truth, separator="\t")
    data = {
        "genome_length": fasta_length(args.reference),
        "window": WINDOW,
        "first": {"rank": -1, "pos": int(shared["pos"][0])},
        "last": {"rank": shared.height - 1, "pos": int(shared["pos"][-1])},
        "mismatch_pos": mismatches["pos"].to_list(),
        "mismatch_rank": mismatches["rank"].to_list(),
        "window_shared": [counts[i]["shared"] if i in counts else 0 for i in range(n_windows)],
        "window_mismatches": [counts[i]["mismatches"] if i in counts else 0 for i in range(n_windows)],
        "planted": [{"name": r["block"], "kind": r["kind"], "start": r["start"], "end": r["end"],
                     "rate": r["substitution_rate"]} for r in truth.iter_rows(named=True)],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("// Written by scripts/05_export_explorer_data.py\nwindow.TUTORIAL_DATA = "
                        + json.dumps(data, separators=(",", ":")) + ";\n")
    print(f"{args.out}: {shared.height:,} shared positions, {mismatches.height:,} mismatches")


if __name__ == "__main__":
    main()
