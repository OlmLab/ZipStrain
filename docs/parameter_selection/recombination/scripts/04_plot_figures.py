"""Draw the tutorial figures from the tables written by 03_find_mosaics.py.

  figures/fig1_genome_overview.png   planted vs detected blocks and windowed divergence
  figures/fig2_block_zooms.png       cumulative mismatches around each planted block
  figures/fig3_identical_stretches.png  stretch length vs significance
  figures/fig4_mixed_population.png  why popANI hides mosaics in a mixed sample
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

# ── Palette (validated: blue/orange pass CVD and contrast checks on the light surface) ──
SURFACE = "#fcfcfb"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
STATE_COLOR = {"low_divergence": "#2a78d6", "high_divergence": "#eb6834"}
STATE_LABEL = {"low_divergence": "Low divergence (recombination)",
               "high_divergence": "High divergence (import)"}
FACTOR = 20  # mosaic_regions default

plt.rcParams.update({
    "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "svg.fonttype": "none",
})


def read(path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t")


def panel_letter(ax, letter, pad=8):
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=13, pad=pad)


def rounded_bar(ax, start, width, y, height, *, color, hatch=False, radius_px=3.5, min_width_px=4):
    """Draw a horizontal bar with rounded ends in data coordinates (x and y may differ in scale)."""
    fig = ax.figure
    fig.canvas.draw_idle()
    box = ax.get_window_extent()
    x_per_px = np.diff(ax.get_xlim())[0] / box.width
    y_per_px = np.diff(ax.get_ylim())[0] / box.height
    if width < min_width_px * x_per_px:  # keep very short blocks visible, centred on their span
        start, width = start + width / 2 - min_width_px * x_per_px / 2, min_width_px * x_per_px
    radius = min(radius_px * x_per_px, width / 2)
    aspect = y_per_px / x_per_px
    radius = min(radius, height / 2 / aspect)
    patch = FancyBboxPatch((start, y - height / 2), width, height,
                           boxstyle=f"round,pad=0,rounding_size={radius}", mutation_aspect=aspect,
                           facecolor=SURFACE if hatch else color, edgecolor=color,
                           linewidth=1.2 if hatch else 0, hatch="////" if hatch else None, zorder=3)
    if hatch:
        patch.set_hatch_linewidth(0.8) if hasattr(patch, "set_hatch_linewidth") else None
    ax.add_patch(patch)


def state_rates(background):
    return {"low_divergence": 1 - (1 - background) ** (1 / FACTOR), "background": background,
            "high_divergence": 1 - (1 - background) ** FACTOR}


def figure_1(truth, results, out):
    mosaics = read(results / "A_vs_B_popani" / "mosaic_regions.tsv")
    windows = read(results / "A_vs_B_popani" / "windows.tsv")
    background = mosaics["background_divergence"][0]
    rates = state_rates(background)
    genome_mb = (windows["window_start"].max() + windows["shared_positions"][-1]) / 1e6

    fig = plt.figure(figsize=(12, 5.8))
    grid = fig.add_gridspec(2, 1, height_ratios=[1.15, 3], hspace=0.32, left=0.085, right=0.86, top=0.86, bottom=0.11)
    tracks, curve = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])

    # ── Tracks: planted (hatched) and detected (solid) ──
    tracks.set_xlim(0, genome_mb)
    tracks.set_ylim(-0.6, 1.75)
    tracks.set_yticks([1, 0], ["Planted", "Detected"])
    tracks.tick_params(axis="y", length=0, labelsize=10, labelcolor=INK)
    tracks.tick_params(axis="x", labelbottom=False, length=0)
    tracks.spines[["left", "bottom"]].set_visible(False)
    for y in (1, 0):
        tracks.axhline(y, color=GRID, linewidth=1, zorder=1)
    for row in truth.iter_rows(named=True):
        rounded_bar(tracks, row["start"] / 1e6, row["length"] / 1e6, 1, 0.5,
                    color=STATE_COLOR[row["kind"]], hatch=True)
        tracks.text((row["start"] + row["length"] / 2) / 1e6, 1.42, row["block"], ha="center",
                    va="center", fontsize=9, color=INK, fontweight="bold")
    for row in mosaics.iter_rows(named=True):
        rounded_bar(tracks, row["start"] / 1e6, (row["end"] - row["start"] + 1) / 1e6, 0, 0.5,
                    color=STATE_COLOR[row["state"]])
    panel_letter(tracks, "a")

    # ── Divergence in 10 kb windows on a log scale; the three model states are equally spaced ──
    floor = 3e-5
    divergence = (windows["mismatches"] / windows["shared_positions"]).to_numpy()
    starts = windows["window_start"].to_numpy() / 1e6
    width = (windows["shared_positions"].to_numpy() / 1e6)
    shown = np.where(divergence > 0, divergence, floor)
    for row in mosaics.iter_rows(named=True):
        curve.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=STATE_COLOR[row["state"]], alpha=0.13,
                      linewidth=0, zorder=0)
    curve.stairs(shown, np.r_[starts, starts[-1] + width[-1]], baseline=None, color=INK_2,
                 linewidth=1.3, zorder=3)
    zero = divergence == 0
    curve.scatter(starts[zero] + width[zero] / 2, np.full(zero.sum(), floor), s=10, color=INK_2,
                  zorder=4, linewidths=0)
    for state, rate in rates.items():
        color = STATE_COLOR.get(state, MUTED)
        curve.axhline(rate, color=color, linewidth=1.2, linestyle=(0, (5, 3)), zorder=2)
        label = {"background": f"Background\n{rate:.3%}",
                 "low_divergence": f"Low\n{rate:.4%}",
                 "high_divergence": f"High\n{rate:.2%}"}[state]
        curve.text(genome_mb * 1.012, rate, label, va="center", ha="left", fontsize=9, color=INK_2,
                   clip_on=False)
    curve.set_yscale("log")
    curve.set_ylim(floor / 1.6, 0.06)
    curve.set_yticks([floor, 1e-4, 1e-3, 1e-2], ["0", "0.01%", "0.1%", "1%"])
    curve.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    curve.set_xlim(0, genome_mb)
    curve.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    curve.set_xlabel("Genome position (Mb)")
    curve.set_ylabel("Divergence per 10 kb window")
    panel_letter(curve, "b")

    handles = [Patch(facecolor=SURFACE, edgecolor=MUTED, hatch="////", label="Planted (truth)"),
               Patch(facecolor=STATE_COLOR["low_divergence"], label=STATE_LABEL["low_divergence"]),
               Patch(facecolor=STATE_COLOR["high_divergence"], label=STATE_LABEL["high_divergence"]),
               Line2D([], [], color=INK_2, linewidth=1.3, label="Observed divergence")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.077, 0.985), ncol=4, fontsize=9.5,
               handlelength=1.6, columnspacing=1.6, labelcolor=INK)
    fig.savefig(out / "fig1_genome_overview.png", dpi=200)
    plt.close(fig)


def figure_2(truth, results, out):
    mosaics = read(results / "A_vs_B_popani" / "mosaic_regions.tsv")
    mismatches = read(results / "A_vs_B_popani" / "mismatches.tsv")["pos"].to_numpy()
    background = mosaics["background_divergence"][0]

    fig, axes = plt.subplots(1, len(truth), figsize=(13, 3.9), gridspec_kw={"wspace": 0.38})
    fig.subplots_adjust(left=0.06, right=0.985, top=0.78, bottom=0.18)
    for ax, row in zip(axes, truth.iter_rows(named=True)):
        flank = max(15_000, int(row["length"] * 0.6))
        lo, hi = row["start"] - flank, row["end"] + flank
        inside = mismatches[(mismatches >= lo) & (mismatches <= hi)]
        x = np.r_[lo, inside, hi] / 1e3
        y = np.r_[0, np.arange(1, inside.size + 1), inside.size]
        for block in mosaics.filter((pl.col("end") >= lo) & (pl.col("start") <= hi)).iter_rows(named=True):
            ax.axvspan(max(block["start"], lo) / 1e3, min(block["end"], hi) / 1e3,
                       color=STATE_COLOR[block["state"]], alpha=0.16, linewidth=0, zorder=0)
        for edge in (row["start"], row["end"]):
            ax.axvline(edge / 1e3, color=INK_2, linewidth=1, linestyle=(0, (2, 2)), zorder=1)
        # Background expectation, anchored at the observed count where the planted block starts.
        at_start = int((inside < row["start"]).sum())
        expected = at_start + (np.array([lo, hi]) - row["start"]) * background
        ax.plot([lo / 1e3, hi / 1e3], expected, color=MUTED, linewidth=1, linestyle=(0, (5, 3)), zorder=2)
        ax.step(x, y, where="post", color=INK, linewidth=1.6, zorder=3)
        ax.set_xlim(lo / 1e3, hi / 1e3)
        ax.set_ylim(0, max(y[-1], expected[-1]) * 1.08)
        ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
        ax.set_xlabel("Position (kb)")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        ax.set_title(row["block"], loc="left", fontweight="bold", fontsize=11, pad=8)
    axes[0].set_ylabel("Cumulative mismatches")

    handles = [Line2D([], [], color=INK, linewidth=1.6, label="Observed mismatches (cumulative)"),
               Line2D([], [], color=MUTED, linewidth=1, linestyle=(0, (5, 3)), label="Expected at background rate"),
               Line2D([], [], color=INK_2, linewidth=1, linestyle=(0, (2, 2)), label="Planted block edges"),
               Patch(facecolor=STATE_COLOR["low_divergence"], alpha=0.35, label="Detected low divergence"),
               Patch(facecolor=STATE_COLOR["high_divergence"], alpha=0.35, label="Detected high divergence")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.052, 0.99), ncol=5, fontsize=9.5,
               handlelength=1.8, columnspacing=1.6, labelcolor=INK)
    fig.savefig(out / "fig2_block_zooms.png", dpi=200)
    plt.close(fig)


def figure_3(truth, results, out):
    stretches = read(results / "A_vs_B_popani" / "identical_stretches.tsv")
    low = truth.filter(pl.col("kind") == "low_divergence")
    labels = {}
    for block in low.iter_rows(named=True):
        hit = stretches.filter((pl.col("start") <= block["end"]) & (pl.col("end") >= block["start"])) \
            .sort("length", descending=True)
        if hit.height:
            labels[hit["stretch_id"][0]] = block
    score = -np.log10(np.clip(stretches["p_value"].to_numpy(), 1e-300, None))
    length = stretches["length"].to_numpy() / 1e3
    planted = np.isin(stretches["stretch_id"].to_numpy(), list(labels))

    fig, ax = plt.subplots(figsize=(9, 5))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.9, bottom=0.12)
    threshold = 3  # p < 0.001
    ax.axhline(threshold, color=INK_2, linewidth=1, linestyle=(0, (5, 3)), zorder=1)
    ax.text(length.max() * 1.02, threshold + 0.45, "p = 0.001", ha="right", fontsize=9, color=INK_2)
    ax.scatter(length[~planted], score[~planted], s=18, color=MUTED, alpha=0.55, linewidths=0, zorder=2,
               label="Other stretches")
    ax.scatter(length[planted], score[planted], s=70, color=STATE_COLOR["low_divergence"],
               edgecolors=SURFACE, linewidths=1.5, zorder=3, label="Overlaps a planted identical block")
    ids = stretches["stretch_id"].to_numpy()
    arrow = {"arrowstyle": "-", "color": BASELINE, "linewidth": 1, "shrinkA": 2, "shrinkB": 5}
    for stretch_id, block in labels.items():
        index = int(np.flatnonzero(ids == stretch_id)[0])
        significant = score[index] > threshold
        ax.annotate(block["block"], (length[index], score[index]),
                    xytext=(-12, 10) if significant else (-16, 14), textcoords="offset points",
                    ha="right" if significant else "left", fontsize=10, fontweight="bold", color=INK,
                    arrowprops=arrow)
    ax.set_xlabel("Identical stretch length (kb)")
    ax.set_ylabel("−log10 p-value")
    ax.set_xlim(0, length.max() * 1.06)
    ax.set_ylim(-0.6, score.max() * 1.12)
    ax.grid(color=GRID, linewidth=0.8, zorder=0)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, fontsize=9.5, labelcolor=INK,
              borderaxespad=0, handletextpad=0.4, columnspacing=1.6)
    fig.savefig(out / "fig3_identical_stretches.png", dpi=200)
    plt.close(fig)


def figure_4(truth, results, profiles, out):
    comparisons = [("A_vs_B_popani", "A vs B", "popANI"), ("A_vs_M_popani", "A vs mix", "popANI"),
                   ("A_vs_M_conani", "A vs mix", "conANI")]
    # A real site inside the 30 kb import, where strain B differs from strain A.
    snps = read(Path(str(truth_path)).parent / "snps.tsv")
    block = truth.filter(pl.col("block") == "I1").row(0, named=True)
    site = snps.filter((pl.col("pos") >= block["start"]) & (pl.col("pos") <= block["end"])).row(0, named=True)
    counts = {}
    for sample in ("A", "B", "M"):
        row = pl.scan_parquet(profiles / f"sample_{sample}" / f"sample_{sample}_profile.parquet") \
            .filter(pl.col("pos") == site["pos"]).select("A", "C", "G", "T").collect().row(0, named=True)
        counts[sample] = (row[site["ref"]], row[site["alt"]])

    fig = plt.figure(figsize=(13.5, 4.9))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 2.25], wspace=0.42, left=0.115, right=0.895, top=0.76, bottom=0.14)
    alleles, tracks = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])

    # ── a: allele frequencies at one site ──
    names = {"A": "Sample A", "B": "Sample B", "M": "Mix\n70% B + 30% A"}
    ref_color, alt_color = INK_2, BASELINE
    for y, sample in enumerate(("M", "B", "A")):
        ref, alt = counts[sample]
        total = ref + alt
        parts = [(ref / total, ref_color, f"ref {site['ref']}"), (alt / total, alt_color, f"alt {site['alt']}")]
        left = 0
        for share, color, base in parts:
            if share == 0:
                continue
            gap = 0.006 if left > 0 else 0
            alleles.barh(y, share - gap, left=left + gap, height=0.56, color=color, linewidth=0)
            if share > 0.12:
                alleles.text(left + share / 2, y, f"{base}  {share:.0%}", ha="center", va="center",
                             fontsize=9.5, color="#ffffff" if color == ref_color else INK, fontweight="bold")
            left += share
    labels = [f"{names[s]}\n{sum(counts[s])}× depth" if s != "M" else f"{names[s]}, {sum(counts[s])}×"
              for s in ("M", "B", "A")]
    alleles.set_yticks(range(3), labels)
    alleles.tick_params(axis="y", length=0, labelcolor=INK)
    alleles.set_xlim(0, 1)
    alleles.set_xticks([0, 0.5, 1], ["0%", "50%", "100%"])
    alleles.spines["left"].set_visible(False)
    alleles.set_xlabel("Allele frequency")
    panel_letter(alleles, "a", pad=34)

    # ── b: blocks found by each comparison ──
    genome_mb = 2.160842
    rows = [("Planted", None, None)] + [(f"{label} · {method}", name, method) for name, label, method in comparisons]
    tracks.set_xlim(0, genome_mb)
    tracks.set_ylim(-0.6, len(rows) - 0.4)
    ys = list(range(len(rows) - 1, -1, -1))
    tracks.set_yticks(ys, [r[0] for r in rows])
    tracks.tick_params(axis="y", length=0, labelcolor=INK)
    tracks.spines["left"].set_visible(False)
    for y in ys:
        tracks.axhline(y, color=GRID, linewidth=1, zorder=1)
    for row in truth.iter_rows(named=True):
        rounded_bar(tracks, row["start"] / 1e6, row["length"] / 1e6, ys[0], 0.5,
                    color=STATE_COLOR[row["kind"]], hatch=True)
    for y, (_, name, _) in zip(ys[1:], rows[1:]):
        found = read(results / name / "mosaic_regions.tsv")
        mismatches = read(results / name / "windows.tsv")["mismatches"].sum()
        for block in found.iter_rows(named=True):
            rounded_bar(tracks, block["start"] / 1e6, (block["end"] - block["start"] + 1) / 1e6, y, 0.5,
                        color=STATE_COLOR[block["state"]])
        tracks.text(genome_mb * 1.015, y, f"{mismatches:,} mismatches", va="center", fontsize=9,
                    color=INK_2, clip_on=False)
    tracks.set_xlabel("Genome position (Mb)")
    panel_letter(tracks, "b", pad=34)

    handles = [Patch(facecolor=SURFACE, edgecolor=MUTED, hatch="////", label="Planted"),
               Patch(facecolor=STATE_COLOR["low_divergence"], label=STATE_LABEL["low_divergence"]),
               Patch(facecolor=STATE_COLOR["high_divergence"], label=STATE_LABEL["high_divergence"])]
    tracks.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, fontsize=9.2,
                  labelcolor=INK, borderaxespad=0, handlelength=1.4)
    fig.savefig(out / "fig4_mixed_population.png", dpi=200)
    plt.close(fig)


truth_path = None


def main():
    global truth_path
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth", default="data/truth/planted_blocks.tsv", type=Path)
    parser.add_argument("--results", default="results", type=Path)
    parser.add_argument("--profiles", default="profiles", type=Path)
    parser.add_argument("--outdir", default="figures", type=Path)
    args = parser.parse_args()
    truth_path = args.truth
    args.outdir.mkdir(parents=True, exist_ok=True)
    truth = read(args.truth)
    figure_1(truth, args.results, args.outdir)
    figure_2(truth, args.results, args.outdir)
    figure_3(truth, args.results, args.outdir)
    figure_4(truth, args.results, args.profiles, args.outdir)
    print(f"Wrote figures to {args.outdir}/")


if __name__ == "__main__":
    main()
