#!/usr/bin/env python3
"""Render comparable ZipStrain benchmark histories as standalone SVG charts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from packaging.version import Version


def plot_selection_history(records: list[dict], output: Path, group: str, ticks: list[int]) -> None:
    """Plot both selection integrations and methods separately for each genome layout."""
    releases = [record["version"] for record in records]
    first = records[0]["cases"]
    layouts = list(dict.fromkeys(case["layout"] for case in first if case["metric_set"] == "selection"))
    if not layouts:
        return
    colors = ["#206C87", "#D78336", "#163E62", "#993D24"]
    markers = ["o", "s", "^", "D"]
    for filename, metric, unit, scale in [
        ("selection-runtime.svg", "seconds", "Wall time (s)", 1),
        ("selection-memory.svg", "peak_rss_bytes", "Peak RSS (MiB)", 1024**2),
    ]:
        fig, axes = plt.subplots(2, len(layouts), figsize=(9.2, 6.6), sharex="col",
                                 sharey="row", layout="constrained", squeeze=False)
        title = f"Selection comparisons · {group}"
        if len(releases) == 1:
            title += f" · {releases[0]}"
        fig.suptitle(title, fontsize=15, x=0.03, ha="left")
        maxima = [0.0, 0.0]
        for column, layout in enumerate(layouts):
            indices = [index for index, case in enumerate(first)
                       if case["layout"] == layout and case["metric_set"] == "selection"]
            labels = []
            for position, index in enumerate(indices):
                cases = [record["cases"][index] for record in records]
                label = f"{cases[0]['method']} · {cases[0]['allele_integration']}"
                labels.append(label.replace(" · ", "\n"))
                collection = [median(run["collection"][metric] for run in case["runs"]) / scale for case in cases]
                pair = [median(pair[metric] for run in case["runs"] for pair in run["isolated_pairs"]) / scale
                        for case in cases]
                for row, values in enumerate((collection, pair)):
                    axis = axes[row, column]
                    maxima[row] = max(maxima[row], max(values))
                    if len(releases) == 1:
                        bars = axis.bar(position, values[0], width=0.62, color=colors[position])
                        axis.bar_label(bars, fmt="%.2f" if metric == "seconds" else "%.0f", padding=3)
                    else:
                        axis.plot(range(len(releases)), values, color=colors[position], marker=markers[position],
                                  linewidth=1.6, markersize=5, label=label)
            genome_count = first[indices[0]]["genomes"]
            genome_label = f"{genome_count} genome{'s' if genome_count != 1 else ''}"
            for row, description in enumerate(("Entire collection", "One isolated sample pair")):
                axis = axes[row, column]
                axis.set_title(f"{genome_label} · {description}", loc="left", fontsize=10)
                if column == 0:
                    axis.set_ylabel(unit)
                axis.yaxis.set_major_locator(MaxNLocator(5))
                axis.grid(axis="y", alpha=0.25)
                axis.set_axisbelow(True)
                axis.spines[["top", "right"]].set_visible(False)
            if len(releases) == 1:
                axes[1, column].set_xticks(range(len(labels)), labels, fontsize=8)
            else:
                axes[1, column].set_xticks(ticks, [releases[i] for i in ticks], rotation=30 if len(ticks) > 4 else 0)
                axes[1, column].set_xlabel("ZipStrain version")
        for row, maximum in enumerate(maxima):
            axes[row, 0].set_ylim(0, maximum * 1.16)
        if len(releases) > 1:
            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="outside lower center", frameon=False, ncol=4, fontsize=9)
        fig.savefig(output / filename, format="svg", metadata={"Date": None})
        plt.close(fig)


def draw_history(records: list[dict], output: Path) -> None:
    """Plot median total and isolated-pair measurements for each case and release."""
    records = sorted(records, key=lambda record: Version(record["version"]))
    if len({(record["suite_version"], record["series"], record["quick"]) for record in records}) != 1:
        raise ValueError("Plot a single benchmark suite and environment series at a time.")
    configurations = []
    for record in records:
        configurations.append([(
            case["layout"], case["metric_set"], case.get("method", "table"), case.get("allele_integration"),
            case["genomes"], case["bases_per_genome"], case["samples"], case["pairs"]
        ) for case in record["cases"]])
    if any(configuration != configurations[0] for configuration in configurations):
        raise ValueError("Benchmark cases or fixture sizes changed within a series.")

    output.mkdir(parents=True, exist_ok=True)
    releases = [record["version"] for record in records]
    x = list(range(len(releases)))
    tick_step = max(1, (len(x) + 7) // 8)
    ticks = sorted(set(x[::tick_step]) | {x[-1]})
    colors = ["#206C87", "#163E62", "#D78336", "#993D24"]
    markers = ["o", "s", "^", "D"]
    group = "GitHub release runner" if records[0]["series"] == "release" else "Local development baseline"
    for filename, metric, unit, scale in [
        ("runtime.svg", "seconds", "Wall time (s)", 1),
        ("memory.svg", "peak_rss_bytes", "Peak RSS (MiB)", 1024**2),
    ]:
        fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True, layout="constrained")
        title = f"{unit} · {group}"
        if len(records) == 1:
            title += f" · {releases[0]}"
        fig.suptitle(title, fontsize=15, x=0.03, ha="left")
        case_labels = []
        for index, case in enumerate(records[0]["cases"]):
            layout, metrics = case["layout"], case["metric_set"]
            if metrics == "selection":
                continue
            cases = [record["cases"][index] for record in records]
            total = [median(run["collection"][metric] for run in case["runs"]) / scale for case in cases]
            pair = [median(p[metric] for run in case["runs"] for p in run["isolated_pairs"]) / scale
                    for case in cases]
            label = f"{cases[0]['genomes']} genome{'s' if layout == 'multiple' else ''}"
            position = len(case_labels)
            case_labels.append(f"{label}\n{metrics} metrics")
            for axis, values in zip(axes, (total, pair), strict=True):
                color = colors[position % len(colors)]
                if len(records) == 1:
                    bars = axis.bar(position, values[0], width=0.56, color=color)
                    axis.bar_label(bars, fmt="%.2f" if metric == "seconds" else "%.0f", padding=3)
                else:
                    axis.plot(x, values, marker=markers[position % len(markers)], markersize=5,
                              linewidth=1.6, color=color, label=f"{label} · {metrics} metrics")
        for axis, title in zip(axes, ("Entire collection", "One isolated sample pair"), strict=True):
            axis.set_title(title, loc="left", fontsize=11)
            axis.set_ylabel(unit)
            axis.set_ylim(bottom=0)
            axis.margins(y=0.16)
            axis.yaxis.set_major_locator(MaxNLocator(5))
            axis.grid(axis="y", alpha=0.25)
            axis.set_axisbelow(True)
            axis.spines[["top", "right"]].set_visible(False)
        if len(records) == 1:
            axes[1].set_xticks(range(len(case_labels)), case_labels)
        else:
            axes[0].legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.18), frameon=False, fontsize=9)
            axes[1].set_xticks(ticks, [releases[index] for index in ticks], rotation=30 if len(ticks) > 4 else 0)
            axes[1].set_xlabel("ZipStrain version")
        fig.savefig(output / filename, format="svg", metadata={"Date": None})
        plt.close(fig)
    plot_selection_history(records, output, group, ticks)


def main() -> None:
    """Load JSON measurements and render the documentation figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("measurements", type=Path, nargs="+")
    args = parser.parse_args()
    records = [json.loads(path.read_text()) for path in args.measurements]
    draw_history(records, args.output_dir)


if __name__ == "__main__":
    main()
