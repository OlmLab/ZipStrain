import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import polars as pl
import pytest

from zipstrain import visualize as vz


def _comparison_frame() -> pl.LazyFrame:
    rows = [
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_b", "genome_pop_ani": 99.99, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_c", "genome_pop_ani": 99.20, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_d", "genome_pop_ani": 99.18, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_b", "sample_2": "sample_c", "genome_pop_ani": 99.19, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_b", "sample_2": "sample_d", "genome_pop_ani": 99.17, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_c", "sample_2": "sample_d", "genome_pop_ani": 99.98, "total_positions": 20000},
    ]
    return pl.DataFrame(rows).lazy()


def _multi_genome_comparison_frame() -> pl.LazyFrame:
    frame = _comparison_frame().collect()
    other = frame.with_columns(pl.lit("genome2").alias("genome"))
    return pl.concat([frame, other]).lazy()


def _sample_to_population() -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "sample": ["sample_a", "sample_b", "sample_c", "sample_d"],
            "population": ["pop_1", "pop_1", "pop_2", "pop_2"],
        }
    ).lazy()


def test_get_silhouette_plot_returns_dense_trace():
    fig = vz.get_silhouette_plot(_comparison_frame(), genome="genome1", min_comp_len=10000)

    assert len(fig.data) == 1
    assert len(fig.data[0].x) == 500
    assert len(fig.data[0].y) == 500
    assert min(fig.data[0].x) >= 99.0
    assert max(fig.data[0].x) <= 100.0


def test_get_cluster_assignments_splits_two_clusters():
    cluster_df = vz.get_cluster_assignments(_comparison_frame())

    assert cluster_df.columns == ["sample", "clonal_cluster", "strain_cluster"]
    assert cluster_df.get_column("clonal_cluster").n_unique() == 2
    assert cluster_df.get_column("strain_cluster").n_unique() == 2


def test_get_cluster_assignments_requires_one_genome_scope():
    with pytest.raises(ValueError, match="Multiple genomes are present"):
        vz.get_cluster_assignments(_multi_genome_comparison_frame())


def test_plot_dendo_returns_figure_with_fraction_bar_axis():
    fig = vz.plot_dendo(
        _comparison_frame(),
        genome="genome1",
        sample_to_population=_sample_to_population(),
        include_fraction_null=True,
    )
    try:
        assert len(fig.axes) == 2
    finally:
        plt.close(fig)


def test_get_clustermap_and_wrapper_return_clustergrid():
    grid = vz.get_clustermap(
        _comparison_frame(),
        genome="genome1",
        sample_to_population=_sample_to_population(),
    )
    wrapped = vz.plot_clustermap(
        _comparison_frame(),
        genome="genome1",
        sample_to_population=_sample_to_population(),
    )
    try:
        assert hasattr(grid, "ax_heatmap")
        assert hasattr(wrapped, "ax_heatmap")
    finally:
        plt.close(grid.fig)
        plt.close(wrapped.fig)
