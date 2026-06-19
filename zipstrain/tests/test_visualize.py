import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pytest

from zipstrain import visualize as vz


def _comparison_frame() -> pl.LazyFrame:
    rows = [
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_b", "genome_pop_ani": 99.99, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_c", "genome_pop_ani": 97.0, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_d", "genome_pop_ani": 97.0, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_b", "sample_2": "sample_c", "genome_pop_ani": 97.0, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_b", "sample_2": "sample_d", "genome_pop_ani": 97.0, "total_positions": 20000},
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


def _sparse_comparison_frame_with_dense_core() -> pl.LazyFrame:
    rows = [
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_b", "genome_pop_ani": 99.9, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_c", "genome_pop_ani": 99.8, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_b", "sample_2": "sample_c", "genome_pop_ani": 99.7, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_d", "genome_pop_ani": 99.0, "total_positions": 20000},
        {"genome": "genome1", "sample_1": "sample_a", "sample_2": "sample_e", "genome_pop_ani": 99.0, "total_positions": 20000},
    ]
    return pl.DataFrame(rows).lazy()


def test_get_silhouette_plot_returns_dense_trace():
    fig = vz.get_silhouette_plot(_comparison_frame(), genome="genome1", min_comp_len=10000)

    assert len(fig.data) == 1
    assert len(fig.data[0].x) == 500
    assert len(fig.data[0].y) == 500
    assert min(fig.data[0].x) >= 99.0
    assert max(fig.data[0].x) <= 100.0


def test_get_silhouette_plot_warns_when_sklearn_is_unavailable(monkeypatch):
    monkeypatch.setattr(vz, "_sklearn_silhouette_score", None)

    with pytest.warns(RuntimeWarning, match="results might not be accurate"):
        fig = vz.get_silhouette_plot(
            _comparison_frame(),
            genome="genome1",
            min_comp_len=10000,
        )

    assert len(fig.data) == 1


def test_get_cluster_assignments_splits_two_clusters():
    cluster_df = vz.get_cluster_assignments(
        _comparison_frame(),
        clonal_cluster_threshold=99.8,
        strain_cluster_threshold=99.8,
    )
    cluster_map = {
        row["sample"]: (row["clonal_cluster"], row["strain_cluster"])
        for row in cluster_df.iter_rows(named=True)
    }

    assert cluster_df.columns == ["sample", "clonal_cluster", "strain_cluster"]
    assert cluster_df.get_column("clonal_cluster").n_unique() == 2
    assert cluster_df.get_column("strain_cluster").n_unique() == 2
    assert cluster_map["sample_a"] == cluster_map["sample_b"]
    assert cluster_map["sample_c"] == cluster_map["sample_d"]
    assert cluster_map["sample_a"] != cluster_map["sample_c"]


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


def test_get_clustermap_returns_clustergrid():
    grid = vz.get_clustermap(
        _comparison_frame(),
        genome="genome1",
        sample_to_population=_sample_to_population(),
    )
    try:
        assert hasattr(grid, "ax_heatmap")
    finally:
        plt.close(grid.fig)


def test_prepare_similarity_matrix_excludes_sparse_samples_before_matrix_build():
    bundle = vz._prepare_similarity_matrix(
        _sparse_comparison_frame_with_dense_core(),
        genome="genome1",
        min_comp_len=10000,
        max_null_samples=2,
    )

    assert bundle.samples == ["sample_a", "sample_b", "sample_c"]
    assert bundle.clustermap_data.columns == ["sample_1", "sample_a", "sample_b", "sample_c"]
    assert bundle.similarity_matrix.shape == (3, 3)


def test_prepare_similarity_matrix_requires_two_connected_samples_after_filtering():
    with pytest.raises(ValueError, match="At least two sufficiently connected samples"):
        vz._prepare_similarity_matrix(
            _sparse_comparison_frame_with_dense_core(),
            genome="genome1",
            min_comp_len=10000,
            max_null_samples=1,
        )


def test_prepare_similarity_matrix_builds_numpy_matrix_with_imputation_and_null_fraction():
    bundle = vz._prepare_similarity_matrix(
        _sparse_comparison_frame_with_dense_core(),
        genome="genome1",
        min_comp_len=10000,
        max_null_samples=10,
        impute_method=97.0,
    )

    assert bundle.samples == ["sample_a", "sample_b", "sample_c", "sample_d", "sample_e"]
    assert bundle.similarity_matrix.shape == (5, 5)
    assert bundle.similarity_matrix[0, 1] == pytest.approx(99.9)
    assert bundle.similarity_matrix[1, 0] == pytest.approx(99.9)
    assert bundle.similarity_matrix[1, 3] == pytest.approx(97.0)
    assert bundle.similarity_matrix[3, 1] == pytest.approx(97.0)

    null_fraction = dict(bundle.null_fraction.iter_rows())
    assert null_fraction["sample_a"] == pytest.approx(0.0)
    assert null_fraction["sample_d"] == pytest.approx(3 / 5)
    assert null_fraction["sample_e"] == pytest.approx(3 / 5)


def test_silhouette_score_precomputed_uses_sklearn_when_available(monkeypatch):
    calls = {}

    def _fake_sklearn(distance_matrix, labels, metric):
        calls["metric"] = metric
        calls["shape"] = distance_matrix.shape
        calls["labels"] = tuple(labels.tolist())
        return 0.321

    monkeypatch.setattr(vz, "_sklearn_silhouette_score", _fake_sklearn)

    distance_matrix = np.array(
        [
            [0.0, 0.1, 0.5],
            [0.1, 0.0, 0.4],
            [0.5, 0.4, 0.0],
        ]
    )
    labels = np.array([1, 1, 2])

    score = vz._silhouette_score_precomputed(distance_matrix, labels)

    assert score == pytest.approx(0.321)
    assert calls == {
        "metric": "precomputed",
        "shape": (3, 3),
        "labels": (1, 1, 2),
    }


def test_silhouette_score_precomputed_falls_back_to_manual(monkeypatch):
    monkeypatch.setattr(vz, "_sklearn_silhouette_score", None)

    distance_matrix = np.array(
        [
            [0.0, 0.1, 0.5, 0.6],
            [0.1, 0.0, 0.4, 0.5],
            [0.5, 0.4, 0.0, 0.1],
            [0.6, 0.5, 0.1, 0.0],
        ]
    )
    labels = np.array([1, 1, 2, 2])

    score = vz._silhouette_score_precomputed(distance_matrix, labels)
    manual_score = vz._silhouette_score_precomputed_manual(distance_matrix, labels)

    assert score == pytest.approx(manual_score)
