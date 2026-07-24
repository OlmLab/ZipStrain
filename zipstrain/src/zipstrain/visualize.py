"""zipstrain.visualize
========================
This module provides statistical analysis and visualization functions for profiling and compare operations.
"""

from dataclasses import dataclass
import warnings
import polars as pl
import plotly.graph_objects as go
import seaborn as sns
import numpy as np
from itertools import chain, combinations
from collections import defaultdict
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colormaps
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.signal import find_peaks
from scipy.spatial.distance import squareform

try:
    from sklearn.metrics import silhouette_score as _sklearn_silhouette_score
except ImportError:
    _sklearn_silhouette_score = None


@dataclass(frozen=True)
class _SimilarityMatrixBundle:
    """Prepared similarity matrix and clustering artefacts for one genome scope."""

    clustermap_data: pl.DataFrame
    null_fraction: pl.DataFrame
    samples: list[str]
    similarity_matrix: np.ndarray
    distance_matrix: np.ndarray
    linkage_matrix: np.ndarray


@dataclass(frozen=True)
class SilhouetteCurveResult:
    """Computed silhouette sweep plus peak statistics for one genome scope."""

    thresholds: np.ndarray
    scores: np.ndarray
    curve: pl.DataFrame
    candidate_peaks: pl.DataFrame
    best_peak: pl.DataFrame
    min_threshold: float
    peak_prominence: float
    peak_distance: int
    used_sklearn: bool
    sample_count: int


def _silhouette_score_precomputed_manual(distance_matrix: np.ndarray, labels: np.ndarray) -> float:
    """Compute an average silhouette score from a precomputed distance matrix."""
    unique_labels = np.unique(labels)
    sample_count = distance_matrix.shape[0]
    if sample_count < 2 or len(unique_labels) < 2 or len(unique_labels) >= sample_count:
        return float("nan")

    scores: list[float] = []
    for idx in range(sample_count):
        same_cluster = labels == labels[idx]
        same_cluster[idx] = False
        if np.any(same_cluster):
            a_i = float(distance_matrix[idx, same_cluster].mean())
        else:
            scores.append(0.0)
            continue

        b_i = min(
            float(distance_matrix[idx, labels == other_label].mean())
            for other_label in unique_labels
            if other_label != labels[idx]
        )
        denom = max(a_i, b_i)
        scores.append(0.0 if denom == 0 else (b_i - a_i) / denom)
    return float(np.mean(scores))


def _silhouette_score_precomputed(distance_matrix: np.ndarray, labels: np.ndarray) -> float:
    """Compute an average silhouette score from a precomputed distance matrix."""
    unique_labels = np.unique(labels)
    sample_count = distance_matrix.shape[0]
    if sample_count < 2 or len(unique_labels) < 2 or len(unique_labels) >= sample_count:
        return float("nan")
    if _sklearn_silhouette_score is not None:
        return float(_sklearn_silhouette_score(distance_matrix, labels, metric="precomputed"))
    return _silhouette_score_precomputed_manual(distance_matrix, labels)


def _summarize_silhouette_curve(
    thresholds: np.ndarray,
    scores: np.ndarray,
    *,
    min_threshold: float = 99.8,
    peak_prominence: float = 0.001,
    peak_distance: int = 3,
) -> tuple[np.ndarray, np.ndarray, pl.DataFrame, pl.DataFrame]:
    """Clean a silhouette curve and extract candidate and best peaks."""
    threshold_array = np.asarray(thresholds, dtype=float)
    score_array = np.asarray(scores, dtype=float)
    mask = np.isfinite(threshold_array) & np.isfinite(score_array)
    threshold_array = threshold_array[mask]
    score_array = score_array[mask]
    order = np.argsort(threshold_array)
    threshold_array = threshold_array[order]
    score_array = score_array[order]
    if threshold_array.size == 0:
        raise ValueError("No finite silhouette scores are available for peak finding.")

    peaks, properties = find_peaks(
        score_array,
        prominence=peak_prominence,
        distance=peak_distance,
    )
    allowed_peak_mask = threshold_array[peaks] > min_threshold
    peaks = peaks[allowed_peak_mask]
    peak_prominences = properties.get("prominences", np.array([], dtype=float))[allowed_peak_mask]

    if len(peaks) > 0:
        best_idx = int(peaks[np.argmax(score_array[peaks])])
        best_source = "detected_peak"
    else:
        allowed = threshold_array > min_threshold
        if allowed.sum() == 0:
            raise ValueError(f"No thresholds above {min_threshold}")
        allowed_idx = np.where(allowed)[0]
        best_idx = int(allowed_idx[np.argmax(score_array[allowed_idx])])
        best_source = "fallback_global_max"

    candidate_peaks = pl.DataFrame(
        {
            "index": peaks.astype(int).tolist(),
            "threshold": threshold_array[peaks].astype(float).tolist(),
            "silhouette": score_array[peaks].astype(float).tolist(),
            "prominence": peak_prominences.astype(float).tolist(),
        }
    )
    if candidate_peaks.height > 0:
        candidate_peaks = candidate_peaks.sort("silhouette", descending=True)

    best_peak = pl.DataFrame(
        {
            "index": [best_idx],
            "threshold": [float(threshold_array[best_idx])],
            "silhouette": [float(score_array[best_idx])],
            "source": [best_source],
        }
    )
    return threshold_array, score_array, candidate_peaks, best_peak


def _prepare_similarity_matrix(
    comps_lf: pl.LazyFrame,
    *,
    genome: str | None = None,
    min_comp_len: int = 10000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    linkage_method: str = "average",
) -> _SimilarityMatrixBundle:
    """Build a square ANI similarity matrix and clustering inputs."""
    schema_names = set(comps_lf.collect_schema().names())
    if genome is not None and "genome" not in schema_names:
        raise ValueError("A genome was requested, but the comparison table has no 'genome' column.")
    if genome is None and "genome" in schema_names:
        genome_count = (
            comps_lf.select(pl.col("genome"))
            .unique()
            .collect(engine="streaming")
            .height
        )
        if genome_count > 1:
            raise ValueError("Multiple genomes are present. Pass `genome=...` or pre-filter to one genome.")

    if not isinstance(impute_method, (int, float)):
        raise NotImplementedError("Only numeric ANI imputation is currently supported.")

    filters = [pl.col("total_positions") > min_comp_len]
    if genome is not None:
        filters.append(pl.col("genome") == genome)
    filtered = (
        comps_lf.filter(pl.all_horizontal(*filters))
        .select("sample_1", "sample_2", "genome_ani")
        .group_by(["sample_1", "sample_2"])
        .agg(pl.col("genome_ani").mean().alias("genome_ani"))
    )
    filtered_pairs = filtered.collect(engine="streaming")
    sample_names = (
        pl.concat(
            [
                filtered_pairs.select(pl.col("sample_1").alias("sample")),
                filtered_pairs.select(pl.col("sample_2").alias("sample")),
            ],
            how="vertical",
        )
        .get_column("sample")
        .unique()
        .sort()
        .to_list()
    )
    if len(sample_names) < 2:
        raise ValueError("At least two samples are required to build a similarity matrix.")

    comparable_counts = (
        pl.concat(
            [
                filtered_pairs.select(pl.col("sample_1").alias("sample")),
                filtered_pairs.select(pl.col("sample_2").alias("sample")),
            ],
            how="vertical",
        )
        .group_by("sample")
        .len()
        .with_columns((pl.col("len") + 1).alias("comparable_count"))
        .select("sample", "comparable_count")
    )
    total_sample_count = len(sample_names)
    exclude_samples = (
        comparable_counts
        .filter((pl.lit(total_sample_count) - pl.col("comparable_count")) > max_null_samples)
        .get_column("sample")
        .to_list()
    )
    if exclude_samples:
        exclude_set = set(exclude_samples)
        sample_names = [sample for sample in sample_names if sample not in exclude_set]
        if len(sample_names) < 2:
            raise ValueError(
                "At least two sufficiently connected samples are required to build a similarity matrix. "
                "Relax `max_null_samples` or provide a denser comparison table."
            )
        filtered_pairs = filtered_pairs.filter(
            (~pl.col("sample_1").is_in(exclude_samples))
            & (~pl.col("sample_2").is_in(exclude_samples))
        )

    samples = sample_names
    sample_index = pl.DataFrame(
        {
            "sample": samples,
            "sample_idx": np.arange(len(samples), dtype=np.int32),
        }
    )
    sample_index_lf = sample_index.lazy()
    indexed_pairs = (
        filtered_pairs.lazy()
        .join(sample_index_lf, left_on="sample_1", right_on="sample", how="inner")
        .rename({"sample_idx": "sample_idx_1"})
        .join(sample_index_lf, left_on="sample_2", right_on="sample", how="inner")
        .rename({"sample_idx": "sample_idx_2"})
        .select("sample_idx_1", "sample_idx_2", "genome_ani")
        .collect(engine="streaming")
    )

    similarity_matrix_raw = np.full((len(samples), len(samples)), np.nan, dtype=np.float64)
    np.fill_diagonal(similarity_matrix_raw, 100.0)
    if indexed_pairs.height > 0:
        sample_idx_1 = indexed_pairs.get_column("sample_idx_1").to_numpy()
        sample_idx_2 = indexed_pairs.get_column("sample_idx_2").to_numpy()
        genome_ani = indexed_pairs.get_column("genome_ani").to_numpy()
        similarity_matrix_raw[sample_idx_1, sample_idx_2] = genome_ani
        similarity_matrix_raw[sample_idx_2, sample_idx_1] = genome_ani

    comparable_column_count = max(len(samples), 1)
    null_fraction = pl.DataFrame(
        {
            "sample_1": samples,
            "null_fraction": np.isnan(similarity_matrix_raw).sum(axis=1) / comparable_column_count,
        }
    )

    similarity_matrix = np.where(
        np.isnan(similarity_matrix_raw),
        float(impute_method),
        similarity_matrix_raw,
    )
    clustermap_data = pl.DataFrame(
        {
            "sample_1": samples,
            **{sample: similarity_matrix[:, idx] for idx, sample in enumerate(samples)},
        }
    )
    distance_matrix = 1 - (similarity_matrix / 100.0)
    np.fill_diagonal(distance_matrix, 0.0)
    linkage_matrix = linkage(squareform(distance_matrix, checks=False), method=linkage_method)
    return _SimilarityMatrixBundle(
        clustermap_data=clustermap_data,
        null_fraction=null_fraction,
        samples=samples,
        similarity_matrix=similarity_matrix,
        distance_matrix=distance_matrix,
        linkage_matrix=linkage_matrix,
    )


def _resolve_population_mapping(
    sample_to_population: pl.LazyFrame,
    samples: list[str],
) -> pl.DataFrame:
    """Return a sample/population frame aligned to a provided sample order."""
    return (
        pl.DataFrame({"sample_1": samples})
        .join(
            sample_to_population.collect(engine="streaming"),
            left_on="sample_1",
            right_on="sample",
            how="left",
        )
        .with_columns(pl.col("population").fill_null("Not Assigned"))
        .select("sample_1", "population")
    )


def _resolve_population_colors(
    sample_to_population: pl.DataFrame,
    *,
    color_map: dict | None = None,
) -> tuple[dict, list]:
    """Build a stable population -> color mapping and row colors."""
    if color_map is None:
        unique_pops = sample_to_population.get_column("population").unique().sort().to_list()
        palette = colormaps["tab20b"]
        divisor = max(len(unique_pops), 1)
        color_map = {population: palette(idx / divisor) for idx, population in enumerate(unique_pops)}
    row_colors = [color_map.get(population, "black") for population in sample_to_population["population"]]
    return color_map, row_colors


def get_cdf(data, num_bins=10000):
    """Calculate the cumulative distribution function (CDF) of the given data."""
    if data[0] == -1:
        return [-1], [-1]
    counts, bin_edges = np.histogram(data, bins=np.linspace(0, 50000, num_bins))
    counts = counts[::-1]
    bin_edges = bin_edges[::-1]
    cummulative_counts = np.cumsum(counts)
    cdf= cummulative_counts / cummulative_counts[-1]
    return bin_edges, cdf

def calculate_strainsharing(
                            comps_lf:pl.LazyFrame,
                            breadth_lf:pl.LazyFrame,
                            sample_to_population:pl.LazyFrame,
                            min_breadth:float=0.5,
                            strain_similarity_threshold:float=99.9,
                            min_total_positions:int=10000
                            )->dict[str, list[float]]:


    """
    Calculate strain sharing between populations based on genome ANI between genomes in their profiles.
    Strain sharing between two samples is defined as the ratio of genomes passing a strain similarity threshold over the total number of genomes in each sample.
    So, for two samples A and B, the strain sharing is defined as (Note the assymetric nature of the calculation):
    strain_sharing(A, B) = (number of genomes in A and B passing the strain similarity threshold) / (number of genomes in A)
    strain_sharing(B, A) = (number of genomes in A and B passing the strain similarity threshold) / (number of genomes in B)
    
    Args:
        comps_lf (pl.LazyFrame): LazyFrame containing the gene profiles of the samples.
        breadth_lf (pl.LazyFrame): LazyFrame containing the genome breadth information.
        sample_to_population (pl.LazyFrame): LazyFrame containing the sample to population mapping.
        min_breadth (float, optional): Minimum genome breadth to consider a genome for strain sharing. Defaults to 0.5.
        strain_similarity_threshold (float, optional): Threshold for strain similarity. Defaults to 0.99.
        min_total_positions (int, optional): Minimum total positions to consider a genome for strain sharing. Defaults to 10000.
    Returns:
        pl.LazyFrame: LazyFrame containing the strain sharing information between populations. It will be in the following form [Sample A, Sample B, Strain Sharing, Relationship]
    """
    comps_lf=comps_lf.filter(
        (pl.col("total_positions")>min_total_positions)
    ).collect(engine="streaming").lazy()
    breadth_lf=breadth_lf.fill_null(0.0)
    breadth_lf_long=(
        breadth_lf.unpivot(
            index=["genome"],
            variable_name="sample",
            value_name="breadth"
        )
    )
    breadth_lf=breadth_lf_long.group_by("sample").agg(num_genomes=(pl.col("breadth")>=min_breadth).sum())
    comps_lf=comps_lf.join(breadth_lf,
        left_on='sample_1',
        right_on='sample',
        how='left',
    ).rename(
        {"num_genomes":"num_genomes_1"}
    ).join(
        breadth_lf,
        left_on='sample_2',
        right_on='sample',
        how='left',
    ).rename(
        {"num_genomes":"num_genomes_2"}
    )
    comps_lf = comps_lf.join(
        sample_to_population,
        left_on='sample_1',
        right_on='sample',
        how='left',
    ).rename(
        {"population":"population_1"}
    ).join(
        sample_to_population,
        left_on='sample_2',
        right_on='sample',
        how='left',
    ).rename(
        {"population":"population_2"}
    )
    comps_lf=comps_lf.join(
        breadth_lf_long,
        left_on=["genome","sample_1"],
        right_on=['genome','sample'],
        how='left',
    ).rename(
        {"breadth":"breadth_1"}
    ).join(
        breadth_lf_long,
        left_on=["genome","sample_2"],
        right_on=['genome','sample'],
        how='left',
    ).rename(
        {"breadth":"breadth_2"}
    )
    comps_lf=comps_lf.filter(
        (pl.col("breadth_1") >= min_breadth) &
        (pl.col("breadth_2") >= min_breadth) &
        (pl.col("genome_ani") >= strain_similarity_threshold)
    )

    comps_lf=comps_lf.group_by(
        ["sample_1", "sample_2"]
    ).agg(
        pl.col("genome").count().alias("shared_strain_count"),
        pl.col("num_genomes_1").first().alias("num_genomes_1"),
        pl.col("num_genomes_2").first().alias("num_genomes_2"),
        pl.col("population_1").first().alias("population_1"),
        pl.col("population_2").first().alias("population_2"),
    ).collect(engine="streaming")
    strainsharingrates=defaultdict(list)
    for row in comps_lf.iter_rows(named=True):
        strainsharingrates[row["population_1"]+"_"+ row["population_2"]].append(row["shared_strain_count"] / row["num_genomes_1"])
        strainsharingrates[row["population_2"]+"_"+ row["population_1"]].append(row["shared_strain_count"] / row["num_genomes_2"])
    return strainsharingrates

def calculate_ibs(
    sample_to_population:pl.LazyFrame, 
    comps_lf:pl.LazyFrame,
    max_perc_id_genes:float=15,
    min_total_positions:int=10000,
)->pl.DataFrame:
    """
    Calculate the Identity By State (IBS) between two populations for a given genome.
    The IBS is defined as the percentage of genes that are identical between two populations for a given genome.
    Args:
        sample_to_population (pl.LazyFrame): LazyFrame containing the sample to population mapping.
        comps_lf (pl.LazyFrame): LazyFrame containing the gene profiles of the samples.
        max_perc_id_genes (float, optional): Maximum percentage of identical genes to consider. Defaults to 0.15.
    Returns:
        pl.LazyFrame: LazyFrame containing the IBS information for the given genome and populations.
    """
    comps_lf_filtered = comps_lf.filter(
        (pl.col('perc_id_genes') <= max_perc_id_genes) &
        (pl.col('total_positions')>min_total_positions)
    )
    comps_lf_filtered=comps_lf_filtered.join(
        sample_to_population,
        left_on='sample_1',
        right_on='sample',
        how='inner',
    ).rename(
        {"population":"population_1"}
    ).join(
        sample_to_population,
        left_on='sample_2',
        right_on='sample',
        how='inner',
        suffix='_2'
    ).rename(
        {"population":"population_2"}
    )
    comps_lf_filtered = comps_lf_filtered.with_columns(
    pl.when(pl.col("population_1") == pl.col("population_2"))
    .then(
        pl.lit("within_population_")
        + pl.col("population_1")
        + pl.lit("|")
        + pl.col("population_2")
    )
    .otherwise(
        pl.concat_str(
            [
                pl.lit("between_population_"),
                pl.concat_str(
                    [
                        pl.min_horizontal("population_1", "population_2"),
                        pl.lit("|"),
                        pl.max_horizontal("population_1", "population_2"),
                    ]
                ),
            ]
        )
    )
    .alias("comparison_type")
    ).fill_null(-1)

    return comps_lf_filtered.group_by(["genome","comparison_type"]).agg(
        pl.col("max_consecutive_length"),
    ).collect(engine="streaming").pivot(
        index="genome",
        columns="comparison_type",
        values="max_consecutive_length",
    ).with_columns(
        pl.col("*").exclude("genome").fill_null([-1])
    )

def plot_ibs_heatmap(
    df:pl.DataFrame,
    vert_thresh:float=0.001,
    populations:list[str]|None=None,
    num_bins:int=10000,
    min_member:int=50,
    title:str="IBS Heatmap",
    xaxis_title:str="Population Pair",
    yaxis_title:str="Genome",
    
):
    """
    Plot the Identity By State (IBS) heatmap for a given genome and two populations.
    Args:
        df (pl.DataFrame): DataFrame containing the IBS information.
        title (str, optional): Title of the plot. Defaults to "IBS Heatmap".
        xaxis_title (str, optional): Title of the x-axis. Defaults to "Population Pair".
        yaxis_title (str, optional): Title of the y-axis. Defaults to "Genome".
    Returns:
        go.Figure: Plotly figure containing the IBS heatmap.
    """
    df = df.with_columns(
    [
        pl.when(pl.col(c).list.len() < min_member)
        .then(pl.lit([-1]))
        .otherwise(pl.col(c))
        .alias(c)
        for c in df.columns if c != "genome"
    ]
)
    if populations is None:
        populations=set(chain.from_iterable(i.replace("within_population_","").replace("between_population_","").split("|") for i in df.columns if i!="genome"))
        populations=sorted(populations)
    heatmap_data = df.rows_by_key("genome", unique=True,include_key=False,named=True)
    fig_data={}
    for genome, genome_data in heatmap_data.items():
        fig_data[genome]={}
        for pop1,pop2 in combinations(populations,2):
            key_between=f"between_population_{min(pop1,pop2)}|{max(pop1,pop2)}"
            key_within_1=f"within_population_{pop1}|{pop1}"
            key_within_2=f"within_population_{pop2}|{pop2}"
            if genome_data.get(key_between, [-1])==[-1] or genome_data.get(key_within_1, [-1])==[-1] or genome_data.get(key_within_2, [-1])==[-1]:
                fig_data[genome][f"{min(pop1,pop2)}-{max(pop1,pop2)}"]=-1
                continue
            between=get_cdf(genome_data[key_between], num_bins=num_bins)
            within=get_cdf(genome_data[key_within_1]+genome_data[key_within_2], num_bins=num_bins)

            between_intersect=between[0][np.where(between[1]>=vert_thresh)[0][0]]
            within_intersect=within[0][np.where(within[1]>=vert_thresh)[0][0]]
            distance=within_intersect-between_intersect
            fig_data[genome][f"{min(pop1,pop2)}-{max(pop1,pop2)}"]=distance
    ###Filter the dataframe to only have useful information
    heatmap_df = pd.DataFrame(fig_data).T
    heatmap_df=heatmap_df.mask(heatmap_df < 0, 0)
    heatmap_df=heatmap_df[heatmap_df.sum(axis=1)>0]
    heatmap_df=heatmap_df[[col for col in heatmap_df.columns if heatmap_df[col].sum()>0]]
    heatmap_df_sorted = heatmap_df.assign(row_sum=heatmap_df.sum(axis=1)).sort_values("row_sum", ascending=True).drop(columns="row_sum")
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_df_sorted.values,
        x=heatmap_df_sorted.columns,
        y=heatmap_df_sorted.index
    ))
    return fig

def plot_strainsharing(
    strainsharingrates:dict[str, list[float]],
    sample_frac:float=1,
    title:str="Strain Sharing Rates",
    xaxis_title:str="Population Pair",
    yaxis_title:str="Strain Sharing Rate",
):
    """
    Plot the strain sharing rates between populations.
    Args:
        strainsharingrates (dict[str, list[float]]): Dictionary containing the strain sharing rates between populations.
        title (str, optional): Title of the plot. Defaults to "Strain Sharing".
        xaxis_title (str, optional): Title of the x-axis. Defaults to "Population Pair".
        yaxis_title (str, optional): Title of the y-axis. Defaults to "Strain Sharing Rate".
    Returns:
        go.Figure: Plotly figure containing the strain sharing plot.
    """
    for key in strainsharingrates.keys():
        strainsharingrates[key] = np.random.choice(strainsharingrates[key], size=int(len(strainsharingrates[key]) * sample_frac), replace=False)
    fig = go.Figure()
    for pair, rates in strainsharingrates.items():
        fig.add_trace(go.Box(
            y=rates,
            name=pair,
            boxpoints='all',
            jitter=0.3,
            pointpos=0
        ))
    fig.update_layout(
        title={"text": title, "x": 0.5},
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title
    )
    return fig
def plot_ibs(df:pl.DataFrame,
            genome:str,
            population_1:str,
            population_2:str,
            vert_thresh_hor_distance:float=0.001,
            num_bins:int=10000,
            title:str="IBS for <GENOME>: <POPULATION_1> vs <POPULATION_2>",
            xaxis_title:str="Max Consecutive Length",
            yaxis_title:str="CDF"
            ):
    """
    Plot the Identity By State (IBS) for a given genome and two populations.
    Args:
        df (pl.DataFrame): DataFrame containing the IBS information.
        genome (str): The genome to plot the IBS for.
        population_1 (str): The first population to plot the IBS for.
        population_2 (str): The second population to plot the IBS for.
        title (str, optional): Title of the plot. Defaults to "IBS for <GENOME>".
        xaxis_title (str, optional): Title of the x-axis. Defaults to "Membership".
        yaxis_title (str, optional): Title of the y-axis. Defaults to "Max Consecutive Length".
    Returns:
        go.Figure: Plotly figure containing the IBS plot.
    """
    df_filtered = df.filter(pl.col("genome") == genome)
    if df_filtered.is_empty():
        raise ValueError(f"Genome {genome} not found in the dataframe.")
    plot_data = {}
    key_within_1=f"within_population_{population_1}|{population_1}"
    key_within_2=f"within_population_{population_2}|{population_2}"
    key_between=f"between_population_{min(population_1,population_2)}|{max(population_1,population_2)}"
    if df_filtered.get_column(key_within_1).list.len()[0]==0 or df_filtered.get_column(key_within_2).list.len()[0]==0 or df_filtered.get_column(key_between).list.len()[0]==0:
        raise ValueError(f"Not enough data for populations {population_1} and {population_2} in genome {genome}.")
    plot_data["within_population"]=df_filtered.get_column(key_within_1)[0].to_list()+df_filtered.get_column(key_within_2)[0].to_list()
    plot_data["between_population"]=df_filtered.get_column(key_between)[0].to_list()
    fig = go.Figure()
    between_pop_cdf=get_cdf(plot_data["between_population"], num_bins=num_bins)
    fig.add_trace(go.Scatter(
        x=between_pop_cdf[0][1:].copy(),
        y=between_pop_cdf[1][1:].copy(),
        mode='lines',
        name='between_population',
        line=dict(color='blue')
    ))
    within_pop_cdf=get_cdf(plot_data["within_population"], num_bins=num_bins)
    fig.add_trace(go.Scatter(
        x=within_pop_cdf[0][1:].copy(),
        y=within_pop_cdf[1][1:].copy(),
        mode='lines',
        name='within_population',
        line=dict(color='green')
    ))

    bin_edges=within_pop_cdf[0]
    cdf=within_pop_cdf[1]
    within_intersect=bin_edges[np.where(cdf>=vert_thresh_hor_distance)[0][0]]
    bin_edges=between_pop_cdf[0]
    cdf=between_pop_cdf[1]
    between_intersect=bin_edges[np.where(cdf>=vert_thresh_hor_distance)[0][0]]  
    distance=within_intersect-between_intersect

    fig.update_layout(
        title={"text": title.replace("<GENOME>", genome).replace("<POPULATION_1>", population_1).replace("<POPULATION_2>", population_2), "x": 0.5},
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        
    )
    ###Add a horizontal line from (between_intersect, vert_thresh_hor_distance) to (within_intersect, vert_thresh_hor_distance)
    fig.add_trace(go.Scatter(
        x=[between_intersect, within_intersect],
        y=[vert_thresh_hor_distance, vert_thresh_hor_distance],
        mode='lines+markers',
        line=dict(color='black'),
        showlegend=False
    ))
    ###Add a text annotation at the middle of the horizontal line with the distance
    fig.add_trace(go.Scatter(
        x=[(between_intersect+within_intersect)/2],
        y=[vert_thresh_hor_distance],
        mode="text",
        text=int(distance),
        textposition="top center",
        showlegend=False
    ))
    ##make both axes logarithmic
    fig.update_xaxes(type='log')
    fig.update_yaxes(type='log')

    return fig
def calculate_identical_frac_vs_popani(
    genome:str,
    population_1:str,
    population_2:str,
    sample_to_population:pl.LazyFrame,
    comps_lf:pl.LazyFrame,
    min_shared_genes_count:int=100,
    min_total_positions:int=10000
    ):
    """
    Calculate the fraction of identical genes vs genome ANI for a given genome and two samples in any possible combination of populations.
    Args:
        genome (str): The genome to calculate the fraction of identical genes vs genome ANI for.
        population_1 (str): The first population to compare.
        population_2 (str): The second population to compare.
        sample_to_population (pl.LazyFrame): LazyFrame containing the sample to population mapping.
        comps_lf (pl.LazyFrame): LazyFrame containing the gene profiles of the samples
    Returns:
        pl.LazyFrame: LazyFrame containing the fraction of identical genes vs genome ANI information for
    """
    comps_lf_filtered=comps_lf.filter(
        (pl.col('genome') == genome) &
        (pl.col("shared_genes_count")>min_shared_genes_count) &
        (pl.col("total_positions")>min_total_positions)
    ).collect(engine="streaming").lazy()

    comps_lf_filtered=comps_lf_filtered.join(
        sample_to_population,
        left_on='sample_1',
        right_on='sample',
        how='left',
    ).rename(
        {"population":"population_1"}
    ).join(
        sample_to_population,
        left_on='sample_2',
        right_on='sample',
        how='left',
        suffix='_2'
    ).rename(
        {"population":"population_2"}
    )
    comps_lf_filtered = comps_lf_filtered.filter(
        (pl.col("population_1").is_in({population_1, population_2})) &
        (pl.col("population_2").is_in({population_1, population_2}))
    ).collect(engine="streaming").lazy()
    groups={
        "same_1":f"{population_1}-{population_1}",
        "same_2":f"{population_2}-{population_2}",
        "diff":f"{population_1}-{population_2}",
    }
    comps_lf_filtered=comps_lf_filtered.with_columns(
        pl.when((pl.col("population_1")==population_1) & (pl.col("population_2")==population_1))
        .then(pl.lit(groups["same_1"]))
        .when((pl.col("population_1")==population_2) & (pl.col("population_2")==population_2))
        .then(pl.lit(groups["same_2"]))
        .otherwise(pl.lit(groups["diff"]))
        .alias("relationship")
    )
    return comps_lf_filtered.group_by("relationship").agg(
        pl.col("perc_id_genes"),
        pl.col("genome_ani")
    ).collect(engine="streaming")

def plot_identical_frac_vs_popani(df:pl.DataFrame,
                                  genome:str,
                                  title:str="Fraction of Identical Genes vs Genome ANI for <GENOME>",
                                  xaxis_title:str="Genome-Wide ANI",
                                  yaxis_title:str="Fraction of Identical Genes",
                                  ):
    """
    Plot the fraction of identical genes vs genome ANI for a given genome and two samples in any possible combination of populations.
    Args:
        df (pl.DataFrame): DataFrame containing the fraction of identical genes vs genome ANI information.
        title (str, optional): Title of the plot. Defaults to "Fraction of Identical Genes vs Genome ANI".
        xaxis_title (str, optional): Title of the x-axis. Defaults to "Genome ANI".
        yaxis_title (str, optional): Title of the y-axis. Defaults to "Fraction of Identical Genes".
    Returns:
        go.Figure: Plotly figure containing the fraction of identical genes vs genome ANI plot.
    """
    fig = go.Figure()
    for group, perc_id_genes, genome_ani in zip(df["relationship"], df["perc_id_genes"], df["genome_ani"]):
        fig.add_trace(go.Scatter(
            x=genome_ani,
            y=perc_id_genes,
            mode='markers',
            name=group
        ))
    fig.update_layout(
        title=title.replace("<GENOME>", genome),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title
    )
    return fig

def compute_silhouette_curve(
    comps_lf: pl.LazyFrame,
    genome: str,
    min_comp_len: int = 100000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    linkage_method: str = "average",
    min_threshold: float = 99.8,
    peak_prominence: float = 0.001,
    peak_distance: int = 3,
) -> SilhouetteCurveResult:
    """Compute a silhouette sweep and peak summary for one genome."""
    bundle = _prepare_similarity_matrix(
        comps_lf,
        genome=genome,
        min_comp_len=min_comp_len,
        impute_method=impute_method,
        max_null_samples=max_null_samples,
        linkage_method=linkage_method,
    )
    distances = np.linspace(0.01, 0.0, 500)
    scores: list[float] = []
    use_sklearn = _sklearn_silhouette_score is not None
    if not use_sklearn:
        warnings.warn(
            "scikit-learn is not installed; falling back to the manual silhouette implementation, "
            "and results might not be accurate.",
            RuntimeWarning,
            stacklevel=2,
        )
    for distance_threshold in distances:
        labels = fcluster(bundle.linkage_matrix, t=distance_threshold, criterion="distance")
        scores.append(_silhouette_score_precomputed(bundle.distance_matrix, labels))

    thresholds = 100 * (1 - distances)
    clean_thresholds, clean_scores, candidate_peaks, best_peak = _summarize_silhouette_curve(
        thresholds,
        np.asarray(scores, dtype=float),
        min_threshold=min_threshold,
        peak_prominence=peak_prominence,
        peak_distance=peak_distance,
    )
    curve = pl.DataFrame(
        {
            "threshold": clean_thresholds.astype(float).tolist(),
            "silhouette": clean_scores.astype(float).tolist(),
        }
    )
    return SilhouetteCurveResult(
        thresholds=clean_thresholds,
        scores=clean_scores,
        curve=curve,
        candidate_peaks=candidate_peaks,
        best_peak=best_peak,
        min_threshold=min_threshold,
        peak_prominence=peak_prominence,
        peak_distance=peak_distance,
        used_sklearn=use_sklearn,
        sample_count=len(bundle.samples),
    )


def plot_silhouette_curve(
    result: SilhouetteCurveResult,
    *,
    title: str = "Simple peak finding on silhouette curve",
    xaxis_title: str = "Clustering threshold",
    yaxis_title: str = "Silhouette score",
) -> go.Figure:
    """Plot a computed silhouette sweep and its peak summary."""
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=result.thresholds,
            y=result.scores,
            mode="lines+markers",
            name="Silhouette score",
            marker=dict(size=5),
            line=dict(width=2),
        )
    )

    if result.candidate_peaks.height > 0:
        peak_indices = result.candidate_peaks.get_column("index").to_numpy()
        fig.add_trace(
            go.Scatter(
                x=result.candidate_peaks.get_column("threshold").to_list(),
                y=result.candidate_peaks.get_column("silhouette").to_list(),
                mode="markers",
                name=f"Detected peaks above {result.min_threshold}",
                marker=dict(size=10, symbol="circle"),
                customdata=np.column_stack([peak_indices]),
                hovertemplate=(
                    "threshold=%{x:.6f}<br>"
                    "silhouette=%{y:.6f}<br>"
                    "index=%{customdata[0]}<extra></extra>"
                ),
            )
        )

    best_peak_row = result.best_peak.row(0, named=True)
    best_threshold = float(best_peak_row["threshold"])
    best_silhouette = float(best_peak_row["silhouette"])
    best_index = int(best_peak_row["index"])
    fig.add_trace(
        go.Scatter(
            x=[best_threshold],
            y=[best_silhouette],
            mode="markers+text",
            name="Best peak",
            marker=dict(size=18, symbol="star"),
            text=[f"best = {best_threshold:.5g}"],
            textposition="top center",
            customdata=[[best_index]],
            hovertemplate=(
                "BEST<br>"
                "threshold=%{x:.6f}<br>"
                "silhouette=%{y:.6f}<br>"
                "index=%{customdata[0]}<extra></extra>"
            ),
        )
    )
    fig.add_vline(
        x=result.min_threshold,
        line_dash="dot",
        opacity=0.7,
        annotation_text=f"min threshold = {result.min_threshold}",
        annotation_position="top left",
    )
    fig.add_vline(
        x=best_threshold,
        line_dash="dash",
        opacity=0.8,
        annotation_text=f"best = {best_threshold:.5g}",
        annotation_position="top right",
    )
    fig.update_layout(
        title={"text": title, "x": 0.5},
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        template="plotly_white",
        width=1000,
        height=550,
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
    )
    fig.update_xaxes(range=[99, 100], autorange="reversed")
    return fig


def get_silhouette_plot(
    comps_lf: pl.LazyFrame,
    genome: str,
    min_comp_len: int = 100000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    linkage_method: str = "average",
    min_threshold: float = 99.8,
    peak_prominence: float = 0.001,
    peak_distance: int = 3,
):
    """Plot silhouette score as a function of ANI threshold for one genome."""
    return plot_silhouette_curve(
        compute_silhouette_curve(
            comps_lf,
            genome=genome,
            min_comp_len=min_comp_len,
            impute_method=impute_method,
            max_null_samples=max_null_samples,
            linkage_method=linkage_method,
            min_threshold=min_threshold,
            peak_prominence=peak_prominence,
            peak_distance=peak_distance,
        )
    )


def get_cluster_assignments(
    comps_lf: pl.LazyFrame,
    min_comp_len: int = 10000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    clonal_cluster_threshold: float = 99.93,
    strain_cluster_threshold: float = 99.8,
    linkage_method: str = "average",
):
    """Get clonal and strain-level cluster assignments from a genome-scoped comparison table."""
    bundle = _prepare_similarity_matrix(
        comps_lf,
        genome=None,
        min_comp_len=min_comp_len,
        impute_method=impute_method,
        max_null_samples=max_null_samples,
        linkage_method=linkage_method,
    )
    clonal_clusters = fcluster(
        bundle.linkage_matrix,
        t=1 - clonal_cluster_threshold / 100,
        criterion="distance",
    )
    strain_clusters = fcluster(
        bundle.linkage_matrix,
        t=1 - strain_cluster_threshold / 100,
        criterion="distance",
    )
    return pl.DataFrame(
        {
            "sample": bundle.samples,
            "clonal_cluster": clonal_clusters,
            "strain_cluster": strain_clusters,
        }
    )


def plot_dendo(
    comps_lf: pl.LazyFrame,
    genome: str,
    sample_to_population: pl.LazyFrame,
    min_comp_len: int = 10000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    linkage_method: str = "average",
    color_map: dict | None = None,
    inches_per_sample: float = 0.15,
    font_size: int = 8,
    color_threshold: float = 0.03,
    clonal_cluster_threshold: float = 99.93,
    strain_cluster_threshold: float = 99.8,
    title: str | None = None,
    include_fraction_null: bool = False,
):
    """Plot a left-oriented dendrogram for one genome with optional null-fraction bars."""
    bundle = _prepare_similarity_matrix(
        comps_lf,
        genome=genome,
        min_comp_len=min_comp_len,
        impute_method=impute_method,
        max_null_samples=max_null_samples,
        linkage_method=linkage_method,
    )
    sample_population = _resolve_population_mapping(sample_to_population, bundle.samples)
    sample_population_dict = dict(zip(sample_population["sample_1"], sample_population["population"]))
    color_map, _ = _resolve_population_colors(sample_population, color_map=color_map)

    fig_height = max(20, len(bundle.samples) * inches_per_sample)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    dendro = dendrogram(
        bundle.linkage_matrix,
        ax=ax,
        labels=[f"{sample}_{sample_population_dict.get(sample, 'Not Assigned')}" for sample in bundle.samples],
        orientation="left",
        color_threshold=color_threshold,
        above_threshold_color="gray",
        leaf_font_size=font_size,
        distance_sort="descending",
    )

    ordered_samples = [bundle.samples[idx] for idx in dendro["leaves"]]
    tick_colors = [color_map.get(sample_population_dict.get(sample, "Not Assigned"), "black") for sample in ordered_samples]
    for tick_label, color in zip(ax.get_yticklabels(), tick_colors):
        tick_label.set_color(color)

    ax.axvline(1 - clonal_cluster_threshold / 100, color="red", linestyle="--", linewidth=1)
    ax.axvline(1 - strain_cluster_threshold / 100, color="blue", linestyle="--", linewidth=1)
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label=population, markerfacecolor=color, markersize=8)
        for population, color in color_map.items()
    ]
    line_legend_handles = [
        Line2D([0], [0], color="red", linestyle="--", label="Clonal cluster threshold"),
        Line2D([0], [0], color="blue", linestyle="--", label="Strain cluster threshold"),
    ]
    ax.legend(
        handles=legend_handles + line_legend_handles,
        title="Populations",
        loc="upper left",
        bbox_to_anchor=(0.1, 1),
        frameon=False,
        fontsize=font_size,
        title_fontsize=font_size + 1,
    )

    ax.set_title(title or genome, fontsize=14, pad=20)
    ax.set_xlabel("Pop-ANI", fontsize=12)
    current_ticks = ax.get_xticks()
    ani_ticks = (1 - current_ticks) * 100
    ax.set_xticks(current_ticks)
    ax.set_xticklabels([f"{tick:.2f}" for tick in ani_ticks])
    ax.set_ylabel("Samples", fontsize=12)
    ax.invert_yaxis()
    fig.tight_layout()

    if include_fraction_null:
        leaf_positions = np.array([tick.get_position()[1] for tick in ax.get_yticklabels()])
        bar_values = (
            pl.DataFrame({"sample_1": ordered_samples}).join(
                bundle.null_fraction,
                on="sample_1",
                how="left",
            )
            .get_column("null_fraction")
            .fill_null(0.0)
            .to_list()
        )
        divider = make_axes_locatable(ax)
        ax_bar = divider.append_axes("right", size="30%", pad=2.5)
        bar_height = (np.diff(leaf_positions).mean() * 0.8) if len(leaf_positions) > 1 else 8.0
        ax_bar.barh(leaf_positions, bar_values, height=bar_height)
        ax_bar.set_ylim(ax.get_ylim())
        ax_bar.set_yticks([])
        ax_bar.set_xlim(0, 1)
        ax_bar.set_xlabel("Null fraction", fontsize=font_size)
        fig.tight_layout()

    return fig


def get_clustermap(
    comps_lf: pl.LazyFrame,
    genome: str,
    sample_to_population: pl.LazyFrame,
    min_comp_len: int = 10000,
    impute_method: float = 97.0,
    max_null_samples: int = 500,
    linkage_method: str = "average",
    color_map: dict | None = None,
):
    """Return a seaborn clustermap for one genome."""
    bundle = _prepare_similarity_matrix(
        comps_lf,
        genome=genome,
        min_comp_len=min_comp_len,
        impute_method=impute_method,
        max_null_samples=max_null_samples,
        linkage_method=linkage_method,
    )
    sample_population = _resolve_population_mapping(sample_to_population, bundle.samples)
    color_map, row_colors = _resolve_population_colors(sample_population, color_map=color_map)
    groups = sample_population.get_column("population").unique().sort().to_list()
    qualitative_palette = [color_map[group] for group in groups]

    grid = sns.clustermap(
        bundle.clustermap_data.to_pandas().set_index("sample_1"),
        figsize=(30, 30),
        row_linkage=bundle.linkage_matrix,
        col_linkage=bundle.linkage_matrix,
        xticklabels=True,
        yticklabels=True,
        row_colors=row_colors,
        col_colors=row_colors,
    )
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xmajorticklabels(), fontsize=0.1)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_ymajorticklabels(), fontsize=0.1)
    legend_handles = [
        mpatches.Patch(color=color, label=label)
        for label, color in zip(groups, qualitative_palette)
    ]
    grid.ax_heatmap.legend(
        handles=legend_handles,
        title="Population",
        title_fontsize=16,
        fontsize=14,
        handlelength=2.5,
        handleheight=2,
        bbox_to_anchor=(-0.15, 1),
        loc="lower left",
    )
    return grid
