"""zipstrain.visualize
========================
This module provides statistical analysis and visualization functions for profiling and compare operations.
"""

import polars as pl
import plotly.graph_objects as go
import seaborn as sns
import numpy as np
from itertools import chain, combinations
from collections import defaultdict
import matplotlib.patches as mpatches
import pandas as pd



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
    Calculate strain sharing between populations based on popANI between genomes in their profiles.
    Strain sharing between two samples is defined as the ratio of genomes passing a strain similarity threshold over the total number of genomes in each sample.
    So, for two samples A and B, the strain sharing is defined as (Note the assymetric nature of the calculation):
    strain_sharing(A, B) = (number of genomes in A and B passing the strain similarity threshold) / (number of genomes in A)
    strain_sharing(B, A) = (number of genomes in A and B passing the strain similarity threshold) / (number of genomes in B)
    
    Args:
        comps_lf (pl.LazyFrame): LazyFrame containing the gene profiles of the samples.
        breadth_lf (pl.LazyFrame): LazyFrame containing the genome breadth information.
        sample_grouping (pl.LazyFrame): LazyFrame containing the sample grouping information.
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
        (pl.col("genome_pop_ani") >= strain_similarity_threshold)
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
    Calculate the fraction of identical genes vs  popANI for a given genome and two samples in any possible combination of populations.
    Args:
        genome (str): The genome to calculate the fraction of identical genes vs popANI for.
        population_1 (str): The first population to compare.
        population_2 (str): The second population to compare.
        sample_to_population (pl.LazyFrame): LazyFrame containing the sample to population mapping.
        comps_lf (pl.LazyFrame): LazyFrame containing the gene profiles of the samples
    Returns:
        pl.LazyFrame: LazyFrame containing the fraction of identical genes vs popANI information for
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
        pl.col("genome_pop_ani")
    ).collect(engine="streaming")

def plot_identical_frac_vs_popani(df:pl.DataFrame,
                                  genome:str,
                                  title:str="Fraction of Identical Genes vs popANI for <GENOME>",
                                  xaxis_title:str="Genome-Wide popANI",
                                  yaxis_title:str="Fraction of Identical Genes",
                                  ):
    """
    Plot the fraction of identical genes vs popANI for a given genome and two samples in any possible combination of populations.
    Args:
        df (pl.DataFrame): DataFrame containing the fraction of identical genes vs popANI information.
        title (str, optional): Title of the plot. Defaults to "Fraction of Identical Genes vs popANI".
        xaxis_title (str, optional): Title of the x-axis. Defaults to "popANI".
        yaxis_title (str, optional): Title of the y-axis. Defaults to "Fraction of Identical Genes".
    Returns:
        go.Figure: Plotly figure containing the fraction of identical genes vs popANI plot.
    """
    fig = go.Figure()
    for group, perc_id_genes, genome_pop_ani in zip(df["relationship"], df["perc_id_genes"], df["genome_pop_ani"]):
        fig.add_trace(go.Scatter(
            x=genome_pop_ani,
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

def plot_clustermap(
    comps_lf:pl.LazyFrame,
    genome:str,
    sample_to_population:pl.LazyFrame,
    min_comp_len:int=10000,
    impute_method:str|float=97.0,
    max_null_samples:int=200,
    color_map:dict|None=None,
):
    """
    Plot a clustermap for the given genome and its associated samples.
    Args:
        comps_lf (pl.LazyFrame): LazyFrame containing the comparison data.
        genome (str): The genome to plot.
        sample_to_population (pl.LazyFrame): LazyFrame containing the sample to population mapping.
    Returns:
        go.Figure: Plotly figure containing the clustermap.
    """
    # Filter the comparison data for the specific genome
    comps_lf_filtered = comps_lf.filter(
        (pl.col("genome") == genome) & (pl.col("total_positions") > min_comp_len)
    ).select(
        pl.col("sample_1"),
        pl.col("sample_2"),
        pl.col("genome_pop_ani"),
    )
    comps_lf_filtered_oposite = comps_lf_filtered.select(
        pl.col("sample_2").alias("sample_1"),
        pl.col("sample_1").alias("sample_2"),
        pl.col("genome_pop_ani"),
    )
    # Combine the filtered data with its opposite pairs
    comps_lf_filtered = pl.concat([comps_lf_filtered, comps_lf_filtered_oposite])
    # Make a synthetic table for similarity of samples with themselves of all samples in sample_1 and sample_2 but each sample exists only once
    self_similarity =(
    pl.concat([
        comps_lf_filtered.select(pl.col("sample_1").alias("sample_1")),
        comps_lf_filtered.select(pl.col("sample_2").alias("sample_1"))
    ])
    .unique()
    .sort("sample_1").with_columns(
        pl.col("sample_1").alias("sample_2"),
        pl.lit(100.0).alias("genome_pop_ani"),
    )
    )
    
    # Combine the self similarity with the filtered data
    comps_lf_filtered = pl.concat([self_similarity, comps_lf_filtered]).collect()
    # Pivot the data for the clustermap
    clustermap_data = comps_lf_filtered.pivot(
        index="sample_1",
        columns="sample_2",
        values="genome_pop_ani"
    )
    # We want to make this a similarity matrix, so we need to frop null values, have sample_1 and sample_2 as index and columns as we
    # Create the clustermap
    exclude_samples=clustermap_data.null_count().transpose(include_header=True, header_name="column", column_names=["null_count"]).filter(pl.col("null_count")>max_null_samples)["column"].to_list()
    # Only include rows and cols not in exclude_samples
    clustermap_data = clustermap_data.filter(~pl.col("sample_1").is_in(exclude_samples))
    clustermap_data = clustermap_data.select(*[col for col in clustermap_data.columns if col not in exclude_samples])
    if isinstance(impute_method, str):
        pass # To be implemented later
    elif isinstance(impute_method, (int, float)):
        clustermap_data = clustermap_data.fill_null(impute_method)
    sample_to_population = clustermap_data.select(pl.col("sample_1")).join(
        sample_to_population.collect(),
        left_on="sample_1",
        right_on="sample",
        how="left")
    sample_to_population_dict = dict(zip(sample_to_population["sample_1"], sample_to_population["population"]))
    if color_map is None:

        num_categories = sample_to_population["population"].n_unique()
        groups= sample_to_population["population"].unique().sort().to_list()
        qualitative_palette = sns.color_palette("hls", num_categories)
        row_colors = [qualitative_palette[groups.index(sample_to_population_dict[sample])] for sample in clustermap_data["sample_1"]]
        col_colors = [qualitative_palette[groups.index(sample_to_population_dict[sample])] for sample in clustermap_data.columns if sample != "sample_1"]
    else:
        groups= list(color_map.keys())
        qualitative_palette= list(color_map.values())
        row_colors = [color_map[sample_to_population_dict[sample]] for sample in clustermap_data["sample_1"]]
        col_colors = [color_map[sample_to_population_dict[sample]] for sample in clustermap_data.columns if sample != "sample_1"]
    fig = sns.clustermap(
        clustermap_data.to_pandas().set_index("sample_1"),
        figsize=(30, 30),
        xticklabels=True, 
        yticklabels=True,
        row_colors=row_colors,
        col_colors=col_colors
    )
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xmajorticklabels(), fontsize=0.1)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_ymajorticklabels(), fontsize=0.1)
    legend_handles = [mpatches.Patch(color=color, label=label)
                  for label, color in zip(groups, qualitative_palette)]
    fig.ax_heatmap.legend(handles=legend_handles,
                          title='Population',
                        title_fontsize=16,   # bigger title
                        fontsize=14,         # bigger labels
                        handlelength=2.5,    # wider color boxes
                        handleheight=2,
                    bbox_to_anchor=(-0.15, 1), loc="lower left")
    return fig