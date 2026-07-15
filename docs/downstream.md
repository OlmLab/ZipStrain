# Downstream Analyses

This page documents the downstream analysis tools available in `zipstrain.visualize`. These functions help you go from raw genome-comparison outputs to interpretable figures and tables — things like strain-sharing rates between populations, identity-by-state distributions, dendrograms, and cluster assignments. All examples are in Python. If you're not a programmer, you can point an AI assistant (like Claude or ChatGPT) at this page and ask it to write the code for you — just describe what you want to visualize and paste in your file paths.

## Input conventions

Most comparison-driven functions expect a `pl.LazyFrame` with some or all of these columns:

- `genome`
- `sample_1`
- `sample_2`
- `genome_ani`
- `total_positions`
- `max_consecutive_length`
- `shared_genes_count`
- `identical_gene_count`
- `perc_id_genes`

Population-aware functions also expect a `sample_to_population` lazy frame with:

- `sample`
- `population`

## Common setup

For most downstream analysis, start from a genome-comparison parquet and a sample metadata table:

```python
import polars as pl
from zipstrain import visualize

comps_lf = pl.scan_parquet("example_genome_compare.parquet")
sample_to_population = pl.scan_csv("sample_to_population.csv")
```

Where:

- `example_genome_compare.parquet` is a genome-compare output parquet
- `sample_to_population.csv` contains at least:
  - `sample`
  - `population`

Example metadata file:

```text
sample,population
sample_a,Population_1
sample_b,Population_1
sample_c,Population_2
sample_d,Population_2
```

Use `pl.scan_parquet(...)` and `pl.scan_csv(...)` when a function expects a `LazyFrame`.
Use `pl.read_parquet(...)` when a function expects an eager `DataFrame`.

## Working with outputs

The downstream helpers return a mix of tables and figures.

### If a function returns a `pl.DataFrame`

You can inspect it directly:

```python
print(df.head())
```

You can save it for later:

```python
df.write_parquet("output.parquet")
df.write_csv("output.csv")
```

### If a function returns a Plotly figure

You can open it interactively:

```python
fig.show()
```

You can save it as an interactive HTML file:

```python
fig.write_html("figure.html")
```

Or save a static image if you have the required Plotly image backend installed:

```python
fig.write_image("figure.png")
```

### If a function returns a Matplotlib figure

You can save it as:

```python
fig.savefig("figure.png", dpi=200, bbox_inches="tight")
```

### If a function returns a Seaborn `ClusterGrid`

Save the underlying Matplotlib figure:

```python
grid.fig.savefig("clustermap.png", dpi=200, bbox_inches="tight")
```

## Function catalogue

### `calculate_strainsharing`

- Purpose: quantify asymmetric strain-sharing rates between populations.
- Inputs:
  - `comps_lf`
  - `breadth_lf`
  - `sample_to_population`
  - `min_breadth`
  - `strain_similarity_threshold`
  - `min_total_positions`
- Output:
  - `dict[str, list[float]]`
- Notes:
  - Uses breadth and ANI filters before computing pairwise sharing rates.

What it does:

- This function asks: for each pair of populations, how often do samples appear to share the same strain?
- It returns a Python dictionary where each key is a population pair and each value is a list of sharing rates across sample pairs.
- You will usually pass the result into `plot_strainsharing(...)`.

Example:

```python
breadth_lf = pl.scan_parquet("example_breadth_table.parquet")

strainsharing = visualize.calculate_strainsharing(
    comps_lf=comps_lf,
    breadth_lf=breadth_lf,
    sample_to_population=sample_to_population,
    min_breadth=0.5,
    strain_similarity_threshold=99.9,
    min_total_positions=10000,
)

print(strainsharing.keys())
```

### `plot_strainsharing`

- Purpose: plot strain-sharing rate distributions.
- Inputs:
  - `strainsharingrates`
  - `sample_frac`
  - `title`
  - `xaxis_title`
  - `yaxis_title`
- Output:
  - `plotly.graph_objects.Figure`

What it does:

- This takes the dictionary from `calculate_strainsharing(...)` and turns it into a boxplot.
- It is useful for comparing how much strain sharing you see within and between groups.

Example:

```python
fig = visualize.plot_strainsharing(
    strainsharingrates=strainsharing,
    sample_frac=1.0,
)
fig.show()
fig.write_html("strainsharing.html")
```

### `calculate_ibs`

- Purpose: aggregate IBS-like distributions across within- and between-population comparisons.
- Inputs:
  - `sample_to_population`
  - `comps_lf`
  - `max_perc_id_genes`
  - `min_total_positions`
- Output:
  - `pl.DataFrame`
- Notes:
  - Produces one row per genome and comparison type.

What it does:

- This function reorganizes IBS-like information into a format that is easier to plot.
- It groups comparisons into:
  - within-population
  - between-population
- The result is a `DataFrame` that you can pass to `plot_ibs(...)` or `plot_ibs_heatmap(...)`.

Example:

```python
ibs_df = visualize.calculate_ibs(
    sample_to_population=sample_to_population,
    comps_lf=comps_lf,
    max_perc_id_genes=15,
    min_total_positions=10000,
)

print(ibs_df.head())
ibs_df.write_parquet("ibs_summary.parquet")
```

### `plot_ibs`

- Purpose: plot IBS cumulative distributions for two populations within one genome.
- Inputs:
  - `df`
  - `genome`
  - `population_1`
  - `population_2`
  - `vert_thresh_hor_distance`
  - `num_bins`
  - `title`
  - `xaxis_title`
  - `yaxis_title`
- Output:
  - `plotly.graph_objects.Figure`

What it does:

- This draws the IBS cumulative distributions for two populations within one genome.
- It is most useful when you want to compare whether within-population pairs and between-population pairs have clearly different IBS behavior.

Example:

```python
fig = visualize.plot_ibs(
    df=ibs_df,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    population_1="Population_1",
    population_2="Population_2",
)
fig.show()
fig.write_html("ibs_curve.html")
```

### `plot_ibs_heatmap`

- Purpose: summarize IBS separation across genomes and population pairs.
- Inputs:
  - `df`
  - `vert_thresh`
  - `populations`
  - `num_bins`
  - `min_member`
  - `title`
  - `xaxis_title`
  - `yaxis_title`
- Output:
  - `plotly.graph_objects.Figure`

What it does:

- This compresses the IBS comparison across many genomes into one heatmap.
- Each cell represents how separated the IBS distributions are for one genome and one population pair.

Example:

```python
fig = visualize.plot_ibs_heatmap(
    df=ibs_df,
    min_member=10,
)
fig.show()
fig.write_html("ibs_heatmap.html")
```

### `calculate_identical_frac_vs_popani`

- Purpose: collect genome-wide popANI and identical-gene fractions for two populations.
- Inputs:
  - `genome`
  - `population_1`
  - `population_2`
  - `sample_to_population`
  - `comps_lf`
  - `min_shared_genes_count`
  - `min_total_positions`
- Output:
  - `pl.DataFrame`

What it does:

- This collects genome-wide popANI values and identical-gene fractions for two chosen populations.
- It is helpful when you want to see whether near-identical genomes also tend to have high identical-gene fractions.

Example:

```python
identical_df = visualize.calculate_identical_frac_vs_popani(
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    population_1="Population_1",
    population_2="Population_2",
    sample_to_population=sample_to_population,
    comps_lf=comps_lf,
)

print(identical_df.head())
identical_df.write_csv("identical_vs_popani.csv")
```

### `plot_identical_frac_vs_popani`

- Purpose: scatterplot identical-gene fraction against genome-wide popANI.
- Inputs:
  - `df`
  - `genome`
  - `title`
  - `xaxis_title`
  - `yaxis_title`
- Output:
  - `plotly.graph_objects.Figure`

What it does:

- This turns the output of `calculate_identical_frac_vs_popani(...)` into a scatterplot.
- It is useful for seeing whether samples form obvious within-group and between-group patterns.

Example:

```python
fig = visualize.plot_identical_frac_vs_popani(
    df=identical_df,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
)
fig.show()
fig.write_html("identical_vs_popani.html")
```

### `compute_silhouette_curve`

- Purpose: compute the silhouette curve and peak summary for one genome.
- Inputs:
  - `comps_lf`
  - `genome`
  - `min_comp_len`
  - `impute_method`
  - `max_null_samples`
  - `min_threshold`
  - `peak_prominence`
  - `peak_distance`
- Output:
  - `SilhouetteCurveResult`
- Notes:
  - Expects one genome scope at a time.
  - Null ANI entries are imputed with a numeric ANI value.

What it does:

- This tries many ANI cutoffs and asks: “At which threshold do the resulting clusters look most coherent?”
- The result includes:
  - the full threshold vector
  - the full silhouette vector
  - `curve`: a table with columns `threshold` and `silhouette`
  - `candidate_peaks`: a table with columns `index`, `threshold`, `silhouette`, and `prominence`
  - `best_peak`: a one-row table with columns `index`, `threshold`, `silhouette`, and `source`
- A higher silhouette score usually means cleaner separation between clusters.

Example:

```python
result = visualize.compute_silhouette_curve(
    comps_lf=comps_lf,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    min_comp_len=100000,
    linkage_method="average",
)
print(result.best_peak)
print(result.candidate_peaks)
```

### `plot_silhouette_curve`

- Purpose: plot a precomputed silhouette curve result.
- Inputs:
  - `result`
- Output:
  - `plotly.graph_objects.Figure`

Example:

```python
result = visualize.compute_silhouette_curve(
    comps_lf=comps_lf,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    min_comp_len=100000,
)

fig = visualize.plot_silhouette_curve(result)
fig.show()
fig.write_html("silhouette.html")
```

### `get_silhouette_plot`

- Purpose: convenience wrapper that computes and plots the silhouette curve in one call.
- Notes:
  - Internally this calls `compute_silhouette_curve(...)` and then `plot_silhouette_curve(...)`.
  - It accepts the same `linkage_method` argument, defaulting to `"average"`.

### `get_cluster_assignments`

- Purpose: assign clonal and strain-level clusters from hierarchical clustering.
- Inputs:
  - `comps_lf`
  - `min_comp_len`
  - `impute_method`
  - `max_null_samples`
  - `clonal_cluster_threshold`
  - `strain_cluster_threshold`
- Output:
  - `pl.DataFrame` with columns:
    - `sample`
    - `clonal_cluster`
    - `strain_cluster`
- Notes:
  - This helper expects `comps_lf` to already represent a single genome scope.

What it does:

- This performs hierarchical clustering on one genome’s ANI comparison matrix.
- It then assigns:
  - `clonal_cluster`
  - `strain_cluster`
- This is useful when you want cluster labels that you can join back to metadata and use downstream.

Example:

```python
single_genome_lf = comps_lf.filter(
    pl.col("genome") == "GCF_000269965.1_ASM26996v1_genomic.fna"
)

cluster_df = visualize.get_cluster_assignments(
    comps_lf=single_genome_lf,
    clonal_cluster_threshold=99.93,
    strain_cluster_threshold=99.8,
    linkage_method="average",
)

print(cluster_df)
cluster_df.write_csv("cluster_assignments.csv")
```

Notes:

- `linkage_method` defaults to `"average"`
- this is passed directly to SciPy hierarchical clustering
- if you want stricter cluster semantics, you can also try `"complete"`

### `plot_dendo`

- Purpose: draw a population-colored dendrogram for one genome.
- Inputs:
  - `comps_lf`
  - `genome`
  - `sample_to_population`
  - `min_comp_len`
  - `impute_method`
  - `max_null_samples`
  - `color_map`
  - `inches_per_sample`
  - `font_size`
  - `color_threshold`
  - `clonal_cluster_threshold`
  - `strain_cluster_threshold`
  - `title`
  - `include_fraction_null`
- Output:
  - `matplotlib.figure.Figure`

What it does:

- This plots a dendrogram for one genome.
- Samples are clustered by ANI similarity.
- Leaf labels can be colored by population.
- If `include_fraction_null=True`, it also adds a side bar showing how much missing pairwise ANI information each sample had before imputation.

Example:

```python
fig = visualize.plot_dendo(
    comps_lf=comps_lf,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    sample_to_population=sample_to_population,
    linkage_method="average",
    include_fraction_null=True,
)
fig.savefig("dendrogram.png", dpi=200, bbox_inches="tight")
```

### `get_clustermap`

- Purpose: build a seaborn clustermap for one genome using ANI as similarity.
- Inputs:
  - `comps_lf`
  - `genome`
  - `sample_to_population`
  - `min_comp_len`
  - `impute_method`
  - `max_null_samples`
  - `color_map`
- Output:
  - `seaborn.matrix.ClusterGrid`

What it does:

- This creates a clustered heatmap of ANI values for one genome.
- It is a good “overview plot” when you want to see sample blocks, clusters, and group structure all at once.

Example:

```python
grid = visualize.get_clustermap(
    comps_lf=comps_lf,
    genome="GCF_000269965.1_ASM26996v1_genomic.fna",
    sample_to_population=sample_to_population,
    linkage_method="average",
)
grid.fig.suptitle("ANI clustermap", y=1.02)
grid.fig.savefig("clustermap.png", dpi=200, bbox_inches="tight")
```
