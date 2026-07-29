# Expected output

This page describes every file ZipStrain writes, the columns in each, and how to read them. Profile outputs are written as parquet. Convert any parquet table to CSV explicitly with `zipstrain utilities parquet-to-csv` when text output is needed. The examples below come from [Tutorial #3](./Tutorial.md) (three replicates of the ZymoBIOMICS standard, Sylph route).

---

## `zipstrain profile`

A profile run writes one folder per sample directly under the run directory, plus shared `profiling_assets/` (auto-generated inputs + logs):

```text
out_profile/
├── SRR12324251/
│   ├── SRR12324251_profile.parquet        # per-position base counts
│   ├── SRR12324251_genome_stats.parquet   # per-genome summary
│   ├── SRR12324251_gene_stats.parquet     # per-gene summary
│   ├── SRR12324251_SNVs.parquet           # SNPs vs. reference
│   └── intermediate_files/                # scratch — safe to ignore
├── SRR12324252/ …
└── profiling_assets/                      # null model, bed, contract, logs
```

### `<sample>_genome_stats.parquet`

One row per genome: how well it was covered, whether it is present, and its identity.

| Column | Meaning |
|--------|---------|
| `genome` | Genome identifier (a GTDB accession on the Sylph route) |
| `coverage` | Mean read depth across the genome (over its full length) |
| `coverage_median`, `coverage_std` | Median and std. dev. of depth across covered positions |
| `breadth` | Fraction of the genome covered by ≥ 1 read |
| `genome_length` | Genome length (bp) |
| `gap_mean`, `gap_std` | Mean / std. dev. of the gap length between consecutive mapped reads |
| `5x_cov_sites` | Number of positions covered ≥ 5× |
| `heterogeneity` | Mean per-site nucleotide diversity (within-sample microdiversity) |
| `ber` | Breadth-error ratio: observed breadth ÷ breadth expected at this coverage (~1 for an evenly-covered genome) |
| `fug` | Read-distribution uniformity (fraction of non-excess gaps); ~0.632 under a random distribution |
| `reads_mapped` | Reads mapped to the genome |
| `ref_ani` | **Population** ANI of the sample to the reference genome (%) — a position counts as matching if the reference allele is present among the reads |
| `conANI_reference` | **Consensus** ANI to the reference (%) — a position counts as matching only if the reference base is the sample's consensus (majority) base. Lower than `ref_ani` at polymorphic sites. Both are reported only when `--reference-fasta` is used |
| `SNS_count`, `SNV_count` | Number of fixed substitutions (`SNS`) and polymorphic sites (`SNV`+`con_SNV`+`pop_SNV`) — the per-genome totals of the SNV table. Reference route only |
| `presence` | Automated `present` / `absent` call (see [Tutorial #3](./Tutorial.md#interpreting-the-benchmarks)) |
| `genome_taxonomy` | GTDB lineage — **Sylph route only** |

Example (one present genome):

```text
genome           coverage  breadth  ber     fug     reads_mapped  ref_ani  presence  genome_taxonomy
GCF_000006945.2  380.34    0.9211   0.9211  0.6467  6579418       99.09    present   d__Bacteria;…;s__Salmonella enterica
```

Here breadth ≈ 0.92 at 380× coverage and `presence = present` — a genuine member of the community. Trace genomes appear with breadth < 0.05 and `presence = absent`.

### `<sample>_gene_stats.parquet`

One row per gene, populated only when gene annotations are supplied (`profile --gene-fasta`, or a gene FASTA on the Sylph route). Empty otherwise.

| Column | Meaning |
|--------|---------|
| `genome` | Genome the gene belongs to |
| `gene` | Gene identifier |
| `length` | Gene length (bp) |
| `breadth` | Fraction of the gene covered by ≥ 1 read |
| `coverage` | Mean read depth across the gene |
| `ref_ani` | ANI of the sample to the reference for this gene (%) |

### `<sample>_profile.parquet`

The core output: one row per reference position with the observed base counts. This is what comparisons are computed from. It is large (tens of millions of rows), and profile runs do not write CSV sidecars.

| Column | Meaning |
|--------|---------|
| `chrom` | Scaffold / contig name |
| `genome` | Genome the scaffold belongs to (`NA` if unbinned) |
| `gene` | Gene at this position (`NA` if none / unannotated) |
| `pos` | Position on the scaffold (**1-based**) |
| `A`, `C`, `G`, `T` | Read counts that pass both the Poisson error ceiling and `--min-freq` allele-frequency filter |
| `ref_base_bitmask` | Reference base as a bitmask: **A = 1, C = 2, G = 4, T = 8** (`0` = reference base unknown). Only present when `--reference-fasta` was used |

Example:

```text
chrom                    genome           gene  pos      A   C   G   T   ref_base_bitmask
GCA_016925055.1__CP0704… GCA_016925055.1  NA    3892095  0   51  0   0   0
GCA_016925055.1__CP0704… GCA_016925055.1  NA    3892097  51  0   0   0   0
```

The null model stores the largest count still compatible with sequencing error, so an allele survives only when its count is strictly greater than that ceiling. `--min-freq` is applied against the original A+C+G+T depth before any allele is removed and defaults to `0.01` (1%); set it to `0` to disable frequency filtering. The generated model supports coverage through 50,000 by default; profiling raises an error instead of silently altering a position when its depth exceeds the model.

### `<sample>_SNVs.parquet`

Every covered position that is **divergent** from the reference — either polymorphic or a substitution — classified in the style of inStrain's `SNVs.tsv`. Skipped without `--reference-fasta` (or with `--no-snvs`). Because the profile counts are sequence-error adjusted, a base with `count > 0` is a "passing" allele.

| Column | Meaning |
|--------|---------|
| `chrom`, `genome`, `gene`, `pos` | Location (`pos` is 1-based) |
| `position_coverage` | Total (error-adjusted) depth at the position |
| `allele_count` | Number of passing alleles (bases with count > 0) |
| `ref_base`, `con_base`, `var_base` | Reference, consensus (majority), and top variant base |
| `ref_freq`, `con_freq`, `var_freq` | Frequencies of those bases over `position_coverage` |
| `A`, `C`, `G`, `T` | Per-base counts |
| `class` | SNP classification (below) |
| `ref_base_bitmask` | Reference base bitmask (A=1, C=2, G=4, T=8) |

`class` values (with `alleles` = passing bases, `con` = consensus, `ref` = reference):

- **`SNS`** — one allele, `con ≠ ref`: a fixed substitution.
- **`SNV`** — ≥2 alleles, `con == ref`: reference is the majority, with a minor variant.
- **`con_SNV`** — ≥2 alleles, `con ≠ ref`, reference **still present** among the alleles.
- **`pop_SNV`** — ≥2 alleles, `con ≠ ref`, reference **absent**: population-divergent (these are what pull `ref_ani`/popANI below 100%).

Monomorphic reference matches are not emitted. The per-genome `SNS_count` and `SNV_count` (= `SNV`+`con_SNV`+`pop_SNV`) in `genome_stats` are the totals of these.

!!! note "Coordinates are 1-based (comparing with inStrain)"
    ZipStrain reports `pos` as **1-based** (the samtools/mpileup convention), whereas inStrain's `SNVs.tsv` uses **0-based** `position`. When intersecting the two tables, add 1 to inStrain's `position` (or subtract 1 from ZipStrain's `pos`) before matching, or every position will appear to disagree.

---

## `zipstrain compare`

Both methods (standard and matrix) write a single comparison table at the top of the run directory: `all_comparisons.parquet` (+ `.csv`), or `all_gene_comparisons.parquet` with `--compare-genes`.

### `all_comparisons.parquet`

One row per genome, per pair of samples.

| Column | Meaning |
|--------|---------|
| `sample_1`, `sample_2` | The two samples being compared |
| `genome` | Genome compared between them |
| `total_positions` | Positions covered in **both** samples (≥ `--min-cov`) — the basis of the comparison |
| `share_allele_pos` | Of those, positions that match under the selected ANI method |
| `genome_ani` | Genome-wide ANI (%) = `share_allele_pos` / `total_positions` × 100. The parquet metadata key `zipstrain_compare_ani_method` records whether this is `popani`, `conani`, or `cosani_<threshold>` |
| `max_consecutive_length` | Longest run of consecutive shared positions (an identical-by-state / IBS measure) |

With gene ranges available, genome-level gene identity summary columns are also added by the default `--calculate all`: `shared_genes_count`, `identical_gene_count`, `perc_id_genes`.

!!! note "ANI method is stored in the parquet header"
    `--ani-method conani` or `--ani-method cosani_<threshold>` changes how the ANI match indicator is computed, but it does not add method-specific output columns. Check `zipstrain_compare_ani_method` in the parquet metadata to know what `genome_ani` represents. Both standard and matrix comparison support `popani`, `conani`, and `cosani_<threshold>`.

Example (same-strain replicates → 100% popANI):

```text
sample_1     sample_2     genome           total_positions  share_allele_pos  genome_ani  max_consecutive_length
SRR12324251  SRR12324252  GCA_031316495.1  3983156          3983156           100.0           768237
```

**Interpreting popANI.** `genome_ani` near **100.0** means the same strain is present in both samples. A common same-strain threshold is **≥ 99.999%** — stringent enough to distinguish strains that diverged only a few years apart. Because ZipStrain uses *population* ANI (a position matches if any read supports the shared base), co-existing minor variants confirm rather than break a match; see [Tutorial #3](./Tutorial.md#interpreting-the-benchmarks) for the popANI-vs-consensus-ANI discussion.
