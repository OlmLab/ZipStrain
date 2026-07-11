# Expected output

This page describes every file ZipStrain writes, the columns in each, and how to read them. All tables are written as parquet, and — when small enough — a matching `.csv` (see `--no-csv` / `--force-csv`). The examples below come from [Tutorial #3](./Tutorial.md) (three replicates of the ZymoBIOMICS standard, Sylph route).

---

## `zipstrain profile`

A profile run writes one folder per sample directly under the run directory, plus shared `profiling_assets/` (auto-generated inputs + logs):

```text
out_profile/
├── SRR12324251/
│   ├── SRR12324251_profile.parquet        # per-position base counts
│   ├── SRR12324251_genome_stats.parquet   # per-genome summary (+ .csv)
│   ├── SRR12324251_gene_stats.parquet     # per-gene summary (+ .csv)
│   ├── SRR12324251_SNVs.parquet           # SNPs vs. reference (+ .csv)
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

The core output: one row per reference position with the observed base counts. This is what comparisons are computed from. It is large (tens of millions of rows) and has no CSV companion by default.

| Column | Meaning |
|--------|---------|
| `chrom` | Scaffold / contig name |
| `genome` | Genome the scaffold belongs to (`NA` if unbinned) |
| `gene` | Gene at this position (`NA` if none / unannotated) |
| `pos` | Position on the scaffold (**1-based**) |
| `A`, `C`, `G`, `T` | Read counts of each base at this position |
| `ref_base_bitmask` | Reference base as a bitmask: **A = 1, C = 2, G = 4, T = 8** (`0` = reference base unknown). Only present when `--reference-fasta` was used |

Example:

```text
chrom                    genome           gene  pos      A   C   G   T   ref_base_bitmask
GCA_016925055.1__CP0704… GCA_016925055.1  NA    3892095  0   51  0   0   0
GCA_016925055.1__CP0704… GCA_016925055.1  NA    3892097  51  0   0   0   0
```

### `<sample>_SNVs.parquet`

SNPs relative to the reference: the subset of `profile` positions that pass coverage (`--snv-min-cov`) and no longer retain the reference allele after sequence-error adjustment. Same columns as `profile`. Skipped without `--reference-fasta` (or with `--no-snvs`).

Read an SNV like this: `ref_base_bitmask = 2` (reference C) with `T = 5` and `C = 0` means the reference base C was replaced by T — a SNP.

```text
chrom                    genome           gene  pos  A  C  G  T  ref_base_bitmask
GCA_018376995.1__JAHAAF… GCA_018376995.1  NA    632  0  0  0  5  2      # ref C -> T
GCA_018376995.1__JAHAAF… GCA_018376995.1  NA    634  0  0  5  0  1      # ref A -> G
```

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
| `share_allele_pos` | Of those, positions where the two samples share an allele (the popANI numerator) |
| `genome_pop_ani` | **Population** ANI (%) = `share_allele_pos` / `total_positions` × 100 — a position matches if any read supports a shared base |
| `share_consensus_pos` | Positions where the two samples' **consensus** (majority) base agrees |
| `consensus_SNPs` | Positions where the consensus bases differ (= `total_positions` − `share_consensus_pos`) |
| `genome_con_ani` | **Consensus** ANI (%) = `share_consensus_pos` / `total_positions` × 100. Stricter than `genome_pop_ani`; the gap between them flags shared populations with a different dominant strain |
| `max_consecutive_length` | Longest run of consecutive shared positions (an identical-by-state / IBS measure) |

The consensus columns come from the `conani` calculation, which is included in the default `--calculate all`. With gene ranges available (or `--compare-genes`), gene-level columns are added: `shared_genes_count`, `identical_gene_count`, `perc_id_genes`.

!!! note "conANI is standard-method only"
    `genome_con_ani` / `share_consensus_pos` / `consensus_SNPs` are produced by `--method standard`. The `--method matrix` store keeps only per-base **presence** (which base is present, not its count) to make popANI a fast boolean operation, and consensus can't be recovered from presence — so matrix mode reports popANI + IBS (+ gene ANI) but not conANI.

Example (same-strain replicates → 100% popANI):

```text
sample_1     sample_2     genome           total_positions  share_allele_pos  genome_pop_ani  max_consecutive_length
SRR12324251  SRR12324252  GCA_031316495.1  3983156          3983156           100.0           768237
```

**Interpreting popANI.** `genome_pop_ani` near **100.0** means the same strain is present in both samples. A common same-strain threshold is **≥ 99.999%** — stringent enough to distinguish strains that diverged only a few years apart. Because ZipStrain uses *population* ANI (a position matches if any read supports the shared base), co-existing minor variants confirm rather than break a match; see [Tutorial #3](./Tutorial.md#interpreting-the-benchmarks) for the popANI-vs-consensus-ANI discussion.
