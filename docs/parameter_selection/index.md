# How sequencing depth affects popANI

When comparing a genome between two samples, popANI calls a position shared when both samples have enough coverage and retain at least one allele in common. At greater sequencing depth, rare alleles become easier to detect. However, this also increases the chance of observing sequencing errors.

This guide explores how sequencing depth affects the accuracy of popANI. The final goal is to select the best parameters for accurate ANI estimation. First, we compare the same sample pair at several read
depths to see how its ANI estimate changes. Then we use simulated reads with known differences to
check whether settings that steady the estimate also preserve real variation. The settings we
examine are the null model's assumed error rate and `p_threshold`, the profile-time allele-frequency
cutoff (`min_freq`), and the ANI method (`popANI`, `cosANI`, or `conANI`).

## How ZipStrain's null model works

Profiling first counts the A, C, G, and T bases observed at each position after read and base
quality filtering. ZipStrain then asks whether each base count could plausibly be explained by
sequencing error alone. This step treats all four bases alike; it does not need to know the
reference base.

At a position with original A+C+G+T depth **N**, the model expects **N × `error_rate` / 3**
errors of any *one* incorrect base. The division by three assumes that an error is equally likely
to produce each of the other three bases. A Poisson distribution gives the largest count still
compatible with errors at the chosen `p_threshold`. A base is retained only when its observed
count is **strictly greater** than that ceiling **and** its count divided by the original depth
is at least `min_freq`. Counts that fail either test become zero; surviving counts are not reduced
by an estimated number of errors.

For example, at depth 100 with the defaults below, the model expects about 0.033 errors of any
one wrong base and sets the error ceiling to 3. A count of 4 survives both filters (4% frequency),
while a count of 3 does not survive the Poisson test.

| Profiling option | Default | What changing it does |
|---|---:|---|
| `--error-rate` | `0.001` (0.1%) | A higher assumed residual substitution-error rate raises the error ceiling and removes more low-count alleles. |
| `--p-threshold` | `0.000001` | A smaller value demands stronger evidence that a count is not just error, so filtering becomes stricter. |
| `--min-freq` | `0.01` (1%) | A higher value removes more low-frequency alleles. Set it to `0` to disable this frequency test, but not the Poisson test. |
| `--max-total-reads` | `50,000` | Sets the highest **per-position depth** covered by the precomputed null-model table. It does not change the cutoff below that depth; profiling errors if a position exceeds it. |

The model assumes residual substitution errors in the counted reads are approximately independent,
have a common rate, and are spread evenly across the three incorrect bases. It does not estimate
that rate from the sample or explicitly account for context-dependent or correlated errors,
mapping mistakes, or insertions and deletions. Real data can violate these assumptions, and a
strict filter can remove genuine minor alleles along with errors. The defaults are therefore a
starting point for strain comparison, not a calibration for every sequencing technology.

ZipStrain currently supports the Poisson model only. If you provide a prebuilt `--null-model`,
the `--error-rate` and `--p-threshold` options do not change that file; rebuild the model and
rerun profiling to use new settings.

!!! note "Test data"

    | | |
    |---|---|
    | **Pair** | SRR30025063 × SRR30025566 |
    | **Genome** | GCF_900475675.1 |
    | **Relationship** | An example genome with popANI similarity near 99.99% cutoff |
    | **Reads retained** | 5%, 25%, 50%, or 100% (`f005`, `f025`, `f050`, `f100` in the figures) |

    We downsampled one or both samples from a pair known to carry the same strain. The source
    samples stayed the same, but downsampling can change which positions have enough coverage and
    which rare alleles are observed. This experiment measures sensitivity to depth; it cannot tell
    us which ANI estimate is closest to the true sequence difference.

## The problem


If a position has different alleles in two samples, a single surviving error that matches the
other sample can make popANI count it as shared. More reads provide more chances to observe such
an allele. Genuine low-frequency alleles can also appear as depth increases, so a rising score is
not proof that every additional shared allele is an error.

With Sample 2 held at full depth and Sample 1 increased from 5% to 100% of its reads, **unfiltered
popANI rose from 99.9227% to 99.9920%**. That crosses a 99.99% cutoff without changing the source
samples. When both samples were downsampled together, the corresponding curve started at **99.9304%**
and ended at the same full-depth value. The figures below keep these two comparisons separate.

## Seeing it: ANI across depth

Flatter lines change less over the depths tested, but are not necessarily more accurate. The y-axis
is zoomed in so small differences are visible. Dashed and dotted lines show **experimental
comparison-time** frequency filters for context; they are not a ZipStrain CLI option. The other
frequency cutoffs were applied **during profiling**. Each chart shows all ten profiling settings
at the same y-axis scale.

<div class="zs-tabs" id="depthtabs">
  <div class="zs-tabbar">
    <button class="zs-tab active" data-t="t1a">Sample 1 varies</button>
    <button class="zs-tab" data-t="t1b">Sample 2 varies</button>
    <button class="zs-tab" data-t="t1c">Both samples vary</button>
  </div>
  <div class="zs-panel active" id="t1a">
    <div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 1A: Sample 1 varies" src="fig1a_sample1_varies.html" style="width:100%;height:990px" loading="lazy"></iframe></div>
    <p class="zs-cap">Sample 2 stays at full depth while Sample 1 varies. The unfiltered popANI curve rises from 99.9227% to 99.9920%.</p>
  </div>
  <div class="zs-panel" id="t1b">
    <div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 1B: Sample 2 varies" src="fig1b_sample2_varies.html" style="width:100%;height:990px" loading="lazy"></iframe></div>
    <p class="zs-cap">Sample 1 stays at full depth while Sample 2 varies. The unfiltered curve follows a different pattern; the effect is not symmetric for this pair.</p>
  </div>
  <div class="zs-panel" id="t1c">
    <div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 1C: Both samples vary" src="fig1c_matched_depths.html" style="width:100%;height:990px" loading="lazy"></iframe></div>
    <p class="zs-cap">Both samples use the same read fraction at each point. Unfiltered popANI rises from 99.9304% to 99.9920%; the other ANI methods change less, but they are not perfectly flat.</p>
  </div>
</div>

## Reducing the depth effect in the real samples

Figure 2 compares profiling settings (rows) with ANI methods (columns). Each number is the
difference between the highest and lowest ANI values as **Sample 1** moves through the four read
fractions while Sample 2 stays at full depth. It is measured in **ANI percentage points**: a lower
number, shown with a lighter color, means less change with depth. A row labeled
"0.1% error, 1% cutoff" uses that assumed error rate and a 1% profile-time allele-frequency
cutoff. "Unfiltered profiles" means no profiling error filter was applied. The p-value used for these real-sample
plots was not recorded, so it should not be assumed to match today's default.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 2: ANI range by profile filters and ANI method" src="fig2_residual_bias_heatmap.html" style="width:100%;height:490px" loading="lazy"></iframe></div>
<p class="zs-cap">ANI range across the four Sample 1 read fractions. Lighter cells and smaller numbers mean less variation with depth. The rows combine an assumed error rate and a profile-time frequency cutoff; the columns compare ANI methods.</p>

For this pair, unfiltered popANI changed by **0.069 percentage points**. Using cosANI with a
threshold of 0.01 reduced that range to **0.002 points**. Keeping popANI but using a null model
with a **1% assumed error rate** and no frequency cutoff reduced it to **0.009 points**. These are
different approaches, and similar stability does not mean they recover true differences equally
well. The simulation below tests that question directly.

## What happens to the number of compared positions?

The settings in Figure 3 leave the number of positions compared in this pair close to **1.71
million**. That is useful, but it does not mean the same alleles survive filtering at those
positions. The simulated mixtures below show what happens to genuine minor alleles.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 3: ANI range and positions compared" src="fig3_fixes_bias_vs_breadth.html" style="width:100%;height:620px" loading="lazy"></iframe></div>
<p class="zs-cap">Left: ANI range across depths (smaller is steadier). Right: positions compared (larger means more positions). The chart includes changes to the ANI method, experimental comparison-time frequency filters, and profiling settings; some profiling settings combine an error-rate assumption with a frequency cutoff.</p>

## Checking accuracy with simulated reads

A steadier ANI estimate is not automatically a better one: a filter could steady the curve by
discarding real minor alleles. In a separate simulation, we know exactly where two sequences
differ. That lets us count differences that ZipStrain recovers, differences hidden by surviving
error alleles, and positions without enough coverage. This checks something the real-sample
experiment cannot.

!!! note "Simulation design"

    | | |
    |---|---|
    | **Reference** | Three 10-kb windows from E. coli K12 U00096.3, spanning 47–54% GC |
    | **Known differences** | 180 substitutions in 18,000 evaluated positions: **99% true ANI** |
    | **Reads** | ART ([Huang et al., 2012](https://doi.org/10.1093/bioinformatics/btr708)) HiSeq2500 paired 150-bp; Badread ([Wick, 2019](https://doi.org/10.21105/joss.01316)) nanopore2023 with Q20 and Q25 read-quality distributions |
    | **Repeated simulations** | Two independently simulated sample pairs per technology; lower depths reuse subsets of the same reads |
    | **Depths** | 10×, 30×, 100×, 300×; also 1,000× for Illumina and one ONT Q20 pair |
    | **Settings held fixed** | Read ANI ≥ 0.95; base quality ≥ 13; mapping quality ≥ 0; all mapped reads included; reference-based BAQ enabled; comparison coverage ≥ 5 after allele filtering |
    | **Range of settings tested** | Six assumed error rates, five minimum allele frequencies, and two `p_threshold` values: **60 combinations** |

Reads were aligned with Bowtie2 ([Langmead & Salzberg, 2012](https://doi.org/10.1038/nmeth.1923))
for Illumina and minimap2 `map-ont` ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191))
for ONT, then processed with SAMtools ([Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008))
and ZipStrain's production pileup parser, allele filtering and popANI expression.
We excluded the first and last 2 kb of each window from the evaluation. These are small,
controlled simulations, not a full metagenomic benchmark. Q20 and Q25 describe simulated
read-quality distributions, not specific ONT chemistries.

### A stricter null-model cutoff reveals hidden differences

Here a **masked difference** is a known difference that is masked because of presence of sequencing errors. An **uncallable** position
does not have enough retained coverage to compare. We held `error_rate=0.001` and
`min_freq=0.01` fixed and changed only `p_threshold` from `0.05` to `0.000001` (one in a million).
The smaller value removes more low-count alleles that could be explained by sequencing error.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 4: Masked true differences across sequencing depths" src="fig4_simulation_masking.html" style="width:100%;height:475px" loading="lazy"></iframe></div>
<p class="zs-cap">Average number of known differences hidden among 180 substitutions. Error bars show the range of two simulated pairs, not a confidence interval. At 1,000×, ONT Q20 has only one pair and ONT Q25 was not tested. Some differences at 10× lack enough coverage; they are not counted as masked.</p>

At **100×**, where all 18,000 positions were callable:

| Simulation | Masked with p = 0.05 | popANI with p = 0.05 | Masked with p = 1e-6 | popANI with p = 1e-6 |
|---|---:|---:|---:|---:|
| Illumina HiSeq2500 | 4 / 180 | 99.022% | 0 / 180 | 99.000% |
| ONT Q20 | 24 / 180 | 99.133% | 0 / 180 | 99.000% |
| ONT Q25 | 9 / 180 | 99.050% | 0 / 180 | 99.000% |

The table shows averages across the two simulated pairs. With the smaller `p_threshold`, both
Illumina pairs recovered all 180 differences at every tested depth from 30× through 1,000×.
Identical-template controls scored 100% popANI at positions with enough coverage, and we observed
no false differences at unchanged positions in these pure-strain simulations. That does not imply
a zero false-positive rate in real data.

!!! example "How one erroneous read hides a difference"

    At one simulated C-versus-T position, Sample A had **78 C reads and one erroneous T read**;
    Sample B had **72 T reads**. With `p_threshold=0.05`, the erroneous T remained, so popANI
    counted a match. At `p_threshold=1e-6`, it was removed and the true difference was recovered.
    The 1% frequency cutoff alone could not remove that T: one read out of 79 is **1.27%**.

### A minimum allele frequency matters at high depth

<div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 5: Effect of a frequency cutoff at 1000-fold coverage" src="fig5_simulation_frequency_floor.html" style="width:100%;height:415px" loading="lazy"></iframe></div>
<p class="zs-cap">One simulated ONT Q20 pair at 1,000×, with `error_rate=0.001` and `p_threshold=0.05`. Without a minimum allele frequency, 94 of 180 known differences were hidden; a 1% cutoff reduced that to three. All 18,000 evaluated positions had enough coverage.</p>

Without the cutoff, this pair reported **99.522% popANI** instead of the known 99% ANI. With
`min_freq=0.01`, it reported **99.017%**. For the tested ONT Q20 pairs with that cutoff, the
average number of hidden differences was highest at 100× rather than increasing steadily with
depth. The depth effect depends on filtering as well as on how popANI defines a match.

The error rate reaching the null model is not the simulator's overall error rate. At 300×,
substitution errors remaining after alignment, read and base filtering, and BAQ were approximately
**0.078% for Illumina, 0.161% for ONT Q20, and 0.055% for ONT Q25**. These values precede
null-model filtering and exclude insertions and deletions. A single average also misses errors
that vary by sequence context.

### The trade-off: genuine rare alleles can be lost

We mixed reads from the two templates to create samples with 1%, 5%, or 10% of a minor strain
at a nominal depth of 300×. Each mixture has 180 positions where the minor allele is known.
The mixtures reuse reads from the first simulation; they are not additional independent repeats.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Figure 6: True minor-allele retention under two null-model thresholds" src="fig6_simulation_minor_alleles.html" style="width:100%;height:475px" loading="lazy"></iframe></div>
<p class="zs-cap">Percentage of 180 known minor-allele positions retained, with `error_rate=0.001` and `min_freq=0.01`. Hover for exact counts and for unchanged positions that retained an erroneous allele. Sampling alone can push a real 1% allele below the 1% cutoff; stricter null filtering can remove still more.</p>

For the **1% ONT Q20 mixture**, tightening the p-value reduced genuine minor-allele retention
from **84/180 to 26/180**, while reducing invariant sites with retained error alleles from
**125 to 15 out of 17,820**. At **5%**, both settings retained all 180 ONT Q20 minor alleles;
Illumina retained 180 with p = 0.05 and 179 with p = 1e-6. All 10% minor alleles
survived both settings in these simulations.

The earlier figures counted how many positions remained available for comparison. This result
shows why that count is not enough: a position can remain while a real minor allele disappears.
When two mixed populations genuinely share an allele, popANI is *supposed* to count it as a match;
not every shared allele is an error to remove.

## Recommendation

!!! tip "Defaults for strain comparison"

    For routine strain comparison, use the current defaults: **`error_rate=0.001`**,
    **`p_threshold=0.000001`**, and **`min_freq=0.01`**. In these simulations, they recovered
    fixed differences well, but they can remove genuine rare alleles. They are not a calibration
    for every sequencing technology or sample.

    Choose **popANI, conANI, or cosANI** for the biological question, not just for the flattest
    curve. The real-sample experiment showed that a 1% assumed error rate or a cosANI threshold
    of 0.01 can reduce depth dependence; it did not show that they preserve real alleles equally.

    Filtering at profiling time is irreversible in the saved counts. If you change null-model
    parameters, rebuild the model and rerun profiling. Use `--force-prepare` to regenerate
    auto-prepared assets.

    For rare-variant work, check sensitivity at the allele frequencies and depths that matter to
    your study. A 1% cutoff does **not** guarantee that true 1% alleles will be detected.

!!! note "Reading the real-sample figures"

    Figures 1–3 show the real-sample pair. Their `p_threshold` value was not recorded and should
    not be inferred from today's defaults. The dashed and dotted comparison-time frequency curves
    are experimental references, not a supported CLI option. The committed
    [pair results](data/simulation_pairs.csv) and [mixture results](data/simulation_mixtures.csv)
    underlie Figures 4–6.

## References

- Huang W, Li L, Myers JR, Marth GT. (2012). **ART: a next-generation sequencing read simulator.** *Bioinformatics*, 28(4), 593–594. [doi:10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).
- Wick RR. (2019). **Badread: simulation of error-prone long reads.** *Journal of Open Source Software*, 4(36), 1316. [doi:10.21105/joss.01316](https://doi.org/10.21105/joss.01316).
- Langmead B, Salzberg SL. (2012). **Fast gapped-read alignment with Bowtie 2.** *Nature Methods*, 9(4), 357–359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923).
- Li H. (2018). **Minimap2: pairwise alignment for nucleotide sequences.** *Bioinformatics*, 34(18), 3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191).
- Danecek P, et al. (2021). **Twelve years of SAMtools and BCFtools.** *GigaScience*, 10(2), giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).
