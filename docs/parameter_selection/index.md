# Keeping popANI honest across depth

popANI can creep upward as sequencing depth goes up, far enough to push a genuine same-strain pair over
the "same clone" line for no biological reason.

This is one of the **profiling parameter** guides. It covers two choices you make together: how
strict to build the **sequencing-error null model** (its error rate, plus an optional `min_freq`
post-filter), and which **ANI metric** to pair it with (`popANI`, `cosANI`, or `conANI`). These choices affect both depth stability and sensitivity to real differences. The rest of the page shows the effect,
measures three ways to remove it, and recommends what ZipStrain ships with.

!!! note "Test data"

    | | |
    |---|---|
    | **Pair** | SRR30025063 × SRR30025566 |
    | **Genome** | GCF_900475675.1 |
    | **Relationship** | Known same strain |
    | **Depth axis** | f005 / f025 / f050 / f100 |

    A single reference pair that we know carries the same strain is downsampled to four read
    fractions. This tests depth sensitivity while holding the underlying samples fixed. Sampling can still
    change callable positions and detection of genuine minor alleles; this pair alone does not
    establish which ANI estimate is closest to biological truth.

## The problem


popANI calls a position "shared" when the two samples have **at least one allele in common**. Deeper
sequencing turns up more alleles, both real rare variants and leftover sequencing error, so more
positions clear that low bar and the score creeps toward 100.

To pin this down, we ran a simple experiment. We took **one pair of samples that we know carry the
same strain** and subsampled their raw reads to four increasing depths (`f005 → f100`). That gives us
a depth axis while the biology stays completely fixed: the same pair, the same genome, just more or
fewer reads. This measures depth dependence, but cannot by itself separate error alleles from genuine rare
alleles detected only at higher depth. On this pair, popANI climbs steadily from **99.9227** at the shallowest depth to **99.9920**
at full depth, crossing the `99.99` "same clone" line for no reason other than the reads we added.
Everything that follows is about removing that artefact: flattening the depth curve while leaving the
real same-versus-different-strain signal alone.

## Seeing it: ANI across depth

Each line is one ANI method. A flat line means the method is depth-invariant; a rising line is the
inflation. The dashed and dotted curves are the compare-time `min_freq` references. The y-axis is
zoomed into the top fraction of a percent so the small differences are actually visible.

<div class="zs-tabs" id="depthtabs"><div class="zs-tabbar"><button class="zs-tab active" data-t="t1c">Both samples deepen together</button><button class="zs-tab" data-t="t1a">Sample 1 varies</button><button class="zs-tab" data-t="t1b">Sample 2 varies</button></div><div class="zs-panel active" id="t1c"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1c_matched_depths.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Matched downsampling. As both samples get deeper (f005 to f100), popANI climbs while cosANI, conANI and the compare-time min_freq references stay flat. Each small panel is one null-model setting.</p></div><div class="zs-panel" id="t1a"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1a_sample1_varies.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Only Sample 1 is downsampled, with Sample 2 held at full depth.</p></div><div class="zs-panel" id="t1b"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1b_sample2_varies.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Only Sample 2 is downsampled while Sample 1 stays at full depth. It is the mirror image, which confirms the effect is symmetric.</p></div></div>

## The real-sample study: two routes to flatter curves

You can attack the bias at two stages: a stricter sequencing-error **null model** (the rows) or a
stricter **ANI metric** (the columns). The grid below shows the residual depth bias for every
combination, and darker means flatter.

<div class="zs-plotwrap"><iframe class="zs-plot" src="fig2_residual_bias_heatmap.html" style="width:860px;height:700px" loading="lazy"></iframe></div>
<p class="zs-cap">Residual depth bias, measured as max minus min ANI across depth. Reading down a column shows what the null model buys you at a fixed metric; reading across a row shows what the metric buys you at a fixed null model. Both routes reach a near-flat floor on this pair. Similar depth stability does not establish equivalent biological accuracy or prove that combining filters is unnecessary.</p>

Two settings are enough on their own: a null model built at a **1% error rate** (span 0.009) and
**cosANI at a threshold of 0.01** (span 0.002). Either one pulls the bias down from 0.069 into the
near-flat range on this pair; no independent measurement-noise baseline was established.

## The real-sample study: little loss of compared positions

A fix that flattens the curve by throwing away data would be a hollow win. That is not what happens
here. Every correction leaves the number of compared positions basically untouched (around 1.71M),
but retained positions do not measure retained minor-allele information. The simulation below tests
that additional cost explicitly.

<div class="zs-plotwrap"><iframe class="zs-plot" src="fig3_fixes_bias_vs_breadth.html" style="width:1260px;height:770px" loading="lazy"></iframe></div>
<p class="zs-cap">Each fix on its own: bias removed on the left (lower is better) against positions kept on the right (higher is better). The right-hand bars are all about the same height, so none of the fixes costs meaningful breadth on this pair.</p>

## Checking accuracy: simulated reads with known differences

A flat depth curve is useful, but it is not enough: a filter could flatten the curve by removing
real minor alleles. We therefore added a second experiment with **known sequence differences**, so
we can distinguish recovered differences, differences hidden by shared errors, and sites that lack
sufficient coverage. This complements the real-sample study above rather than replacing it.

!!! note "Simulation design"

    | | |
    |---|---|
    | **Reference** | Three 10-kb windows from E. coli K12 U00096.3, spanning 47–54% GC |
    | **Truth** | 180 substitutions in 18,000 evaluated positions: **99% ANI** |
    | **Reads** | ART ([Huang et al., 2012](https://doi.org/10.1093/bioinformatics/btr708)) HiSeq2500 paired 150-bp; Badread ([Wick, 2019](https://doi.org/10.21105/joss.01316)) nanopore2023 with Q20 and Q25 read-quality distributions |
    | **Replication** | Two independent simulated sample pairs per technology; lower depths are nested read subsets |
    | **Depths** | 10×, 30×, 100×, 300×; also 1,000× for Illumina and one ONT Q20 pair |
    | **Fixed filters** | Read ANI ≥ 0.95; baseQ ≥ 13; mapQ ≥ 0; all-mapped inclusion; reference-based BAQ enabled; comparison coverage ≥ 5 after allele filtering |
    | **Sweep** | Six error rates × five frequency floors × two p-values: **60 parameter combinations** |

Reads were aligned with Bowtie2 ([Langmead & Salzberg, 2012](https://doi.org/10.1038/nmeth.1923))
for Illumina and minimap2 `map-ont` ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191))
for ONT, then processed with SAMtools ([Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008))
and ZipStrain's production pileup parser, allele filtering and popANI expression.
The evaluation excludes 2 kb at either end of each window. These are small controlled simulations,
not a full metagenomic benchmark. Q20/Q25 are simulator settings, not named ONT chemistries.

### Stricter p-values recover hidden differences

Here **masked** means a true difference has sufficient coverage in both samples but is called
shared because an erroneous allele survives filtering. It is distinct from an **uncallable** site.
The comparison holds `error_rate=0.001` and `min_freq=0.01` fixed and changes only the null-model
p-value between `0.05` and `0.000001`.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Simulation: masked true differences across sequencing depths" src="fig4_simulation_masking.html" style="width:1100px;height:550px" loading="lazy"></iframe></div>
<p class="zs-cap">Mean masked differences per 180 known substitutions. Error bars span the two simulation replicates, not a confidence interval. ONT Q20 at 1,000× has only one replicate; ONT Q25 was not tested at 1,000×. Hover for replicate counts and uncallable differences. At 10×, additional differences are uncallable and are not counted as masked.</p>

At **100×**, where all 18,000 positions were callable:

| Simulation | Masked with p = 0.05 | popANI with p = 0.05 | Masked with p = 1e-6 | popANI with p = 1e-6 |
|---|---:|---:|---:|---:|
| Illumina HiSeq2500 | 4 / 180 | 99.022% | 0 / 180 | 99.000% |
| ONT Q20 | 24 / 180 | 99.133% | 0 / 180 | 99.000% |
| ONT Q25 | 9 / 180 | 99.050% | 0 / 180 | 99.000% |

These are replicate averages. For Illumina, the stricter setting recovered all 180 differences
at every tested depth from 30× through 1,000×, in both pairs. Identical-template controls produced
100% popANI on callable positions across the parameter sweep; there were also no false differences
at unchanged sites in the pure-strain comparisons. This finite pilot does not establish a zero
false-positive rate in real samples.

!!! example "How one erroneous read hides a difference"

    At one simulated C-versus-T position, sample A had **78 C + 1 erroneous T**, while sample B
    had **72 T**. With p = 0.05, the erroneous T survived and popANI called a match. With p = 1e-6,
    it was removed and the true difference was recovered. The 1% frequency floor alone cannot
    remove that T: one read out of 79 is **1.27%**.

### The frequency floor matters at high depth

<div class="zs-plotwrap"><iframe class="zs-plot" title="Simulation: effect of a frequency floor at 1000-fold coverage" src="fig5_simulation_frequency_floor.html" style="width:1100px;height:490px" loading="lazy"></iframe></div>
<p class="zs-cap">One ONT Q20 pair at 1,000×, holding error_rate = 0.001 and p = 0.05 fixed. Removing the frequency floor hid 94 of 180 differences; a 1% floor reduced this to three. All evaluation positions were callable.</p>

Without a frequency floor, this pair reported **99.522%** ANI instead of 99%. With `min_freq=0.01`,
it reported **99.017%**. With the floor enabled, masking did **not** progressively worsen with depth;
in this pilot the largest mean masking occurred around 100×. The coverage-dependent effect therefore
depends on the filtering settings, not just on the definition of popANI.

The error rate seen by the null model is also not the simulator's total error rate. At 300×,
residual substitution fractions after alignment, read/base filtering and BAQ were approximately
**0.078% for Illumina, 0.161% for ONT Q20, and 0.055% for ONT Q25**. These are before null-model
filtering and exclude indels; a global mean alone does not capture context-specific errors.

### The cost: minor alleles, even when breadth is unchanged

We mixed reads from the two templates to create 1%, 5% and 10% minor strains at nominal 300×.
Each mixture has 180 known minor-allele sites. These mixtures reuse the first simulation's reads
and are not additional independent replicates.

<div class="zs-plotwrap"><iframe class="zs-plot" title="Simulation: true minor-allele retention under two null-model thresholds" src="fig6_simulation_minor_alleles.html" style="width:1100px;height:550px" loading="lazy"></iframe></div>
<p class="zs-cap">Retained true minor-allele sites divided by 180, with error_rate = 0.001 and min_freq = 0.01 fixed. Hover for exact counts and the number of invariant positions retaining an erroneous allele. Sampling around a 1% cutoff already removes some genuine 1%-frequency alleles; stricter null filtering removes more.</p>

For the **1% ONT Q20 mixture**, tightening the p-value reduced genuine minor-allele retention
from **84/180 to 26/180**, while reducing invariant sites with retained error alleles from
**125 to 15 out of 17,820**. At **5%**, both settings retained all 180 ONT Q20 minor alleles;
Illumina retained 180 with p = 0.05 and 179 with p = 1e-6. All 10% minor alleles
survived both settings in these simulations.

This is why unchanged comparison breadth is not proof that a correction is free. It can preserve
positions while removing biologically meaningful alleles. Genuine shared alleles in mixed populations
are intentional popANI matches, not necessarily errors to eliminate.

## Recommendation

!!! tip "Defaults for strain comparison"

    Use **error_rate = 0.001**, **p_threshold = 0.000001**, and **min_freq = 0.01**.
    These defaults are supported by the controlled pilot for recovering fixed differences, but
    can discard genuine rare alleles. They are not a universal technology-specific calibration.

    Choose **popANI, conANI or cosANI** according to the biological question, not only which makes
    the depth curve flattest. The real-pair study's 1% error-rate and cosANI observations
    remain useful, but do not establish that the two approaches are interchangeable.

    Filtering at profiling time is irreversible in the saved counts. If you change null-model
    parameters, rebuild the model and rerun profiling. Use `--force-prepare` to regenerate
    auto-prepared assets.

    For rare-variant work, validate sensitivity at the allele frequencies and depths of interest
    before loosening the filters. A 1% frequency floor is not a guarantee of reliable 1%-allele detection.

!!! note "Reading the real-sample figures"

    Figures 1–3 show the real-pair results. Their p-value was not recorded on this page
    and should not be inferred from the recommended defaults.
    Compare-time frequency-filter curves are experimental references, not a supported CLI option.

## References

- Huang W, Li L, Myers JR, Marth GT. (2012). **ART: a next-generation sequencing read simulator.** *Bioinformatics*, 28(4), 593–594. [doi:10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).
- Wick RR. (2019). **Badread: simulation of error-prone long reads.** *Journal of Open Source Software*, 4(36), 1316. [doi:10.21105/joss.01316](https://doi.org/10.21105/joss.01316).
- Langmead B, Salzberg SL. (2012). **Fast gapped-read alignment with Bowtie 2.** *Nature Methods*, 9(4), 357–359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923).
- Li H. (2018). **Minimap2: pairwise alignment for nucleotide sequences.** *Bioinformatics*, 34(18), 3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191).
- Danecek P, et al. (2021). **Twelve years of SAMtools and BCFtools.** *GigaScience*, 10(2), giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).
