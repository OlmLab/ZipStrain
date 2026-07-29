# Keeping popANI honest across depth

popANI creeps upward as sequencing depth goes up, far enough to push a genuine same-strain pair over
the "same clone" line for no biological reason.

This is one of the **profiling parameter** guides. It covers two choices you make together: how
strict to build the **sequencing-error null model** (its error rate, plus an optional `min_freq`
post-filter), and which **ANI metric** to pair it with (`popANI`, `cosANI`, or `conANI`). Get these
right and the score stays stable across sequencing depth. The rest of the page shows the effect,
measures three ways to remove it, and recommends what ZipStrain ships with.

!!! note "Test data"

    | | |
    |---|---|
    | **Pair** | SRR30025063 × SRR30025566 |
    | **Genome** | GCF_900475675.1 |
    | **Relationship** | Known same strain |
    | **Depth axis** | f005 / f025 / f050 / f100 |

    A single reference pair that we know carries the same strain is downsampled to four read
    fractions. Because the biology never changes, any movement in ANI across those fractions is pure
    depth bias, which makes this pair an ideal ruler.

## The problem

popANI calls a position "shared" when the two samples have **at least one allele in common**. Deeper
sequencing turns up more alleles, both real rare variants and leftover sequencing error, so more
positions clear that low bar and the score creeps toward 100.

To pin this down, we ran a simple experiment. We took **one pair of samples that we know carry the
same strain** and subsampled their raw reads to four increasing depths (`f005 → f100`). That gives us
a depth axis while the biology stays completely fixed: the same pair, the same genome, just more or
fewer reads. So anything that moves in popANI across those depths is a depth artefact, not real
signal. On this pair, popANI climbs steadily from **99.9227** at the shallowest depth to **99.9920**
at full depth, crossing the `99.99` "same clone" line for no reason other than the reads we added.
Everything that follows is about removing that artefact: flattening the depth curve while leaving the
real same-versus-different-strain signal alone.

## Seeing it: ANI across depth

Each line is one ANI method. A flat line means the method is depth-invariant; a rising line is the
inflation. The dashed and dotted curves are the compare-time `min_freq` references. The y-axis is
zoomed into the top fraction of a percent so the small differences are actually visible.

<div class="zs-tabs" id="depthtabs"><div class="zs-tabbar"><button class="zs-tab active" data-t="t1c">Both samples deepen together</button><button class="zs-tab" data-t="t1a">Sample 1 varies</button><button class="zs-tab" data-t="t1b">Sample 2 varies</button></div><div class="zs-panel active" id="t1c"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1c_matched_depths.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Matched downsampling. As both samples get deeper (f005 to f100), popANI climbs while cosANI, conANI and the compare-time min_freq references stay flat. Each small panel is one null-model setting.</p></div><div class="zs-panel" id="t1a"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1a_sample1_varies.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Only Sample 1 is downsampled, with Sample 2 held at full depth.</p></div><div class="zs-panel" id="t1b"><div class="zs-plotwrap"><iframe class="zs-plot" src="fig1b_sample2_varies.html" style="width:1410px;height:1600px" loading="lazy"></iframe></div><p class="zs-cap">Only Sample 2 is downsampled while Sample 1 stays at full depth. It is the mirror image, which confirms the effect is symmetric.</p></div></div>

## The headline: two knobs, the same cure

You can attack the bias at two stages: a stricter sequencing-error **null model** (the rows) or a
stricter **ANI metric** (the columns). The grid below shows the residual depth bias for every
combination, and darker means flatter.

<div class="zs-plotwrap"><iframe class="zs-plot" src="fig2_residual_bias_heatmap.html" style="width:860px;height:700px" loading="lazy"></iframe></div>
<p class="zs-cap">Residual depth bias, measured as max minus min ANI across depth. Reading down a column shows what the null model buys you at a fixed metric; reading across a row shows what the metric buys you at a fixed null model. Both routes reach the same near-flat floor, and they overlap: they are correcting the <em>same</em> error, so stacking them buys you very little.</p>

Two settings are enough on their own: a null model built at a **1% error rate** (span 0.009) and
**cosANI at a threshold of 0.01** (span 0.002). Either one pulls the bias down from 0.069 into the
range of ordinary measurement noise.

## The cost: essentially none

A fix that flattens the curve by throwing away data would be a hollow win. That is not what happens
here. Every correction leaves the number of compared positions basically untouched (around 1.71M),
so the choice between them comes down to *where* the correction lives, not to lost data.

<div class="zs-plotwrap"><iframe class="zs-plot" src="fig3_fixes_bias_vs_breadth.html" style="width:1260px;height:770px" loading="lazy"></iframe></div>
<p class="zs-cap">Each fix on its own: bias removed on the left (lower is better) against positions kept on the right (higher is better). The right-hand bars are all about the same height, so none of the fixes costs meaningful breadth on this pair.</p>

## Recommendation

!!! tip "Recommended parameters"

    **Pick one: a 1% null model, or cosANI at 0.01.** Either one removes the depth bias on its own.
    They correct the same error, so there is no reason to apply both.

    **If you want your profiles to stay general, keep the null model loose and use cosANI.** A loose
    null model at profiling time keeps all of the profile's allele detail, so the same profiles work
    for any downstream comparison. You then apply the depth correction later, at compare time, with
    cosANI, and you can change the cosANI threshold whenever you like without re-profiling. Baking a
    strict null model into the profile does the opposite: it throws away information up front and
    locks one level of strictness into the file.

    **Why there is no compare-time `min_freq` option.** A frequency filter applied only at compare
    time reproduces cosANI almost exactly (a depth span of 0.003 versus 0.002 at the 0.01 level). It
    adds nothing new, and because it runs as a custom per-comparison step it slows things down. So its
    job lives in the two places that already do it efficiently: **null-model generation** at profile
    time, and **cosANI** at compare time.

---

<p class="zs-foot">Depth-bias parameter study. One same-strain pair (SRR30025063 × SRR30025566,
GCF_900475675.1), four read fractions, swept over null-model strictness, ANI metric, and a
compare-time frequency filter. The figures are interactive, so hover for exact values.</p>
