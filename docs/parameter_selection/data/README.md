# Simulation result provenance

These CSVs are a selected subset of the sequencing-error calibration pilot. They contain measured
outputs, not synthetic chart values. Each pair row identifies its technology, nominal depth, seed,
truth class and filter settings. Counts are integers; figures average independent pairs at each depth.

## Inputs and software

- ZipStrain 1.0.4, commit `e3d04ec90667dd73029f99fe743a3a0991b8ada6` (before the default p-value change).
- Reference: [E. coli K12 U00096.3](https://www.ebi.ac.uk/ena/browser/view/U00096.3).
- Three 10,000-base windows starting at zero-based offsets 4,260,000; 3,170,000; 270,000.
  Named gc0, gc1, gc2; GC fractions 0.4706, 0.5115, 0.5356. Only positions 2001–8000
  (one-based, inclusive) of each window are evaluated.
- A reference FASTA SHA256: `bbd7ecd8b7004cb1317ed6ebefa6c3959f4dc7a1f018a0ec59eb04e4d858cf43`.
- B: 60 random substitutions per evaluated window (mutation seed81273), 180 total / 18,000 positions.
- [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) 2.5.8: HS25,
  paired150bp, fragment mean350/SD50.
- [Badread](https://github.com/rrwick/Badread) commit `e2ccccc608c8132fe22844fcc6c9a544fbb44a77`:
  nanopore2023 error and qscore models; identity20,3 or25,3; length4000,2000.
  Junk/random reads, chimeras, glitches and adapters disabled.
- Independent library seeds11,29 for A and100011,100029 for B. Lower depths reuse subsets
  of each source library; they are not independent replicates across depth.
- Bowtie2 2.5.5 default end-to-end; minimap2 2.31-r1302 map-ont; samtools1.21.
- Both samples mapped to A; no competing references. Read ANI0.95, mapQ0, all-mapped inclusion,
  baseQ13 and reference BAQ. Comparison requires retained coverage5 in both samples.

## Settings tested and included data

We tested six assumed error rates (0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), five minimum
allele frequencies (0, 0.005, 0.01, 0.02, 0.05), and two `p_threshold` values (0.05, 1e-6).
That makes 60 combinations. The null-model table covers depth up to 12,000, beyond the highest
depth observed here. The full experiment produced 54 raw profile datasets, 3,180 pair-and-setting
rows, and 540 mixture-and-setting rows.

`simulation_pairs.csv` includes an assumed error rate of 0.001, minimum allele frequencies of 0
or 0.01, both p-values, and all available depths. The `different`, `identical_A`, and `identical_B`
labels distinguish true A/B differences from same-template controls. Seed 0 denotes a control
comparing libraries 11 and 29, not a third independent simulation. At 1,000×, Illumina has two
A/B pairs; ONT Q20 has only seed 11; ONT Q25 has no data.

Important columns:

- `sites`, `callable_sites`: evaluated positions and positions passing pairwise retained coverage.
- `raw_joint_cov5`: positions meeting coverage5 before allele filtering.
- `true_differences`:180 for A/B;0 for identical controls.
- `recovered_differences`, `masked_differences`, `uncallable_differences`: partition true differences.
- `true_diff_callable`: recovered + masked; `same_callable`: callable unchanged positions.
- `false_differences`: unchanged sites incorrectly called different.
- `popani`, `true_ani_callable`: observed ANI and truth restricted to the callable positions, in percent.

`simulation_mixtures.csv` includes an assumed error rate of 0.001, a minimum allele frequency of
0.01, and both p-values. Mixtures at a nominal depth of 300× combine A/B reads from source seed 11
at B fractions of 0.01, 0.05, and 0.10. Read-name hashing chooses
whole reads/pairs, not independently sampled position counts. These are reused libraries, not new
independent simulations. `retained_minor_sites` is out of180; `raw_minor_present` records minor-base
presence before allele filtering. `error_allele_sites` counts callable invariant sites retaining any
incorrect allele, out of `nonvariant_callable` (17,820 here). It is not a pairwise false-difference count.
`median_observed_minor_fraction` is the median raw minor-base fraction at true variant positions;
finite read sampling means it need not equal the nominal mixture proportion.

## Validation and regeneration

Production ZipStrain Poisson thresholds and allele filters were checked against an independent
indexed threshold calculation for all 60 settings. Production popANI was checked against presence-mask
intersection. Duplicate positions, output-key uniqueness, truth partitions and mixture count bounds
were validated. No false differences occurred in this finite pilot; this does not imply zero real-world risk.

From the repository root:

```bash
python Scripts/plot_parameter_selection_simulation.py
# Optional visual-QA PNGs; requires Kaleido and Chrome:
python Scripts/plot_parameter_selection_simulation.py --preview-dir /tmp/parameter-selection-previews
```

The script regenerates Figures 4–6 from these CSVs; it does not rerun simulation or mapping. Raw
BAMs, results from the full range of settings, and simulation scripts are not bundled with the
tutorial. The original experiment artifacts were kept outside the production repository in
`/tmp/zipstrain-error-calibration-20260913`; that is a temporary local directory, not a durable
public download. The committed CSVs are the durable source for the displayed figures. Figures 1–3
have separate provenance in the original tutorial.
