# Ambient correction: attempted at scale, rejected by its own acceptance test

**Outcome: the correction does not work and its output must not be used.**
Recorded in full because this is the first step in the project that could
certify itself, and it certified itself as failed.

---

## 1. CellBender is not affordable here — measured, not assumed

Installed (0.3.2, CPU torch) and run on one library. Two facts:

- Its checkpointing is incompatible with torch 2.13
  (`cannot pickle 'weakref.ReferenceType'`); `--checkpoint-mins 100000`
  works around it.
- **322 seconds per epoch on 4 CPU cores.** At 50 epochs that is ~4.5 h
  for one library and **~300 h for all 73**.

One library is being run to completion as a benchmark. Full-dataset
CellBender is out of reach in this container.

## 2. What was run instead

SoupX's non-expressed-gene estimator, built on this project's own
finding. The soup profile is measured directly from empty droplets
(this dataset ships the full 737,280-barcode whitelist, so they are
available), and the contamination fraction ρ is read off genes a lineage
cannot express — immunoglobulin for NK/CD8T/CD4T/Myeloid, granzymes for
B. Applied at the (library × lineage × condition) level, which is exact
because summation is linear.

### The measurement is sound and worth keeping

| lineage | ρ (median) |
|---|---|
| **B** | **0.205** |
| **NK** | **0.151** |
| Myeloid | 0.076 |
| CD4T | 0.064 |
| CD8T | 0.058 |

**NK carries 2.6× the soup fraction of CD8 T.** That is independent
confirmation of the RNA-content result: NK has ~25% less endogenous RNA
on the tumour side, so the same soup occupies more of its profile.

## 3. The acceptance test, and how it rejected the method

Two versions of the test were run, and **they disagree — which is the
lesson.**

**Weak version (mine): does the marker's *within-sample fraction* fall?**
Held-out markers, never used for fitting, fell 95–100% in every lineage.
**Passed.**

**Strong version: does the marker's *tumour-vs-normal log2FC* go to
zero?** Its true value is exactly zero, so this is the criterion that
matters for a differential analysis.

| | uncorrected | corrected | corrected + re-filtered | target |
|---|---|---|---|---|
| \|log2FC\| immunoglobulin **in NK** | 3.985 | **4.216** | (genes dropped) | → 0 |
| \|log2FC\| granzyme/perforin **in B** | 0.862 | **1.314** | **1.440** | → 0 |
| NK median SE | 0.120 | 0.199 | **0.239** | unchanged |
| B median SE | 0.088 | 0.286 | **0.290** | unchanged |
| stop rule S3 | PASS | **STOP (2.12×)** | **STOP (2.59×)** | PASS |

**Both known-zero controls got worse, and S3 broke.** Re-deriving the
detection floor on corrected counts (arm 3) correctly removed all eight
immunoglobulin genes, but did not rescue the granzymes and made the
variance worse.

## 4. Diagnosis

Two mechanisms, both intrinsic to the method as implemented:

**(a) ρ is estimated from three marker genes per sample, so ρ is noisy.**
Subtracting `ρ × total × soup` from every gene injects that per-sample
noise into every gene at once — a sample-level random effect that
inflates residual variance across the board. NK's SE doubles; B's
triples. CellBender avoids this by modelling the whole droplet
population jointly rather than fitting one scalar per sample.

**(b) Near-zero residuals give unstable ratios.** For genes that are
almost entirely soup, subtraction leaves a small residual floored at
zero, and the log2 ratio of two noisy near-zero values is wild — which
is why immunoglobulin's |log2FC| in NK went *up*.

## 5. What this leaves standing

- **The ambient measurement stands** (ρ by lineage above, and the
  shared-emulsion calibration: Ig |log2FC| in NK 5.67 → 1.68, in myeloid
  4.68 → 0.21). Ambient is real, large, and NK-heavy.
- **The corrected pseudobulk must not be used.** `results/pseudobulk_
  decontaminated.npz` and the `quadrants_decontaminated*.csv` files are
  retained as the record of a failed attempt, not as analysis inputs.
- **Q2 remains unevaluable.** The blockade was ambient; the ambient
  correction available here is worse than the disease.

## 6. What would actually work

1. **CellBender with adequate compute** (GPU, or many CPU-hours). The
   acceptance test is already written and would validate it.
2. **A better ρ estimator** — many more marker genes, or ρ pooled across
   samples within a library so it stops being a per-sample random
   effect. This is the cheap fix and was not attempted, because two
   failed arms is the point at which iterating blindly becomes the
   error this project keeps making.
3. **Modelling ambient as an offset in the DE model** rather than
   subtracting counts, which avoids the near-zero-residual problem
   entirely.

## 7. The methodological point

This is the first analysis step in the project with a built-in,
known-truth acceptance test. It failed that test, twice, and the failure
was visible immediately rather than after external argument.

Note also which version of the test caught it. The within-sample
marker-fraction check **passed** while the differential check **failed**.
A decontamination method can remove the soup's *level* convincingly and
still corrupt the *contrast*, which is the only thing a differential
analysis reads. **The acceptance test has to be stated in the same units
as the claim.**

---

# Addendum: the two authorised follow-up runs

Two runs were pre-registered and executed after the failure above, on
the explicit condition that the line closes afterwards regardless of
outcome. Both are recorded here in full.

## 8. Run 1 — is ρ's per-sample noise the binding constraint?

Section 4 offered diagnosis (a): ρ is fitted per sample from three
markers, so subtraction injects a sample-level random effect into every
gene. The test written before the run:

> ρ pooled across conditions within a library, so ρ is no longer a
> per-sample quantity.
> **SE returns to baseline** → (a) is fixed, (b) is binding, move to the
> regression approach.
> **SE does not return** → diagnosis (a) was wrong, stop the line.

`scripts/rho_variants.py` computed three estimators in one pass over the
73 libraries; `scripts/rho_variant_se.py` read out the standard errors
on a gene panel fixed on the uncorrected counts, so the arms differ only
in the counts.

**Median SE** (`results/rho_variant_se.csv`):

| arm | NK | CD8T | CD4T | B | Myeloid |
|---|---|---|---|---|---|
| baseline (no correction) | 0.120 | 0.071 | 0.059 | 0.088 | 0.071 |
| ρ per sample (the failed arm) | 0.199 | 0.096 | 0.091 | 0.286 | 0.075 |
| **ρ pooled across conditions** | **0.193** | 0.095 | 0.091 | **0.270** | 0.075 |
| ρ from ~20 markers | 0.238 | 0.133 | 0.079 | 0.516 | 0.234 |

**Pooling ρ changed essentially nothing** — NK 0.199 → 0.193, B 0.286 →
0.270, against a baseline of 0.120 and 0.088. The pre-registered branch
that fires is the second one: **diagnosis (a) was wrong.**

The real mechanism is arithmetic, not statistical. Subtracting
`ρ × total × soup` moves soup-dominated genes into a lower-expression
regime, and voom assigns variance as a function of expression level, so
those genes get larger precision-weighted standard errors. The variance
inflation is not noise imported from ρ; it is the honest consequence of
the counts being smaller. No ρ estimator fixes that.

Three further things the run showed, none of which were the question
asked:

- **ρ estimated from ~20 lineage-foreign genes is roughly double the
  three-marker estimate** (NK 0.151 → 0.295, B 0.205 → 0.476). Two
  defensible marker panels disagree by 2×, so **a single scalar ρ does
  not describe this contamination.** Some "foreign" genes are more
  soup-loaded than others, which is exactly what the regression in §9
  assumes and the subtraction model denies.
- **The wide-marker arm passes S3 while being the worst arm.** Its
  NK:CD4T SE ratio is 1.30 against baseline's 1.69 — better-looking —
  because every lineage degraded together. A ratio-based stop rule is
  blind to uniform degradation. Recorded as a limitation of S3.
- **NK's ρ varies between the two conditions within a library far more
  than any other lineage does**: within-(library × lineage) SD of ρ
  across conditions is 0.120 for NK against a median ρ of 0.151, versus
  0.005 for CD4T. Ambient's *differential* component — the only part
  that distorts a paired contrast — is concentrated in NK. This is the
  third NK-specific technical effect (see §10).

## 9. Run 2 — regress ambient out of log2FC instead of subtracting counts

The alternative, and the better idea: never modify a count. Model the
observed log2FC as a function of each gene's **soup fraction** in that
lineage,

    f_g = ρ × total × soup_g / observed_g   (clipped to [0,1])

fit a monotone binned curve across genes, and take residuals. Structural
advantages: nothing is floored at zero, no per-sample quantity enters,
and the standard errors are untouched by construction.

The acceptance test is stated in the units of the claim — differential,
out of sample. The known-zero genes (immunoglobulin in NK/CD8T/CD4T/
Myeloid, granzymes in B) are **excluded from the fit**, so their
residuals are a genuine held-out check.

**SE, as predicted, is identical to baseline**: NK 0.122, CD8T 0.071,
CD4T 0.058, B 0.088, Myeloid 0.071. The variance problem that killed
subtraction is gone.

**The controls** (`results/ambient_regression_summary.csv`):

| lineage | estimated soup fraction of the control genes | \|log2FC\| before | after | change |
|---|---|---|---|---|
| **B** (granzymes) | **1.000** | 0.861 | **0.200** | **−77%** |
| Myeloid (Ig) | 0.767 | 3.422 | 3.025 | −12% |
| CD4T (Ig) | 0.408 | 2.928 | 3.128 | **+7%** |
| **NK** (Ig) | **0.291** | 3.982 | 3.608 | −9% |
| CD8T (Ig) | 0.286 | 3.192 | 3.338 | **+5%** |

Read the first column against the last. These genes are known to be
**pure soup — their true soup fraction is 1.0 in every row.** The method
corrects them almost perfectly in the one lineage where the estimate is
right (B, f = 1.000, −77%) and does nothing or slightly harms where the
estimate is badly low (NK f = 0.291, CD8T f = 0.286). The correction is
doing what it should; it is being told the wrong thing.

**So the binding constraint is neither the regression nor ρ's noise: it
is the per-gene soup fraction.** `f_g` is computed from the empty-droplet
profile, and that profile under-represents immunoglobulin's share of
what leaks into cell-containing droplets — the empty droplets are not a
faithful sample of the soup the cells actually swim in. This is precisely
the quantity CellBender infers jointly instead of assuming, and it is
now the single identified blocker.

Background |log2FC| across all other genes moved down modestly in every
lineage (NK −11%, B −19%, CD8T −2%), which is consistent with removing a
real component rather than shrinking everything.

**Verdict: the regression is the right shape and cannot be validated
here.** `results/ambient_regression_genes.csv` is retained as a method
record. It is **not** used to recompute quadrants — a correction that
demonstrably works in one lineage out of five would silently
re-differentiate the lineage comparison, which is the exact failure
mode this project is trying to detect.

## 10. Closing the ambient line — and what it did to the framework

CellBender was killed at epoch 7/50 (322 s/epoch → ~300 h for the
cohort). No GPU is waited for: the ambient line closes here.

The line closes, but it did not end empty, and its finding is not about
ambient. Three **independent, NK-specific** technical effects were
measured in this dataset:

1. **Immunoglobulin soup fraction** — NK carries the second-highest ρ of
   any lineage (0.151, behind B's 0.205), and its estimated soup fraction
   at known-zero genes is the *lowest* (0.291), i.e. NK's contamination
   is both large and the worst-characterised.
2. **Differential dissociation stress** — NK's ρ differs between tumour
   and adjacent-normal within the same library by an SD of 0.120 on a
   mean of 0.151, versus 0.005 for CD4T. Only NK's ambient burden is
   condition-dependent, and only the condition-dependent part biases a
   paired contrast.
3. **RNA content loss on the tumour side, in NK alone** — median UMI
   1611 → 1156, log2 ratio −0.323, against +0.061 CD8T, +0.612 B,
   +0.212 myeloid. NK vs the others, **p = 6.0 × 10⁻⁸**.

The four-quadrant framework rests on one assumption: that technical
effects hit all lineages alike, land in Q3, and are therefore
discardable — which is what licenses reading Q4 as biology. **All three
effects above are lineage-specific and NK-specific. They are Q4-shaped
by construction.** On this dataset the assumption is not merely
unverified; it is falsified in the direction that manufactures the
project's target class.

That is the project's actual result. It is a negative result about the
method, obtained from the method's own controls, and it is worth more
than the 46-gene Q4 list it invalidates.
