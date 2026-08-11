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
