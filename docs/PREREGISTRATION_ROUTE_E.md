# PREREGISTRATION — X1 Route E (within-compartment NK clone structure)

**Committed before any Route-E statistic is computed.** Read-only after first
commit; amendment requires a version bump in `docs/CHANGELOG.md` with reason and
timestamp.

**Why this exists.** Tier 1 was stopped because the deposited variant set is
ascertained to be blood-present, so any *cross*-compartment statistic answers its
own question (`docs/X1_TIER1_RESULT.md`, DEC-21). Route E uses **only each
compartment's own variants**, so no cross-compartment variant selection enters and
that confound cannot arise. It is the one X1 computation that is currently
runnable.

**What it is and is not.** It compares the **clone structure** of NK *inside*
tumour, inside NILT and inside blood. Structure is not ancestry. This is a
**hint, not a verdict**, and the asymmetric interpretation rule in §4 is fixed
here precisely so that it cannot be read as one later.

---

## 1. Estimator (fixed)

Per donor `d`, per compartment `c`, over NK cells only, using variants selected
**within that compartment's libraries alone**.

**Primary estimator — clonal concentration.** Let clone `i` in compartment `c`
hold `n_i` NK cells out of `N` assigned. Report:

```
Simpson-style concentration   D_c = Σ_i ( n_i (n_i − 1) ) / ( N (N − 1) )
```

`D_c` is the probability that two NK cells drawn without replacement from that
compartment share a clone. It is chosen over Shannon entropy and over "fraction of
cells in expanded clones" for three reasons, all fixed here:

1. it is **unbiased at fixed N** under sampling without replacement, which is what
   makes the downsampling in §2 sufficient rather than approximate;
2. it depends on the **large** clones, which is what "local expansion" means, and
   is insensitive to the singleton tail, which is where clone-calling error lives;
3. it needs **no threshold** for what counts as "expanded".

**Secondary estimators**, reported alongside and never substituted for the primary:
- fraction of assigned NK cells in clones of size ≥ 3
- size of the largest clone / N

**Predictions, fixed before computation:**

| | prediction |
|---|---|
| tumour NK **recruited** from blood | `D_tumour ≈ D_blood` — flat, blood-like |
| tumour NK from a **resident** pool | `D_tumour > D_blood`, and closer to `D_NILT` |

`D_NILT` is the resident-side reference. `D_blood` is the recruited-side
reference. **Both references are required**; a tumour value reported without both
is uninterpretable and is not to be reported.

---

## 2. Downsampling (fixed)

Clonal concentration is strongly N-dependent, and the compartments differ by more
than an order of magnitude (tumour NK 281–509; blood NK up to 7,713). Therefore:

1. Per donor, set `n* = min(NK cells assigned to a clone in tumour, NILT, blood)`.
2. Draw `n*` NK cells **without replacement** from each compartment; recompute `D`.
3. Repeat **200 times** (fixed). Report **median and 2.5–97.5 percentile
   interval** per compartment per donor.
4. **Report per donor. Never pool donors.** Cross-donor comparison of clone
   structure is not meaningful and pooling lets the deepest donor set the answer.
5. Undownsampled values are reported as secondary. Any qualitative disagreement
   between downsampled and raw is reported, not resolved silently.

**Admissibility.** A donor-compartment cell with `n* < 30` clone-assigned NK cells
is **not analysed** — the same 30-cell floor already preregistered for power. If
fewer than 3 donors have all three compartments admissible, **Route E fails power
and is reported as failing power.** No threshold is lowered to produce output.

---

## 3. Clone calling (fixed — inherited unchanged)

Identical to `docs/PREREGISTRATION_TIER1.md` §3 **except** that variants are taken
**per compartment**, not by donor union — which is the entire point of Route E:

- mgatk defaults: strand correlation ≥ 0.65, VMR ≥ 0.01, detected in ≥ 5 cells,
  **within that compartment's libraries only**
- chrM:307–314 excluded
- variants in > 20% of that compartment's cells removed
- heteroplasmy binarised at 0.07
- cells carrying no admitted variant are **excluded and counted**, not dropped silently

**Sensitivity analysis, fixed in advance and run regardless of outcome:** repeat at
≥ 3 cells, and at binarisation 0.05 and 0.10. If the ordering of `D_tumour`
relative to `D_blood` and `D_NILT` is not stable across these, the result is
reported as **unstable**, not as the primary value.

---

## 4. Interpretation rule (fixed — deliberately asymmetric)

This is the clause that keeps Route E honest, and it is written before any number
exists.

**A positive is informative. A negative is not.**

- **`D_tumour` clearly above `D_blood` and comparable to `D_NILT`** → consistent
  with a locally expanded, resident-like pool. Reported as **a hint consistent
  with X1, not as evidence of ancestry.**
- **`D_tumour ≈ D_blood`** → **reported as uninformative, NOT as evidence against
  X1.** A tumour can recruit on top of a resident pool and flatten its clonality;
  a resident pool can be diluted below detection by later infiltration; and clone
  assignment in solid tissue runs ~15–26%, so the tumour compartment is the one
  most prone to losing structure to missing assignment.

**No p-value is reported for a cross-compartment difference in `D`.** The
comparison is descriptive, on 4 donors, on a quantity whose null distribution is
not established for this data. Percentile intervals from the resampling are
reported as sampling variability **only**, and explicitly not as inference.

**Nothing in Route E can address the temporal clause of X1.** It cannot
distinguish a pool resident before the tumour from one that became resident after
it. Route E speaks to clause (a) weakly and clause (b) not at all.

---

## 5. Controls (fixed)

The same estimator is computed on **myeloid** and **T** cells in the same donors,
same downsampling. These calibrate the estimator on lineages whose behaviour the
source paper characterised — myeloid is reported there as clonally linked to
circulating monocytes, so myeloid is expected toward the blood-like end.

**If NK is not distinguishable from the myeloid control, the NK result is not
reported as a finding.** The controls are part of the primary output, not a
supplementary check.

---

## 6. Known limits carried in

1. **NK typing is cluster-level and validated only on blood** (CD56⁺ 74.6% vs
   CD56⁻ 5.5%; DEC-08). In tissue, where the ILC1/trNK boundary is exactly what
   this claim is about, it is a power estimate rather than a validated classifier.
   Route E inherits that uncertainty in full, and it bears directly on the tumour
   compartment where the numbers are smallest.
2. **Clone assignment ~15–26% in solid tissue**, higher in blood — so the tumour
   compartment loses more structure to non-assignment than blood does, which
   biases `D_tumour` **downward** relative to `D_blood`. That bias runs *against*
   the resident hypothesis, which is why a positive is informative and a negative
   is not (§4).
3. **n = 4 donors** at most.
4. This is a re-analysis of one dataset by one person; it replicates nothing.
