# PREREGISTRATION v1.4 — the effect scale is fixed, before any further result

**v1.0 (`PREREGISTRATION.md`) remains untouched and read-only.** v1.1–v1.3
stand as written.

**Timestamp:** 2026-08-13, after the v1.3 run.
**Authorised by:** the project owner.

**Direction of effect: neither. This fixes a scale so that later results
cannot be read on whichever scale flatters them.** It is committed *before*
any further result is generated, which is the only time such a choice is
worth anything.

---

## A10 — Percentage points are the primary scale. z is auxiliary and never read alone.

**Why it must be settled now.** v1.3 A9 established that z is not comparable
across genes on a bounded scale near its boundary: a z of 5 buys a 32.1-point
effect at >30% baseline detection and 1.63 points at <0.5%. With two scales
available and no rule, every future result could be reported on whichever one
favours it.

**Why no transformation fixes it.** This is not a statistical defect to be
repaired by a link function — it is a choice about which scale answers the
question. The question is *"does a larger fraction of cells carry this
chemokine transcript?"*, and a fraction difference is naturally in points.
A logit transform would move the reading **further** toward low-baseline
genes, since logit expands exactly the region where detection is compressed
against zero. Changing scale changes which genes look important; it does not
remove the arbitrariness.

**v1.4 rule.**

1. **The primary effect scale is percentage points** of detection rate.
2. **z is auxiliary.** It is a standardised effect against a matched null and
   is (a) never converted to a p, and (b) **never reported, ranked or
   discussed without the percentage-point effect beside it.**
3. No table in `RESULTS.md` is sorted by z.

---

## A11 — The measurable band, and what falls outside it

The design detects subset differences only for genes already detected in a
workable fraction of cells. Measured across the 27 primary targets, pooled
over decidual dNK1–3 at the primary depth floor:

| band | n | genes |
|---|---|---|
| **floor, < 5% — no power in points** | **14** | PPBP 0.000, CCL1 0.001, ACKR2 0.002, CCR2 0.002, CXCL12 0.003, CCR9 0.004, S1PR1 0.004, S1PR5 0.005, CX3CR1 0.005, CXCL10 0.006, CCR5 0.008, CXCR6 0.012, CXCL16 0.024, SELL 0.033 |
| measurable, 5–85% | 12 | CCL20 0.052, CCL2 0.054, CCR1 0.064, CXCR3 0.177, CXCL8 0.235, CCL3L1 0.256, CCL4L2 0.383, CXCR4 0.421, CCL5 0.585, XCL1 0.756, XCL2 0.801, CCL3 0.817 |
| ceiling, > 85% | 1 | CCL4 0.944 |

**More than half the primary panel — 14 of 27 — has no power on the primary
scale at 10x depth in this dataset.** That is a property of the design, and it
is stated here rather than discovered later.

**The double failure closes on the same genes.** Four of those 14 floor-band
genes — **S1PR5, SELL, CCR9, CXCR6** — are exactly the ones that produced the
largest and most misleading z values in v1.3 (S1PR5 z = +10.34 on +1.5 pp;
SELL z = +10.55 on +4.8 pp). So the floor band has **no power in points and
false power in z simultaneously**, in the same genes. That is the strongest
possible argument for A10, and it is why no effect-size filter was added
after the fact: the filter is the scale rule, declared in advance.

### Correction carried into this amendment

The reading "CCL3 does not reproduce, which may mean no power rather than no
difference" is **not supported and is not adopted**. CCL3 sits at **0.817**,
inside the measurable band; its flat, sign-inconsistent result is a finding
with power behind it, not an absence of power. The genuinely
power-limited genes are the 14 listed above, and CCL3 is not among them.
CCL4 (0.944) is limited by the **ceiling**, which is a different failure with
a different remedy (greater depth does not help; it makes it worse).

---

## A12 — Consequence for reading a null result

A null result in this study means one of three things, and which one must be
stated with it:

| detection band | a flat result means |
|---|---|
| < 5% | **untested.** Report as no power, never as no difference |
| 5–85% | **tested and flat.** A real null, within the CI reported |
| > 85% | **uninterpretable.** Saturated; the contrast has no room (R5) |

---

## Unchanged

R1–R11 as amended by v1.1–v1.3; the frozen 68-gene panel; the 30-cell
threshold; the 15% purity stop; `unmeasurable` ≠ `negative`; no FDR (v1.3 A8);
v1.2 A4's narrow re-grouping exception.
