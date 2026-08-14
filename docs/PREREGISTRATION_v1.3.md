# PREREGISTRATION v1.3 — the FDR criterion is withdrawn, not adjusted

**v1.0 (`PREREGISTRATION.md`) remains untouched and read-only** — one commit in
its history. v1.1 and v1.2 stand as written.

**Timestamp:** 2026-08-13, after the v1.2 run.
**Authorised by:** the project owner.

**Direction of effect: this amendment removes a decision rule and adds none.**
Nothing is promoted. The number of "significant" rows becomes undefined
because the criterion that produced it is withdrawn.

---

## A8 — No FDR is reported. The pre-registered criterion is not estimable at this n.

**v1.0/v1.1/v1.2 rule.** BH over the empirical p, family of 81 (v1.2 A5).

**Why it fails, arithmetically.** The empirical p is rank-based, so with N
null genes per target its floor is 1/(N+1). BH needs the smallest p to reach
`0.05 × k / 81` for k rows tied at that floor:

| N per target | rank floor | rows needed at the floor before any can pass |
|---|---|---|
| 112 | 0.00885 | **15** |
| 317 | 0.00314 | **6** |
| 371 | 0.00269 | **5** |
| 399 | 0.00251 | **5** |

There are only 8 candidate rows in the whole study. **So the number of rows
that pass is determined by the null-set size, not by the effects.** Measured:
the same data gave 6 passes, then 0, then 8, purely as N changed.

That is not a rule needing a tuned parameter. It is a rule that does not hold
at this n, in the same way and for the same reason as ruler B's fixed
denominator (v1.2 A7) and the signed-rank BH (v1.1 A2).

**v1.3 rule.** No FDR, no q value, no pass/fail list. Every affected row is
emitted with `fdr_estimable = FALSE` and the reason. Reported instead, as
descriptive quantities:

- effect in **percentage points**
- **z against the matched null**, as a standardised effect size, explicitly
  **never converted to a p**
- **donor sign concordance** (all donors moving the same way)
- **ruler-A margin** (soup fraction against the 0.092 ceiling and 0.021 floor)
- the observed effect's **rank** within its null set, and the null-set size

**No replacement decision rule is invented.** Readers weigh the descriptive
quantities. A study that cannot support a decision rule should say so rather
than substitute a weaker one and keep the appearance of a threshold.

---

## A9 — Null genes are matched on baseline detection rate, not CPM

**Defect.** The readout is a detection rate, bounded in [0,1], whose sampling
variance is p(1−p)/n — set by the rate, not by CPM. Matching nulls on CPM
gives a low-detection target a null set whose variance is compressed against
the zero floor for reasons unrelated to the target. That is what let a
**1.5-percentage-point S1PR5 effect** reach the top of the v1.2 list.

**v1.3 rule.** Each target's null set is matched on the **reference arm's
baseline detection rate**, computed over the whole transcriptome at matched
depth. Since the reference arm differs per contrast, each contrast builds its
own null sets.

**Direction: intended to be stricter. It is not, and that is a finding.**

Measured after the change:

| baseline detection of reference arm | n | median null SD (pp) | median abs effect (pp) |
|---|---|---|---|
| < 0.5% | 29 | 0.33 | 0.17 |
| 0.5–2% | 11 | 0.82 | 0.48 |
| 2–10% | 8 | 2.06 | 2.07 |
| 10–30% | 14 | 5.31 | 3.14 |
| > 30% | 17 | 6.41 | 10.15 |

**A z of 5 means a 32.1-point effect for a gene detected in >30% of cells and
a 1.63-point effect for one detected in <0.5% — a 20-fold difference in what
the same z buys.** Matching on detection rate did not fix that; it made it
sharper, because the nulls for a low-detection target are now themselves
low-detection and their spread is smaller still. S1PR5 went from z = +5.01 to
z = **+10.34** on the same +1.5 pp effect; SELL entered at z = +10.55 on
+4.8 pp.

**Conclusion, recorded as a limit rather than repaired.** The problem is not
the matching axis. **z is not comparable across genes on a bounded scale near
its boundary**, however the nulls are chosen. Therefore:

> **z may never be read without the effect in percentage points beside it.**
> Both are emitted in `null_descriptive.tsv`, and no ranking by z alone
> appears in `RESULTS.md`.

---

## What this leaves standing

The two surviving claims of the study depend on **none** of the machinery
withdrawn in v1.1–v1.3 — no p, no q, no FDR, no empirical null:

1. **XCL1/XCL2 are higher in dNK2 than dNK1**: ≥ +39.3 / +32.9 pp, 6/6 donors
   concordant, survives marker-based regrouping and the within-label test
   (+6.7 / +5.7 pp, 3/3 donors), ruler A 0.021 against a 0.092 ceiling, and
   carries the opposite-sign internal control (XCL1 +6.7 while CXCR4 −3.5 in
   the same contrast on the same cells).
2. **CXCR4's subset gradient reverses within the dNK1 cluster** (−3.5 pp,
   3/3 donors) and does not hold.

That is the methodological result of this project: **what survived is exactly
what never needed the null distribution.**

---

## Unchanged

R1–R11 as amended; the frozen 68-gene panel; the 30-cell threshold; the 15%
purity stop; the n-floor clause; `unmeasurable` ≠ `negative`; v1.2 A4's narrow
re-grouping exception.

## Known and not repaired in this round

The ruler-A ceiling is a **minimum over ambient controls**, so it is set by
whichever control the estimator handles worst, and it falls monotonically as
controls are added. Measured class structure (D23): haematopoietic controls
return a median soup fraction of 0.250 and 1.8% saturation, non-haematopoietic
1.000 and 90% — so the ceiling of 0.092 is set entirely by the haematopoietic
class, which is the class the estimator handles worst because NK cells are
themselves haematopoietic. A lower-prediction-bound construction within
lineage class would be more stable. Deferred.
