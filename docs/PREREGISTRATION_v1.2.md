# PREREGISTRATION v1.2 — amendments to v1.1

**v1.0 (`PREREGISTRATION.md`) remains untouched and read-only** — still one
commit in its history. v1.1 remains as written. This file records four
amendments plus one authorised exception, each with its direction of effect.

**Timestamp:** 2026-08-13, after the v1.1 run.
**Authorised by:** the project owner, in response to the v1.1 report.

**Disclosure of direction.** Unlike v1.1, **three of these four amendments make
the results weaker or narrower** (A5, A6, A7), and one is neutral in verdict
but corrects a precision claim (A6). None was chosen to promote a gene.

---

## A4 — Authorised exception: marker-gated re-grouping, for a circularity test only

v1.0 §7 bars re-clustering and re-labelling. That bar stands for the main
analysis. **One exception is authorised**, for a test the bar itself makes
impossible to run otherwise.

**Why it is needed.** dNK1/dNK2/dNK3 come from Vento-Tormo's own **unsupervised
whole-transcriptome clustering**, with the subsets then characterised by
differential expression (verified against the published methods). CCL5, CXCR4
and XCL1 are among the genes that characterisation named. Measuring, in the
same data under the same labels, that dNK3 is CCL5/CXCR4-high therefore partly
**restates how the clusters were drawn** — not because anyone thresholded on
CCL5, but because the cluster boundary came from a transcriptome that includes
it.

**The exception.** `src/circularity_test.py` may re-group decidual NK cells
using **only** surface-marker genes mapping onto the Huhn CyTOF gating:

| gate | definition | ≈ |
|---|---|---|
| gate1 | ENTPD1⁺ (CD39⁺) | dNK1 |
| gate3 | ENTPD1⁻ ITGAE⁺ (CD39⁻CD103⁺) | dNK3 |
| gate2 | ENTPD1⁻ ITGAE⁻ | dNK2 |

`ITGB2` (CD18) is reported but not gated on. **`CD160` and `KLRB1` are barred
as gating genes** — they are themselves dNK3 characterisation genes, so gating
on them would rebuild the circularity the test exists to break. No target gene
may inform the gate.

**Constraints.** The published annotation is not modified. Gate assignments
live in their own tables. **No output of this script feeds the main analysis.**

---

## A5 — The BH family is every analysis-A test at once

**v1.1 rule.** BH within each pairwise comparison, 27 genes per family.

**Defect.** v1.0 said "BH across the panel" without saying whether the three
pairwise contrasts form one family or three. They are not three independent
hypothesis families — they are three views of the same cells — so correcting
within each understates multiplicity by up to 3×.

**v1.2 rule.** One family: **27 genes × 3 contrasts = 81 tests**, corrected
together (`q_bh_empirical`). The per-comparison version is retained as
`q_bh_empirical_per_comparison` so the difference is visible, not assumed away.

**Direction: stricter.** Measured effect: all 8 rows still clear q ≤ 0.05, but
two of them land at **q = 0.0499** — inside the threshold by 0.0001. That
margin must be reported as the knife-edge it is.

---

## A6 — The empirical p may not be reported below its resolution

**v1.1 defect.** Rows were reported at `q < 0.0001`. A rank-based empirical p
over N null genes cannot resolve below **1/(N+1)**; at N = 405 that is
**0.00246**. Anything smaller could only come from converting z through a
normal tail — and a normal tail at z = +14.6 is exactly what 405 resamples
cannot verify.

**v1.2 rule.** The empirical p uses the standard +1 correction,
`(n_ge + 1)/(N + 1)`, is floored at its resolution, and each row carries
`empirical_p_resolution_floor` and `empirical_p_at_resolution_floor`. A row at
the floor is reported as **"≤ 0.0025, at the resolution limit"**, never as a
smaller number. The z score is still reported, as a standardised effect size,
**not** converted to a p value.

**Direction: neutral in verdict, stricter in claim.** No row changes
pass/fail; six rows stop claiming precision the resampling cannot support.

---

## A7 — Ruler B's denominator is the largest TOTAL contributor, and NK dominance means `unmeasurable`

**v1.1 defect.** The denominator was the lineage with the highest **per-cell
CPM**. The ambient pool's composition is set by **total transcript output —
CPM × cells × depth** — so ranking by CPM lets a small, highly-expressing
population hijack the denominator. Measured: v1.1 assigned XCL1's source to
**ILC3, a 184-cell population** that cannot plausibly dominate any ambient
pool.

**v1.2 rule.** The denominator is the lineage contributing the most **total
counts** of that gene. And:

> If **NK itself** is the largest contributor of a gene to the pool, there is
> no competing source, the ambient argument does not apply, and ruler B
> records **`unmeasurable`** — never `positive`.

This is the same principle already in T3: `unmeasurable` ≠ `negative`. Two
reasons are now distinguished in the supplement
(`scale_b_unmeasurable_reason`):

- `no_competing_source_NK_is_largest_contributor` — **supports** genuine
  expression; ruler A carries the call
- `source_lineage_does_not_express_it_enough` — no power in either direction

**Direction: strictly weaker on paper.** Under v1.2, `XCL1`, `XCL2`, `CCL5`,
`CXCR4` and `CCL4` all move from ruler-B `above_pickup_band` to
`unmeasurable`, so their `combined_verdict` moves from `positive` to
`unmeasurable`. **The "both rulers positive" criterion can no longer be met by
any of the headline genes.** Their soup evidence now rests on ruler A alone
(soup fraction 0.021–0.048 against a 0.092 ceiling and a 0.021 true-NK floor)
plus the reason code above.

Band recalibrated under v1.2 (decidua): **0.0538**, with every ambient control
against a sensible source — LYZ→Myeloid 0.054, C1QA→Myeloid 0.044,
IGFBP1→Stromal 0.036, COL1A1→Stromal 0.016, DCN→Stromal 0.011,
HLA-G→Trophoblast 0.002, IGKC→Plasma 0.002, CSH1→Trophoblast 0.0005.

---

## Reading the headline genes under v1.2

The pre-registered "positive requires both rulers" rule is retained unchanged,
and under v1.2 **no headline gene satisfies it**. That is not a downgrade of
the evidence but a change in which instrument carries it:

| instrument | XCL1 / XCL2 / CCL5 / CXCR4 |
|---|---|
| ruler A | `above_ambient`, 0.021–0.048 vs ceiling 0.092 |
| ruler B | `unmeasurable` — NK is the largest contributor |
| empirical null | q ≤ 0.05 in the 81-test family |
| circularity (A4) | see `RESULTS.md` §3 — it separates them |

A gene whose largest pool contributor is NK itself cannot be argued to be
ambient pickup; but neither can ruler B be used to argue it is real. **The
verdict column reports `unmeasurable` and the reason code carries the
interpretation.**

---

## Unchanged

R1–R11 as amended by v1.1; the frozen 68-gene panel; the 30-cell threshold;
the 15% purity stop; the n-floor clause; the `unmeasurable` ≠ `negative` rule;
everything in v1.0 §7 except as narrowly authorised in A4.
