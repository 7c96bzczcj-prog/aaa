# PREREGISTRATION v1.1 — amendments to v1.0

**v1.0 (`PREREGISTRATION.md`) is unchanged and remains read-only.** It was
frozen at commit `175713a` before any expression value was read, and it has
exactly one commit in its history. This file records three amendments, the
reason for each, and — importantly — **which direction each one moves the
results**, so that a reader can judge whether an amendment was made to rescue
an unwelcome outcome.

**Timestamp:** 2026-08-13, after the v1.0 run completed.
**Authorised by:** the project owner, in response to the v1.0 report.

**Disclosure of direction.** Amendments A1 and A3 both *increase* the number
of positive calls. That is the pattern of a rule loosened to manufacture
results, so each one below states the pre-existing methodological defect it
repairs, independent of any gene. Amendment A2 does not loosen a threshold at
all — it replaces a statistic that provably has no power at this n with one
that does. None of the three was chosen after inspecting which genes it would
promote; A1 and A3 were both flagged as defects in the v1.0 report itself,
before the amended results existed.

---

## A1 — R5 (ceiling/floor) becomes a per-arm rule

**v1.0 rule.** A contrast is `uninterpretable_ceiling` if **either** arm has
detection rate > 0.85.

**Defect.** The rule was written from a case where a *null* was misread at
93.75% detection — both arms saturated, no room to move, and "no difference"
was wrongly concluded. It does not follow that a **one-sided** ceiling is
uninterpretable. When only the high arm saturates, saturation compresses the
observed difference **toward zero**. A large positive difference measured
against a saturated high arm is therefore a **lower bound** on the true
difference, not an unreadable number. v1.0 discarded exactly the contrasts
where the effect was largest, which is the opposite of conservative.

**v1.1 rule.**

| condition | verdict |
|---|---|
| both arms > 0.85 | `uninterpretable_ceiling` — unchanged, still excluded |
| one arm > 0.85 **and** the effect points toward that arm | direction retained; `magnitude_is_lower_bound = TRUE` |
| one arm > 0.85 and the effect points away from it | `uninterpretable_ceiling` (implies both arms are high) |

The symmetric case at the floor is treated the same way
(`magnitude_is_lower_bound` for a one-sided floor).

**Direction of effect on results: increases positives.** Under v1.0, 86
contrasts were excluded as ceiling-limited; under v1.1, 5 are excluded (both
arms saturated) and 81 are retained as direction-only with the magnitude
marked a lower bound. This is the amendment most in need of scrutiny, which
is why the new columns `ceiling_arm_a`, `ceiling_arm_b` and
`magnitude_is_lower_bound` are emitted per row: every reclassified contrast
can be audited individually.

**Not changed.** The 0.85 and 0.05 thresholds themselves, and the rule that a
magnitude may never be quoted from a saturated arm.

---

## A2 — BH over the signed-rank p is reported, not adjudicating. The empirical null adjudicates.

**v1.0 rule.** BH across the 27-gene primary-target family, on the exact
Wilcoxon signed-rank p.

**Defect, established before the data.** v1.0 §2.2 already computed that at
n = 6 a single gene needs ≥ 17 of 27 genes simultaneously at the p = 0.03125
floor for q ≤ 0.05, and at n = 5 the requirement (≥ 34 of 27) is impossible.
**BH over a floor-limited p can never reject at these donor counts, whatever
the biology.** A test statistic that cannot reject is not a conservative
test; it is an absent one.

**v1.1 rule.** The primary inference is the **empirical null** — the
distribution of the same donor-level effect over ≥ 400 expression-matched
random genes, run through the identical pipeline (v1.0 §6.5, unchanged).
Multiplicity is controlled by **BH over the empirical p** within the same
27-gene family (`q_bh_empirical` in T4). `q_bh` over the signed-rank p is
retained in T2 and reported truthfully, but does not decide anything.

**This is not a loosened threshold.** α stays at 0.05 and the family stays at
27 genes. What changes is the statistic the correction is applied to: from
one whose resolution is 2/2ⁿ (0.031 at best here) to one whose resolution is
1/n_null (0.0025 at 405 genes). The empirical p is also the stricter test in
a substantive sense — it asks whether an effect exceeds what *matched random
genes* produce in this same dataset, with the baseline measured rather than
assumed (measured: +0.47 ± 4.05 pp, not zero).

**Direction of effect on results: increases positives**, from 0 rows at
q ≤ 0.05 to 8. Every one of those 8 sits at |z| ≥ 5.6 against the matched
null, so none is a marginal call.

**Unchanged:** R1 (n is donors), R2 (`min_achievable_p` and `on_test_floor`
travel with every row, including the empirical rows).

---

## A3 — Ruler B is per-gene, against that gene's own dominant source lineage

**v1.0 rule.** Ruler B = NK / **myeloid** CPM ratio, one band for all genes.

**Defect, measured in the v1.0 run.** Ruler B asks whether NK-gate signal is
explainable as pickup from a competing lineage. A **fixed** denominator only
tests that when the fixed lineage is the gene's actual ambient source. For a
stroma- or trophoblast-sourced gene, the NK gate and the myeloid gate pick up
the ambient pool **equally**, so the ratio sits near 1 regardless of the
truth. Measured in v1.0: the band calibrated over all ambient controls ran to
**1.253**, versus **0.054** over myeloid-sourced controls only — a 23-fold
difference that is the concept error made visible, not noise.

**v1.1 rule.** For each gene, the denominator is the **lineage with the
highest measured CPM for that gene** in that compartment, chosen from the
manifest's mapped lineages (≥ 30 cells, ≥ 10 CPM, ≥ 5% detection). The band
is calibrated the same way: each ambient control is judged against its own
dominant source. This is data-driven and generic — no dataset branch, no
hand-assigned gene→lineage table.

Choosing the *highest-expressing* lineage is the conservative choice: it
maximises the denominator, minimises the ratio, and therefore makes a gene
harder to call positive.

Measured under v1.1 in decidua, every ambient control now lands in a low
band against a sensible source, which is what "pure pickup" should look like:

| control | dominant source | source CPM | NK/source |
|---|---|---|---|
| CSH1 | Trophoblast | 13,408 | 0.0005 |
| HLA-G | Trophoblast | 4,877 | 0.0020 |
| IGKC | Plasma | 10,428 | 0.0019 |
| COL1A1 | Perivascular | 2,080 | 0.0030 |
| DCN | Stromal | 5,535 | 0.0106 |
| HBB | Epithelial | 138 | 0.0123 |
| IGFBP1 | Stromal | 2,500 | 0.0361 |
| LYZ | cDC1 | 2,124 | 0.0428 |
| C1QA | Myeloid | 3,962 | 0.0441 |

**Band = 0.0441**, against v1.0's 1.253 and prior work's 0.118.

**Direction of effect on results: increases positives** (a narrower band is
easier to clear). The v1.0 verdict is retained per gene as
`scale_b_verdict_v10_myeloid_only` in the supplement, so both readings are on
the record; the frozen `nk_myeloid_cpm_ratio` column is also unchanged.

`unmeasurable` still means what it meant: the source lineage does not express
the gene enough for the ratio to have power. It never means `negative`.

---

## Unchanged from v1.0

R1–R4, R6–R11; the 68-gene frozen panel; the 30-cell admission threshold; the
15% purity stop; the three-outcome structure of §3 and the n-floor clause;
everything in §7 (barred practices). No result was reclassified except through
the three amendments above, and no threshold in v1.0 was moved.
