# RESULTS — dNK chemokine profile, VT2018

Dataset: `VT2018` (E-MTAB-6701, Vento-Tormo et al., Nature 2018).
Panel: `chemokine_panel_v1` (68 genes, 98.5% matched).
Pre-registration: [`PREREGISTRATION.md`](PREREGISTRATION.md), frozen at commit
`175713a` **before** any expression value was read.

Primary depth floor **2137 UMI** (10th percentile of NK-cell depth);
sensitivity floors 1808 and 2784. Every rate below is a detection rate at
matched depth (R3). Every test consumes one number per donor (R1).

---

## 0. The two things that bound everything else

**The pre-registered stop condition fired**, and the pre-registered power
ceiling held exactly as written. Both were established before the results
were looked at, and neither can be argued away afterwards.

| | pre-registered expectation | measured |
|---|---|---|
| donors, analysis A | ≤ 6 | 6 (dNK1/dNK2), 5 (dNK3 contrasts) |
| min achievable p, dNK1 vs dNK2 | 0.03125 | 0.03125 |
| min achievable p, dNK3 contrasts | 0.0625 — **p < 0.05 unreachable** | 0.0625 |
| donors, analysis B vs CD56-dim | 4 — **p < 0.05 unreachable** | 4 |
| analysis B vs CD56-bright | `unmeasurable` at n = 2 | n = 2, min p = 0.500 |
| genes needed at the floor for BH q ≤ 0.05 | ≥ 17 of 27 | 6 of 27 reached it |

**No gene in the study reaches q ≤ 0.05, and the pre-registration says why
before the fact: at n = 6 a single gene cannot, however large its effect.**
Everything below is therefore reported as effect size, sign concordance,
empirical-null z, and dual-ruler verdict — not as a star.

---

## 1. Conclusions

Each conclusion carries: evidence rows, n donors, effect in percentage
points, both rulers, and what it does **not** support.

---

### C1. dNK2 cells carry CCL5 transcript at a higher rate than dNK1 from the same donor — in all six donors.

- **Evidence:** `out/VT2018/donor_level_tests.tsv` line 5;
  `soup_calibration.tsv` line 72; `null_distribution.tsv` line 5.
- **n = 6 donors** (D6, D7, D8, D9, D10, D12).
- **Effect: +23.3 percentage points** (pooled detection dNK1 33.5% →
  dNK2 56.8%), 95% bootstrap CI [17.1, 30.1]. p_raw = 0.03125 = the exact test floor
  (`on_test_floor = TRUE`, i.e. all six donors moved the same way),
  q_BH = 0.141.
- **Ruler A:** soup fraction 0.025, against an ambient ceiling of 0.092 and a
  true-NK floor of 0.021 → **above ambient**.
  **Ruler B:** NK/myeloid CPM ratio 6.26, band upper edge 1.25 →
  **above pickup band**. `combined_verdict = positive`.
- **Empirical null:** +23.3 pp against a matched-expression baseline of
  +0.47 ± 4.05 pp over 405 genes → **z = +5.6, empirical p = 0.0025** (the
  floor at 405 null genes).
- **Stable across all three depth floors** (+20.0 / +23.3 / +21.3 pp,
  p = 0.03125 at each) and across the T-cell-purged arm (+22.5 pp).
- Neither group is at the ceiling or the floor — this is the **only**
  primary-target contrast in the study that clears every pre-registered
  filter simultaneously.

**This does not support:** that dNK2 secretes more CCL5 protein; that CCL5
protein reaches any receiving cell; that this difference matters for any
process in the decidua. It is a transcript-detection difference measured at
one depth in one cohort of six donors, and it is not significant after
multiplicity correction.

---

### C2. CCL5 rises monotonically dNK1 < dNK2 < dNK3, but only the dNK1-vs-dNK2 step is interpretable.

- **Evidence:** `donor_level_tests.tsv` lines 5, 43, 81;
  `null_distribution.tsv` lines 5, 43, 81.
- Pooled detection: dNK1 33.5%, dNK2 56.8%, **dNK3 90.4%**.
- dNK1 vs dNK3 **+57.3 pp** (n = 5, z = +10.9); dNK2 vs dNK3 **+35.5 pp**
  (n = 5, z = +9.8). Both sit on their test floor (all donors concordant).
- **Both dNK3 contrasts are flagged `uninterpretable_ceiling`** — dNK3 is at
  90.4% detection, above the 0.85 bar, at every depth floor tested
  (0.907 / 0.904 / 0.900). R5 removes them from the conclusions.

**This does not support:** a magnitude for the dNK3 difference. At 90%
detection the measurement is compressed and the effect is a lower bound
only. It also does **not** license the opposite error — the historical one
this rule exists to prevent — of reading a saturated rate as "no difference".

---

### C3. XCL1 and XCL2 show the largest subset differences in the panel, in the published direction, and both are ceiling-limited.

- **Evidence:** `donor_level_tests.tsv` lines 7, 8, 45, 46;
  `soup_calibration.tsv` lines 74, 75; `null_distribution.tsv` lines 7, 45.
- dNK1 vs dNK2: XCL1 **+39.3 pp** (n = 6, p = 0.03125 on floor, z = +9.6,
  empirical p = 0.000); XCL2 **+32.9 pp** (z = +8.0).
  dNK1 vs dNK3: XCL1 +31.2 pp, XCL2 +23.0 pp (n = 5, both on floor).
- **Ruler A** 0.021 (identical to the true-NK gene floor);
  **ruler B** 13.5 and 12.8 versus a band edge of 1.25 → `positive` on both.
  Of every gene in the panel these sit furthest from ambient.
- **Flagged `uninterpretable_ceiling`**: dNK2 reaches 91.1% and dNK3 85.7%
  detection. At the lowest depth floor dNK2 is still 89.7%, so this is not a
  depth artefact that a different floor would dodge.
- dNK2 vs dNK3 is flat (−5.0 pp, z = −1.4): the XCL1 difference is
  **dNK1 versus the other two**, not a graded trend.

**This does not support:** any statement about how much more XCL1 dNK2 makes.
Above 85% detection the readout is saturated and only the direction survives.

---

### C4. CCL3 and CCL4 do not reproduce a subset difference. CCL4 cannot be tested at this depth at all.

- **Evidence:** `donor_level_tests.tsv` lines 3, 4, 41, 42, 79, 80.
- CCL4 detection is 92.9–96.2% across all four dNK subsets —
  `uninterpretable_ceiling` in every contrast. Effects are +0.5 to +2.9 pp
  with p = 0.44–0.81.
- CCL3 is 77.7–88.4%; the largest effect is +10.1 pp (dNK2 vs dNK3,
  p = 0.19, z = +2.8), and its sign is inconsistent across contrasts
  (−2.3 pp for dNK1 vs dNK2).
- CCL3 ruler verdict `weak` (ruler A 0.043 above ambient, ruler B 0.39 inside
  the 1.25 band); CCL4 `positive` on both rulers but uninterpretable on
  detection.

**This does not support:** "dNK subsets do not differ in CCL3/CCL4". For CCL4
this dataset has no resolving power — near-universal detection means the
question is unanswerable here, which is the pre-registered
**A-uninterpretable** outcome for that gene, not a negative result.

---

### C5. PTN and OGN in the NK gate are entirely stromal pickup.

- **Evidence:** `soup_calibration.tsv` lines 126, 127;
  `soup_calibration_supplement.tsv` (same genes).
- **Ruler A soup fraction = 1.000 for both** — the value a *pure ambient*
  gene returns.
- Forced parallel stromal readout: NK/stromal CPM ratio **0.007 (PTN)** and
  **0.016 (OGN)**; stromal CPM 118.3 and 90.1 against NK CPM 0.86 and 1.43.
- Ruler B returns `unmeasurable`, correctly: myeloid is not these genes'
  source, so a NK/myeloid ratio has no power over them. This is exactly the
  case the `unmeasurable`-is-not-`negative` rule exists for.

**This does not support:** any claim about PTN or OGN in NK cells in either
direction. It says the observed signal is indistinguishable from
stroma-derived ambient RNA, which is what the panel flagged these two genes
for in advance.

---

### C6. The XCR1 negative control confirms a low, compartment-comparable ambient background.

- **Evidence:** `detection_rates.tsv`, gene `XCR1`.
- Decidual dNK: 0.20% (dNK1), 0.22% (dNK2), 0.53% (dNK3), 0.00% (dNKp).
  Blood NK: 0.00%.
- Decidua−blood difference **0.22 pp < the pre-registered 2 pp gate**, so
  analysis B is formally readable (`analysis_b_gate.json`).

**This does not support:** the ambient background being negligible for every
gene. XCR1 is a low-abundance transcript; abundant tissue genes (PAEP-class,
DCN, IGFBP1) are picked up at soup fraction 1.000 in the same cells.

---

### C7. Analysis B is descriptive only, exactly as pre-registered.

- **Evidence:** `donor_level_tests.tsv`, all `B_*` rows (n = 4 against
  CD56-dim, n = 2 against CD56-bright).
- **Every B row carries `cross_compartment_soup_caveat = TRUE`** and a
  minimum achievable p of 0.125 (or 0.500 at n = 2). **No B row can be
  significant.**
- Direction, for the record only (`mean_diff_pp` is blood minus decidua, so a
  positive number means blood NK are higher). Against CD56-dim, blood NK
  exceed dNK on `SELL` (+25.6 to +30.6 pp), `S1PR5` (+45.9 to +47.9) and
  `CX3CR1` (+20.9 to +21.2) — the egress/recirculation set. dNK exceed blood
  NK on `XCL1` (−36.9 pp for dNK1, −74.3 for dNK2, −70.4 for dNK3). `CCL5`
  changes sign by subset (+38.5 pp for dNK1, +15.8 for dNK2, −12.5 for dNK3),
  which is the same dNK1 < dNK2 < dNK3 gradient as C2 seen from the blood
  side. Many B rows are additionally ceiling-flagged.
- Against **CD56-bright**, the biologically matched comparator, only 2 donors
  clear 30 cells. Recorded as **`unmeasurable`**, not as a null.

**This does not support:** any decidua-versus-blood conclusion. The decidual
soup is stroma/trophoblast-rich and the blood soup is platelet/erythroid-rich,
so the two backgrounds differ in composition regardless of the XCR1 gate.

---

## 2. Verdict against the pre-registered outcomes

The pre-registration allows A-positive, A-negative, A-uninterpretable, and
A-underpowered. **The honest answer is per-gene, not one label:**

| gene | outcome | reason |
|---|---|---|
| **CCL5** (dNK1 vs dNK2) | **A-positive in substance, A-underpowered formally** | direction matches the literature, both rulers positive, all 6 donors concordant, z = +5.6; but q_BH = 0.141, and the pre-registration established beforehand that a single gene cannot reach q ≤ 0.05 at n = 6 |
| **XCL1 / XCL2** | **A-positive in direction, A-uninterpretable in magnitude** | largest effects in the panel and furthest from ambient, but at the detection ceiling at every depth floor |
| **CCL5** (dNK3 contrasts) | **A-uninterpretable** | dNK3 saturated at 90% |
| **CCL4** | **A-uninterpretable** | 93–96% detection everywhere |
| **CCL3** | **A-negative** | small, sign-inconsistent, ruler B inside the band |

Read together: **the subset-specific chemokine pattern reported from CyTOF
does appear at transcript level for XCL1/XCL2 and CCL5, with unanimous
donor-level direction and effects far outside an expression-matched null —
but this dataset cannot certify it at conventional significance, and cannot
measure its size for the genes where the effect is biggest.**

That is neither a confirmation nor a refutation of the published protein
work. It is a statement about what a 6-donor 10x cohort can carry.

---

## 3. Caveats that would change these readings

1. **The purity stop fired.** Raw triple-positive `TRBC2⁺CD3E⁺CD3D⁺` in dNK1
   is **17.4%**, above the 15% stop. Diagnosis in
   [`DECISIONS.md` D10](DECISIONS.md): it is driven by sequencing depth and
   by NK–T doublets, not by mislabelled T cells, and every conclusion above
   is reproduced in the T-purged arm
   (`out/VT2018/robustness_tpurged/`). C1 moves from +23.3 to +22.5 pp.
2. **Differential cell loss is severe (R6).** At the primary floor, donor D8
   retains 85.9% of dNK1 but only **42.8%** of dNK3; D12 retains 56.0% vs
   33.4%. Three donors exceed the 10 pp trigger. The mandated three-floor
   sensitivity is in `donor_level_tests_depth_sensitivity.tsv`; the dNK1-vs-
   dNK2 findings are stable across floors, the dNK3 findings weaken at the
   deepest floor as D8's dNK3 falls away.
3. **Ruler B has little power in decidua.** Calibrated on all ambient
   controls the pickup band runs to 1.25, because stroma- and
   trophoblast-sourced genes are picked up about equally by NK and by
   myeloid cells. Restricted to controls myeloid actually sources
   (LYZ, C1QA) the band is 0.054 — a **23-fold** difference. The
   conservative (wider) band was used for all verdicts, so ruler B is if
   anything under-calling positives.
4. **PAEP is absent from this matrix entirely** (see D11). The specification
   named it the decidua-specific ambient ceiling gene; the published matrix
   does not contain it, so the ceiling rests on DCN, COL1A1, IGFBP1, HLA-G
   and CSH1 — all of which do return soup fraction 1.000.
5. **`CCL3L1` matched only via Ensembl** and its symbol maps to three
   Ensembl gene IDs. Its rates are reported alone and were never summed into
   CCL3.
6. **Ruler A's dynamic range is narrower here than in prior work.** Ambient
   ceiling 0.092 versus a true-NK floor of 0.021 is a 4.4× separation, where
   the earlier calibration had 0.402 versus 0.012 (33×). The ruler still
   separates the headline genes cleanly (0.021–0.025), but the margin is
   smaller than the prior numbers would suggest.

---

## 4. What was not done

No pathway scoring, no pseudotime, no ligand–receptor inference, no
re-clustering, no re-labelling, no cross-dataset matrix merging, and no
cell-level p value anywhere. `NETSKAR2024` was closed without promotion
(D1). `HECA` was not reached.
