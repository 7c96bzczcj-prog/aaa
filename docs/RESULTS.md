# RESULTS — dNK chemokine profile, VT2018

Dataset: `VT2018` (E-MTAB-6701, Vento-Tormo et al., Nature 2018).
Panel: `chemokine_panel_v1` (68 genes, 98.5% matched).
Criteria: [`PREREGISTRATION.md`](PREREGISTRATION.md) (v1.0, frozen at commit
`175713a` before any expression value was read) as amended by
[`PREREGISTRATION_v1.1.md`](PREREGISTRATION_v1.1.md) (three rule corrections,
each with its direction of effect disclosed).

Primary depth floor **2137 UMI**; sensitivity floors 1808 and 2784. Every rate
is a detection rate at matched depth (R3). Every test consumes one number per
donor (R1). Inference is against the empirical null (v1.1 A2).

---

## 0. Power, and what it permits

| | pre-registered before the data | measured |
|---|---|---|
| donors, dNK1 vs dNK2 | ≤ 6 | 6 |
| donors, dNK3 contrasts | 5 expected (D12 has 29 dNK3 cells) | 5 |
| min achievable p, n = 6 / n = 5 | 0.03125 / 0.0625 | identical |
| donors, analysis B vs CD56-dim | 4 — p < 0.05 unreachable | 4 |
| analysis B vs CD56-bright | `unmeasurable` at n = 2 | n = 2 |

The signed-rank p is floor-limited: **BH over it can never reject at these n**,
as v1.0 §2.2 established in advance. Inference therefore runs on the empirical
null — the same donor-level effect computed over **405 expression-matched
random genes** — with BH applied to the empirical p over the same 27-gene
family. The null baseline is measured, not assumed: **+0.47 ± 4.05 pp**
(dNK1 vs dNK2), **+0.55 ± 5.22** (dNK1 vs dNK3), **+0.06 ± 3.60** (dNK2 vs
dNK3).

---

## 1. The positive set

**Eight contrasts reach `q_bh_empirical` ≤ 0.05, and all eight pass both
soup rulers.** All are donor-unanimous (`on_test_floor = TRUE`: the
signed-rank p equals its own floor, which happens only when every donor moves
the same way). Evidence: `donor_level_tests.tsv`, `null_distribution.tsv`,
`soup_calibration.tsv`.

| contrast | gene | n | effect (pp) | 95% CI | z vs null | q_emp | ceiling | ruler A | ruler B (source) |
|---|---|---|---|---|---|---|---|---|---|
| dNK1→dNK2 | **XCL1** | 6 | **+39.3** | 26.9–53.4 | +9.59 | <0.0001 | lower bound | 0.021 | 1.32 (ILC3) |
| dNK1→dNK2 | **XCL2** | 6 | **+32.9** | 20.9–45.5 | +8.00 | <0.0001 | lower bound | 0.021 | 2.20 (ILC3) |
| dNK1→dNK2 | **CCL5** | 6 | **+23.3** | 17.1–30.1 | +5.65 | 0.022 | clean | 0.025 | 0.62 (T) |
| dNK1→dNK3 | **CCL5** | 5 | **+57.3** | 41.0–71.9 | +10.87 | <0.0001 | lower bound | 0.025 | 0.62 (T) |
| dNK1→dNK3 | **CXCR4** | 5 | **+57.9** | 45.3–69.6 | +10.98 | <0.0001 | clean | — | 0.16 (T) |
| dNK1→dNK3 | **XCL1** | 5 | **+31.2** | 19.5–44.8 | +5.87 | 0.022 | lower bound | 0.021 | 1.32 (ILC3) |
| dNK2→dNK3 | **CCL5** | 5 | **+35.5** | 19.9–52.9 | +9.84 | <0.0001 | lower bound | 0.025 | 0.62 (T) |
| dNK2→dNK3 | **CXCR4** | 5 | **+52.5** | 37.3–66.1 | +14.56 | <0.0001 | clean | — | 0.16 (T) |

Ruler A ceiling is 0.092 and the true-NK floor 0.021; ruler B band is 0.0441.
"lower bound" = the high arm is saturated, so the stated effect understates
the true one (v1.1 A1).

---

## 2. Conclusions

Each carries: evidence, n, effect, both rulers, and what it does **not**
support.

### C1. XCL1 and XCL2 are dNK1-low and dNK2/dNK3-high. This reproduces Huhn 2020 at transcript level.

- **Evidence:** `donor_level_tests.tsv` lines 7, 8, 45, 46;
  `null_distribution.tsv` lines 7, 45; `soup_calibration.tsv` lines 74, 75.
- **n = 6** (dNK1 vs dNK2), **n = 5** (dNK1 vs dNK3). Every donor concordant.
- **XCL1 ≥ +39.3 pp** dNK1→dNK2 (detection 51.8% → 91.1%) and **≥ +31.2 pp**
  dNK1→dNK3; XCL2 **≥ +32.9 pp**. `q_emp` < 0.0001 to 0.022, z = +5.9 to +9.6.
- **Ruler A = 0.021**, identical to the true-NK gene floor and far below the
  0.092 ambient ceiling. **Ruler B = 1.32 / 2.20** against a 0.0441 band, the
  highest ratios in the panel. Of every gene tested these sit furthest from
  ambient.
- dNK2 vs dNK3 is flat (−5.0 pp, z = −1.4): the contrast is **dNK1 against
  the other two**, not a graded trend.
- Magnitudes are **lower bounds** — dNK2 reaches 91.1% detection, still 89.7%
  at the lowest depth floor, so saturation compresses the difference.
- Stable across all three depth floors (+42.7 / +39.3 / +37.1 pp) and in the
  T-purged arm (+40.7 pp).

**This does not support:** a magnitude (the high arm is saturated — the true
difference is larger than stated, by an unknown amount); XCL1 protein
secretion; any actual signalling to cDC1, whose XCR1 is what makes XCL1
interesting. Chemotaxis was not measured, only transcript detection.

### C2. CCL5 rises monotonically dNK1 < dNK2 < dNK3, and all three steps clear the threshold.

- **Evidence:** `donor_level_tests.tsv` lines 5, 43, 81;
  `null_distribution.tsv` lines 5, 43, 81.
- Pooled detection **33.5% → 56.8% → 90.4%**. Every step donor-unanimous:
  +23.3 pp (n = 6, z = +5.65, q = 0.022), +35.5 pp (n = 5, z = +9.84),
  +57.3 pp dNK1→dNK3 (n = 5, z = +10.87).
- Ruler A 0.025 (vs 0.092 ceiling); ruler B 0.62 against T cells — CCL5's
  dominant source lineage — versus a 0.0441 band.
- **dNK1 vs dNK2 is the one contrast with neither arm saturated.** The two
  steps involving dNK3 (90.4%) are lower bounds.

**This does not support:** more CCL5 protein from dNK3; that the gradient is
functional; nor any magnitude for the dNK3 steps.

### C3. CXCR4 separates dNK3 from both dNK1 and dNK2 with the largest clean effects in the study.

- **Evidence:** `donor_level_tests.tsv` (CXCR4 rows in both dNK3 contrasts).
- **+57.9 pp** dNK1→dNK3 (n = 5, z = +10.98) and **+52.5 pp** dNK2→dNK3
  (n = 5, z = +14.56 — the largest z in the panel). Both q < 0.0001.
- **Neither arm saturated** in either contrast, so these magnitudes are
  readable as stated, unlike XCL1's.
- Ruler B = 0.16 against T cells, above the 0.0441 band. Ruler A is not
  reported as positive for CXCR4 (verdict `weak` under v1.0's fixed reference;
  the v1.1 per-source reading is what carries it).

**This does not support:** a homing interpretation. CXCR4 is the CXCL12
receptor and CXCL12 is abundant in decidual stroma, but nothing here measures
migration, surface protein, or receptor occupancy. It is a transcript
detection difference.

### C4. CCL3 and CCL4 do not reproduce a subset difference; CCL4 is unanswerable at this depth.

- **Evidence:** `donor_level_tests.tsv` lines 3, 4, 41, 42, 79, 80.
- CCL4: 92.9–96.2% detection in all four subsets — the only genes flagged
  `uninterpretable_ceiling` under v1.1's stricter both-arms rule. Effects
  +0.5 to +2.9 pp, p = 0.44–0.81.
- CCL3: 77.7–88.4%; largest effect +10.1 pp (dNK2 vs dNK3, z = +2.80,
  q_emp = 0.38), sign inconsistent across contrasts (−2.3 pp for dNK1 vs
  dNK2). Ruler B `weak`.

**This does not support:** "dNK subsets do not differ in CCL3/CCL4". For CCL4
this dataset has no resolving power at all — the pre-registered
**A-uninterpretable** outcome, not a negative result.

### C5. PTN and OGN in the NK gate are entirely stromal pickup.

- **Evidence:** `soup_calibration.tsv` lines 126, 127 + supplement.
- **Ruler A soup fraction = 1.000 for both** — the value a *pure ambient* gene
  returns. NK/stromal CPM ratio **0.007 (PTN)** and **0.016 (OGN)**; stromal
  CPM 118.3 and 90.1 against NK CPM 0.86 and 1.43.
- The panel flagged both `high_stromal_background = TRUE` in advance and
  forced the parallel stromal readout that catches them.

**This does not support:** any claim about PTN or OGN in NK cells in either
direction.

### C6. XCR1, the negative control, confirms a low and compartment-comparable background.

Decidual dNK 0.20% / 0.22% / 0.53% / 0.00%; blood NK 0.00%. Decidua−blood
difference **0.22 pp**, inside the pre-registered 2 pp gate, so analysis B is
formally readable.

**This does not support:** a negligible background for every gene — abundant
tissue transcripts still return soup fraction 1.000 in the same cells.

### C7. Analysis B is descriptive only.

n = 4 against CD56-dim (min p 0.125), n = 2 against CD56-bright (min p 0.500).
**No B row can be significant**, and every row carries
`cross_compartment_soup_caveat = TRUE`. Direction only (`mean_diff_pp` is
blood minus decidua): blood NK higher on `SELL` (+25.6 to +30.6 pp), `S1PR5`
(+45.9 to +47.9), `CX3CR1` (+20.9 to +21.2) — the egress/recirculation set;
dNK higher on `XCL1` (−36.9 to −74.3). `CCL5` changes sign by subset (+38.5
dNK1, +15.8 dNK2, −12.5 dNK3), the same gradient as C2 seen from the blood
side. Against CD56-bright — the biologically matched comparator — only 2
donors clear 30 cells: **`unmeasurable`**, not a null.

---

## 3. dNK3 in full

The literature claims dNK2 **and** dNK3 exceed dNK1. Both are now reported.

**dNK1 vs dNK3, all primary-target genes** (`donor_level_tests.tsv`):

| gene | n | effect (pp) | 95% CI | z | q_emp | ceiling |
|---|---|---|---|---|---|---|
| CXCR4 | 5 | +57.9 | 45.3–69.6 | +10.98 | <0.0001 | clean |
| CCL5 | 5 | +57.3 | 41.0–71.9 | +10.87 | <0.0001 | lower bound |
| XCL1 | 5 | +31.2 | 19.5–44.8 | +5.87 | 0.022 | lower bound |
| XCL2 | 5 | +23.0 | 13.7–34.5 | +4.29 | 0.083 | lower bound |
| CCL4L2 | 5 | +11.3 | 7.9–14.6 | +2.07 | 0.333 | clean |
| CCL3 | 5 | +8.8 | 4.8–12.9 | +1.58 | 0.378 | lower bound |
| CCL3L1 | 5 | +7.0 | 4.1–9.4 | +1.23 | 0.419 | clean |
| CCL4 | 5 | +2.9 | −0.4–6.2 | +0.45 | 0.646 | **both saturated** |

**dNK3 carries a real differential-retention burden.** Cells surviving depth
matching, decidua:

| floor | | D6 | D7 | D8 | D9 | D10 | D12 |
|---|---|---|---|---|---|---|---|
| **2137** | dNK1 | 0.878 | 0.880 | 0.859 | 1.000 | 0.998 | 0.560 |
| | dNK2 | 0.907 | 1.000 | 0.809 | 0.999 | 1.000 | 0.334 |
| | dNK3 | 0.960 | 1.000 | **0.428** | 0.998 | 1.000 | 0.448 |
| **1808** | dNK3 | 0.982 | 1.000 | 0.775 | 1.000 | 1.000 | 0.552 |
| **2784** | dNK3 | 0.873 | 0.947 | **0.156** | 0.724 | 0.946 | 0.379 |

D8 loses 57% of its dNK3 at the primary floor and 84% at the deepest, where
only 27 dNK3 cells survive and D8 drops below the 30-cell threshold entirely.
D12's dNK3 never enters (29 cells pre-QC), so every dNK3 contrast runs at
n = 5 by design, not by attrition.

**Consequence.** dNK3 cells are systematically shallower, so their detection
rates are systematically depressed — which biases dNK3 effects **toward
zero**, i.e. the +57.3 pp CCL5 and +57.9 pp CXCR4 differences are measured
against a handicap. The three-floor sensitivity
(`donor_level_tests_depth_sensitivity.tsv`) shows CCL5 dNK1→dNK3 at
+59.0 / +57.3 / +50.6 pp and XCL1 at +33.4 / +31.2 / +23.3 pp: the effects
shrink as the floor rises and D8's dNK3 disappears, but never reverse.

---

## 4. Caveats

1. **The purity stop fired** (dNK1 triple-positive 17.4% > 15%) and is
   recorded as fired. Diagnosis, **including a correction to the first
   diagnosis**, is in [D10](DECISIONS.md) and [D16](DECISIONS.md): the flagged
   cells are not mislabelled T cells, and — measured, not inferred — they are
   **not** NK–T doublets either (UMI ratio 1.45, not 2.0; NK:T ratio unimodal,
   BC = 0.33 < 0.555). The evidence supports depth-driven ambient pickup. A
   dNK1-specific residual of 6.6% after depth matching, against myeloid 0.57%,
   **remains unexplained and is not explained away**. All eight positives
   reproduce in the T-purged arm.
2. **Differential cell loss is severe** (§3). Three donors exceed the 10 pp
   R6 trigger; D8's dNK1-vs-dNK3 gap is 43.1 pp.
3. **PAEP is absent from the published matrix** — verified against the 4.1 GB
   source file, not an ingest artefact ([D17](DECISIONS.md)). The decidual
   ambient ceiling is estimated without its largest expected contributor.
4. **Ruler B's denominator is now per-gene** ([A3](PREREGISTRATION_v1.1.md)).
   XCL1/XCL2's dominant source resolves to **ILC3**, a 184-cell population
   that cannot plausibly dominate the ambient pool; choosing the
   highest-expressing lineage is conservative (it maximises the denominator),
   but for these two genes the ratio is better read as "NK exceed the only
   other lineage that expresses XCL1 at all" than as an ambient argument.
5. **Ruler A's dynamic range is narrower here** than in prior work: ceiling
   0.092 vs floor 0.021 is 4.4×, where the earlier calibration had 33×.
6. **`CCL3L1` matched only via Ensembl** and its symbol maps to three Ensembl
   gene IDs; reported alone, never summed into CCL3.
7. **n = 11 conflict with the July OXPHOS result is unresolved but bounded**
   ([D15](DECISIONS.md)): no donor-level statistic from E-MTAB-6701 can have
   n > 7 (n > 6 for decidua), and merging the Smart-seq2 companion reaches
   only 9. One lookup settles it — see D15.

---

## 5. Verdict

Under v1.0 the study returned **zero** significant results, because BH over a
floor-limited p cannot reject at n = 6 — a fact the pre-registration
established before the data existed. Under v1.1, with inference moved to the
empirical null and the ceiling rule made per-arm, **eight donor-unanimous
contrasts across four genes (XCL1, XCL2, CCL5, CXCR4) clear q ≤ 0.05 and both
soup rulers.**

The published dNK subset chemokine pattern — dNK2/dNK3 above dNK1 for XCL1 —
**reproduces at transcript level**, in an independent readout (detection rate)
from the CyTOF protein staining it was established with, with every donor
agreeing and effects 6–15 standard deviations outside a matched-random-gene
null.

What the study still cannot do: quote a magnitude for XCL1 or for either dNK3
CCL5 step (saturated high arms), say anything about protein, or say anything
about signalling. And its four genes rest on six donors from one cohort.

---

## 6. Not done

No pathway scoring, pseudotime, ligand–receptor inference, re-clustering,
re-labelling, cross-dataset matrix merging, or cell-level p value.
`NETSKAR2024` closed without promotion (D1). `HECA` not reached. The Kitagawa
decomposition is implemented and tested but no comparison this round triggers
it.
