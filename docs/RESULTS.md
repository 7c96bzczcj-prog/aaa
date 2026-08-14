# RESULTS — dNK chemokine profile, VT2018

Dataset: `VT2018` (E-MTAB-6701, Vento-Tormo et al., Nature 2018).
Panel: `chemokine_panel_v1` (68 genes, 98.5% matched).
Criteria: [`PREREGISTRATION.md`](PREREGISTRATION.md) (v1.0, frozen at commit
`175713a` before any expression value was read), as amended by
[`v1.1`](PREREGISTRATION_v1.1.md) and [`v1.2`](PREREGISTRATION_v1.2.md). Every
amendment discloses its direction of effect; v1.1's three loosened the reading,
v1.2's four tightened it.

> **The methodological result of this project, stated first because it frames
> everything below.** Four statistical instruments were withdrawn across
> v1.1–v1.3 — BH over the signed-rank p, ruler B's fixed denominator, the FDR
> criterion, and z as a standalone ranking. Each failed the same way: its
> output was determined by a nuisance parameter (n, the denominator lineage,
> the null-set size, the baseline detection rate) rather than by the effect.
> **What survived is exactly what never needed a null distribution:**
> donor-unanimous direction, effect size in percentage points, and
> reproduction inside a fixed cluster. No p value or q value in this document
> adjudicates anything.

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

## 1. What the primary targets do, with no decision rule applied

The pre-registered FDR criterion is **withdrawn** (v1.3 A8): the number of rows
it passes is set by the null-set size, not by the effects — the same data gave
6, then 0, then 8 passes as N alone changed. No q value appears below.

Descriptive quantities only, sorted by effect. `z` is against each gene's own
null set, matched on the **reference arm's baseline detection rate** (v1.3 A9),
and is a standardised effect size that is **never converted to a p**.

| contrast | gene | n | effect (pp) | z | donors concordant | baseline det. (ref arm) | ceiling | ruler A |
|---|---|---|---|---|---|---|---|---|
| dNK1→dNK3 | **CXCR4** | 5 | **+57.9** | +7.97 | 5/5 | 0.176 | clean | 0.048 |
| dNK1→dNK3 | **CCL5** | 5 | **+57.3** | +6.72 | 5/5 | 0.413 | lower bound | 0.025 |
| dNK2→dNK3 | **CXCR4** | 5 | **+52.5** | +10.04 | 5/5 | 0.194 | clean | 0.048 |
| dNK1→dNK2 | **XCL1** | 6 | **+39.3** | +7.74 | 6/6 | 0.405 | lower bound | 0.021 |
| dNK2→dNK3 | **CCL5** | 5 | **+35.5** | +5.72 | 5/5 | 0.483 | lower bound | 0.025 |
| dNK1→dNK2 | **XCL2** | 6 | **+32.9** | +5.75 | 6/6 | 0.583 | lower bound | 0.021 |
| dNK1→dNK3 | **XCL1** | 5 | **+31.2** | +4.09 | 5/5 | 0.405 | lower bound | 0.021 |
| dNK1→dNK2 | **CCL5** | 6 | **+23.3** | +4.36 | 6/6 | 0.413 | clean | 0.025 |
| dNK1→dNK2 | *SELL* | 6 | *+4.8* | *+10.55* | 6/6 | **0.005** | floor | — |
| dNK2→dNK3 | *CXCR6* | 5 | *+2.3* | *+4.78* | 5/5 | **0.002** | floor | — |
| dNK2→dNK3 | *S1PR5* | 5 | *+1.5* | *+10.34* | 5/5 | **0.000** | floor | — |
| dNK1→dNK3 | *S1PR5* | 5 | *+1.5* | *+5.81* | 5/5 | **0.001** | floor | — |
| dNK2→dNK3 | *CCR9* | 5 | *+1.0* | *+5.45* | 5/5 | **0.001** | floor | — |

**The italicised rows are why z may not be read alone.** SELL and S1PR5 carry
the two largest z values in the study on effects of 4.8 and 1.5 percentage
points, because their reference-arm baseline detection is 0.5% and 0.0% and
their null sets are correspondingly narrow. Measured across all targets:

| baseline detection | median null SD (pp) | what z = 5 buys |
|---|---|---|
| > 30% | 6.41 | **32.1 pp** |
| 2–10% | 2.06 | 10.3 pp |
| < 0.5% | 0.33 | **1.63 pp** |

**A 20-fold difference in the effect that the same z represents.** Matching the
nulls on detection rate rather than CPM — the methodologically correct axis —
made this sharper rather than milder (S1PR5 went from z = +5.01 to +10.34), so
it is a property of a bounded scale near its boundary, not of the matching
(D25). No effect-size filter was added afterwards to remove those rows.

### Ruler B is `unmeasurable` for every headline gene

Not for lack of signal: **NK is itself the largest total contributor** of each
to the decidual ambient pool (v1.2 A7), so an ambient argument has no competing
source to run against. That is supportive — but the pre-registered "both rulers
positive" bar **is not met by any headline gene**, and the alternative reading
(ruler A plus the circularity test) is **post hoc** and labelled as such.

"Largest contributor to the pool" is a statement about **total counts**; NK is
roughly 70% of decidual leukocytes here, so it does **not** mean highest
per-cell expression.

Ruler A ceiling 0.092 — the **minimum** over ambient controls, and set entirely
by the haematopoietic class, which is the class the estimator handles worst
(D23). True-NK floor 0.021; ruler B band 0.0538.

---

## 2. Conclusions

Each carries: evidence, n, effect, both rulers, and what it does **not**
support.

### C1. XCL1 and XCL2 are dNK1-low and dNK2/dNK3-high. This reproduces Huhn 2020 at transcript level.

- **Evidence:** `donor_level_tests.tsv` lines 7, 8, 45, 46;
  `null_distribution.tsv` lines 7, 45; `soup_calibration.tsv` lines 74, 75.
- **n = 6** (dNK1 vs dNK2), **n = 5** (dNK1 vs dNK3). Every donor concordant.
- **XCL1 ≥ +39.3 pp** dNK1→dNK2 (detection 51.8% → 91.1%) and **≥ +31.2 pp**
  dNK1→dNK3; XCL2 **≥ +32.9 pp**. q = 0.033 and 0.050, z = +5.9 to +9.6.
- **Ruler A = 0.021**, identical to the true-NK gene floor and far below the
  0.092 ambient ceiling — the furthest from ambient of any gene tested.
  **Ruler B = `unmeasurable`**: NK is the largest total contributor of XCL1 to
  the decidual pool, so no competing source exists to test pickup against.
- **Survives the circularity test** (§3), including within the dNK1 cluster.
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
  +23.3 pp (n = 6, z = +5.65, q = **0.0499**, marginal), +35.5 pp (n = 5,
  z = +9.84, q = 0.033), +57.3 pp dNK1→dNK3 (n = 5, z = +10.87, q = 0.033).
- Ruler A 0.025 (vs 0.092 ceiling); ruler B `unmeasurable` — NK is the largest
  contributor of CCL5 to the pool, ahead of T cells.
- **Partly survives the circularity test** (§3): the CD103 axis reproduces it
  inside the dNK2 cluster (+6.0 pp, 3/3 donors), but the between-cluster
  magnitude does not.
- **dNK1 vs dNK2 is the one contrast with neither arm saturated.** The two
  steps involving dNK3 (90.4%) are lower bounds.

**This does not support:** more CCL5 protein from dNK3; that the gradient is
functional; nor any magnitude for the dNK3 steps.

### C3. CXCR4 carries the study's largest z and its weakest claim.

- **Evidence:** `donor_level_tests.tsv` (CXCR4 rows in both dNK3 contrasts);
  `circularity_within_label_tests.tsv`.
- **+57.9 pp** dNK1→dNK3 (n = 5, z = +10.98) and **+52.5 pp** dNK2→dNK3
  (n = 5, **z = +14.56, the largest in the panel**), both q = 0.033. Neither
  arm saturated, so the magnitudes are readable as stated.
- Ruler A 0.048 (below the 0.092 ceiling); ruler B `unmeasurable`.
- **But it fails the circularity test worst of the four** (§3): only 0.30× of
  the effect survives marker-based regrouping, and **inside the dNK1 cluster
  the effect reverses sign** (−3.5 pp, 3/3 donors concordant).

**This does not support:** a homing interpretation — nothing here measures
migration, surface protein or receptor occupancy. And it does not support the
between-cluster magnitude as biology: **a large part of it is a restatement of
how the dNK3 cluster was drawn.** A z of +14.56 measured across labels is not
evidence against that; the within-label reversal is evidence for it.

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

## 3. Circularity: the credibility ordering inverts the z ordering

dNK1/dNK2/dNK3 are Vento-Tormo's own **unsupervised whole-transcriptome**
clusters, characterised afterwards by differential expression — and CCL5,
CXCR4 and XCL1 are among the genes that characterisation named. Measuring them
across the same labels partly restates the clustering. (Verified against the
published methods: the clustering was unsupervised, not marker-guided, so the
contamination is partial rather than total.)

`src/circularity_test.py` (authorised by v1.2 A4) discards the labels and
re-groups the same cells on **CD39/CD103 equivalents (ENTPD1/ITGAE)** — markers
no target gene informs. `CD160` and `KLRB1` are barred from gating: they are
themselves dNK3 characterisation genes.

**Effect retained after marker regrouping:**

| gene | gate1→2 | gate1→3 | gate2→3 | median retained |
|---|---|---|---|---|
| XCL1 | +24.3 (4/4 donors) | +27.7 (4/4) | +5.5 | **0.62×** |
| XCL2 | +17.6 (4/4) | +22.9 (4/4) | +4.5 | **0.53×** |
| CCL5 | +17.1 (4/4) | +26.3 (4/4) | +6.8 | **0.46×** |
| CXCR4 | +8.6 (4/4) | +17.6 (4/4) | +8.7 (5/5) | **0.30×** |

n falls to 4–5, so `min_achievable_p` is 0.125 and **no p value in this section
is interpretable**. Direction and effect size are all that is claimed.

Shrinkage alone does **not** settle it. The gates are impure (gate3 is only
35% dNK3) and dilution shrinks everything. It is tempting to argue that
dilution is gene-agnostic, so a 0.62×–0.30× spread must mean circularity —
**that argument is too weak to use.** Detection rates are bounded in [0,1], so
dilution compresses each gene by an amount depending on where the
contaminating cells' rate sits between the two arms; it is not proportional.
The within-label test is what carries the conclusion.

**The within-label test settles it.** A gate contrast computed *inside one
published cluster* cannot restate how that cluster was drawn, and is immune to
between-gate dilution:

| inside | contrast | gene | n | effect | donors concordant |
|---|---|---|---|---|---|
| dNK1 | CD39⁺→CD39⁻ | **XCL1** | 3 | **+6.7 pp** | 3/3 |
| dNK1 | CD39⁺→CD39⁻ | **XCL2** | 3 | **+5.7 pp** | 3/3 |
| dNK1 | CD39⁺→CD39⁻ | CCL5 | 3 | +2.5 pp | 3/3 |
| dNK1 | CD39⁺→CD39⁻ | **CXCR4** | 3 | **−3.5 pp** | 3/3 |
| dNK2 | CD103⁻→CD103⁺ | **CCL5** | 3 | **+6.0 pp** | 3/3 |
| dNK2 | CD103⁻→CD103⁺ | XCL1 | 3 | −5.2 pp | 3/3 |

n = 3, `min_achievable_p` = 0.25 — again, **no p value here means anything**.

**The single strongest observation in this study is in that table's first and
fourth rows.** Within the *same* contrast, on the *same* cells, XCL1 moves
**+6.7 pp** while CXCR4 moves **−3.5 pp**, both with all three donors
concordant. A depth or dropout artefact in the gate pushes every gene the same
way — if CD39⁻ cells were merely shallower, every detection rate would fall
together. **Opposite signs within one contrast exclude that common artefact**,
which no single effect size and no z score can do.

**Resulting tiers:**

| tier | rows | why |
|---|---|---|
| **cleanest** | **XCL1/XCL2, dNK2 > dNK1** | survives regrouping (0.53–0.62×) *and* survives inside dNK1 (+6.7/+5.7, 3/3); and it is a **cross-modality reproduction of Huhn 2020**, which gated on CD39/CD103 protein by CyTOF |
| middle | CCL5, dNK2 > dNK1 | CCL5 is a dNK3 characterisation gene, but neither arm of this contrast is defined by it; survives inside dNK2 along CD103 |
| **most circular** | all dNK3 rows (CCL5/CXCR4/XCL1) | these are the genes that characterise dNK3; CXCR4 retains only 0.30× and **reverses inside dNK1** |

**So the study's largest z (CXCR4, +14.56) is its least trustworthy row, and
XCL1 dNK2 > dNK1 — a smaller z — is its strongest.** Direction is established
as non-circular for XCL1/XCL2 and CCL5; **magnitude is not established for
anything**, since within-label effects are ~17% of between-cluster ones.

---

## 4. dNK3 in full

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

## 5. Caveats

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
   source file ([D17](DECISIONS.md)). The v1.2 CHANGELOG called this harmless
   because five tissue controls already return 1.000; **that argument was
   wrong and is retracted** ([D21](DECISIONS.md)). The ceiling is the *minimum* over ambient
   controls (IGKC, 0.092), not the maximum. If PAEP would return anything
   between 0.048 and 0.092, the ceiling drops below CXCR4's 0.0484 and CXCR4's
   ruler-A verdict flips to `indistinguishable_from_ambient`. XCL1 (0.021) and
   CCL5 (0.025) keep a wide margin.
4. **Ruler B now ranks sources by total pool contribution**
   ([A7](PREREGISTRATION_v1.2.md)), after v1.1's per-cell-CPM ranking assigned
   XCL1's source to a 184-cell ILC3 population. Consequence: for every headline
   gene NK is itself the largest contributor, so ruler B returns
   `unmeasurable` and **the "both rulers positive" bar is met by none of them**.
5. **Ruler A's dynamic range is narrower here** than in prior work: ceiling
   0.092 vs floor 0.021 is 4.4×, where the earlier calibration had 33×.
6. **`CCL3L1` matched only via Ensembl** and its symbol maps to three Ensembl
   gene IDs; reported alone, never summed into CCL3.
7. **The dNK1 residual is not library composition** (D19): triple-positive rate
   against library T fraction gives Pearson r = +0.061 (p = 0.90) over 7
   libraries. That excludes a strong relationship, not a weak one. The other
   two candidates remain undecidable here.
8. **n = 11 conflict with the July OXPHOS result is unresolved but bounded**
   ([D15](DECISIONS.md)): no donor-level statistic from E-MTAB-6701 can have
   n > 7 (n > 6 for decidua), and merging the Smart-seq2 companion reaches
   only 9. One lookup settles it — see D15.

---

## 6. Verdict

Under v1.0 the study returned **zero** significant results, because BH over a
floor-limited p cannot reject at n = 6 — a fact the pre-registration
established before the data existed. Under v1.1–v1.2, with inference on the
empirical null, the ceiling rule per-arm, and multiplicity over all 81
analysis-A tests, **eight donor-unanimous contrasts across four genes clear
q ≤ 0.05**.

But the eight are not equal, and the circularity test (§3) reorders them
against their z scores:

**The claim that holds.** *XCL1 and XCL2 are detected in a larger fraction of
dNK2 than of dNK1 cells from the same donor, in all six donors, by at least
39.3 and 32.9 percentage points.* It survives marker-based regrouping, it
survives **inside** the dNK1 cluster where circularity is impossible, it sits
at the true-NK-gene floor on ruler A, and it is a **cross-modality
reproduction of Huhn 2020** — a CyTOF protein result recovered from transcript
detection with an independent grouping.

**The claim that is weaker than its z suggests.** Everything involving dNK3.
CCL5 and CXCR4 are the genes that *define* dNK3 in the source annotation;
CXCR4 retains only 0.30× under regrouping and **reverses inside dNK1**. Its
z = +14.56 measures agreement with the cluster definition as much as biology.

**What the study cannot do**, unchanged: quote a magnitude for any saturated
contrast, or for anything — within-label effects are ~17% of between-cluster
ones; say anything about protein; say anything about signalling; or meet the
pre-registered "both rulers positive" bar, since ruler B is `unmeasurable`
wherever NK is the pool's largest contributor.

**And the two rows at q = 0.0499** are inside the threshold by 0.0001. They
are not established.

## 7. Not done

No pathway scoring, pseudotime, ligand–receptor inference, re-clustering,
re-labelling, cross-dataset matrix merging, or cell-level p value.
`NETSKAR2024` closed without promotion (D1). `HECA` not reached. The Kitagawa
decomposition is implemented and tested but no comparison this round triggers
it.
