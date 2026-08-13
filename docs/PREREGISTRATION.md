# PREREGISTRATION — dNK chemokine profile, v1.0

**Status: FROZEN.** Committed before any expression value was read from any
count matrix. Changes require a new file `PREREGISTRATION_v1.1.md` plus a
CHANGELOG entry giving the reason and timestamp. Editing this file after
results exist invalidates the study.

Source specification: "dNK 趋化因子谱 — 执行规范 v1.0", 2026-08-13.
Frozen panel: [`panel/chemokine_panel_v1.tsv`](../panel/chemokine_panel_v1.tsv), 68 genes.

---

## 0. What is being asked

At **donor** level, after ambient-RNA (soup) calibration, measure the
chemokine ligand and receptor profile of human decidual NK subsets
(dNK1 / dNK2 / dNK3 / dNKp), using a dataset-agnostic pipeline driven
entirely by manifests.

Explicitly **not** asked: pathway scoring, pseudotime, cell–cell
communication inference (no CellPhoneDB, no CellChat), cross-dataset
matrix integration.

---

## 1. Primary readout

`detection_rate` = fraction of cells in a `donor × compartment × subset`
unit with **count ≥ 1** for a gene, measured **after** downsampling every
cell to a common UMI depth. `mean_cpm` is reported alongside but carries
no conclusion.

---

## 2. Power, established before analysis

The donor structure of the anchor dataset was read from
`meta_10x.txt` (cell annotations only — no expression values). It differs
materially from the values the specification recalled, and that changes
what the study can conclude. **These numbers are pre-analysis facts and
are recorded here so that no result can be reinterpreted against a
remembered n.**

| quantity | spec recalled | measured in E-MTAB-6701 |
|---|---|---|
| total cells | ~70,000 | **64,734** |
| dNK1+dNK2+dNK3 cells | 11,927 | **11,202** (11,932 including dNKp) |
| usable donors | 11 | **7 total; 6 with decidua; 4 with paired blood** |

The 10x experiment (E-MTAB-6701) contains donors D6–D12. The paper's
larger donor count includes the Smart-seq2 experiment (E-MTAB-6678),
which is a different platform and is barred from matrix-level merging by
**R9**.

### 2.1 Minimum achievable p (R2)

Two-sided Wilcoxon signed-rank on n non-tied donor pairs cannot return a
p below `2 / 2^n`:

| n donors | min achievable p | is p < 0.05 reachable? |
|---|---|---|
| 3 | 0.250 | no |
| 4 | 0.125 | **no** |
| 5 | 0.0625 | **no** |
| 6 | 0.03125 | yes, only if every donor agrees in sign |
| 7 | 0.015625 | yes |

### 2.2 Consequences that are fixed in advance

**Main analysis A (within-decidua, across subsets).** At most **6**
donors (D6, D7, D8, D9, D10, D12). Ceiling n = 6 → best possible raw
p = 0.03125, attainable only when all six donors move the same way.
Pre-downsampling cell counts show dNK3 is the limiting subset: D12 has
29 dNK3 cells, below the 30-cell admission threshold, and D10 has 37,
which is unlikely to survive QC plus depth-matching. **The dNK1-vs-dNK3
and dNK2-vs-dNK3 contrasts are therefore expected to run at n = 4–5,
where p < 0.05 is mathematically unreachable.** Only dNK1-vs-dNK2 is
expected to retain n = 6.

**Multiplicity.** BH correction is applied across the primary-target set
(27 genes: 15 ligands + 12 receptors; XCR1 is a control, not a
hypothesis). For a q ≤ 0.05 at the n = 6 floor of p = 0.03125, BH
requires at least `ceil(27 × 0.03125 / 0.05)` = **17 of 27 genes** to sit
at that floor simultaneously. A single strong gene cannot reach q ≤ 0.05
at this n, no matter how large its effect.

**Secondary analysis B (decidua vs paired blood).** 4 donors
(D6, D7, D8, D9) → n = 4 → **min achievable p = 0.125. p < 0.05 is
mathematically impossible for analysis B, before any data are seen.**
Analysis B is therefore pre-declared **descriptive with effect sizes and
confidence intervals only**; no significance claim from it will be made
or accepted.

**Analysis B against CD56-bright blood NK is pre-declared unmeasurable.**
Blood `NK CD16-` counts are D6 = 20, D7 = 3, D8 = 72, D9 = 125 — only 2
donors clear 30 cells. Since **R-6.3 forbids merging CD56-bright with
CD56-dim**, the comparison that would be biologically correct (dNK vs
pbNK_CD56bright) runs at n = 2 and will be reported as
`unmeasurable`, not as a null result.

**This is a power statement, not a result.** It commits the study to
reporting effect sizes and floors rather than mining for stars.

---

## 3. Decision rules — main analysis A

Outcome is one of exactly three, decided by the rules below and by
nothing else.

### A-positive
Donor-level differences among dNK subsets in `XCL1 / CCL5 / CCL3 / CCL4`
that are (a) significant at the n available, (b) `combined_verdict =
positive` on both soup scales, and (c) directionally consistent with the
literature (dNK2/dNK3 > dNK1).

*Means:* subset-specific chemokine output holds **at transcript level**.
*Does not mean:* protein secretion; actual signalling to EVT or cDC1.

### A-negative
Differences not significant, or not passing both soup scales.

*Means:* the published subset differences rest largely on
post-stimulation intracellular protein staining (CyTOF) and do not
reproduce at transcript level. **This does not refute the literature** —
it says transcript abundance is the wrong readout, and any future design
using transcripts as the readout must be reconsidered.

### A-uninterpretable
Target genes fall largely on the ceiling/floor (detection rate > 0.85 or
< 0.05) or return `unmeasurable`.

*Means:* this dataset lacks the depth to answer the question. Move to
HECA or a deeper dataset. **Forcing an interpretation is barred.**

### Pre-committed n-floor clause
If a contrast's n makes p < 0.05 unreachable (§2.1), that contrast may
**not** be reported as A-negative on the basis of its p value. It is
reported as **A-underpowered** with `on_test_floor = TRUE`, its effect
size, and its CI. Absence of significance at n = 4 is not evidence of
absence.

---

## 4. Decision rules — secondary analysis B

**Readability gate.** B's cross-compartment comparison is interpretable
only if `XCR1` detection in the NK gate differs by **< 2 percentage
points** between decidua and blood. Otherwise the ambient backgrounds are
not comparable and B is **reported without conclusions**.

Every B row carries `cross_compartment_soup_caveat = TRUE`
unconditionally: decidual soup is stroma/trophoblast-rich, blood soup is
platelet/erythroid-rich, so the ambient contribution to CCL3/CCL4/CCL5
has entirely different provenance on the two sides.

---

## 5. Soup calibration — the two rulers

- **Ruler A (conservative ceiling):** soup fraction of pure-ambient genes
  measured inside the NK gate. Prior work put this at 0.402–1.000 against
  a true value of 1.0; **0.402** is the conservative ceiling. A target
  gene whose soup fraction approaches the ceiling is indistinguishable
  from pure ambient.
- **Ruler B (pickup band):** NK/myeloid CPM ratio. Prior work put the
  pure-ambient pickup band at **≤ 0.118**. A target gene must sit
  materially above that band to count as really expressed.
- **Floor:** true-NK genes previously sat at soup fraction 0.012–0.013.

Both rulers are **re-calibrated inside this dataset**; the prior numbers
are order-of-magnitude reference only. Verdicts:

| both rulers pass | one ruler passes | neither | reference lineage not enriched over NK |
|---|---|---|---|
| `positive` | `weak` | `indistinguishable_from_ambient` | `unmeasurable` |

`unmeasurable` is **not** `negative`. When the reference lineage does not
express the gene far above NK, the ratio test has no power and must say
so.

---

## 6. Stop conditions, committed in advance

1. **Purity stop.** If triple-positive (`TRBC2⁺CD3E⁺CD3D⁺`) T-cell
   contamination exceeds **15%** in any dNK subset, stop and resolve the
   annotation before running main analysis.
2. **Ceiling/floor stop (R5).** Genes with detection rate > 0.85 or
   < 0.05 are flagged `uninterpretable_ceiling` / `uninterpretable_floor`
   and removed from conclusions.
3. **Differential loss stop (R6).** If downsampling retention differs by
   > 10 percentage points across compared groups, sensitivity analysis at
   three depth floors is mandatory.
4. **Match-rate stop (R11).** Panel match rate < 90% → hard exception,
   no silent continuation.
5. **Admission threshold.** `donor × subset` units with < 30 cells after
   downsampling are excluded and listed in `out/<id>/excluded_units.tsv`
   with the reason.

---

## 7. Barred in advance

- Any cell-level p value, anywhere, including appendices.
- Merging count matrices across datasets before the main analysis.
- Pathway scoring / GSEA as evidence of expression.
- Any ligand–receptor prediction tool.
- Re-clustering or re-labelling cell types. A contradiction with the
  published annotation is **reported, not repaired**.
- Summing `CCL3L1`/`CCL4L2` into `CCL3`/`CCL4`.
- Any statement of the form "X is expressed in dNK" without soup
  calibration attached.
- Signal-to-noise metrics whose denominator is a single-sample
  difference (R7).
- Amending this document because a result is unwelcome.

---

## 8. Registration

Frozen at commit time, before `detection_rates.py`, `donor_tests.py`,
`soup_calibration.py`, or `null_distribution.py` had been run on any
expression data. The only data read beforehand were cell-annotation
metadata (`meta_10x.txt`, 4 columns: donor, location, cluster,
annotation) and the remote `obs` table of the NETSKAR2024 atlas — neither
contains expression values.
