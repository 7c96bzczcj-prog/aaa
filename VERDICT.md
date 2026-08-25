# VERDICT — BM-PB-NK preflight

spec `BM-PB-NK-preflight-v1.0`, graded and worded per the v1.1 addendum.
Every number traces to `results/` via `results/PROVENANCE.md`.

---

## 0. Limits of this preflight

**Material.** Not all admitted material is aspirate, and the addendum's premise
is corrected here rather than reproduced: HCA Census of Immune Cells and
GSE120221 are iliac **aspirate** mononuclear cells; GSE233304's marrow arm is
**surgically taken calvarial bone** for ten donors plus iliac BMMC for three.
No public dataset carries vertebral flush.

**Non-transferability.** Blood dilution is present in aspirate and scRNA-seq has
no purity criterion equivalent to IGRA/N, so task A cannot separate ambient-RNA
confounding from real blood NK carried in by dilution; both compress the BM–PB
difference in the same direction. **The task A negative therefore does not
reverse onto any other material. Only a positive would have been a transferable
warning.** Sampling method is not testable here at all.

Nothing below is a statement about any cohort other than the public datasets
named. All three tasks stop at the transcriptome.

**Evidence grading (v1.1 §1).** The question is the difference between two
compartments *of one person*. Grading follows design, not n:

| dataset | design | usable donors | grade |
|---|---|---|---|
| GSE233304 | **same-donor paired** BM/PB | 4 (P1, P15, P16 glioblastoma; c9 non-neoplastic control) | **higher** |
| HCA Census of Immune Cells | **unpaired**, different individuals per compartment | 8 BM vs 8 PB | lower — magnitude reference only |

HCA's ρ_BM vs ρ_PB carries between-donor variance as well as between-compartment
variance and cannot answer the question directly. Where the two conflict,
GSE233304 governs.

---

## 1. Where each task lands in its pre-written table

| task | cell of the pre-written table | one sentence |
|---|---|---|
| **A** | **row 1 — 通过** (both datasets) | Δ_ambient is near zero on every testable marker and ρ_BM ≈ ρ_PB, so compartment soup is not a systematic confounder in these datasets. |
| **B** | **row 3 — 不可测** | CXCR6 detection is insufficient, so ltNK cannot be reliably defined at the transcriptome layer — the direct implication is that ltNK must be read by flow, not by transcriptome. |
| **C** | **row 1 — P3 成立** | No intermediate state before or after QC and X ≤ 0.0093%, so in marrow of the kind specified in §4 the marrow is not where NK output is read, and a marrow arm there should ask about residency rather than development. |

---

## 2. Task A — detail

Instrument acceptance passes in every compartment of both datasets (27–42
testable foreign genes, 0.70–1.00 of them within 3× of ρ against a required
2/3), so the numbers may be read.

| dataset | design | ρ_BM | ρ_PB | ρ ratio | p | max abs Δ_ambient (testable genes) | genes over the pre-set 0.10 effect |
|---|---|---|---|---|---|---|---|
| GSE233304 | paired, n=4 | 0.0090 | 0.0028 | 3.23 | 0.32 | 0.0081 (11 genes) | **0** |
| HCA_ICA | unpaired, 8 vs 8 | 0.0175 | 0.0154 | 1.13 | 0.39 | 0.0076 (14 genes) | **0** |

Three HCA genes clear FDR (CD69, CD96, KLRC1) at effect sizes near 0.005 —
about a twelfth of what the frozen table treats as a difference — so they do not
move the verdict.

**Recorded conflict (v1.1 §1).** The two datasets disagree on sign: all 11
testable genes run BM-higher in GSE233304, 10 of 14 run PB-higher in HCA. Both
magnitudes sit near a twelfth of the threshold, so the disagreement lies inside
the noise rather than between the verdicts. The higher-grade dataset
(GSE233304) is the one reported as governing; it also carries the larger ρ ratio
(3.23) at n = 4, which the frozen rule does not treat as fired because it is not
significant.

**One premise of the task-A argument does not reproduce here.** The claim that
NK carries a uniquely large soup fraction (ρ_NK 0.151 vs ρ_CD8T 0.058, ratio
2.60) was established on dissociated tumour tissue. In marrow and blood, on the
same estimator and the same marker set, ρ_NK/ρ_T is 0.44–1.86 across five
dataset × compartment cells and ρ itself is 0.004–0.019, an order of magnitude
below the tumour figure (`results/rho_by_lineage.tsv`).

---

## 3. Task B — detail

The section 5 gate is checked in six dataset × compartment cells and passes in
none.

| dataset | compartment | donors | NK cells | CXCR6 in NK | ambient floor | CD69 in NK | ambient floor |
|---|---|---|---|---|---|---|---|
| GSE120221 | BM | 8 | 2,576 | 0.62% | 0.09% | 60.5% | **30.7%** |
| GSE233304 | BM | 7 | 726 | 0.83% | 1.95% | 87.5% | 10.9% |
| GSE233304 | BMMC | 3 | 2,337 | 1.28% | 1.37% | 85.9% | 6.98% |
| GSE233304 | PB | 8 | 1,212 | 1.32% | — | 81.4% | — |
| HCA_ICA | BM | 8 | 22,362 | 0.44% | 0.03% | 32.9% | 4.24% |
| HCA_ICA | PB | 8 | 32,780 | 0.29% | — | 25.3% | — |

Floor = the same gene's detection rate in erythroid cells, which express
neither. CXCR6 spans 0.29%–1.32% against a 5% admission floor and is at or
below its own floor in two cells. The age regression is present in
`results/ltnk_age_regression.tsv` only because the script computes it
unconditionally; its `gate_passed_in` column reads NONE and it is not read.

---

## 4. Task C — detail, and the boundary of the claim

| quantity | value |
|---|---|
| bone marrow cells examined | 374,305 post-QC, 441,171 pre-QC droplets, 8 donors, 71 libraries |
| early NK gate (CD34/KIT⁺ IL2RB⁺, no cytotoxic programme) | **25 before QC, 25 after QC** |
| exclusion bound X (95% upper) | **0.0093%**; **0.031%** after correcting for IL2RB detection (0.299 in mature NK) |
| assumption-light superset (CD34/KIT⁺ IL7R⁺, must contain any NK precursor) | 331 cells, X ≤ **0.0968%** |
| positive control, pro-B, structurally identical gate | 1,284 cells (0.343%) |
| positive control, erythroid progenitor, structurally identical gate | 2,981 cells (0.796%) |
| mature NK, same machinery | 13,437 cells (3.59%) |
| at Melsen 2022's own QC floor (≥1,000 genes/cell) | 20 of the 25 early-NK cells survive, while only 18.1% of all cells do |

**Scope of the claim (v1.1 §2).** This bound holds for **adult (26–52 y) iliac
aspirate bone marrow mononuclear cells, unsorted, 10x 3′ v2, from the HCA
Census of Immune Cells**. It is not a statement that early NK stages are absent
from human bone marrow. The admissible wording is: *in bone marrow scRNA-seq
data of the conditions above, early NK stages are detected at no more than X*.
Age range outside 26–52 y, other skeletal sites, CD45⁺ or other enrichment, and
other chemistries are each capable of changing the bound and none of them is
covered here. GSE233304's marrow gives X ≤ 0.112%, consistent but underpowered
and drawn from calvarial bone.

**What the source actually says.** Checked against the full text of Melsen 2022
(Front Immunol 13:1044398), not the abstract: its own scRNA-seq is **one
donor** (6,525 cells); QC removed every cell with **fewer than 1,000 expressed
genes**; and the "no developmental connection between CD34⁺ progenitors and NK"
claim rests on integration with the HCA marrow dataset **as downloaded through
HCAData, i.e. the filtered 378,000-cell matrix**, so the primary source never
examined an unfiltered droplet matrix. **It reports no frequency bound.** The
firm product of this preflight is therefore: the negative holds, the source
never bounded it, this run bounds it, and it is not a QC false negative — the
early-NK count is identical before and after QC, and those cells are *more*
likely to survive a 1,000-gene floor than an average marrow cell.

---

## 5. Failure modes produced by this run

Full write-ups in [`docs/TRANSFERABLE.md`](docs/TRANSFERABLE.md). None is
specific to NK or to marrow.

**T1 — a gating rule that is itself compartment-dependent.** The inherited rule
"an NK call requires zero CD3/TCR counts" dropped 46.5% of NK-cluster cells in
marrow and 73.7% in blood, because blood soup is T-dominated and blood NK are
sequenced deeper. Any gate of the form *define a population by a gene being
zero* has this problem whenever it is applied across two compartments with
different ambient levels: the gate's stringency is set by the nuisance
parameter under study. The rule adopted instead is **depth-symmetric**: keep a
cell if its own lineage panel outscores the competing panel on log1p(CP10K).
Both quantities are measured in the same cell, so sequencing depth and the
ambient level cancel in the comparison, which an absolute zero-count test
cannot do. Applying it moved blood NK from 1.6% to 11.3% of mononuclear cells
and marrow NK to 2.8%, both in the range flow cytometry reports.

**T2 — a per-gene ambient statistic on a low-count gene is a doublet detector.**
The first acceptance test failed on marrow (MS4A1 ambient share 0.24) not
because the soup model was wrong but because 5 of 8 MS4A1 counts sat in one
cell. Doublets are now removed before any gate is formed, and acceptance is
judged on the spread of the per-gene ratios around ρ over genes whose
*expected* ambient counts clear a floor.

**T3 — the detection-rate criterion is void without an ambient floor, and this
is now a revision to criterion ②.** In GSE120221, HBB is detected in 100% of
cells, GNLY in 52% and NKG7 in 51%, erythroid and myeloid included: the 5%
floor is cleared by soup alone. The revised rule, to be applied from here on:

> Criterion ② (detection rate inside 5%–85%) is not a criterion on its own. It
> holds only when reported together with the same gene's detection rate in a
> population that certainly does not express it, measured in the same
> libraries; the admissible quantity is the gap, not the rate. Where no
> unfiltered matrix exists, the floor cannot be estimated and the criterion
> cannot be applied at all.

Earlier conclusions that used the bare criterion need re-examination on this
point. Every detection rate quoted above carries its floor.

---

## 6. Datasets

`manifest.yaml` carries all twelve, each with its admission decision and the
evidence for it. Four of the five accessions named in the instruction fail the
hard gate: GSE120221 (2,994 barcodes per sample), GSE185381 (27,170),
GSE130430 (~1,000–2,500) and Melsen's own GSE199411 (7,000) all deposit called
cells only. GSE231946 — same-subject BM and PB for eight subjects including two
healthy controls, the best design found anywhere — fails the same gate. Three
of the instruction's descriptions were wrong and are corrected in the manifest,
including that GSE130430 does contain two same-individual BM/PB pairs.

---

## 7. E1 — adaptive-NK BM–PB magnitude on the paired dataset (spec v1.2)

**Falls in row 2 of the E1 table: the public data give no stable magnitude, so
the distribution is reported and no point estimate is offered.** Six of eight
candidate genes clear the T3 gate (KLRC2, B3GAT1, CX3CR1 high; FCER1G, KLRC1,
GZMK low); ZEB2 is excluded because its ambient floor is 0.547 against a
detection of 0.848 (1.55×, below the required 3×) and IL7R because its
detection of 0.076 sits on a floor of 0.063. The four per-donor Δ (marrow −
blood, adaptive share of NK, depth-symmetric score) are **+0.010 (P1), −0.110
(P15), −0.081 (P16), −0.010 (c9)** — median −0.046, IQR [−0.088, −0.005], but
only three of four share a sign, and the presence/absence sensitivity arm
returns all four positive, so the two definitions disagree even in direction.
GSE233304 is calvarial marrow from glioblastoma patients plus one
non-neoplastic intracranial control, n = 4, and this is a magnitude reading on
public data only. Full gate and per-donor table:
`results/e1_adaptive_delta_gse233304.tsv`.

## 8. E2 — donor age structure of the admitted datasets (spec v1.2)

**Falls in row 2 of the E2 table — the middle band exists — with row 3 applying
alongside it.** Across the admitted human bone-marrow datasets there are 46
donors, of whom **28 (61%) carry a deposited age** and 18 do not; GSE233304 and
GSE181543 deposit none at all. Ages run 24–84 y, and **17 donors fall in
25–55 y** (GSE120221 n = 9, HCA Census of Immune Cells n = 8), so a continuous
marrow curve across the middle band is in principle drawable from public data
and the question of why it has not been drawn is a separate one. The structural
absence is at the young end, not the middle: **0 donors under 18** and 1 under
25 in any admitted marrow dataset. This is an inventory of public data.
Per-donor table: `results/e2_age_distribution.tsv`; decade counts:
`results/e2_age_summary.tsv`.
