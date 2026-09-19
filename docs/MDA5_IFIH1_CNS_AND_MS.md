# MDA5 (IFIH1) in the human nervous system, and what changes in multiple sclerosis

Date: 2026-09-19. MDA5 is encoded by **IFIH1** (aliases MDA-5, Helicard, IDDM19);
it is the cytoplasmic long-dsRNA sensor that signals through MAVS to induce type I
interferon. Its paralogues are `DDX58` (RIG-I) and `DHX58` (LGP2). IFIH1 is itself
an interferon-stimulated gene, so "baseline expression" and "inducible response"
are reported separately throughout.

All numbers below were computed here from primary data. Tables are in
[`results/mda5/`](../results/mda5/); the GSE180759 run log is
[`log_GSE180759_analysis.txt`](../results/mda5/log_GSE180759_analysis.txt).

---

## Answer 1 — MDA5 is a vascular and myeloid gene in the CNS, not a neuronal one

Four independent resources agree on the same ranking, and the spread is about
**60-fold from top to bottom**.

| cell class | Census CP10k | HPA nCPM | relative to neurons |
|---|---|---|---|
| **endothelial cell** | **1.095** | **118.9** | **63×** |
| microglia / CNS macrophage | 0.338 | 35.5 | 19.5× |
| fibroblast | 0.230 | 23.2 | 13.3× |
| pericyte / smooth muscle / leptomeningeal | 0.219 | 14.1 | 12.6× |
| lymphocyte / leukocyte | 0.141 | 30.2 | 8.1× |
| ependymal / choroid plexus | 0.098 | 11.2 | 5.6× |
| astrocyte | 0.086 | 8.5 | 5.0× |
| OPC | 0.084 | 9.0 | 4.9× |
| oligodendrocyte | 0.061 | 7.5 | 3.5× |
| **neuron** | **0.017** | **0.3** | **1.0×** |

**Two bulk, antibody-free technologies give the same ranking**, which matters
because single-nucleus data under-detect. Immunopanned human brain cells (Zhang
2016, GSE73721, FPKM): endothelial 6.5, myeloid/microglia 5.13, astrocytes 1.70,
oligodendrocytes 1.18, neuron 0.8. Immunopanned mouse cells (GSE52564, FPKM):
endothelial 2.06, microglia 0.54, neuron 0.42, OPC 0.22, astrocyte 0.21, with both
oligodendrocyte fractions at the 0.1 detection floor. Astrocyte IFIH1 is not fixed:
astrocytes purified from the core of human brain tumours reach 6.5–11.4 FPKM, four
to seven times normal cortical astrocytes.

Sources: **CELLxGENE Census** (2025-11-08 release), 18,551,076 human brain cells,
783 donors, 144 datasets, single-nucleus subset, CP10k = IFIH1 counts per 10,000
total counts. **Human Protein Atlas** single-nuclei brain panel (Siletti atlas),
nCPM. Full tables: [`out_census_baseline_nucleus.csv`](../results/mda5/out_census_baseline_nucleus.csv),
[`out_hpa_brain_ifih1.csv`](../results/mda5/out_hpa_brain_ifih1.csv),
[`out_baseline_consolidated.csv`](../results/mda5/out_baseline_consolidated.csv).

**The endothelial result is not an artefact of pooling.** Computed per dataset,
endothelium exceeds microglia, astrocytes, oligodendrocytes and OPCs in **14 of 15**
datasets that measure all five, and exceeds microglia in **15 of 17**. Median
across datasets: endothelial 0.946, microglia 0.298, mural 0.170, astrocyte 0.100,
OPC 0.058, oligodendrocyte 0.055, neuron 0.006
([`out_census_per_dataset.csv`](../results/mda5/out_census_per_dataset.csv)).

**Neurons are the floor, consistently.** In the HPA panel every neuronal subtype
sits below 5 nCPM and most below 1.5 (hippocampal dentate gyrus 0.0, medium spiny
neuron 0.1, Purkinje-adjacent populations 0.1–0.2), against 118.9 for endothelium.
Within Census, granule cells (0.003) and Purkinje cells (0.005) are the lowest
measured populations in the brain.

### Confirmed in both MS cohorts' own control tissue

| cell type | GSE180759 controls, CP10k | GSE279180 controls, CP10k |
|---|---|---|
| lymphocyte / T cell | 0.352 | 0.785 |
| endothelial / vascular | 0.069 (n=65 nuclei) | 0.560 |
| microglia | 0.118 | 0.283 |
| astrocyte | 0.101 | 0.103 |
| OPC | 0.020 | 0.052 |
| oligodendrocyte | 0.029 | 0.047 |
| neuron | 0.000 (n=34 nuclei) | 0.017 |

### Tissue level, for orientation

GTEx v8 median TPM: IFIH1 is a **peripheral-immune** gene first (EBV-transformed
lymphocytes 87.0, spleen 15.9, lung 13.4) and the brain is among its lowest
tissues. Within the nervous system, **peripheral nerve is ~5× any brain region**
(tibial nerve 11.7, spinal cord 4.3, substantia nigra 3.0, cortex 2.2, cerebellum
1.0) — consistent with the cell-level result, since nerve is rich in endothelium,
fibroblasts and resident immune cells while brain is dominated by neurons and
oligodendrocytes ([`out_gtex_ifih1_tissues.csv`](../results/mda5/out_gtex_ifih1_tissues.csv)).

**So: MDA5's CNS expression tracks the vascular/immune compartment and the barrier
interfaces. The brain parenchyma proper — neurons above all, then oligodendrocytes,
OPCs and astrocytes — is where MDA5 is scarcest.**

---

## Answer 2 — in MS, MDA5 rises in microglia/macrophages and in oligodendrocytes

Three human cohorts were used. GSE180759 (Absinta 2021) is the dataset in
question; the other two were added because GSE180759 has only **three** control
brains, which is too few to settle a donor-level question on its own.

| cohort | accession | design | donors |
|---|---|---|---|
| Absinta 2021 *Nature* | GSE180759 | snRNA-seq, MRI-staged white-matter lesions | 5 MS + 3 control |
| Lerma-Martin 2024 *Nat Neurosci* | GSE279180 | snRNA-seq, subcortical MS lesions (CA/CI) | 7 MS + 6 control |
| Jäkel 2019 *Nature* | GSE118257 (via Census) | MS white matter | 4 MS + 5 control |

### Combined donor-level meta-analysis (21 brains)

Each donor contributes one depth-normalised CP10k per cell type; values are scaled
by their own cohort's control median, then pooled. The donor, not the nucleus, is
the unit — cell-level p-values here are pseudoreplicated and not used for the verdict.

| cell type | control donors | MS donors | MS / control | Mann-Whitney p |
|---|---|---|---|---|
| **microglia / macrophage** | 8 | 12 | **2.07×** | **0.0061** |
| **oligodendrocyte** | 9 | 12 | **2.42×** | **0.0077** |
| astrocyte | 9 | 12 | 0.87× | 0.37 |
| OPC | 9 | 12 | 0.82× | 0.80 |
| endothelial / vascular | 5 | 12 | 1.32× | 0.51 |

[`out_meta_donor_result.csv`](../results/mda5/out_meta_donor_result.csv),
per-donor values in [`out_meta_donor_cp10k.csv`](../results/mda5/out_meta_donor_cp10k.csv).

**Astrocytes and OPCs do not change. Neither does endothelium** — despite carrying
by far the highest baseline. The change is confined to the myeloid and
oligodendroglial compartments.

A second, more conservative model agrees on the effect size and is less sure of
the p-value. Negative-binomial GLM on donor pseudobulk counts, `IFIH1 ~ condition
+ cohort` with `log(total counts)` as offset
([`out_meta_nb_glm.csv`](../results/mda5/out_meta_nb_glm.csv)):

| cell type | fold change | 95% CI | p |
|---|---|---|---|
| microglia / macrophage | 1.99 | 0.97 – 4.11 | 0.062 |
| oligodendrocyte | 1.97 | 1.05 – 3.70 | 0.034 |
| astrocyte | 0.78 | 0.39 – 1.55 | 0.47 |
| OPC | 0.97 | 0.43 – 2.18 | 0.95 |
| endothelial | 0.97 | 0.42 – 2.26 | 0.95 |

**The fold-change estimate is robust and the p-value is not.** Varying the
dispersion parameter from 0.1 to 2.0 moves microglia only between 1.79× and 2.18×
and oligodendrocytes between 1.96× and 2.01×, while p ranges from 0.004 to 0.26.
The defensible claim is **a roughly two-fold increase in both cell types, with
21 brains giving borderline significance** — not a precisely quantified p.

### How stable are the two positive results?

Re-run under progressively stricter donor-inclusion rules:

| threshold | microglia Fisher p | oligodendrocyte Fisher p |
|---|---|---|
| n ≥ 30 nuclei (primary) | 0.018 | 0.007 |
| n ≥ 100 | 0.061 | 0.007 |
| n ≥ 100 and median depth ≥ 1000 | 0.050 | 0.011 |
| n ≥ 200 and median depth ≥ 1000 | (too few control donors) | 0.011 |

**The oligodendrocyte result is stable; the microglial result sits on the
significance boundary** and depends on retaining small, shallow control samples.
Per cohort, one-sided: microglia GSE180759 ratio 2.71 p=0.095 (only 2 usable
control donors), GSE279180 ratio 1.88 p=0.027. Oligodendrocyte GSE180759 ratio
2.25 p=0.125, GSE279180 ratio 2.51 p=0.007.

### Within MS, MDA5 tracks lesion activity and the disease-associated microglial states

Lerma-Martin microglia, depth-normalised CP10k by annotated subtype:

| microglial state | n nuclei | CP10k |
|---|---|---|
| **MG_CA** (chronic active) | 547 | **0.690** |
| **MG_Rim** (lesion rim) | 916 | **0.589** |
| MG_Dis (disease-associated) | 1281 | 0.508 |
| MG_Homeo3 | 1243 | 0.419 |
| MG_Phago1 | 1008 | 0.357 |
| MG_Homeo1 | 1825 | 0.353 |
| MG_Homeo2 | 2077 | 0.330 |

By lesion type: chronic active 0.467, chronic inactive 0.458, control 0.283
([`out_lerma_microglia_subtype.csv`](../results/mda5/out_lerma_microglia_subtype.csv)).

In GSE180759 the same gradient appears across MRI-staged regions — immune-cell
IFIH1 CP10k: control WM 0.118 → periplaque 0.253 → chronic inactive edge 0.275 →
chronic active edge 0.282 → lesion core 0.386.

**This places MDA5 on the chronic-active rim**, which is the compartment Absinta
2021 identified as the lymphocyte–microglia–astrocyte axis driving lesion
expansion.

### The rise is inside genuine microglia, not a macrophage influx

MS recruits CD163⁺/MRC1⁺ macrophages, so a shift in myeloid composition is the
obvious confounder. Restricting GSE180759 to **P2RY12⁺ nuclei** (homeostatic
microglia, a marker macrophages lack):

| population | control | MS | odds ratio | p |
|---|---|---|---|---|
| P2RY12⁺ microglia only | 0.147 CP10k (1.49%) | 0.375 CP10k (8.02%) | 5.75 | 0.0038 |
| CD163⁺/MRC1⁺ macrophages | 0.183 CP10k (4.00%) | 0.408 CP10k (8.81%) | 2.32 | 0.72 |
| all immune nuclei | 0.118 CP10k (1.38%) | 0.290 CP10k (6.11%) | 4.63 | 0.00024 |

Baseline IFIH1 is similar in the two myeloid populations (microglia 0.354 vs
macrophage 0.405 CP10k), so substituting one for the other cannot manufacture the
effect. The induction is cell-intrinsic to microglia.

### Is it MDA5 specifically, or the whole interferon module?

Judged on raw detection rates, Lerma-Martin looks like a blanket interferon
response — IFIH1 OR 4.13, but DDX58 2.69, STAT1 2.77, MX1 2.19, OAS1 3.19, BST2
6.42, ISG15 14.13 as well. **That reading is an artefact of the 3.84× depth
difference.** Re-run depth-normalised and at donor level, the picture is much more
specific ([`out_lerma_isg_panel_donor.csv`](../results/mda5/out_lerma_isg_panel_donor.csv)):

| gene | oligodendrocyte | microglia | astrocyte |
|---|---|---|---|
| **IFIH1** (MDA5) | **2.45× p=0.014** | 1.64× p=0.054 | 1.12× n.s. |
| DDX58 (RIG-I) | **1.42× p=0.022** | 1.12× n.s. | 1.11× n.s. |
| EIF2AK2 (PKR) | **1.34× p=0.035** | 0.99× n.s. | 1.01× n.s. |
| IRF7 | **4.96× p=0.012** | 1.32× n.s. | **1.93× p=0.008** |
| MAVS | 0.85× p=0.051 | **0.76× p=0.014** | **0.71× p=0.002** |
| BST2 | 5.64× n.s. (near zero) | **2.64× p=0.005** | 1.18× n.s. |
| STAT1, MX1, OAS1, ISG15, IFIT3, XAF1 | n.s. | n.s. | n.s. |

**In oligodendrocytes the cytoplasmic RNA-sensing arm moves as a unit** — MDA5,
RIG-I and PKR all rise together with IRF7 — while the downstream effectors (MX1,
OAS1, BST2, ISG15) stay at essentially zero (0.000–0.003 CP10k). Oligodendrocytes
in MS look **primed to sense dsRNA without executing an interferon response**.

In microglia, MDA5 and BST2 rise while the classic effectors do not. In
GSE180759's microglia the same selectivity appears from the other direction:
IFIH1 2.45× and IRF7 2.59× and BST2 1.49× up, while RIG-I (0.65×), STAT1 (0.75×),
MX1 (0.34×), OAS1 (0.27×) and IFIT3 (0.46×) fall. RIG-I staying flat while MDA5
rises in the same nuclei is a clean internal control against a depth artefact.

**`MAVS`, the adaptor MDA5 signals through, falls** in microglia (0.76×, p=0.014)
and astrocytes (0.71×, p=0.002) — so more sensor is being made while less of its
obligate adaptor is. That is worth following up rather than asserting.

`ADAR`, which edits endogenous dsRNA and is the brake on MDA5, does **not** rise
with it in GSE180759 microglia (1.342 → 1.023, 0.76×, n.s.), while ADAR2/`ADARB1`
rises 1.61×. Too weak to build on, but it is the natural next thing to test.

---

---

## Answer 3 — split by compartment: the MS effect is a white-matter effect

**No single public dataset with accessible processed data sampled both compartments.**
GSE180759, GSE279180 and GSE118257 are all subcortical **white matter**. The only
study designed for this question, Macnair 2025 (*Neuron*, 632,000 nuclei, 54 MS +
26 controls, WM lesion / normal-appearing WM / GM lesion / normal-appearing GM),
deposits its data at EGA under controlled access (EGAD00001009169) and could not
be used. So the grey-matter arm comes from a separate study:
**Schirmer 2019 (*Nature*)**, cortical snRNA-seq, 48,919 nuclei, 9 control and 12
MS donors, obtained from the UCSC Cell Browser. The compartment comparison is
therefore **cross-study**, and that is its main weakness.

### Baseline: the cell-type hierarchy is the same in grey and white matter

Absolute units differ between studies, so each cell type is expressed relative to
**astrocytes in its own dataset** ([`out_gm_vs_wm_baseline.csv`](../results/mda5/out_gm_vs_wm_baseline.csv)):

| cell type | cortex / GM | WM (Absinta) | WM (Lerma-Martin) |
|---|---|---|---|
| endothelial | 10.40 | 0.68 | 5.44 |
| microglia | 2.80 | 1.18 | 2.75 |
| astrocyte | 1.00 | 1.00 | 1.00 |
| oligodendrocyte | 0.70 | 0.29 | 0.46 |
| OPC | 0.29 | 0.20 | 0.51 |
| neuron | 0.16 | 0.00 | 0.17 |

Same ordering in both compartments: endothelium highest, then microglia, then
astrocytes, then oligodendroglia, with neurons lowest. **Compartment does not
change which cells express MDA5.** (The Absinta endothelial figure rests on 65
control nuclei and should be ignored.)

### Disease effect: depth-matched, same statistic in both compartments

[`out_gm_wm_depthmatched.csv`](../results/mda5/out_gm_wm_depthmatched.csv)

| cell type | WM Absinta | WM Lerma-Martin | cortex / GM (Schirmer) |
|---|---|---|---|
| **oligodendrocyte** | **OR 1.84, p=0.0019** | **OR 2.40, p=7×10⁻⁴²** | OR 1.14, p=0.77 |
| **microglia** | OR 2.55, p=0.17 | **OR 1.64, p=0.0017** | OR 1.78, p=0.54 |
| astrocyte | 0.49, n.s. | 1.14, n.s. | 0.96, n.s. |
| OPC | 2.01, n.s. | 0.72, n.s. | 2.16, p=0.13 |
| endothelial | 4.24, n.s. (n=54) | 0.83, n.s. | 1.92, p=0.16 (n=90) |
| neuron | — | 0.88, n.s. | **1.12, p=0.43** |

Three things follow.

**The oligodendrocyte effect is white-matter-specific, and this is not a power
problem.** Cortex contributes 3,070 oligodendrocytes per matched group and returns
OR 1.14 with p=0.77; the donor-level test agrees (1.56×, p=0.37). In white matter
the same statistic gives 1.84 and 2.40 in two independent cohorts. MDA5 goes up in
oligodendrocytes where the lesions are demyelinating white-matter tracts, not in
cortical oligodendrocytes.

**Neurons do not change, and grey matter is what finally allows that to be tested.**
White-matter datasets contain almost no neurons (34 control neuronal nuclei in
GSE180759). Cortex gives 10,134 per matched group, and the answer is a clean
negative: OR 1.12, p=0.43. Neurons have the lowest MDA5 in the brain and MS does
not raise it.

**Microglia point the same way in all three cohorts** (OR 2.55, 1.64, 1.78) but
cortex cannot confirm it: control cortex yields only **159** microglial nuclei,
and just 2 of 9 control donors reach 20. The microglial result stands on white
matter.

### A separate cortical observation, not about MDA5

Microglia are **0.91%** of control cortical nuclei and **4.84%** of MS cortical
nuclei, a roughly five-fold expansion, while neurons fall from 63.2% to 40.3%.
That is a composition change in MS cortex, and it is the reason the control
microglial baseline there is so thin.

### One nominal hit that does not survive inspection

Cortical OPCs scored 16.9× with p=0.041 at donor level. Five of nine control
donors have **zero** IFIH1-positive OPC nuclei, so the control median is 0 and the
ratio is division by roughly zero. Depth-matched, the same comparison gives OR
2.16, p=0.13. It is not a finding.

---

## Prior art: what was already published, and what these analyses add

A systematic literature and database sweep was run alongside the computation, with
every claim adversarially re-verified against primary sources. It changes the
standing of several results above.

### The oligodendrocyte result was already published — by Lerma-Martin themselves

Supplementary Table 6 of Lerma-Martin 2024 (pseudobulk DESeq2, per cell type)
already reports IFIH1 as significantly up in oligodendrocytes:

| contrast | cell type | log2FC | adjusted p |
|---|---|---|---|
| chronic active vs control | oligodendrocyte | +0.905 | 0.0122 |
| chronic inactive vs control | oligodendrocyte | +1.095 | 0.0047 |
| chronic inactive vs control | myeloid | +0.950 | 0.116 (not significant) |

IFIH1 does not appear for astrocytes, OPCs, neurons, endothelium, B cells or
stromal cells in any contrast. **So the oligodendrocyte finding is not new, and my
independent reanalysis reproduces the authors' own conclusion by a different
route.** Their myeloid result is likewise non-significant, which matches the
borderline status of the microglial claim here rather than contradicting it.

### Two bulk cohorts confirm it, and confirm the RIG-I selectivity

- **GSE138614** (Elkjaer 2019, 73 MS white-matter areas from 10 progressive-MS
  donors vs 25 control WM areas). IFIH1 up in active lesions (log2FC +0.558,
  FDR 0.021) and chronic active lesions (+0.478, FDR 0.0073); trend in
  remyelinating lesions; not significant in inactive lesions or normal-appearing
  WM. `DDX58` and `DHX58` are **never** significant in any lesion type. That is
  the same MDA5-over-RIG-I selectivity seen in GSE180759 microglia, now in bulk
  tissue and an independent cohort.
- **GSE283092** (Chen 2025, *Brain*, 96 dissected samples). IFIH1 is named in the
  paper among the most significant genes up in foamy active lesions. Paired
  against the same donor's normal-appearing WM, foamy active lesions give
  fold-change 1.79 (p=0.0015) and remyelinated lesions 1.36 (p=0.0009).

### What is genuinely unreported

- **IFIH1 appears nowhere in Absinta 2021's own published output** — not in any
  marker list, DE list, WGCNA module or drug table across all 13 supplementary
  tables. Neither do DDX58, MX1, ISG15, STAT1, IRF7, OAS1 or EIF2AK2. The
  GSE180759 analysis here is new for that dataset.
- **IFIH1 is absent from every Schirmer 2019 supplementary table**, having fallen
  below their detection filter. Their lesion-stage list does contain DDX58, DHX58,
  IFIT1, IFIT3, STAT1 and ISG15, all *declining* from early to late lesion. So the
  cortical null reported above is an original measurement, not a restatement.
- **No primary study has ever measured MDA5 in oligodendrocytes or OPCs**, in
  any species, by any method. Nor in primary or iPSC-derived human microglia; the
  "microglial MDA5" statements in reviews trace back to THP-1 and monocyte-derived
  macrophage experiments.

### The protein-level evidence is absent, and what exists points the other way

**There is no immunohistochemistry, immunofluorescence, western blot or proteomic
measurement of MDA5 in human MS brain tissue.** Everything above, and everything in
the literature, is transcript-level.

Worse for the naive reading, the available protein data in *normal* brain are
**discordant with the mRNA ranking**. The Human Protein Atlas IHC (antibody
HPA002656, reliability Approved) scores cerebral cortex as neuronal cells Medium,
endothelial cells Low, and **glial cells Not detected**; cerebellar Purkinje cells
Medium with granular and molecular layer cells not detected. A macaque SIV
encephalitis study (Co et al., *J Infect Dis* 2011) likewise found MDA5 protein
constitutively expressed in neurons and undetectable in perivascular macrophages
and astrocytes in uninfected brain, becoming induced in those glial compartments
only on infection.

Two readings are possible and the data do not separate them: either MDA5 protein
really is neuronal while the transcript is vascular/myeloid, or the single
antibody plus semi-quantitative four-level scoring misses thin glial processes
that are easy to miss and large neuronal somata that are easy to score. Either
way, **no transcript-level claim in this document should be restated as a claim
about MDA5 protein.**

### An independent parallel for the white-matter specificity

In that same macaque SIV model, MDA5 mRNA rose about 21-fold **in white matter** by
day 7, while **grey matter had higher basal levels and minimal infection-driven
change**. That is the same compartment asymmetry found here in MS: the inducible
compartment is white matter, and grey matter moves little.

### IFIH1 is not a multiple sclerosis susceptibility gene

This is worth stating because it is easy to assume otherwise from the T1D
literature. The IMSGC 2019 genomic map (47,429 cases, 68,374 controls) contains no
associated variant in the IFIH1 region. The GWAS Catalog has **zero** MS entries
for rs1990760 (A946T), which is genome-wide significant for type 1 diabetes,
psoriasis, SLE, autoimmune thyroid disease, inflammatory bowel disease and
selective IgA deficiency. Open Targets scores IFIH1–MS at 0.016 on literature
evidence alone, with a genetic-association score of exactly zero. The single
positive report (Martinez 2008, Spanish, n=412) was driven by rs2068330, which
lies in KCNH7 outside IFIH1, and failed replication in 591 French trio families.
The T1D-protective loss-of-function variant rs35667974 has no MS effect at 99%
power in 3,037 cases.

**So any MDA5 change in MS lesions is a downstream tissue response, not an
inherited driver of susceptibility.**

---

## Methods

- **GSE180759**: the authors' filtered matrix (29,432 genes × 66,432 nuclei) and
  their annotation were downloaded from GEO. Per-nucleus library sizes were computed
  by a full pass over the matrix (median 2,852 counts; 197.6M total UMI).
  Normalisation is CP10k. Cell types are the authors' own labels.
- **GSE279180**: the nine per-cell-type `.h5ad` files (103,780 nuclei, 13 donors),
  raw integer counts, chunked extraction of the same gene panel.
- **Census / Jäkel**: `cellxgene-census` 1.18.0, stable release, `raw` layer.
- **Schirmer 2019 (cortex / grey matter)**: the authors' normalised matrix and
  metadata from the UCSC Cell Browser (`cells.ucsc.edu/ms`), 48,919 nuclei, 21
  donors, four cortical regions. The matrix is per-cell library-size normalised and
  log-transformed rather than raw counts, so the reported statistics there are mean
  normalised expression and decile-depth-matched detection using the metadata's own
  per-nucleus UMI totals. Raw reads for that study are SRA PRJNA544731.
- **Statistics**: donor-level Mann-Whitney on depth-normalised CP10k is the primary
  test. Cell-level Fisher tests are reported only as descriptive, because nuclei
  within a donor are not independent. Cross-cohort combination is Fisher's method
  on one-sided per-cohort p-values.

## Caveats that materially limit these conclusions

1. **Depth is confounded with disease status in both cohorts.** In Lerma-Martin,
   MS microglia carry 3.84× the counts of control microglia (3,564 vs 927); in
   GSE180759, 1.48×. Raw detection rates are therefore badly inflated: the naive
   Lerma microglial effect is OR 4.13, but decile-depth-matched it is **OR 1.64**.
   Every headline number above is depth-normalised or depth-matched; any published
   detection-rate comparison on these data that is not will overstate the effect.
2. **GSE180759 has a complete platform confound.** All four periplaque libraries
   are NovaSeq 6000; all eleven lesion-edge libraries are HiSeq 2500. The apparent
   oligodendrocyte increase in that cohort sits entirely in the periplaque stratum,
   and a platform-matched cell-level test nulls it (OR 1.18, p=0.38). The
   oligodendrocyte conclusion therefore rests on Lerma-Martin.
3. **Three control brains, one with 10 immune nuclei.** GSE180759 can supply only
   two usable control donors for any myeloid test, which is why its microglial
   p-value never falls below 0.095 however strong the cell-level signal looks.
4. **Ambient RNA is not excluded.** Across GSE180759 libraries, microglial and
   oligodendrocyte IFIH1 detection correlate (Spearman ρ=0.63, p=0.009), which is
   equally consistent with shared soup and with a shared tissue-level interferon
   tone. Arguing against soup: microglial IFIH1 runs 2–9× the oligodendrocyte rate
   in the same library, and that ratio rises in MS.
5. **snRNA-seq under-detects.** IFIH1 is present in only 2.18% of GSE180759 nuclei
   (1,645 UMI total). Everything here is a low-count measurement.
6. **Jäkel adds direction but no power** — 78 control microglia, 0% positive vs
   3.14% in MS, p=0.23.
7. **The grey-matter versus white-matter comparison is cross-study**, because the
   one cohort that sampled both compartments (Macnair 2025) is under EGA controlled
   access. Cortex comes from Schirmer 2019 and white matter from three other
   studies, so a compartment difference is confounded with study, tissue handling
   and normalisation. The two conclusions that survive this are the ones resting on
   a *negative* in a well-powered cortical sample — oligodendrocytes (3,070 per
   matched group) and neurons (10,134) — since a cross-study artefact is a poor
   explanation for an effect that is strong in two independent WM cohorts and
   absent in an adequately sized GM one.

## What would settle it

- **Protein.** There is no MDA5 immunostaining of MS lesions anywhere in the
  literature, and the normal-brain antibody data point to neurons rather than glia.
  Until someone stains an MS lesion, every statement here is about transcript.
  This is now the single largest gap.
- A depth-balanced snRNA-seq design, or CITE-seq, in which control and lesion
  libraries are sequenced to matched depth.
- Testing the ADAR1 axis directly: if MDA5 induction on the rim reflects loss of
  the editing brake, A-to-I editing rates should fall in the same nuclei, which is
  measurable from the existing FASTQs (SRA SRP329679).
