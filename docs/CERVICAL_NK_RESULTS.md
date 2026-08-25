# NK proportion in cervical tissue — results

Full per-donor table: `results/cervical/NK_TABLE.csv`
Primary paired result: `results/cervical/NK_PAIRED.csv`
Method and marker validation: `docs/CERVICAL_NK_METHODS.md`

37 libraries, 4 datasets, 4,374 gated NK cells. No author annotation used.

---

## Primary result — three same-donor pairs (n=3)

This is the only comparison that holds donor, batch and dissociation protocol
constant. Everything else in this document is context.

| donor | normal | tumour | NK/lymph normal | NK/lymph tumour | fold | Fisher p |
|---|---|---|---|---|---|---|
| P45 | N1 · 15/896 | T1 · 6/1850 | 1.67% (1.02–2.74) | 0.32% (0.15–0.71) | 0.19 | **0.0005** |
| P50 | N2 · 23/722 | T2 · 21/828 | 3.19% (2.13–4.73) | 2.54% (1.66–3.85) | 0.80 | 0.45 |
| P51 | N3 · 1/54 | T3 · 2/335 | 1.85% (0.33–9.77) | 0.60% (0.16–2.15) | 0.32 | 0.36 |

All three point down. One reaches significance. The sign test over three pairs
does not (p=0.25), and Mantel–Haenszel across strata gives OR 1.96.

**n=3, one adequately powered.** P51 rests on 1 NK cell against 2 and is
uninformative at any gating threshold. This is a direction, not an effect size.

## Cervical cancer, all 15 tumour libraries

| denominator | median | range |
|---|---|---|
| of lymphocytes | **2.76%** | 0.27 – 12.09% |
| of CD45+ | **1.99%** | 0.19 – 7.12% |
| of all cells | **1.22%** | 0.02 – 4.66% |
| CD56bright (of NK) | 21.5% | 0 – 50.4% |

1,172 NK cells across 15 libraries from 3 datasets. These are pooled across
datasets **only to state a range**, not to compare against anything.

## Why the cross-donor comparison is dead

Between-donor spread in **healthy cervix alone** is 60-fold: 0.38% to 23.23% of
lymphocytes across 8 donors (GSE173231). That exceeds every between-group
difference in the study. No normal-vs-tumour, normal-vs-HSIL, or cross-dataset
statement survives it.

The spread is *not* measurement noise. Two donors contributed two independent
biopsies each, and the replicates agree closely:

| donor | biopsy A | biopsy B |
|---|---|---|
| CX6 | 23.23% | 23.13% |
| CX7 | 8.79% | 12.23% |

So the measurement is reproducible; the variance is real between-donor and
between-site variation — plausibly ectocervix vs endocervix vs transformation
zone, which none of these submissions record. Reproducible variance is still
fatal to a cross-donor comparison; it just means the problem is sampling, not
instrumentation.

Within-dataset the same problem recurs. In E-MTAB-12305 the immune fraction
ranges from 1.3% (N3, effectively a stromal biopsy: 4,321 muscle and 3,332
fibroblast cells against 54 T/NK) to 85.5% (H1). Composition is dominated by
what the biopsy caught and how it dissociated.

## External calibration

Converted to the CD45+ denominator, the healthy-donor median is **2.97%**
(donor-level, n=8) against a flow-cytometry expectation of ~2.7% for ectocervix.
The gate lands on the independent expectation, which is the main evidence that
it is calibrated rather than merely self-consistent.

## Numbers that should not be trusted, and why

**Rows resting on <20 NK cells** — the percentage is dominated by Poisson noise:

| library | NK | issue |
|---|---|---|
| N3 (E-MTAB) | 1 | only 54 lymphocytes in the whole library; 1.85% is one cell |
| T3 (E-MTAB) | 2 | 9,936 of 10,401 cells are epithelial; 4.4% immune |
| ADC_3 (GSE197461) | 1 | 0.27% is one cell |
| T1 (E-MTAB) | 6 | drives the one significant pair; estimate rests on 6 cells |
| CX1, CX5 (GSE173231) | 5, 11 | |
| ADC_2 (GSE197461) | 12 | |
| N1, CX4 (E-MTAB, GSE173231) | 15, 15 | |

**The HSIL numbers look elevated and are not interpretable.** E-MTAB HSIL is
16–22% of lymphocytes against 1.7–3.2% for that dataset's adjacent-normal, and
GSE208653 HSIL is 5.2–7.5% against 1.1–10.2% for its normals. But every HSIL
library comes from a *different donor* than every normal library, so this
comparison sits inside the 60-fold between-donor spread. It is not evidence of
NK enrichment in HSIL.

**Cross-dataset absolute values are not comparable at all**, per the dissociation
caveat below. GSE208653 tumours median 5.96% and E-MTAB tumours median 1.57% of
lymphocytes is a protocol difference, not a biological one.

---

## Required limitations

### 1. The capture confound points the same way as the expected answer

The positive half of the gate requires *detecting* NKG7/GNLY/NCAM1/KLRF1.
Detection is depth-limited, so NK cells carrying less endogenous RNA are
preferentially missed. If tumour-side NK carry ~25% less RNA, the gate loses
tumour NK first, and "NK fall in cancer" is partly manufactured by the
measurement. **The bias direction is identical to the reported direction**, so
the paired result above cannot be taken at face value.

Measured here (`results/cervical/NK_CAPTURE.csv`), as the NK/T median-UMI ratio
within each library, which absorbs that library's own depth:

| tissue | NK/T UMI |
|---|---|
| normal (healthy) | 1.01 |
| tumour | 1.15 |
| normal (adjacent) | 1.22 |
| HSIL | 1.23 |
| normal (GSE208653) | 1.31 |

**The ~25% tumour-side deficit is not reproduced in these data.** NK carry
slightly *more* RNA than the T cells beside them in every tissue, and tumour
(1.15) sits inside the normal range (1.01–1.31).

This does **not** clear the bias, and the reason matters: this ratio is measured
only on NK cells that survived QC and passed a detection-based gate. Cells lost
to low RNA are by construction absent from it. The measurement is conditioned on
survival and can only ever underestimate the effect. It says the surviving
tumour NK are not RNA-poor; it cannot say how many were never counted.

### 2. Dissociation dominates composition

Single-cell proportions are set largely by dissociation protocol and by whether
a library was immune-enriched. Cross-dataset absolute values must not be
compared. Only within-dataset, and ideally within-donor, comparisons carry
information — which is why the paired n=3 is the primary result despite being
the smallest number in this document.

### 3. Ambient correction is a model, not a measurement

decontX was validated (foreign keratins driven to ~0 in immune cells while PTPRC
holds) but it is still an inference. Ambient is highest for exactly the NK marker
genes: uncorrected GNLY/NKG7 soup inflates NK-like signal in non-NK cells, while
ambient CD3 in true NK pushes them past a negative gate. Both are corrected here;
neither is eliminated.

### 4. The best-graded dataset is missing

Sci Adv `add8977` / PMC10277926 — 5 normal, 4 HSIL, 5 microinvasive, 6 invasive,
the only properly graded series — is deposited in GSA-Human (`subHRA005765`),
controlled access, requiring an approved application. Its absence is the main
reason the HSIL question stays open.

### 5. CD3E is expressed by NK cells in these data

Established by discriminant test, not assumption: ambient is refuted by sign
(CD3E detection *falls* with a sample's T fraction, rho −0.90 / −0.71 / −0.36),
doublets by a flat doublet score and a CD3-module lift of only 1.17–1.36.
Anyone applying a CD3E-negative or TRBC2-negative definition of NK to cervical
scRNA-seq will discard most real NK cells. See `docs/CERVICAL_NK_METHODS.md`.
