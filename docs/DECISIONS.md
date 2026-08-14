# DECISIONS

Every judgement call: what was done, why, and what the alternative was.

---

## D1 — NETSKAR2024 contains no decidual sample. Closed as a dataset.

**Question (spec §9, cheap task 1):** does the Netskar 2024 pan-cancer NK
atlas (Zenodo 10.5281/zenodo.8434224, 89,216 NK cells) include decidua or
uterus?

**What was done.** The atlas `obs` table was read *without downloading the
matrix*: `all_nk_cells.h5ad` is 1.77 GB, and only its HDF5 `obs` group is
needed. A small HTTP range-request file object was handed to `h5py`, which
fetched **8.8 MB** of the 1,765 MB file.

**Answer — no.** Verified against three independent obs columns and 89,216
cells:

- `source` (14 levels): PBMC, brain_normal, breast_normal, breast_tumor,
  glioblastoma_tumor, lung_normal, lung_tumor, melanoma_tumor,
  pancreas_normal, pancreas_tumor, prostate_normal, prostate_tumor,
  sarcoma_tumor, skin_normal
- `tumor_type` (8 levels): breast, glioblastoma, lung, melanoma, none,
  pancreas, prostate, sarcoma
- `dataset` (47 studies), `sample` (61 levels) — no reproductive-tract entry

**Decision.** NETSKAR2024 is **not** promoted to a dataset. Recorded as: *the
principal NK-centric reference atlas contains no decidual or uterine sample*,
which is itself worth stating — the reference map that a decidual NK study
would naturally be projected onto has never seen this tissue.

**Alternative considered.** Downloading the per-tissue h5ads to search for a
mislabelled reproductive sample. Rejected: three orthogonal metadata columns
agree, and the 47 contributing studies are all cancer/normal-tissue cohorts.

---

## D2 — VT2018 has 7 donors, not 11. Everything downstream is bounded by this.

**Question (spec §9, cheap task 2):** verify the local E-MTAB-6701 against
the remembered "~70k cells / dNK 11,927 / 11 usable donors".

**What was done.** The dataset was **not** on local disk — this is a fresh
container, and `data/` is gitignored. Re-fetched from BioStudies
(`scripts_fetch_vt2018.sh`). `meta_10x.txt` (64,734 annotated cells) was
tabulated before any expression value was read.

**Measured:**

| quantity | recalled | measured |
|---|---|---|
| total cells | ~70,000 | **64,734** |
| dNK1+dNK2+dNK3 | 11,927 | **11,202** (11,932 including dNKp) |
| usable donors | 11 | **7** (D6–D12) |

Compartments: Decidua 36,186 / Placenta 18,547 / Blood 10,001.

Decidual dNK by donor:

| subset | D6 | D7 | D8 | D9 | D10 | D12 |
|---|---|---|---|---|---|---|
| dNK1 | 1388 | 108 | 1518 | 161 | 482 | 141 |
| dNK2 | 634 | 99 | 554 | 3409 | 326 | 386 |
| dNK3 | 1227 | 75 | 173 | 446 | 37 | **29** |
| dNKp | 49 | 24 | 263 | 216 | 102 | 73 |

Donors with decidua: D6 D7 D8 D9 D10 D12 (**6**).
Donors with blood: D6 D7 D8 D9 (**4**). D11 is placenta-only.

**Reading of the discrepancy.** The recalled cell count is close to
dNK1-3 + dNKp (11,932 vs 11,927) rather than to dNK1-3, so the remembered
figure most likely included dNKp. The 11-donor figure spans the companion
**Smart-seq2** experiment E-MTAB-6678; E-MTAB-6701 is the 10x experiment
alone. **R9 bars merging the two matrices**, so the extra donors are not
recoverable by combining them.

**Decision.** Proceed on measured values. The consequence is recorded in the
pre-registration *before* any test was run: analysis A is capped at n = 6 and
analysis B is pinned at n = 4, where a two-sided signed-rank test cannot
return p < 0.05 at all. See `docs/PREREGISTRATION.md` §2.

**Alternative considered.** Adding E-MTAB-6678 for donor count. Rejected
under R9 — different platform, and merging is exactly the confound the design
excludes.

---

## D3 — Analysis B against CD56-bright blood NK is `unmeasurable`, not negative.

Blood `NK CD16-` (CD56-bright) counts: D6 = 20, D7 = 3, D8 = 72, D9 = 125.
Only two donors clear the 30-cell threshold.

dNK is predominantly CD56-bright, so the biologically correct blood
comparator is CD56-bright — but **§6.3 forbids merging bright with dim**,
since that would fold an NK-subtype difference into what is reported as a
tissue difference.

**Decision.** Run and report dNK vs CD56-dim (n = 4) with its caveats, and
record dNK vs CD56-bright as **`unmeasurable` at n = 2** rather than
quietly substituting the dim population or pooling the two.

---

## D4 — Library id is taken from the barcode prefix.

`meta_10x.txt` publishes no per-cell library/batch column, but R6 requires
retention to be reportable against some library grouping. Cell barcodes carry
an `FCA*` run prefix, which is the 10x run. Ingest uses that prefix as the
library surrogate **only when** the manifest supplies no library column, and
only when more than one prefix exists.

**Alternative considered.** Reconstructing libraries from the SDRF's
assay-to-donor table. Kept in reserve; the prefix already separates runs and
needs no join that could silently mis-key.

---

## D5 — Ruler A is estimated leave-one-out.

With no empty droplets in a processed matrix, the ambient profile is
estimated as the abundance-weighted average transcriptome of the compartment,
and the per-cell ambient load `rho` from lineage-foreign genes.

The estimate of `rho` for gene *g* **excludes gene *g***. Without that, every
ambient control returns a soup fraction of exactly 1.0 by construction and
the ruler is calibrating against itself — the spread over held-out ambient
genes (the thing that produced the prior 0.402–1.000 range) would collapse
to a point and the ceiling would be meaningless.

---

## D6 — `unmeasurable` on ruler B keys on the reference lineage's absolute expression.

Ruler B asks whether NK signal is explainable as pickup from the reference
lineage. If the reference barely expresses the gene, no pickup is possible
*and* the ratio is uninformative — so the verdict must be `unmeasurable`.

Keying `unmeasurable` on the **ratio** instead would be circular: a gene
genuinely expressed by NK produces a high ratio, which would then be
mislabelled unmeasurable. So the gate is on the reference's own CPM
(>= 10) and detection (>= 5%), not on the ratio.

---

## D7 — dNKp is reported but never tested against dNK1/2/3.

Per §6.2: proliferating cells differ in RNA content and transcriptome
structure, so a detection-rate contrast against them is confounded by cell
state rather than by chemokine biology. dNKp detection rates appear in T1;
no T2 row tests them against the trio.

---

## D8 — The fetch failure was mawk, not the endpoint.

Worth recording because the first diagnosis was wrong and the wrong fix would
have been plausible. The parallel fetcher produced parts of wildly incorrect
size (one 3.3 GB part for a 515 MB range), which looked like the server
ignoring `Range` under concurrency.

The actual cause: the size was parsed with
`awk 'BEGIN{IGNORECASE=1}/^content-length/...'`. The system awk is **mawk**,
which accepts `IGNORECASE` and then ignores it, so the header `Content-Length`
never matched, `TOTAL` was empty, every computed byte range was nonsense —
and curl still exited 0. Range requests were working the whole time.

**Decision.** Parse headers with `grep -i`, verify each part's size against
the range requested, and verify the concatenated total. The size checks, not
the diagnosis, are what make the script safe against this class of failure.

---

## D9 — Panel Ensembl IDs come from Ensembl REST, not from the dataset.

Resolving panel symbols against the dataset's own gene list would guarantee a
100% match rate and make R11 vacuous. IDs were resolved against Ensembl REST
(GRCh38) independently, so the match rate measures something.

`CCL3L1` resolves to **three** Ensembl gene IDs (ENSG00000276085,
ENSG00000277336, ENSG00000277768). That is not an annotation nuisance to be
cleaned up — it is precisely the multi-mapping hazard §4.1 names. The primary
ID is recorded; if the gene misses in a dataset, that stands as a result and
must not be repaired by folding its reads into `CCL3`.

---

## D10 — The purity stop fired. Diagnosis: depth and doublets, not mislabelled T cells.

**What fired.** `qc_report.py` exits 3: raw triple-positive
`TRBC2⁺CD3E⁺CD3D⁺` is **17.35% in dNK1** and 16.92% in dNKp, above the
pre-registered 15% stop (`PREREGISTRATION.md` §6.1). dNK2 4.90%, dNK3 3.32%.

Per the stop, the main analysis was halted and `src/purity_diagnosis.py` was
run first. It reports; it re-labels nothing.

**Finding 1 — the raw metric is depth-confounded.** Requiring three genes at
≥1 count each is strongly depth-dependent, and step 4 of the pipeline runs
*before* depth matching. Median UMI: dNK1 4,159, dNK2 3,903, dNK3 3,563 —
dNK1 is the deepest of the trio. Re-measured at a common 2,348 UMI floor:

| subset | median UMI | raw triple-pos | depth-matched |
|---|---|---|---|
| T | 3,608 | 63.5% | 44.7% |
| **dNK1** | 4,159 | **17.4%** | **6.6%** |
| dNKp | 9,103 | 16.9% | 3.2% |
| dNK2 | 3,903 | 4.9% | 1.9% |
| dNK3 | 3,563 | 3.3% | 1.3% |
| Myeloid | 3,502 | 2.5% | 0.6% |
| Stromal | 12,784 | 0.08% | 0.01% |

Depth explains most of the level but not the ranking: dNK1 stays ~10× above
myeloid and ~500× above stromal after matching.

**Finding 2 — the flagged cells are not T cells.** In flagged dNK1 cells the
T-gate signal is substantial (median 8 counts; 79% carry ≥6), so this is not
one-count ambient pickup. But those same cells carry **more** NK signal than
the unflagged dNK1 cells (median 618 vs 407 counts over NKG7/KLRD1/GNLY/PRF1)
and higher total UMI (5,605 vs 3,873). A mislabelled T cell would be NK-low.
Both programmes at once, at elevated depth, is the **doublet** signature.
The same pattern holds in every subset with flagged cells, including Myeloid
(7,742 vs 3,448 UMI) and Stromal (15,852 vs 12,784).

**Finding 3 — not donor-specific.** dNK1 ranges 10.6%–34.3% across all six
donors, so it is not one bad library.

**Decision.** The published annotation is not wrong, so there is nothing to
"resolve" in the sense of re-labelling — and re-labelling is barred anyway.
The pre-specified remedy for contamination already exists in the plan
(§5 step 4: remove the flagged cells, re-run, keep both versions), so both
arms were run:

- `out/VT2018/` — all cells (primary)
- `out/VT2018/robustness_tpurged/` — triple-positive cells removed

Every conclusion in `RESULTS.md` holds in both. The headline CCL5 effect moves
from +23.3 to +22.5 pp with p unchanged at 0.03125.

**The stop is recorded as FIRED, not waived.** It is reported at the top of
`RESULTS.md`, and the differential doublet rate across the trio
(17.4% / 4.9% / 3.3%) is itself a caveat on any dNK1-versus-other contrast,
which is why the purged arm is a required companion rather than an appendix.

**Alternative considered.** Running a doublet caller (Scrublet, DoubletFinder)
to remove doublets properly. Rejected for this round: both depend on numba or
on re-clustering, which R10 and the no-re-clustering rule exclude. The
triple-positive purge is a blunter instrument that removes strictly more than
the doublets, which is the conservative direction.

---

## D11 — PAEP is absent from the published matrix.

The specification requires `PAEP` (glycodelin) as the decidua-specific
ambient ceiling gene, warning that without it the ceiling is underestimated.
It is **not in E-MTAB-6701's matrix at all**: no `gene_symbol` containing
"PAEP" and no `ENSG00000122133`, among 31,764 genes. This is the single panel
miss on the symbol route (match rate 98.5%, above the 90% bar).

**Decision.** Proceed with the remaining decidual ceiling genes — `DCN`,
`COL1A1`, `IGFBP1`, `HLA-G`, `CSH1`. All five return soup fraction **1.000**,
the correct value for a pure ambient gene, so the ceiling is anchored by five
independent tissue-specific transcripts rather than one. Recorded because the
specification's stated reason for adding PAEP (ceiling underestimation)
cannot be checked against PAEP itself here.

**Alternative considered.** Substituting another glycodelin-class transcript.
Rejected — the panel is frozen, and adding a gene after seeing the data is
exactly what the freeze exists to prevent.

---

## D12 — Ruler B has almost no power in decidua, and the conservative reading was kept.

Re-calibrated inside this dataset as the specification requires, ruler B's
pickup band runs to **1.253** in decidua (prior work: ≤ 0.118).

The reason is visible per gene. Ruler B asks whether NK signal is explainable
as pickup from the reference lineage (myeloid). That question only has an
answer for genes myeloid actually sources:

| ambient control | NK/myeloid ratio | myeloid CPM | compartment CPM |
|---|---|---|---|
| LYZ | 0.054 | 1,690 | 187 |
| C1QA | 0.044 | 3,962 | 387 |
| DCN | 1.078 | 54 | 3,290 |
| IGFBP1 | 1.136 | 79 | 1,510 |
| COL1A1 | 1.226 | 5 | 337 |
| IGKC | 1.253 | 16 | 23 |

Stroma- and trophoblast-sourced genes are picked up about equally by NK and
by myeloid cells, so their ratio sits near 1 and carries no information about
pickup. Restricted to controls the reference lineage actually sources
(myeloid CPM above the compartment average — a generic, data-driven test, not
a dataset branch), the band is **0.0538**, close to the prior 0.118.

**Decision.** Verdicts use the **wider** band (1.253), which is the
conservative choice: a wider band makes it *harder* for a gene to be called
positive. Both numbers are reported in `soup_ruler_calibration.tsv`
(`ruler_b_band_upper` and `ruler_b_band_upper_ref_sourced`). Switching to the
narrower band would only add positives, so no conclusion in `RESULTS.md`
depends on the choice.

---

## D13 — Differential cell loss is severe, and the sensitivity analysis is mandatory here.

R6 triggers. At the primary 2,137 UMI floor, retention within the tested trio:

| donor | dNK1 | dNK2 | dNK3 | gap |
|---|---|---|---|---|
| D6 | 0.878 | 0.907 | 0.960 | 8.3 pp |
| D7 | 0.880 | 1.000 | 1.000 | 12.0 pp |
| **D8** | 0.859 | 0.809 | **0.428** | **43.1 pp** |
| D9 | 1.000 | 0.999 | 0.998 | 0.2 pp |
| D10 | 0.998 | 1.000 | 1.000 | 0.2 pp |
| **D12** | 0.560 | 0.334 | — | **22.6 pp** |

Three donors exceed the 10 pp trigger, and D8 loses 57% of its dNK3 cells.
This is the exact failure mode the rule was written for — differential loss
whose direction can align with the conclusion.

**Decision.** `donor_tests.py` now repeats every comparison at all three
depth floors and writes `donor_level_tests_depth_sensitivity.tsv`. Result:
the dNK1-vs-dNK2 findings are stable (CCL5 +20.0/+23.3/+21.3 pp, p = 0.03125
at every floor; XCL1 +42.7/+39.3/+37.1 pp, p = 0.03125 at every floor), while
the dNK3 contrasts weaken at the deepest floor as D8's dNK3 falls below the
admission threshold. Reported in `RESULTS.md` §3.2.

---

## D14 — No gene reaches q ≤ 0.05, and that was known before the data.

6 of 27 primary-target genes sit exactly on the test floor in dNK1-vs-dNK2
(expected by chance: 0.84), 15 of 27 in dNK1-vs-dNK3, 10 of 27 in
dNK2-vs-dNK3. BH would need ≥ 17 of 27 simultaneously at the n = 6 floor, and
≥ 34 of 27 — an impossibility — at the n = 5 floor.

**Decision.** Report effect sizes, sign concordance, empirical-null z, and
dual-ruler verdicts; report `q_bh` truthfully as failing; and do **not**
retitle any result as "significant". The pre-registration fixed this reading
in advance (§2.2, §3 n-floor clause) precisely so that the outcome could not
be renegotiated once the effects turned out to be large.

The empirical null is what carries the weight instead: CCL5 at z = +5.6 and
XCL1 at z = +9.6 against 405 expression-matched genes, with the baseline
measured rather than assumed (+0.47 ± 4.05 pp — not zero, as §6.5 warned).

---

## D15 — n = 11 is unreachable in VT2018 by any route. The 10x donor parse is correct.

**Challenge.** A July metabolic result recorded `S_oxphos B−A −0.041, 11/11
donors negative, p = 0.001`. A two-sided sign test at n = 11 has a floor of
2×(1/2)¹¹ = 0.00098 ≈ 0.001, so that p value is the n = 11 floor exactly — it
implies 11 donors. Against this study's n = 6, that is impossible. Two
candidate explanations were put: a donor-column parse error here, or a
cross-platform merge there.

**Audit — the parse is not the problem.**

`E-MTAB-6701.sdrf.txt` has 30 rows (= 30 libraries) and exactly **one**
donor-like column, `Characteristics[individual]`, with **7** values (D6–D12).
There is no second identity column that could have been missed.

| donor | libraries | compartments | cells in matrix | runs in matrix |
|---|---|---|---|---|
| D6 | 5 | blood, decidua | 8,340 | 5 |
| D7 | 4 | blood, decidua | 6,820 | 4 |
| D8 | 4 | blood, decidua, placenta | 17,591 | 4 |
| D9 | 4 | blood, decidua, placenta | 13,327 | 4 |
| D10 | 5 | decidua, placenta | 10,430 | 4 |
| D11 | 2 | placenta | 2,502 | 1 |
| D12 | 6 | decidua, placenta | 5,724 | 3 |

Decisively: **libraries mapping to more than one donor = 0**, and **matrix
runs mapping to more than one donor = 0**. Neither collapsing several
libraries into one donor nor splitting one donor across runs is possible
here. 25 of the 30 libraries appear in the published matrix; no run in the
matrix is absent from the SDRF.

**Audit — a cross-platform merge cannot produce 11 either.**

`E-MTAB-6678` (Smart-seq2, decidua only) carries donors **D3, D5, D6, D7, D8,
D9** plus a `not available` group (4,136 rows).

| set | n | donors |
|---|---|---|
| 10x | 7 | D6–D12 |
| Smart-seq2 | 6 | D3, D5, D6–D9 |
| intersection | 4 | D6, D7, D8, D9 |
| **union** | **9** | D3, D5, D6–D12 |
| union, decidua only | 8 | D3, D5, D6–D10, D12 |

**Merging both platforms reaches 9 donors, not 11** (8 for decidua). So the
"it was a cross-platform analysis" explanation does not reconstruct 11 either.

No library-level count reproduces 11 as cleanly: decidual libraries = 15,
decidual runs in the matrix = 14, decidual runs containing any dNK1-3 = 12,
runs with ≥ 30 dNK1-3 cells = 8, donor × compartment units = 10.

**What this does and does not settle.** It settles that **if the July result
used E-MTAB-6701, its n = 11 is wrong** — under every grouping available in
that accession, alone or merged with its Smart-seq2 companion. It does not
identify where 11 came from, because **this repository contains no record of
that analysis**: a full-text search for `oxphos`, `S_oxphos`, `metabolic`,
`glycoly`, `n = 11` and `11 donors` across every `.md`, `.py` and `.csv` here
returns nothing. That work lives in another workspace and cannot be audited
from here.

**Recommended next check, and it is one lookup.** Print the distinct values of
the donor column used in the July analysis. `D1…D11` or anything containing
`D3`/`D5` means a cross-platform merge; `FCA…` values mean libraries were
counted as donors (an R1 violation, and the reading that best fits a
suspiciously round 11/11); values from another accession entirely mean the
conflict dissolves and both results stand.

**Addendum.** No library-level count listed above *equals* 11 — but 11 sits
between "decidual runs containing any dNK1-3" (12) and "runs with >= 30
dNK1-3 cells" (8), so a different cell-count threshold lands on it exactly.
The libraries-as-donors reading is therefore **not** excluded by the absence
of a literal 11 in that list, and it remains the most economical explanation
of a clean 11/11 with a p at the n=11 sign-test floor.

**Until then the OXPHOS result is neither confirmed nor withdrawn here, and
should not be cited.** What is established is the constraint it must satisfy:
no donor-level statistic from E-MTAB-6701 can have n > 7, or n > 6 restricted
to decidua.

---

## D16 — CORRECTION: the doublet explanation in D10 is not supported by measurement.

**D10 concluded** that the fired purity stop was "depth confounding plus NK–T
doublets". The doublet half of that was an inference from a pattern
(flagged cells carry more NK signal *and* more UMI), not a measurement. It
was challenged on exactly that ground, correctly: depth-driven ambient pickup
produces the same pattern.

**The two measurements that separate them** (`src/doublet_evidence.py`):

| prediction | doublets | depth-driven pickup |
|---|---|---|
| total UMI, flagged / unflagged | ≈ 2.0 (two cells in one droplet) | ≈ 1.0–1.4 (one deeper cell) |
| NK:T count-ratio distribution | bimodal (NK mode + doublet mode) | unimodal with a tail |

**Measured, decidua:**

| subset | flagged | UMI ratio | NK:T bimodality coef. |
|---|---|---|---|
| dNK1 | 659 | **1.45** | 0.332 |
| dNK2 | 265 | **1.31** | 0.236 |
| dNK3 | 66 | **1.67** | 0.466 |
| dNKp | 123 | **1.36** | 0.387 |
| T | 1,173 | 1.39 | 0.300 |

Every UMI ratio is 1.31–1.67, far from 2.0. Every bimodality coefficient is
below the 5/9 = 0.555 bimodality threshold, i.e. **unimodal**. The figure
(`out/VT2018/doublet_evidence.png`) shows the flagged distribution sitting
inside the unflagged one and shifted slightly right — not a second population.

**Corrected conclusion.** The flagged cells are **not** mislabelled T cells
(D10's first finding stands: they carry more NK signal, not less), and they
are **not** predominantly NK–T doublets either. The evidence supports
**depth-driven ambient pickup**: deeper cells are more likely to register
≥ 1 count on each of three genes. This is consistent with the raw 17.4% →
6.6% collapse under depth matching.

**What remains unexplained, and is not being explained away.** After depth
matching dNK1 is still 6.6% against myeloid 0.57% and stromal 0.01%. Depth
alone does not account for a dNK1-specific residual an order of magnitude
above other lineages. Candidate causes — annotation-boundary cells between
dNK1 and T, genuine low-level CD3D/CD3E transcription in NK, or dNK1 being
enriched in T-rich libraries — are **not** adjudicated by anything measured
here, and no further inference is offered.

The `modes` column in `doublet_evidence.tsv` is retained but is not evidence:
smoothed-histogram peak counting is noise-sensitive at these sample sizes
(66 flagged dNK3 cells yield 9 spurious "modes"). The bimodality coefficient
is the statistic to read.

**Consequence for the results: none.** Both arms were run and both are
retained. The v1.1 positive set is identical in the purged arm.

---

## D17 — PAEP verified absent in the source file, not lost at ingest.

D11 recorded PAEP as missing. That was checked against the ingested h5ad
only, so it could have been an ingest bug. Verified directly against the
4.1 GB published matrix:

```
grep -o -E "^(PAEP|PP14|CSH1|DCN|IGFBP1|GNLY)_ENSG[0-9]+" raw_data_10x.txt | sort -u
  CSH1_ENSG00000136488
  DCN_ENSG00000011465
  GNLY_ENSG00000115523
  IGFBP1_ENSG00000146678
```

PAEP is absent from the source; the four comparators are present. No alias
(`PP14`, `GdA`, `glycodelin`) and no `ENSG00000122133` appears anywhere in the
gene list.

The matrix carries 31,764 genes against the 33,694 of the CellRanger
GRCh38-1.2.0 reference — a filtered subset, missing 1,930 genes. It is not a
regional dropout: PAEP's chromosome-9 neighbours (GLIPR2, CNTFR, DCTN3, SIT1,
RMRP, KIF24, NPR2, SPAG8, HINT2) are all present.

**Consequence.** The decidual ambient ceiling is estimated without its single
largest expected contributor. Since ruler A returns 1.000 for five other
tissue-specific controls, the ceiling is anchored — but a genuinely higher
ambient plateau cannot be ruled out from this matrix. Recovering PAEP would
require re-quantifying from FASTQ, which is outside this round.

---

## D18 — Circularity test: the ranking of the eight positives is inverted by it.

**The problem.** dNK1/dNK2/dNK3 are Vento-Tormo's own **unsupervised
whole-transcriptome** clusters, characterised afterwards by differential
expression — verified against the published methods, so the contamination is
*partial*, not total. Nobody thresholded on CCL5. But the cluster boundary came
from a transcriptome containing CCL5, CXCR4 and XCL1, and those are among the
genes the characterisation named. Measuring them across the same labels partly
restates the clustering.

**Test** (`src/circularity_test.py`, authorised by v1.2 A4). Discard the
labels; re-group the same cells on CD39/CD103 equivalents (ENTPD1/ITGAE) —
markers no target gene informs — and re-run. CD160 and KLRB1 are barred as
gating genes, being dNK3 characterisation genes themselves.

**Gate composition** (this matters for reading the result):

| published | gate1 CD39⁺ | gate2 CD39⁻CD103⁻ | gate3 CD39⁻CD103⁺ |
|---|---|---|---|
| dNK1 | 1218 | 1965 | 155 |
| dNK2 | 99 | 4500 | 384 |
| dNK3 | 48 | 1483 | 291 |

Concordance 59.2%. **gate1 is 89.2% dNK1 — clean. gate3 is only 35.1% dNK3 —
badly diluted**, because single-gene dropout makes CD103⁺ a sparse gate in
scRNA where it is not in CyTOF protein.

**Result — effect retained versus the label-based effect:**

| gene | gate1→2 | gate1→3 | gate2→3 | median retained |
|---|---|---|---|---|
| XCL1 | +24.3 (4/4 donors) | +27.7 (4/4) | +5.5 | **0.62×** |
| XCL2 | +17.6 (4/4) | +22.9 (4/4) | +4.5 | **0.53×** |
| CCL5 | +17.1 (4/4) | +26.3 (4/4) | +6.8 | **0.46×** |
| CXCR4 | +8.6 (4/4) | +17.6 (4/4) | +8.7 (5/5) | **0.30×** |

n falls to 4–5, so `min_achievable_p` is 0.125 and **no p value here is
interpretable**. Direction and effect size are the readable quantities.

**Dilution is a competing explanation for shrinkage, and it is not sufficient.**
Impure gates shrink every effect. But dilution is gene-agnostic: it should
shrink all four genes by a similar factor. Observed shrinkage spans **0.62×
to 0.30×, a 2-fold spread**, ordered exactly as the circularity prior predicts
— XCL1 least affected, CXCR4 most.

**The within-label test settles it.** A gate contrast computed *inside a single
published cluster* cannot restate how that cluster was drawn, and is immune to
between-gate dilution:

| inside | contrast | gene | n | effect (pp) | donors concordant |
|---|---|---|---|---|---|
| dNK1 | CD39⁺ → CD39⁻ | **XCL1** | 3 | **+6.7** | 3/3 |
| dNK1 | CD39⁺ → CD39⁻ | **XCL2** | 3 | **+5.7** | 3/3 |
| dNK1 | CD39⁺ → CD39⁻ | CCL5 | 3 | +2.5 | 3/3 |
| dNK1 | CD39⁺ → CD39⁻ | **CXCR4** | 3 | **−3.5** | 3/3 |
| dNK2 | CD103⁻ → CD103⁺ | **CCL5** | 3 | **+6.0** | 3/3 |
| dNK2 | CD103⁻ → CD103⁺ | XCL1 | 3 | −5.2 | 3/3 |
| dNK2 | CD103⁻ → CD103⁺ | CXCR4 | 3 | +3.1 | 0/3 |

n = 3 throughout, so `min_achievable_p` = 0.25 — **no p value is meaningful
here either**. Sign concordance and magnitude are all that is claimed.

**Conclusions.**

1. **XCL1/XCL2 survive inside dNK1** (+6.7/+5.7 pp, all 3 donors), so the
   CD39⁻-over-CD39⁺ direction is real structure, not a restatement of the
   clustering.
2. **CCL5 survives inside dNK2** (+6.0 pp, all 3 donors) along the CD103 axis.
3. **CXCR4 REVERSES inside dNK1** (−3.5 pp, all 3 donors concordant) while the
   between-cluster effect is strongly positive. Its between-cluster gradient
   does not reproduce when the cluster is held fixed. **CXCR4 is the most
   circular of the four**, which is exactly the row that carried the study's
   largest z (+14.56).
4. Within-label effects are ~17% of between-cluster ones. **The direction is
   established as non-circular; the magnitude is not.**

**Consequence for RESULTS.** The eight positives are re-ranked. The strongest
claim is **XCL1/XCL2, dNK2 > dNK1** — it survives both the marker regrouping
and the within-label test, and it is a cross-modality reproduction of Huhn
2020's CyTOF protein result. The weakest is **CXCR4 in the dNK3 contrasts**,
despite its z. The z ordering and the credibility ordering are inverted.

---

## D19 — The dNK1 residual is not explained by library T-cell content.

D16 left three candidates for the 6.6% depth-matched dNK1 triple-positive
residual. Only one is decidable from this data: that dNK1 is enriched in
T-rich libraries.

Per decidual library with ≥ 30 dNK1 cells (n = 7), triple-positive rate at
matched depth against the library's T-cell fraction:

| library | T fraction | dNK1 triple-pos |
|---|---|---|
| FCA7474062 | 4.7% | 4.2% |
| FCA7196218 | 4.9% | 10.0% |
| FCA7167219 | 9.6% | 4.7% |
| FCA7196224 | 10.6% | 4.5% |
| FCA7167223 | 12.0% | 3.7% |
| FCA7167221 | 12.8% | 2.1% |
| FCA7511881 | 14.5% | 13.2% |

**Pearson r = +0.061 (p = 0.90); Spearman r = −0.036 (p = 0.94).** No
relationship. Figure: `out/VT2018/residual_vs_library_Tfrac.png`.

**Decision.** The library-composition candidate is **excluded**. At n = 7
libraries this excludes a strong relationship, not a weak one. The remaining
two candidates — annotation-boundary cells, and genuine low-level CD3D/CD3E
transcription in NK — are not decidable from this dataset, and the residual
stays flagged as unexplained. It does not block anything: the purged arm gives
the same eight positives.

---

## D20 — Reporting corrections carried into v1.2.

Three items raised against the v1.1 report, all fixed in code rather than
in prose:

- **BH family (v1.2 A5).** v1.0's "correct across the panel" did not say
  whether three pairwise contrasts are one family or three. Now one family of
  81. Effect: all 8 rows still clear, but two land at **q = 0.0499**, inside
  by 0.0001 — reported as the knife-edge it is.
- **Empirical-p resolution (v1.2 A6).** `q < 0.0001` was unsupportable: a
  rank-based p over 405 nulls cannot resolve below 1/406 = 0.0025. Now
  computed as `(n_ge+1)/(N+1)`, floored, with the floor emitted per row. Six
  rows are reported as "≤ 0.0025, at the resolution limit". z is kept as a
  standardised effect size and **not** converted to a p.
- **Ruler B denominator (v1.2 A7).** Ranking sources by per-cell CPM let a
  184-cell ILC3 population own XCL1's denominator. Now ranked by **total
  counts contributed**, and NK being the largest contributor forces
  `unmeasurable`. Effect: XCL1, XCL2, CCL5, CXCR4 and CCL4 all move from
  ruler-B positive to `unmeasurable`, so **no headline gene meets the
  "both rulers positive" bar any more**. Their ambient evidence rests on
  ruler A alone plus the `no_competing_source_NK_is_largest_contributor`
  reason code.

PAEP's absence (D17) is downgraded to a limitation: ruler A is already
saturated at 1.000 on five tissue controls, so a sixth cannot raise the
ceiling. The residual risk is confined to ruler B's band upper edge, now set
by LYZ/Myeloid at 0.0538, which PAEP could in principle raise.

---

## D21 — Ruler A did not drift. Three different quantities were reported as one.

**Challenge.** The v1.1 report said the tissue ambient controls "all return
soup fraction 1.000"; the v1.2 report said ruler A is "0.021–0.048 against a
0.092 ceiling". Same estimator or a silent definition change?

**Same estimator, verified by history, not by memory.** The ruler-A block and
the `soup_fraction` formula appear in exactly one commit — `bcf1947`, the
original pipeline commit — and have never been edited. v1.1 and v1.2 changed
ruler B only.

**Three quantities, not one:**

| quantity | value | what it is |
|---|---|---|
| ambient controls' soup fraction | DCN/COL1A1/IGFBP1/HLA-G/CSH1 = **1.000**; IGKC 0.092, LYZ 0.164, C1QA 0.176, HBB 0.408 | pure-ambient genes, true value 1.0 |
| **conservative ceiling** | **0.092** | the **minimum** of the row above |
| true-NK floor | 0.021 | KLRD1/NKG7/GNLY/PRF1/GZMB/NCR1 |
| target genes | 0.021–0.048 | XCL1 0.021, CCL5 0.025, CXCR4 0.048 |

**But checking it exposed an error in the v1.2 CHANGELOG, retracted here.**
That entry argued PAEP's absence was harmless because "ruler A is already
saturated at 1.000 on five tissue controls, so a sixth cannot raise the
ceiling". **The ceiling is the MINIMUM, not the maximum.** Those five 1.000s
do not set it; IGKC's 0.092 does. The argument was written into the record
without checking the one line of code that defines the ceiling.

The correct risk statement: PAEP matters only if it would return a value
**below** 0.092, and then it lowers the ceiling. **If PAEP would return
anything between 0.048 and 0.092, the ceiling drops below CXCR4's 0.0484 and
CXCR4's ruler-A verdict flips to `indistinguishable_from_ambient`.** XCL1
(0.021) and CCL5 (0.025) have far more margin. This is a third independent
line pointing at CXCR4 as the weakest of the headline genes, and it was
concealed by the earlier wrong argument.

---

## D22 — The empirical null is not robust to how it is constructed. This is the round's most consequential negative.

Raising the null from 405 to ~10,000 genes was meant to be a precision fix.
Implementing it exposed two defects in the v1.1 null and then a limit that no
implementation removes.

**Defect 1 — one pooled null for all targets.** v1.1 pooled every null gene
into a single distribution and judged all 27 targets against the same mu/sd.
That is not expression matching: a target at 2,315 CPM was being judged
against a pool with median 5.5 CPM, whose variance is far smaller, inflating
every z. Fixed: each target now gets its own matched null set
(`null_set_is_expression_matched`).

**Defect 2 — null genes were consumed, not shared.** Each target banned the
genes it took, starving whichever targets came later — worst exactly at the
extremes of expression, where the headline genes live. Measured before the
fix: **XCL1 received only 98 null genes** (so its empirical p of 0.0101 *was*
its resolution floor, 1/99), and **XCL2 and CXCR4 received none at all** and
silently fell back to the pooled distribution. Fixed: null sets are shared
between targets, since each is evaluated independently.

**The limit that remains.** With per-target matching, the null size per target
is bounded by how many genes of comparable expression exist. Requesting 10,000
yields 8,537 unique genes but only ~317 per target — **resolution 1/318 =
0.0031, not 1/10001**. Widening the window to get more genes degrades the
matching it exists to provide. Precision and matching trade off, and this
dataset cannot have both.

**Consequence — the pass list moves:**

| version | matched? | resolution | rows at q ≤ 0.05 |
|---|---|---|---|
| pooled, 399 nulls | no | 0.0025 | 6 |
| per-target, 112 nulls | yes | 0.0088 | **0** |
| per-target, 317 nulls | yes | 0.0031 | 8 |

The middle row has zero passes purely because 112 nulls cannot resolve a p
small enough to survive BH over 81 tests — not because the effects weakened.

**And a new failure mode appeared.** Under per-target matching, `S1PR5` passes
at q = 0.032 on an effect of **+1.5 pp**, because its own null set consists of
low-expression genes whose detection rates are compressed against the zero
floor, giving a null SD of 0.30. **Statistically significant, biologically
meaningless.** No effect-size threshold was added after the fact to remove it;
it is reported as it stands, as a warning about the method.

**Decision.** Report the per-target, 317-null version as the methodologically
correct one, publish the sensitivity table beside it
(`out/VT2018/null_sensitivity/`), and state plainly: **no row passes under
every construction tried.** The empirical null is the best instrument
available at n = 5–6, and it is not a stable one. Weight the
construction-independent evidence — the within-label circularity test —
accordingly.

Under the correct version, XCL1 dNK1→dNK2 sits at **q = 0.0463**, and CCL5
dNK1→dNK2 at 0.0318.

---

## D23 — Ruler A's soup fraction splits by lineage CLASS. PAEP's absence is low-risk after all.

The nine panel ambient controls split perfectly: the four that do **not**
saturate (IGKC 0.092, LYZ 0.164, C1QA 0.176, HBB 0.408) are all
haematopoietic; the five that saturate at 1.000 (DCN, COL1A1, IGFBP1, HLA-G,
CSH1) are all non-haematopoietic. Tested by pulling marker genes for every
mapped lineage **straight from the data** and running them through the same
leave-one-out estimator (`src/ceiling_class_structure.py`).

| class | n | median soup fraction | saturated at 1.000 |
|---|---|---|---|
| haematopoietic | 56 | 0.250 | **1.8%** |
| non-haematopoietic | 60 | **1.000** | **90%** |

Mann–Whitney p = 3.1e-19. Gene-level, not donor-level: it characterises the
estimator, not a biological claim, so R1 does not apply and nothing from it
enters T2.

**A first pass got this wrong and the fix matters.** Ranking candidate markers
by lineage CPM alone made Stromal look like an exception (0% saturated). The
cause was selection, not biology: the lineage with the most cells (Stromal,
12,583) owns the total counts of every housekeeping gene too, so the top-CPM
list filled with ribosomal genes that NK also expresses, dragging their soup
fraction down. Requiring 5× specificity over the next-highest lineage put
Stromal at 100% saturated, in line with the panel's own DCN/COL1A1/IGFBP1.

**Consequence for PAEP.** Glandular-epithelial/stromal markers are 92%
saturated with median 1.000. PAEP is glandular epithelial / decidualised
stroma, so it would be expected to return 1.000 and **leave the ceiling
untouched** — the ceiling being a *minimum*. **The CXCR4 risk from PAEP's
absence, raised in D21, is therefore small.** D21's risk statement is
downgraded accordingly; its correction of the *mechanism* (ceiling = minimum,
not maximum) still stands.

**And it exposes a deeper defect.** The ceiling of 0.092 is set entirely by
the haematopoietic class — which is the class the estimator handles worst,
because NK cells are themselves haematopoietic and rho's single-load
assumption fails for same-class genes. A minimum over controls is therefore
set by the estimator's worst case and falls monotonically as controls are
added. Recorded in `PREREGISTRATION_v1.3.md` as known and deferred.

---

## D24 — The FDR criterion is withdrawn, not adjusted.

Walk the arithmetic: with N nulls per target the rank floor is 1/(N+1), and BH
at m = 81 needs `0.05 × k / 81` for k rows tied at that floor — 15 rows at
N = 112, 6 at N = 317, 5 at N = 371. **Only 8 rows are candidates in the whole
study.** The same data gave 6 passes, then 0, then 8, as N alone changed.

The pass count is a function of the null-set size, not of the effects. That is
not a rule needing a tuned parameter; it is a rule that does not hold at this
n — the same failure as ruler B's fixed denominator (v1.2 A7) and BH over the
signed-rank p (v1.1 A2), both of which were also withdrawn rather than tuned.

**Decision.** No q value, no pass/fail list, no replacement rule. Rows carry
`fdr_estimable = FALSE` with the reason. Descriptive quantities are reported:
effect in points, z as a standardised effect size never converted to a p,
donor sign concordance, ruler-A margin, and the observed effect's rank within
its own null set.

Retaining BH as a "reported but non-adjudicating" column, as v1.1 A2 did, was
a half-measure: a number printed beside a result is read as a verdict whatever
the caption says.

---

## D25 — Matching nulls on detection rate was the right diagnosis and did not fix it.

The readout is a detection rate, bounded in [0,1], with sampling variance
p(1−p)/n set by the rate rather than by CPM. Matching nulls on CPM therefore
gave low-detection targets a null set compressed against the zero floor for
unrelated reasons — which is how a **+1.5 pp** S1PR5 effect reached the top of
the v1.2 list. Matching was moved to the reference arm's **baseline detection
rate**, computed over the whole transcriptome at matched depth, per contrast.

**It made the problem sharper, not smaller.**

| baseline detection | n | median null SD (pp) | median abs effect (pp) |
|---|---|---|---|
| < 0.5% | 29 | 0.33 | 0.17 |
| 2–10% | 8 | 2.06 | 2.07 |
| > 30% | 17 | 6.41 | 10.15 |

**A z of 5 buys a 32.1-point effect at >30% baseline detection and a
1.63-point effect at <0.5% — 20-fold.** S1PR5 went from z = +5.01 to
**+10.34** on the same +1.5 pp; SELL entered at z = +10.55 on +4.8 pp.

**Conclusion, recorded as a limit rather than patched.** The defect is not the
matching axis. **z is not comparable across genes on a bounded scale near its
boundary, however the nulls are chosen.** So z may never be read without the
percentage-point effect beside it; both are emitted in `null_descriptive.tsv`,
and no z-alone ranking appears in `RESULTS.md`. No effect-size threshold was
added after the fact to suppress S1PR5 or SELL — they are left in as evidence
about the instrument.

**The pattern across v1.1–v1.3.** Four instruments have now been withdrawn
rather than tuned: BH over the signed-rank p, ruler B's fixed denominator,
the FDR criterion, and z as a standalone ranking. In each case the failure was
that the output was determined by a nuisance parameter — n, the denominator
lineage, the null-set size, the baseline rate — instead of by the effect.
**What survived all four withdrawals is the evidence that never used any of
them:** donor-unanimous direction, effect size in points, and reproduction
within a fixed cluster.

---

## D26 — CXCR4's ceiling margin is now a probability, and it is negative on every axis tested.

> **SUPERSEDED by [D29](#d29) and [D31](#d31).** The probability table below is
> withdrawn (controls are specific genes, not exchangeable draws) and the
> "every axis" claim is retracted: CXCR4 is negative on **two** axes, not four.
> Kept unedited as the record of what was concluded and when.

The ruler-A ceiling is a **minimum** over ambient controls, and D23 showed it
is set entirely by the haematopoietic class, whose median soup fraction is
0.250. So 0.092 is an extreme lower tail of 56 draws, and adding controls can
only lower it.

Empirical distribution of the 56 data-derived haematopoietic controls:

| threshold | P(a single haematopoietic control falls below) |
|---|---|
| 0.092 — the current ceiling (IGKC) | 0.107 (6/56) |
| **0.048 — CXCR4** | **0.054 (3/56)** |
| **0.021 — XCL1** | 0/56 — see the correction below |

Probability the ceiling drops below a gene after adding k more haematopoietic
controls:

| k | below CXCR4 (0.048) | below XCL1 (0.021) |
|---|---|---|
| 1 | 5.4% | 0.0% |
| 3 | 15.2% | 0.0% |
| 5 | **24.1%** | 0.0% |
| 10 | **42.3%** | 0.0% |

Only 4 of the panel's 9 ambient controls are haematopoietic, so "a few more"
is an ordinary scenario, not a contrived one.

**Correction to the 0.0% figure.** 0/56 is a point estimate, not zero. By the
rule of three the 95% upper bound on a single control falling below 0.021 is
about 3/56 = 5.4%, so the worst case over five added controls is ~24% — the
same as CXCR4's. **XCL1's real defence is not that probability; it is the
4.4-fold gap between 0.021 and the observed minimum of 0.043**, against
CXCR4's 1.1-fold gap. Stating it as "0.0%" overstated the asymmetry.

**And treating "add a control" as random sampling is itself a heuristic.**
Controls are specific genes, not exchangeable draws. The probability framing
communicates the instability of a minimum-over-controls ceiling; it is not a
fix. The lower-prediction-bound construction (v1.3, deferred) remains the
actual repair.

**CXCR4 is now negative on every axis this project has measured:**

| axis | result |
|---|---|
| circularity — marker regrouping | retains **0.30×**, the worst of the four |
| circularity — within-label | **reverses**, −3.5 pp, 3/3 donors |
| ruler-A ceiling margin | ~24% chance of flipping after 5 more controls |

Three independent axes agreeing in direction is also a positive result *about
the pipeline*: the axes were built to answer different questions and were not
tuned to agree.

**XCL1 passes all three**: survives regrouping at 0.62×, survives within dNK1
at +6.7 pp (3/3), and no haematopoietic control in 56 falls below its 0.021.

---

## D27 — More than half the primary panel has no power on the primary scale.

Measured, pooled over decidual dNK1–3 at the primary depth floor: **14 of 27**
primary targets sit below 5% detection, 12 in the measurable 5–85% band, and
1 above 85%. Full lists in `PREREGISTRATION_v1.4.md` A11.

**The two failure modes close on the same genes.** Four of the 14 floor-band
genes — S1PR5, SELL, CCR9, CXCR6 — are precisely those that produced the
largest z values in v1.3 (S1PR5 z = +10.34 on +1.5 pp). The floor band has
**no power in percentage points and false power in z, simultaneously**. That
is why the remedy is a declared scale rule (v1.4 A10) rather than an
after-the-fact effect-size filter.

**A correction.** The suggested reading that CCL3's flat result "may be no
power rather than no difference" is **not supported**: CCL3 is at 0.817,
inside the measurable band, so its flat and sign-inconsistent result carries
power behind it. CCL4 (0.944) is ceiling-limited, which is a different problem
with the opposite remedy — more depth makes saturation worse, not better. The
power-limited genes are the 14 listed, and neither CCL3 nor CCL4 is among
them.

---

## D28 — CCR5 is in the floor band. The check transfers; the number does not.

`CCR5` sits at **68 detected of 10,130** depth-matched decidual dNK1–3 cells =
**0.67%**, deep inside the floor band declared in v1.4 A11 — the band shown to
have *no power in percentage points and false power in z simultaneously*.

Per-subset, at the 2,137-UMI floor: dNK1 0.15%, dNK2 0.91%, dNK3 1.43%,
Myeloid 2.51%, T 4.00%. At p ≈ 0.007 and n ≈ 2,000 cells the binomial
resolution is about 0.4 percentage points, so a detection-rate contrast on
CCR5 in this dataset can only see effects far larger than any plausible
biological difference at that abundance.

**This matters beyond this project** because any archived analysis whose
conclusion rests on a CCR5/Ccr5 detection rate from a droplet platform was
operating in that band, and would carry both failure modes at once.

**But the number does not transfer, and saying otherwise would repeat the
exact error this project has now made twice.** 0.67% is a fact about human
first-trimester decidual NK on 10x v2 at a 2,137-UMI floor. Chemokine-receptor
expression is not conserved across species, tissue or activation state, so it
implies nothing about *Ccr5* in mouse tumour-infiltrating NK. Both prior
errors here — the PAEP ceiling argument and the CCL3 band assignment — had the
same shape: **a correct general principle applied straight to a specific
number without checking that number.**

**What transfers is the check**, and it is one number:

1. compute the gene's mean detection rate in the **reference** condition at
   matched depth;
2. **< 5%** → floor band: no power in points, inflated standardised
   statistics; any published effect from that band needs re-reading;
3. **> 85%** → ceiling: compressed, and a null there is uninterpretable;
4. in between → usable; report in points.

Written up for reuse in [`TRANSFERABLE.md`](TRANSFERABLE.md) §4. **The data
needed to run it on the MC38/Kaede cohort is not in this workspace** — the
same workspace that holds the donor column blocking D15. Both are one lookup
each and can be done in the same sitting.

---

## D29 — The ceiling definition set is declared. Two margins had been computed against two different sets.

**The error.** D26 stated "XCL1's real defence is the 4.4-fold gap to the
observed minimum of 0.043, against CXCR4's 1.1-fold." Those two multiples come
from **different minima**: 0.021 × 4.4 = 0.092 (the panel controls' minimum),
while 0.048 ÷ 0.043 = 1.1 (the class-structure markers' minimum). One sentence,
two definition sets. Same failure class as the two errors before it.

**And the consequence was larger than a wrong multiple.** If 0.043 were the
operative ceiling, CXCR4 (0.0484) would already sit **above** it — its verdict
would be `indistinguishable_from_ambient` *today*, not at some future 24%
probability. The probability framing in D26 was also invalid on its own terms:
those 56 controls already exist, so "adding five more" describes nothing.

**The decision, which had never been made explicitly.**

| definition set | ceiling | XCL1 (0.021) | CXCR4 (0.048) | CXCR4 verdict |
|---|---|---|---|---|
| **A — the 9 pre-registered panel ambient controls** | **0.0920** (IGKC) | **4.4×** | **1.9×** | `above_ambient` |
| B — the 56 data-derived haematopoietic markers | 0.0430 (**CCL3**) | 2.1× | 0.89× | `indistinguishable_from_ambient` |

**Set A is operative, and the reason is structural, not convenient.** The 56
markers of set B were auto-selected *after the fact* to characterise the
estimator's class behaviour (D23), and the three that fall below CXCR4 are
`CCL3` (0.0430), `LINC00861` (0.0438) and `IFNG` (0.0452). **`CCL3` is a
primary target of this study.** Letting set B define the ceiling would let a
gene under test set the threshold that judges genes under test — circular in
exactly the way the whole project has been guarding against. Set B is a
diagnostic set; it was never an ambient-control set and must not become one.

Note the direction: choosing A is the choice that lets CXCR4 *pass* ruler A.
It is not being chosen to strengthen the case against CXCR4.

**What this changes in the record.**

- Both margins are now quoted against set A: **XCL1 4.4×, CXCR4 1.9×.**
- **D26's probability table is withdrawn.** "P(the ceiling drops below gene X
  after adding k controls)" was never a well-posed quantity: controls are
  specific genes, not exchangeable draws, and the draws in question already
  existed.
- Replacing it is a sharper and simpler statement: **CXCR4's ruler-A verdict
  depends on which set defines the ceiling, and XCL1's does not.** XCL1 clears
  both (4.4× and 2.1×); CXCR4 clears one and fails the other (1.9× and 0.89×).

**So CXCR4 is negative on a fourth axis** — not "at risk in future", but
*conditional on a definition that was only settled today*:

> **RETRACTED by [D31](#d31).** The three set-B markers below CXCR4 are all
> genes NK expresses itself, so set B never carried a verdict here. The table
> below keeps only its first two rows; the ruler-A rows are withdrawn.

| axis | CXCR4 | XCL1 |
|---|---|---|
| circularity, marker regrouping | 0.30× retained, worst of four | 0.62×, best |
| circularity, within-label | **reverses**, −3.5 pp, 3/3 donors | survives, +6.7 pp, 3/3 |
| ruler-A margin, set A | 1.9× | 4.4× |
| ruler-A verdict under set B | **flips to indistinguishable** | holds |

The deferred repair is unchanged and now better motivated: a minimum over
controls is set by whichever control the estimator handles worst, and the
lower-prediction-bound construction within lineage class replaces it.

---

## D30 — On the decidual side, CCR5 is not detectably expressed at all.

The per-subset breakdown says more than "CCR5 is in the floor band". At the
2,137-UMI floor: dNK1 **0.15%**, dNK2 0.91%, dNK3 1.43%, dNKp 0.64% — against
**T cells 4.00%** and **Myeloid 2.51%**. Every dNK subset is below both
lineages that actually express it.

Ruler B agrees on the source: CCR5's dominant pool contributor is **Myeloid**
(57.5% of the non-NK pool), NK is *not* the largest contributor, and the
NK/Myeloid ratio is 0.388 — but the verdict is `unmeasurable`, because 0.388
is above the 0.0538 pickup band while ruler A returns 0.0575, just below the
0.092 ceiling. The two rulers disagree, which at this abundance is what
"no signal to adjudicate" looks like.

**The stronger and better-supported statement for the decidual side is
therefore: in this dataset, CCR5 expression is not detectable in dNK at
transcript level** — the signal present is at or below what the two
higher-expressing lineages would deposit as pickup.

**This does not transfer to mouse tumour NK** (D28), and the point of stating
it is narrower: if a CCR5-based line were ever to be extended to the decidual
side, this dataset offers no starting point for it.

---

## D31 — Set B never carried a verdict on CXCR4. Retracting the fourth axis.

D29 excluded set B (the 56 data-derived haematopoietic markers) because it
contained `CCL3`, a primary target. **That reason was too narrow and the
conclusion drawn from it was wrong in CXCR4's disfavour.**

**All three markers below CXCR4 are genes NK expresses itself:**

| gene | NK CPM | NK detection | NK / source lineage | |
|---|---|---|---|---|
| CCL3 | 2223.6 | **90.9%** | 0.38 | NK expresses it heavily |
| IFNG | 27.3 | 6.5% | 0.11 | an NK effector gene |
| LINC00861 | 7.0 | 3.0% | 0.11 | a T/NK-lineage lncRNA |
| — the panel's real controls — | | | | |
| IGKC | 19.9 | 5.2% | **0.002** | plasma-cell Ig |
| LYZ | 91.0 | 19.4% | **0.05** | myeloid |
| C1QA | 174.1 | 42.1% | **0.04** | macrophage |

**Their low soup fractions are correct readings, not estimator failures.** NK
really does express them, so a low "ambient fraction" is the right answer —
and it carries **no information about an ambient ceiling.** Set B never
delivered an adverse verdict on CXCR4 that could be discounted; it delivered
none at all.

**Confirmed by recomputation.** Imposing NK/source ≤ 0.05, the criterion the
three real controls meet, leaves 49 of 56 markers, and the set's minimum rises
from 0.043 to **0.0821 — above CXCR4's 0.0484, with zero markers below it.**
Both readings agree: the original set is inadmissible, and the filtered set is
not adverse.

**Why this is not repairable by a better filter.** That NK/source screen uses
`nk_cpm`, and `soup_fraction = rho x amb_cpm / nk_cpm` has the same
denominator, so screening on low NK expression necessarily selects high soup
fractions. **In the haematopoietic class there is no data-driven way to
separate "NK expresses it" from "NK picks it up" — that separation is
precisely what ruler A exists to measure.** Selecting ambient controls from
data works on the non-haematopoietic side only, where prior lineage knowledge
guarantees NK cannot transcribe stromal or trophoblast genes.

**Corrections that follow.**

1. **The fourth axis against CXCR4 is withdrawn.** The ruler-A statement is
   set A alone: ceiling 0.0920, **CXCR4 margin 1.9×, XCL1 4.4×**.
2. **CXCR4 is negative on TWO axes, not four:** circularity retention 0.30×
   (worst of four genes), and within-label sign reversal (−3.5 pp, 3/3
   donors). The conclusion is unchanged — the within-label reversal alone is
   decisive — but "negative on every axis measured" overstated it and must not
   be written that way.
3. **D23's haematopoietic quantification is qualified.** Median 0.250, 1.8%
   saturated and p = 3.1e-19 were computed over a set containing NK-expressed
   genes, so they are not a clean estimate of the estimator's behaviour on
   true haematopoietic ambient genes. Under the (circular) NK/source ≤ 0.05
   screen the median moves only 0.250 → 0.270 and saturation 1.8% → 2%, so the
   quantities are insensitive to it — but insensitivity is not cleanliness.
   **The class effect itself stands**, resting on the panel's three
   prior-selected controls (IGKC 0.092, LYZ 0.164, C1QA 0.176, none saturated)
   against five non-haematopoietic controls (all 1.000).
4. **HBB should not have been counted as a fourth haematopoietic control.**
   Its nominal source, Myeloid, has a CCR5-like problem: myeloid CPM for HBB is
   **1.9**, and decidua contains no erythroid population, so HBB has no source
   lineage here and its 0.408 is uninterpretable.
5. **The PAEP downgrade (D23) stands**, since it rests entirely on the
   non-haematopoietic side, which this contamination does not touch.

**This correction moves in CXCR4's favour, and is made for that reason:** the
whole review had been accumulating evidence against CXCR4, and this piece of it
does not survive checking.

---

## D32 — Two corrections to the CCR5 statements.

**The "rulers disagree" claim was wrong.** Ruler A returns 0.0575, below the
0.092 ceiling → `above_ambient`. Ruler B returns 0.388, above the 0.0538 band.
**Both point the same way**; there is no contradiction. The `unmeasurable`
verdict comes from the power gate, not from a conflict: CCR5's dominant pool
contributor is Myeloid, whose own CCR5 detection is **4.0%**, below the 5%
gate, and whose CPM is 9.37, below the 10 CPM gate. The verdict means *the
source lineage does not express it enough for the ratio to have power* — which
is the gate working, not the rulers fighting.

**And "CCR5 expression is not detectable in dNK at transcript level" is
withdrawn as overstated.** The floor band means the **readout has no
comparative power**, not that the gene is absent; 0.67% is entirely compatible
with genuine expression in a rare subpopulation.

The defensible version: **CCR5 detection in dNK is 0.67%, below myeloid (2.51%)
and T (4.00%) in the same tissue, and inside the band where this study has no
power. This dataset establishes neither expression nor its absence.** The
lineage ordering is suggestive of pickup, and that is all it is.

---

## D33 — Scope correction to D31, and two numeric-provenance stamps.

**D31 over-corrected.** It said "in the haematopoietic class there is no
data-driven way to separate NK-expresses from NK-picks-up", and then let that
read as though ruler A were unusable on the haematopoietic side. **The second
part is wrong.**

What fails is **auto-expanding the control set with a statistic**. What does
*not* fail is a control chosen from prior exclusivity knowledge, and the
panel's three are exactly that: "NK does not transcribe immunoglobulin"
(IGKC), "…lysozyme" (LYZ), "…complement C1q" (C1QA) are as hard a prior as
"…collagen". **Ruler A remains usable for haematopoietic target genes**, on
those three controls, ceiling 0.0920.

The practical consequence: XCL1's ruler-A reading (0.021 against 0.092) is not
withdrawn. XCL1's conclusion never depended on it — it rests on donor
unanimity, marker regrouping and the within-label test — but there is no
reason to discard a layer of evidence that stands.

**HBB, secondary correction.** D31 said HBB should not have counted as a
control because decidua contains no erythroid population. That is right for
**ruler B**, which needs a source lineage — and the current calibration
confirms it, assigning HBB's dominant source to Stromal at 8.3 CPM and
excluding it below the 10 CPM power gate. But it does **not** automatically
void its **ruler A** reading: in a well-perfused tissue, free haemoglobin
transcript comes from red cells never captured as cells at all, which is a
different contamination route from RNA shed by captured cells. HBB's
non-saturated 0.408 may be reporting that route rather than failing. Low
priority; recorded so it is not silently treated as settled.

**Two numeric-provenance stamps.** Both quantities below appeared at two
values across these documents. Both readings are correct at their own stage;
the drift was in not stamping them.

| quantity | value | reading |
|---|---|---|
| CCL3 detection in dNK | **0.817** | dNK1–3, depth-matched to 2137 UMI, per-subset mean — **the v1.4 A11 band assignment uses this** |
| | 0.909 | all dNK incl. dNKp, raw counts, no depth matching — the figure quoted in D31 |
| ruler B pickup band | **0.0538** | v1.2 onward, source ranked by **total counts** (LYZ→Myeloid) — **current** |
| | 0.0441 | v1.1, source ranked by **per-cell CPM** (LYZ→cDC1) — superseded; `PREREGISTRATION_v1.1.md` keeps it as the record and should not be edited |

Neither affects a conclusion: CCL3 is inside the measurable band on both
readings, and every ruler-B verdict in the current output was computed against
0.0538.

---

## D34 — "CCL3 does not reproduce" was too broad. The pipeline had it right; the prose did not.

**Trigger.** v1.4 A11 assigned CCL3 to the measurable band on a **cross-subset
mean of 0.817**, only 3.3 points below the 0.85 ceiling threshold. A mean is
compatible with one arm sitting above it, and R5 (as amended by v1.1 A1) is a
**per-arm** rule. Checked:

| subset | per-donor mean | cell-weighted |
|---|---|---|
| dNK1 | 0.8003 | 0.8161 |
| dNK2 | 0.7772 | 0.7536 |
| **dNK3** | **0.8837** | 0.8646 |

**dNK3 is above the ceiling.** The band assignment used a statistic the
decision rule does not use.

**The pipeline was already correct.** T2 flags exactly what v1.1 A1 requires:

| contrast | effect | ceiling arm | flag |
|---|---|---|---|
| dNK1 vs dNK2 | −2.31 pp | none | clean |
| dNK1 vs dNK3 | **+8.79 pp** | dNK3 | `magnitude_is_lower_bound = TRUE` |
| dNK2 vs dNK3 | **+10.15 pp** | dNK3 | `magnitude_is_lower_bound = TRUE` |

**So the claim must be split, and half of it withdrawn.**

- **dNK1 vs dNK2: a real null.** −2.3 pp with both arms inside the measurable
  band. This is a powered true negative, and D27's wording holds here.
- **Both dNK3 contrasts: not a null at all.** +8.8 and +10.1 pp, both pointing
  the same way (dNK3 higher), against a saturated dNK3 arm that compresses the
  difference toward zero. **A one-sided-ceiling null is not evidence of
  absence** — direction is readable, magnitude is a lower bound, and "no
  difference" cannot be adjudicated.

**"CCL3 does not reproduce a subset difference" is therefore withdrawn as
stated.** The supported version: *CCL3 shows no dNK1-vs-dNK2 difference
(−2.3 pp, both arms measurable), and its dNK3 contrasts are ceiling-limited
lower bounds of +8.8 and +10.1 pp that this dataset cannot adjudicate.*

**Not a revival of the earlier claim.** The suggestion rejected in D27 — that
CCL3's flat result might be a power problem — was argued from the **floor**
band and was wrong; 0.817 refuted it. This is the opposite end: **the ceiling**,
and it was triggered by this project's own R5 line rather than by judgement.

**Process fix.** Band assignment must be reported **per arm**, matching the
rule that consumes it. A cross-subset mean can place a gene in the measurable
band while one arm is saturated, which is exactly what happened. The
per-contrast flags in T2 remain authoritative.

---

## D35 — The same power gate caught two unrelated failures.

Worth one line as evidence of internal consistency. The `>= 10 CPM` gate on a
source lineage caught, independently:

- **CCR5** (D32): dominant source Myeloid at **9.37 CPM** → ruler B
  `unmeasurable`, correctly, since myeloid barely expresses it either;
- **HBB** (D33): dominant source assigned to Stromal at **8.3 CPM** →
  excluded from ruler B calibration, correctly, since decidua has no erythroid
  population for it to come from.

Two different scenarios — a genuinely low-abundance receptor, and a gene whose
source lineage is absent from the tissue — caught by one rule that was not
tuned for either.

---

## D36 — CCL3's dNK3 effect does not survive the circularity test either, and circularity cannot explain it.

The follow-up was already computed: `CCL3` has been in `circularity_test.py`'s
target list from the start, and only the printout filtered it out. Read from
the existing outputs, no rerun.

**Prior expectation, and why it was reasonable.** `CCL3` is **not** among
Vento-Tormo's dNK3 characterisation genes (CCL5, CXCR4, XCL1, CD160, KLRB1,
GZMK, IFNG, TNFSF14, TIGIT), so its dNK3 effect sat in the *low*-circularity
tier alongside XCL1, not with CCL5/CXCR4. The test was expected to support it.

**It does not.**

| gene | retained under marker regrouping | within-label |
|---|---|---|
| XCL1 | **0.62×**, 4/4 donors concordant | dNK1: +6.7 pp, 3/3 |
| XCL2 | 0.53×, 4/4 | dNK1: +5.7 pp, 3/3 |
| CCL5 | 0.46×, 4/4 | dNK2: +6.0 pp, 3/3 |
| CXCR4 | 0.30× | dNK1: **−3.5 pp**, 3/3 (reversal) |
| **CCL3** | **0.28× — the lowest of the five** | dNK1: −2.3 pp; dNK2: +0.25 pp, **0/3 concordant** |

All three of CCL3's gate contrasts return `all_donors_same_sign = FALSE`,
which none of the other four do in the contrasts that carry them.

**Circularity cannot be the explanation**, since CCL3 is not a dNK3
characterisation gene. The available explanations are ordinary ones: gate
dilution bites hardest on a gene whose arms are *both* high (0.78–0.88), where
mixing compresses an already small difference; and the +8.8 pp label-based
effect was itself a ceiling-limited lower bound of unknown size.

**Status, which is a third category and not either of the earlier two.**
D34 corrected "CCL3 does not reproduce" to "the dNK3 contrasts are not nulls —
direction readable, magnitude a lower bound". The gate test now removes the
*direction* claim as well: **CCL3's dNK3 contrasts are unresolved.** The
ceiling prevents adjudicating a null, and the marker-gated test does not
support the effect. Neither established nor excluded.

**Consequence for a reading that is now withdrawn before it was made.** The
attractive inference — that CCL3 (dNK1 ≈ dNK2 < dNK3), CCL5 (monotonic) and
XCL1 (dNK1 < dNK2 ≈ dNK3) show three *different* patterns, implying subsets
carry distinct chemokine combinations rather than different amounts of one
programme — **requires CCL3's pattern to be established, and it is not.**
The pattern claim rests on CCL5 and XCL1 alone, which differ from each other,
and that much is unchanged.

This was the last open scientific follow-up in the run. It closes negative.

---

## D37 — Provenance stamp: 0.817 vs 0.8204.

Same quantity, two weightings, neither wrong:

| value | weighting |
|---|---|
| **0.8167** | every `donor × subset` row equally weighted (17 rows) — **what v1.4 A11 used** |
| 0.8204 | each subset averaged first, then the three subsets equally weighted |

The gap exists because **dNK3 clears the admission threshold in 5 donors while
dNK1 and dNK2 clear it in 6**, so row-weighting gives dNK3 — the highest and
the saturated arm — slightly less weight. Stamped rather than reconciled; the
per-arm figures (0.800 / 0.777 / 0.884) are what the ceiling rule consumes and
neither summary enters a decision.
