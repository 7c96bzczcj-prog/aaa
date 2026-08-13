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
