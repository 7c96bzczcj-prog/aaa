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
