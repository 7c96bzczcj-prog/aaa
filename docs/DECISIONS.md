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
