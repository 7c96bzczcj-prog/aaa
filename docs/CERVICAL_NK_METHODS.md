# Recounting NK cells in cervical tissue from public single-cell data

## Why recount at all

Published cervical single-cell studies almost uniformly report a merged
`T/NK`, `NK/T` or `T&NK` cluster. An NK-only percentage cannot be read off those
figures at any denominator, so the number has to be recomputed from the count
matrices. Everything below is computed from raw counts; no author cell-type
annotation is used at any point.

## Datasets

### Used

| Accession | Composition | Why included |
|---|---|---|
| **E-MTAB-12305** | 3 tumour-adjacent normal, 2 HSIL, 4 tumour, 1 metastatic LN (10 libraries) | The only source covering the full normal→HSIL→tumour→node series, and it is **paired** |
| **GSE208653** | 2 normal HPV-neg, 2 normal HPV-pos, 2 HSIL, 2 SCC, 1 ADC (9 libraries) | Independent normal→HSIL→cancer series |
| **GSE197461** | 3 SCC, 5 ADC (8 libraries) | Adds adenocarcinoma breadth |
| **GSE173231** | 10 ectocervix biopsies from 8 healthy donors | A *healthy-donor* normal baseline, not tumour-adjacent normal |

`E-MTAB-12305` is the accession given in the 2024 corrigendum
(Front Immunol 15:1386072); the accession printed in the original
Li & Hua 2022 paper (Front Immunol 13:897366) is wrong.

**Donor pairing in E-MTAB-12305** is not stated in the paper text but is
recoverable from the SDRF: `Characteristics[age]` ties N1/T1 (45 y), N2/T2 (50 y)
and N3/T3 (51 y) into three adjacent-normal/tumour pairs from the same donor, and
T4/L1 (48 y) into one tumour plus its own metastatic node. Those pairs carry
almost all of the interpretable signal, because they are the only comparisons
that hold dissociation protocol, batch and donor constant.

### Considered and rejected, with the reason

| Source | Status |
|---|---|
| eLife RP97335 (11 ADC + 4 SCC) | Deposited as **raw FASTQ only** (SRA `PRJNA1231266`); no count matrices. Also contains **no normal tissue**, so it cannot contribute a comparison. |
| Sci Adv `add8977` / PMC10277926 (5 normal, 4 HSIL, 5 MIC, 6 ICC) | Deposited in **GSA-Human** (`subHRA005765`), which is controlled-access and requires an approved data-access application. Not obtainable in this environment. This is the single most costly exclusion — it is the largest properly graded series. |
| Front Immunol 2025;16:1658705 / `GSE308792` (3 paired CSCC + adjacent) | The only supplementary file is an **integrated, normalized** expression CSV (values such as `0.399408…`, not counts). Ambient correction and doublet detection both require counts, so it cannot meet the analysis rules. Raw data is FASTQ-only. |

No accession here is guessed. Where data could not be obtained, that is stated
rather than substituted.

### Independence check

`E-MTAB-12305` (Li & Hua, Fudan), `GSE197461` (Qiu … Hua K, Fudan) and
`GSE208653` (same associated PMID 37794698) come from overlapping author groups,
so a re-deposit of the same library would fake replication. Two tests are run in
`scripts/cervical/check_overlap.py`: cell-barcode overlap against the ~737k 10x
whitelist expectation, and pseudobulk log-CPM correlation. Results in
`results/cervical/donor_overlap_*.csv`.

## Pipeline

Per library, in this order:

1. **Load** cellranger-filtered counts.
2. **QC**: ≥200 genes, ≥500 UMI, <20% mitochondrial reads.
3. **Doublet removal** — Scrublet (the scDblFinder/DoubletFinder equivalent
   available here). T–NK doublets are the largest single source of false NK
   calls, so a failure of this step is treated as fatal rather than skipped.
4. **Ambient RNA correction** — decontX (Yang et al., *Genome Biology* 2020).

   SoupX proper needs the raw empty-droplet matrix to build the soup profile.
   Every cervical dataset here is distributed as filtered matrices only, so the
   empty droplets do not exist to be read. decontX is the standard equivalent
   that estimates the contamination distribution from the cell matrix itself,
   using the other clusters as the ambient pool. It is implemented directly in
   `cerv_common.decontx` and validated before use (below).

Then per dataset (pooling only to define a shared clustering, never to report a
number):

5. Normalize, HVG, PCA, Harmony over libraries, Leiden → major lineages assigned
   by per-cluster argmax of lineage signature scores.
6. The **T/NK compartment is re-clustered from scratch** at higher resolution.
7. **NK gating** at cluster level, on positive *and* negative markers.
8. Per-donor tallies under three denominators.

### Ambient correction: acceptance test

decontX was not trusted on assertion. The test is that correction must strip
lineage-foreign transcripts from cells that cannot be making them, while leaving
those cells' own defining transcripts intact. Measured inside immune cells of
sample H1 (detection rate before → after):

```
KRT13   0.108 -> 0.000     foreign, should drop
KRT5    0.039 -> 0.006     foreign, should drop
EPCAM   0.009 -> 0.002     foreign, should drop
PTPRC   0.865 -> 0.863     own, should hold
CD3D    0.346 -> 0.333     own, should hold
GNLY    0.275 -> 0.186     own, but heavily soup-contaminated
```

The GNLY drop is the point of doing this at all: GNLY/NKG7/GZMB are among the
most prominent genes in the ambient pool, so without correction a substantial
fraction of non-NK cells carry NK-marker counts they did not transcribe.

### NK gating rule

A cluster in the re-clustered T/NK compartment is called NK when **all** hold,
on decontaminated counts:

- at least **3 of 4** T markers negative (detection rate < 0.25):
  `CD3D`, `CD3E`, `CD3G`, `TRBC2`
- `NKG7` detected in ≥60% and `GNLY` in ≥40% of cells
- at least one NK receptor above threshold: `KLRD1` ≥0.40, `KLRF1` ≥0.20,
  `NCAM1` ≥0.10

Requiring cytotoxic-granule genes (`GNLY`, `NKG7`, plus `PRF1` reported) as well
as a receptor keeps ILC1s, which lack them, from being counted as NK.

### TRBC1 was tested and dropped

The task asks whether `TRBC1` is reliable in human data before using it.
Measured in `CD3D+` cells of E-MTAB-12305 sample H1:

```
TRBC1  0.616      TRBC2  0.822      TRAC  0.778      CD3E  0.793
```

`TRBC1` and `TRBC2` are allelically excluded — an individual T cell uses one
constant region or the other — so a genuine `TRBC1` detection rate should be
roughly 0.35–0.40, not 0.62. A rate that high in `CD3D+` cells indicates read
cross-mapping between the two highly similar constant regions. `TRBC1` is
therefore **not** used as a negative marker. `TRBC2`, at a detection rate
comparable to `CD3E`, is retained.

### Denominators

Three, because the literature mixes them and they differ several-fold:

- **of all cells** — every cell surviving QC and doublet removal
- **of CD45+ immune cells** — T/NK, B, plasma, myeloid, mast, pDC lineages
- **of lymphocytes** — T/NK, B, plasma lineages

## Reported limitations

These are properties of the data, not of the pipeline, and they bound how far
the numbers can be pushed.

- **Cross-dataset absolute values are not comparable.** Single-cell composition
  is dominated by dissociation protocol and by whether a library was
  immune-enriched. Only within-dataset — ideally within-donor — comparisons mean
  anything.
- **Capture efficiency is confounded with biology.** NK cells in tumour tissue
  carry roughly 25% less endogenous RNA than in matched normal, so part of any
  cross-tissue proportion difference reflects differential capture and QC
  survival, not differential abundance.
- **Ambient contamination is highest for exactly the NK marker genes**, and it
  cuts both ways: uncorrected `GNLY`/`NKG7` soup inflates NK-like signal in
  non-NK cells, while ambient `CD3` in true NK cells pushes them past the
  negative-marker gate. Both are corrected here, but correction is a model, not
  a measurement.
- **Small denominators.** Any donor contributing few T/NK cells produces a
  percentage with wide binomial error; those rows are flagged in the table.
