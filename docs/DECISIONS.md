# DECISIONS — NK lineage-assertion evidence screen

Every judgement call: what was done, why, and what the alternative was.
Appended in order. Entries written before retrieval are marked *(pre-retrieval)*.

---

## DEC-01 — The spec says 24 claims; §1 lists 25 *(pre-retrieval)*

**What.** Spec §0/§2 says "全部 24 条". Counting §1: D1–D4 (4), P1–P3 (3),
T1–T6 (6), S1–S3 (3), X1–X5 (5), R1–R4 (4) = **25**.

**Decision.** Audit all 25. The count in the prose is treated as a
typo, not as an instruction to drop a row.

**Alternative rejected.** Guessing which row was meant to be excluded.
There is no principled way to pick, and dropping a row silently is
exactly the failure mode §5's "report incomplete" rule exists to prevent.

---

## DEC-02 — Repository layout coexists with the prior project *(pre-retrieval)*

**What.** `docs/` already contains the previous quadrant-mining project's
record (`PHASE0_LITERATURE.md`, `DEVIATIONS.md`, …).

**Decision.** The new deliverables are added alongside under the names
the spec fixes (`docs/PREREGISTRATION.md`, `DECISIONS.md`, `RESULTS.md`,
`CHANGELOG.md`, `out/`, `refs/`). Nothing from the prior project is
modified or deleted.

**Alternative rejected.** Nesting this round under a subdirectory. The
spec fixes the paths explicitly.

---

## DEC-03 — Preregistration committed before the first search *(pre-retrieval)*

**What.** Spec §5 requires the decision rule to be committed before any
retrieval begins.

**Decision.** `docs/PREREGISTRATION.md` and this file were written and
committed in a dedicated commit that contains **no results**, before any
web search, PubMed query, or full-text fetch was issued in this session.
The commit hash of that commit is the timestamp of record.

---

## DEC-04 — chrM coverage threshold fixed in advance *(pre-retrieval)*

**What.** Spec C.2 requires the chrM coverage threshold to be set from
the variant caller's documentation and written into DECISIONS before the
power table is computed.

**Decision.** Pre-registered at **median per-cell chrM coverage ≥ 20×**,
with informative variants defined by mgatk defaults (strand concordance
≥ 0.65, variance–mean ratio ≥ 0.01, detected in ≥ 5 cells). mgatk is the
standard mtscATAC-seq heteroplasmy caller; Lareau et al. 2021 operate at
~40–100× and low-coverage cells are treated as uncallable.

**Alternative considered.** Setting the bar at the caller's absolute
floor (~10×) to admit more datasets. Rejected: C.2 forbids lowering a
threshold to make an analysis runnable, and admitting near-floor coverage
would produce clone assignments whose error rate is the dominant term in
any compartment-sharing statistic.

---

## DEC-05 — What counts as "the same claim" when sources disagree *(pre-retrieval)*

**What.** Several rows (T1, T4, P3) are contested in the literature: the
question is not only how strong the evidence is but whether the two
sides are even asserting the same transition.

**Decision.** The claim text in spec §1 is authoritative for what is
being audited. Where a paper addresses a *neighbouring but different*
transition, it is recorded in `strongest_alternative_explanation` or
`counter_evidence`, not used to grade the row's primary evidence.

---
## DEC-06 — GSE302113 profiled; the spec's description of it needed two corrections

**What.** Spec C.1 nominates `GSE302113` (Liu et al., *Cancer Cell* 2026) as the
known starting point, described as "约 4 对配对病人" NSCLC mtscATAC-seq.

**Measured.** The series is 38 libraries / 218,715 cells over **10 donors**:
5 NSCLC (SU-L-001…005: lung tumour + lung normal + PBMC each) and 5 ovarian
(SU-O-001…005: tumour + PBMC, plus one omentum metastasis). Published as
Liu VV et al., *Cancer Cell* 2026;44(7):1509-1521.e4, PMID 42242233.

**Two corrections to the spec.** (1) It is 5 paired NSCLC donors, not ~4, and the
ovarian arm was not mentioned at all. (2) The compartment structure is better than
the spec assumed for X1: lung tumour **and** matched non-involved lung **and** blood
are present for all 5 NSCLC donors, and the NK-containing fraction (CD45+CD3−) was
separately sorted in both tissue compartments.

**chrM retention (the spec's explicit trap).** Every library retains chrM: the
fragments are aligned to a NUMT-masked reference (`hg38_v20-mtMask`,
cellranger-atac 2.0.0), and mgatk output is deposited per library
(`variant_stats.tsv.gz` over all 16,566 chrM positions, `cell_heteroplasmic_df.tsv.gz`
per cell). **All 38/38 libraries clear the preregistered 20× floor**; median chrM
coverage ranges 23.8× – 158.9×. Table: `out/T4_power_GSE302113_libraries.tsv`.

**What is NOT deposited: cell type annotations.** This is the binding constraint on
the power check, not coverage — see DEC-07.

---

## DEC-07 — First NK-typing attempt FAILED its own control; reported, not patched

**What.** With no deposited annotations, NK counts had to be derived. First attempt:
count fragments over marker gene bodies + 2 kb promoters (59 genes, 7 lineages) per
cell, convert to enrichment over background, z-score each lineage across cells, and
assign each cell the argmax lineage.

**The control.** GSE302113 contains its own positive/negative control pair:
GSM9096509 is **CD56+ sorted** PBMC and GSM9096510 is the matched **CD56−** fraction
from the same donor (SU-L-003). The CD56+ library must come out NK-high.

**Result: the control failed.** CD56+ = 2.9% NK, CD56− = 5.3% NK — a ratio of
**0.5×, i.e. backwards** — and 70% of cells were unclassifiable. The non-ambiguous
calls were near-uniform across the seven lineages (134/133/128/111/108/90), which is
the signature of noise rather than of biology.

**Diagnosis (measured, not guessed).** Median total fragments per cell is 7,538, but
the marker windows are 4.26 Mb = 0.14% of the genome, so a cell carries a **median of
4 NK-marker fragments and 4 T-marker fragments**. Per-cell argmax over counts of ~4 is
noise. The instrument was underpowered, and z-scoring guaranteed each lineage would
claim ~1/7 of the tail regardless of content.

**The panel itself is fine — the same data proves it in aggregate.** Pooling the same
counts to library level, CD56+ vs CD56− gives NK **+1.70 log2FC** with the eight
top-ranked genes all NK (SH2D1B +3.18, KIR2DL4 +3.05, NCR1 +3.04, KLRF1 +2.58,
GNLY +2.56, NKG7 +2.43, PRF1 +2.24, KLRD1 +2.12) and T/B markers depleted
(IL7R −1.90, CD4 −1.59, CD5 −1.51, CD3G −1.40, MS4A1 −1.27). Marker choice and
fragment counting are correct; only per-cell assignment was.

**Decision.** Do not tune thresholds until the control passes — that would be fitting
the classifier to the answer. Replace the instrument with the standard scATAC
approach (genome-tiled LSI → clustering → annotate *clusters*, where marker counts
aggregate over hundreds of cells and are reliable), and re-run the same CD56+/CD56−
control as the acceptance test. If the control fails again, the NK counts are
reported as underivable rather than reported anyway.

**Alternative rejected.** Estimating an NK *fraction* per library by bulk
deconvolution and multiplying by cell count. It would have produced a number for
every cell of the power table, but the preregistered threshold (≥30 NK cells per
donor per compartment) consumes a count of cells that can actually be carried into a
clone-sharing analysis, and G6 requires the summary quantity to be computed on the
same basis as the rule that consumes it. A fraction estimate is not that.

---

## DEC-08 — Second attempt also failed; the fix was the baseline, not the algorithm

**What.** Replacing per-cell argmax with the standard scATAC route (5 kb tiles →
TF-IDF → LSI → KMeans → annotate clusters) still failed the control:
CD56+ 23.2% NK vs CD56− 18.1% NK, ratio 1.3×.

**Diagnosis.** Same statistical error as DEC-07, one level up. Cluster annotation
z-scored each lineage's enrichment **across the clusters of its own library**. That
silently assumes every library is a mixture of lineages. For a CD56+ *sorted*
library most clusters really are NK, so centring on the library's own mean forces
about half of them below zero and re-labels them as something else. The procedure
could not report a homogeneous library as homogeneous.

**Fix.** Separate clustering from annotation, and z-score each lineage across the
**pooled cluster set from all libraries** instead. The pool genuinely spans all
lineages, so the baseline is no longer a function of the library being scored.
Enrichment is log-transformed before z-scoring, since enrichment is multiplicative
and a few very pure clusters would otherwise set the scale for their lineage.

**Result: control PASSES.** CD56+ = **62.1% NK**, CD56− = **0.0% NK** (that library
resolves as T-dominated, 1197/2109 cells, which is what depleting CD56 from PBMC
should leave). Acceptance test was fixed in the script before the run:
`NK_fraction(CD56+) > 0.5 AND > 3 × NK_fraction(CD56−)`.

**What this control does and does not license.** It shows the pipeline can tell
NK-rich from NK-free *blood* libraries at this depth. It does not validate NK calls
in tumour or lung tissue, where NK is rarer and ILC1/tissue-NK boundaries are
exactly what claims T1/T4/X1 are about. Tissue NK counts are therefore reported as
what they are — a power estimate for admitting or refusing a dataset cell, not a
biological result — and the residual 62.1% (rather than ~90%) is a reminder that
these counts are conservative-to-noisy, not exact.

---
