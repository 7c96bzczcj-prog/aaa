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

## DEC-09 — The mtDNA substrate for a cross-compartment question exists; NK counts are the open question

**What.** Before counting NK cells it is worth knowing whether GSE302113 has
informative mtDNA variants *shared between compartments within a donor* at all. A
dataset can pass a coverage threshold and still be useless for clone sharing if each
compartment's informative variants are disjoint.

**Measured** (preregistered mgatk filters, per donor by union across libraries,
chrM:307-314 excluded as the authors do):

| donor | compartments | union informative | present in >= 2 compartments |
|---|---|---|---|
| SU-L-001 | normal, tumour, PBMC | 199 | 54 |
| SU-L-002 | normal, tumour, PBMC | 642 | 178 |
| SU-L-003 | normal, tumour, PBMC | 215 | 32 |
| SU-L-004 | normal, tumour, PBMC | 232 | 52 |
| SU-L-005 | normal, tumour, PBMC | 1067 | 319 |
| SU-O-002 | tumour, PBMC | 632 | 43 |
| SU-O-004 | tumour, PBMC | 573 | 109 |
| SU-O-005 | omentum met, tumour, PBMC | 733 | 250 |
| SU-O-001 | tumour only | 204 | 0 |
| SU-O-003 | tumour only | 51 | 0 |

**Reading.** All 5 NSCLC donors carry the three compartments X1 needs and 32–319
cross-compartment informative variants each, so the barcode substrate is not the
limiting factor. SU-O-001 and SU-O-003 have a single compartment and cannot
contribute to any cross-compartment claim; they are excluded from those rows rather
than counted as donors.

**Caveat carried forward.** Variant sharing across compartments is necessary, not
sufficient: these counts are over *all* cells, not NK cells. Whether enough NK cells
carry them is exactly what the power table decides.

---

## DEC-10 — Early-termination gate does NOT fire; Phase A continues

**Rule** (preregistered §5.2): if more than 3 of the first six claims (D1, D2, D3, D4, P1, P2)
come back `fulltext_available = FALSE`, Phase A halts and a literature-access report is
delivered instead of a grading table.

**Result.** Of the first four returned (D1, D2, D3, D4): **zero** are FALSE; all four are
PARTIAL, meaning the principal primary source was read in full and one or more secondary or
corroborating sources were not. Under the preregistered wording only FALSE counts toward the
threshold, so the gate cannot fire on the remaining two. Phase A proceeds to all 25 claims.

**Note on PARTIAL.** PARTIAL is doing real work here and is not a fudge: for these claims the
literature is a chain of papers, and the rule that matters is G1 — the classification must rest
on a methods section actually read. Every PARTIAL row carries the methods quote it was graded
from. Rows where the *decisive* paper was unreadable are marked FALSE, not PARTIAL.

---

## DEC-11 — Remaining Phase A batches launched in parallel before the gate fully closed

**What.** Batches for T/P3, S/X and R claims were launched while the gate batch was still
finishing.

**Why this does not violate §5.2.** The rule's purpose is to prevent *delivering* a grading
table built on abstracts. Running audits concurrently does not deliver anything; had the gate
fired, every parallel row would have been discarded and the access report delivered instead.
The cost of being wrong was wasted compute, not a compromised table. With four of six gate rows
already back and none FALSE, the risk was already resolved when the last two batches launched.

---

## DEC-12 — Power check result, and what the zeros do and do not mean

**Acceptance control re-run on the final pooled annotation** (mandatory, because the pooled
baseline changes as libraries are added): CD56+ = **74.6% NK** (1948/2612), CD56− = **5.5%**
(117/2109), ratio **13.4×**. PASS, and better than the 2-library pool.

**Independent corroboration.** The typing reproduces, without being told, a qualitative claim
the source paper makes in prose: NK is scarce in most tumours, and one ovarian tumour
(SU-O-005) is the exception with substantial NK infiltration. Measured here: SU-O-005 tumour
**2,476 NK**, while SU-O-002/003/004 tumours return none at cluster resolution. That is a second
validation, on a different axis from the sorted control.

**Preregistered verdicts** (>= 30 NK per donor per compartment; >= 3 donors with two comparable
compartments). Full table: `out/T4_power.tsv`.

| question | compartments | donors passing | verdict |
|---|---|---|---|
| **X1** | adjacent normal lung vs tumour | **4** (SU-L-001, -002, -004, -005) | **PASS** |
| X2 (lung) | blood vs tumour | **4** (same) | **PASS** |
| X2 (ovarian) | blood vs tumour | 1 (SU-O-005) | **FAIL** |

Coverage passes everywhere (23.8×–158.9× median, floor 20×), and informative variants per donor
range 51–1067, so neither is limiting.

**The zeros are a detection floor, not a measurement of zero.** Annotation is at cluster level,
so an NK population smaller than roughly one cluster can be missed entirely. For the libraries
returning no NK the floor is ~324 cells (SU-L-003 tumour), ~331 (SU-O-002), ~124 (SU-O-003) and
~379/~614 (SU-O-004). Every one of those floors is **above** the 30-cell threshold the rule
tests. So those rows must be read as "no NK population detectable at cluster resolution", and
they are **not** evidence that fewer than 30 NK cells are present. The FAIL verdicts are
therefore conservative in the direction of *excluding* donors that might in fact be usable —
they may under-count available power, never over-count it.

**Why this does not weaken the X1 PASS.** The passing rows clear the threshold by one to three
orders of magnitude (tumour NK 281–509; adjacent-normal NK 394–6,971), so the verdict survives
even large classification error. The counts are a power estimate for admitting or refusing a
dataset cell, exactly as DEC-08 scoped them, and are not offered as biological quantities.

**Not fixed post hoc.** Raising cluster count to lower the floor would change the instrument
after seeing which donors fail, which is precisely what C.2 forbids. The floor is reported
instead.

---

## DEC-13 — Rule G8: where the spec's initial calls were wrong, and how

Phase A ran 25 audits plus 25 independent adversarial re-examinations (48 agents
across four batches, 0 errors, ~1,690 tool calls). The adversarial pass changed
**10 of 25** classifications, in both directions, so it was not a rubber stamp.

**Citations first.** Every citation the spec flagged as possibly misremembered turned
out to be **real and correctly attributed to a real paper**:

| spec citation | verdict |
|---|---|
| Picant, *Nat Commun* 2025 (T6) | real — PMID 40610398 |
| Serger 2026 snRNA+snATAC (X4) | real — *Sci Immunol* 2026, PMID 42247486 |
| Schmid/Wiedemann ATAC (R3) | real — bioRxiv 2026.02.11.705354 |
| Barahona/Yokoyama (T2) | real — *eLife* 2026, PMID 42417504 |
| Gamliel 2018 (R4) | real — *Immunity* 2018 |

**But three of those four carried the wrong evidence class**, which is precisely what
G1 exists to catch:

- **T6** — spec said TRANSFER. There is no adoptive transfer anywhere in Picant et al.;
  a search of the archived full text returns zero hits for "adoptive transfer",
  "congenic", "CD45.1", "NSG" or "NOD scid". Human in vitro throughout. → `INF_INVITRO`.
- **R3** — spec said OBSERVATION. TGF-β genuinely was withdrawn and chromatin and
  function genuinely were measured separately, exactly as the spec described, but the
  entire withdrawal experiment is in vitro. → `INF_INVITRO`, no `OBS_*` class.
- **T2** — classification upheld (`OBS_TRANSFER`), attribution corrected: the author is
  Josselyn D. Barahona (Yokoyama lab), not "Barahona Ponce", who is a different
  researcher working on gallbladder cancer genetics.
- **X4** — spec's call confirmed on every element; the audit sharpened it to
  `INF_TRAJECTORY` under G2.

**The largest single overturn is S1.** The spec called it "pseudotime only = pure
inference". The audit found something stronger: **Vento-Tormo 2018, the paper that
defines dNK1/2/3, makes no ordering claim at all** — its methods state verbatim that
"only cells that were identified as trophoblast were considered for trajectory
analysis". The ordering is a later accretion, and the three papers that do assert one
give three **mutually incompatible topologies** (Wang 2021: dNK1→dNK2→dNK3 with dNK1
*immature*; Huhn 2020: dNK3→dNK2→dNK1 with dNK1 the *mature* endpoint — opposite
polarity for the same subset; Guo 2021: three parallel sibling branches with no
ordering among dNK1/2/3 at all).

**Two claims changed kind, not just class.**

- **P3** — the issue is not a false transition but an **entity misassignment**. Paust
  2010 gated CD45+NK1.1+CD3− with no CD49a and no DX5, a gate that pools conventional
  NK with liver ILC1. Wang 2018, sorting cNK, IL-7Rα− LrNK and IL-7Rα+ ILC1 from the
  same donors, found only the ILC1 fraction conferred hapten recall. The memory cell
  may never have been an NK cell, and so never "became" anything.
- **X5** — downgraded `OBS_TRANSFER` → `INF_KINETIC`, consistent with the spec's own
  suspicion that the end state is undefined.

**Where an audit was wrong in the *negative* direction.** X1's audit asserted a
universal negative — that no parabiosis, fate map, photoconversion or barcode exists
"anywhere in this or any other study for this transition". The adversarial pass
falsified it by producing Dadi et al., *Cell* 2016. G3 is usually invoked against
overclaiming absence of evidence; here it caught an overconfident absence claim inside
our own audit. That is the rule working in the direction it is least often applied.

---

## DEC-14 — P1 occupancy: found independently, and it also tests my own threshold

**What.** While building the C.1 dataset inventory I found, independently of the Phase B
occupancy agent, that the exact instrument P1 would need has already been applied to P1's
population: **Rückert T, Lareau CA, et al., "Clonal expansion and epigenetic inheritance
of long-lasting NK cell memory", *Nat Immunol* 2022, PMID 36289449** (GEO GSE197037 /
GSE197008; code at github.com/timorueckert/Clonal_NK).

**Methods, read not skimmed.** Human, 4 HCMV+ and 3 HCMV− donors, mtscATAC-seq/ASAP-seq
with ADT and HTO capture; mgatk in tenx mode; informative variants filtered at strand
concordance > 65%, variance-to-mean ratio > 0.01, detected in ≥ 3 cells; clonotypes built
by clustering a neighbourhood graph on mtDNA mutation frequency. Result: substantial clonal
expansion of adaptive NK cells, with clonotypes associated not only with the adaptive
compartment as a whole but with specific adaptive subclusters, and mutations "specifically
enriched in the adaptive NK cell compartment" in HCMV+ donors.

**Why this is recorded here rather than left to the agent.** B4 is the criterion most
likely to be answered from memory, and a wrong FALSE there is what corrupts a shortlist
(the preregistration says so explicitly). Finding the occupant independently gives a check
on the agent's answer rather than a substitute for it; the two are compared in RESULTS.

**Scope caution before this is used to kill P1.** Rückert establishes clonal expansion and
persistence *within* the adaptive compartment. P1 as written asserts a *descent* relation,
cNK → adaptive NKG2C+. Whether adaptive-restricted clonotypes settle that descent question,
or only the expansion question, is a real distinction and is left to the Phase B occupancy
verdict rather than resolved by assertion here.

**It also tests my own preregistered threshold, and I am not moving it.** Rückert operate at
a **median chrM coverage of 11–20×**, below the 20× floor I preregistered in DEC-04, arguing
that sensitivity for high-heteroplasmy mutations is relatively coverage-independent. So a
published, peer-reviewed, field-defining application of exactly this method would be
**refused by my own admissibility rule**. That is worth stating plainly: my threshold is
conservative relative to accepted practice. It stays as preregistered — loosening it now,
after seeing that it excludes a dataset I would like to use, is precisely the post-hoc
adjustment C.2 forbids. The consequence is recorded as a limitation, not repaired.

---

## DEC-15 — Normalising the adversarial verdicts; two conditional cases adjudicated

**Problem.** The Phase A adversarial schema left `revised_evidence_type_strongest`
as a free string. Several reviewers returned a reasoned paragraph instead of a
token, so the field was not machine-readable.

**Rule applied** (`scripts/normalize_verdicts.py`):
1. If the reviewer led with a class token, lift it.
2. If the reviewer challenged the class but supplied no token, **keep the audit
   class** rather than guess.
3. If the verdict was genuinely conditional on a reading of the claim, adjudicate
   explicitly under DEC-05 (spec §1 claim text is authoritative) and record why.

**Two conditional cases:**

- **P1** — reviewer gave `INF_KINETIC` for the human NKG2C/HCMV claim and
  `OBS_DIRECT` for the mouse Ly49H/MCMV homolog. NKG2C and HCMV are human-specific;
  the mouse work uses a different receptor, virus and species. Neighbouring but
  different transition → `INF_KINETIC`.
- **X1** — reviewer gave `OBS_TRANSFER` under an operational trNK definition and
  `INF_MARKER` under a strict Eomes+ conventional-NK definition, where the same
  figure becomes counter-evidence. Resolved against the registry's own ontology:
  T4 ("NK → ILC1") and P3 (NK-or-ILC1 identity dispute) both presuppose ILC1 ≠ NK,
  so trNK means NK-lineage tissue-resident NK → `INF_MARKER`.

**Alternative rejected.** Re-running the adversarial agents with a constrained
enum. It would have produced cleaner fields but would have re-rolled verdicts
already returned, after I had seen which ones were inconvenient. Normalising the
existing output with a stated rule keeps the record auditable.

---

## DEC-16 — MEDIUM-confidence FALSE on occupancy is surfaced, not reclassified

**What.** Both shortlist entries (S2, S3) carry `already_traced = FALSE` at
`search_confidence = MEDIUM`. The preregistration told B4 agents to prefer
`NOT_SEARCHED` over a thin `FALSE`.

**Decision.** Apply §5 criterion 4 exactly as written — a `FALSE` admits the claim
without the `occupancy_unverified` flag — and add `occupancy_search_confidence` as
a column in T5 plus an explicit paragraph in RESULTS §1.

**Why not reclassify to NOT_SEARCHED.** That would be me overriding the agent's own
assessment of its search in the direction that weakens my only two shortlist
entries, after seeing the result. The preregistered rule keys on the recorded
value; the confidence is reported alongside so a reader can discount it. Changing
the value would be a post-hoc edit of §5 by the back door.

---

## DEC-17 — v1.1 re-grade of X1 and S1 under the bidirectional criterion

**Trigger.** The spec's author judged the v1.0 `load_bearing` definition too
narrow: it asked only what presupposes a claim is TRUE. Amendment recorded in
`docs/PREREGISTRATION.md` §A1 and `docs/CHANGELOG.md` before the re-grade was run.

**X1: MEDIUM → HIGH, carried by the negation side.** The v1.0 finding is not
withdrawn — Horowitz 2026's ctrNK ACT platform really is origin-agnostic and really
does survive X1 being false. What v1.0 could not see is that four primary papers
plus an engineered-cell therapy presuppose X1's negation: Serger 2026
(PMID 42247486), Gao 2017 (PMID 28759001), Cortez 2017 (PMID 28759002), Horowitz
2026's *mechanism* as distinct from its platform (PMID 42090477), and
DN-TGF-βRII NK constructs (PMID 28109751, PMID 36524207). Full dependency list with
directions: `phaseA/load_bearing_v11.json`.

**Consequence: X1 enters the shortlist**, having already passed criteria 2, 3 and 4.

**Registry-internal note.** X4 (Serger 2026) is itself one of X1's negation-side
dependencies, so X1 and X4 are mutually constraining. X4 stays `INF_TRAJECTORY`
under G2 regardless — the bidirectional amendment touches criterion 1 only.

**S1: MEDIUM → HIGH, status unchanged.** S1 turns out to carry dependencies in
*both* directions — Guo 2021's RPL disease model requires an ordering, while the
larger compositional literature (representative: Huhn 2020) presupposes none, and
Vento-Tormo 2018 itself frames dNK1/2/3 as co-existing "states". Under the union
that is HIGH. **It changes nothing**: S1 fails criterion 3 on both decidability
gates independently of load-bearing, because the asserted ordering is a directed
arrow over a reversible, milieu-instructed state and clonal structure cannot order
a reversible state. Re-graded for the register only, and recorded as such so that
nobody later reads HIGH as a near-miss.

**Scope discipline.** Applied to X1 and S1 only, on instruction. The other 23 rows
keep their v1.0 grade and are flagged `load_bearing_basis = v1.0_positive_only` in
T5 so a stale grade is never mistaken for a bidirectional one. Re-grading all 25
under a changed criterion after seeing the v1.0 outcome would let the amendment
reshape the result set retrospectively; that is deliberately not done.

**Not touched in the same motion:** the 20× chrM floor, known to be uncalibrated
(DEC-14), is left exactly as preregistered. Recalibrating a threshold in the same
change that admits a claim it bears on is the failure mode the preregistration
exists to prevent, whichever direction the recalibration would go.

---

## DEC-18 — A decidual-NK fate-mapping paper found during the re-run; S2/S3 unchanged

**What.** While searching S1's negation side I found Zhang et al. (?), "T-bet Fate
Mapping Reveals Gestational Stage-Specific Transcriptional Adaptation of Decidual
NK", 2026, PMID 41782473 / PMC12961421 — a genuine fate-mapping study of uterine
and decidual NK.

**Methods, read not skimmed.** `Rosa26^RFP × Tbx21^Cre` T-bet fate-mapping **mouse**
model, tracking NK in uterus, decidua and placenta through pregnancy; flow
cytometry, bulk RNA-seq of fate-mapped cells, scRNA-seq of CD45⁺Lin⁻ cells at mid
and late gestation.

**Effect on S2/S3: none, and I checked rather than assumed.** S2 is a *human*
subset-correspondence claim (uNK1/2/3 → dNK1/2/3); those subsets are defined by
human scRNA-seq and have no mouse counterpart. S3 is the *cycle-to-pregnancy*
transition; mice do not menstruate, so there is no cycling-endometrium-to-decidua
transition to trace, and a T-bet fate map labels a lineage rather than a temporal
transition from a pre-pregnancy compartment. Neither `already_traced` value changes.

**But it is recorded, for two reasons.** It is now the nearest occupant on
method-plus-compartment for both rows, and it reinforces that S2 and S3 rest on
`occupancy_search_confidence = MEDIUM` (RESULTS §1) — this paper is recent enough
that the occupancy agents may not have weighted it. It is a reason to re-run B4 on
S2/S3 with the mouse fate-mapping literature explicitly in scope before either is
acted on, not a reason to change the value now.

---

## DEC-19 — Only X1 could ever have moved; T1 re-graded; five rows opened and left ungraded

**Checked before spending effort, and it changes what the re-run is for.** v1.1 is
a *union* (positive OR negation), so it can only raise a grade — verified
mechanically. The criterion-defined re-grade set is therefore exactly the 14 rows
graded below HIGH under v1.0, and it is complete and unbiased by construction, as
the spec author argued.

**But of those 14, X1 is the only row whose shortlist membership a load-bearing
re-grade could change.** Every other row is blocked by an independent criterion:
T1/T3/T5 by `already_traced = TRUE`; R2/T2/X2 because their evidence is already
`OBS_*`; D3/R1/R3/S1/T3/T5/T6/X3/X4 on decidability. Encoded as a check in
`scripts/verify_results.py`. The remaining re-grades are register consistency,
not membership.

**T1: MEDIUM → HIGH, carried by the negation side, status unchanged.** The
negation side is unusually clean — Sojka et al., eLife 2014 (PMID 24714492)
asserts T1's negation *in its title*: "Tissue-resident natural killer (NK) cells
are cell lineages **distinct from** thymic and conventional splenic NK cells".
Klose 2014 (PMID 24725403) supplies the separate-progenitor architecture and
Nixon 2022 (PMID 35394814) reports the ILC1 lineage "did not interconvert with NK
cells". The ILC classification framework presupposes T1 false. T1 still fails
criterion 4, so nothing moves; the re-grade removes an internal inconsistency —
a MEDIUM sitting next to X1's HIGH on an adjacent claim.

**A structural finding worth keeping.** T1 and X1 are near-complementary, and the
DN-TGFBR2 engineered-NK construct demonstrates it: the *same* therapeutic premise
requires **T1 true** (tumour TGF-β converts infused cNK) and **X1 false**
(intratumoural trNK are infiltrating-and-converted rather than pre-resident). One
construct, opposite directions on two adjacent claims.

**T3, T5, T6, R1, R3: opened, NOT graded.** Targeted searches for their negation
side returned reviews rather than primary sources, so no negation dependency was
*established*. Under G3 that is NOT_SEARCHED-equivalent, not NONE; under G7 a
grade may not be raised on an unverified dependency. Their v1.0 grades stand and
stay flagged `v1.0_positive_only`. Recorded as unfinished in
`phaseA/load_bearing_v11.json` under `_unfinished` rather than left to look
complete. D3, R2, T2, X2, X3 and X4 are also part of the criterion-defined set and
were never opened.

---

## DEC-20 — X1 design memo kept outside the screening record

**What.** `docs/X1_DESIGN.md` sets out how X1 would actually be adjudicated.

**Why it is a separate file.** Spec §0 fixed "no new directions" as a non-goal of
this round, and RESULTS.md honours that. The design question was asked separately,
after the screen closed. Putting the design into RESULTS.md would retroactively
violate a non-goal the preregistration fixed before any retrieval; keeping it in
its own document leaves the screening record clean and lets the memo be revised
without touching a preregistered result.

**Nothing in the memo feeds back.** No grade, table row or verdict is derived from
it.

---

## DEC-21 — Tier 1 stopped before computing: an ascertainment confound, not a power failure

**What.** The X1 Tier-1 re-analysis was preregistered (`4aaa154`) and then stopped
before any statistic was computed. Full record: `docs/X1_TIER1_RESULT.md`.

**Why.** The preregistered primary statistic asks whether tumour–NILT shared clones
are under-represented in blood. GSE302113 deposits **per-library** heteroplasmy
matrices, not the donor-union matrix the authors' own methods describe building. A
variant can only be evaluated in a compartment where it was deposited, and it is
deposited only if it passed mgatk's ≥5-cell filter there — so the variant set
evaluable across all three compartments is **ascertained to be blood-present**:
4.5% of SU-L-001's 220 tri-compartment variants and **0.0%** of SU-L-005's 507 have
zero blood carriers. The statistic would answer its own question by construction.

**Decision: report the stop; run no substitute.** Tier-1 preregistration §5 fixes
that no threshold is lowered to produce a runnable result. A weaker pairwise
statistic was available and inherits the same bias, so it was not run either.

**This is the GSE221064 trap in a new place.** The spec warned that archives strip
chrM. Here chrM is fully retained and clears every threshold — what is stripped is
the *variant union*, one level up. Checking that files exist and clear a coverage
floor is not the same as checking that their ascertainment can carry the question.

**Correction issued.** `docs/X1_DESIGN.md` §3 called Tier 1 "nearly free" because
the data was "already on disk". That was wrong for exactly this reason. The section
is corrected in place with the error left visible rather than rewritten.

**Routes forward, costed.** (A) Ask the authors for the union matrices or their
Mitotrek clone assignments — Mitotrek is public, the matrices are not. (B) Rebuild
from SRA FASTQs (per-sample SRX accessions exist) with cellranger-atac + mgatk over
38 libraries — TB-scale, days of compute. (C) Abandon Tier 1. Tier 1 is
**expensive, not blocked**.

**One published sentence recorded because it cuts against X1**, not despite it:
Liu 2026 groups NK with monocytes/macrophages/DCs as showing "high levels of clone
sharing, suggesting that they originated from recent hematopoietic output without
substantial clonal bottlenecks prior to tissue infiltration" — the blood-origin
side. Qualitative, never given NK-specific numbers, attached to heatmaps whose NK
entries were never separately reported.

---
