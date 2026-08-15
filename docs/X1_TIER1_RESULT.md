# X1 Tier 1 — STOPPED before computing, on an ascertainment confound

**Preregistered**: `docs/PREREGISTRATION_TIER1.md`, committed at `4aaa154` before
any Tier-1 statistic was computed. Nothing in that file was amended.

**Outcome: the preregistered statistic is not executable on the deposited data.**
Not a power failure — an **ascertainment** failure, which the preregistration did
not anticipate. No substitute statistic was run.

---

## 1. Step 1 as fixed: recover the published Figure 2F value first

**Attempted, and it is not recoverable from the published record.**

- Fig 2F legend, verbatim: *"Heatmap showing the fraction of all cell pairs
  belonging to the same clone and consisting of a lung tumor cell type and a NILT
  cell type."*
- The three-way comparison is in fact **all three plotted**: Fig 2E (within
  tumour / within NILT), Fig 2F (tumour × NILT), Fig S6C (PBMC × solid tissue).
- The only supplementary file deposited with the preprint is Table S1, which is
  **patient metadata** (donor ID, diagnosis, age, stage, PD-L1). No source data
  for the heatmaps; the panels are images.
- Result: the tumour-NK × NILT-NK cell exists in three published figures and its
  numeric value is **not recoverable** without the authors' underlying matrices.

**One published sentence does bear on X1, and it cuts against it.** The text groups
NK with myeloid: *"innate immune cells, including monocytes, macrophages, DCs, and
NK cells, consistently exhibited high levels of clone sharing, suggesting that they
originated from recent hematopoietic output without substantial clonal bottlenecks
prior to tissue infiltration."* That is the authors' own reading and it is the
blood-origin side. It is qualitative, NK-specific numbers are never given, and it
is the interpretation attached to heatmaps whose NK entries were never separately
reported — but it is on the record and should not be ignored because it is
inconvenient.

---

## 2. Why the preregistered statistic cannot be computed

The primary statistic asks whether clones shared between tumour and NILT are
**under-represented in blood**. That requires evaluating, for a given mtDNA
variant, whether blood cells carry it. Three structural facts make that impossible
on the deposited files.

**(a) What is deposited is per-library, not the donor union.** The authors' own
methods describe building "a cell by variant heteroplasmy matrix combining all
samples for each donor" — the union procedure. That union matrix is **not
deposited**. What is deposited is one heteroplasmy matrix per library, with that
library's own mgatk-selected variants as columns. Measured column overlap between
libraries of the same donor:

| donor | libraries | union of variants | intersection | ratio |
|---|---|---|---|---|
| SU-L-001 | 5 | 1,516 | 119 | 0.08 |
| SU-L-002 | 5 | 2,183 | 73 | 0.03 |
| SU-L-004 | 5 | 1,146 | 35 | 0.03 |
| SU-L-005 | 7 | 3,388 | 89 | 0.03 |

Even two libraries of the *same compartment* in SU-L-005 share only 0.06–0.23 of
their columns. Per-library selection, confirmed.

**(b) A variant can only be evaluated where it was deposited.** For a variant
selected in tumour and NILT but not in blood, the deposited data contains no
per-cell heteroplasmy in blood, so "is this clone absent from blood?" is
unanswerable — absence of the column is not absence of the variant.

**(c) Restricting to variants evaluable in all three compartments creates exactly
the confound being measured.** A variant appears in a compartment's file only if it
passed mgatk's filter there, which requires ≥5 confidently-detected cells in that
compartment. So the usable variant set is **ascertained to be blood-present**:

| donor | variants evaluable in all 3 compartments | carried by **zero** blood cells (≥0.07) |
|---|---|---|
| SU-L-001 | 220 | 10 (**4.5%**) |
| SU-L-005 | 507 | 0 (**0.0%**) |

The statistic would return "tumour–NILT shared clones are present in blood" for
**95.5–100% of the usable variants by construction of the variant set**, not by
biology. Running it would produce a confident, publishable-looking, meaningless
number pointing away from X1.

The secondary pairwise statistic inherits the same bias and was not run either.

---

## 3. What this does and does not support

*Does not support*: that X1 is false, or that tumour NK are blood-derived. Nothing
biological was measured.

*Does not support*: that GSE302113 is unusable for X1. The **power** check still
passes — 4 donors, ≥30 NK in both compartments, 23.8–158.9× chrM coverage, 32–319
informative variants shared across ≥2 compartments per donor (DEC-09). The
limitation is in what was **deposited**, not in what was **sequenced**.

*Supports*: that the cheap tier is **not cheap**. See §4.

---

## 4. Three routes, with honest costs

| route | what it needs | cost |
|---|---|---|
| **A. Ask the authors** | the donor-union heteroplasmy matrices, or the Mitotrek clone assignments they already computed | one email; Mitotrek is public at `github.com/vincent6liu/mitotrek`, the matrices are not |
| **B. Rebuild from raw** | FASTQs are in SRA (per-sample SRX accessions, e.g. `SRX29602347`); re-run cellranger-atac + mgatk over 38 libraries, then Mitotrek | large — TB-scale download, days of compute. Technically unblocked, not cheap |
| **C. Abandon Tier 1** | — | loses the only human-side evidence, leaving X1 entirely on Tier 2 |
| **D. Switch dataset** | a different mtDNA dataset with the right design | **CLOSED — no alternative exists.** See §6 |
| **E. Within-compartment statistic** | compare NK clone structure *inside* each compartment, using only that compartment's own variants | escapes the confound, but is a hint not a verdict. Not preregistered; would need its own preregistration. See §7 |

Route A is obviously first. Route B is what makes the honest statement possible:
**Tier 1 is not blocked, it is expensive** — which is a materially different claim
from the one in `docs/X1_DESIGN.md` §3, where Tier 1 was described as "nearly
free" because the data was "already on disk". That description was wrong, and it
was wrong because I checked that the files existed and cleared the power
thresholds without checking that their *variant ascertainment* could carry the
question. Corrected here rather than quietly.

---

## 5. What was actually gained

1. The preregistration did its job: the confound was found **before** a number was
   produced, so there is no result to retract.
2. A concrete, checkable statement about the published record: the tumour-NK ×
   NILT-NK clone-sharing value exists in three figures of Liu 2026 and has never
   been quantified or interpreted, and its numeric value is not in the public
   supplementary material.
3. X1's cost estimate is corrected upward. The claim survives the screen unchanged
   — `INF_MARKER`, load-bearing HIGH, decidable, both answers actionable,
   `already_traced = FALSE` — but its cheap human tier requires either the
   authors' cooperation or a full reprocessing run.


---

## 6. Route D — checked, and it closes

The T3 dataset table classified candidates at the **file** level ("fragments or
mgatk output present → chrM present"). That is precisely the check this episode
showed to be the wrong one, so the candidates' status was **unknown, not
second-best**. Re-audited at the variant level
(`out/T3b_datasets_variant_level.tsv`, `scripts/routeD_variant_level.py`):

| dataset | mgatk object deposited | design | verdict |
|---|---|---|---|
| **GSE302113** (Liu 2026) | **no** — derived per-library summaries only | tumour + non-involved organ + blood ✓ | the failure documented above |
| **GSE197037 / GSE197008** (Rückert 2022) | **yes** — `*_mgatk.rds.gz` per library | **PBMC only** ✗ | no ascertainment confound, but no tissue compartment to compare |
| GSE148508 | no | CD34 + PBMC | fragments only |
| GSE173936 | no | mtDNA-disorder samples | fragments only |
| GSE216911 | no | cultured T cells | nothing usable |
| GSE276283 | filelist unretrievable | glioma myeloid | unknown |

**Route D closes: no alternative dataset exists.** The only candidate with the
right compartment design is the one with the lossy deposit, and the only candidate
with the un-lossy deposit has no tissue compartment. That strengthens Route A
rather than substituting for it.

**Fragments cannot substitute.** A fragments file records intervals and a barcode,
not alleles. mtDNA variants are not recoverable from fragments at any effort, so
"fragments deposited" is not partial credit for this question — it is zero.

**A finding that makes the Route A ask easy to grant.** Rückert deposits the mgatk
object, which is the un-lossy form. Depositing it is therefore normal practice in
this exact subfield, and the request is for an intermediate the published pipeline
already produces.

**One correction to this re-audit, recorded because it is the same error class.**
The first pass reported GSE197037 as depositing nothing — the filename regex did
not match `*_mgatk.rds.gz`. Caught by printing the actual filenames instead of
trusting the pattern. Fixed in the script with the reason in a comment.

---

## 7. Route E — noted, not run, and it would need its own preregistration

A within-compartment statistic escapes the confound: compare the clone structure
of NK **inside** tumour, inside NILT and inside blood, each using only that
compartment's own variants, so no cross-compartment variant selection enters.

- tumour NK showing a **blood-like, flat** clonal structure → leans recruitment
- tumour NK showing **resident-like local expansion** → leans in-situ origin

**Why it is not run here.** It is a hint, not a verdict — it measures structure,
not ancestry, and structure has other causes (a tumour can flatten a resident
pool's clonality by recruiting on top of it). It also demands the same sampling
discipline as the original: diversity and expansion estimators are extremely
sensitive to cell number, and tumour carries only 281–509 NK against up to 7,713
in blood, so it must be downsampled to equal n per donor and reported per donor.

Running it would require preregistering its estimator, its downsampling scheme and
its interpretation rule **before** computation, exactly as Tier 1 was. It is
recorded here as an available option, not executed.
