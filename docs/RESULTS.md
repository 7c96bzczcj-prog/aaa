# RESULTS — NK lineage-assertion evidence screen v1.0

Preregistration committed at `9692011`, **before any retrieval**. Decision rule
unchanged since. Every conclusion below carries its `out/T1_evidence_audit.tsv`
row, evidence class, n with its unit, and **what it does not support**.

**Scale.** 25 claims. Phase A = 25 audits + 25 independent adversarial
re-examinations (48 agents, 0 errors, ~1,690 tool calls). Phase B = 17 claims ×
(decidability + occupancy) = 34 agents, 0 errors. Phase C = 38 mtscATAC libraries
downloaded and parsed. Primary-source texts archived under `refs/`.

---

## 0. Headline

**The shortlist is non-empty and contains two claims: S2 and S3.**

Per the preregistration (§5.1), this does **not** mean they are worth doing. Moat,
feasibility and time cost are a separate layer of judgement, explicitly outside
this round. This round says only: their evidence is inferential, lineage can
decide them, both answers change what you do next, and nobody has done them.

**The user-designated highest-priority claim, X1, is not on the shortlist.** It
fails on load-bearing, not on decidability or occupancy, and it fails after
passing the Phase C power check. Details in §3.

---

## 1. The shortlist (out/T5_shortlist.tsv)

| claim | T1 row | evidence | n (unit) | LB | decidable | both actionable | traced |
|---|---|---|---|---|---|---|---|
| **S2** endometrial uNK1/2/3 → decidual dNK1/2/3 | 14 | `INF_MARKER` | see below (donors/women) | HIGH | TRUE | TRUE | FALSE |
| **S3** eNK → dNK across the cycle-to-pregnancy transition | 15 | `INF_INVITRO` | 26 (women) | HIGH | TRUE | TRUE | FALSE |

**S2** — the uNK→dNK subset correspondence rests on cluster-label matching across
two separate datasets, i.e. marker similarity between clusters computed in
different studies. n is donors (individual women), not cells; the cell counts in
those atlases are not independent replicates (G5).
*Does not support*: that the correspondence is wrong. It supports only that the
correspondence has never been established by anything stronger than shared
markers.

**S3** — the cycle-to-pregnancy transition, primary source Feyaerts 2024, n = 26
women, and critically **0 endometrium–decidua pairs from the same woman**.
*Does not support*: that decidual NK are not carried over from endometrium. It
supports only that no human dataset pairs the two timepoints in one individual,
so the transition has never been observed in the same person.

**Both carry `occupancy_search_confidence = MEDIUM`, not HIGH.** The
preregistration instructed agents to record `NOT_SEARCHED` rather than a thin
`FALSE`; both returned `FALSE` but at medium confidence. Under §5 criterion 4 a
`FALSE` admits them without the `occupancy_unverified` flag, so that is how they
are recorded — but the confidence is surfaced here rather than smoothed over,
because a `FALSE` that should have been `NOT_SEARCHED` is exactly what corrupts a
shortlist.

**A structural observation, not a value judgement.** S2 and S3 are two framings of
one anatomical transition (endometrium → decidua): S2 asks whether the *subsets*
correspond, S3 whether the *cells* carry over. They would be adjudicated by
overlapping material, and one compartment they need — the same woman before
conception and while pregnant — is separated by months by construction, not by
logistics.

---

## 2. Why the other 23 failed (out/T5_shortlist.tsv, `failed_criteria`)

| failure mode | claims | n |
|---|---|---|
| criterion 2 — evidence is already `OBS_*` | D1, D2, D4, P3, R2, T2, T4, X2 | 8 |
| criterion 3 — not lineage-decidable, or both answers identical | D3, P2, R1, R3, R4, S1, T3, T5, T6, X3, X4, X5 | 12 |
| criterion 4 — already traced | P1, T1, T3, T5 | 4 |
| criterion 1 — load-bearing below HIGH | D3, R1, R2, R3, S1, T1, T2, T3, T5, T6, X1, X2, X3, X4 | 14 |

(Claims can fail several criteria; the table lists all.)

**B2 was the decisive filter, as the spec predicted.** Twelve of 17 Phase B claims
died there. The recurring reason is that the claim is functional rather than
ancestral: P2, R1, R3, T3, T5 and T6 all ask whether a state can be induced or
restored, and lineage does not measure function.

**One claim failed B2 for the opposite reason** — D3 (CD56^bright^→CD56^dim^),
where the reviewer judged the branches asymmetric rather than identical: a NO
would overturn the textbook arrow, void a live cell-manufacturing rationale
(iPSC/CD34-derived NK products are infused on the assumption they mature to
CD56^dim^ in the patient), and turn GATA2-deficiency patients into a test bed;
a YES would confirm the consensus and change nobody's protocol.

---

## 3. X1 — the priority claim (T1 row 22)

**Evidence** `INF_MARKER`. **n = 10 patients** (5 NSCLC with matched tumour +
non-involved lung + blood; 5 ovarian, 3 with matched blood) — donors, not cells.
**Load-bearing MEDIUM.** **`already_traced = FALSE` at HIGH confidence.**
**Lineage-decidable TRUE, both answers actionable TRUE.** **`mouse_only = TRUE`.**

X1 fails the shortlist on **criterion 1 alone**. Three things are worth separating.

**(a) Its own audit was wrong in the negative direction.** The audit asserted that
"no pre-tumour timepoint, no fate mapping, no photoconversion, no parabiosis, no
intravascular labelling, no exogenous barcode exists anywhere in this or any other
study for this transition." The adversarial pass falsified that by producing
Dadi et al., *Cell* 2016 (PMID 26806130) — congenic CD45.1/CD45.2 parabiosis in
tumour-bearing MMTV-PyMT mice with pre-transformation enumeration of the same
population. Rule G3 is normally invoked against claiming absence of evidence; here
it caught an overconfident absence claim inside our own audit.

**(b) The class turned on a nomenclature fork, resolved against the registry.**
Under an operational trNK definition (CD49a⁺CD103⁺NK1.1⁺), Dadi's population is
the claim's subject and the class upgrades to `OBS_TRANSFER`. Under a strict
Eomes⁺ conventional-NK definition, the *same figure* is counter-evidence, because
cNK showed **high** non-host chimerism — conventional NK in tumour-bearing tissue
are continuously exchanged with the circulation. This registry's own ontology
treats ILC1 as distinct from NK (claim T4 is "NK → ILC1"; claim P3 turns on
whether liver memory cells are NK or ILC1), so trNK here means NK-lineage
tissue-resident NK, and Dadi's population — explicitly Nfil3-independent and
"distinct from conventional NK cells" — addresses a neighbouring transition.
X1 keeps `INF_MARKER` (DEC-05, DEC-15).

**(c) Load-bearing HIGH → MEDIUM on a specific argument.** The headline dependency
was Horowitz 2026's ctrNK adoptive-cell-therapy platform. That platform
differentiates ctrNK **ex vivo from peripheral NK**, so it is origin-agnostic and
survives X1 being false. A dependency that survives the claim being false is not
evidence of HIGH load-bearing. What does survive is interpretive: the CD49a⁺ NK
prognostic literature in HCC reads causality off a residency pedigree, and the
choice between preserve-the-resident-pool and promote-recruitment therapeutic
levers depends on the answer.

**What X1's failure does not support.** It does **not** mean the question is
settled, uninteresting, or unanswerable. It means that on this registry's
preregistered criteria the downstream structure resting on it is interpretive
rather than interventional. Under a different load-bearing rubric it would qualify.

---

## 4. Phase C — X1 was run regardless, and it passes power

Run unconditionally per §4. Dataset `GSE302113` (Liu VV et al., *Cancer Cell*
2026, PMID 42242233). Tables: `out/T3_datasets.tsv`, `out/T4_power.tsv`.

**Corrections to the spec's description.** It is **5 paired NSCLC donors, not ~4**,
plus an unmentioned 5-donor ovarian arm; 38 libraries, 218,715 cells.

**chrM is retained in all 38/38 libraries** — NUMT-masked reference
(`hg38_v20-mtMask`, cellranger-atac 2.0.0), per-library mgatk output deposited.
Median chrM coverage **23.8×–158.9×**; every library clears the preregistered 20×
floor. The GSE221064 trap does not recur here.

**The binding constraint was that no cell-type annotations were deposited**, so NK
counts had to be derived. Two attempts failed the dataset's own CD56⁺/CD56⁻ sorted
control and are recorded as failures rather than tuned (DEC-07, DEC-08). The third
passes: **CD56⁺ 74.6% NK vs CD56⁻ 5.5%, ratio 13.4×**, against an acceptance test
fixed before the run.

| question | compartments | donors ≥30 NK in both | verdict |
|---|---|---|---|
| **X1** | adjacent normal lung vs tumour | **4** (SU-L-001/002/004/005) | **PASS** |
| X2 | blood vs tumour (lung) | 4 | PASS |
| X2 | blood vs tumour (ovarian) | 1 | FAIL |

**The zeros are a detection floor, not a measurement.** Cluster-level annotation
cannot see an NK pool smaller than roughly one cluster: 124–614 cells in the
libraries that returned none — every floor **above** the 30-cell threshold being
tested. Those rows mean "not detectable at cluster resolution", not "fewer than
30", so the FAIL verdicts may under-count available power and cannot over-count it.

**Two independent lines converged on the same obstacle.** My measurement found
that lung donors have matched normal tissue but thin tumour NK (281–509) while the
one NK-rich tumour (ovarian SU-O-005, 2,476 NK) has no matched normal ovary. The
Phase B agent, reasoning only from the literature, wrote: "where compartment A
exists, compartment B is thin; where B is rich, A is missing."

**What the Phase C pass does not support.** It does **not** mean X1 can be answered
in humans. Phase B returned `mouse_only = TRUE` for a reason that is logical, not
logistical: **no human sample exists of an organ from before its own tumour arose.**
Contemporaneous uninvolved lung is a proxy that assumes the uninvolved region was
not remodelled by the tumour-bearing state. Passing power means a
compartment-sharing analysis is *runnable* on GSE302113; it does not mean that
analysis answers the temporal clause of X1.

---

## 5. Rule G8 — where the spec's own initial calls were overturned

The adversarial pass changed **9 of 25** classifications, in both directions, and
re-examined three more (P1, T2, X1) without changing the class — T2's weight was
downgraded within `OBS_TRANSFER`, and P1 and X1 were conditional verdicts resolved
by adjudication (DEC-15).

**Every citation the spec flagged as possibly misremembered is real**: Picant
*Nat Commun* 2025 (PMID 40610398), Serger *Sci Immunol* 2026 (PMID 42247486),
Schmid/Wiedemann (bioRxiv 2026.02.11.705354), Barahona/Yokoyama *eLife* 2026
(PMID 42417504), Gamliel *Immunity* 2018. The errors were **wrong evidence classes
attached to real papers**, which is the failure mode G1 targets rather than the
fabrication G8 anticipated.

| claim | spec said | audit found |
|---|---|---|
| T6 (row 21) | TRANSFER | `INF_INVITRO` — no transfer exists in the paper; zero hits for "adoptive transfer", "congenic", "CD45.1", "NSG" in the full text. n = 3 (human blood donors) |
| R3 (row 11) | OBSERVATION | `INF_INVITRO` — TGF-β genuinely withdrawn and ATAC and function genuinely measured separately, but entirely in vitro. n = 7 (human blood donors) |
| D3 (row 4) | mostly inference | `INF_INVITRO` — restored to the spec's call by a route it did not anticipate: the audit graded `OBS_TRANSFER` on Chan 2007's NOD-SCID transfer, and the adversarial pass refused to credit an `OBS_*` class to a design known only from an abstract (403 on all hosts, no PMC). n not stated in any accessible text |
| P1 (row 6) | "best evidence, near clonal" | `INF_KINETIC` — human evidence is longitudinal population phenotyping. n not extractable from the founding source (Gumá 2004, paywalled); unit is donors |
| S1 (row 13) | pseudotime only | worse — see below. n = 24 women (scRNA-seq), 43–47 (flow) |
| T2 (row 17) | "Barahona Ponce" | class upheld (`OBS_TRANSFER`), author corrected to Josselyn D. Barahona; "Barahona Ponce" is a different researcher (gallbladder cancer genetics). n = 4 pregnant dams |

**S1 is the largest overturn.** The spec called it "pseudotime only = pure
inference." It is weaker than that: **Vento-Tormo 2018, the paper that defines
dNK1/2/3, makes no ordering claim at all** — its methods state that "only cells
that were identified as trophoblast were considered for trajectory analysis." The
ordering is a later accretion, and the three papers asserting one give three
**mutually incompatible topologies**: Wang 2021 orders dNK1→dNK2→dNK3 with dNK1
*immature*; Huhn 2020 places dNK1 as the *mature* endpoint — opposite polarity for
the same subset; Guo 2021 recovers three parallel sibling branches with no
ordering among dNK1/2/3 at all.
*Does not support*: that dNK subsets are unordered. It supports only that no
published ordering is anchored in anything but transcriptome adjacency, and that
the published orderings contradict each other.

**P3 (row 8) changed kind, not just class.** The issue is not a false transition
but an **entity misassignment**: Paust 2010 gated CD45⁺NK1.1⁺CD3⁻ with no CD49a
and no DX5, pooling conventional NK with liver ILC1; Wang 2018, sorting all three
populations from the same donors, found only the IL-7Rα⁺ ILC1 fraction conferred
hapten recall. n = 4–6 mice per group.
*Does not support*: that liver memory does not exist. Transferable
antigen-specific liver memory in T/B-deficient mice is replicated across four
labs. What is unsupported is the word "NK" in "liver memory NK cell".

---

## 6. Occupancy (B4) — and a cross-check on the criterion most easily faked

`already_traced = TRUE` for **P1, T1, T3, T5**.

**P1's occupant was found twice, independently.** While building the C.1 dataset
inventory I found Rückert T, Lareau CA, et al., *Nat Immunol* 2022 (PMID 36289449,
GSE197037/GSE197008): human mtscATAC-seq, 4 HCMV⁺ and 3 HCMV⁻ donors, mgatk,
variants filtered at strand concordance > 65% / VMR > 0.01 / ≥ 3 cells, clonotypes
from a neighbourhood graph on mtDNA mutation frequency. The Phase B occupancy
agent then returned `already_traced = TRUE` at HIGH confidence citing the same
paper. B4 is the criterion most likely to be answered from memory, and a wrong
`FALSE` there is what corrupts a shortlist, so the independent agreement matters.

*Does not support*: that P1's descent claim is settled. Rückert establishes clonal
expansion and persistence **within** the adaptive compartment; P1 asserts descent
from conventional NK. The occupancy verdict reflects that the instrument has been
applied to this population, not that the arrow has been drawn.

---

## 7. Limitations, stated rather than repaired

1. **My own coverage threshold would refuse the field's defining dataset.**
   Rückert 2022 operates at **11–20× median chrM coverage**, below the 20× floor
   preregistered in DEC-04. A peer-reviewed application of exactly this method is
   inadmissible under my own rule. The threshold stays as preregistered; loosening
   it after seeing what it excludes is the post-hoc adjustment C.2 forbids.
2. **Cluster-level NK typing has a 124–614-cell detection floor** (§4), above the
   30-cell threshold it feeds.
3. **NK counts are a power estimate, not a biological result.** The CD56 control
   validates NK-rich vs NK-free *blood*; it does not validate NK calls in tumour
   or lung tissue, where the ILC1/tissue-NK boundary is exactly what T1, T4 and X1
   dispute. The passing rows clear threshold by 1–3 orders of magnitude, so the
   verdicts survive large classification error.
4. **`fulltext_available = FALSE` on P3 only** (O'Leary 2006, no PMC deposit,
   Nature paywall). 10 claims are TRUE, 14 PARTIAL. No row's class was asserted
   from an abstract; where the decisive paper was unreadable the row says so
   (D3 is the worked example).
5. **Two shortlist entries rest on MEDIUM-confidence occupancy searches** (§1).
6. **Phase A adversarial schema was under-constrained**: several reviewers
   returned reasoned prose where a bare enum was expected. Normalisation lifted a
   leading token where given, kept the audit class where a reviewer challenged
   without supplying one, and routed the two genuinely conditional verdicts (P1,
   X1) through an explicit adjudication table rather than string-matching
   (`scripts/normalize_verdicts.py`, DEC-05/DEC-15).

---

## 8. What this round did not do

No review, no summary of biological significance, no new directions, and **no
opinion on the importance or scientific value of any claim** — including the two
on the shortlist and including X1. Those were non-goals fixed in the
preregistration before any retrieval, and they remain unaddressed by design.
