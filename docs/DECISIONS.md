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
