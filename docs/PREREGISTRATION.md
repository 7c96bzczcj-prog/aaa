# PREREGISTRATION — NK lineage-assertion evidence screen

**Spec**: NK 谱系断言筛查规范 v1.0
**Committed**: 2026-08-14, **before any literature retrieval was performed.**
**Status after first commit: READ-ONLY.** Any change requires a v1.1 bump
recorded in `docs/CHANGELOG.md` with reason and timestamp. Post-hoc
relaxation of any threshold below invalidates the affected rows.

This file fixes, in advance: (§5) the shortlist admission criteria, (§2)
the evidence-type taxonomy and the abstract-vs-methods rule, (§3) the
decidability gates, (§4/C.2) the power thresholds and the chrM
variant-calling threshold, and (§6) the hard rules whose violation voids
a row.

---

## 1. Scope fixed in advance

The claim list is the 25 rows of spec §1 (D1–D4, P1–P3, T1–T6, S1–S3,
X1–X5, R1–R4). Additions or deletions are permitted but must be recorded
in `docs/DECISIONS.md` with a reason, and any added claim is audited
under exactly these criteria.

The "Claude 初判" column of spec §1 is entered as a **hypothesis to be
overturned**, not as evidence. Every row is re-derived from primary
sources; the initial call is compared against the audit result only
afterwards, and disagreement is recorded in `audit_overturns_claude`.

**Non-goals, fixed here so they cannot be drifted into**: no review, no
summary of biological significance, no new directions, and **no opinion
on the importance or scientific value of any claim.** This round produces
an evidence grading and a decidability verdict, nothing else.

---

## 2. Evidence-type taxonomy (fixed)

Each claim is assigned to **exactly one** strongest class; all classes
present are also listed.

| code | definition | strength |
|---|---|---|
| `OBS_DIRECT` | the same cell or clone is observed to convert (live imaging, barcoding, genetic lineage tracing) | strongest |
| `OBS_TRANSFER` | adoptive transfer: A is transferred, B is recovered from the recipient | strong |
| `OBS_GENETIC` | genetic necessity: gene knocked out, B fails to appear | medium — proves necessity, not lineage |
| `INF_INVITRO` | A cultured into a B-like phenotype in vitro | weak |
| `INF_MARKER` | A and B share/overlap markers | weak |
| `INF_TRAJECTORY` | pseudotime / RNA velocity / trajectory inference | **weakest** |
| `INF_KINETIC` | A appears before B in time | weak |

**Fixed rule (G2): `INF_TRAJECTORY` is never an observation.** Pseudotime
orders cells by transcriptome similarity and then calls that order time.
It cannot distinguish "A becomes B" from "A and B are two independent
states that happen to be adjacent in expression space." No trajectory
result may be promoted to `OBS_*` for any reason.

**Fixed rule (G1): methods, not abstracts.** `what_was_measured` must be
written in the form "under condition X, at timepoint Y, marker Z was
measured", never "showed that A becomes B". If full text is unobtainable,
`fulltext_available = FALSE` with the reason (paywall / 403 / no PMC),
and the row may not carry a methods-derived classification asserted from
an abstract.

**Fixed rule (G3): not-found ≠ absent.** `counter_evidence` and
`already_traced` take the value `NOT_SEARCHED` whenever a search was not
run or returned nothing that was actually read. `NONE` may be recorded
only for a search that was run and read.

**Fixed rule (G5)**: `n` is always accompanied by `n_unit` ∈ {donors,
animals, wells, cells, embryos, samples}. An `n` without a unit is void.

---

## 3. Phase B gates (fixed)

Applied to every claim whose Phase A strongest evidence type is `INF_*`.
All four must be answered; **all four must pass** to reach Phase C.

- **B1 `lineage_decidable`** — clonal barcoding supplies exactly one
  primitive: whether two cells share an ancestor, and how recently. A
  claim of the form "the progeny of A are B" is decidable. A claim of the
  form "A becomes functionally B-like" is **not** decidable — lineage
  does not measure function.
- **B2 `both_answers_actionable`** — the next experiment must be written
  out separately for the "yes, lineage-related" and "no, not
  lineage-related" answers. **If the two are the same, the claim is
  out.** This is the primary filter.
- **B3 `compartments_obtainable`** — the compartments containing A and B
  must be obtainable *simultaneously*, preferably from the same
  individual. The required compartments are enumerated and their
  pairability stated.
- **B4 `already_traced`** — search must cover at least: `mtDNA` +
  lineage / clonal + the claim's cell types; `mtscATAC`; `MAESTER`;
  `barcode` + NK; and for mouse, `fate mapping` / `confetti` /
  `lineage tracing`. Records `TRUE` / `FALSE` / `NOT_SEARCHED`.

---

## 4. Phase C thresholds (fixed — no post-hoc adjustment)

Phase C is run only for claims passing Phase B, **plus X1 unconditionally**
(user-designated highest priority, run regardless of its score).

### 4.1 Dataset admissibility

Modality preference, fixed: `mtscATAC-seq` > `MAESTER` > plain scATAC
with retained chrM reads > scRNA (mtDNA variant detection generally
insufficient). **Whether the archived files retain chrM fragments must be
recorded for every dataset**; a chrM-stripped archive is unusable and is
recorded as such rather than silently dropped.

### 4.2 Power thresholds

| threshold | value | unit |
|---|---|---|
| NK cells per donor per compartment | **< 30 → that cell is unusable** | cells |
| donors with ≥ 2 comparable compartments each | **< 3 → dataset unusable** | donors |
| median per-cell chrM coverage | **< 20× → dataset unusable** | mean reads per chrM base per cell |

The 20× figure is the pre-registered operating point for `mgatk`-style
heteroplasmy calling (the standard mtscATAC-seq pipeline; Lareau et al.
2021 report ~40–100× and mgatk's own guidance treats low-coverage cells
as uncallable). An **informative variant** is pre-defined as one passing
mgatk defaults: strand concordance ≥ 0.65, variance-mean ratio ≥ 0.01,
detected in ≥ 5 cells. If a different caller is used, its documented
threshold is substituted and the substitution is recorded in
`docs/DECISIONS.md` **before** the numbers are computed.

**If the power check fails, the result is "fails power". Thresholds are
not lowered to manufacture a runnable analysis.**

### 4.3 Species constraint

mtDNA barcoding works in humans and is **weak in young inbred mice**:
somatic mtDNA variants accumulate with age, and laboratory mice are both
young and isogenic, so informative variants are scarce. Any
mouse-dominant lineage plan must argue its feasibility explicitly in
`docs/DECISIONS.md` or switch to a non-mtDNA recorder (Confetti /
Polylox / CRISPR recorder / genetic fate mapping). Claims answerable only
in mouse are flagged `mouse_only = TRUE`.

---

## 5. SHORTLIST ADMISSION CRITERIA (the pre-registered decision rule)

Verbatim from spec §5:

> **进入最终短名单的条件，四条全满足：**
> 1. `load_bearing = HIGH`
> 2. Phase A 最强证据类型属于 `INF_*`
> 3. `lineage_decidable = TRUE` 且 `both_answers_actionable = TRUE`
> 4. `already_traced = FALSE`（若为 `NOT_SEARCHED`，该条留在短名单但标 `occupancy_unverified = TRUE`）

Operationally:

A claim enters the final shortlist **iff all four hold**:

1. `load_bearing == HIGH`
2. `evidence_type_strongest` ∈ {`INF_INVITRO`, `INF_MARKER`,
   `INF_TRAJECTORY`, `INF_KINETIC`}
3. `lineage_decidable == TRUE` **and** `both_answers_actionable == TRUE`
4. `already_traced == FALSE`; if `already_traced == NOT_SEARCHED`, the
   claim **remains** on the shortlist but is flagged
   `occupancy_unverified = TRUE`

Note the asymmetry, fixed here in advance: `OBS_GENETIC` is **not** an
`INF_*` code, so a claim whose strongest evidence is genetic necessity
fails criterion 2 and does not enter the shortlist, even though genetic
necessity does not establish lineage. This is a deliberate consequence of
the taxonomy as written and is not to be reinterpreted after seeing which
claims it excludes.

### 5.1 The three possible outcomes and their meanings (written in advance)

- **Non-empty shortlist** — a set of questions whose evidence is
  inferential, which lineage can decide, where both answers are
  actionable, and which have not been done. **This does not mean they are
  worth doing.** Moat, feasibility and time cost are a separate layer of
  judgement, explicitly outside this round.
- **Empty shortlist** — these claims are either already hard enough,
  or not lineage-decidable, or already done. **This is a valuable
  conclusion, not a failure**, and is reported as such.
- **Many rows stuck at `NOT_SEARCHED`** — search depth was insufficient
  to support a verdict. Reported as **"incomplete"**. Unsearched is never
  reported as blank or as negative.

### 5.2 Early-termination condition (fixed)

Verbatim from spec §5:

> **提前终止条件**：若 Phase A 的前 6 条里有 3 条以上拿不到全文，停下来先报告文献获取问题，不要继续硬跑——一份建立在摘要上的证据分级表是有害的。

Operationally: the first six claims audited are **D1, D2, D3, D4, P1,
P2**, in that batch. If **more than 3** of those six come back
`fulltext_available = FALSE`, Phase A halts and a literature-access
report is delivered instead of a completed grading table. A grading table
built on abstracts is harmful and is not to be produced.

---

## 6. Hard rules that void a row (fixed)

Violating any of these voids the affected result.

- **G1** Read the methods, not the abstract. "We show A becomes B" in an
  abstract routinely corresponds to marker similarity in the methods.
- **G2** `INF_TRAJECTORY` is never an observation.
- **G3** Not found ≠ does not exist. Record `NOT_SEARCHED`.
- **G4** **Search the mechanism's home field, not the NK field.** Any
  claim invoking a general process (differentiation, plasticity,
  residency, exhaustion) must additionally be searched once in the T
  cell / ILC / myeloid literature. A borrowed mechanism's fatal collision
  lives in the field it was borrowed from and is invisible to NK
  keywords.
- **G5** `n` is donors/animals, not cells. Record the unit of every `n`.
- **G6** Any summary quantity is computed on the same basis as the rule
  that consumes it.
- **G7** No loosening of a classification to make a claim interesting.
  Classification follows only what the methods did.
- **G8** Where the §1 initial call is wrong, say so explicitly and change
  it. That table is a hypothesis to be overturned.

## 7. Explicitly forbidden (fixed)

- treating an abstract as methods
- treating pseudotime / RNA velocity as lineage evidence
- writing "not searched" as "blank" or as "nothing there"
- modifying §5 after the fact
- opining on the scientific value of any claim
- loosening any criterion in order to make the shortlist non-empty
