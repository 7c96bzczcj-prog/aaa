# NK 趋化因子跨区室调查规范 — the bridge from the closed 四格 record

Status: **dataset-side groundwork only.** The protocol text itself is
not in this repository, so nothing here implements it. What is here is
the work that any version of it needs first, and that the closed
four-quadrant record already licenses: *which chemokine axes this
dataset can measure, and between which compartments.* Everything is
computed from files already on disk by
[`scripts/chemokine_measurability.py`](../scripts/chemokine_measurability.py)
— no new download, no re-gating.

The reason for doing this before reading the spec is D10. There, the
protocol's own positive control (`TOX`) was carried through five phases
and reported as a result before anyone checked that it cleared the
detection floor; it did not, and the headline claim had to be retracted.
Measurability is therefore checked first here, against the panel rather
than against a conclusion.

---

## 1. What carries over from the closed project

| asset | state | use to a chemokine protocol |
|---|---|---|
| balanced pseudobulks, 5 lineages × tumour/normal | on disk, 29 individuals | the substrate; no re-gating needed |
| ambient soup profile + per-lineage ρ | on disk | the only calibrated way to tell NK expression from myeloid pickup in this dataset |
| paired voom+limma, TOST equivalence, permutation nulls | implemented, tested | reusable estimators |
| detection floor / saturation guards (2.5, 2.6) | implemented; 2.6 needs raw cells | the floor gate below is the pseudobulk half only |
| synthetic-recall harness | implemented | measures sensitivity by construction, for any new call rule |

**Dead on arrival, and worth stating so no phase is written that depends
on them:**

- **Anything that reads "NK changes, the others do not" as biology.**
  Three NK-specific technical effects were measured on this dataset —
  ambient burden 2.6× CD8T's, *differential* ambient between the two
  sides of the same library, and a tumour-side NK RNA-content loss
  (log2 −0.323, p = 6.0 × 10⁻⁸). All three are Q4-shaped by
  construction. A chemokine protocol inherits every one of them: they
  are properties of the dataset, not of the four-quadrant rule.
- **Replication in a second cohort.** The GEO rescan found exactly one
  public CD45⁺ + paired + ≥20-individual dataset, and it is the one in
  use.
- **Decontamination as a way out.** Both affordable corrections failed
  their own acceptance tests; CellBender is ~300 h in this container.

---

## 2. The compartment inventory — the blood arm is n = 2

GSE154826 has 108 human libraries over 35 individuals. A cross-compartment
contrast is paired, so its n is the number of *individuals holding both
compartments*:

| contrast | individuals |
|---|---|
| tumour vs adjacent normal | **29** |
| tumour vs blood | **2** |
| normal vs blood | **2** |
| all three | **2** |

Blood exists in this cohort — 11 PBMC libraries — but in two donors.
Whatever "跨区室" is meant to span, **a blood↔tissue arm cannot be run
here**: n = 2 fails the protocol family's own admission rule (A4, 30
cells per lineage per group, ≥ 20 individuals) by an order of magnitude,
and no amount of analysis repairs it.

Two further facts about those PBMC libraries, both of which would bite a
protocol that tried anyway:

- Their `prep` field is **unrecorded**, while the tissue libraries are
  CD45⁺ bead / CD45⁺ FACS / dead-cell-depleted. A blood-vs-tissue
  contrast would therefore confound compartment with enrichment
  protocol — and this project has already measured that enrichment, not
  cohort size, is what decides whether a rare lineage is readable.
- They were never processed. The Phase 2 gate admitted tumour/normal
  only, so no PBMC pseudobulk exists, and `data/raw/` is not in this
  container. Adding them means re-fetching the tarballs
  (`scripts_fetch.sh`) and a re-run of Phase 2 — cheap in code, not in
  time.

So on this dataset "cross-compartment" can only mean **tumour vs
adjacent normal lung**, at n = 29.

---

## 3. Measurability — the output is readable, the receptor repertoire is not

Full chemokine system (18 receptors + 5 atypical + 37 ligands) plus the
retention/egress module, against the protocol's own floor (≥ 10 counts in
≥ 50 % of NK pseudobulk samples; at the mean NK library of 637,948
counts that is ≈ 16 CPM) and against the ambient scales calibrated inside
NK from genes of known truth (genuinely-NK genes read soup fraction
0.013; known-ambient genes read 0.402, their true value being 1.0;
NK/Myeloid for myeloid-only genes 0.070–0.118).

| group | measurable in NK |
|---|---|
| receptors, inflammatory | **1 / 9** — only `CX3CR1` |
| receptors, homing | **1 / 9** — only `CXCR4` |
| receptors, atypical | **0 / 5** |
| retention / egress | **6 / 8** — `S1PR5`, `SELL`, `SELPLG`, `CD69`, `ITGAE`, `KLRG1` |
| ligands, lymphoid | **7 / 7** — `CCL3`, `CCL3L3`, `CCL4`, `CCL4L2`, `CCL5`, `XCL1`, `XCL2` |
| ligands, myeloid/stromal | **3 / 30** |
| **panel total** | **18 / 68** |

**The asymmetry is the finding.** NK's chemokine *output* is fully
measurable here — all seven lymphoid ligands clear both gates, at 128 to
5,295 CPM, with soup fractions of 0.02–0.12 against a 0.402 ambient
reading. NK's chemokine *input* is not: 2 of 23 receptors survive. A
protocol phrased as "which chemokine receptors position NK in the
tumour" is unanswerable on this dataset; one phrased as "what does NK
secrete, and how does that differ between compartments" is answerable.

Three qualifications, because the table is easy to over-read:

- **Below the floor is not absence.** 18 of the 46 below-floor genes sit
  within 2× of the floor (8–32 CPM), including `CXCR3`, `CXCR6`,
  `CCR7`, `S1PR1` and `ITGA1` — all real NK receptors. `CXCR3` misses on
  detection breadth (0.483 against the 0.50 gate), which is a knife-edge.
  This is 3′ V2 sampling depth, not biology, and it names its own fix: 5′
  chemistry, deeper libraries, a targeted panel, or protein (the 20 ADT
  libraries in this cohort are unexploited).
- **The floor gate here is the pseudobulk half only.** The cell-level
  saturation guard (2.6) — the one that caught `TOX` — needs per-cell
  counts, i.e. `data/raw/`. "MEASURABLE" in the table means *not
  excluded by the pseudobulk gates*, and must be re-checked against 2.6
  before any gene is read.
- **`AMBIENT-INDISTINGUISHABLE` and `BELOW FLOOR` are different deaths.**
  `CXCL8`, `CXCL16`, `CCL20`, `CXCL3` are abundant in NK pseudobulk and
  fail on the soup scale — they are myeloid transcripts arriving as
  pickup. That is a contamination result, not a detection one, and it is
  the reason the panel cannot be read on CPM alone.

`S1PR5` — the canonical NK egress receptor, and the single most relevant
gene to a cross-compartment question — is measurable (148 CPM, soup
0.009, NK/Myeloid 264) but **falls out of the tested universe**, because
the relaxed filter demands a witness lineage and `S1PR5` is NK-restricted
(D7). A chemokine protocol must set its own gene universe rather than
inherit this one, and state the witness rule it uses.

---

## 4. Effect sizes: computed, and not readable

`results/chemokine_panel_verdict.csv` carries the paired tumour-vs-normal
statistics for the 18 survivors, joined from the existing Phase 5 table.
They are printed by the script as a feasibility readout and are recorded
here as **descriptive only**: the largest NK effects in the panel are
`CCL2` (+0.82), `CD69` (+0.69), `CXCL2` (+0.60), `XCL1` (+0.57) and
`SELPLG` (−0.52), and every one of them is a candidate for the same
NK-specific artefacts that falsified the four-quadrant reading. Two
genes (`SELL`, `XCL2`) carry a Q4 call from the retracted analysis; that
call is not evidence.

Nothing in section 4 should enter a hypothesis without the Phase 8
B1–B3 controls that were never run.

---

## 5. What a workable protocol has to pre-register

Written as constraints the dataset imposes, so the spec can be checked
against them rather than discovering them in phase 5:

1. **The measurability gate runs before the hypothesis**, and the panel
   is fixed in writing before the table is read. Section 3 is that gate,
   already run.
2. **The compartment axis is tumour vs adjacent normal, n = 29.** Any
   blood arm is a study-design request, not an analysis.
3. **Every NK-vs-other-lineage claim carries the three NK-specific
   confounders**, and needs a control that is not "the other lineages
   did not move" — the closed record shows why that control is void
   here.
4. **The ambient scale is calibrated inside NK from genes of known
   truth**, per gene, per lineage. Raw CPM is not admissible evidence
   for a myeloid-dominant chemokine.
5. **The cell-level saturation guard must be armed** (needs `data/raw/`),
   or every "MEASURABLE" call stays provisional.
6. **Sensitivity is measured by construction, not assumed** — the
   synthetic-recall harness applies to whatever call rule the protocol
   defines, and reporting a rule without its recall is what D14 exists
   to prevent.

---

## 6. Open, and needed from the spec

- The protocol text ("NK 趋化因子跨区室调查规范") is not in the repo.
  The phases, stop rules, and the intended meaning of 区室 —
  tumour/normal only, or blood-inclusive, or spatial — decide which of
  the above is a constraint and which is a blocker.
- If 区室 is blood-inclusive, this dataset is the wrong instrument
  (n = 2) and the choice is a different cohort or a redefinition.
- If the receptor repertoire is the object of study, 3′ pseudobulk is
  the wrong readout and the 20 unexploited ADT libraries are the first
  place to look.
