# 四格挖掘协议 — implementation and Phase 0–2 record

Mining "NK resistance" genes (Q4) and the shared cytotoxic-lymphocyte
programme (Q2) from paired tumour / adjacent-normal single-cell data.

This repository contains a validated implementation of Phases 2–8, plus
the completed record for Phase 0 (prior art), Phase 1 (dataset
admission) and Phase 2 (gating and cell counts) on GSE154826.

**No Q4 gene list exists yet, and none should be quoted from this
repository.** What is finished is the part that decides whether a Q4
list would mean anything: the prior-art check, the dataset admission
gate, and a statistical core tested against planted ground truth.

---

## Headline results

### Phase 0 — stop rule S1 does not fire

No prior work defines or systematically reports the Q4 class, and no
prior work formalises cross-lineage effect sharing using equivalence
testing. Details and the query battery: [`docs/PHASE0_LITERATURE.md`](docs/PHASE0_LITERATURE.md).

One finding sharpens the project. The pan-cancer NK atlas (*Cell* 2023,
716 patients) reports tumour-infiltrating NK cells as marked by
**DNAJB1, HSPA1A, FOS, JUN** — which is the Marsh dissociation-stress
set, a textbook Q3 candidate, being read as NK biology. Nobody runs the
cross-lineage control that would separate the two readings. Even if Q4
comes back empty, showing that a widely-cited NK signature is Q3 rather
than Q1 is a result.

### Phase 1 — GSE154826 verified in full; stop rule S2 does not fire

Everything the protocol remembered about GSE154826 was checked directly
against GEO and the raw files, and all of it holds:

| Criterion | Verified value |
|---|---|
| A1 paired individuals | **29** patients with both tumour and adjacent normal |
| A2 condition ⟂ library | **Pass** — 20 libraries carry *both* conditions in one droplet emulsion |
| A3 all lineages in one pool | **Pass** — CD45⁺ enrichment; 4 CD2⁺-enriched libraries excluded as instructed |
| A4 NK ≥ 30 per group | **Pass — 27/29 patients** (measured, see below) |
| Raw droplet matrices | **Yes** — `barcodes.tsv` = 737,280 lines = full 10x v2 whitelist |
| HTO shared soup | **Yes** — 20 libraries, 8 patients, tumour + normal in one soup |
| ADT | **Partial — 20 of 73 libraries only** (see correction below) |

**Correction to the protocol's assumptions.** The ADT panel is present
in only the 20 hashed libraries, not dataset-wide, and even there
CD56/CD8/CD19/CD16 median counts are 2–3 per cell — too sparse to
threshold. Protein-based gating, which the protocol prefers precisely to
avoid circularity, is therefore available for 8 of 29 patients and
cannot carry the main analysis. Lineage assignment is RNA-based, with
every gating gene barred from testing via `excluded_genes.txt`.

GSE178341 (Pelka CRC) is admitted as the Phase 7 replication cohort:
**36 paired patients**, and the protocol's memory that it has no raw
droplet matrices is confirmed. GSE131907 passes A1 with 10 paired
patients. GSE176078 (Wu breast) is **rejected** — tumour-only, no paired
normal, fails A1 outright. Full table: [`phase1_registry/registry.csv`](phase1_registry/registry.csv).

### Phase 2 — A4 measured on real data

73 libraries processed, 463,518 cells called from 737,280-barcode raw
droplet matrices.

```
NK >= 20 in both conditions : 28/29 patients
NK >= 30 in both conditions : 27/29 patients   <- protocol threshold
NK >= 50 in both conditions : 27/29 patients
NK >= 100 in both conditions: 23/29 patients
```

NK is the limiting lineage in **53 of 64** (patient × condition) groups,
which is exactly why Phase 2.4 power balancing is not optional. Median
limiting-lineage size is 242 cells, and 61 of 64 groups clear 30.

### The protocol's central premise, quantified

Planted-truth simulation where **every gene is a true Q3** (all five
lineages move, NK included), with NK starved of cells. Every Q4 call is
false by construction:

| NK cells/sample | false Q4 by "p > 0.05 ⇒ unchanged" | false Q4 by TOST |
|---|---|---|
| 3 | **49 / 58** | **0** |
| 5 | 36 / 60 | 0 |
| 8 | 22 / 60 | 1 |
| 15 | 8 / 60 | 3 (= 5%, TOST's nominal rate) |
| 150 | 0 | 0 |

Without equivalence testing, a rare lineage produces an ~84% false Q4
rate. With it, the error rate is at its nominal level. The protocol's
core design decision is correct and now has a number attached.

---

## Two specification bugs found during validation

Both were found by running planted ground truth through the classifier,
and both are documented with evidence in [`docs/DEVIATIONS.md`](docs/DEVIATIONS.md).

**1. Q3 absorbed every Q4 gene.** The Q3 rule as written
(`>= 4 lineages changed`) is satisfied by any Q4 gene, since a Q4 gene
has four lineages moving. On planted truth this misfiled **25 of 25**
true Q4 genes as Q3 — the protocol's primary endpoint had a
true-positive rate of zero. Q3 now additionally requires that NK is not
TOST-equivalent. Q4 recovery: 0/25 → **24/25**.

**2. The `p50` equivalence criterion discards most of the signal.** A
gene with a genuinely zero effect sits below its own null *median* only
half the time. Requiring that across four lineages retains ~6% of true
Q1 genes:

| equivalence bound | Q1 | Q2 | Q3 | Q4 | false Q4 |
|---|---|---|---|---|---|
| `p50` (as written) | 0/25 | 3/25 | 25/25 | 18/25 | 0 |
| `p95` (recommended) | 18/25 | 20/25 | 25/25 | 24/25 | 0 |

Loosening it costs **nothing** in false positives, because TOST already
controls false equivalence at 5%. The bound is now a parameter; the
protocol default still runs.

A third issue is documented but not a bug: the Phase 3 split-half null
is *unpaired* while Phase 4 is *paired*, so its percentiles are wider
than the estimator's true noise — and that bias runs toward
manufacturing Q4 calls. A structurally matched sign-flip permutation
null is implemented and used by default.

---

## What is NOT done

Being explicit, because a half-run protocol that looks finished is worse
than one that looks unfinished:

- **Phases 3–5 have not been run on real data.** They are implemented
  and pass planted-truth validation; they have not been executed on
  GSE154826 pseudobulks.
- **Phase 6 (ambient) not run.** CellBender has not been run. The
  shared-soup HTO contrast — the protocol's cleanest control, and the
  strongest single asset of this dataset — is implemented but not
  executed.
- **Phase 7 (replication) not run.** GSE178341 is admitted but not
  downloaded.
- **Phase 8 not run.**
- **NK gate purity is unresolved.** The independent T-probe panel
  (LCK/CD2/THEMIS/SKAP1/TRAT1/ITK) shows 57.6% positivity in the NK gate
  against 25% in the B gate (the soup floor). This is *not*
  interpretable as contamination, because LCK, CD2 and SKAP1 are
  genuinely expressed by NK cells — the probe panel was poorly chosen.
  The sub-cluster purity check that protocol 2.1 says actually catches
  hidden T cells has **not** been run. Treat NK purity as an open
  question.
- **The limma reimplementation has not been diffed against R limma.**
  It is validated against planted truth and closed-form references, not
  against the reference implementation.

---

## Layout

```
src/nkmine/
  de.py                voom + limma empirical-Bayes moderation (Phase 4)
  quadrant.py          TOST equivalence + four-quadrant rules (Phase 5)
  pseudobulk.py        power balancing, detection/saturation gates (Phase 2.3-2.6)
  null_calibration.py  per-lineage per-gene permutation nulls (Phase 3)
  gating.py            cell calling, HTO demux, lineage assignment (Phase 2.1-2.2)
  artifacts.py         Kitagawa 3-term decomposition, ambient (Phase 6)
  controls.py          built-in positive/negative controls, stop rule S4 (Phase 5.3)
  replication.py       cross-dataset replication, mundane screen (Phases 7-8)
  simulate.py          planted-truth generator
scripts/phase2_gate_all.py   Phase 2 driver over all 73 libraries
phase1_registry/registry.csv Phase 1 admission table
results/                     A4 counts, purity, library metadata
docs/DEVIATIONS.md           every departure from the protocol, with evidence
docs/PHASE0_LITERATURE.md    prior-art record
tests/                       33 tests, all passing
```

## Running it

```bash
python3 -m pip install numpy scipy pandas statsmodels pytest
python3 -m pytest tests/ -q

# Phase 2 (needs the GEO tarballs in data/raw/)
python3 scripts/phase2_gate_all.py
```

Stop rules are implemented as functions that return an explicit verdict
rather than as comments: `pseudobulk.se_balance_report` (S3),
`controls.evaluate_stop_rule_s4` (S4),
`replication.evaluate_stop_rules_s5_s6` (S5/S6).

---

## Suggested next step

The single highest-value remaining action is the **shared-soup HTO
contrast** (Phase 6.1) on the 20 hashed libraries covering 8 patients.
It needs no decontamination tool, no parameter choices, and no external
data: within those libraries both conditions share one ambient soup, so
ambient contributes the same offset to both and largely cancels. Running
Phases 4–5 separately on the hashed and unhashed sets and comparing
cross-lineage concordance gives a direct estimate of ambient's
contribution — and it is a cheap, decisive check that can be done before
committing to the full pipeline.
