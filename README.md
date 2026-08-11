# 四格挖掘协议 — implementation and Phase 0–8 record

Mining "NK resistance" genes (Q4) and the shared cytotoxic-lymphocyte
programme (Q2) from paired tumour / adjacent-normal single-cell data.

This repository contains an implementation of Phases 2–8, the completed
record for Phase 0 (prior art) and Phase 1 (dataset admission), and
executed results for Phases 2–7 plus part of Phase 8 on GSE154826, with
replication attempted against GSE131907.

The work was then audited adversarially, and the audit reversed a
headline claim. **Stop rule S4 fires: the
pipeline has no valid Q4 positive control, so its 46 Q4 candidates are
not cleared for reading.**

An earlier revision claimed this was demonstrated by `TOX` reaching Q4.
That rested on the Phase 2.6 guard, which was implemented and documented
but **never actually called**; armed, `TOX` is floor-saturated (4.3% /
4.9% NK detection) and excluded. That claim is retracted
([D10](docs/DEVIATIONS.md)) and the present one rests on the synthetic
truth set instead.

### Stop-rule status

| rule | condition | status |
|---|---|---|
| S1 | someone already reported Q4-type results | **does not fire** — the class is undefined in the literature, and equivalence testing has never been applied to cross-lineage sharing |
| S2 | fewer than 2 datasets clear A1–A4 | **FIRES** at the protocol's 30-cell threshold — only GSE154826 qualifies. Two qualify at a relaxed 20 cells, the second at exactly n = 8 |
| S3 | NK SE > 2× the other lineages after balancing | does not fire — 0.113 vs 0.079, ratio 1.43 |
| S4 | TOX misses Q4, or IEGs land in Q1/Q4 | **fires as written** — TOX is floor-saturated and can never qualify (D10). But the question S4 exists to ask is now answered by a valid instrument: **76% recall on a constructed truth set** (D14) |
| S5 | fewer than 10 Q4 genes replicate | **fires, and is uninformative** — 0/46 replicate, but so do Q1 (0/10) and Q3 (0/43), and the replication cohort yields only **6** genes from which a synthetic Q4 case could even be built. That cohort cannot test Q4 (D14) |
| S6 | Q4 replication rate below the Q1/Q3 baseline | does not fire — the baseline is itself zero |
| S7 | all survivors explained by B1–B6 | not reachable — Phase 8 incomplete |
| S8 | Phase 0–5 exceeds one week on dataset 1 | not applicable |

**Three stop rules fire, but they no longer mean the same thing.** S2 is
a dataset-availability fact. S4 fires only in its literal wording — its
substance is answered. S5 fires on a cohort now shown incapable of
testing Q4 at all.

The primary endpoint therefore reads: **Q4 is tested and detectable in
the discovery cohort (76% recall), and unreplicated rather than refuted**
— because no available second cohort can test it. That is a materially
different record from "Q4 is rare or absent", which is what the
pre-registered rule would have concluded had recall been measured in
only one of the two designs.

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

### Phase 1 — GSE154826 verified in full

Everything the protocol remembered about GSE154826 was checked directly
against GEO and the raw files, and all of it holds:

| Criterion | Verified value |
|---|---|
| A1 paired individuals | **29** patients with both tumour and adjacent normal |
| A2 condition ⟂ library | **Pass** — 20 libraries carry *both* conditions in one droplet emulsion |
| A3 all lineages in one pool | **Pass** — CD45⁺ enrichment; 4 CD2⁺-enriched libraries excluded as instructed |
| A4 NK ≥ 30 per group | **Pass — 27/29 patients** (measured, see below) |
| Stop rule S2 | **Fires** — only this dataset clears A1–A4 as written |
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

### Stop rule S2 fires at the protocol's stated threshold

A4 was then measured on the candidate replication cohorts too, and it
changes the outlook. "Balanced" below means all five lineages clear the
threshold in both conditions, which is what Phase 2.4 actually needs:

| dataset | paired patients | NK ≥ 30 | balanced ≥ 30 | balanced ≥ 20 |
|---|---|---|---|---|
| GSE154826 (Leader, **CD45⁺**) | 29 | 27 | **27** | **28** |
| GSE131907 (Kim, unsorted) | 10 | 8 | 6 | **8** |
| GSE178341 (Pelka, unsorted) | 36 | 5 | 5 | 7 |

**Only one dataset clears A1–A4 as written, so S2 fires.** At a relaxed
20-cell threshold two qualify, but GSE131907 sits at exactly n = 8.

The mechanism generalises: Pelka has the *most* paired patients (36) and
is still the worst, because unsorted tissue leaves NK at 1.1% of cells.
**CD45⁺ enrichment, not cohort size, is what delivers A4 for a rare
lineage** — so a dataset search should filter on enrichment protocol
first, which is the opposite of the intuitive order and the reason the
protocol's own nomination of GSE178341 does not survive contact with the
data. GSE176078 (Wu breast) is rejected outright — tumour only.
Full table: [`phase1_registry/registry.csv`](phase1_registry/registry.csv).

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

### Phase 6.1 — the shared-soup contrast: ~40% of cross-lineage sharing is technical

GSE154826 contains both designs at once, which makes this measurable
with no decontamination tool. 8 patients have tumour and adjacent-normal
cells pooled into **one droplet emulsion** and separated afterwards by
hashing, so both conditions sit in the same ambient soup and it cancels
from the paired contrast. 19 patients have their conditions built as
**separate libraries**, where it does not cancel — and ambient pushes
every lineage the same way at once.

Power-matched (the separate group repeatedly subsampled to the shared
group's size, since correlation rises with precision and an unmatched
comparison would find the predicted gap whether or not ambient exists).
Patient 581 is excluded from both groups — it contributes libraries of
*both* designs, so it cannot cleanly represent either:

| | mean pairwise cross-lineage r | ≥4 lineages concordant |
|---|---|---|
| shared soup (n=7) | **0.291** | 0.38% |
| separate libraries (matched n) | **0.490** [0.393, 0.583] | 2.29% |

0 of 25 matched draws fell at or below the shared-soup value. So roughly
**40% of the apparent cross-lineage sharing in a conventional design is
technical, not biological**, and the ≥4-lineage concordance rate — the
Q3 rate — is inflated about six-fold.

A depth confounder was checked and runs the *wrong way* for an artefact
explanation: shared-soup pseudobulks are about twice as deep (median
1.13M vs 0.58M counts), and more depth means less attenuation and so
*higher* correlation. The shared group is deeper and still less
correlated, which strengthens rather than weakens the reading.

Two caveats remain, both real. The shared emulsion removes *all*
technical differences between conditions, not ambient alone, so this is
an upper bound on ambient specifically. And the two groups are
**different patients**, so biology is confounded with design; the single
patient carrying both designs is not enough for a within-patient test.

### Phases 3–5 on real data, and why S4 fires

Balanced pseudobulks: **4,525 genes × 270 samples from 27 patients**
(27 × 2 conditions × 5 lineages). **Stop rule S3 passes** — after Phase
2.4 balancing, NK's median moderated SE is 0.113 against 0.079 for the
other lineages (ratio 1.43, limit 2.0). "NK did not change" is
falsifiable on this dataset.

Corrected quadrant counts at δ = 0.5, `p95`: **Q4 46, Q4_attenuated 7,
Q3 43, Q3_NK_indeterminate 10, Q1 10, unclassified 4,409.**

**S4 fires, and the reason is conceptual rather than mechanical.** The
protocol's Q4 positive control is self-defeating. TCR-proximal genes are
T-restricted, so the Phase 2.5 detection floor removes them ([D7](docs/DEVIATIONS.md));
and they are barely expressed in NK, so the Phase 2.6 saturation guard
excludes them ([D10](docs/DEVIATIONS.md)). Both guards are right. The
conflict is that "NK does not respond to a TCR-driven programme" is not
evidence of resistance — NK cells have no TCR, so it is mundane
explanation **B3 (receptor not expressed)**, exactly what Phase 8 exists
to rule out. **The protocol nominates as its positive control an
instance of the artefact it is designed to reject.**

A usable Q4 control has to be a gene NK demonstrably expresses and could
in principle regulate, whose driver is shared with the witness lineages.
Until one exists and passes, no Q4 list from this pipeline is biology.

### Phase 7 — categorical replication fails, the pattern replicates

Replication was attempted against GSE131907 (n = 8, the only cohort
available, at the relaxed threshold). Using the author's own published
cell-type labels, so the replication is independent of this repo's
gating.

**Every quadrant fails categorical replication:** Q1 0/10, Q3 0/43,
Q4 0/46. Read naively that kills the project. It should not be read
naively — the replication cohort's standard errors are 2–4× larger, and
for CD8 T the 90% CI half-width is 0.507 against an equivalence margin
of 0.5, so **a CD8T equivalence call is mathematically unattainable
there**. Q3, the easiest category, replicating at 0/53 is the tell: this
is global power failure, not a verdict on Q4. This is exactly why the
protocol requires the Q1/Q3 baselines alongside Q4 (S6) — without them
the obvious and wrong conclusion is "Q4 is noise".

A power-appropriate continuous test tells a different story. In the
independent cohort:

| group | n | median \|NK log2FC\| | median \|witness log2FC\| | NK/witness |
|---|---|---|---|---|
| **Q4 candidates** | 46 | **0.146** | 0.498 | **0.29** |
| Q3 genes | 43 | 0.573 | 0.900 | 0.64 |
| background | 4,409 | 0.193 | 0.253 | 0.76 |

Q4 candidates carry the pattern into an independent cohort: NK moves
~29% as much as the witness lineages, against 64% for Q3 and 76% for
background. `|NK log2FC|` for Q4 vs Q3: **p = 7.2 × 10⁻⁶** (ratio-free
primary test; Q3 is the right comparator because its witness effects are
*larger*, so this is not a witness-magnitude artefact).

One honest qualification: Q4's absolute NK effect is **not** smaller than
background (0.146 vs 0.193, p = 0.14). What replicates is the *contrast*
— large witness effects alongside a small NK effect — not a uniquely
quiet NK.

The obvious mundane explanation is excluded: Q4 candidates are not
merely underpowered in NK — their NK standard error in the replication
cohort is *smaller* than Q3's (0.249 vs 0.294), and their effects are
not small everywhere (myeloid 0.665, B 0.495, NK 0.145).

The protocol bans ratios in the decision path and that ban is respected:
the ratio above is descriptive, and the primary test is the ratio-free
comparison of `|NK log2FC|` between Q4 and Q3.

### Q4 sensitivity — measured by construction (D14)

The protocol's positive control was unavailable in principle, so the
truth set was built: take real genes whose witness lineages move *and*
whose NK moves, then divide out NK's fitted effect symmetrically,
leaving its variance, cell counts and depth untouched. Those genes are
Q4 by construction.

| λ (NK effect retained) | ≈ residual NK log2FC | recall (Q4) |
|---|---|---|
| **0.00** | 0 | **76.0%** |
| 0.25 | ≈ 0.19 | 49.3% |
| 0.50 | ≈ 0.37 | 0.0% |
| 0.75 | ≈ 0.56 | 0.0% |

n = 75 constructed genes. **76% clears the pre-registered 50% bar: the
design has Q4 power.** The dose curve gives the detection floor as an
effect size — Q4 stays callable while NK's residual effect is under
about **0.19 log2FC**, which matches the arithmetic (δ = 0.5, NK median
SE 0.113 → 90% CI half-width 0.19).

**The construction matters.** Removing the effect by *permuting* the
condition label within individual — the obvious approach — gives **0%
recall**, because sign-flipping a real effect *d* absorbs *d*² into the
residual variance and inflates the SE. Had that construction been used,
this project would have concluded "no Q4 power", the opposite of the
truth, from an artefact of the test.

**And the same test in the replication cohort cannot be run at all**:
GSE131907 yields only **6** genes with moving witnesses alongside a
moving NK. So its zero replication is uninformative — **Q4 is
unreplicated, not refuted.**

### A negative control, on the false-positive side

S4 fires because the Q4 positive control cannot work, which leaves
*sensitivity* unproven. *Specificity* is still testable: flip the
tumour/normal label within each individual — preserving pairing, depth,
composition and ambient structure — and re-run everything.

| quadrant | observed | null (12 permutations) |
|---|---|---|
| Q4 | 40 | **0.0** [0, 0] |
| Q3 | 50 | **0.0** [0, 0] |
| Q1 | 8 | 0.2 [0, 2] |

**0 of 12 permutations produced a single Q4 call.** The pipeline does
not manufacture Q4 from noise. That bounds the false-positive side only
— a pipeline that found nothing under the null *and* nothing under real
signal would look the same — which is why S4 still stands.

Note also that the Q4 count is soft: 40 here, 46 in `s4_diagnosis.py`,
the only difference being 60 vs 100 Phase 3 calibration permutations.
The pattern-level results carry the weight, not the count.

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

## Specification bugs and implementation defects found

All are documented with evidence in [`docs/DEVIATIONS.md`](docs/DEVIATIONS.md).
The first two were found by running planted ground truth through the
classifier; the third only appeared on real data.

**0. The Phase 2.6 guard was never armed (D10).** `saturation_flags()`
was implemented, documented and unit-tested — and never called, so no
published result had ever been checked for the ceiling/floor
conflation. Arming it retracts the TOX result and re-fires S4. Found by
an adversarial audit of this repo, not by me. Four further defects
(Q3 guarded on the noise-screened value rather than raw TOST; Q3 firing
when NK was merely indeterminate; "equivalent" reported as "equivalent
to zero"; voom using `log2(mean(lib))` where limma uses
`mean(log2(lib))`) are in [D10–D11](docs/DEVIATIONS.md).

**0b. The detection floor deletes the positive control (D7).** Phase 2.5
keeps only genes detected in *all five* lineages. TCR-proximal genes are
T-restricted, so `TOX`, `PDCD1`, `CTLA4`, `LAG3`, `TIGIT`, `ZAP70`,
`CD28` and `TNFRSF9` — 11 of the 14 Q4 controls — are removed before
testing, and stop rule S4 then fires for the only reason left: no
control survives to reach Q4. **A healthy pipeline fails S4 as written.**
Relaxing the filter to "detected in NK plus ≥ 2 lineages" (4,159 →
4,525 genes) brings them back into the panel. This was invisible in
simulation, where every gene is expressed in every lineage by
construction. It does not rescue S4, because of D10.

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

- **The 46 Q4 candidates are candidates, not findings.** The design is
  now shown to detect Q4 (76% recall), so they come from a working
  pipeline — but no individual gene is validated, and no second cohort
  capable of testing them exists. Phases 6–8 are
  the gates that decide whether any of them survive, and none has run.
- **Phase 6 (ambient) not run.** CellBender has not been run. The
  shared-soup HTO contrast — the protocol's cleanest control, and the
  strongest single asset of this dataset — is implemented but not
  executed.
- **Phase 7 (replication) not run.** GSE178341 is admitted but not
  downloaded.
- **Phase 8 is only partly run.** B4 (ceiling/floor) is now armed, B5
  (power) is excluded for the Q4 set, and B6 (ambient) is partly
  addressed — the Q4 pattern survives where ambient cancels
  (p = 2.2 × 10⁻⁸, [D13](docs/DEVIATIONS.md)), though that test is not
  independent of selection. **B1 (spatial exposure), B2 (NK turnover)
  and B3 (missing receptor) are unaddressed**, and B3 is the one that
  matters most: D10 showed the protocol's own positive control fails
  precisely because it is a B3 case.
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

## What actually holds up

Ordered by how much weight the evidence bears.

1. **The protocol's central statistical premise is correct and now
   quantified.** Without equivalence testing, a rare lineage yields an
   ~84% false Q4 rate; with TOST, the nominal 5%.
2. **~40% of cross-lineage sharing in a conventional design is
   technical** (shared-soup contrast, depth confounder ruled out).
   This is a real, reusable measurement about scRNA-seq design, largely
   independent of whether Q4 exists.
3. **CD45⁺ enrichment, not cohort size, is what makes a dataset usable
   for a rare lineage** — the 36-patient cohort is unusable and the
   10-patient one is not.
4. **Q4 sensitivity is measured, not assumed: 76% recall** on a
   constructed truth set, with a detection floor of ≈ 0.19 log2FC
   residual NK effect. Paired with the null control (0/12 permutations
   produce any Q4), the pipeline is shown to be both sensitive and
   specific for Q4.
5. **The pipeline does not manufacture Q4 from noise**, and the Q4
   pattern replicates continuously in an independent cohort
   (p = 7 × 10⁻⁶) and survives ambient cancellation (p = 2 × 10⁻⁸).
6. **Specification defects in the protocol**, each with the evidence
   that exposed it: D1 (Q3 absorbs Q4), D2 (`p50` costs the power),
   D3 (null/estimator mismatch), D7 + D14 (the Q4 positive control is
   unsatisfiable in principle), D10 (guard never armed), and the Q3
   positive control being mis-specified for a paired design
   ([CELL2023_Q3_CHECK.md](docs/CELL2023_Q3_CHECK.md)).

Against that: **no individual gene is validated**, no second cohort can
test the candidates, and B1–B3 of Phase 8 are untouched. The 46 Q4
candidates are a lead, not a finding — but the design behind them is now
measured rather than assumed.

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

**Find a Q4 positive control that can actually work.** S4 fires, and
until it clears, nothing this pipeline reports about Q4 is readable as
biology. D10 shows why the protocol's own choice cannot serve: a
TCR-driven gene is both T-restricted and near-absent from NK, so NK's
non-response is mundane explanation B3 rather than resistance.

A valid control needs three properties at once:

1. **NK demonstrably expresses it** — clear of the Phase 2.6 floor, so
   there is room to move.
2. **Its driver reaches NK** — a receptor or pathway NK actually
   carries, so a non-response is informative rather than trivial.
3. **The witness lineages measurably respond to that same driver** in
   this contrast.

Candidate families worth testing against those criteria: TGF-β target
genes (NK carries TGFBR2 and is known to respond), type I/II interferon
targets, and hypoxia/HIF targets — all shared drivers with NK-expressed
receptors, unlike the TCR module. The test is cheap: each is a named
gene set, and `scripts/s4_diagnosis.py` already has the machinery to
check where a control set lands.

If no such control passes, that is itself the answer — it would mean
this design cannot demonstrate Q4 sensitivity on this data, and the
protocol's stop rules should be honoured rather than worked around.

**Second priority, if a control clears:** Phase 8's B1–B3, which are
untouched. B3 (receptor not expressed) matters most, since D10 showed
the protocol's own positive control failed precisely as a B3 case — the
same trap will be waiting for individual candidates.

**Not a priority:** more discovery. The bottleneck is not candidate
count, it is that nothing validates the candidates already in hand.
