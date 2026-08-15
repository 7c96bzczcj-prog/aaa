# X1 — how it would actually be done

**Scope note.** This document is **outside** the screening round. Spec §0 fixed
"no new directions" as a non-goal, and `docs/RESULTS.md` honours that. This file
exists because the design question was asked separately, after the screen closed.
Nothing here feeds back into T1–T5 or changes any grade.

**What the screen established, and is assumed here.** X1 = organ-resident trNK →
functional intratumoural NK. Evidence `INF_MARKER`; load-bearing HIGH under v1.1;
lineage-decidable TRUE; both answers actionable TRUE; `already_traced = FALSE` at
HIGH confidence; `mouse_only = TRUE`; Phase C power PASS on 4 donors.

---

## 1. Decompose the claim before designing anything

X1 bundles three assertions that need different instruments:

| clause | statement | instrument |
|---|---|---|
| **(a) provenance** | intratumoural NK descend from the NK pool of *that organ*, not from the circulating pool | clonal barcoding — decidable |
| **(b) temporal** | that ancestral pool was in the organ *before the tumour* | requires a pre-tumour label — **not obtainable in humans** |
| **(c) functional** | the intratumoural cells are *functional* | must be a **sort gate**, not an inference |

Clause (c) is the one most easily lost. "Functional" has to be imposed at sort
(GzmB⁺ / degranulating / IFN-γ⁺), otherwise the design silently answers a
different question — where all intratumoural NK came from, functional or not.

---

## 2. The discriminating statistic is coalescence *depth*, not sharing

Every NK cell in one person shares an HSPC ancestor, so bare clone-sharing is
trivially positive under every hypothesis and is not a discriminator. What
separates the hypotheses is **how recently** two cells coalesce:

- **X1 true** — tumour NK coalesce *recently* with uninvolved-organ NK, and only
  at HSPC depth with blood NK.
- **Recruitment true** — the reverse: tumour NK sit inside the contemporaneous
  circulating repertoire.

This forces four compartments, not two. Blood is a **mandatory discriminator**,
not a nice-to-have; and a marrow/HSPC-level baseline is needed to know what
"deep" means for this individual. Without blood, shared clones between organ and
tumour cannot be separated from common marrow output and *neither* answer is
reachable.

---

## 3. Tier 1 — the human analysis that is runnable now

> **CORRECTION (2026-08-15).** This section was written before the deposited data
> was checked for variant ascertainment. It is **wrong** that Tier 1 is nearly
> free. The deposited heteroplasmy matrices are per-library, not the donor union,
> and the variant set evaluable across all three compartments is ascertained to be
> blood-present (95.5–100%), which is the very confound the statistic measures.
> Tier 1 is not blocked but it is **expensive** — see `docs/X1_TIER1_RESULT.md`.
> The section is left standing rather than rewritten, because the error is part of
> the record.

`GSE302113` is already downloaded and parsed in this repo. Four NSCLC donors
(SU-L-001, -002, -004, -005) pass the preregistered power check with ≥30 NK in
both tumour and matched non-involved lung; all five have blood.

**The statistic**, using the authors' own definition so it is comparable to their
published myeloid results (G6 — same basis as the rule consuming it): fraction of
same-clone cell pairs between two cell sets, normalised by all possible
within-donor pairs. Compute for NK:

```
sharing(tumour-NK, NILT-NK)   vs   sharing(tumour-NK, blood-NK)
```

with within-compartment sharing as the shallow baseline, and with **myeloid and
CD8 T as concurrent controls** — those are the lineages the paper already
characterised, so they calibrate the statistic on the same donors.

**Clone calling must follow the source pipeline**: mgatk union across the donor's
libraries (a variant passing filters in any one sample counts in all samples of
that donor), heteroplasmy binarised at 0.07, variants in >20% of a donor's cells
or <3 cells removed, and chrM:307–314 excluded as a homopolymeric artefact
region. Substrate already verified: 32–319 informative variants shared across ≥2
compartments per lung donor (DEC-09).

**One number already exists and was never reported.** The X1 adversarial pass
found that Liu 2026 Figure 2F plots tumour-cell-type × NILT-cell-type clone
sharing, and NK/ILC is an annotated cell type in that atlas. The tumour-NK ×
NILT-NK entry is therefore **plotted but never quantified, tested, or
interpreted** — the accompanying text discusses only CD8 T, B/plasma and myeloid.
Quantifying it is the cheapest possible first move.

### What Tier 1 can and cannot deliver

**Can**: adjudicate clause (a), provenance, against a blood baseline.
**Cannot**: touch clause (b). Non-involved lung is *contemporaneous*, not
pre-tumour, and treating it as a proxy assumes the uninvolved region was not
remodelled by the tumour-bearing state — an assumption the same dataset gives
reason to doubt, since it reports NK depleted in tumour relative to NILT.

**Three limits to state up front rather than discover later:**

1. **NK typing.** The cell-type annotations are not deposited. My cluster-level
   calls (DEC-08) were validated only on *blood* (CD56⁺ 74.6% vs CD56⁻ 5.5%) and
   are explicitly a power estimate, not a tissue-validated classifier. A real
   analysis needs better NK calls in tissue — either the authors' annotations, or
   a classifier validated where ILC1/trNK boundaries actually bite.
2. **Assignment rate.** mtDNA clone assignment runs ~15–26% in solid tissue vs
   ~20–32% in PBMC, so per-donor NK clone counts are small and a null result at
   n = 4 would be uninformative rather than negative.
3. **The anti-correlation.** Where matched normal tissue exists, tumour NK are
   thin (281–509 cells); where tumour NK are abundant (ovarian SU-O-005, 2,476)
   there is no matched normal organ. Two independent routes found this — my power
   table and the Phase B agent reasoning from literature alone.

**So Tier 1 is the cheap move, not the decisive one.** It is worth running because
it is nearly free and because a strong negative would matter; it cannot settle X1.

---

## 4. Tier 2 — the decisive experiment, and why it is mouse

Clause (b) needs a label applied to the organ **before the tumour exists**. No
human sample exists of an organ from before its own tumour arose — not
post-mortem, not at surgery. That is a logical barrier, not a logistical one, and
it is why Phase B returned `mouse_only = TRUE`.

**mtDNA will not work here.** Preregistration §4.3: somatic mtDNA variants
accumulate with age, and laboratory mice are young and isogenic, so informative
variants are scarce. Tier 2 must use a non-mtDNA recorder.

**Design.**

- **Model**: autochthonous, so the tumour arises *in the labelled organ* — e.g. a
  KP lung model. A transplanted flank tumour cannot address X1 at all, because
  the organ-resident pool it would need was never there.
- **Label**: pre-tumour marking of the organ-resident pool. Either
  photoconversion (Kaede / KikGR) of the target lobe before induction, or
  inducible `Ncr1-CreERT2 × Confetti` with tamoxifen given **after** the resident
  pool is established and **before** tumour induction. The label must be applied
  below the progenitor stage, or both hypotheses predict the same labelling.
- **Readout**: fraction of functionally-gated intratumoural NK carrying the
  pre-tumour organ label, against a concurrently labelled circulating pool.
- **Mandatory control**: intravascular anti-CD45 at harvest, to exclude
  vascular-contaminating cells being scored as intratumoural. Its absence is a
  standard failure mode in this literature.
- **Registry-internal**: the same design reads out X4. If the labelled,
  pre-existing pool contains both GZMK⁺ and ENTPD1⁺ cells at the outset, Serger's
  trajectory is pre-existing resident heterogeneity read as a trajectory.

**Then the branch**, taken from the Phase B record:

- **If YES** — convert the label into a functional handle: DTR restricted to the
  marked resident pool, deplete at t₀ before tumour induction, and read tumour
  incidence and growth against unlabelled controls and against a
  recruitment-blocked arm that leaves the resident pool intact.
- **If NO** — the target compartment moves to the circulation: serially label
  blood/spleen at staggered times relative to tumour onset and ask whether a
  *single* recruited clone gives rise to both trNK-like and conventional
  intratumoural progeny. That question is not asked on the YES branch at all, and
  it is the selection-versus-conversion fork.

---

## 5. What would falsify X1

Stated in advance, so the result is not read off after the fact:

- Intratumoural NK coalescing with blood NK at the same depth as blood NK
  coalesce with each other, while showing no shallower coalescence with
  uninvolved-organ NK than the HSPC baseline predicts.
- In mouse: pre-tumour organ label essentially absent from functionally-gated
  intratumoural NK, while the concurrently labelled circulating pool is well
  represented.

Note what neither result would show: that resident NK are irrelevant to tumour
control. Ablating them could still change tumour growth even if they contribute
no progeny to the intratumoural pool.

---

## 6. Honest summary

Tier 1 is cheap, uses data already on disk, and quantifies a number that exists
in a published figure but was never reported. It addresses provenance only.

Tier 2 is the experiment X1 actually needs, is mouse-only for a reason that no
amount of human sample collection can fix, and is the one the field currently
states it lacks the tools to do cleanly.
