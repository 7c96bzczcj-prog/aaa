# What transfers out of this project

The scientific output here is one reproduction, one negative result, and a
reusable pipeline. **This file is the part with the widest reach**, and it is
written down because it would otherwise be lost in five rounds of review
threads.

Read it as: *structural claims transfer; specific numbers do not.*

---

## 1. The failure form that recurred four times

Four instruments were withdrawn rather than tuned across v1.1–v1.4:

| instrument | its output was actually determined by |
|---|---|
| BH over the signed-rank p | **n** — the floor 2/2ⁿ made rejection impossible at n ≤ 6 |
| ruler B's fixed myeloid denominator | **the denominator lineage** — ~1 for any gene myeloid doesn't source |
| the FDR criterion | **the null-set size** — the same data gave 6, then 0, then 8 passes |
| z as a standalone ranking | **the baseline rate** — z = 5 buys 32.1 pp at >30% detection, 1.63 pp at <0.5% |

**The common form: the output was a function of a nuisance parameter rather
than of the effect.** That is diagnosable *before* running anything — vary the
nuisance parameter over its plausible range and see whether the verdict moves.
If it does, the instrument is not measuring the effect.

**What survived all four withdrawals** is exactly what never used any of them:

- direction agreed by **every** donor (a sign count, no distributional model)
- effect size in the **natural units of the readout**
- **reproduction under a grouping that could not have produced the effect**
  (here: within a fixed cluster)

That triple is the thing to design *toward*, not the thing to fall back on.

---

## 2. The bounded-readout pathology, and how to pre-empt it

A readout bounded in [0,1] — detection rate, %positive, fraction of anything —
has variance p(1−p)/n, which collapses at both ends. Consequences:

- near **0**: tiny absolute changes are enormous relative changes. Any
  standardised statistic (z, fold-change, SNR) inflates. **False power.**
- near **1**: differences compress toward zero and a null is unreadable.
  **No power, and a null that looks like a finding.**
- **both happen in the same panel at once**, and in this study four of the
  fourteen floor-band genes were exactly the ones producing the largest and
  most misleading z values.

**No transformation fixes this.** A logit pushes the reading *further* toward
low-baseline items, because it expands precisely the compressed region. It is
not a statistical defect to repair; it is a choice of scale to declare.

### The pre-emption, in three steps, all before seeing results

1. **Declare the measurable band.** Compute each item's baseline level in the
   *reference* condition and bin it: floor (no power), measurable, ceiling
   (uninterpretable). Publish the lists.
2. **Fix the effect scale in advance.** Points? Fold? Something else? Pick one,
   commit it, and forbid any table sorted by the auxiliary statistic.
3. **State what a flat result means per band** — *untested* at the floor, a
   real null in the middle, *uninterpretable* at the ceiling. These are three
   different findings and must never share a word.

---

## 3. The same pathology in wet-lab readouts

**Flow cytometry %positive is the identical object.** A marker at 2% in the
control condition has minute variance, so any change is a large multiple —
that is S1PR5's shape in a FACS plot. Everything in §2 applies unchanged.

Checklist, all fixed **before unblinding**:

- [ ] which planned readouts sit near 0% or 100% in the **control** condition?
      Declare those unmeasurable in advance, exactly as the floor/ceiling
      bands are declared here
- [ ] effect scale committed in writing — percentage-point difference? fold?
      MFI? — before the first batch is looked at
- [ ] gating strategy frozen before unblinding
- [ ] **n is animals or donors, never wells.** The wet-lab form of R1. A well
      is a technical replicate; it inflates df without adding information

The nuisance parameters to stress-test in that setting are **replicate count,
normalisation baseline, and readout scale**. Vary each over its plausible
range and check whether the verdict moves.

---

## 4. Numbers do not transfer. A worked example.

`CCR5` sits in this study's floor band: **68 detected of 10,130 depth-matched
decidual dNK1–3 cells = 0.67%**, against a 5% band threshold.

That is a fact about **human first-trimester decidual NK cells on 10x v2 at a
2,137-UMI floor.** It says nothing about *Ccr5* in mouse tumour-infiltrating NK
cells, and it must not be used as if it did. Expression of a chemokine
receptor is not conserved across species, tissue, or activation state, and
inferring one from the other is the same move — *correct general principle
applied directly to a specific number without checking that number* — that
produced two errors in this project's own review (the PAEP ceiling argument
and the CCL3 band assignment).

**What transfers is the check, not the value.** For any dataset where a
conclusion rests on a detection-rate readout:

| step | action |
|---|---|
| 1 | compute the gene's mean detection rate in the **reference** condition, at matched depth |
| 2 | **< 5%** → floor band: the readout has no power in points, and any standardised statistic on it is inflated. A published effect from this band needs re-reading before it is relied on |
| 3 | **> 85%** → ceiling band: differences are compressed; a null there is uninterpretable, not negative |
| 4 | in between → the readout is usable; report the effect in points |

This check is one number and one comparison. **Run it on any archived result
whose conclusion depends on a detection rate before that result is relied on
again** — the cost is minutes and the outcome is binary.

---

## 5. Ambient controls must come from prior knowledge, not from a statistic

Late addition, and the most portable single rule in the file.

An ambient control must be a gene the target population **cannot transcribe**,
so that its entire signal in the target gate is pickup.

**What fails is data-driven selection, not a region of the lineage tree.**
Prior-selected controls remain valid wherever the exclusivity claim is solid,
including on the target's own side: "NK does not transcribe immunoglobulin"
(IGKC), "NK does not transcribe lysozyme" (LYZ), "NK does not transcribe
complement C1q" (C1QA) are as hard as "NK does not transcribe collagen". Those
three carry ruler A on the haematopoietic side and it remains usable there.

**What fails is auto-expanding the control set with a statistic**, and on the
target's own side it fails *always*, because separating *"the target expresses
it"* from *"the target picks it up"* is exactly the quantity the ambient
estimator measures. Any screen for the first shares a term with the second:
here `soup_fraction = rho × amb_cpm / nk_cpm`, and the obvious screen
"NK/source is low" divides by that same `nk_cpm`, so screening for low NK
expression necessarily selects high soup fractions.

Concretely: auto-selecting haematopoietic ambient controls for an **NK** gate
returned `IFNG` (an NK effector), `CCL3` (a study target, ~90% detected in NK)
and a T/NK lncRNA. Their soup fractions were low, correctly, because NK
expresses them — and they were briefly allowed to set an ambient ceiling. The
non-haematopoietic side had no such problem, because no prior belief permits an
NK cell to transcribe collagen: there, auto-selection crosses a boundary the
target cannot cross, and is admissible.

**Rule.** Ambient controls come from prior exclusivity knowledge. Data-driven
selection is admissible only across a boundary the target cannot cross. Do not
over-correct this into "the estimator is unusable on the target's side" — the
prior-selected controls are unaffected.

---

## 6. The same rule in the wet lab: negative controls are chosen by biology

The isotype or fluorescence-minus-one control is chosen because the biology
says the signal cannot be there — **never because a channel happened to be
quiet in the data.**

**The loading-control case is the one most likely to bite in metabolic work.**
`GAPDH` is a glycolytic enzyme. Using it to normalise an experiment that
manipulates the glycolysis/OXPHOS balance is this error in its purest form:
it is chosen *because it looks stable*, and whether it is stable is the very
quantity being perturbed. `ACTB` moves under metabolic and hypoxic conditions
too. The same applies to the **normalisation basis** of a flux assay — if the
perturbation changes cell size or protein content, normalising per protein
subtracts the effect. Commit the basis before the first plate is read, not
after comparing which one looks cleaner.

Carry that straight into an experiment: **the negative control must be a
pathway that should not respond, or a readout that should not move, named in
advance for a biological reason.** Choosing "whichever readout came out flat"
is the wet-lab form of the circularity in §5 — the selection statistic is the
measurement, so the control is guaranteed to look clean and guarantees
nothing.

Same shape, five settings:

| setting | the wrong move | the right move |
|---|---|---|
| ambient RNA | pick controls a statistic says are quiet in the target gate | pick genes the target lineage cannot transcribe |
| flow cytometry | pick the channel with least signal as "background" | isotype/FMO chosen from the staining biology |
| perturbation | call whichever readout stayed flat the negative control | name a should-not-respond pathway before unblinding |
| **Western / qPCR loading control** | pick the housekeeper that happens to be **flattest across your lanes** | pick one the perturbation **cannot** affect, on biological grounds |
| **Seahorse / assay normaliser** | choose cells vs protein vs DNA after seeing which looks cleanest | commit the basis before the first plate is read |

---

## 7. A reported summary must share its definition with the rule that consumes it

Second instance of one failure form, so it earns its own rule.

> **Any summary you report must be computed on the same definition as the rule
> that acts on it — otherwise it will be read *as* that rule.**

Both instances in this project:

| what was reported | what the rule actually used | result |
|---|---|---|
| BH q "reported but not adjudicating" | nothing — the criterion had been withdrawn | the number was still read as a verdict, so it had to be deleted, not demoted |
| a gene's **cross-subset mean** detection (0.817) placing it in the measurable band | the ceiling rule is **per arm**, and one arm was at 0.884 | a conclusion was archived as a "powered true negative" when one arm was saturated |

The two look unrelated and are the same mistake: **a reporting layer and a
decision layer computed on different definitions.** A summary printed beside a
result acquires that result's authority whatever the caption says.

Two practical forms:

- if a rule is withdrawn, **delete its number**, do not demote it to
  "informational";
- if a rule is per-arm / per-group / per-replicate, **report it that way**, not
  as an average across the very units it distinguishes.

---

## 8. Weak observations combine; report the set, not the line

When each observation is individually weak — small n, a sign count, a bounded
readout — the evidence lives in **whether the observations are mutually
consistent**, not in any one of them.

Worked instance. Each within-label contrast here had n = 3, so a 3/3 sign
concordance has null probability 2/2³ = 0.25 — nothing on its own. But the
paralogue pair counts once, leaving three independent observations, two
surviving and one *reversing* in the predicted place: 0.25³ ≈ 1.6%.

Three practical consequences:

1. **Count independence honestly.** Paralogues, correlated readouts and
   repeated measures on the same subject are one observation, not several.
   Inflating the count here is the same error as counting wells as n.
2. **Never quote a member as if it were the set.** "3/3 donors agreed" invites
   the reader to weigh one line at 0.25. Quote the pattern and name its
   exceptions.
3. **An exception is informative only under the same premise.** If the set is
   what carries weight, a member that breaks the pattern is meaningful only
   when it had the power to conform. Check that before reading a departure as
   a negative — here the one exception turned out to be ceiling-limited, so it
   carried no information at all.

The stronger version of the same idea is a *relationship* between observations
rather than a tally: two genes moving in **opposite** directions in one
contrast on one set of cells excludes any artefact that pushes everything one
way. Design for relationships that a nuisance parameter cannot produce.

---

## 9. Provenance and versioning of the numbers above

Derived from `docs/DECISIONS.md` D22–D32 and `PREREGISTRATION_v1.1`–`v1.4`.

**Two quantities appear at more than one value across those documents. Both
are correct at their own stage; check the stamp before quoting.**

| quantity | value | which reading |
|---|---|---|
| CCL3 detection in dNK | **0.817** | dNK1–3, depth-matched to 2137 UMI, per-subset mean — **the band assignment uses this** |
| | 0.909 | all dNK including dNKp, raw counts, no depth matching — quoted in D31 |
| ruler B pickup band | **0.0538** | v1.2 onward: source lineage ranked by **total counts** (LYZ→Myeloid) — **current** |
| | 0.0441 | v1.1: source ranked by **per-cell CPM** (LYZ→cDC1) — superseded, retained in `PREREGISTRATION_v1.1.md` as the record |
| CCL3 cross-subset summary | 0.8167 | every donor×subset row equally weighted — what v1.4 A11 used |
| | 0.8204 | each subset averaged first, then subsets equally weighted (dNK3 clears in 5 donors, dNK1/dNK2 in 6) |

Neither CCL3 summary enters a decision: the ceiling rule consumes the per-arm
figures (0.800 / 0.777 / **0.884**). See §7 — that is exactly why the summary
should not have been quoted as if it did.

Depth matching lowers a detection rate and it is not optional for anything
comparative; the raw figure is only ever a sanity check.
