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

## 5. Provenance

Derived from `docs/DECISIONS.md` D22–D27 and `PREREGISTRATION_v1.1`–`v1.4`,
which carry the measurements behind every claim above.
