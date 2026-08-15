# TRANSFERABLE — rules that generalise beyond this round

**Note on numbering.** The project's canonical TRANSFERABLE list, including §5
(controls are chosen by the biology) and §7b (a summary quantity is computed on
the basis of the rule that consumes it), lives outside this repository. Those two
are referenced here by the numbering used there; their text is **not** reproduced,
because reproducing it from memory would risk drifting from the original. This
file records the rules this round produced, for merging into that list.

---

## Data availability must be checked against the *statistic*, not against the *modality*

**Rule.** "The files exist and clear the threshold" and "their ascertainment
structure can carry this question" are two different checks. Only the second one
licenses an analysis. Run it before costing the work, not after.

**Where it came from.** GSE302113, X1 Tier 1 (`docs/X1_TIER1_RESULT.md`, DEC-21).

Everything modality-level passed:

- chrM retained in **38/38** libraries, against a NUMT-masked reference
- median chrM coverage **23.8×–158.9×**, every library over the preregistered 20× floor
- **32–319** informative variants shared across ≥2 compartments per donor
- power check passed on **4 donors**

And the analysis was still impossible. What was missing sat one level up from the
modality: the deposit carries **per-library selected** heteroplasmy matrices, not
the **per-donor union** the authors' own methods describe building. A variant is
evaluable only where it was deposited, and it is deposited only where it passed a
≥5-cell filter — so the variant set evaluable across all three compartments is
**ascertained to be blood-present** (0.0% of one donor's 507 tri-compartment
variants have zero blood carriers). The statistic asked "are these clones absent
from blood?" and the usable variant set answered "no" by construction.

**The specific trap.** This is the chrM-stripping trap one level up. The prior
round was burned by archives that strip chrM; here chrM is fully present and the
*variant union* is what is missing. The general form: **an archive can preserve
the measurement and discard the ascertainment**, and only the ascertainment
determines whether a comparison is possible.

**How to run the check.** Take the statistic's discriminating quantity, and ask
what would have to be true of the deposit for both of its answers to be reachable.
Here: "clone absent from blood" is only reachable if the variant can be evaluated
in blood cells that do not carry it. Then verify that on the actual files —
counting, not reasoning from filenames.

**Corollary — a lossy deposit is not visible from a file listing.** Two deposits of
the same assay differ decisively:

| deposit | what it permits |
|---|---|
| mgatk object (per-cell allele counts at all positions) — e.g. GSE197037 | any variant evaluable in any cell; no ascertainment confound |
| derived per-library heteroplasmy at *selected* variants — GSE302113 | only the selected variants, in only the libraries that selected them |

Both list as "mtDNA data deposited". Only the first supports a cross-compartment
question.

**Second corollary — this applies to the pattern-matching too.** The first pass of
the Route D re-audit classified GSE197037 as depositing nothing, because the
filename regex did not include `mgatk*.rds`. That is the same error one more level
up: trusting a pattern over the actual list. The fix was to print the filenames and
read them. Recorded in the script itself
(`scripts/routeD_variant_level.py`) so it is not lost.

---

## Load-bearing is bidirectional

**Rule.** A claim is load-bearing when the field presupposes it TRUE **or** when
the field presupposes it FALSE. Work built on a claim being false must be re-read
if it turns out true, and that is load in the same sense.

**Where it came from.** `docs/PREREGISTRATION.md` §A1, DEC-17. The one-directional
version graded X1 MEDIUM on the ground that its only positive dependency was
origin-agnostic; the negation side turned out to carry four primary papers plus an
engineered-cell therapeutic strategy.

**Two properties worth keeping.**

1. **It is monotone.** The bidirectional criterion is a union, so it can only raise
   a grade. Therefore the set of rows a one-directional pass could have got wrong
   is exactly *the rows it graded below HIGH* — a set defined by the criterion, not
   by the desired outcome, and so complete and unbiased.
2. **A single dependency can sit on both sides of two adjacent claims.** The
   dominant-negative TGFBR2 NK construct presupposes T1 **true** (TGF-β converts
   infused conventional NK) and X1 **false** (intratumoural trNK are
   infiltrating-and-converted rather than pre-resident). When claims are
   near-complementary, grading them independently will double-count or miss.

**Limiting clause, learned the same day.** When a claim's negation underwrites a
therapeutic construct, what the claim turning true collapses is usually the
construct's *narrative about the endogenous compartment*, not the construct
itself — infused cells are infiltrators by construction regardless.

---

## Don't change an admission threshold and an admission result in the same motion

**Rule.** Recalibrating a threshold between rounds is legitimate — an uncalibrated
instrument should be fixed. Recalibrating it *in the same change that admits or
rejects a result it bears on* is not, in either direction.

**Where it came from.** The 20× chrM floor (DEC-04, DEC-14) is known to be
uncalibrated: it would exclude Rückert 2022, a field-defining application of the
method, which operates at 11–20×. It was left untouched through the v1.1
amendment precisely because that amendment was admitting X1, and X1's Phase C
result depends on that floor.
