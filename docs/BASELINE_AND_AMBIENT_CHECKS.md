# Three checks before CellBender, and what they changed

Ordered as specified. Check 1 could have invalidated the Q4 list; it
invalidated one reading of it instead. Checks 2 and 3 both came back
against the stated hypothesis and in a more useful direction.

---

## Check 1 — the ceiling objection: right about the genes, wrong about the category

**The objection.** `TIGIT`, `CD38` and `PRDM1` are constitutively high on
NK cells. So "NK +0.08 while CD8T +1.12" may mean *NK started where the
others finished*, not *NK resisted*. Structurally, Q4 would then be
biased toward genes with a high NK baseline, since a high baseline
leaves less room to move and makes a TOST-equivalent call easier.

### On the three named genes: confirmed, quantitatively

Absolute levels (median log2 CPM):

| gene | NK normal | NK tumour | witness normal | **witness tumour** |
|---|---|---|---|---|
| TIGIT | **6.17** | 5.82 | 3.73 | **4.68** |
| CD38 | **6.02** | 5.72 | 3.54 | **4.72** |
| PRDM1 | **8.00** | 7.47 | 6.53 | **6.93** |

**NK's starting point is 1.1–1.5 log2 units above where the witnesses
end up.** The earlier reading — "NK does not upregulate the exhaustion
programme" — is withdrawn. The parsimonious reading is **"NK already
expresses these constitutively; there is nothing to upregulate."**

### On Q4 as a category: not supported, by three tests

| test | Q4 | background | p |
|---|---|---|---|
| NK baseline (log2 CPM) | 6.33 | 5.82 | **0.00038** |
| NK baseline **minus witness baseline** | +0.329 | +0.332 | **0.64** |
| "NK normal ≥ witness tumour" | 87% | 83% | — |

Genome-wide correlation between NK baseline and |NK log2FC|:
**r = −0.017**.

That last number is the direct test of the proposed mechanism — high
baseline ⇒ small NK effect ⇒ spurious equivalence. **It does not operate
genome-wide.** And Q4 genes are *not* specifically genes where NK is
high relative to the other lineages (p = 0.64); the +0.33 NK-over-witness
offset is the same in the background, so it is a property of NK, not of
Q4.

Q4 genes *are* more highly expressed overall (p = 0.0004). That is a
**power** effect, not a ceiling effect: better-measured genes have
tighter CIs, so equivalence is easier to establish. It biases Q4 toward
*under*-ascertaining lowly-expressed genes, not toward false calls on
highly-expressed ones.

**Verdict: the top three genes are baseline artefacts and their
interpretation is withdrawn. The count of 46 does not need discounting.**

---

## Check 2 — the detection floor is not leaking; ambient defeats it

**The objection.** `PRF1` reaching a "B cells not equivalent" verdict
implies the §2.5 floor failed, since B cells do not express perforin.

**Measured — the floor is working.** B-cell detection rates, ≥10 counts
in ≥50% of pseudobulks:

| gene | B detection rate | passes floor? |
|---|---|---|
| GZMB | **0.93** | yes |
| GZMA | 0.91 | yes |
| CTSW | 0.87 | yes |
| GZMK | 0.87 | yes |
| PRF1 | **0.57** | yes |
| KLRG1 | 0.56 | yes |

These genes clear the floor legitimately. The relaxed filter (D7) does
admit some genes with an undetected lineage — 2.3% for B, 5.0% for
CD8T — but **not these**.

**The real mechanism is worse than a leak.** B-cell pseudobulks contain
≥10 granzyme counts in 93% of samples *because ambient puts them there*.
Soup carries abundant NK/T-derived granzyme transcripts; every B droplet
takes some; summed over hundreds of balanced cells it clears any
count-based floor.

> **A detection floor cannot separate "expressed" from "contaminated to
> a detectable level."** The filter is functioning exactly as specified
> and still admits pure soup.

**This reverses the recommended ordering.** Since the route is ambient
rather than filter logic, **CellBender addresses it directly** — remove
the soup, detection rates fall, and these genes get filtered out
correctly. There is no separate filter bug to fix first.

---

## Check 3 — NK loses RNA on the tumour side, and only NK

Measured in the 20 shared-emulsion libraries, where both conditions sit
in the *same* droplet pool, so any difference cannot be a soup
difference.

| lineage | median UMI, tumour | median UMI, normal | log2 ratio |
|---|---|---|---|
| **NK** | 1,156 | 1,611 | **−0.323** |
| CD8T | 2,089 | 2,069 | +0.061 |
| B | 2,917 | 2,038 | +0.612 |
| Myeloid | 4,016 | 3,423 | +0.212 |

**NK is the only lineage whose RNA content falls, and it falls
significantly (NK vs others, p = 6.0 × 10⁻⁸).**

This explains the residual that check 3 was written to explain. With
identical soup on both sides, ~25% less endogenous RNA on the tumour
side means **soup occupies a larger fraction of the tumour-side NK
transcriptome**, producing an apparent ambient-driven effect even though
the soup itself cancels. The 1.68 residual needs no "incomplete
cancellation" explanation.

### Two consequences beyond the calibration

**(a) It is NK-specific, so it is a Q1-shaped artefact source.** This is
the third independent instance in this project of a technical effect
that lands only on NK, and the protocol has no control for that shape.

**(b) Phase 2.4 balances cells, not RNA.** Power balancing equalises
*cell numbers* across lineages. Two groups with equal cell counts but
different UMI per cell have different effective depth, and that
difference falls specifically on tumour-side NK. The balancing is
incomplete in exactly the direction that matters.

**Direction for Q4:** a rising soup fraction pulls tumour-side NK
*toward* the soup, i.e. makes NK move, which breaks equivalence and
**destroys** Q4 rather than creating it. Q4 is conservative with respect
to this artefact. Q2 and the Q1/Q3 boundary are not.

---

## Revised order

1. ~~Check the detection floor~~ — done; the floor is sound, ambient
   defeats it, and CellBender is the right instrument.
2. **CellBender on the raw droplet matrices**, with the two-sided
   acceptance test: immunoglobulin → 0 in NK *and* granzyme/perforin
   → 0 in B. Both true values are known.
3. Re-run the quadrant assignment; only then is Q2 evaluable.
4. Re-examine the Q4 list for baseline artefacts **gene by gene** using
   the four absolute levels (`results/Q4_with_absolute_levels.csv`),
   since the category-level test cleared it but the three headline genes
   did not survive.
5. Over-dispersion — after 2, for the reason already established: soup
   is additive and depresses apparent dispersion, so pre-correction
   variance measurements would track the contamination gradient. Note
   that check 3 adds a second reason: NK's RNA content differs between
   conditions, which moves dispersion estimates on its own.
