# Closing the three deferred leads

Deferred items consume attention silently. All three are resolved here:
two closed, one downgraded to a candidate node inside a live project.
Total cost, as scoped: about twenty minutes. Two of the three came back
differently than expected, and both differences are recorded.

---

## 1. NK-CHIP — closed, but not for the stated reason

The lead was deferred as blocked by a **methodological wall**: that
calling an NK "clone" needs single-cell lineage identity, that T cells
get this by expanding into a clone before sequencing, that NK cells do
not expand, that bulk duplex-seq gives mutational burden without clonal
structure, and that mtDNA tracing gives structure without nuclear
drivers. On that reading the blocker was that **no method exists**, so
waiting is not a plan.

**One search was authorised, with a pre-registered rule: if someone has
solved the wall, the direction dies. It fires.**

### The wall does not exist

**Rodriguez-Sevilla et al., *Nature Communications* 2025**
(doi:10.1038/s41467-025-58662-0, PMC11992119), on NK dysfunction and
immune escape in early-stage MDS/CCUS, ran **Mission Bio Tapestri
single-cell targeted DNA + surface-protein sequencing** and reports:

> "high-throughput single-cell targeted DNA and surface protein
> sequencing analysis using the Mission Bio Tapestri Platform showed
> that, independently of the genetic alterations, NK cells isolated from
> five CCUS patients were part of the mutant clone population and had a
> clonal burden similar to that detected in myeloid cells"

That is NK single-cell nuclear-driver genotyping with immunophenotype
attached — the thing the wall said does not exist. **The premise that
made it a wall is a T-cell-era assumption**: droplet single-cell DNA
genotyping reads the genome of individual cells directly and never needs
clonal expansion, so "NK cells do not expand" stops being an obstacle.
The platform is commercial and the NK result is published.

### What the deferral got right, and what it got wrong

The **Savola** paper is often cited as the "someone with the tools tried
and could not." Its actual sentence, from the supplementary methods:

> "However, the yield of the fraction containing NK cells was so small
> that NK cells were not sequenced."

They FACS-sorted CD3⁻CD19⁺, CD3⁻CD14⁺ and **CD56⁺CD16⁺** fractions and
ran **amplicon sequencing** on bulk sorted DNA. So the failure is
**sorted-cell yield in a bulk amplicon design** — a mundane input
problem — not the clonal-structure wall, and not duplex-seq. The
outcome ("NK went unsequenced") was read correctly; the cause was not.

### Verdict

**Closed — and this is a stronger death than the one proposed.** Not
"infeasible until a method is invented" but **"the method exists and the
headline result is already published by someone else."** A residual
exists (that paper establishes mutant-vs-wild-type membership and clonal
burden in NK, not NK-specific clone-size hierarchy or phylogeny), but
that is an incremental extension of another group's result, on their
platform, in the aging/haematology domain that is already off-command.
No further work.

---

## 2. NK complosome — premise checked on data already on disk

The premise "does NK transcribe intracellular C3, and carry `C3AR1`"
had never been run despite being answerable in ten minutes from
`results/pseudobulk_uncorrected_matched.npz`.

**It cannot be answered by a detection call**, because `C3` is a
canonical myeloid/stromal transcript and this project has measured that
NK carries the second-highest ambient burden of any lineage. Detection
in NK is exactly what contamination produces. Two scales are therefore
calibrated **inside NK, from genes of known truth**
(`scripts/complosome_premise.py`, `results/complosome_premise.csv`):

| scale | genuinely-NK genes (KLRD1/NKG7/GNLY) | fully-ambient genes (IGKC/C1QA/LYZ) |
|---|---|---|
| estimated soup fraction | 0.012 – 0.013 | **0.402 – 1.000** (true value 1.0) |
| NK/Myeloid CPM | — | **0.070 – 0.118** |

The soup-fraction ceiling is 0.402, not 1.0, because run 2 of the
ambient work established that this estimator reads ~0.29 in NK for genes
whose true soup fraction is 1.0. **A raw "0.19" therefore does not mean
"19% ambient"; it means "below where a fully-ambient gene reads."**

### Result

| gene | NK detection | NK CPM | soup frac (ceiling 0.402) | NK/Mye (pickup band ≤0.118) | call |
|---|---|---|---|---|---|
| **C3** | 0.69 | 16.2 | **0.191** | **0.437** | clears both scales |
| **C3AR1** | 0.91 | 20.7 | 0.388 | 0.208 | clears both, **marginally** |
| CD46 | 0.98 | 68.5 | 0.148 | 0.929 | clears, but uninformative |
| C5AR1 | 0.88 | 21.7 | 0.471 | 0.132 | not separable |
| CTSL | 0.90 | 90.6 | 0.476 | 0.188 | not separable |
| CFB | 0.28 | 0.7 | — | — | below detection |

**The premise is supported, weakly.** `C3` clears the ambient
calibration on both independent scales — soup fraction 0.191 against a
fully-ambient reading of 0.402, and NK/Myeloid 0.437 against a
pure-pickup band of 0.070–0.118, i.e. ~4× more C3 in NK than ambient
pickup accounts for. `C3AR1` clears too but only just (0.388 against a
0.402 ceiling).

Three limits, stated rather than buried:

- **Both are low-abundance.** 16 and 21 CPM, against KLRD1's 1050. This
  is a detectable transcript, not an abundant one.
- **`CTSL` and `C5AR1` must not be read as negatives.** The NK/Myeloid
  ratio test only discriminates for genes myeloid expresses far more
  than NK; for broadly-expressed genes it flags ambient by construction.
  Those two are *untestable this way*, not refuted. This matters,
  because cathepsin L is the protease in the proposed mechanism.
- **`CD46` proves nothing.** It is expressed on all nucleated human
  cells; it is included only to show the scales behave.
- Transcript presence is also not the whole premise: complosome C3 can
  be taken up from serum rather than transcribed, and intracellular C3
  detection has a documented antibody-specificity dispute.

### Verdict

**Not revived as a direction** — it is a wholesale turn into complement
immunology, in a paradigm the Kemper group owns, and off-command.
**Downgraded to a candidate node**, as proposed, on a real anatomical
overlap: NK lytic granules *are* secretory lysosomes, and complosome C3
sits in lysosomes and is cleaved by cathepsin L. `C3` and `C3AR1` join
the Aim 2 candidate panel alongside DGKα, Cbl-b and Cdc42–Par6. The
premise check above is what licenses that, and it is now on file rather
than assumed.

---

## 3. The over-dispersion axis — dead, and the death certificate

Previously deferred as "prior-art check not deep enough, not invested
in." No prior-art check is needed: **this project's own data kills it.**

Dispersion estimates in NK carry two lineage-specific confounders, and
both push toward the false-positive direction:

| confounder | measurement | effect on NK's apparent dispersion |
|---|---|---|
| ambient is **additive**, so it inflates the mean while adding little variance | NK ρ = 0.151, **2.6× CD8T's 0.058** | depresses apparent dispersion **specifically in NK** |
| NK's **RNA content differs between the two conditions**, and sequencing depth drives dispersion estimation directly | median UMI 1611 → 1156, log2 −0.323; the only lineage that falls (CD8T +0.061, B +0.612, Myeloid +0.212), **p = 6.0 × 10⁻⁸** | makes the tumour-vs-normal dispersion contrast a depth contrast |

**This is the same death as the four-quadrant framework's.** The
technical effect is not shared across lineages; it lands on NK, and it
points in the direction that manufactures the finding being looked for.
An axis whose readout is systematically biased in exactly the lineage
and exactly the direction of interest cannot be rescued by careful
analysis.

Nor is there a route out:

- **Decontamination cannot fix it.** CellBender is unaffordable here
  (322 s/epoch, ~300 h for 73 libraries), and the two affordable
  corrections both failed their own acceptance tests
  ([AMBIENT_CORRECTION_ATTEMPT.md](AMBIENT_CORRECTION_ATTEMPT.md)).
- **Even a perfect decontamination would not fix it.** The RNA-content
  confounder is biology and technique entangled — a real difference in
  NK's transcriptional output on the tumour side is indistinguishable
  from a capture-efficiency difference, and removing ambient does not
  touch it.
- **A different dataset cannot fix it.** The GEO rescan found exactly
  one cohort meeting CD45⁺ + paired + ≥20 individuals, and it is the one
  in use.

**Verdict: dead. Cause of death recorded above.**

---

## Net

| lead | prior status | now |
|---|---|---|
| NK-CHIP | deferred (platform) | **closed** — the wall does not exist; the result is published |
| complosome | deferred | **closed as a direction; `C3`/`C3AR1` added to the Aim 2 candidate panel** on a passed premise check |
| over-dispersion | deferred (prior art) | **dead** — killed by this project's own NK-specific confounders |

Zero deferred items remain.
