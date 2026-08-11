# The "tumour-NK signature is really Q3" hypothesis — refuted

**Status: dead.** Killed by route 3, this repository's own data. Recorded
in full because it was the most attractive lead in the project and it did
not survive.

## The hypothesis

Pan-cancer single-cell NK papers report tumour-infiltrating NK cells as
marked by **DNAJB1, HSPA1A, FOS, JUN** — which is, essentially exactly,
the Marsh et al. (2022) conserved dissociation-stress geneset. If that
signature is a universal (Q3) artefact shared by every lineage, then
reading it as NK biology is a category error, and demonstrating so would
need only the Q3 quadrant: the best-powered cell, already working, with a
shared-soup control in hand and no TOST or equivalence testing anywhere
in the argument.

The question reduces to one sentence: **do these genes move in B cells
and myeloid cells from the same libraries?**

## Route 1 — read the methods, not the abstract

*Pan-cancer profiling of tumour-infiltrating NK cells through
transcriptional reference mapping* (Nat Immunol 2024, PMC11291284):

- **No cross-lineage control.** The methods describe CellTypist-based
  annotation and give no check of whether NK-enriched genes also move in
  B, myeloid or T cells from the same libraries. The gap the hypothesis
  assumed is real.
- **But a dissociation control does exist, and it is a good one.** The
  authors state they "took several measures to rule out digestion
  artefacts", removed ambient RNA with decontX, and — decisively —
  showed the stress signature in **spatial transcriptomics data "that
  have not undergone any upstream tissue dissociation/digestion"**.
  Tissue that was never dissociated cannot carry a dissociation
  artefact. This is a different defence from the one the hypothesis
  attacked, and it is a legitimate one.
- **No within-patient paired tumour/normal comparison.**

The *Cell* 2023 panorama's control was different again: it compared
mitochondrial content between the `c6-DNAJB1` and `c7-NR4A3` clusters to
argue the stress phenotype was not a cell-quality effect. That addresses
quality, not cross-lineage sharing.

**Verdict on route 1:** the specific gap is real, but these papers are
less naive than the hypothesis assumed. The spatial control in
particular is a serious answer.

## Route 2 — has anyone already run this check on NK atlases?

Not found. Marsh et al. established the shared dissociation geneset in
microglia and brain tissue, not in NK or tumour. No publication was
found applying a cross-lineage artefact control to a tumour NK
signature.

## Route 3 — run it on GSE154826. **This is where it dies.**

Paired tumour vs adjacent normal, 27 individuals, all five lineages from
the same droplet pools, δ = 0.5:

| gene | NK | CD8T | CD4T | B | Myeloid |
|---|---|---|---|---|---|
| DNAJB1 | **+0.74 changed** | +0.33 indet | +0.39 indet | +0.11 indet | +0.18 equiv |
| HSPA1A | **+1.11 changed** | +0.33 indet | +0.09 indet | −0.14 indet | +0.06 equiv |
| FOS | **+0.87 changed** | −0.28 indet | +0.10 equiv | −0.12 equiv | +0.31 indet |
| JUN | **+1.19 changed** | +0.58 changed | +0.64 changed | +0.49 indet | +0.36 indet |
| HSPA1B | **+1.04 changed** | −0.00 indet | −0.09 indet | +0.03 indet | −0.04 equiv |
| JUNB | **+0.82 changed** | +0.21 indet | +0.35 indet | −0.04 equiv | +0.31 indet |
| EGR1 | **+1.39 changed** | −0.08 equiv | −0.05 equiv | **−0.85 changed** | **+0.92 changed** |
| HSPH1 | **+1.15 changed** | +0.41 indet | +0.41 indet | −0.01 equiv | +0.35 indet |

**NK moves for all eight. B and myeloid mostly do not.** Not one is
classified Q3. The signature is NK-predominant in this dataset, which is
the opposite of the hypothesis and is consistent with the published
interpretation.

The effect is not a power artefact in the convenient direction: NK
carries the *largest* standard error of the five lineages (0.113 vs
0.070–0.094), so it is the lineage where "changed" is hardest to reach,
and it is the one that reaches it.

## Why the protocol expected otherwise, and what that costs it

Phase 5.3 nominates the immediate-early genes as the **Q3 positive
control** — "解离应激/即刻早期基因，所有谱系共享". On this data that
expectation is simply wrong, and the reason is structural rather than
incidental:

**In a paired tumour vs adjacent-normal contrast, both sides were
dissociated by the same protocol, so the shared dissociation component
cancels in the within-individual difference.** What survives the
subtraction is the *differential* stress response between the two
tissues — and that difference can be, and here is, lineage-specific.

Two consequences, both awkward for the protocol:

1. **The Q3 positive control is mis-specified for this design.** It
   would be appropriate for a fresh-vs-dissociated comparison, not for a
   paired contrast where dissociation is common to both arms. Every IEG
   landed in `unclassified` in the Phase 5.3 control check; that was
   recorded at the time as "not fatal" without being diagnosed. This is
   the diagnosis.

2. **The protocol has no control for a Q1-shaped artefact.** Its whole
   defensive architecture assumes technical effects are universal, hence
   Q3, hence discardable. A lineage-specific differential artefact —
   tumour tissue needing harsher dissociation, NK being more fragile
   than B or myeloid cells — lands in Q1 and is indistinguishable from
   real NK-specific biology by any test in the protocol. That is a gap,
   and this gene set is a live instance of it: **nothing here
   distinguishes "NK mounts a genuine tumour stress programme" from "NK
   is differentially damaged by tumour dissociation".**

## Caveat on my own gating

An alternative reading is that my RNA-based NK gate preferentially
captures stressed cells on the tumour side, manufacturing the
NK-specific signal. This was not tested. It would need the sub-cluster
purity check that Phase 2.1 asks for and that has not been run.

## What remains usable

The falsification is clean and cheap, and the reasoning transfers: any
claim that a lineage-restricted tumour signature is real should be
checked against the other lineages in the same droplet pool, which is
one query once the pseudobulks exist. Here it argued *for* the published
interpretation rather than against it.
