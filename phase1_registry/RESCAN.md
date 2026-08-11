# Dataset rescan under the corrected criterion

Phase 1 originally screened by intuition — big cohort first, protocol
second — and looked at three datasets. The measured result of this
project inverts that ordering:

> **What determines whether a dataset can test a rare lineage is CD45⁺ /
> immune enrichment, not cohort size.** Pelka (GSE178341) has the most
> paired patients of any candidate (36) and yields **5** usable ones,
> because unsorted tissue leaves NK at 1.1% of cells.

So the landscape was rescanned with enrichment as the first filter.

## Criterion

1. CD45⁺ sorted or otherwise immune-enriched (or unsorted with a
   demonstrated high immune fraction)
2. Tumour **and** matched normal from the **same individual**
3. **≥ 20 individuals** meeting both

## Method

NCBI eutils search of GEO DataSets, three keyword strategies covering
CD45/immune-enrichment, paired/adjacent-normal, and tumour terms,
restricted to *Homo sapiens* expression profiling — **114 unique GSE
series**. Series were screened on title and summary, then the promising
ones were verified by pulling sample-level titles and counting
**individuals**, not samples. Raw query output and the candidate list
are in `phase1_registry/rescan/`.

## Result — no second qualifying cohort exists

| dataset | enrichment | paired | **individuals** | verdict |
|---|---|---|---|---|
| **GSE154826** Leader NSCLC | **CD45⁺** | ✓ | **29** | **the only one that qualifies** |
| GSE140228 Zhang HCC | CD45⁺ | ✓ | ~10 | too few |
| GSE131907 Kim NSCLC | unsorted | ✓ | 10 | measured: 8 usable at a relaxed threshold; cannot test Q4 (D14) |
| GSE114725 Azizi breast | CD45⁺ | ✓ | 8 | too few |
| GSE178341 Pelka CRC | unsorted | ✓ | 36 | measured: 5 usable; NK = 1.1% of cells |
| GSE181061 RCC | CD45⁺ | ✓ | **2** | too few |
| GSE139555 | T cells only | ✓ | — | fails A3: no B or myeloid comparator |
| E-MTAB-6149 Qian pan-cancer | unsorted | ✓ | 36 across **4 cancer types** (~9 each) | inferred to fail on both counts; **not measured** |

The combination of CD45⁺ enrichment, within-individual pairing, and
≥ 20 individuals is met by **exactly one** public dataset, and it is the
one already in use.

## What this licenses, and what it does not

**Licensed:** the statement *"no currently available public dataset can
replicate the Q4 candidates"* — and, more usefully, the reason. The
binding constraint is not cohort size but enrichment protocol, and it is
a property of how the data were generated rather than of the analysis.
It also explains, rather than merely reporting, why Phase 7 is stuck:
GSE131907 could not test Q4 (only 6 genes there admit a synthetic Q4
construction, D14), and no better cohort exists to replace it.

**Not licensed:** a claim that no such dataset exists anywhere. This
scan covered **GEO only**, through keyword matching on titles and
summaries. It cannot see:

- ArrayExpress, EGA, dbGaP, the Human Cell Atlas portal, or CNGB
- controlled-access clinical cohorts
- series whose summaries never mention CD45 or enrichment
- datasets published after the search

Qian 2020 is recorded as *inferred*, not measured — the same class of
assumption that cost this project three weeks in an earlier round, so it
is flagged rather than relied upon.

## Consequence for the record

Q4 remains **unreplicated rather than refuted**, and that status is now
attributable: not to an analysis failure, but to the absence of a second
dataset with the design required. The generative fix is a study design
question — CD45⁺ enrichment, paired sampling, ≥ 20 individuals — not an
analysis question.
