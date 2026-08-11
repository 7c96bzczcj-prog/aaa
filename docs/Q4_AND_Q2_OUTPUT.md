# The output nobody read: the Q4 gene list, and why Q2 is empty

The instrument was audited to exhaustion — sensitivity measured, null
rate measured, construction compared, retraction corrected, 42 agents of
adversarial review — and **the output was never printed once**. The 46
Q4 genes appear in no report in this repository's history. Neither does
the fact that **Q2 = 0**.

Both are below. Full list: `results/Q4_gene_list.csv`.

---

## 1. The 46 Q4 genes

Ranked by mean witness effect. NK is flat in all of them by
construction; the question is what the other lineages are doing.

| rank | gene | NK | CD8T | CD4T | B | Myeloid |
|---|---|---|---|---|---|---|
| 1 | **TIGIT** | +0.08 | **+1.12** | **+1.77** | +0.52 | **+1.04** |
| 2 | **CD38** | +0.01 | **+1.41** | **+1.08** | **+1.52** | +0.17 |
| 3 | SEPT11 | −0.16 | −0.84 | −1.08 | −0.66 | −1.07 |
| 4 | MYADM | −0.06 | −0.84 | −0.75 | −1.17 | −0.85 |
| 5 | S100A10 | +0.16 | −0.77 | −0.71 | −1.23 | −0.59 |
| 6 | PDLIM1 | −0.13 | −0.39 | −0.69 | −0.61 | −1.53 |
| 7 | S100A6 | +0.10 | −0.61 | −0.61 | −1.33 | −0.59 |
| … | (40 more) | | | | | |

**Read with discipline — this is eyeballing, not enrichment.** No
GO/pathway database is available offline, so no formal enrichment test
was run. What follows is pattern-recognition on gene symbols and should
be treated as a hypothesis about the list, not a result.

Three groupings are visible:

**(a) Inhibitory / activation-exhaustion, at the very top.**
`TIGIT`, `CD38`, and `PRDM1` (Blimp-1) rise sharply in T and myeloid
cells while NK does not move. This answers the first branch of the
question directly: **the strongest Q4 genes are exhaustion-programme
genes, and the reading is "NK does not upregulate the exhaustion
programme that the other lymphocytes do."** Note this is the same shape
as `TOX` and `TIGIT` in the positive-control work — and `TIGIT`'s
appearance here is not independent evidence, since it was already the
one control gene that survived.

**(b) A cytoskeleton / membrane-scaffold block, going DOWN in the
witnesses while NK holds steady.** `SEPT11`, `PDLIM1`, `TAGLN2`,
`ANXA2`, `S100A10`, `S100A6`, `LGALS1`, `MYADM`, `TUBA4A`, `PECAM1`,
`STX11`. If real, this reads as *the other lineages dismantle
cytoskeletal/trafficking machinery on tumour entry and NK does not* —
which is the **opposite** direction from "NK acquires a brake".

**(c) Interferon-stimulated genes**: `GBP1`, `GBP2`, `ISG20`, `EPSTI1`,
and `IFNGR1` (down in witnesses).

Direction split: 19 genes up in witnesses, 27 down.

`SELL` (CD62L) is present — up in CD8T/CD4T/B, flat in NK — which bears
directly on mundane explanation **B2 (NK turnover / recent immigration)**
and should be checked there rather than read as a finding.

---

## 2. Q2 = 0, and the reason is not what was expected

Q2 was never reported. It is **zero** in every configuration run except
δ = 0.25 / p95, where it is 2 (`CORO1B`, `SPSB3` — both with tiny
effects, ≈ +0.3).

### The funnel (δ = 0.5, p95)

| stage | genes |
|---|---|
| total | 4,525 |
| NK changed | 345 |
| + CD8T also changed | 66 |
| + same direction | 63 |
| + B equivalent | 4 |
| + Myeloid equivalent = **Q2** | **0** |

### The expectation that Q2 is better-powered than Q4 is wrong

The argument was: Q4 puts the hard *equivalence* claim on NK (rare, worst
power), while Q2 puts it on B and myeloid (abundant, good power). That
reasoning is sound but **it is not the binding constraint.** Among the 63
candidates:

| lineage | genuine effect (\|log2FC\| ≥ δ) | small effect, failed equivalence | equivalent |
|---|---|---|---|
| B | **49 / 63** | 10 | 4 |
| Myeloid | **39 / 63** | 16 | 8 |

**B and myeloid are not failing the equivalence test for lack of power.
They are failing it because they genuinely move.** Nearly every gene
that NK and CD8 T share is also moving in B and myeloid — i.e. it is
Q3-shaped, not Q2-shaped.

---

## 3. Why they move: ambient, proved by immunoglobulin

The candidate pool contains `IGKC +3.72`, `IGHG1 +4.02`, `IGLC2 +4.37`,
`MZB1 +2.92` — **in NK cells**.

| gene | NK | CD8T | CD4T | B | Myeloid |
|---|---|---|---|---|---|
| IGKC | +3.72 | +3.24 | +2.84 | +2.27 | +3.23 |
| IGHG1 | +4.02 | +3.04 | +3.01 | +3.31 | +3.99 |
| IGLC3 | +4.68 | +3.69 | +3.41 | +3.50 | +3.49 |
| MZB1 | +2.92 | +2.48 | +2.12 | +2.07 | +2.49 |

NK cells and myeloid cells do not produce immunoglobulin. A +3 to +4.7
log2FC in those lineages can only be soup — plasma cells expand in
tumour, and their transcripts land in every droplet.

**This is a better ambient positive control than anything the protocol
specified**, because its true within-lineage value is known to be exactly
zero, so its magnitude *is* the ambient contribution on an absolute
scale.

### Calibrated against the shared-emulsion design

Median |log2FC| across the 8 immunoglobulin/plasma-cell genes:

| | in NK | in Myeloid |
|---|---|---|
| shared soup (n = 7) | **1.68** | **0.21** |
| separate libraries (n = 19) | **5.67** | **4.68** |

Ambient is **3.4× smaller in NK and 22× smaller in myeloid** when the two
conditions share one droplet emulsion. This is a cleaner and more
interpretable measurement than the correlation-based estimate reported
earlier (0.291 vs 0.490), because it is in log2FC units against a known
zero.

*(Open question, not explained away: why NK's residual is 1.68 rather
than ~0 in shared soup. Candidates are n = 7 noise, HTO demultiplexing
error, or genuine B/plasma contamination of the NK gate — the sub-cluster
purity check that Phase 2.1 asks for has still not been run.)*

### The structural consequence: ambient specifically destroys Q2

Ambient pushes **every** lineage off zero, in the same direction, at once.
Q2 is the only quadrant that requires **two abundant lineages to be
TOST-equivalent**. So ambient attacks Q2 precisely where Q2 is defined.

Suggestive but not conclusive support, using CI-based criteria in the two
designs:

| | candidates (NK+CD8T concordant) | Q2 | conversion |
|---|---|---|---|
| shared soup (n = 7) | 29 | **2** | 6.9% |
| separate libraries (n = 19) | 105 | 1 | 1.0% |

A ~7× higher conversion rate where ambient cancels, despite a third of
the individuals. **n = 7, so this is a direction, not a number.**

### The protocol's own Q2 positive control fails the same way

Phase 5.3 nominates the cytotoxic programme as the Q2 control. It does
not land there:

| gene | quadrant | NK | CD8T | **B** |
|---|---|---|---|---|
| GZMB | unclassified | −0.66 | +0.16 | **−0.98 (not equivalent)** |
| PRF1 | unclassified | −0.87 | −0.47 | **−1.41 (not equivalent)** |
| GZMA | unclassified | −0.34 | +0.19 | **−0.68 (not equivalent)** |
| CTSW | **Q3** | −0.65 | −0.29 | **−0.96 (not equivalent)** |
| GZMH | **Q3** | −0.88 | −0.55 | **−0.93 (not equivalent)** |

In every case the blocker is **B cells failing equivalence for perforin
and granzymes** — genes B cells do not meaningfully express. That is
ambient again.

**So Q2 is in the same position Q4 was: its positive control is
destroyed before the quadrant can be evaluated.** Same diagnosis,
different quadrant, different mechanism — Q4's control was killed by the
detection floor and the saturation guard (D7/D10), Q2's is killed by
ambient. **Q2 is currently untestable in separate-library data, not
empty of biology.**

---

## 4. The polarization-machinery panel: no prioritisation obtained

Checked as requested, with its stated status: prioritisation only, never
falsification.

20 of 24 panel genes cleared the detection floor. **None is in Q2** —
which follows, since Q2 is empty. Two are in Q4:

| gene | quadrant | NK | CD8T | B | Myeloid |
|---|---|---|---|---|---|
| STX11 | **Q4** | −0.08 | −0.07 | −0.73 | −0.78 |
| TUBA4A | **Q4** | +0.23 | −0.19 | +0.58 | +0.58 |

`STX11` (syntaxin-11, the FHL4 gene, required for lytic granule fusion)
is the more interesting of the two — but note its witness pattern is
B and myeloid moving while CD8T does **not** (−0.07), which is not the
lymphocyte-shared shape the hypothesis predicts, and is the pattern
ambient produces. **One gene, from an unreplicated list, with an
ambient-suspicious signature. It does not prioritise anything.**

The intended interface — "if the brake is lymphocyte-shared it lands in
Q2" — cannot be used, because Q2 cannot currently be evaluated.

---

## 5. What this changes

- **Q4's top genes are exhaustion genes** (`TIGIT`, `CD38`, `PRDM1`),
  supporting the "NK does not enter the exhaustion programme" reading
  over the "tissue-adaptation" reading — with the caveat that this is
  eyeballed, unreplicated, and `TIGIT` is not independent of the
  control work.
- **Q2 is not empty of biology; it is untestable under ambient**, and
  that is a specific, mechanistic, fixable diagnosis rather than a null
  result.
- **Immunoglobulin genes are a free, unambiguous, absolute-scale ambient
  control** for any droplet-based multi-lineage experiment. This should
  be a standard QC panel and, as far as this project's searches went, is
  not used as one.
- The correct next step for Q2 is **ambient correction** (CellBender on
  the raw droplet matrices, which are already on disk), then re-running
  the quadrant assignment — not more discovery.
