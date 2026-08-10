# Deviations from the protocol, and why

Every entry is a place where the protocol as written could not be
followed literally, or where following it literally would defeat its own
stated purpose. Each records what was specified, what was done instead,
and the evidence for the change. Nothing here was changed silently.

---

## D1. Q3 absorbs every Q4 gene (specification conflict)

**Specified.** The summary table defines Q3 as "all lineages change in
the same direction". The operational rule in Phase 5.2 says
`>= 4 lineages changed and concordant`.

**Problem.** With five lineages these disagree exactly on the case the
protocol exists to find. A textbook Q4 gene — NK equivalent to zero,
CD8T/CD4T/B/Myeloid all moving together — satisfies "≥ 4 lineages
changed" and is therefore filed as Q3. On the planted-truth simulation
this misclassified **25 of 25** true Q4 genes as Q3, i.e. the primary
endpoint of the protocol had a true-positive rate of zero.

**Change.** Q3 additionally requires that NK is *not* TOST-equivalent.
A gene whose NK effect is demonstrably inside the equivalence margin is
by definition not a universal shift.

**Evidence.** `tests/test_pipeline.py::test_q3_does_not_swallow_q4` and
`::test_q3_still_fires_when_nk_also_moves`. After the fix, Q4 recovery
goes 0/25 → 24/25 with no loss of Q3 recovery (25/25).

---

## D2. The `p50` equivalence criterion costs most of the protocol's power

**Specified.** Phase 5.1: a lineage counts as unchanged when it is
TOST-equivalent **and** `|log2FC| <` that gene's 50th-percentile null
from Phase 3.

**Problem.** A gene with a genuinely zero effect falls below its own
null *median* only 50% of the time — that is what a median is. Q1
requires this simultaneously in four lineages (0.5⁴ ≈ 6% retention) and
Q2 in two (25%). The criterion is not a noise filter; it is a coin flip
applied once per lineage.

**Measured cost** (planted truth, δ = 0.5, 25 genes per class):

| equivalence bound | Q1 | Q2 | Q3 | Q4 | false Q4 |
|---|---|---|---|---|---|
| `p50` (as written) | 0/25 | 3/25 | 25/25 | 18/25 | 0 |
| `p95` | 18/25 | 20/25 | 25/25 | 24/25 | 0 |
| `p99` | 24/25 | 21/25 | 25/25 | 24/25 | 0 |

Loosening the bound cost **nothing** in false Q4 calls, because TOST is
already doing the work of controlling false equivalence at its nominal
5%.

**Change.** The bound is a parameter (`equiv_key`). The protocol default
`p50` is preserved and still runs; `p95` is the recommended setting and
what the reported analyses use. The "changed" side keeps `p99`
unaltered — an effect *exceeding* the empirical noise ceiling is a
sensible one-sided guard and costs little.

---

## D3. The Phase 3 null does not match the Phase 4 estimator

**Specified.** Phase 3 builds the null by taking adjacent-normal samples
only and splitting individuals into random halves.

**Problem.** That contrast is **unpaired**, while Phase 4 estimates a
**paired** `~ patient + condition` effect. The unpaired split-half null
therefore contains between-individual variance that the paired estimator
removes by construction, so its percentiles are systematically wider
than the noise the real estimator actually has. Feeding those wide
percentiles into Phase 5 makes "changed" harder to reach and — through
the `p50` criterion of D2 — makes "unchanged" *easier*. The bias runs
toward manufacturing exactly the Q4 calls the protocol is trying to
protect against.

**Change.** Both nulls are implemented.
`null_calibration.sign_flip_null` permutes the condition label *within*
each individual, which is the exact permutation null for the paired
estimator, and is what Phase 5 consumes by default.
`null_calibration.split_half_null` implements the protocol text and is
reported alongside for comparison.

**Trade-off, stated plainly.** The sign-flip null uses both conditions,
so a real effect is randomised into the null and slightly inflates it
(conservative). The split-half null avoids real-effect contamination but
mismatches the estimator's variance structure. The first error is
conservative for Q4; the second is anti-conservative. That asymmetry is
why the sign-flip null is the default.

---

## D4. The purity criterion is unreachable at this dataset's depth

**Specified.** Phase 2.1: the CD8 T gate should be > 90% CD3
triple-positive (`CD3D⁺ AND CD3E⁺ AND CD3G⁺` on raw counts).

**Problem.** Measured on GSE154826 library 48, median UMI per called
cell is **2,161**. At that depth, among cells that are unambiguously
CD3-protein-high, detection is CD3D 69%, CD3E 42%, CD3G 42%. Even if the
gate were perfectly pure, the expected triple-positive rate is roughly
0.69 × 0.42 × 0.42 ≈ **12%**. The observed 3% is low, but a 90%
threshold cannot be met by *any* gate on data this shallow — the
criterion silently assumes deep sequencing.

**Change.** The triple-positive rate is still computed and reported, but
as a *relative* quantity: each gate is compared against the B-cell gate,
which the protocol already designates as the ambient/soup floor
reference, and against the depth-adjusted expectation computed from that
library's own per-gene detection rates. The protocol's intent — catch T
cells hiding in the NK gate — is preserved by the discrete sub-cluster
check, which is what it says actually caught the problem last time.

---

## D5. limma-voom reimplemented in Python

**Specified.** Phase 4 recommends limma-voom or DESeq2.

**Constraint.** No R toolchain in this environment, and no network route
to install one reliably.

**Change.** `src/nkmine/de.py` reimplements TMM normalisation
(Robinson & Oshlack 2010), voom precision weights (Law et al. 2014) and
limma's empirical-Bayes variance moderation (Smyth 2004) directly in
NumPy/SciPy.

**Verification.** The batched weighted least squares is checked against
an explicit per-gene reference to ~1e-15; `fit_f_dist` recovers a known
`(d0, s0²)` prior from simulated scaled-inverse-χ² variances; effect
estimates are unbiased against planted truth (|bias| < 0.03, r > 0.99
per lineage). It has **not** been diffed against R limma itself, which
remains the right check before any result is published.

---

## D6. Cell calling uses a UMI/gene threshold, not EmptyDrops

**Specified.** Implied by Phase 6.1's use of raw droplet matrices.

**Change.** `read_10x_tar` calls cells at ≥ 500 UMI and ≥ 200 genes.
This is cruder than EmptyDrops. It is adequate for the use it is put to
here — counting cells per lineage to settle admission criterion A4, and
assigning lineage — and the ambient modelling that would justify
EmptyDrops' extra machinery is Phase 6's job, on the full droplet matrix
that is retained on disk either way.

---

## D7. The Phase 2.5 detection floor deletes the Phase 5.3 positive control

**This is the finding that mattered most on real data.** It was invisible
in simulation, because simulated genes are expressed in every lineage by
construction.

**Specified.** Phase 2.5 keeps only genes that clear a detection floor
in **all** lineages, so that "NK did not change" cannot be a restatement
of "NK does not express this gene". Phase 5.3 then requires `TOX` and
the TCR-proximal panel to land in Q4, as the positive control proving
the pipeline can detect Q4 at all.

**Problem.** These two requirements are not jointly satisfiable as
written. TCR-proximal genes are T-restricted; they are not expressed by
B or myeloid cells, so they cannot clear a detection floor imposed in
all five lineages. On GSE154826 the strict filter cut 33,723 genes to
**4,159**, and removed `TOX`, `PDCD1`, `CTLA4`, `LAG3`, `TIGIT`,
`ZAP70`, `CD28` and `TNFRSF9` — 11 of the 14 Q4 controls. Stop rule S4
then fired for the only possible reason: no positive control was left in
the panel to reach Q4.

A pipeline can therefore fail S4 while being perfectly healthy. The
first real-data run produced 45 Q4 candidates and *correctly* refused to
let them be read, for a reason that had nothing to do with the
candidates.

**Change.** `detection_filter` gains `require_nk` + `min_lineages`,
expressing the weaker condition that actually protects the inference:
the gene must be measurable **in NK** — whose non-response is the claim
being made — and in enough other lineages to supply witnesses. It need
not be measurable in lineages that are not being asked to witness
anything.

**Evidence** (`scripts/s4_diagnosis.py`, δ = 0.5, `p95`):

| configuration | genes | Q4 controls in panel | S4 |
|---|---|---|---|
| strict filter (all 5 lineages) | 4,159 | 3/14 | **STOP** |
| relaxed (NK + ≥ 2 lineages) | 4,525 | 6/14 | **PASS** |

Under the relaxed filter, `TOX` and `TIGIT` both land in Q4 with the
textbook pattern:

```
TOX     NK +0.054 (TOST-equivalent)  CD8T +0.794  CD4T +0.691  Myeloid +0.689  B +0.158 (equiv)
TIGIT   NK +0.083 (TOST-equivalent)  CD8T +1.116  CD4T +1.772  Myeloid +1.043
```

TCR-driven exhaustion genes rise in T cells and do not move in NK, which
is exactly what the protocol predicted for its own positive control.
**The pipeline has demonstrated Q4 power.**

**A prediction that was wrong, recorded because it was wrong.** Before
running this, the expectation was that the Q4 witness rule would *also*
block: witnesses are drawn from `{CD8T, B, Myeloid}`, and a TCR-driven
gene should move only in T lineages, supplying at most one witness. The
data refuted that — `TOX` and `TIGIT` also move in **myeloid** cells, so
they reach two witnesses without needing CD4T. Adding CD4T to the
witness set changes the Q4 count from 58 to 64 but was never the
blocker. The witness set is left at the protocol's specification.

(That myeloid cells move for `TOX` and `TIGIT` at all is worth a second
look later — it is either real biology or ambient contamination from the
abundant T compartment, which is precisely the question Phase 6.1 exists
to answer, and it has not been run.)

---

## D8. Stop rule S2 fires at the protocol's stated threshold

**Specified.** A1 requires ≥ 8 paired individuals; A4 requires ≥ 30 NK
cells per (individual × condition); S2 says stop if fewer than **two**
datasets clear A1–A4, because Phase 7 makes cross-dataset replication a
hard gate.

**Measured** (`results/A4_cross_dataset_scan.csv`; "balanced" = all five
lineages clear the threshold in both conditions, which is what Phase 2.4
actually requires):

| dataset | paired patients | NK ≥ 30 | balanced ≥ 30 | balanced ≥ 20 |
|---|---|---|---|---|
| GSE154826 (Leader, CD45⁺) | 29 | 27 | **27** | **28** |
| GSE131907 (Kim, unsorted) | 10 | 8 | 6 | **8** |
| GSE178341 (Pelka, unsorted) | 36 | 5 | 5 | 7 |

**At the protocol's own threshold, exactly one dataset qualifies, so S2
fires.** At a relaxed 20-cell threshold two qualify — but GSE131907
lands at exactly n = 8, the bare A1 minimum.

**The mechanism, and it generalises.** Pelka has the *most* paired
patients of the three (36) and is still the worst: it is unsorted, so NK
is 3,924 of 370,115 cells (1.1%), and NK is the limiting lineage in 53
of its 72 groups. **CD45⁺ enrichment, not patient count, is what
delivers A4 for a rare lineage.** Any future dataset search should
filter on enrichment protocol first and cohort size second — the
opposite of the intuitive ordering, and the reason the protocol's own
nomination of GSE178341 as a discovery cohort does not survive contact
with the data.

**Consequence for Phase 7.** Replication is possible only against
GSE131907, at n = 8 and a relaxed NK threshold. That is enough to
attempt, but not enough to interpret a *negative* result: a Q4 gene
failing to replicate there would be ambiguous between a false original
finding and insufficient power. Any replication analysis must report the
Q1/Q3 baseline replication rates alongside Q4 (which the protocol
already requires via S6) precisely so that this ambiguity is visible.

---

## D9. Categorical replication is the wrong instrument at n = 8

**Specified.** Phase 7 defines replication as: same quadrant, same sign,
`|log2FC|` within 2×.

**Problem.** A quadrant call is a conjunction of several threshold
crossings, at least one of which is an *equivalence* claim. In
GSE131907 the median moderated SE for CD8 T is 0.308, so the 90% CI
half-width is 1.645 × 0.308 = 0.507 — wider than the equivalence margin
δ = 0.5. A CD8T equivalence call is therefore **arithmetically
impossible** in that cohort, and any quadrant requiring one cannot
replicate no matter how real the underlying biology.

Measured: Q1 0/10, Q3 0/53, Q4 1/64 replicate categorically. Q3 is the
easiest category — large concordant effects in every lineage — and it
replicates at zero.

**Change.** Categorical replication is still computed and reported (it
is what the protocol asks for, and S5 fires on it). Alongside it, a
continuous criterion is reported that does not depend on any threshold:
the magnitude of the NK effect for candidate genes, compared against Q3
genes and against background, in the replication cohort. On that
criterion the Q4 set replicates at p = 1.3 × 10⁻⁶.

**Why this is not a lowered bar.** The continuous test is *more*
demanding in the way that matters: it compares Q4 against Q3, and Q3
genes have larger witness effects, so a gene set that merely had
"big effects in other lineages" would fail it. It also survives the
power confounder — Q4 candidates' NK SE in the replication cohort is
smaller than Q3's, so their small NK effects are not a precision
artefact.

**What it does not license.** One cohort at n = 8 is not the two
independent replications the protocol demands, and a continuous signal
across a 64-gene set says nothing about which individual genes are real.
Phase 8 has not been run.

---

## D10. The Phase 2.6 guard was never armed, and arming it retracts the TOX result

Found by an independent adversarial audit of this repository, not by the
author. Recorded in full because the correction reverses a headline claim.

**The defect.** `pseudobulk.saturation_flags()` was implemented,
documented, and unit-tested — and never called. Every `classify_table`
call site omitted `saturated=`, so the argument defaulted to all-False
and **no gene in any published result had ever been tested for the
ceiling/floor conflation**. `is_equivalent`'s docstring asserted a
protection that was not in force.

**What arming it changed.** `TOX` is detected in **4.3%** of NK cells in
tumour and **4.9%** in adjacent normal. It is floor-saturated: there is
no room for it to move in NK. With the guard armed it is barred from
equivalence, lands in `unclassified`, and **stop rule S4 fires again**.

**The claim being retracted.** An earlier revision of the README stated
that the pipeline "has demonstrated Q4 power" because `TOX` and `TIGIT`
reached Q4. That rested on a guard that was not running. `TOX` does not
reach Q4 once Phase 2.6 is enforced. The correct status is that **S4
fires and the Q4 candidate list is not cleared for reading.**

**The deeper problem, which is the protocol's and not the code's.** The
Q4 positive control is self-defeating. TCR-proximal genes are
T-restricted, so the Phase 2.5 detection floor removes them (D7); and
they are barely expressed in NK, so the Phase 2.6 saturation guard
excludes them. Both guards are correct. The conflict is that "NK does
not respond to a TCR-driven programme" is not evidence of NK
*resistance* — NK cells have no TCR, so it is mundane explanation **B3
(receptor not expressed)**, which Phase 8 exists to rule out. The
protocol nominates as its positive control an instance of the artefact
it is designed to reject.

**Consequence.** There is currently **no valid Q4 positive control**. A
usable one must be a gene that NK cells demonstrably express and could
in principle regulate, whose driver is nonetheless shared with the
witness lineages. Until such a control exists and passes, S4 cannot be
satisfied and no Q4 list from this pipeline should be read as biology.

## D11. Other defects found in the same audit

- **Q3 guarded on the noise-screened `equiv["NK"]` rather than raw
  TOST**, so a gene demonstrably TOST-equivalent in NK — the exact Q4
  signature — fell into Q3 whenever the Phase 3 screen happened to fail.
  Observed on real data at the protocol's own `p50` default. Now guarded
  on the raw TOST result.
- **Q3 fired when NK was merely *indeterminate*** (neither changed nor
  equivalent): 10 of 53 Q3 genes. Four positive demonstrations plus one
  non-result is the "p > 0.05 therefore unchanged" inference this
  protocol exists to refuse, applied at the Q3 end. Now labelled
  `Q3_NK_indeterminate` and kept separate.
- **"Equivalent" was reported as "equivalent to zero".** TOST
  establishes only that an effect is smaller than δ; 8 of 64 Q4 calls
  had an NK 90% CI *excluding* zero, and one (`RAB11FIP1`) moved
  opposite to its witnesses. Split into `Q4` (indistinguishable from no
  change) and `Q4_attenuated` (moves, but by less than δ).
- **voom's mean-variance trend used `log2(mean(lib_size))` where limma
  uses `mean(log2(lib_size))`.** The two differ by Jensen's inequality
  and shift the trend x-axis, feeding into every precision weight.
  Corrected.
- **CI90 and the effective null thresholds were not written to the
  output**, so an equivalence call could not be audited from the
  artefact. Now emitted, along with `equiv_key`/`change_key`.
- **A docstring pointed at `scripts/sensitivity_equiv_bound.py`, which
  does not exist.** Corrected to point at D2 and `s4_diagnosis.py`.

---

## D12. A negative control, since the positive control is unavailable

S4 fires because the protocol's Q4 positive control cannot work (D10),
which leaves the pipeline's sensitivity unproven. Its *specificity* is
still testable with the data in hand, and that is worth having on
record.

`scripts/q4_null_rate.py` flips the tumour/normal label **within each
individual**, for all five lineages together, and re-runs the whole
Phase 4–5 stack. Pairing, library sizes, lineage composition, gene
abundance and the ambient structure are all preserved; only the
condition assignment is randomised. Any quadrant call surviving that is
manufactured.

| quadrant | observed | null (12 permutations) |
|---|---|---|
| Q4 | 40 | **0.0** [0, 0] |
| Q3 | 50 | **0.0** [0, 0] |
| Q1 | 8 | 0.2 [0, 2] |

0 of 12 permutations produced a single Q4 call. The pipeline does not
generate Q4 calls out of noise.

**What this does and does not establish.** It bounds the false-positive
side only. A pipeline that reported nothing under the null *and* nothing
under real signal would look identical here, which is exactly why a
positive control is not optional and why S4 still stands. Combined with
the independent continuous replication (D9), the aggregate evidence is
that the Q4 set is not noise — but neither result validates any
individual gene, and Phase 8 has not been run.

**An instability worth recording.** The observed Q4 count is 40 here and
46 in `s4_diagnosis.py`. The only difference is the number of Phase 3
calibration permutations (60 vs 100), which moves the per-gene noise
percentiles and therefore the classification boundary. A ±15% swing in
the headline count from a calibration parameter is a reminder that the
Q4 *count* is soft; the pattern-level results are what carry weight.
