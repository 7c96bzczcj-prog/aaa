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
