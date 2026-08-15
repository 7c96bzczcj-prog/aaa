# The archived OXPHOS "11/11 donors, p = 0.001" — it does not stand

**Verdict: the claim as archived is void.** It fails on arithmetic, before any
re-analysis, and the failure is not recoverable by recomputation.

**What this rests on.** Two inputs, no data access:

1. **Given** (spec author): that dataset has **at most 7 donors**, and **6** if
   restricted to decidua.
2. **Archived**: the claim reads **"11/11 donors, p = 0.001"**.

Everything below follows from those two.

---

## 1. The arithmetic contradiction

**11 > 7.** If the dataset has at most 7 donors, the 11 units in "11/11 donors"
**cannot be donors.** That is decisive on its own and needs no test theory: the
claim's stated denominator exceeds the number of independent units in existence.

Whatever the 11 are, they are one of:

| decomposition | what the 11 really are |
|---|---|
| 6 decidua donors × 2 conditions = 12, minus one dropped | **donor × condition**, not donors |
| ≤ 7 donors contributing 11 libraries | **libraries**, not donors |
| 6 decidua + 5 other-compartment | **cross-compartment samples** pooled as if independent |

All three are the same error class: **the unit of replication is not the unit that
was counted.** This is `n = wells` with a different label.

---

## 2. The p-value is a fingerprint of the count, not of the effect

`p = 0.001` is almost certainly **2 / 2¹¹ = 0.0009766**, the two-sided sign /
binomial test on 11 units all pointing the same way.

**k = 11 is the *unique* k for which 2/2^k rounds to 0.001:**

| k | 2/2^k | rounds to |
|---|---|---|
| 10 | 0.0019531 | 0.002 |
| **11** | **0.0009766** | **0.001** |
| 12 | 0.0004883 | 0.000 |

This matters more than it looks. **The p-value is a function of the count itself**,
not of effect sizes. So correcting the unit does not require re-running anything —
it rewrites p directly, and the ceiling can be computed now:

| n used | best achievable two-sided p (all same direction) |
|---|---|
| 11 (as archived) | **0.001** |
| **7 (dataset maximum)** | **0.016** |
| **6 (decidua only)** | **0.031** |

**Even in the most favourable case — every real donor pointing the same way —
p = 0.001 is unreachable.** The ceiling at the true n is 0.016, and 0.031 if the
claim is about decidua. Neither is 0.001, and 0.031 is marginal by any convention.

---

## 3. The one reading under which the ceiling would not apply, and why it does not rescue the claim

The ceiling above is exact for **sign and signed-rank tests** (Wilcoxon signed-rank
also bottoms out at 2/2⁷ = 0.016 at n = 7). A **paired t-test** could in principle
reach p < 0.001 at n = 7 given a large, consistent effect — so if the archived p
came from a t-test, the ceiling argument is not airtight.

**But the phrasing rules that out.** "11/11 donors" is a *count-of-directions*
statement — it is the sign-test framing by construction. A t-test result would not
be reported as "11/11".

**And even under the t-test reading the claim still fails**, because §1 is
independent of test choice: the units still are not donors, so whatever test was
run was run on pseudoreplicated units and must be recomputed at the donor level
regardless. Reading A kills it arithmetically; reading B kills it by
pseudoreplication. There is no reading in which it survives untouched.

---

## 4. What to check in the workspace, in this order

1. **`donor` column unique values.** Confirms which decomposition in §1 applies.
   Expect ≤ 7 distinct real donors and 11 distinct labels.
2. **Grep archived p-values equal to 2/2^k** — this is the signature, and it will
   cluster in the analyses sharing the same unit error:
   ```
   grep -rnE "0\.0009766|0\.00098|0\.001953|0\.0039|0\.0078|0\.0156|0\.0312" .
   grep -rnE "\b(11|12|10)/(11|12|10)\b" .
   ```
   Any p in `{0.031, 0.016, 0.008, 0.004, 0.002, 0.001}` reported alongside an
   "all k of k" statement is the same construction.
3. **Whether the 11 mixes compartments.** If it is 6 decidua + 5 other, the error
   is worse than pseudoreplication — it is a between-compartment comparison
   reported as within.

---

## 5. Consequence for the main line

The spec author's own framing was that this is **one of three independent supports
for the main line direction**. On the above it is not a support at all:

- it cannot be repaired by recomputation — the ceiling at the true n is 0.016/0.031;
- **it was never independent** if the other two supports are computed on the same
  donor set with the same unit convention. That is the second thing to check, and
  it is more important than this one result: a shared unit error does not produce
  three independent supports, it produces one error reported three times.

**So the proposal section resting on it has to be rewritten**, and the rewrite
should not simply substitute p = 0.016. At n = 6–7 with a directional consistency
argument, the honest statement is *"consistent in 6/6 (or 7/7) donors, p = 0.031
(0.016) by sign test"* — reported as a small, directionally consistent
observation, which is what it is.

**What this does not support.** It does **not** show the underlying biology is
wrong. Every real donor may well point the same way, and at n = 6 that is still
worth reporting. What is void is the *strength* claimed for it, and its status as
an independent pillar.

---

## 6. Confidence, stated plainly

- §1 (11 cannot be donors): **certain**, given the stated n ≤ 7.
- §2 (p = 2/2¹¹ sign test): **strong inference** — k = 11 is the unique k rounding
  to 0.001, and the "11/11" phrasing is the sign-test form. Not verified against
  the code.
- §5 (the three supports may not be independent): **hypothesis to check**, not a
  finding.

The one number that settles §2 and §4 at once is the `donor` column's unique
values.
