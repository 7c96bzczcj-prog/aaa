# PREREGISTRATION — X1 Tier 1 (human clone-sharing re-analysis of GSE302113)

**Committed before any Tier-1 statistic was computed.** Read-only after first
commit; amendment requires a version bump recorded in `docs/CHANGELOG.md` with
reason and timestamp, exactly as for the screening preregistration.

**Scope.** Tier 1 addresses **clause (a) of X1 only — provenance**. It cannot
address clause (b), the temporal clause, because non-involved lung tissue (NILT)
is contemporaneous with the tumour, not pre-tumour, and no human sample exists of
an organ from before its own tumour arose (`docs/X1_DESIGN.md` §4).

**Why clause (a) alone is worth running.** If intratumoural NK are clonally
related to NILT-resident NK and *not* to blood NK, then "TGF-β converted
infiltrating NK" fails as an account of that population **regardless of when those
cells arrived**. The cheap tier carries most of the value; Tier 2 serves only the
temporal half.

---

## 1. The statistic is three-way, not pairwise (fixed)

Pairwise tumour∩NILT sharing is **insufficient**: tumour and NILT NK could both
have arrived recently from blood and both been converted, which produces tumour–
NILT sharing without any residency. And because every NK cell in one person shares
an HSPC ancestor, bare sharing is trivially positive at sufficient depth under
every hypothesis.

**The discriminator is whether the clones shared between tumour and NILT are
under-represented in blood.** Rationale: a resident pool should carry a stronger
clonal structure than blood — local expansion, slower turnover.

Fixed predictions, written before computation:

| hypothesis | prediction |
|---|---|
| **resident origin** (X1) | clones shared between tumour and NILT are **rare or absent in blood** |
| **blood origin** (recruitment) | those same clones appear in blood at **comparable frequency** |

### 1.1 Primary statistic

Per donor `d`, per clone `c`, let `t_c`, `n_c`, `b_c` be the number of **NK** cells
of clone `c` recovered in tumour, NILT and blood.

Define the **tumour–NILT shared set** `S_d = { c : t_c > 0 AND n_c > 0 }`.

**Primary readout** — blood representation of `S_d`:

```
f_blood(S_d) = ( # clones in S_d with b_c > 0 ) / |S_d|
```

compared against a **size-matched null**: clones not in `S_d`, sampled to match
the distribution of total tissue size `(t_c + n_c)`, giving `f_blood(null_d)`.

Report `Δ_d = f_blood(S_d) − f_blood(null_d)` with a permutation CI
(**1000 permutations**, fixed).

- `Δ_d` significantly **negative** → tumour–NILT shared clones are
  under-represented in blood → **supports resident origin**.
- `Δ_d ≈ 0` or **positive** → **supports blood origin**.

Size-matching is mandatory because larger clones are more likely to be detected in
any compartment, which would otherwise manufacture the effect.

### 1.2 Secondary statistic (also fixed in advance)

All three pairwise sharing fractions — `s(T,N)`, `s(T,B)`, `s(N,B)` — computed with
**the source paper's own definition** (fraction of same-clone cell pairs between
two sets, normalised by all within-donor cross-set pairs), so the numbers are
comparable to that paper's published myeloid results (G6: a summary computed on
the same basis as the rule consuming it).

### 1.3 Concurrent controls (fixed)

The identical statistic is computed for **myeloid** and **CD8 T / T** cells in the
same donors. These calibrate the statistic on lineages whose behaviour the source
paper already characterised: myeloid clones are reported as clonally linked to
circulating monocytes, so myeloid should sit toward the blood-origin end. A NK
result that does not differ from the myeloid control is not evidence of residency.

---

## 2. Sampling correction (fixed)

Clone-sharing statistics are strongly sensitive to the number of cells sampled,
and the tumour compartment carries only **281–509 NK cells** in the passing donors
while blood carries up to 7,713.

Fixed procedure:

1. **Downsample each compartment to equal `n` per donor**, where
   `n = min(NK cells in tumour, NILT, blood)` for that donor.
2. Repeat **200 times** (fixed); report the **median and 2.5–97.5 percentile
   interval** across resamples.
3. **Report per donor. Never pool donors into a single number.** Cross-donor cell
   pairs are not valid comparisons (the source pipeline excludes them for the same
   reason), and pooling would let one deep-sequenced donor set the result.
4. The undownsampled values are reported alongside as a **secondary** readout, and
   any qualitative disagreement between downsampled and raw is reported, not
   resolved silently.

---

## 3. "Same clone" thresholds (fixed — this is an admission threshold)

Treated as the same class of object as the 20× chrM floor: **fixed before looking
at the data, and not adjusted afterwards in either direction.**

**Variant admission** (primary):
- mgatk defaults as already preregistered in DEC-04: **strand correlation ≥ 0.65,
  variance–mean ratio ≥ 0.01, confidently detected in ≥ 5 cells**.
- **Union across the donor's libraries** — a variant passing in any library of a
  donor is admitted for all libraries of that donor. This is both the source
  pipeline's procedure and the only way to compare compartments at all.
- **chrM:307–314 excluded** (homopolymeric region; the source authors exclude it
  because it creates spurious clonal links).
- Variants detected in **> 20% of a donor's cells** removed as likely homoplasmic
  or technical.

**Cell-level clone assignment**:
- Per-cell heteroplasmy **binarised at 0.07**, the source pipeline's benchmarked
  cutoff, on the stated rationale that exact heteroplasmy is unreliable given
  stochastic mtDNA partitioning at division.
- A cell is assigned to a clone only if it carries **≥ 1 admitted variant** after
  binarisation; cells carrying none are **excluded and counted**, not silently
  dropped.

**Sensitivity analysis, also fixed in advance** (run regardless of the primary
result): repeat at the source authors' looser variant threshold (**≥ 3 cells**
instead of ≥ 5) and at binarisation **0.05 and 0.10**. If the sign of `Δ_d` is not
stable across these, the result is reported as **unstable**, not as the primary
value.

---

## 4. Order of operations (fixed)

1. **Extract the published number first.** The source paper's Figure 2F plots
   tumour-cell-type × NILT-cell-type clone sharing, and NK/ILC is an annotated
   cell type in that atlas, so the tumour-NK × NILT-NK cell **already exists in a
   published figure and was never quantified, tested or interpreted in the text**.
   Recover it from the figure or supplementary data before running anything.
2. **Then run the analysis above.**
3. **If the two disagree, that disagreement is the first thing to explain** — it
   is reported before any biological reading, and no biological conclusion is
   drawn until it is resolved or explicitly recorded as unresolved.

---

## 5. Admissibility and stopping (fixed)

- Donors admitted: those passing the screening power check — **SU-L-001, SU-L-002,
  SU-L-004, SU-L-005**. SU-L-003 is excluded (no NK detectable in tumour at
  cluster resolution).
- **If fewer than 3 donors yield ≥ 1 clone in `S_d`, the analysis fails power and
  is reported as failing power.** No threshold is lowered to produce a runnable
  result.
- A null result is reported as a null result. Given the known ~15–26% mtDNA clone
  assignment rate in solid tissue, **a null is expected to be uninformative rather
  than evidence against X1**, and will be reported that way.

## 6. Known limits carried in from the screen

1. NK typing is **cluster-level and validated only on blood** (CD56⁺ 74.6% vs
   CD56⁻ 5.5%). It is a power estimate, not a tissue-validated classifier, and the
   ILC1/trNK boundary is exactly where it is weakest — which is where this claim
   lives. Any Tier-1 result inherits that uncertainty.
2. NILT is contemporaneous, not pre-tumour (§0).
3. n = 4 donors.
