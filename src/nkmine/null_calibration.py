"""Per-lineage, per-gene noise calibration (Protocol Phase 3).

The protocol's motivating observation: replicate noise is both lineage-
specific and gene-specific (one lineage's detection rate can swing 15
percentage points between two same-condition replicates while another
gene in another lineage swings 0.2).  A single global threshold is
therefore meaningless, and Phase 5 consumes per-(lineage, gene)
percentiles instead.

Two nulls are provided, and they are not interchangeable:

`sign_flip_null` (default)
    Permutes the condition label *within* each individual, i.e. flips
    the sign of each within-pair difference.  This is the exact
    permutation null for the paired ``~ patient + condition`` estimator
    that Phase 4 actually uses, so its spread is directly comparable to
    the real effect estimates.

`split_half_null`
    Literally what the protocol text specifies: take one condition only
    (adjacent-normal), split individuals into halves, contrast them.
    This is an *unpaired* contrast, so it carries between-individual
    variance that the paired Phase 4 estimator removes.  Its percentiles
    are consequently wider than the paired estimator's noise.  Using
    them as the Phase 5 floors makes "changed" harder to reach and, via
    the p50 criterion, makes "equivalent" *easier* -- which biases
    toward exactly the Q4 calls the protocol is trying to protect.

Both are computed and reported; the sign-flip percentiles are the ones
wired into Phase 5 by default.  See docs/DEVIATIONS.md.
"""

from __future__ import annotations

import numpy as np

from .de import _wls, fit_f_dist, voom


def _design(patient, condition, case_label):
    pats = sorted(set(np.asarray(patient).tolist()))
    cols = [np.ones(len(condition))]
    for p in pats[1:]:
        cols.append((np.asarray(patient) == p).astype(float))
    cols.append((np.asarray(condition) == case_label).astype(float))
    return np.column_stack(cols)


def _fit_lfc(log_cpm, weights, X, coef_index):
    """One WLS + eBayes pass, returning just the coefficient of interest."""
    beta, sigma2, df_resid, stdev_unscaled, _ = _wls(log_cpm, X, weights)
    d0, s02 = fit_f_dist(sigma2, df_resid)
    if np.isinf(d0):
        s2_post = np.full_like(sigma2, s02)
    else:
        s2_post = (d0 * s02 + df_resid * sigma2) / (d0 + df_resid)
    return beta[:, coef_index], stdev_unscaled[:, coef_index] * np.sqrt(s2_post)


def sign_flip_null(
    counts: np.ndarray,
    patient: np.ndarray,
    condition: np.ndarray,
    case_label: str,
    n_perm: int = 200,
    rng: np.random.Generator | None = None,
):
    """Exact paired permutation null: randomise condition within individual.

    voom weights are estimated once under the real design and reused
    across permutations (standard practice -- the mean-variance trend is
    a property of the count data, not of the label assignment, and
    re-estimating it 200x is the dominant cost).
    """
    rng = rng or np.random.default_rng(0)
    counts = np.asarray(counts, dtype=float)
    patient = np.asarray(patient)
    condition = np.asarray(condition)

    X_real = _design(patient, condition, case_label)
    log_cpm, weights, _ = voom(counts, X_real)
    coef = X_real.shape[1] - 1

    pats = sorted(set(patient.tolist()))
    null = np.empty((n_perm, counts.shape[0]))
    for k in range(n_perm):
        perm_cond = condition.copy()
        for p in pats:
            m = np.flatnonzero(patient == p)
            if len(m) == 2 and rng.random() < 0.5:
                perm_cond[m] = perm_cond[m[::-1]]  # swap the two labels
        X = _design(patient, perm_cond, case_label)
        lfc, _ = _fit_lfc(log_cpm, weights, X, coef)
        null[k] = lfc
    return null


def split_half_null(
    counts: np.ndarray,
    patient: np.ndarray,
    n_perm: int = 200,
    rng: np.random.Generator | None = None,
):
    """Protocol-literal null: split individuals of ONE condition in half.

    `counts` must already be restricted to a single condition.
    """
    rng = rng or np.random.default_rng(0)
    counts = np.asarray(counts, dtype=float)
    patient = np.asarray(patient)
    pats = np.array(sorted(set(patient.tolist())))
    if len(pats) < 4:
        raise ValueError("need >= 4 individuals to form split halves")

    X_dummy = np.column_stack([np.ones(len(patient)), rng.integers(0, 2, len(patient)).astype(float)])
    log_cpm, weights, _ = voom(counts, X_dummy)

    null = np.empty((n_perm, counts.shape[0]))
    for k in range(n_perm):
        shuffled = rng.permutation(pats)
        group_a = set(shuffled[: len(pats) // 2].tolist())
        grp = np.array([1.0 if p in group_a else 0.0 for p in patient])
        X = np.column_stack([np.ones(len(patient)), grp])
        lfc, _ = _fit_lfc(log_cpm, weights, X, 1)
        null[k] = lfc
    return null


def null_percentiles(null: np.ndarray) -> dict:
    """Collapse a permutation matrix into the percentiles Phase 5 needs.

    Percentiles are taken on |log2FC|, because both Phase 5 criteria are
    stated on the magnitude of the effect.
    """
    a = np.abs(null)
    return {
        "p50": np.percentile(a, 50, axis=0),
        "p95": np.percentile(a, 95, axis=0),
        "p99": np.percentile(a, 99, axis=0),
        "max": a.max(axis=0),
    }


def calibrate(
    ps,
    lineages=("NK", "CD8T", "CD4T", "B", "Myeloid"),
    case_label: str = "Tumor",
    n_perm: int = 200,
    rng: np.random.Generator | None = None,
) -> dict:
    """Run the sign-flip null for every lineage of a PseudobulkSet."""
    rng = rng or np.random.default_rng(0)
    out = {}
    for l in lineages:
        sub = ps.subset_lineage(l)
        null = sign_flip_null(
            sub.counts, sub.patient, sub.condition, case_label,
            n_perm=n_perm, rng=rng,
        )
        out[l] = null_percentiles(null)
        out[l]["_null"] = null
    return out


def cross_lineage_noise_summary(null_stats: dict) -> "object":
    """Phase 3 use #2: show that noise levels differ *between* lineages.

    Justifies why Phase 5 uses per-lineage floors rather than one global
    threshold.
    """
    import pandas as pd

    rows = []
    for l, s in null_stats.items():
        rows.append(
            {
                "lineage": l,
                "median_null_p50": float(np.median(s["p50"])),
                "median_null_p99": float(np.median(s["p99"])),
                "max_null_p99": float(np.max(s["p99"])),
                "genes": len(s["p99"]),
            }
        )
    df = pd.DataFrame(rows)
    ref = df["median_null_p99"].median()
    df["p99_vs_lineage_median"] = df["median_null_p99"] / ref
    return df
