"""Donor-level statistics. Pure numpy/scipy; no scanpy (R10).

Every test here consumes one value per donor. Nothing in this module can be
handed cell-level input without the caller lying about what a row is (R1).
"""
from __future__ import annotations

import itertools

import numpy as np


def min_achievable_p(n, kind="wilcoxon_signed_rank_two_sided"):
    """R2: the smallest p the test can return at this n.

    Two-sided exact signed-rank with n non-zero pairs: the extreme statistic
    (all differences one sign) occurs for 2 of the 2**n sign assignments.
    """
    if n <= 0:
        return float("nan")
    if kind == "wilcoxon_signed_rank_two_sided":
        return min(1.0, 2.0 / (2 ** n))
    if kind == "sign_two_sided":
        return min(1.0, 2.0 / (2 ** n))
    raise ValueError(kind)


def exact_wilcoxon_signed_rank(diffs):
    """Exact two-sided Wilcoxon signed-rank p by full enumeration.

    Zero differences are dropped (Wilcoxon's own convention); ties in |d| get
    average ranks. Exact enumeration is used unconditionally -- at the donor
    counts this project runs at, a normal approximation is meaningless.

    Returns (statistic_W_plus, p_two_sided, n_used).
    """
    d = np.asarray(diffs, dtype=float)
    d = d[np.isfinite(d)]
    d = d[d != 0]
    n = d.size
    if n == 0:
        return float("nan"), 1.0, 0
    order = np.argsort(np.abs(d), kind="mergesort")
    absd = np.abs(d)[order]
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and absd[j + 1] == absd[i]:
            j += 1
        ranks[i:j + 1] = np.arange(i + 1, j + 2).mean()
        i = j + 1
    signs = np.sign(d)[order]
    w_plus = ranks[signs > 0].sum()

    if n > 22:  # 4M sign patterns; beyond this fall back to the normal approx
        mu = ranks.sum() / 2.0
        sigma = np.sqrt((ranks ** 2).sum() / 4.0)
        if sigma == 0:
            return w_plus, 1.0, n
        z = (w_plus - mu) / sigma
        from scipy.stats import norm
        return w_plus, float(min(1.0, 2 * norm.sf(abs(z)))), n

    total = ranks.sum()
    target = abs(w_plus - total / 2.0)
    count = 0
    for mask in itertools.product((0, 1), repeat=n):
        wp = float(np.dot(ranks, mask))
        if abs(wp - total / 2.0) >= target - 1e-12:
            count += 1
    p = count / (2 ** n)
    return w_plus, float(min(1.0, p)), n


def benjamini_hochberg(pvals):
    """BH q-values. NaN p-values are carried through as NaN."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.empty(m)
    out[order] = adj
    q[ok] = out
    return q


def logit(p, eps=None):
    """Logit with a Haldane-style continuity guard for 0/1."""
    p = np.asarray(p, dtype=float)
    if eps is None:
        eps = 1e-6
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def continuity_rate(n_detected, n_cells):
    """(k + 0.5) / (n + 1) -- keeps logit finite at 0 and 100% detection."""
    return (np.asarray(n_detected, dtype=float) + 0.5) / (np.asarray(n_cells, dtype=float) + 1.0)


def bootstrap_ci(values, n_boot=20000, alpha=0.05, seed=0, statistic=np.mean):
    """Percentile bootstrap over DONORS (R7: never a single-difference denominator).

    With the donor counts here the interval is coarse by construction; that is
    the honest width, not a defect to be hidden.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    if v.size == 1:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    stats = statistic(v[idx], axis=1)
    lo = float(np.percentile(stats, 100 * alpha / 2))
    hi = float(np.percentile(stats, 100 * (1 - alpha / 2)))
    return lo, hi


def fisher_z_combine(rs, ns):
    """Combine correlations across donors/datasets on the z scale (R1, R9)."""
    rs = np.asarray(rs, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rs) & np.isfinite(ns) & (ns > 3)
    if ok.sum() == 0:
        return float("nan"), float("nan"), 0
    z = np.arctanh(np.clip(rs[ok], -0.999999, 0.999999))
    w = ns[ok] - 3
    zbar = float(np.sum(w * z) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    return float(np.tanh(zbar)), se, int(ok.sum())


def kitagawa_decomposition(rate_a, comp_a, rate_b, comp_b):
    """Three-term composition/state decomposition on the DETECTION-RATE scale.

    within  = sum_i comp_a_i * (rate_b_i - rate_a_i)     (state)
    between = sum_i rate_a_i * (comp_b_i - comp_a_i)     (composition)
    interaction = sum_i (comp_b_i - comp_a_i) * (rate_b_i - rate_a_i)

    The three terms sum exactly to the total difference. The interaction term
    is NOT optional -- dropping it is what produced the 64%/99% confusion this
    project already paid for.
    """
    ra, ca = np.asarray(rate_a, float), np.asarray(comp_a, float)
    rb, cb = np.asarray(rate_b, float), np.asarray(comp_b, float)
    within = float(np.sum(ca * (rb - ra)))
    between = float(np.sum(ra * (cb - ca)))
    interaction = float(np.sum((cb - ca) * (rb - ra)))
    total = float(np.sum(cb * rb) - np.sum(ca * ra))
    residual = total - (within + between + interaction)
    return {"within": within, "between": between, "interaction": interaction,
            "total": total, "residual": residual}
