"""Paired pseudobulk differential expression (Protocol Phase 4).

A faithful NumPy/SciPy reimplementation of the limma-voom pipeline
(Law et al. 2014, Genome Biology; Smyth 2004, SAGMB), because no R
toolchain is available in this environment.

Why voom rather than DESeq2: the protocol needs a *standard error* per
gene per lineage as a first-class output (Phase 5 runs TOST on the
confidence interval).  voom's precision weights plus limma's moderated
variance give exactly that, on the log2 scale the quadrant rules are
written in.

The protocol forbids ratios (Phase 4: "report the vector, not the
ratio"), so everything downstream consumes (log2FC, SE, df) triples.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.special import digamma, polygamma


# --------------------------------------------------------------------------
# Normalisation
# --------------------------------------------------------------------------
def tmm_norm_factors(
    counts: np.ndarray,
    ref_index: int | None = None,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
) -> np.ndarray:
    """TMM normalisation factors (Robinson & Oshlack 2010).

    `counts` is genes x samples.  Returns one factor per sample,
    normalised to geometric mean 1 (as edgeR's calcNormFactors does).
    """
    counts = np.asarray(counts, dtype=float)
    lib = counts.sum(axis=0)
    if np.any(lib <= 0):
        raise ValueError("every sample needs a non-zero library size")

    if ref_index is None:
        # edgeR picks the sample whose upper-quartile is closest to the mean
        f75 = np.array(
            [np.percentile(counts[:, j][counts[:, j] > 0], 75) if np.any(counts[:, j] > 0) else 0.0
             for j in range(counts.shape[1])]
        )
        f75 = f75 / lib
        ref_index = int(np.argmin(np.abs(f75 - f75.mean())))

    ref = counts[:, ref_index]
    ref_lib = lib[ref_index]
    factors = np.ones(counts.shape[1])

    for j in range(counts.shape[1]):
        obs = counts[:, j]
        keep = (obs > 0) & (ref > 0)
        if keep.sum() < 2:
            continue
        o, r = obs[keep], ref[keep]
        # log ratio (M) and mean abundance (A)
        m = np.log2((o / lib[j]) / (r / ref_lib))
        a = 0.5 * np.log2((o / lib[j]) * (r / ref_lib))
        # approximate asymptotic variance, used as the weight
        v = (lib[j] - o) / lib[j] / o + (ref_lib - r) / ref_lib / r

        finite = np.isfinite(m) & np.isfinite(a) & (a > -1e10)
        m, a, v = m[finite], a[finite], v[finite]
        if m.size < 2:
            continue

        n = m.size
        lo_m, hi_m = np.floor(n * logratio_trim) + 1, n - np.floor(n * logratio_trim)
        lo_a, hi_a = np.floor(n * sum_trim) + 1, n - np.floor(n * sum_trim)
        rank_m = stats.rankdata(m)
        rank_a = stats.rankdata(a)
        keep2 = (rank_m >= lo_m) & (rank_m <= hi_m) & (rank_a >= lo_a) & (rank_a <= hi_a)
        if keep2.sum() == 0:
            continue
        w = 1.0 / v[keep2]
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            continue
        factors[j] = 2.0 ** (np.sum(w * m[keep2]) / np.sum(w))

    return factors / np.exp(np.mean(np.log(factors)))


# --------------------------------------------------------------------------
# Weighted least squares
# --------------------------------------------------------------------------
def _wls(y: np.ndarray, X: np.ndarray, w: np.ndarray):
    """Per-gene weighted least squares.

    y, w are genes x samples; X is samples x coefficients.
    Returns (coefficients, residual variance, residual df, stdev_unscaled,
    fitted values).  Ranks are checked once on X, so a rank-deficient
    design fails loudly instead of silently returning garbage.
    """
    n_genes, n_samp = y.shape
    n_coef = X.shape[1]
    rank = np.linalg.matrix_rank(X)
    if rank < n_coef:
        raise ValueError(
            f"design matrix is rank deficient (rank {rank} < {n_coef} columns): "
            "condition is collinear with the blocking factor"
        )
    df_resid = n_samp - rank
    if df_resid < 1:
        raise ValueError("no residual degrees of freedom left after fitting the design")

    # Batched normal equations: one (n_coef x n_coef) system per gene.
    # Weights are gene-specific (voom), so this cannot collapse to a
    # single shared solve, but it does vectorise cleanly.
    XtWX = np.einsum("ni,gn,nj->gij", X, w, X, optimize=True)
    XtWy = np.einsum("ni,gn,gn->gi", X, w, y, optimize=True)

    try:
        beta = np.linalg.solve(XtWX, XtWy[..., None])[..., 0]
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        # a gene with all-zero weights can make its own system singular
        pinv = np.linalg.pinv(XtWX)
        beta = np.einsum("gij,gj->gi", pinv, XtWy)
        XtWX_inv = pinv

    fitted = beta @ X.T
    resid = y - fitted
    sigma2 = np.einsum("gn,gn->g", w, resid**2) / df_resid
    stdev_unscaled = np.sqrt(
        np.maximum(np.einsum("gii->gi", XtWX_inv), 0.0)
    )
    return beta, sigma2, df_resid, stdev_unscaled, fitted


# --------------------------------------------------------------------------
# voom
# --------------------------------------------------------------------------
def _lowess(x: np.ndarray, y: np.ndarray, frac: float = 0.5, n_iter: int = 0) -> np.ndarray:
    """LOWESS smoother returning fitted values at the input points."""
    import statsmodels.nonparametric.smoothers_lowess as sl

    out = sl.lowess(y, x, frac=frac, it=n_iter, return_sorted=True)
    return out  # (n, 2) sorted by x


def voom(
    counts: np.ndarray,
    design: np.ndarray,
    lib_size: np.ndarray | None = None,
    span: float = 0.5,
):
    """Estimate voom precision weights for log-CPM values.

    Returns (log_cpm, weights, lib_size).  Mirrors limma::voom: fit an
    unweighted model, regress sqrt(residual sd) on average log-count,
    then read the trend back at each observation's fitted log-count.
    """
    counts = np.asarray(counts, dtype=float)
    if lib_size is None:
        lib_size = counts.sum(axis=0) * tmm_norm_factors(counts)

    log_cpm = np.log2((counts + 0.5) / (lib_size + 1.0) * 1e6)

    ones = np.ones_like(counts)
    beta, sigma2, df_resid, _, fitted = _wls(log_cpm, design, ones)

    sy = np.sqrt(np.sqrt(np.maximum(sigma2, 0.0)))  # sqrt of residual sd
    # mean log2 count (not CPM) — the x-axis of the mean-variance trend
    sx = log_cpm.mean(axis=1) + np.log2(np.mean(lib_size + 1.0)) - np.log2(1e6)

    keep = np.isfinite(sx) & np.isfinite(sy) & (sy > 0)
    if keep.sum() < 10:
        # too few genes to fit a trend; fall back to unit weights
        return log_cpm, np.ones_like(counts), lib_size

    sm = _lowess(sx[keep], sy[keep], frac=span)
    xs, ys = sm[:, 0], sm[:, 1]

    # fitted log2 count for every observation
    lam = fitted + np.log2(lib_size + 1.0)[None, :] - np.log2(1e6)
    lam = np.clip(lam, xs.min(), xs.max())
    sy_hat = np.interp(lam, xs, ys)
    sy_hat = np.maximum(sy_hat, 1e-8)
    weights = 1.0 / sy_hat**4
    return log_cpm, weights, lib_size


# --------------------------------------------------------------------------
# Empirical Bayes moderation
# --------------------------------------------------------------------------
def _trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """Solve trigamma(y) = x for y (limma's trigammaInverse)."""
    x = np.asarray(x, dtype=float)
    y = np.where(x > 1e7, 1.0 / np.sqrt(np.maximum(x, 1e-300)), np.nan)
    y = np.where(x < 1e-6, 1.0 / np.maximum(x, 1e-300), y)
    todo = np.isnan(y)
    if np.any(todo):
        y0 = 0.5 + 1.0 / x[todo]
        for _ in range(50):
            tri = polygamma(1, y0)
            dif = tri * (1 - tri / x[todo]) / polygamma(2, y0)
            y0 = y0 + dif
            if np.max(np.abs(dif / y0)) < 1e-8:
                break
        y[todo] = y0
    return y


def fit_f_dist(sigma2: np.ndarray, df: float) -> tuple[float, float]:
    """Moment estimation of the prior (d0, s0^2) — limma's fitFDist."""
    s2 = np.asarray(sigma2, dtype=float)
    ok = np.isfinite(s2) & (s2 > 0)
    if ok.sum() < 2:
        return 0.0, float(np.nanmean(s2)) if ok.any() else 1.0
    z = np.log(s2[ok])
    e = z - digamma(df / 2.0) + np.log(df / 2.0)
    ebar = e.mean()
    n = e.size
    evar = np.sum((e - ebar) ** 2) / (n - 1) - polygamma(1, df / 2.0)
    if evar > 0:
        d0 = 2.0 * float(_trigamma_inverse(np.array([evar]))[0])
        s02 = float(np.exp(ebar + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0 = np.inf
        s02 = float(np.exp(ebar))
    return d0, s02


@dataclass
class DEResult:
    """Per-gene result for one lineage, one contrast."""

    genes: np.ndarray
    log2fc: np.ndarray
    se: np.ndarray          # moderated standard error
    t: np.ndarray
    pvalue: np.ndarray
    df_total: float
    ci_lo_95: np.ndarray
    ci_hi_95: np.ndarray
    ci_lo_90: np.ndarray
    ci_hi_90: np.ndarray
    ave_expr: np.ndarray
    sigma2_post: np.ndarray

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(
            {
                "gene": self.genes,
                "log2FC": self.log2fc,
                "SE": self.se,
                "t": self.t,
                "P": self.pvalue,
                "df_total": self.df_total,
                "CI95_lo": self.ci_lo_95,
                "CI95_hi": self.ci_hi_95,
                "CI90_lo": self.ci_lo_90,
                "CI90_hi": self.ci_hi_90,
                "AveExpr": self.ave_expr,
            }
        )


def paired_de(
    counts: np.ndarray,
    genes,
    patient: np.ndarray,
    condition: np.ndarray,
    case_label: str,
    robust_span: float = 0.5,
) -> DEResult:
    """Paired differential expression: ``~ patient + condition``.

    n is the number of *individuals*, not cells (protocol Phase 2.3).
    `condition` must be two-valued; `case_label` names the level whose
    effect is reported (case minus the other level).
    """
    counts = np.asarray(counts, dtype=float)
    genes = np.asarray(genes)
    patient = np.asarray(patient)
    condition = np.asarray(condition)

    levels = sorted(set(condition.tolist()))
    if len(levels) != 2:
        raise ValueError(f"condition must have exactly 2 levels, got {levels}")
    if case_label not in levels:
        raise ValueError(f"case_label {case_label!r} not among {levels}")

    # design: intercept + patient dummies (first patient absorbed) + condition
    pats = sorted(set(patient.tolist()))
    cols = [np.ones(len(condition))]
    for p in pats[1:]:
        cols.append((patient == p).astype(float))
    cond_col = (condition == case_label).astype(float)
    cols.append(cond_col)
    X = np.column_stack(cols)
    coef_index = X.shape[1] - 1

    log_cpm, weights, lib = voom(counts, X, span=robust_span)
    beta, sigma2, df_resid, stdev_unscaled, _ = _wls(log_cpm, X, weights)

    d0, s02 = fit_f_dist(sigma2, df_resid)
    if np.isinf(d0):
        s2_post = np.full_like(sigma2, s02)
        df_total = np.inf
    else:
        s2_post = (d0 * s02 + df_resid * sigma2) / (d0 + df_resid)
        df_total = df_resid + d0

    lfc = beta[:, coef_index]
    se = stdev_unscaled[:, coef_index] * np.sqrt(s2_post)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, lfc / se, 0.0)
    if np.isinf(df_total):
        p = 2 * stats.norm.sf(np.abs(t))
        q95, q90 = stats.norm.ppf(0.975), stats.norm.ppf(0.95)
    else:
        p = 2 * stats.t.sf(np.abs(t), df_total)
        q95, q90 = stats.t.ppf(0.975, df_total), stats.t.ppf(0.95, df_total)

    return DEResult(
        genes=genes,
        log2fc=lfc,
        se=se,
        t=t,
        pvalue=p,
        df_total=float(df_total),
        ci_lo_95=lfc - q95 * se,
        ci_hi_95=lfc + q95 * se,
        ci_lo_90=lfc - q90 * se,
        ci_hi_90=lfc + q90 * se,
        ave_expr=log_cpm.mean(axis=1),
        sigma2_post=s2_post,
    )
