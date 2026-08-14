"""Cell QC, depth matching and detection-rate computation.

Pure numpy/scipy/anndata (R10 -- scanpy's numba path is unavailable on the
target machine, so normalisation and depth matching are implemented here).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def as_csr(X):
    if sp.issparse(X):
        return X.tocsr()
    return sp.csr_matrix(X)


def verify_integer_counts(X, sample_rows=2000, seed=0):
    """Return (is_integer, max_value, n_checked). Raw counts are required."""
    X = as_csr(X)
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    rows = np.arange(n) if n <= sample_rows else rng.choice(n, sample_rows, replace=False)
    sub = X[rows]
    data = sub.data
    if data.size == 0:
        return True, 0.0, int(rows.size)
    frac_ok = np.allclose(data, np.rint(data))
    return bool(frac_ok), float(data.max()), int(rows.size)


def cell_qc(X, gene_symbols, min_genes=200, max_mito_frac=0.10,
            mito_prefix="MT-"):
    """Cell-level QC mask plus the per-cell quantities the QC table needs."""
    X = as_csr(X)
    n_genes = np.diff(X.indptr)
    total = np.asarray(X.sum(axis=1)).ravel()
    mito_cols = np.array([i for i, s in enumerate(gene_symbols)
                          if str(s).upper().startswith(mito_prefix)], dtype=int)
    if mito_cols.size:
        mito = np.asarray(X[:, mito_cols].sum(axis=1)).ravel()
    else:
        mito = np.zeros(X.shape[0])
    with np.errstate(invalid="ignore", divide="ignore"):
        mito_frac = np.where(total > 0, mito / np.maximum(total, 1), 0.0)
    keep = (n_genes >= min_genes) & (mito_frac < max_mito_frac)
    return keep, {"n_genes": n_genes, "total_umi": total, "mito_frac": mito_frac,
                  "n_mito_genes_in_matrix": int(mito_cols.size)}


def downsample_columns(X, cols, depth, seed=0):
    """Exact multivariate-hypergeometric downsampling to a fixed UMI depth.

    Every cell's whole transcriptome is thinned to `depth` reads without
    replacement; the counts returned are those of `cols` under that thinning.
    Implemented as a sequential hypergeometric over the requested columns
    (vectorised across cells), with all remaining genes absorbed into the
    "other" pool -- which is exactly the marginal of full-transcriptome
    downsampling, without materialising the full thinned matrix.

    Cells whose total UMI is below `depth` cannot be thinned to it and are
    reported by the caller as differential loss (R6).

    Returns (counts[n_cells_kept, len(cols)], keep_mask over input rows).
    """
    X = as_csr(X)
    total = np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    keep = total >= depth
    if keep.sum() == 0:
        return np.zeros((0, len(cols)), dtype=np.int64), keep

    Xk = X[keep]
    remaining_total = np.asarray(Xk.sum(axis=1)).ravel().astype(np.int64)
    remaining_draws = np.full(remaining_total.shape, int(depth), dtype=np.int64)
    rng = np.random.default_rng(seed)

    sub = Xk[:, np.asarray(cols, dtype=int)]
    sub = np.asarray(sub.todense(), dtype=np.int64) if sp.issparse(sub) else np.asarray(sub, np.int64)

    out = np.zeros_like(sub)
    for j in range(sub.shape[1]):
        ngood = sub[:, j]
        nbad = remaining_total - ngood
        nsample = remaining_draws
        # guard the degenerate cells: nothing left to draw, or nothing to draw from
        active = (nsample > 0) & (ngood > 0) & (remaining_total > 0)
        k = np.zeros_like(ngood)
        if active.any():
            ns = np.minimum(nsample[active], remaining_total[active])
            k[active] = rng.hypergeometric(ngood[active], np.maximum(nbad[active], 0), ns)
        out[:, j] = k
        remaining_draws -= k
        remaining_total -= ngood
        np.maximum(remaining_draws, 0, out=remaining_draws)
        np.maximum(remaining_total, 0, out=remaining_total)
    return out, keep


def downsample_marginal(X, cols, depth, seed=0, block=2000):
    """Per-gene MARGINAL of full-transcriptome downsampling to `depth`.

    For a cell with total T and c_g copies of gene g, thinning the whole cell
    to D reads gives g the marginal Hypergeometric(c_g, T - c_g, D). That
    marginal is the same for every gene regardless of ordering, so when only
    per-gene detection rates are needed -- as in the null distribution, where
    each gene is evaluated independently -- the genes can be drawn
    independently instead of sequentially.

    Identical per-gene distribution to `downsample_columns`, but blockable and
    without the O(n_cols) sequential state, which is what makes a 10,000-gene
    null affordable. It does NOT preserve the joint constraint that the drawn
    counts sum to `depth`, so it must NOT be used where cross-gene structure
    matters.

    Returns (counts[n_kept, len(cols)], keep_mask).
    """
    X = as_csr(X)
    total = np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    keep = total >= depth
    if keep.sum() == 0:
        return np.zeros((0, len(cols)), dtype=np.int64), keep
    Xk = X[keep]
    T = total[keep]
    rng = np.random.default_rng(seed)
    cols = np.asarray(cols, dtype=int)
    out = np.zeros((Xk.shape[0], cols.size), dtype=np.int64)
    for start in range(0, cols.size, block):
        cc = cols[start:start + block]
        sub = Xk[:, cc]
        sub = np.asarray(sub.todense(), dtype=np.int64) if sp.issparse(sub) else \
            np.asarray(sub, dtype=np.int64)
        ngood = sub
        nbad = T[:, None] - ngood
        active = (ngood > 0) & (nbad >= 0)
        blk = np.zeros_like(ngood)
        if active.any():
            blk[active] = rng.hypergeometric(
                ngood[active], nbad[active],
                np.broadcast_to(np.minimum(depth, T)[:, None], ngood.shape)[active])
        out[:, start:start + block] = blk
    return out, keep


def detection_and_cpm(counts, depth):
    """Detection rate (count >= 1) and mean CPM from depth-matched counts."""
    if counts.shape[0] == 0:
        z = np.full(counts.shape[1], np.nan)
        return np.zeros(counts.shape[1], dtype=int), z, z
    n_detected = (counts >= 1).sum(axis=0).astype(int)
    detection_rate = n_detected / counts.shape[0]
    mean_cpm = (counts / float(depth) * 1e6).mean(axis=0)
    return n_detected, detection_rate, mean_cpm


def choose_depth_floor(total_umi, quantile=0.10):
    """Depth floor from the observed distribution of the target cells."""
    t = np.asarray(total_umi)
    t = t[np.isfinite(t) & (t > 0)]
    if t.size == 0:
        return 0
    return int(np.floor(np.quantile(t, quantile)))


def triple_positive_fraction(X, cols):
    """Fraction of cells with count >= 1 for EVERY column in `cols`.

    Used for the TRBC2+CD3E+CD3D+ T-cell gate. Single-gene lineage calls are
    barred: TRBC1 alone has previously called 95.6% of NK cells positive.
    """
    X = as_csr(X)
    if X.shape[0] == 0 or len(cols) == 0:
        return float("nan"), 0
    sub = X[:, np.asarray(cols, dtype=int)]
    sub = np.asarray(sub.todense()) if sp.issparse(sub) else np.asarray(sub)
    hit = (sub >= 1).all(axis=1)
    return float(hit.mean()), int(hit.sum())
