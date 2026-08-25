"""
Per-NK-cell effector-gene scoring, shared verbatim by the cervical and lung arms.

Both arms import this module, so the score is computed by the same code on the
same gene list; nothing about the scoring differs between cancer types.

Depth normalisation
-------------------
Three scores are produced per cell, and the figures state which is used:

  score_cp10k   PRIMARY. Counts are scaled to a fixed library size of 10,000
                (counts per 10k, "CP10K"), log1p-transformed, and averaged over
                the six effector genes. This is the standard scRNA-seq depth
                normalisation.

  frac_effector Compositional control: effector UMIs as a fraction of the cell's
                total UMIs, no log. Immune to the log's variance-stabilising
                distortion, and reads directly as "share of the transcriptome
                spent on effector genes".

  score_ds      Depth-matched control. Every cell is first downsampled by
                binomial thinning to a common total UMI count, then scored as
                CP10K. CP10K rescales but does not remove depth: dropout is
                depth-dependent, so a deeply sequenced cell detects more of a
                low-abundance transcript than a shallow one at identical true
                expression. Since the two cancer types come from different
                studies at different depths, this is the control that decides
                whether a shape difference is biology or sequencing.

Counts must be the ambient-corrected (decontX) counts, matching the cervical
pipeline: GZMB, PRF1 and NKG7 are among the most prominent genes in the ambient
pool, so scoring uncorrected counts measures the soup as much as the cell.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

# the gene list, fixed, identical for every dataset
EFFECTOR_GENES = ["GZMB", "PRF1", "GZMA", "IFNG", "NKG7", "KLRD1"]

# common depth for the downsampling control
DOWNSAMPLE_TO = 1500


def _col(M, j):
    c = M[:, j]
    return c.toarray().ravel() if sp.issparse(c) else np.asarray(c).ravel()


def downsample_counts(X, target, seed=0):
    """
    Binomial thinning of each cell to `target` total counts.

    Cells already below target are left alone and flagged, because thinning
    upward is not possible and silently keeping them would mix depths in the
    'depth-matched' comparison.
    """
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(X).astype(np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    keep = tot >= target
    out = X.copy()
    for i in np.where(keep)[0]:
        s, e = out.indptr[i], out.indptr[i + 1]
        d = out.data[s:e]
        n = int(d.sum())
        if n <= target:
            continue
        # multivariate hypergeometric via successive binomials
        p = d / d.sum()
        out.data[s:e] = rng.multinomial(target, p).astype(np.float64)
    return out, keep


def effector_table(counts, var_names, obs, genes=EFFECTOR_GENES,
                   target_sum=1e4, downsample_to=DOWNSAMPLE_TO, seed=0):
    """
    Per-cell effector scores.

    counts    : cells x genes, ambient-corrected counts (CSR)
    var_names : gene symbols aligned to columns
    obs       : DataFrame indexed like the cells; carried through to the output
    """
    X = sp.csr_matrix(counts).astype(np.float64)
    var_names = np.asarray(var_names)
    present = [g for g in genes if g in set(var_names)]
    missing = [g for g in genes if g not in set(var_names)]
    idx = {g: int(np.where(var_names == g)[0][0]) for g in present}

    tot = np.asarray(X.sum(axis=1)).ravel()
    tot_safe = np.where(tot > 0, tot, 1.0)

    # --- primary: CP10K + log1p, averaged over the gene list ----------------
    per_gene = {}
    acc = np.zeros(X.shape[0])
    for g in present:
        raw = _col(X, idx[g])
        ln = np.log1p(raw / tot_safe * target_sum)
        per_gene[f"ln_{g}"] = ln
        acc += ln
    score = acc / max(len(present), 1)

    # --- compositional control ----------------------------------------------
    eff_umi = np.zeros(X.shape[0])
    for g in present:
        eff_umi += _col(X, idx[g])
    frac = eff_umi / tot_safe

    # --- depth-matched control ----------------------------------------------
    Xd, ds_ok = downsample_counts(X, downsample_to, seed=seed)
    totd = np.asarray(Xd.sum(axis=1)).ravel()
    totd_safe = np.where(totd > 0, totd, 1.0)
    accd = np.zeros(X.shape[0])
    for g in present:
        accd += np.log1p(_col(Xd, idx[g]) / totd_safe * target_sum)
    score_ds = accd / max(len(present), 1)

    out = pd.DataFrame(
        {
            "total_counts_decont": tot,
            "n_effector_umi": eff_umi,
            "score_cp10k": score,
            "frac_effector": frac,
            "score_ds": np.where(ds_ok, score_ds, np.nan),
            "depth_ok": ds_ok,
            **per_gene,
        },
        index=obs.index,
    )
    for c in obs.columns:
        out[c] = obs[c].values
    out.attrs["genes_used"] = present
    out.attrs["genes_missing"] = missing
    return out


def summarise(tab, group_cols=("dataset", "cancer", "tissue_group", "donor"),
              score="score_cp10k", min_cells=10):
    """
    Per-donor summary. Never pooled across donors: each donor contributes one
    row, and the group-level statistics in the figures are taken over those
    donor rows, not over cells.
    """
    rows = []
    for key, g in tab.groupby(list(group_cols), observed=True):
        if len(g) < 1:
            continue
        d = dict(zip(group_cols, key if isinstance(key, tuple) else (key,)))
        v = g[score].dropna()
        d.update(
            n_cells=len(g),
            median=float(np.median(v)) if len(v) else np.nan,
            q25=float(np.percentile(v, 25)) if len(v) else np.nan,
            q75=float(np.percentile(v, 75)) if len(v) else np.nan,
            mean=float(v.mean()) if len(v) else np.nan,
            median_depth=float(np.median(g.total_counts_decont)),
            median_frac_effector=float(np.median(g.frac_effector)),
            usable=len(g) >= min_cells,
        )
        rows.append(d)
    return pd.DataFrame(rows)
