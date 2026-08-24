"""
Re-apply the validated NK gate to every dataset, from cache.

QC, doublet removal, decontX and clustering are unchanged and expensive, so this
reuses the cached T/NK objects and only redoes the call. Writes the same
per-sample outputs build_dataset.py does, so make_table.py is unaffected.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

CACHE = "/home/user/cervical_work/cache"
OUT = "/home/user/aaa/results/cervical"
DATASETS = ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]


def cd56_bright(tnk, nk_mask):
    """CD56bright = NCAM1-detected, FCGR3A-absent, among gated NK cells."""
    ln = tnk.layers["lognorm"]

    def vec(g):
        if g not in tnk.var_names:
            return np.zeros(tnk.n_obs)
        c = ln[:, tnk.var_names.get_loc(g)]
        return c.toarray().ravel() if sp.issparse(c) else np.asarray(c).ravel()

    return nk_mask & (vec("NCAM1") > 0) & (vec("FCGR3A") <= 0)


def main():
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        obsp = f"{OUT}/obs_{ds}.csv.gz"
        if not (os.path.exists(p) and os.path.exists(obsp)):
            print(f"[skip] {ds}: cache not ready")
            continue
        tnk = sc.read_h5ad(p)
        tab, nk_clusters = cc.call_nk_clusters(tnk)
        tab.to_csv(f"{OUT}/nkgate_{ds}.csv")

        keys = np.asarray(tnk.obs["tnk_cluster"].values).astype(str)
        nk_mask = np.isin(keys, list(nk_clusters))
        stage1 = int(nk_mask.sum())

        # second stage: drop internally-mixed subclusters lacking NK receptors
        nk_mask, subtab = cc.refine_nk_subclusters(tnk, nk_mask)
        if subtab is not None:
            subtab.to_csv(f"{OUT}/nksub_{ds}.csv")

        bright = cd56_bright(tnk, nk_mask)

        print(f"\n=== {ds} ===")
        print(f"NK clusters: {sorted(nk_clusters)}  "
              f"-> stage1 {stage1} -> refined {int(nk_mask.sum())} "
              f"of {tnk.n_obs} T/NK cells "
              f"({100*(stage1-int(nk_mask.sum()))/max(stage1,1):.0f}% dropped as receptor-negative)")
        if subtab is not None:
            print(subtab.round(3).to_string())
        show = [c for c in ["NCAM1", "KLRF1", "KLRD1", "NKG7", "GNLY", "PRF1",
                            "CD3D", "CD3G", "CD6", "CD3E", "TRAC",
                            "n_T_neg", "cytotoxic", "receptor", "is_NK"] if c in tab]
        print(tab[show].round(3).to_string())

        obs = pd.read_csv(obsp, index_col=0)
        obs["is_NK"] = False
        obs.loc[tnk.obs_names[nk_mask], "is_NK"] = True
        obs["is_bright"] = False
        obs.loc[tnk.obs_names[bright], "is_bright"] = True
        obs["is_immune"] = obs.lineage.isin(cc.IMMUNE_LINEAGES)
        obs["is_lymph"] = obs.lineage.isin(cc.LYMPHOID_LINEAGES)

        rows = []
        for smp, g in obs.groupby("sample", observed=True):
            nk = int(g.is_NK.sum())
            rows.append({
                "dataset": ds, "sample": smp, "donor": g.donor.iloc[0],
                "tissue": g.tissue.iloc[0], "histology": g.histology.iloc[0],
                "n_cells_total": len(g),
                "n_immune": int(g.is_immune.sum()),
                "n_lymphocyte": int(g.is_lymph.sum()),
                "n_TNK": int((g.lineage == "T_NK").sum()),
                "n_NK": nk,
                "pct_of_all": 100 * nk / len(g) if len(g) else np.nan,
                "pct_of_immune": 100 * nk / g.is_immune.sum() if g.is_immune.sum() else np.nan,
                "pct_of_lymphocyte": 100 * nk / g.is_lymph.sum() if g.is_lymph.sum() else np.nan,
                "n_NK_bright": int(g.is_bright.sum()),
                "pct_bright_of_NK": 100 * g.is_bright.sum() / nk if nk else np.nan,
                "decontx_theta_median": float(g.decontx_theta.median()),
                "median_UMI": float(g.total_counts.median()),
            })
        res = pd.DataFrame(rows)
        res.to_csv(f"{OUT}/nk_per_sample_{ds}.csv", index=False)
        print(res[["sample", "tissue", "n_cells_total", "n_TNK", "n_NK",
                   "pct_of_all", "pct_of_immune", "pct_of_lymphocyte",
                   "pct_bright_of_NK"]].round(2).to_string(index=False))

        # keep the per-cell NK flags so make_table's capture check can use them
        tnk.obs["is_NK"] = nk_mask
        tnk.write(p)


if __name__ == "__main__":
    main()
