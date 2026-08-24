"""
Why is CD3E still at ~64% of the T level inside the CD56+ cluster?

Three explanations, with opposite consequences for the NK count:

  ambient   decontX systematically under-removes high-abundance T transcripts.
            CD3E is then not evidence of T identity and must not gate.
  doublet   the CD56+ cluster carries residual T-NK doublets. CD3E IS evidence
            of T identity, and the current gate is passing contaminated cells.
  real      NK cells genuinely transcribe CD3E. It cannot gate, and the
            published CD3E-negative definition of NK is wrong in scRNA data.

Each leaves a different fingerprint, and they are separable without assuming the
answer:

  1. ambient scales with its source. Per sample, the CD3E detection rate inside
     the candidate NK cells should track that sample's T-cell fraction. Real
     expression should not.
  2. doublets carry a size signature. CD3E+ NK cells would show higher UMI,
     higher gene count and higher scrublet score than CD3E- NK cells.
  3. doublets carry the whole T module, not one gene. CD3E+ NK cells would also
     carry CD3D/CD3G far above chance.
  4. within-cluster retest. Cluster-mean differences are partly a product of the
     clustering itself, so the CD56+ cells are re-clustered on their own. A
     contaminated cluster splits into a CD3-high and a CD3-low subcluster; real
     or ambient expression stays spread across subclusters.

The candidate NK cluster is chosen on POSITIVE markers only. No CD3 gene takes
part in selecting the cells whose CD3 is then interrogated.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

CACHE = "/home/user/cervical_work/cache"
OUT = "/home/user/aaa/results/cervical"
DATASETS = ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]


def vec(a, g, layer="lognorm"):
    if g not in a.var_names:
        return np.zeros(a.n_obs)
    M = a.layers[layer] if layer else a.X
    c = M[:, a.var_names.get_loc(g)]
    return c.toarray().ravel() if sp.issparse(c) else np.asarray(c).ravel()


def pick_nk_cluster(a):
    """Positive markers only: NCAM1 anchor, cytotoxic-granule requirement."""
    ex = cc.cluster_expression(a, "tnk_cluster",
                               ["NCAM1", "KLRF1", "KLRD1", "NKG7", "GNLY", "PRF1"])
    cyto = (ex["NKG7"] > ex["NKG7"].median()) & (ex["GNLY"] > ex["GNLY"].median())
    return ex["NCAM1"].where(cyto, -1).idxmax()


def main():
    summary = []
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        if not os.path.exists(p):
            print(f"[skip] {ds}")
            continue
        a = sc.read_h5ad(p)
        keys = np.asarray(a.obs["tnk_cluster"].values).astype(str)
        nk_cl = pick_nk_cluster(a)
        nk = np.where(keys == nk_cl)[0]
        print("\n" + "=" * 78)
        print(f"{ds}: candidate NK cluster {nk_cl}, n={len(nk)}")

        cd3e = vec(a, "CD3E")[nk]
        cd3d = vec(a, "CD3D")[nk]
        cd3g = vec(a, "CD3G")[nk]
        pos_e = cd3e > 0

        # ---- 1. does CD3E track the T-cell fraction of its sample? ----------
        smp = np.asarray(a.obs["sample"].values).astype(str)
        rows = []
        for s in np.unique(smp):
            in_s = smp == s
            n_nk_s = int((keys[in_s] == nk_cl).sum())
            if n_nk_s < 15:
                continue
            t_frac = float((keys[in_s] != nk_cl).sum() / in_s.sum())
            sel = np.where(in_s & (keys == nk_cl))[0]
            rows.append({
                "sample": s, "n_NK": n_nk_s,
                "T_frac_of_TNK": round(t_frac, 3),
                "CD3E_det_in_NK": round(float((vec(a, "CD3E")[sel] > 0).mean()), 3),
                "CD3D_det_in_NK": round(float((vec(a, "CD3D")[sel] > 0).mean()), 3),
                "theta": round(float(a.obs["decontx_theta"].values[sel].mean()), 3),
            })
        st = pd.DataFrame(rows)
        if len(st) >= 4:
            r_e, p_e = stats.spearmanr(st.T_frac_of_TNK, st.CD3E_det_in_NK)
            r_d, p_d = stats.spearmanr(st.T_frac_of_TNK, st.CD3D_det_in_NK)
            print(f"\n  [1] ambient scaling across {len(st)} samples")
            print(st.to_string(index=False))
            print(f"      CD3E vs T-fraction : rho={r_e:+.3f} p={p_e:.3f}")
            print(f"      CD3D vs T-fraction : rho={r_d:+.3f} p={p_d:.3f}   (near-zero gene, control)")
        else:
            r_e = p_e = np.nan
            print("  [1] too few samples with >=15 NK cells")

        # ---- 2. do CD3E+ NK cells look like doublets? -----------------------
        tot = a.obs["total_counts"].values[nk]
        ngen = a.obs["n_genes_by_counts"].values[nk]
        dsc = (a.obs["doublet_score"].values[nk]
               if "doublet_score" in a.obs else np.full(len(nk), np.nan))
        print(f"\n  [2] doublet signature  (CD3E+ n={int(pos_e.sum())} vs "
              f"CD3E- n={int((~pos_e).sum())})")
        for name, arr in [("median UMI", tot), ("median genes", ngen),
                          ("median doublet score", dsc)]:
            if np.all(np.isnan(arr)):
                continue
            hi, lo = np.median(arr[pos_e]), np.median(arr[~pos_e])
            try:
                u_p = stats.mannwhitneyu(arr[pos_e], arr[~pos_e]).pvalue
            except ValueError:
                u_p = np.nan
            print(f"      {name:22s} CD3E+ {hi:9.1f}   CD3E- {lo:9.1f}   "
                  f"ratio {hi/lo if lo else np.nan:5.2f}   p={u_p:.3g}")

        # ---- 3. does CD3E travel with the rest of the CD3 module? -----------
        both = float(((cd3d > 0) & pos_e).mean())
        exp_ind = float((cd3d > 0).mean() * pos_e.mean())
        print(f"\n  [3] CD3 module coherence")
        print(f"      P(CD3D+)={float((cd3d>0).mean()):.3f}  P(CD3E+)={float(pos_e.mean()):.3f}")
        print(f"      observed P(CD3D+ and CD3E+)={both:.4f}  "
              f"independent expectation={exp_ind:.4f}  "
              f"lift={both/exp_ind if exp_ind else np.nan:.2f}")
        print(f"      CD3D+ rate among CD3E+ cells: "
              f"{float((cd3d[pos_e]>0).mean()) if pos_e.sum() else np.nan:.3f}"
              f"   among CD3E- cells: "
              f"{float((cd3d[~pos_e]>0).mean()) if (~pos_e).sum() else np.nan:.3f}")

        # ---- 4. within-cluster retest ---------------------------------------
        sub = a[nk].copy()
        sub.X = sub.layers["lognorm"].copy()
        print(f"\n  [4] within-cluster retest: re-clustering the {sub.n_obs} CD56+ cells")
        try:
            sc.pp.highly_variable_genes(sub, n_top_genes=1500)
            s2 = sub[:, sub.var.highly_variable].copy()
            sc.pp.scale(s2, max_value=10)
            sc.tl.pca(s2, n_comps=min(20, s2.n_obs - 1, s2.n_vars - 1),
                      svd_solver="arpack")
            sc.pp.neighbors(s2, n_neighbors=min(15, s2.n_obs - 1))
            sc.tl.leiden(s2, resolution=0.6, key_added="sub",
                         flavor="igraph", n_iterations=2, directed=False)
            sub.obs["sub"] = s2.obs["sub"].values
            ex = cc.cluster_expression(sub, "sub",
                                       ["CD3E", "CD3D", "CD3G", "CD6", "TRAC",
                                        "NCAM1", "KLRF1", "KLRD1", "NKG7", "GNLY",
                                        "FCGR3A", "IL7R"])
            ex["n"] = sub.obs.groupby("sub", observed=True).size()
            ex["med_UMI"] = sub.obs.groupby("sub", observed=True).total_counts.median()
            if "doublet_score" in sub.obs:
                ex["dbl"] = sub.obs.groupby("sub", observed=True).doublet_score.median()
            print(ex.round(3).to_string())
            spread = ex["CD3E"].max() / max(ex["CD3E"].min(), 1e-6)
            print(f"      CD3E max/min across NK subclusters = {spread:.1f}"
                  f"   (a contaminated cluster splits; ambient/real does not)")
        except Exception as e:
            print(f"      subclustering failed: {e}")
            spread = np.nan

        summary.append({"dataset": ds, "nk_cluster": nk_cl, "n_NK": len(nk),
                        "rho_CD3E_vs_Tfrac": round(float(r_e), 3) if r_e == r_e else np.nan,
                        "p_CD3E_vs_Tfrac": round(float(p_e), 4) if p_e == p_e else np.nan,
                        "CD3E_pos_frac": round(float(pos_e.mean()), 3),
                        "CD3D_lift_in_CD3Epos": round(both / exp_ind, 2) if exp_ind else np.nan,
                        "CD3E_subcluster_spread": round(float(spread), 2) if spread == spread else np.nan})

    if summary:
        df = pd.DataFrame(summary)
        df.to_csv(f"{OUT}/cd3e_discriminant.csv", index=False)
        print("\n\n=== SUMMARY ===")
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
