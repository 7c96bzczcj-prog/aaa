"""
Which negative markers actually discriminate NK from T in these libraries?

The task says to verify TRBC1 in the data before using it, and not to use it if
it is unreliable. The same test has to be applied to every negative marker,
because assuming a marker works is exactly the failure mode that check guards
against.

Method, and it is deliberately not circular: the reference NK cluster is
identified using POSITIVE markers only (KLRF1 / NCAM1 / KLRD1 plus the cytotoxic
granule genes). No negative marker takes part in choosing it. Each candidate
negative marker is then scored by how far it falls in that cluster relative to
the T-cell level in the same dataset:

    discrimination = mean_expression_in_NK_cluster / T_reference_level

where the T reference is the 90th percentile of per-cluster means across the
T/NK compartment, i.e. the level in a genuine T cluster. A usable negative
marker sits far below its T level in NK cells; one that stays near it is
carrying ambient or genuine NK expression and cannot be used to exclude T cells.
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

# candidate T-lineage negative markers, spanning a wide abundance range
CANDIDATES = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2",
              "CD5", "CD6", "TRAT1", "THEMIS", "ITK", "CD28", "LCK", "SKAP1",
              "CD40LG", "IL7R"]
# a positive control: CD2 is expressed by NK as well as T, so it must FAIL
CONTROL = ["CD2"]

USABLE_MAX = 0.25          # NK level must be below a quarter of the T level
PAN_T_MIN = 0.80           # and the marker must be positive in >=80% of T clusters


def cluster_means(a, genes):
    M = a.layers["lognorm"]
    idx = {g: a.var_names.get_loc(g) for g in genes if g in a.var_names}
    rows = {}
    for k in a.obs.tnk_cluster.cat.categories if hasattr(a.obs.tnk_cluster, "cat") \
            else sorted(a.obs.tnk_cluster.unique()):
        sel = np.where(a.obs.tnk_cluster.values == k)[0]
        if len(sel) == 0:
            continue
        d = {"n": len(sel)}
        for g, j in idx.items():
            v = M[sel][:, j]
            v = v.toarray().ravel() if sp.issparse(v) else np.asarray(v).ravel()
            d[g] = float(v.mean())
        rows[k] = d
    return pd.DataFrame(rows).T


def main():
    datasets = [d for d in ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]
                if os.path.exists(f"{CACHE}/tnk_{d}.h5ad")]
    allrows = []
    for ds in datasets:
        a = sc.read_h5ad(f"{CACHE}/tnk_{ds}.h5ad")
        genes = CANDIDATES + CONTROL + cc.NK_POS
        cm = cluster_means(a, genes)

        # ---- pick the NK reference cluster from POSITIVE markers only --------
        def col(g):
            return cm[g] if g in cm else pd.Series(0.0, index=cm.index)

        # Anchor on NCAM1 (CD56), requiring cytotoxic granules.
        #
        # Ranking on KLRF1 instead picks the wrong cluster: CD8 TEMRA / "NK-like
        # T" cells carry KLRF1, FCGR3A, GNLY and PRF1 while remaining genuine
        # CD3D+/CD3G+ T cells. Anchoring on those would make every CD3 gene look
        # non-discriminative, which is an artefact of the anchor, not a property
        # of the marker. CD56 is the canonical NK definition and selects the
        # FCGR3A-low, CD3D-low cluster in every dataset here.
        cyto_ok = (col("NKG7") > col("NKG7").median()) & (col("GNLY") > col("GNLY").median())
        cand = col("NCAM1").where(cyto_ok, -1)
        nk_cl = cand.idxmax()
        print(f"\n=== {ds} ===")
        print(f"NK reference cluster (positive markers only): {nk_cl} "
              f"(n={int(cm.loc[nk_cl,'n'])} cells)")
        print("  KLRF1=%.3f NCAM1=%.3f KLRD1=%.3f NKG7=%.3f GNLY=%.3f PRF1=%.3f" % (
            col("KLRF1")[nk_cl], col("NCAM1")[nk_cl], col("KLRD1")[nk_cl],
            col("NKG7")[nk_cl], col("GNLY")[nk_cl], col("PRF1")[nk_cl]))

        for g in CANDIDATES + CONTROL:
            if g not in cm:
                allrows.append({"dataset": ds, "gene": g, "status": "absent"})
                continue
            ref = float(np.percentile(cm[g].values, 90))
            nkv = float(cm.loc[nk_cl, g])
            # Pan-T consistency. A negative marker must be expressed by *every* T
            # subset, not most of them. CD28, IL7R and CD40LG are all low in CD8
            # effector/TEMRA cells, so using them would let an effector-T cluster
            # accumulate "negative" votes and be miscalled NK -- the exact error
            # the negative panel exists to prevent.
            others = cm.drop(index=nk_cl)[g].values
            pan_t = float((others > USABLE_MAX * ref).mean()) if ref > 1e-6 else np.nan
            allrows.append({
                "dataset": ds, "gene": g,
                "NK_cluster_mean": round(nkv, 3),
                "T_reference_p90": round(ref, 3),
                "ratio_NK_over_T": round(nkv / ref, 3) if ref > 1e-6 else np.nan,
                "frac_T_clusters_positive": round(pan_t, 3),
                "is_control": g in CONTROL,
            })

    df = pd.DataFrame(allrows)
    df.to_csv(f"{OUT}/marker_validation.csv", index=False)

    print("\n\n=== negative-marker validation ===")
    ratio = df.pivot_table(index="gene", columns="dataset",
                           values="ratio_NK_over_T", aggfunc="first")
    pan = df.pivot_table(index="gene", columns="dataset",
                         values="frac_T_clusters_positive", aggfunc="first")
    tab = pd.DataFrame({
        "worst_ratio_NK_over_T": ratio.max(axis=1),
        "worst_frac_T_positive": pan.min(axis=1),
    }).sort_values("worst_ratio_NK_over_T")
    tab["low_in_NK"] = tab.worst_ratio_NK_over_T < USABLE_MAX
    tab["pan_T"] = tab.worst_frac_T_positive >= PAN_T_MIN
    tab["USABLE"] = tab.low_in_NK & tab.pan_T & ~tab.index.isin(CONTROL)
    print(tab.round(3).to_string())

    usable = tab.index[tab.USABLE].tolist()
    print(f"\nUSABLE negative markers "
          f"(NK level < {USABLE_MAX} of T level, AND positive in "
          f">={PAN_T_MIN:.0%} of T clusters):")
    print("   " + (", ".join(usable) if usable else "none"))
    for g in tab.index[~tab.USABLE & ~tab.index.isin(CONTROL)]:
        why = []
        if not tab.loc[g, "low_in_NK"]:
            why.append(f"stays at {tab.loc[g,'worst_ratio_NK_over_T']:.2f} of T level inside NK")
        if not tab.loc[g, "pan_T"]:
            why.append(f"positive in only {tab.loc[g,'worst_frac_T_positive']:.0%} of T clusters")
        print(f"   REJECT {g:8s} - " + "; ".join(why))
    if "CD2" in tab.index:
        print(f"\ncontrol CD2 (NK express it, so must be rejected): "
              f"usable={tab.loc['CD2','USABLE']}")


if __name__ == "__main__":
    main()
