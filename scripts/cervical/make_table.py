"""
Assemble the final per-donor NK table and the summary statistics.

Reads the per-dataset outputs written by build_dataset.py and produces:
  results/cervical/NK_TABLE.csv        the deliverable table
  results/cervical/NK_SUMMARY.csv      medians/ranges by tissue, per dataset
  results/cervical/NK_PAIRED.csv       within-donor normal->tumour deltas
  results/cervical/NK_SENSITIVITY.csv  gating-threshold scan
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

OUT = "/home/user/aaa/results/cervical"

DATASETS = ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]
TISSUE_ORDER = ["normal_healthy", "normal_adj", "HSIL", "tumor", "lymph_node"]
TISSUE_LABEL = {
    "normal_healthy": "normal cervix (healthy donor)",
    "normal_adj": "normal cervix (tumour-adjacent)",
    "HSIL": "HSIL / CIN",
    "tumor": "cervical cancer",
    "lymph_node": "metastatic lymph node",
}

# a donor-level number is flagged unreliable below these
MIN_CELLS_FOR_TRUST = 500
MIN_TNK_FOR_TRUST = 100


def load_all() -> pd.DataFrame:
    frames = []
    for d in DATASETS:
        p = f"{OUT}/nk_per_sample_{d}.csv"
        if os.path.exists(p):
            frames.append(pd.read_csv(p))
        else:
            print(f"  [missing] {p}")
    if not frames:
        raise SystemExit("no per-sample results found; run build_dataset.py first")
    df = pd.concat(frames, ignore_index=True)
    df["tissue_label"] = df.tissue.map(TISSUE_LABEL)
    df["tissue"] = pd.Categorical(df.tissue, TISSUE_ORDER, ordered=True)
    df["flag"] = ""
    df.loc[df.n_cells_total < MIN_CELLS_FOR_TRUST, "flag"] += "few_cells;"
    df.loc[df.n_TNK < MIN_TNK_FOR_TRUST, "flag"] += "few_TNK;"
    df.loc[df.n_NK < 20, "flag"] += "few_NK;"
    return df.sort_values(["dataset", "tissue", "donor", "sample"])


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (ds, tis), g in df.groupby(["dataset", "tissue"], observed=True):
        for col in ["pct_of_all", "pct_of_immune", "pct_of_lymphocyte"]:
            v = g[col].dropna()
            if not len(v):
                continue
            rows.append({
                "dataset": ds, "tissue": tis, "tissue_label": TISSUE_LABEL[tis],
                "denominator": col, "n_samples": len(v),
                "median_pct": round(float(v.median()), 2),
                "min_pct": round(float(v.min()), 2),
                "max_pct": round(float(v.max()), 2),
                "total_NK_cells": int(g.n_NK.sum()),
            })
    return pd.DataFrame(rows)


def paired(df: pd.DataFrame) -> pd.DataFrame:
    """Within-donor normal vs tumour. Only valid inside one dataset."""
    rows = []
    for ds, g in df.groupby("dataset", observed=True):
        for donor, gd in g.groupby("donor", observed=True):
            tis = set(gd.tissue.astype(str))
            if "tumor" in tis and ("normal_adj" in tis or "normal_healthy" in tis):
                nrm = gd[gd.tissue.astype(str).isin(["normal_adj", "normal_healthy"])].iloc[0]
                tum = gd[gd.tissue.astype(str) == "tumor"].iloc[0]
                for col in ["pct_of_all", "pct_of_immune", "pct_of_lymphocyte"]:
                    rows.append({
                        "dataset": ds, "donor": donor, "denominator": col,
                        "normal_sample": nrm["sample"], "tumor_sample": tum["sample"],
                        "normal_pct": round(nrm[col], 2), "tumor_pct": round(tum[col], 2),
                        "delta_pp": round(tum[col] - nrm[col], 2),
                        "ratio_tumor_over_normal":
                            round(tum[col] / nrm[col], 2) if nrm[col] else np.nan,
                    })
    return pd.DataFrame(rows)


def sensitivity() -> pd.DataFrame:
    """Re-gate the cached T/NK objects across a grid of thresholds."""
    import scanpy as sc
    CACHE = "/home/user/cervical_work/cache"
    grid = [
        ("primary",  0.25, 3, 0.60, 0.40),
        ("strict",   0.15, 4, 0.70, 0.50),
        ("lenient",  0.35, 3, 0.50, 0.30),
        ("3neg_only", 0.25, 3, 0.00, 0.00),
        ("4neg",     0.25, 4, 0.60, 0.40),
    ]
    rows = []
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        if not os.path.exists(p):
            continue
        tnk = sc.read_h5ad(p)
        markers = sorted(set(cc.NK_POS + cc.T_NEG))
        tab = cc.cluster_marker_table(tnk, "tnk_cluster", markers, layer="counts")
        sizes = tnk.obs.groupby(["tnk_cluster", "sample"], observed=True).size().unstack(fill_value=0)
        for name, tneg, minneg, nkg7, gnly in grid:
            n_t_neg = sum((tab[g] < tneg).astype(int) for g in cc.T_NEG if g in tab)
            rec = np.zeros(len(tab), dtype=bool)
            for g, thr in {"KLRD1": 0.40, "KLRF1": 0.20, "NCAM1": 0.10}.items():
                if g in tab:
                    rec |= (tab[g] >= thr).values
            is_nk = ((n_t_neg >= minneg).values
                     & (tab.get("NKG7", 0) >= nkg7).values
                     & (tab.get("GNLY", 0) >= gnly).values
                     & rec)
            nk_clusters = tab.index[is_nk]
            per_sample = sizes.loc[sizes.index.isin(nk_clusters)].sum(axis=0)
            tot = sizes.sum(axis=0)
            for smp in tot.index:
                rows.append({
                    "dataset": ds, "setting": name, "sample": smp,
                    "n_NK": int(per_sample.get(smp, 0)),
                    "n_TNK": int(tot[smp]),
                    "pct_of_TNK": round(100 * per_sample.get(smp, 0) / tot[smp], 2)
                    if tot[smp] else np.nan,
                })
    return pd.DataFrame(rows)


def main():
    df = load_all()
    cols = ["dataset", "donor", "sample", "tissue_label", "histology",
            "n_cells_total", "n_immune", "n_lymphocyte", "n_TNK", "n_NK",
            "pct_of_all", "pct_of_immune", "pct_of_lymphocyte",
            "n_NK_bright", "pct_bright_of_NK",
            "decontx_theta_median", "median_UMI", "flag"]
    tab = df[cols].copy()
    for c in ["pct_of_all", "pct_of_immune", "pct_of_lymphocyte", "pct_bright_of_NK"]:
        tab[c] = tab[c].round(2)
    tab.to_csv(f"{OUT}/NK_TABLE.csv", index=False)
    print("=== PER-DONOR TABLE ===")
    print(tab.to_string(index=False))

    s = summarise(df)
    s.to_csv(f"{OUT}/NK_SUMMARY.csv", index=False)
    print("\n=== SUMMARY BY TISSUE ===")
    print(s.to_string(index=False))

    p = paired(df)
    if len(p):
        p.to_csv(f"{OUT}/NK_PAIRED.csv", index=False)
        print("\n=== WITHIN-DONOR PAIRS ===")
        print(p.to_string(index=False))

    try:
        sens = sensitivity()
        if len(sens):
            sens.to_csv(f"{OUT}/NK_SENSITIVITY.csv", index=False)
            piv = (sens.groupby(["dataset", "setting"])
                       .apply(lambda g: pd.Series({
                           "total_NK": g.n_NK.sum(),
                           "median_pct_of_TNK": round(g.pct_of_TNK.median(), 2)}),
                          include_groups=False)
                       .reset_index())
            print("\n=== GATING SENSITIVITY ===")
            print(piv.to_string(index=False))
    except Exception as e:
        print(f"\n[sensitivity scan failed: {e}]")


if __name__ == "__main__":
    main()
