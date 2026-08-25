"""
Cervical arm: effector scores for every gated NK cell.

Reads the existing T/NK objects and their NK flags. Nothing upstream is re-run:
QC, doublet removal, decontX, clustering and the two-stage NK gate are all taken
as they stand.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "cervical"))
import effector_common as ec

CACHE = "/home/user/cervical_work/cache"
OUT = "/home/user/aaa/results/xcancer"
os.makedirs(OUT, exist_ok=True)

DATASETS = ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]

# the four tissue groups asked for. GSE208653's normals are colposcopy biopsies
# from patients, not tumour-adjacent tissue, so they group with healthy rather
# than with the E-MTAB adjacent-normals; the distinction is kept in `tissue_raw`.
TISSUE_GROUP = {
    "normal_healthy": "healthy",
    "normal": "healthy",
    "normal_adj": "adjacent",
    "HSIL": "HSIL",
    "tumor": "cancer",
    "lymph_node": "met_LN",
}


def main():
    parts = []
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        if not os.path.exists(p):
            print(f"[skip] {ds}")
            continue
        a = sc.read_h5ad(p)
        if "is_NK" not in a.obs:
            raise SystemExit(f"{ds}: no is_NK flag; run regate.py first")
        nk = a[a.obs.is_NK.astype(bool)].copy()
        print(f"{ds}: {nk.n_obs} NK cells of {a.n_obs} T/NK")
        if nk.n_obs == 0:
            continue
        obs = nk.obs[["sample", "donor", "tissue", "histology"]].copy()
        obs["dataset"] = ds
        obs["cancer"] = "cervical"
        obs["tissue_raw"] = obs["tissue"]
        obs["tissue_group"] = obs["tissue"].map(TISSUE_GROUP).astype(str)
        tab = ec.effector_table(nk.layers["counts"], nk.var_names, obs)
        print(f"   genes used: {tab.attrs['genes_used']}"
              + (f"  MISSING {tab.attrs['genes_missing']}" if tab.attrs["genes_missing"] else ""))
        parts.append(tab)

    tab = pd.concat(parts, ignore_index=False)
    tab.to_parquet(f"{OUT}/cervical_nk_effector_cells.parquet")
    print(f"\ntotal cervical NK cells scored: {len(tab)}")

    s = ec.summarise(tab)
    s.to_csv(f"{OUT}/cervical_nk_effector_by_donor.csv", index=False)
    print("\n=== per donor x tissue (no pooling) ===")
    print(s.sort_values(["tissue_group", "dataset", "donor"]).to_string(index=False))

    print("\n=== cells per tissue group ===")
    print(tab.groupby("tissue_group", observed=True).size().to_string())


if __name__ == "__main__":
    main()
