"""
Lung arm: effector scores for every gated NK cell, same code as the cervical arm.

The NK cells come from build_dataset.py + regate.py run on the lung datasets with
no parameter changed, and the score comes from effector_common.effector_table,
the same function the cervical arm calls. Nothing about the scoring differs
between cancer types.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import effector_common as ec

CACHE = "/home/user/cervical_work/cache"
OUT = "/home/user/aaa/results/xcancer"
os.makedirs(OUT, exist_ok=True)

DATASETS = ["GSE131907", "GSE154826"]
TISSUE_GROUP = {"tumor": "cancer", "normal_adj": "adjacent"}


def main():
    parts = []
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        if not os.path.exists(p):
            print(f"[skip] {ds}: not built")
            continue
        a = sc.read_h5ad(p)
        if "is_NK" not in a.obs:
            raise SystemExit(f"{ds}: no is_NK flag; run regate.py")
        nk = a[a.obs.is_NK.astype(bool)].copy()
        print(f"{ds}: {nk.n_obs} NK cells of {a.n_obs} T/NK")
        if nk.n_obs == 0:
            continue
        obs = nk.obs[["sample", "donor", "tissue", "histology"]].copy()
        obs["dataset"] = ds
        obs["cancer"] = "lung"
        obs["tissue_raw"] = obs["tissue"]
        obs["tissue_group"] = obs["tissue"].map(TISSUE_GROUP).astype(str)
        tab = ec.effector_table(nk.layers["counts"], nk.var_names, obs)
        print(f"   genes used: {tab.attrs['genes_used']}"
              + (f"  MISSING {tab.attrs['genes_missing']}" if tab.attrs["genes_missing"] else ""))
        parts.append(tab)

    if not parts:
        raise SystemExit("no lung datasets built yet")
    tab = pd.concat(parts)
    tab.to_parquet(f"{OUT}/lung_nk_effector_cells.parquet")
    print(f"\ntotal lung NK cells scored: {len(tab)}")

    s = ec.summarise(tab)
    s.to_csv(f"{OUT}/lung_nk_effector_by_donor.csv", index=False)
    print("\n=== per donor x tissue (no pooling) ===")
    print(s.sort_values(["tissue_group", "dataset", "donor"]).to_string(index=False))
    print("\n=== cells per tissue group ===")
    print(tab.groupby("tissue_group", observed=True).size().to_string())


if __name__ == "__main__":
    main()
