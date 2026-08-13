#!/usr/bin/env python
"""Diagnose a fired NK-gate purity stop.

The stop in docs/PREREGISTRATION.md section 6.1 says: if triple-positive
TRBC2+CD3E+CD3D+ exceeds 15% in any dNK subset, stop and resolve the
annotation before the main analysis. This script does the resolving. It
answers, in order:

  1. is the raw triple-positive fraction confounded by sequencing depth?
     Detecting three genes at >= 1 count each is strongly depth-dependent,
     and the raw metric compares subsets at whatever depth they happen to
     have.
  2. how strong is the T signal in the flagged cells? One count of each gene
     is what ambient pickup looks like; a real T cell carries many.
  3. do the flagged cells look like T cells, or like NK cells with pickup, or
     like NK-T doublets? Judged by whether NK identity genes are retained.
  4. is it concentrated in particular donors or libraries?

It reports. It does not re-label anything (barred by section 7).

    python src/purity_diagnosis.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr, cell_qc, downsample_columns  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

T_GATE = ["TRBC2", "CD3E", "CD3D"]
NK_ID = ["NKG7", "KLRD1", "GNLY", "PRF1"]
FOCUS = ["dNK1", "dNK2", "dNK3", "dNKp", "T", "Myeloid", "Stromal", "cDC1", "Plasma"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--compartment", default="Decidua")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=20260813)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    ds = load_dataset(mf)
    hits, _, _, _ = ds.panel_columns(panel)
    X = as_csr(ds.X)
    keep, qc = cell_qc(X, ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs
    total = qc["total_umi"]

    gate = [hits[g] for g in T_GATE]
    nkid = [hits[g] for g in NK_ID if g in hits]
    comp = (obs["compartment"] == args.compartment).to_numpy()

    # --- 1. depth by subset, and the depth-matched triple-positive rate ----
    sel = keep & comp & obs["subset"].isin(FOCUS).to_numpy()
    depths = []
    for (ss,), idx in unit_indices(obs, sel, ["subset"]):
        depths.append({"subset": ss, "n_cells": len(idx),
                       "median_umi": float(np.median(total[idx])),
                       "q10_umi": float(np.quantile(total[idx], 0.10)),
                       "q90_umi": float(np.quantile(total[idx], 0.90))})
    D = pd.DataFrame(depths).sort_values("median_umi")
    print(f"[1] sequencing depth by subset in {args.compartment}")
    print(D.to_string(index=False))

    floor = int(np.floor(np.quantile(total[sel], 0.10)))
    print(f"\n[1] depth-matched triple-positive rate at a common {floor} UMI floor")
    rows = []
    for (ss,), idx in unit_indices(obs, sel, ["subset"]):
        raw_pos = np.asarray((X[idx][:, gate] >= 1).sum(axis=1)).ravel() == 3
        counts, kept = downsample_columns(X[idx], gate, floor, seed=args.seed)
        matched = (counts >= 1).all(axis=1) if counts.shape[0] else np.array([], bool)
        rows.append({"subset": ss, "n_cells": len(idx),
                     "raw_triple_pos": float(raw_pos.mean()),
                     "n_after_depth_match": int(kept.sum()),
                     "matched_triple_pos": float(matched.mean()) if matched.size else np.nan,
                     "median_umi": float(np.median(total[idx]))})
    P = pd.DataFrame(rows).sort_values("raw_triple_pos", ascending=False)
    P["dataset_id"] = mf.dataset_id
    P["compartment"] = args.compartment
    P["depth_floor"] = floor
    print(P[["subset", "n_cells", "median_umi", "raw_triple_pos",
             "n_after_depth_match", "matched_triple_pos"]].round(4).to_string(index=False))
    P.to_csv(os.path.join(outdir, "purity_depth_matched.tsv"), sep="\t", index=False)

    # --- 2/3. what do the flagged cells look like? -------------------------
    print(f"\n[2/3] character of the triple-positive cells in {args.compartment}")
    ch = []
    for (ss,), idx in unit_indices(obs, sel, ["subset"]):
        g = np.asarray(X[idx][:, gate].todense())
        pos = (g >= 1).all(axis=1)
        if pos.sum() == 0:
            continue
        nk = np.asarray(X[idx][:, nkid].todense())
        tsum = g[pos].sum(axis=1)
        ch.append({
            "subset": ss, "n_triple_pos": int(pos.sum()),
            "median_T_counts_in_flagged": float(np.median(tsum)),
            "frac_flagged_with_T_counts_eq3": float((tsum == 3).mean()),
            "frac_flagged_with_T_counts_ge6": float((tsum >= 6).mean()),
            "median_NK_counts_in_flagged": float(np.median(nk[pos].sum(axis=1))),
            "median_NK_counts_in_unflagged": float(np.median(nk[~pos].sum(axis=1)))
                if (~pos).sum() else np.nan,
            "median_umi_flagged": float(np.median(total[idx][pos])),
            "median_umi_unflagged": float(np.median(total[idx][~pos]))
                if (~pos).sum() else np.nan})
    C = pd.DataFrame(ch).sort_values("n_triple_pos", ascending=False)
    C["dataset_id"] = mf.dataset_id
    C["compartment"] = args.compartment
    print(C.round(3).to_string(index=False))
    C.to_csv(os.path.join(outdir, "purity_character.tsv"), sep="\t", index=False)

    # --- 4. donor / library concentration ----------------------------------
    print(f"\n[4] triple-positive rate by donor, dNK subsets ({args.compartment})")
    dsel = keep & comp & obs["subset"].isin(["dNK1", "dNK2", "dNK3", "dNKp"]).to_numpy()
    dr = []
    for (donor, ss), idx in unit_indices(obs, dsel, ["donor", "subset"]):
        if len(idx) < 30:
            continue
        g = np.asarray(X[idx][:, gate].todense())
        dr.append({"donor": donor, "subset": ss, "n_cells": len(idx),
                   "triple_pos": float((g >= 1).all(axis=1).mean()),
                   "median_umi": float(np.median(total[idx]))})
    R = pd.DataFrame(dr)
    if len(R):
        print(R.pivot_table(index="subset", columns="donor",
                            values="triple_pos").round(3).to_string())
        print("\n     median UMI, same cells")
        print(R.pivot_table(index="subset", columns="donor",
                            values="median_umi").round(0).to_string())
    R["dataset_id"] = mf.dataset_id
    R.to_csv(os.path.join(outdir, "purity_by_donor.tsv"), sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
