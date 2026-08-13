#!/usr/bin/env python
"""T1 detection_rates.tsv — the atomic table everything else joins to.

Pipeline steps 5-6:
  5. depth matching. Every cell is thinned to a common UMI floor, and the
     retained fraction is reported per donor x compartment x subset (R6).
     Two alternate floors are computed alongside the primary for sensitivity.
  6. admission. Units below the cell threshold are recorded in
     excluded_units.tsv rather than silently dropped.

Primary readout is DETECTION RATE (R3). mean_cpm rides along and carries no
conclusion.

    python src/detection_rates.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import (as_csr, cell_qc, choose_depth_floor,  # noqa: E402
                            detection_and_cpm, downsample_columns)
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

CEILING = 0.85
FLOOR = 0.05
RETENTION_GAP_TRIGGER = 0.10


def compute_table(ds, panel, hits, keep_mask, depth, min_cells, seed, subsets=None):
    """Detection rates for every donor x compartment x subset unit at one depth."""
    obs = ds.obs
    X = as_csr(ds.X)
    genes = [g for g in panel["gene_symbol"] if g in hits]
    cols = [hits[g] for g in genes]

    sel = keep_mask & obs["subset"].notna().to_numpy()
    if subsets is not None:
        sel = sel & obs["subset"].isin(subsets).to_numpy()

    rows, retention = [], []
    for (donor, comp, ss), rows_idx in unit_indices(
            obs, sel, ["donor", "compartment", "subset"]):
        Xu = X[rows_idx]
        n_pre = Xu.shape[0]
        counts, kept = downsample_columns(Xu, cols, depth, seed=seed)
        n_post = int(kept.sum())
        frac = n_post / n_pre if n_pre else np.nan
        retention.append({"dataset_id": ds.dataset_id, "donor_id": donor,
                          "compartment": comp, "subset": ss, "depth_floor": depth,
                          "n_cells_pre_downsample": n_pre,
                          "n_cells_post_downsample": n_post,
                          "cells_retained_frac": frac})
        n_det, det, cpm = detection_and_cpm(counts, depth)
        for j, g in enumerate(genes):
            if n_post == 0:
                flag = "below_min_cells"
                d = c = np.nan
                nd = 0
            else:
                d, c, nd = float(det[j]), float(cpm[j]), int(n_det[j])
                if n_post < min_cells:
                    flag = "below_min_cells"
                elif d > CEILING:
                    flag = "uninterpretable_ceiling"
                elif d < FLOOR:
                    flag = "uninterpretable_floor"
                else:
                    flag = "ok"
            rows.append({"dataset_id": ds.dataset_id, "donor_id": donor,
                         "compartment": comp, "subset": ss, "gene": g,
                         "n_cells": n_post, "n_detected": nd, "detection_rate": d,
                         "mean_cpm": c, "depth_floor": depth,
                         "cells_retained_frac": frac, "qc_flag": flag})
    return pd.DataFrame(rows), pd.DataFrame(retention)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260813)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    ds = load_dataset(mf)
    hits, misses, rate, _ = ds.panel_columns(panel)
    keep, qc = cell_qc(as_csr(ds.X), ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs

    # depth floor from the NK cells that carry the analysis
    nk = keep & obs["subset"].fillna("").str.startswith(("dNK", "pbNK")).to_numpy()
    nk_umi = qc["total_umi"][nk]
    primary = choose_depth_floor(nk_umi, 0.10)
    alt_lo = choose_depth_floor(nk_umi, 0.05)
    alt_hi = choose_depth_floor(nk_umi, 0.25)
    floors = sorted({primary, alt_lo, alt_hi})
    print(f"[depth] NK UMI: q05={alt_lo} q10={primary} q25={alt_hi} "
          f"(n={nk.sum()} cells)")
    print(f"[depth] primary floor {primary}; sensitivity floors {floors}")

    all_rows, all_ret = [], []
    for d in floors:
        t, r = compute_table(ds, panel, hits, keep, d, args.min_cells, args.seed)
        all_rows.append(t)
        all_ret.append(r)
        print(f"[depth {d}] {len(t)} rows, {r['n_cells_post_downsample'].sum()} cells retained")

    T1 = pd.concat(all_rows, ignore_index=True)
    T1 = T1[["dataset_id", "donor_id", "compartment", "subset", "gene", "n_cells",
             "n_detected", "detection_rate", "mean_cpm", "depth_floor",
             "cells_retained_frac", "qc_flag"]]
    T1.to_csv(os.path.join(outdir, "detection_rates.tsv"), sep="\t", index=False)
    print(f"[T1] {len(T1)} rows -> {outdir}/detection_rates.tsv")

    RET = pd.concat(all_ret, ignore_index=True)
    RET.to_csv(os.path.join(outdir, "retention_by_unit.tsv"), sep="\t", index=False)

    # --- R6: differential cell loss ---------------------------------------
    prim = RET[RET["depth_floor"] == primary]
    nkr = prim[prim["subset"].str.startswith(("dNK", "pbNK"))]
    gaps = []
    for (donor, comp), g in nkr.groupby(["donor_id", "compartment"], observed=True):
        g = g[g["n_cells_pre_downsample"] >= args.min_cells]
        if len(g) < 2:
            continue
        gap = float(g["cells_retained_frac"].max() - g["cells_retained_frac"].min())
        gaps.append({"donor_id": donor, "compartment": comp, "retention_gap": gap,
                     "worst_subset": g.loc[g["cells_retained_frac"].idxmin(), "subset"],
                     "best_subset": g.loc[g["cells_retained_frac"].idxmax(), "subset"],
                     "triggers_sensitivity": gap > RETENTION_GAP_TRIGGER})
    G = pd.DataFrame(gaps)
    G.to_csv(os.path.join(outdir, "retention_gaps.tsv"), sep="\t", index=False)
    triggered = bool(G["triggers_sensitivity"].any()) if len(G) else False
    print(f"\n[R6] across-subset retention gap at the primary floor:")
    if len(G):
        print(G.round(3).to_string(index=False))
    print(f"[R6] sensitivity analysis {'REQUIRED' if triggered else 'not required'} "
          f"(trigger: any gap > {RETENTION_GAP_TRIGGER:.0%}); "
          f"three floors were computed regardless")

    # --- T6 append: units lost to downsampling ----------------------------
    lost = prim[(prim["n_cells_pre_downsample"] >= args.min_cells) &
                (prim["n_cells_post_downsample"] < args.min_cells)].copy()
    exc_path = os.path.join(outdir, "excluded_units.tsv")
    add = pd.DataFrame({
        "dataset_id": lost["dataset_id"], "donor": lost["donor_id"],
        "compartment": lost["compartment"], "subset": lost["subset"],
        "n_cells_post_qc": lost["n_cells_post_downsample"], "stage": "post_downsample",
        "reason": (f"cleared {args.min_cells} cells after QC but fell below it after "
                   f"depth matching to {primary} UMI")})
    if os.path.exists(exc_path):
        prev = pd.read_csv(exc_path, sep="\t")
        pd.concat([prev, add], ignore_index=True).to_csv(exc_path, sep="\t", index=False)
    else:
        add.to_csv(exc_path, sep="\t", index=False)
    print(f"[T6] {len(add)} further units lost to depth matching")

    meta = {"dataset_id": mf.dataset_id, "primary_depth_floor": int(primary),
            "sensitivity_floors": [int(f) for f in floors], "seed": args.seed,
            "min_cells": args.min_cells, "ceiling": CEILING, "floor": FLOOR,
            "panel_match_rate": rate, "panel_unmatched": misses,
            "retention_sensitivity_triggered": triggered}
    with open(os.path.join(outdir, "detection_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
