#!/usr/bin/env python
"""T5 qc_report.tsv + T6 excluded_units.tsv, plus the NK-gate purity check.

Pipeline steps 2-4 of the specification:
  2. cell QC (n_genes < 200 or mito >= 10% removed), counted per
     donor x compartment x celltype before and after
  3. annotation verification against subset markers -- REPORTED, never repaired
  4. NK-gate purity: TRBC2+CD3E+CD3D+ triple-positive fraction per subset,
     with a non-T lineage carried alongside as the ambient floor

    python src/qc_report.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr, cell_qc, triple_positive_fraction  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

T_GATE = ["TRBC2", "CD3E", "CD3D"]
PURITY_STOP = 0.15


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--min-cells", type=int, default=30)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    ds = load_dataset(mf)
    hits, misses, rate, _ = ds.panel_columns(panel)
    obs = ds.obs
    X = as_csr(ds.X)

    keep, qc = cell_qc(X, ds.symbols, min_genes=args.min_genes,
                       max_mito_frac=args.max_mito)
    obs = obs.assign(qc_pass=keep, n_genes=qc["n_genes"],
                     total_umi=qc["total_umi"], mito_frac=qc["mito_frac"])
    print(f"[QC] mito genes found in matrix: {qc['n_mito_genes_in_matrix']}")
    print(f"[QC] cells passing: {keep.sum()} / {len(keep)} ({keep.mean():.1%})")

    # --- T5 -------------------------------------------------------------
    grp = ["donor", "compartment", "celltype_raw", "subset"]
    rows = (obs.groupby(grp, observed=True, dropna=False)
            .agg(n_cells_pre_qc=("qc_pass", "size"),
                 n_cells_post_qc=("qc_pass", "sum"),
                 median_umi=("total_umi", "median"),
                 median_genes=("n_genes", "median"),
                 median_mito_frac=("mito_frac", "median"))
            .reset_index())
    rows["frac_cells_kept_qc"] = rows["n_cells_post_qc"] / rows["n_cells_pre_qc"]
    rows["dataset_id"] = mf.dataset_id
    rows["min_genes_threshold"] = args.min_genes
    rows["max_mito_threshold"] = args.max_mito
    rows = rows[["dataset_id"] + grp + ["n_cells_pre_qc", "n_cells_post_qc",
                                        "frac_cells_kept_qc", "median_umi",
                                        "median_genes", "median_mito_frac",
                                        "min_genes_threshold", "max_mito_threshold"]]
    rows.to_csv(os.path.join(outdir, "qc_report.tsv"), sep="\t", index=False)
    print(f"[T5] {len(rows)} rows -> {outdir}/qc_report.tsv")

    # --- T6 -------------------------------------------------------------
    post = obs[obs["qc_pass"]]
    unit = (post[post["subset"].notna()]
            .groupby(["donor", "compartment", "subset"], observed=True)
            .size().rename("n_cells_post_qc").reset_index())
    excl = unit[unit["n_cells_post_qc"] < args.min_cells].copy()
    excl["dataset_id"] = mf.dataset_id
    excl["stage"] = "post_cell_qc"
    excl["reason"] = (f"n_cells < {args.min_cells} after cell QC; "
                      f"below the admission threshold for donor-level statistics")
    excl = excl[["dataset_id", "donor", "compartment", "subset",
                 "n_cells_post_qc", "stage", "reason"]]
    excl.to_csv(os.path.join(outdir, "excluded_units.tsv"), sep="\t", index=False)
    print(f"[T6] {len(excl)} units excluded at this stage -> {outdir}/excluded_units.tsv")

    # --- step 3: annotation verification (report only) --------------------
    marker_rows = []
    mk = panel[panel["category"] == "marker_subset"]["gene_symbol"].tolist()
    mk_cols = [hits[g] for g in mk if g in hits]
    mk_names = [g for g in mk if g in hits]
    mapped_keep = keep & obs["subset"].notna().to_numpy()
    if mk_cols:
        det = (X[:, mk_cols] >= 1).astype(np.int8)
        for (comp, ss), rows_idx in unit_indices(obs, mapped_keep,
                                                 ["compartment", "subset"]):
            d = np.asarray(det[rows_idx].sum(axis=0)).ravel() / max(1, len(rows_idx))
            for g, v in zip(mk_names, d):
                marker_rows.append({"dataset_id": mf.dataset_id, "compartment": comp,
                                    "subset": ss, "gene": g, "detection_rate": float(v),
                                    "n_cells": len(rows_idx)})
    mdf = pd.DataFrame(marker_rows)
    mdf.to_csv(os.path.join(outdir, "annotation_check.tsv"), sep="\t", index=False)
    if len(mdf):
        dnk = mdf[mdf["subset"].isin(["dNK1", "dNK2", "dNK3", "dNKp"])]
        if len(dnk):
            print("\n[step 3] subset marker detection (annotation is CHECKED, never repaired):")
            print(dnk.pivot_table(index="gene", columns="subset",
                                  values="detection_rate").round(3).to_string())

    # --- step 4: NK gate purity ------------------------------------------
    gate_cols = [hits[g] for g in T_GATE if g in hits]
    missing_gate = [g for g in T_GATE if g not in hits]
    if missing_gate:
        raise SystemExit(f"purity gate genes missing from matrix: {missing_gate}")
    pur_rows = []
    for (comp, ss), rows_idx in unit_indices(obs, mapped_keep, ["compartment", "subset"]):
        frac, nhit = triple_positive_fraction(X[rows_idx], gate_cols)
        pur_rows.append({"dataset_id": mf.dataset_id, "compartment": comp, "subset": ss,
                         "n_cells": len(rows_idx), "n_triple_positive": nhit,
                         "triple_positive_frac": frac})
    pdf = pd.DataFrame(pur_rows).sort_values(["compartment", "subset"])
    pdf.to_csv(os.path.join(outdir, "nk_gate_purity.tsv"), sep="\t", index=False)

    print("\n[step 4] TRBC2+CD3E+CD3D+ triple-positive fraction "
          "(single-gene calls are barred: TRBC1 alone has previously scored 95.6% of NK)")
    print(pdf[pdf["n_cells"] >= 30].round(4).to_string(index=False))

    dnk_pur = pdf[pdf["subset"].isin(["dNK1", "dNK2", "dNK3", "dNKp"]) & (pdf["n_cells"] >= 30)]
    worst = dnk_pur["triple_positive_frac"].max() if len(dnk_pur) else 0.0
    stop = bool(worst > PURITY_STOP)
    summary = {"dataset_id": mf.dataset_id,
               "worst_dnk_triple_positive_frac": float(worst),
               "purity_stop_threshold": PURITY_STOP,
               "purity_stop_fires": stop,
               "panel_match_rate": rate,
               "panel_unmatched": misses,
               "cells_pre_qc": int(len(keep)), "cells_post_qc": int(keep.sum())}
    with open(os.path.join(outdir, "qc_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\nworst dNK triple-positive fraction: {worst:.3%} "
          f"(stop threshold {PURITY_STOP:.0%})")
    if stop:
        print("STOP CONDITION FIRES — resolve the annotation before the main analysis.")
        return 3
    print("purity stop does not fire.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
