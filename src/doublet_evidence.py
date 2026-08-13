#!/usr/bin/env python
"""Discriminate NK-T doublets from depth-driven ambient pickup.

`purity_diagnosis.py` showed the triple-positive cells are not mislabelled
T cells (they carry MORE NK signal, not less). That rules out one explanation
but does not separate the remaining two, because BOTH produce "flagged cells
are deeper and carry both programmes":

  doublets       -> two cells in one droplet: total UMI roughly DOUBLED, and
                    the NK:T ratio is BIMODAL (a genuine NK mode plus a
                    doublet mode at an intermediate ratio)
  ambient pickup -> one cell that sampled more free RNA: total UMI varies
                    CONTINUOUSLY, and the NK:T ratio is UNIMODAL with a tail

So the two discriminating measurements are:
  1. the total-UMI ratio flagged/unflagged -- near 2.0 argues doublets,
     near 1.0-1.3 argues depth
  2. the shape of the NK:T count-ratio distribution -- bimodal vs unimodal,
     tested with a dip statistic and reported as a histogram

Produces out/<id>/doublet_evidence.tsv and doublet_evidence.png.

    python src/doublet_evidence.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr, cell_qc  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

T_GATE = ["TRBC2", "CD3E", "CD3D"]
NK_ID = ["NKG7", "KLRD1", "GNLY", "PRF1"]
SUBSETS = ["dNK1", "dNK2", "dNK3", "dNKp"]


def dip_like_bimodality(x, nbins=40):
    """Cheap bimodality evidence: Hartigan-style dip is overkill here.

    Returns (bimodality_coefficient, n_modes_in_smoothed_histogram).
    BC = (skew^2 + 1) / kurtosis; BC > 5/9 = 0.555 suggests bimodality.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 30:
        return np.nan, 0
    n = x.size
    m, s = x.mean(), x.std(ddof=1)
    if s == 0:
        return np.nan, 0
    z = (x - m) / s
    skew = float((z ** 3).mean())
    kurt = float((z ** 4).mean())
    bc = (skew ** 2 + 1.0) / kurt if kurt > 0 else np.nan
    hist, edges = np.histogram(x, bins=nbins)
    k = np.ones(3) / 3.0
    sm = np.convolve(hist.astype(float), k, mode="same")
    # NOTE: this peak count is noise-sensitive on small samples (66 flagged
    # dNK3 cells produce spurious peaks). It is reported for completeness; the
    # bimodality coefficient is the statistic to read.
    modes = int(sum(1 for i in range(1, len(sm) - 1)
                    if sm[i] > sm[i - 1] and sm[i] >= sm[i + 1] and sm[i] > 0.05 * sm.max()))
    return bc, modes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--compartment", default="Decidua")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
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
    sel = keep & comp & obs["subset"].isin(SUBSETS + ["T"]).to_numpy()

    rows, panels = [], {}
    for (ss,), idx in unit_indices(obs, sel, ["subset"]):
        g = np.asarray(X[idx][:, gate].todense()).sum(axis=1)
        nk = np.asarray(X[idx][:, nkid].todense()).sum(axis=1)
        pos = np.asarray((X[idx][:, gate] >= 1).todense()).all(axis=1)
        if pos.sum() < 20:
            continue
        u = total[idx]

        # 1. UMI ratio flagged / unflagged
        umi_ratio = float(np.median(u[pos]) / np.median(u[~pos])) if (~pos).sum() else np.nan

        # 2. NK:T ratio shape, on flagged cells (log scale, T floored at 1)
        lr = np.log2((nk[pos] + 1) / (g[pos] + 1))
        bc, modes = dip_like_bimodality(lr)
        lr_all = np.log2((nk + 1) / (g + 1))
        bc_all, modes_all = dip_like_bimodality(lr_all)

        rows.append({
            "dataset_id": mf.dataset_id, "compartment": args.compartment,
            "subset": ss, "n_cells": int(len(idx)), "n_flagged": int(pos.sum()),
            "flagged_frac": float(pos.mean()),
            "median_umi_flagged": float(np.median(u[pos])),
            "median_umi_unflagged": float(np.median(u[~pos])) if (~pos).sum() else np.nan,
            "umi_ratio_flagged_over_unflagged": umi_ratio,
            "nkT_log2ratio_bimodality_coef_flagged": bc,
            "nkT_log2ratio_modes_flagged": modes,
            "nkT_log2ratio_bimodality_coef_all": bc_all,
            "nkT_log2ratio_modes_all": modes_all})
        panels[ss] = (u, pos, lr_all)

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(outdir, "doublet_evidence.tsv"), sep="\t", index=False)
    print("=== discriminating measurements ===")
    print("doublets predict umi_ratio ~2.0 and a bimodal NK:T ratio (BC > 0.555, 2 modes)")
    print("depth-driven pickup predicts umi_ratio ~1.0-1.3 and a unimodal ratio\n")
    print(R[["subset", "n_cells", "n_flagged", "flagged_frac",
             "median_umi_flagged", "median_umi_unflagged",
             "umi_ratio_flagged_over_unflagged",
             "nkT_log2ratio_bimodality_coef_flagged",
             "nkT_log2ratio_modes_flagged"]].round(3).to_string(index=False))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib unavailable; table written, figure skipped)")
        return 0

    ks = [k for k in SUBSETS if k in panels]
    fig, axes = plt.subplots(2, len(ks), figsize=(3.4 * len(ks), 6.4), squeeze=False)
    for j, ss in enumerate(ks):
        u, pos, lr_all = panels[ss]
        ax = axes[0][j]
        bins = np.logspace(np.log10(max(1, u.min())), np.log10(u.max()), 45)
        ax.hist(u[~pos], bins=bins, alpha=.65, label="unflagged", color="#4C78A8")
        ax.hist(u[pos], bins=bins, alpha=.65, label="triple-pos", color="#E45756")
        ax.set_xscale("log")
        r = np.median(u[pos]) / np.median(u[~pos])
        verdict = "doublets" if r > 1.8 else "depth"
        ax.set_title(f"{ss}\nUMI ratio {r:.2f}x  (doublets predict 2.0) -> {verdict}",
                     fontsize=9)
        ax.set_xlabel("total UMI"); ax.set_ylabel("cells")
        if j == 0:
            ax.legend(fontsize=7)
        ax2 = axes[1][j]
        ax2.hist(lr_all[~pos], bins=45, alpha=.65, color="#4C78A8", label="unflagged")
        ax2.hist(lr_all[pos], bins=45, alpha=.65, color="#E45756", label="triple-pos")
        bc, _ = dip_like_bimodality(lr_all[pos])
        shape = "BIMODAL -> doublets" if bc > 0.555 else "unimodal -> depth"
        ax2.set_title(f"NK:T log2 ratio, flagged\nBC={bc:.3f} (>0.555 = bimodal) "
                      f"-> {shape}", fontsize=9)
        ax2.set_xlabel("log2 (NK counts+1)/(T counts+1)"); ax2.set_ylabel("cells")
    fig.suptitle(f"{mf.dataset_id} {args.compartment}: doublets or depth-driven pickup?",
                 fontsize=11)
    fig.tight_layout()
    out = os.path.join(outdir, "doublet_evidence.png")
    fig.savefig(out, dpi=140)
    print(f"\nfigure -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
