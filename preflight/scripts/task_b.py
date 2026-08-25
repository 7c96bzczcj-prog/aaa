#!/usr/bin/env python3
"""Task B -- does the ltNK fraction move with age.

Section 5 of the spec puts a gate in front of the question: unless CD69 AND
CXCR6 both detect inside the 5%-85% window in NK cells of that dataset and
compartment, the answer is "unmeasurable" and nothing further is computed.
The gate is checked and reported whatever it says, because a failed gate is
itself the deliverable ("ltNK cannot be defined from transcriptome; it needs
flow").

Definitions, frozen before the run.  All are transcript-presence gates, so no
threshold is tuned:
    ltNK        CD69 >= 1 count AND CXCR6 >= 1 count
    CD56bright  FCGR3A == 0 AND NCAM1 >= 1
    CD56dim     FCGR3A >= 1
All three are reported together, because a change in ltNK that is really a
change in the CD56bright/dim balance is not an ltNK result.

Statistical unit = donor.  Age is continuous and never binned.  Dataset enters
as a covariate.
"""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

GATE_GENES = ["CD69", "CXCR6"]
SUBSET_GENES = ["CD69", "CXCR6", "NCAM1", "FCGR3A", "SELL", "GZMK", "KLRC1"]


def counts_for(mat, cells, gene):
    g2i = {}
    for i, g in enumerate(mat.genes):
        g2i.setdefault(g, []).append(i)
    rows = g2i.get(gene)
    if not rows:
        return None
    return np.asarray(mat.mat[rows][:, cells].sum(axis=0)).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gate", default="ratio")
    args = ap.parse_args()
    key = {"ratio": "lineage", "cluster": "lineage_cluster",
           "cd3zero": "lineage_cd3zero"}[args.gate]
    entries = []
    for f in args.index:
        entries += json.load(open(f))
    cache, rows = {}, []
    for e in entries:
        if e["labels"] not in cache:
            z = np.load(e["labels"], allow_pickle=True)
            cache[e["labels"]] = (z["barcode"], z["library"], z[key])
        bc, libcol, lin = cache[e["labels"]]
        m = (libcol == e["library"]) & (lin == "NK")
        nk_bc = set(bc[m])
        mat = (pf.read_10x_h5(e["path"]) if e["kind"] == "h5"
               else pf.read_10x_mtx(*e["path"].split("|"), e["library"]))
        sel = np.flatnonzero(np.isin(mat.barcodes, list(nk_bc)))
        r = {"dataset": e["dataset"], "donor": e["donor"], "library": e["library"],
             "compartment": e["compartment"], "age": e.get("age", np.nan),
             "n_nk": len(sel)}
        if len(sel) >= pf.MIN_NK_CELLS:
            c = {g: counts_for(mat, sel, g) for g in SUBSET_GENES}
            for g in SUBSET_GENES:
                r[f"det_{g}"] = float((c[g] > 0).mean()) if c[g] is not None else np.nan
            cd69 = c.get("CD69")
            cxcr6 = c.get("CXCR6")
            ncam1 = c.get("NCAM1")
            fcgr3a = c.get("FCGR3A")
            if cd69 is not None and cxcr6 is not None:
                r["n_ltNK"] = int(((cd69 > 0) & (cxcr6 > 0)).sum())
                r["frac_ltNK"] = r["n_ltNK"] / len(sel)
            if ncam1 is not None and fcgr3a is not None:
                bright = (fcgr3a == 0) & (ncam1 > 0)
                dim = fcgr3a > 0
                r["n_bright"] = int(bright.sum())
                r["frac_bright"] = float(bright.mean())
                r["n_dim"] = int(dim.sum())
                r["frac_dim"] = float(dim.mean())
        rows.append(r)
        print(f"  {e['dataset']:<11} {e['donor']:<12} {e['library']:<24} "
              f"NK={r['n_nk']:<6} detCD69={r.get('det_CD69', float('nan')):.3f} "
              f"detCXCR6={r.get('det_CXCR6', float('nan')):.3f}", flush=True)
        del mat
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
