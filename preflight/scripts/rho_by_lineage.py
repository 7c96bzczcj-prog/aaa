#!/usr/bin/env python3
"""Cross-lineage rho, per library.

The instruction's premise for task A rests on an earlier cross-dataset finding
of ours: NK carries a much larger soup fraction than the other lineages
(rho_NK 0.151 vs rho_CD8T 0.058), and therefore a cross-lineage control cannot
absorb NK-specific contamination.  That premise was established on tumour /
adjacent-normal tissue.  This checks whether it reproduces in marrow and blood,
because if it does not, the "NK is special" argument does not carry over and
should not be quoted as if it did.

Markers are genes each lineage does not express; the estimator is the inherited
one, so these numbers are on the same scale as the 0.151 / 0.058 pair.
"""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

IG = ["IGKC", "IGHG1", "IGHA1"]
CYTO = ["GZMB", "GZMA", "GZMH"]
MARKERS = {"NK": IG, "T": IG, "Myeloid": IG, "B": CYTO, "Erythroid": IG,
           "HSPC": IG, "GranProg": IG}
MIN_CELLS = 30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    idx = json.load(open(a.index))
    cache, rows = {}, []
    for e in idx:
        if e["labels"] not in cache:
            z = np.load(e["labels"], allow_pickle=True)
            cache[e["labels"]] = (z["barcode"], z["library"], z["lineage"])
        bc, libcol, lin = cache[e["labels"]]
        sel_lib = libcol == e["library"]
        mat = (pf.read_10x_h5(e["path"]) if e["kind"] == "h5"
               else pf.read_10x_mtx(*e["path"].split("|"), e["library"]))
        prof, n_empty = pf.soup_profile(mat)
        for lineage, mk in MARKERS.items():
            m = sel_lib & (lin == lineage)
            if m.sum() < MIN_CELLS:
                continue
            sel = np.flatnonzero(np.isin(mat.barcodes, list(set(bc[m]))))
            if len(sel) < MIN_CELLS:
                continue
            obs = np.asarray(mat.mat[:, sel].sum(axis=1)).ravel().astype(float)
            rows.append({"dataset": e["dataset"], "donor": e["donor"],
                         "compartment": e["compartment"], "library": e["library"],
                         "lineage": lineage, "n_cells": len(sel),
                         "rho": pf.estimate_rho(obs, mat.genes, prof, mk),
                         "total_umi": float(obs.sum())})
        del mat
        print(f"  {e['donor']:<12}{e['library'][:26]:<28} done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(a.out, sep="\t", index=False)
    print(df.groupby(["dataset", "compartment", "lineage"])
          .apply(lambda g: pd.Series({
              "n_lib": len(g),
              "rho_mean": float(np.average(g.rho, weights=g.n_cells)),
              "rho_median": float(g.rho.median())}), include_groups=False)
          .to_string())


if __name__ == "__main__":
    main()
