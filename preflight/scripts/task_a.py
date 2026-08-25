#!/usr/bin/env python3
"""Task A -- how much of the BM-PB NK difference is compartment ambient RNA.

Per library (an emulsion = a soup):
  1. ambient profile from empty droplets of the UNFILTERED matrix
  2. NK pseudobulk from the cluster-level NK gate
  3. rho_NK  = observed counts of NK-foreign genes / counts expected if the
     whole NK pseudobulk were soup                      (inherited estimator)
  4. per target gene, the share of that gene's NK counts the ambient model
     attributes to soup

Donor value = NK-count-weighted mean over that donor's libraries.
Statistical unit = donor, never library and never cell (spec section 2.3).

Acceptance test: genes NK cannot express (MS4A1, CD79A, LYZ, HBB) must come
back with an ambient share near 1.  If they do not, the model is wrong and no
Delta_ambient from this run may be read.
"""
from __future__ import annotations

import argparse, glob, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

ACCEPTANCE_GENES = ["MS4A1", "CD79A", "LYZ", "HBB"]
ACCEPTANCE_MIN = 0.50           # ambient share of an NK-foreign gene must clear this
EMPTY_WINDOWS = [(1, 20), (1, 100)]   # primary first, then the sensitivity arm


def rows_for(genes, wanted):
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(g, i)
    return [idx.get(g, -1) for g in wanted]


def one_library(mat, labels_bc, labels_lin, lib_name, gate="ratio"):
    out = {}
    keep = labels_lin == "NK"   # doublets already carry the label "doublet"
    bc_nk = set(labels_bc[keep])
    sel = np.flatnonzero(np.isin(mat.barcodes, list(bc_nk)))
    cells = pf.call_cells(mat)
    res = {"library": lib_name, "n_cells_called": len(cells), "n_nk": len(sel)}
    if len(sel) < pf.MIN_NK_CELLS:
        res["skipped"] = f"NK<{pf.MIN_NK_CELLS}"
        return res, None
    obs = np.asarray(mat.mat[:, sel].sum(axis=1)).ravel().astype(float)
    genes = mat.genes
    tgt_rows = rows_for(genes, pf.TARGET_GENES)
    acc_rows = rows_for(genes, ACCEPTANCE_GENES)
    per_window = {}
    for lo, hi in EMPTY_WINDOWS:
        prof, n_empty = pf.soup_profile(mat, lo, hi)
        rho_n = pf.estimate_rho(obs, genes, prof, pf.RHO_MARKERS_NK)
        rho_w = pf.estimate_rho(obs, genes, prof, pf.RHO_MARKERS_NK_WIDE)
        rho_m = pf.rho_median(obs, genes, prof)
        ok, n_test, frac_in, med_r = pf.acceptance(obs, genes, prof, rho_n)
        ok_m, _, frac_in_m, _ = pf.acceptance(obs, genes, prof, rho_m)
        af = pf.ambient_fraction(obs, prof, rho_n, [r for r in tgt_rows])
        afm = pf.ambient_fraction(obs, prof, rho_m, [r for r in tgt_rows])
        acc = pf.ambient_fraction(obs, prof, rho_n, [r for r in acc_rows])
        per_window[(lo, hi)] = {"rho_narrow": rho_n, "rho_wide": rho_w,
                                "rho_median": rho_m, "n_empty": n_empty,
                                "amb": af, "amb_rhomed": afm, "acc": acc,
                                "accept_pass": ok, "accept_n": n_test,
                                "accept_frac_within": frac_in,
                                "accept_median_ratio": med_r,
                                "accept_pass_rhomed": ok_m,
                                "accept_frac_within_rhomed": frac_in_m}
    lo, hi = EMPTY_WINDOWS[0]
    p = per_window[(lo, hi)]
    res.update({"n_empty_droplets": p["n_empty"], "rho_narrow": p["rho_narrow"],
                "rho_wide": p["rho_wide"], "rho_median": p["rho_median"],
                "accept_pass": p["accept_pass"], "accept_n": p["accept_n"],
                "accept_frac_within": p["accept_frac_within"],
                "accept_median_ratio": p["accept_median_ratio"],
                "accept_pass_rhomed": p["accept_pass_rhomed"],
                "accept_frac_within_rhomed": p["accept_frac_within_rhomed"],
                "rho_narrow_w100": per_window[EMPTY_WINDOWS[1]]["rho_narrow"],
                "nk_total_umi": float(obs.sum())})
    for k, g in enumerate(pf.TARGET_GENES):
        r = tgt_rows[k]
        res[f"amb_{g}"] = p["amb"][k]
        res[f"amb100_{g}"] = per_window[EMPTY_WINDOWS[1]]["amb"][k]
        res[f"ambM_{g}"] = p["amb_rhomed"][k]
        res[f"cnt_{g}"] = float(obs[r]) if r >= 0 else np.nan
        if r >= 0:
            sub = mat.mat[r, sel]
            res[f"det_{g}"] = float((np.asarray(sub.todense()).ravel() > 0).mean())
        else:
            res[f"det_{g}"] = np.nan
    for k, g in enumerate(ACCEPTANCE_GENES):
        res[f"acc_{g}"] = p["acc"][k]
    return res, obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--index", required=True,
                    help="JSON: [{donor, compartment, library, path, kind, labels}]")
    ap.add_argument("--gate", default="ratio", choices=["ratio", "cluster", "cd3zero"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    index = json.load(open(args.index))
    key = {"ratio": "lineage", "cluster": "lineage_cluster", "cd3zero": "lineage_cd3zero"}[args.gate]

    cache = {}
    rows = []
    for e in index:
        lab_path = e["labels"]
        if lab_path not in cache:
            z = np.load(lab_path, allow_pickle=True)
            cache[lab_path] = (z["barcode"], z["library"], z[key])
        bc, libcol, lin = cache[lab_path]
        m = (libcol == e["library"])
        if e["kind"] == "h5":
            mat = pf.read_10x_h5(e["path"])
        else:
            mat = pf.read_10x_mtx(*e["path"].split("|"), e["library"])
        r, _ = one_library(mat, bc[m], lin[m], e["library"], args.gate)
        r.update({"dataset": args.dataset, "donor": e["donor"],
                  "compartment": e["compartment"], "gate": args.gate,
                  "group": e.get("group", "")})
        rows.append(r)
        print(f"  {e['donor']:<12} {e['compartment']:<5} {e['library']:<28} "
              f"NK={r['n_nk']:<6} rho={r.get('rho_narrow', float('nan')):.4f}", flush=True)
        del mat
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out} ({len(df)} libraries)")


if __name__ == "__main__":
    main()
