#!/usr/bin/env python3
"""Task C supplement.

Three things the main pass does not answer:

  * an assumption-free superset bound.  Every NK developmental stage 1-3 is
    CD34+ or CD117(KIT)+, so the whole CD34/KIT+ compartment is an upper bound
    that assumes nothing about NK priming genes at all.  A ladder of
    progressively tighter gates is reported so the reader can see which
    assumption buys which bound.

  * the literature's own QC boundary.  Melsen 2022 removed every cell with
    fewer than 1,000 expressed genes.  How many cells -- and how many gated
    cells -- would that threshold delete here?

  * gate sensitivity.  A bound built on "IL2RB >= 1 count" is only as good as
    IL2RB's detection rate.  Detection is measured inside populations that
    certainly carry the transcript at protein level.
"""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "scripts"))
from task_c import GATES, apply_gate, _rows, PRE_QC_MIN_UMI, POST_QC_MIN_UMI, POST_QC_MIN_GENES  # noqa: E402

LIT_QC_MIN_GENES = 1000     # Melsen 2022's own threshold
EXTRA_GATES = {
    "prog_anchor_any": dict(any_of=[["CD34", "KIT"]], zero=[]),
    "prog_anchor_noCyto": dict(any_of=[["CD34", "KIT"]],
                               zero=["PRF1", "GZMB", "GNLY", "NKG7"]),
}
DETECT_IN_NK = ["IL2RB", "IL7R", "KIT", "CD34", "NCAM1", "KLRF1"]
DETECT_IN_PROG = ["IL2RB", "IL7R", "KIT", "CD34", "DNTT", "GATA1", "TOX", "ID2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    idx = json.load(open(a.index))
    rows = []
    allg = {**GATES, **EXTRA_GATES}
    for e in idx:
        mat = (pf.read_10x_h5(e["path"]) if e["kind"] == "h5"
               else pf.read_10x_mtx(*e["path"].split("|"), e["library"]))
        g2i = _rows(mat.genes)
        umi, ng = mat.umi, mat.n_genes_det
        post = np.flatnonzero((umi >= POST_QC_MIN_UMI) & (ng >= POST_QC_MIN_GENES))
        r = {"dataset": e["dataset"], "donor": e["donor"], "library": e["library"],
             "compartment": e["compartment"], "n_postQC": len(post),
             "n_postQC_ge1000genes": int((ng[post] >= LIT_QC_MIN_GENES).sum())}
        sub = mat.mat[:, post]
        ng_post = ng[post]
        hits = {}
        for gname, spec in allg.items():
            h = apply_gate(sub, g2i, spec)
            hits[gname] = h
            k = int(h.sum()) if h is not None else -1
            r[f"{gname}_n"] = k
            r[f"{gname}_n_ge1000genes"] = (int((h & (ng_post >= LIT_QC_MIN_GENES)).sum())
                                           if h is not None else -1)
            r[f"{gname}_median_genes"] = (float(np.median(ng_post[h])) if h is not None
                                          and h.any() else np.nan)
        for pop, genes in (("NK", DETECT_IN_NK), ("PROG", DETECT_IN_PROG)):
            h = hits["matureNK_marker"] if pop == "NK" else hits["prog_anchor_any"]
            if h is None or not h.any():
                continue
            s = sub[:, np.flatnonzero(h)]
            for g in genes:
                rr = g2i.get(g)
                if not rr:
                    continue
                v = np.asarray(s[rr].sum(axis=0)).ravel()
                r[f"det_{pop}_{g}"] = float((v > 0).mean())
                r[f"n_{pop}"] = int(h.sum())
        rows.append(r)
        print(f"  {e['donor']:<11}{e['library'][:22]:<24} post={len(post):>6} "
              f"ge1000g={r['n_postQC_ge1000genes']:>6} prog={r['prog_anchor_any_n']:>5} "
              f"CLP={r['CLP_like_n']:>4} eNK={r['earlyNK_permissive_n']:>3} "
              f"detNK_IL2RB={r.get('det_NK_IL2RB', float('nan')):.3f}", flush=True)
        del mat
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
