#!/usr/bin/env python3
"""Task C -- is the early NK developmental series really absent from marrow,
or is it QC attrition.

Three design decisions, all frozen before any early-NK count was read:

1.  Every count is made TWICE, before QC (any droplet with >= PRE_QC_MIN_UMI
    counts, no gene and no mitochondrial filter) and after (the pipeline's own
    UMI/gene thresholds).  The difference is reported, because that difference
    IS the question the spec asks.

2.  The early-NK gate is deliberately PERMISSIVE.  The deliverable is an
    exclusion bound, and an over-inclusive gate makes the bound conservative:
    if even a loose definition finds almost nothing, the bound means something.
    A strict variant is reported alongside so the gap is visible.

3.  Two developmental series marrow certainly contains -- pro-B and erythroid
    progenitors -- are gated with STRUCTURALLY IDENTICAL definitions
    (progenitor anchor + lineage-priming gene + no cytotoxic programme) and
    differ from the NK gate in exactly one term: which priming gene is
    required.  If the controls come back present and NK comes back empty, the
    empty result is about NK.  If the controls come back empty too, the
    instrument is blind and no conclusion may be drawn.

A first version of these gates required zero counts across a long exclusion
list (LYZ, HBB, MS4A1 among them).  It returned 1-9 mature NK cells in CD45+
marrow libraries that clustering scores in the hundreds, because those genes
are abundant in the soup and a zero-count requirement on them is a filter on
sequencing depth, not on identity.  Ambient-abundant genes were removed from
every exclusion list for that reason.
"""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

PRE_QC_MIN_UMI = 100
POST_QC_MIN_UMI = pf.CELL_MIN_UMI
POST_QC_MIN_GENES = pf.CELL_MIN_GENES

PROG_ANCHOR = ["CD34", "KIT"]
NO_CYTOTOXIC = ["PRF1", "GZMB", "GNLY", "NKG7"]
STRICT_EXTRA = ["DNTT", "VPREB1", "IGLL1", "MPO", "ELANE", "CD3D", "CD3E", "TRAC"]

GATES = {
    "earlyNK_permissive": dict(any_of=[PROG_ANCHOR, ["IL2RB"]], zero=NO_CYTOTOXIC),
    "earlyNK_strict": dict(any_of=[PROG_ANCHOR, ["IL2RB"]],
                           zero=NO_CYTOTOXIC + STRICT_EXTRA),
    "earlyNK_TFvariant": dict(any_of=[PROG_ANCHOR, ["IL2RB", "IL7R"],
                                      ["TOX", "ID2", "ZBTB16", "GATA3"]],
                              zero=NO_CYTOTOXIC),
    "proB_control": dict(any_of=[PROG_ANCHOR, ["DNTT", "VPREB1"]], zero=NO_CYTOTOXIC),
    "eryProg_control": dict(any_of=[PROG_ANCHOR, ["GATA1", "KLF1"]], zero=NO_CYTOTOXIC),
    "CLP_like": dict(any_of=[PROG_ANCHOR, ["IL7R"]], zero=NO_CYTOTOXIC),
    "CD34_any": dict(any_of=[["CD34"]], zero=[]),
    "matureNK_marker": dict(any_of=[["KLRF1", "NKG7"], ["PRF1", "GNLY", "GZMB"]],
                            zero=["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"]),
}


def _rows(genes):
    d = {}
    for i, g in enumerate(genes):
        d.setdefault(g, []).append(i)
    return d


def apply_gate(sub, g2i, spec):
    n = sub.shape[1]
    ok = np.ones(n, dtype=bool)
    for grp in spec["any_of"]:
        rows = [i for g in grp for i in g2i.get(g, [])]
        if not rows:
            return None
        ok &= np.asarray(sub[rows].sum(axis=0)).ravel() > 0
    if spec["zero"]:
        rows = [i for g in spec["zero"] for i in g2i.get(g, [])]
        if rows:
            ok &= np.asarray(sub[rows].sum(axis=0)).ravel() == 0
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    index = json.load(open(args.index))
    rows = []
    for e in index:
        mat = (pf.read_10x_h5(e["path"]) if e["kind"] == "h5"
               else pf.read_10x_mtx(*e["path"].split("|"), e["library"]))
        g2i = _rows(mat.genes)
        umi, ng = mat.umi, mat.n_genes_det
        pre = np.flatnonzero(umi >= PRE_QC_MIN_UMI)
        post = np.flatnonzero((umi >= POST_QC_MIN_UMI) & (ng >= POST_QC_MIN_GENES))
        r = {"dataset": e["dataset"], "donor": e["donor"], "library": e["library"],
             "compartment": e["compartment"], "n_preQC": len(pre),
             "n_postQC": len(post), "n_lost_to_QC": len(pre) - len(post),
             "median_genes_postQC": float(np.median(ng[post])) if len(post) else np.nan,
             "frac_postQC_under_1000_genes":
                 float((ng[post] < 1000).mean()) if len(post) else np.nan}
        for tag, cols in (("preQC", pre), ("postQC", post)):
            sub = mat.mat[:, cols]
            for gname, spec in GATES.items():
                hit = apply_gate(sub, g2i, spec)
                k = int(hit.sum()) if hit is not None else -1
                r[f"{gname}_{tag}_n"] = k
                r[f"{gname}_{tag}_frac"] = (k / len(cols)) if (hit is not None and len(cols)) else np.nan
        rows.append(r)
        print(f"  {e['dataset'][:9]:<10}{e['donor']:<11}{e['library'][:24]:<25} "
              f"pre={len(pre):>7} post={len(post):>6} | "
              f"eNKperm={r['earlyNK_permissive_postQC_n']:>4} "
              f"eNKstr={r['earlyNK_strict_postQC_n']:>4} "
              f"proB={r['proB_control_postQC_n']:>5} ery={r['eryProg_control_postQC_n']:>5} "
              f"CLP={r['CLP_like_postQC_n']:>4} CD34={r['CD34_any_postQC_n']:>5} "
              f"mNK={r['matureNK_marker_postQC_n']:>5}", flush=True)
        del mat
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
