#!/usr/bin/env python3
"""Task C aggregation and the pre-written verdict lookup (spec section 4)."""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

VERDICT_TABLE = {
    "C_P3_HOLDS": "no intermediate state before or after QC, and X small enough "
                  "(X <= 0.1%): P3 holds; the marrow is not where NK output is "
                  "read, and the marrow arm should ask about residency, not "
                  "development.",
    "C_QC_ARTEFACT": "present before QC and gone after: P3 is a QC false "
                     "negative and the literature claim is drawn too wide. A "
                     "publishable methodological correction.",
    "C_UNMEASURABLE": "no X can be stated: unmeasurable, and it may not be "
                      "written up as either answer.",
}
GATES = ["earlyNK_permissive", "earlyNK_strict", "earlyNK_TFvariant",
         "proB_control", "eryProg_control", "CLP_like", "CD34_any",
         "matureNK_marker"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--libs", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.libs, sep="\t")
    rows = []
    for (ds, donor, comp), g in df.groupby(["dataset", "donor", "compartment"]):
        r = {"dataset": ds, "donor": donor, "compartment": comp,
             "n_libraries": len(g),
             "n_preQC": int(g.n_preQC.sum()), "n_postQC": int(g.n_postQC.sum()),
             "n_lost_to_QC": int(g.n_lost_to_QC.sum()),
             "median_genes_postQC": float(g.median_genes_postQC.median()),
             "frac_postQC_under_1000_genes": float(
                 np.average(g.frac_postQC_under_1000_genes, weights=g.n_postQC))}
        for gate in GATES:
            for tag in ("preQC", "postQC"):
                k = int(g[f"{gate}_{tag}_n"].sum())
                n = r[f"n_{tag}"]
                r[f"{gate}_{tag}_n"] = k
                r[f"{gate}_{tag}_frac"] = k / n if n else np.nan
                r[f"{gate}_{tag}_X_upper95"] = pf.clopper_pearson_upper(k, n)
        rows.append(r)
    dn = pd.DataFrame(rows)
    dn.to_csv(os.path.join(a.outdir, "early_nk_stages.tsv"), sep="\t", index=False)

    # pooled bound per dataset x compartment, and per-donor worst case
    pool = []
    for (ds, comp), g in dn.groupby(["dataset", "compartment"]):
        p = {"dataset": ds, "compartment": comp, "n_donors": g.donor.nunique(),
             "n_preQC": int(g.n_preQC.sum()), "n_postQC": int(g.n_postQC.sum())}
        for gate in GATES:
            for tag in ("preQC", "postQC"):
                k = int(g[f"{gate}_{tag}_n"].sum())
                n = p[f"n_{tag}"]
                p[f"{gate}_{tag}_n"] = k
                p[f"{gate}_{tag}_frac"] = k / n if n else np.nan
                p[f"{gate}_{tag}_X_upper95"] = pf.clopper_pearson_upper(k, n)
                p[f"{gate}_{tag}_X_upper95_maxdonor"] = float(
                    g[f"{gate}_{tag}_X_upper95"].max())
        pool.append(p)
    po = pd.DataFrame(pool)
    po.to_csv(os.path.join(a.outdir, "early_nk_pooled.tsv"), sep="\t", index=False)
    json.dump(VERDICT_TABLE, open(os.path.join(a.outdir, "task_c_verdict_table.json"), "w"),
              indent=1)
    cols = ["dataset", "compartment", "n_donors", "n_preQC", "n_postQC"]
    for gate in GATES:
        cols += [f"{gate}_preQC_n", f"{gate}_postQC_n"]
    print(po[cols].to_string(index=False))
    print()
    show = ["dataset", "compartment", "n_postQC",
            "earlyNK_permissive_postQC_n", "earlyNK_permissive_postQC_X_upper95",
            "earlyNK_permissive_preQC_n", "earlyNK_permissive_preQC_X_upper95",
            "proB_control_postQC_n", "eryProg_control_postQC_n",
            "matureNK_marker_postQC_n"]
    print(po[show].to_string(index=False))


if __name__ == "__main__":
    main()
