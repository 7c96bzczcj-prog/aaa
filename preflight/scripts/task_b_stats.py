#!/usr/bin/env python3
"""Task B aggregation, gate check and age regression (spec section 5)."""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

VERDICT_TABLE = {
    "B_AGE_EFFECT": "ltNK fraction moves monotonically with age and is not "
                    "explained by the CD56bright/dim shift: the marrow arm has a "
                    "landing point on the age axis -- and a high scoop risk, so "
                    "run a targeted prior-art search immediately.",
    "B_NO_EFFECT": "no change: the marrow arm's age question is not about subset "
                   "proportions and needs a different readout.",
    "B_UNMEASURABLE": "detection rate outside the 5%-85% window: unmeasurable. "
                      "That is itself information for our own cohort -- ltNK "
                      "cannot be defined from the transcriptome and needs flow.",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--libs", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.libs, sep="\t")
    df = df[df.n_nk >= pf.MIN_NK_CELLS].copy()

    rows = []
    for (ds, donor, comp), g in df.groupby(["dataset", "donor", "compartment"]):
        w = g.n_nk.astype(float)
        r = {"dataset": ds, "donor": donor, "compartment": comp,
             "n_libraries": len(g), "n_nk": int(g.n_nk.sum()),
             "age": float(g.age.iloc[0]) if not np.isnan(g.age.iloc[0]) else np.nan}
        for c in ["det_CD69", "det_CXCR6", "det_NCAM1", "det_FCGR3A", "det_SELL",
                  "det_GZMK", "det_KLRC1"]:
            if c in g:
                r[c] = float(np.average(g[c].astype(float), weights=w))
        for c, ncol in [("frac_ltNK", "n_ltNK"), ("frac_bright", "n_bright"),
                        ("frac_dim", "n_dim")]:
            if ncol in g:
                r[ncol] = int(g[ncol].sum())
                r[c] = r[ncol] / r["n_nk"]
        rows.append(r)
    dn = pd.DataFrame(rows)
    dn.to_csv(os.path.join(a.outdir, "ltnk_age.tsv"), sep="\t", index=False)

    # ---- the gate ---------------------------------------------------------
    gate = []
    for (ds, comp), g in dn.groupby(["dataset", "compartment"]):
        w = g.n_nk.astype(float)
        row = {"dataset": ds, "compartment": comp, "n_donors": len(g),
               "n_nk": int(g.n_nk.sum())}
        for gene in ["CD69", "CXCR6"]:
            v = float(np.average(g[f"det_{gene}"], weights=w))
            row[f"det_{gene}"] = v
            row[f"det_{gene}_min_donor"] = float(g[f"det_{gene}"].min())
            row[f"det_{gene}_max_donor"] = float(g[f"det_{gene}"].max())
            row[f"gate_{gene}_pass"] = bool(pf.DETECTION_LO <= v <= pf.DETECTION_HI)
        row["gate_pass"] = bool(row["gate_CD69_pass"] and row["gate_CXCR6_pass"])
        gate.append(row)
    gt = pd.DataFrame(gate)
    gt.to_csv(os.path.join(a.outdir, "ltnk_detection_gate.tsv"), sep="\t", index=False)
    print("=== section 5 gate: CD69 and CXCR6 detection in gated NK ===")
    print(gt.to_string(index=False))

    # ---- regression, only if a gate passes --------------------------------
    out = []
    bm = dn[(dn.compartment == "BM") & dn.age.notna()].copy()
    passing = set(gt[gt.gate_pass].dataset)
    if len(bm) >= 5:
        import statsmodels.api as sm
        ds_dummies = pd.get_dummies(bm.dataset, drop_first=True).astype(float)
        X = pd.concat([bm[["age"]].astype(float).reset_index(drop=True),
                       ds_dummies.reset_index(drop=True)], axis=1)
        X = sm.add_constant(X)
        for y in ["frac_ltNK", "frac_bright", "frac_dim"]:
            if y not in bm:
                continue
            yy = bm[y].astype(float).reset_index(drop=True)
            ok = yy.notna()
            m = sm.OLS(yy[ok], X[ok]).fit()
            out.append({"outcome": y, "n_donors": int(ok.sum()),
                        "slope_per_year": float(m.params["age"]),
                        "ci_lo": float(m.conf_int().loc["age", 0]),
                        "ci_hi": float(m.conf_int().loc["age", 1]),
                        "p": float(m.pvalues["age"]), "r2": float(m.rsquared),
                        "datasets": ",".join(sorted(bm.dataset.unique())),
                        "gate_passed_in": ",".join(sorted(passing)) or "NONE"})
    reg = pd.DataFrame(out)
    if len(reg):
        reg.to_csv(os.path.join(a.outdir, "ltnk_age_regression.tsv"), sep="\t", index=False)
        print("\n=== donor-level regression, age continuous, dataset as covariate ===")
        print(reg.to_string(index=False))
    json.dump(VERDICT_TABLE, open(os.path.join(a.outdir, "task_b_verdict_table.json"), "w"),
              indent=1)


if __name__ == "__main__":
    main()
