#!/usr/bin/env python3
"""Task A statistics and the pre-written verdict lookup (spec section 3).

Order of operations, fixed:
  1. drop libraries with fewer than MIN_NK_CELLS NK cells
  2. aggregate library -> donor, weighted by NK cell count (statistical unit is
     the donor, spec 2.3)
  3. admission criterion C2: a target gene is TESTABLE in a dataset only if its
     NK detection rate sits inside 5%-85% in BOTH compartments
  4. instrument acceptance: the soup model must pass in both compartments
  5. contrast BM vs PB -- paired where the donors are shared, unpaired where
     they are not -- with BH-FDR across the testable genes
  6. read the frozen verdict table
"""
from __future__ import annotations

import argparse, json, os, sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

VERDICT_TABLE = {
    "A_PASS": "Delta_ambient near zero on most markers and rho_BM ~ rho_PB: "
              "compartment soup is not a systematic confounder, paired "
              "Delta(BM-PB) is an interpretable quantity, marrow-arm premise P1 "
              "holds.",
    "A_PARTIAL": "Delta_ambient significantly non-zero on SOME genes: those "
                 "genes' BM-PB differences are uninterpretable and go on a "
                 "blacklist; the rest remain usable.",
    "A_FAIL": "Delta_ambient significantly non-zero on MOST genes, or "
              "rho_BM >> rho_PB: a substantial part of the published BM-vs-PB NK "
              "difference may be a product of compartment ambient RNA. That is a "
              "methodological result in its own right, and it means our own "
              "paired design must move to a protein readout.",
    "A_UNMEASURABLE": "fewer than 2 admitted datasets, or the soup model fails "
                      "its acceptance test: no conclusion, and a filtered matrix "
                      "may not be used to make up the number.",
}


def donor_table(df):
    df = df[df["n_nk"] >= pf.MIN_NK_CELLS].copy()
    cols = (["rho_narrow", "rho_wide", "rho_median", "rho_narrow_w100",
             "accept_frac_within", "accept_frac_within_rhomed"]
            + [f"amb_{g}" for g in pf.TARGET_GENES]
            + [f"amb100_{g}" for g in pf.TARGET_GENES]
            + [f"ambM_{g}" for g in pf.TARGET_GENES]
            + [f"det_{g}" for g in pf.TARGET_GENES]
            + [f"acc_{g}" for g in ["MS4A1", "CD79A", "LYZ", "HBB"]])
    out = []
    for (ds, donor, comp, gate), g in df.groupby(["dataset", "donor", "compartment", "gate"]):
        row = {"dataset": ds, "donor": donor, "compartment": comp, "gate": gate,
               "n_libraries": len(g), "n_nk": int(g["n_nk"].sum()),
               "n_cells_called": int(g["n_cells_called"].sum()),
               "group": g["group"].iloc[0],
               "accept_pass_frac": float(g["accept_pass"].mean())}
        w = g["n_nk"].astype(float).values
        for c in cols:
            if c not in g:
                continue
            v = g[c].astype(float).values
            ok = ~np.isnan(v)
            row[c] = float(np.average(v[ok], weights=w[ok])) if ok.any() else np.nan
        out.append(row)
    return pd.DataFrame(out)


def testable_genes(dn, bm="BM", pb="PB"):
    """Admission criterion C2, per dataset and per compartment."""
    rows = []
    for gene in pf.TARGET_GENES:
        c = f"det_{gene}"
        r = {"gene": gene}
        ok = True
        for comp in (bm, pb):
            g = dn[dn.compartment == comp]
            if not len(g) or c not in g:
                r[f"det_{comp}"] = np.nan
                ok = False
                continue
            v = float(np.average(g[c].astype(float), weights=g.n_nk.astype(float)))
            r[f"det_{comp}"] = v
            if not (pf.DETECTION_LO <= v <= pf.DETECTION_HI):
                ok = False
        r["testable"] = ok
        rows.append(r)
    return pd.DataFrame(rows)


def contrast(dn, paired, bm="BM", pb="PB", amb_prefix="amb_"):
    a_df = dn[dn.compartment == bm]
    b_df = dn[dn.compartment == pb]
    if paired:
        shared = sorted(set(a_df.donor) & set(b_df.donor))
        a_df = a_df.set_index("donor").loc[shared]
        b_df = b_df.set_index("donor").loc[shared]
    rows = []
    quantities = ["rho_narrow", "rho_wide", "rho_median"] + \
                 [f"{amb_prefix}{g}" for g in pf.TARGET_GENES]
    for q in quantities:
        if q not in dn:
            continue
        a = a_df[q].astype(float).values
        b = b_df[q].astype(float).values
        if paired:
            ok = ~np.isnan(a) & ~np.isnan(b)
            a, b = a[ok], b[ok]
            n = len(a)
            t, p = stats.ttest_rel(a, b) if n >= 3 else (np.nan, np.nan)
        else:
            a = a[~np.isnan(a)]
            b = b[~np.isnan(b)]
            n = min(len(a), len(b))
            t, p = stats.ttest_ind(a, b, equal_var=False) if n >= 3 else (np.nan, np.nan)
        rows.append({"quantity": q,
                     "gene": q[len(amb_prefix):] if q.startswith(amb_prefix) else "",
                     "n_bm": len(a), "n_pb": len(b),
                     "mean_bm": float(np.mean(a)) if len(a) else np.nan,
                     "mean_pb": float(np.mean(b)) if len(b) else np.nan,
                     "delta": float(np.mean(a) - np.mean(b)) if len(a) and len(b) else np.nan,
                     "sd_bm": float(np.std(a, ddof=1)) if len(a) > 1 else np.nan,
                     "sd_pb": float(np.std(b, ddof=1)) if len(b) > 1 else np.nan,
                     "paired": paired,
                     "t": float(t) if t == t else np.nan,
                     "p": float(p) if p == p else np.nan})
    out = pd.DataFrame(rows)
    m = out.gene != ""
    out.loc[m, "q_bh"] = pf.bh_fdr(out.loc[m, "p"].values)
    return out


def verdict_for(con, tst, dn):
    g = con[(con.gene != "")].merge(tst[["gene", "testable"]], on="gene", how="left")
    g = g[g.testable.fillna(False)]
    n_test = len(g)
    sig = g[(g.q_bh < pf.FDR_ALPHA) & (g.delta.abs() >= pf.DELTA_AMBIENT_EFFECT)]
    rho = con[con.quantity == "rho_narrow"].iloc[0] if (con.quantity == "rho_narrow").any() else None
    rho_ratio = (rho.mean_bm / rho.mean_pb) if (rho is not None and rho.mean_pb) else np.nan
    rho_sig = bool(rho is not None and rho.p == rho.p and rho.p < pf.FDR_ALPHA)
    acc = dn.groupby("compartment").accept_pass_frac.mean().to_dict()
    acc_ok = all(v >= 0.5 for v in acc.values()) and len(acc) >= 2
    if not acc_ok:
        key = "A_UNMEASURABLE"
    elif (rho_ratio >= pf.RHO_RATIO_LARGE and rho_sig) or \
         (n_test and len(sig) / n_test > pf.MAJORITY_FRAC):
        key = "A_FAIL"
    elif len(sig) > 0:
        key = "A_PARTIAL"
    else:
        key = "A_PASS"
    return {"n_testable_genes": n_test, "n_significant": len(sig),
            "significant_genes": ",".join(sorted(sig.gene)),
            "rho_bm": float(rho.mean_bm) if rho is not None else np.nan,
            "rho_pb": float(rho.mean_pb) if rho is not None else np.nan,
            "rho_ratio_bm_over_pb": float(rho_ratio),
            "rho_p": float(rho.p) if rho is not None else np.nan,
            "acceptance_pass_frac": json.dumps({k: round(v, 3) for k, v in acc.items()}),
            "verdict_key": key, "verdict": VERDICT_TABLE[key]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--libs", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    df = pd.concat([pd.read_csv(f, sep="\t") for f in a.libs], ignore_index=True)
    df.to_csv(os.path.join(a.outdir, "soup_by_library.tsv"), sep="\t", index=False)
    dn = donor_table(df)
    dn.to_csv(os.path.join(a.outdir, "soup_by_donor.tsv"), sep="\t", index=False)

    all_con, all_tst, verdicts = [], [], []
    for (ds, gate), g in dn.groupby(["dataset", "gate"]):
        tst = testable_genes(g)
        tst.insert(0, "dataset", ds)
        tst.insert(1, "gate", gate)
        all_tst.append(tst)
        paired = len(set(g[g.compartment == "BM"].donor) &
                     set(g[g.compartment == "PB"].donor)) >= 3
        con = contrast(g, paired)
        con.insert(0, "dataset", ds)
        con.insert(1, "gate", gate)
        all_con.append(con)
        v = verdict_for(con, tst, g)
        v.update({"dataset": ds, "gate": gate, "paired": paired,
                  "n_donors_bm": int((g.compartment == "BM").sum()),
                  "n_donors_pb": int((g.compartment == "PB").sum())})
        verdicts.append(v)
    con = pd.concat(all_con, ignore_index=True)
    tst = pd.concat(all_tst, ignore_index=True)
    vd = pd.DataFrame(verdicts)
    con.to_csv(os.path.join(a.outdir, "soup_by_compartment.tsv"), sep="\t", index=False)
    tst.to_csv(os.path.join(a.outdir, "target_gene_detection.tsv"), sep="\t", index=False)
    vd.to_csv(os.path.join(a.outdir, "task_a_verdict.tsv"), sep="\t", index=False)
    json.dump(VERDICT_TABLE, open(os.path.join(a.outdir, "task_a_verdict_table.json"), "w"),
              indent=1)
    print(tst.to_string(index=False))
    print()
    print(con[con.gene != ""].to_string(index=False))
    print()
    print(con[con.gene == ""].to_string(index=False))
    print()
    print(vd.to_string(index=False))


if __name__ == "__main__":
    main()
