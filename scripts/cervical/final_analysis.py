"""
Final deliverable.

Structure follows what the data can actually support:

  PRIMARY   the three same-donor adjacent-normal/tumour pairs in E-MTAB-12305,
            n=3, stated as n=3, with exact tests and binomial intervals.
  CONTEXT   everything else. Cross-donor and cross-dataset absolute values are
            not comparable at the observed between-donor variance and are
            reported only to bound the range.

Every proportion carries a Wilson interval, because several rest on single-digit
NK counts and a bare percentage would hide that.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

OUT = "/home/user/aaa/results/cervical"
DATASETS = ["E-MTAB-12305", "GSE208653", "GSE197461", "GSE173231"]
TISSUE_LABEL = {
    "normal_healthy": "normal cervix (healthy donor)",
    "normal": "normal cervix",
    "normal_adj": "normal cervix (tumour-adjacent)",
    "HSIL": "HSIL / CIN",
    "tumor": "cervical cancer",
    "lymph_node": "metastatic lymph node",
}
# a donor row is flagged when the denominator cannot support a percentage
MIN_NK_TRUST = 20
MIN_LYMPH_TRUST = 200


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * max(0.0, c - h), 100 * min(1.0, c + h))


def load() -> pd.DataFrame:
    frames = []
    for d in DATASETS:
        p = f"{OUT}/nk_per_sample_{d}.csv"
        if os.path.exists(p):
            frames.append(pd.read_csv(p))
        else:
            print(f"  [missing] {d}")
    df = pd.concat(frames, ignore_index=True)
    df["tissue_label"] = df.tissue.map(TISSUE_LABEL).fillna(df.tissue)
    lo, hi = zip(*[wilson(k, n) for k, n in zip(df.n_NK, df.n_lymphocyte)])
    df["lymph_ci_lo"], df["lymph_ci_hi"] = np.round(lo, 2), np.round(hi, 2)
    df["flag"] = ""
    df.loc[df.n_NK < MIN_NK_TRUST, "flag"] += "few_NK;"
    df.loc[df.n_lymphocyte < MIN_LYMPH_TRUST, "flag"] += "few_lymphocytes;"
    return df


def paired(df):
    """The primary result: same donor, same batch, same dissociation."""
    rows = []
    for ds, g in df.groupby("dataset"):
        for donor, gd in g.groupby("donor"):
            tis = set(gd.tissue)
            if "tumor" in tis and ("normal_adj" in tis or "normal" in tis):
                n = gd[gd.tissue.isin(["normal_adj", "normal"])].iloc[0]
                t = gd[gd.tissue == "tumor"].iloc[0]
                # Fisher on the raw cell counts, not on the percentages
                table = [[int(n.n_NK), int(n.n_lymphocyte - n.n_NK)],
                         [int(t.n_NK), int(t.n_lymphocyte - t.n_NK)]]
                try:
                    orr, pv = stats.fisher_exact(table)
                except Exception:
                    orr, pv = np.nan, np.nan
                rows.append({
                    "dataset": ds, "donor": donor,
                    "normal": n["sample"], "tumor": t["sample"],
                    "NK_normal": int(n.n_NK), "lymph_normal": int(n.n_lymphocyte),
                    "NK_tumor": int(t.n_NK), "lymph_tumor": int(t.n_lymphocyte),
                    "pct_normal": round(n.pct_of_lymphocyte, 2),
                    "ci_normal": f"{n.lymph_ci_lo:.2f}-{n.lymph_ci_hi:.2f}",
                    "pct_tumor": round(t.pct_of_lymphocyte, 2),
                    "ci_tumor": f"{t.lymph_ci_lo:.2f}-{t.lymph_ci_hi:.2f}",
                    "fold_change": round(t.pct_of_lymphocyte / n.pct_of_lymphocyte, 2)
                    if n.pct_of_lymphocyte else np.nan,
                    "fisher_p": round(pv, 4) if pv == pv else np.nan,
                    "informative": "yes" if (n.n_NK >= 10 and t.n_NK >= 10) else "NO - too few NK",
                })
    return pd.DataFrame(rows)


def capture_bias():
    """
    Measure the confound that runs in the same direction as the expected answer.

    The positive side of the gate depends on DETECTING NKG7/GNLY/NCAM1/KLRF1.
    Detection is depth-limited, so if tumour NK carry less endogenous RNA they
    are preferentially missed, and 'fewer NK in tumour' is partly manufactured
    by the gate. Compared against the T cells of the same library, which absorbs
    the library's own depth.
    """
    import scanpy as sc
    CACHE = "/home/user/cervical_work/cache"
    rows = []
    for ds in DATASETS:
        p = f"{CACHE}/tnk_{ds}.h5ad"
        if not os.path.exists(p):
            continue
        obs = sc.read_h5ad(p).obs
        if "is_NK" not in obs:
            continue
        for smp, g in obs.groupby("sample", observed=True):
            nk, tc = g[g.is_NK.astype(bool)], g[~g.is_NK.astype(bool)]
            if len(nk) < 10 or len(tc) < 30:
                continue
            rows.append({
                "dataset": ds, "sample": smp, "tissue": g.tissue.iloc[0],
                "n_NK": len(nk),
                "NK_med_UMI": round(float(nk.total_counts.median())),
                "T_med_UMI": round(float(tc.total_counts.median())),
                "NK_over_T_UMI": round(float(nk.total_counts.median()
                                             / tc.total_counts.median()), 3),
            })
    return pd.DataFrame(rows)


def main():
    df = load()

    print("=" * 100)
    print("PRIMARY RESULT - within-donor pairs (the only comparison holding donor,")
    print("batch and dissociation constant)")
    print("=" * 100)
    pr = paired(df)
    if len(pr):
        pr.to_csv(f"{OUT}/NK_PAIRED.csv", index=False)
        print(pr.to_string(index=False))
        ok = pr[pr.informative == "yes"]
        print(f"\n  pairs available: {len(pr)}   adequately powered: {len(ok)}")
    else:
        print("  none")

    print("\n" + "=" * 100)
    print("CONTEXT - per-donor table, all datasets")
    print("=" * 100)
    cols = ["dataset", "donor", "sample", "tissue_label", "histology",
            "n_cells_total", "n_immune", "n_lymphocyte", "n_NK",
            "pct_of_all", "pct_of_immune", "pct_of_lymphocyte",
            "lymph_ci_lo", "lymph_ci_hi", "pct_bright_of_NK", "flag"]
    tab = df[cols].copy().sort_values(["dataset", "tissue_label", "donor"])
    tab.to_csv(f"{OUT}/NK_TABLE.csv", index=False)
    print(tab.to_string(index=False))

    print("\n" + "=" * 100)
    print("CONTEXT - spread by tissue, WITHIN each dataset only")
    print("=" * 100)
    rows = []
    for (ds, tl), g in df.groupby(["dataset", "tissue_label"]):
        for col, name in [("pct_of_all", "of all cells"),
                          ("pct_of_immune", "of CD45+"),
                          ("pct_of_lymphocyte", "of lymphocytes")]:
            v = g[col].dropna()
            if len(v):
                rows.append({"dataset": ds, "tissue": tl, "denominator": name,
                             "n_samples": len(v), "median": round(v.median(), 2),
                             "min": round(v.min(), 2), "max": round(v.max(), 2),
                             "total_NK": int(g.n_NK.sum())})
    sm = pd.DataFrame(rows)
    sm.to_csv(f"{OUT}/NK_SUMMARY.csv", index=False)
    print(sm.to_string(index=False))

    print("\n" + "=" * 100)
    print("CAPTURE BIAS - runs in the same direction as the expected answer")
    print("=" * 100)
    cap = capture_bias()
    if len(cap):
        cap.to_csv(f"{OUT}/NK_CAPTURE.csv", index=False)
        print(cap.to_string(index=False))
        byt = cap.groupby("tissue")["NK_over_T_UMI"].median().round(3)
        print("\n  median NK/T UMI ratio by tissue:")
        print(byt.to_string())


if __name__ == "__main__":
    main()
