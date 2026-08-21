#!/usr/bin/env python3
"""
S-1a §S1a.2 — build a bright/dim stratifier that excludes every outcome gene.

Defined from FCGR3A and NCAM1 only. XCL1, XCL2, CCL3, CCL4 and CCL5 are
explicitly excluded, so the label cannot encode the outcome.

Two threshold schemes, reported side by side, because NCAM1 drops out
badly on 10x:
  strict : NCAM1 > 0  AND  FCGR3A == 0            -> bright
           FCGR3A > 0                              -> dim
           NCAM1 == 0 AND FCGR3A == 0              -> UNCLASSIFIABLE
  loose  : FCGR3A == 0                             -> bright
           FCGR3A > 0                              -> dim

Thresholds are raw-count based and identical across every organ and
chemistry; they are never tuned per organ.

Reports the cell-level confusion matrix against Manually_curated_celltype
and the per-organ bright-share difference (clean minus original), then
applies §S1a.2's four pre-written verdicts.

Reads only the cached per-cell table. No re-download, no re-scan.
No matching, no stratified estimate, no decomposition (§S1a.5).
"""
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
pd.set_option("display.width", 220)

OUT = "results"
OUTCOME_GENES = ["XCL1", "XCL2", "CCL3", "CCL4", "CCL5"]


def main():
    c = pd.read_parquet(f"{OUT}/organ_cells_cache.parquet")
    print("=" * 78)
    print("S-1a §S1a.2  CLEAN STRATIFIER (outcome genes excluded by construction)")
    print("=" * 78)
    print(f"excluded from the label: {OUTCOME_GENES}")
    print("built from: FCGR3A, NCAM1 only")
    print("thresholds (identical across all organs and chemistries, raw counts):")
    print("   strict: bright = NCAM1>0 & FCGR3A==0 ; dim = FCGR3A>0 ; "
          "else UNCLASSIFIABLE")
    print("   loose : bright = FCGR3A==0           ; dim = FCGR3A>0")

    # restrict to author-annotated NK: we are re-drawing the bright/dim
    # stratifier within the same NK set, so the comparison is like for like.
    nk = c.author_ct.astype(str).isin(["NK_CD16+", "NK_CD56bright_CD16-"])
    d = c[nk].copy()
    print(f"\nauthor-annotated NK cells: {len(d)}")

    ncam, fcgr = d["NCAM1"].values, d["FCGR3A"].values
    strict = np.where(fcgr > 0, "dim",
                      np.where(ncam > 0, "bright", "UNCLASSIFIABLE"))
    loose = np.where(fcgr > 0, "dim", "bright")
    d["strict"], d["loose"] = strict, loose
    d["orig"] = np.where(d.author_ct.astype(str) == "NK_CD56bright_CD16-",
                         "bright", "dim")

    print("\noverall label distribution")
    for k in ("orig", "strict", "loose"):
        vc = d[k].value_counts()
        print(f"   {k:7s} " + "  ".join(f"{a}={b} ({100*b/len(d):.1f}%)"
                                        for a, b in vc.items()))
    unc = float((d.strict == "UNCLASSIFIABLE").mean())
    print(f"\nUNCLASSIFIABLE under strict: {unc:.1%} "
          "(NCAM1 and FCGR3A both zero — the NCAM1 dropout problem)")

    # ---- confusion matrices ------------------------------------------
    print("\n" + "-" * 78)
    print("CONFUSION MATRIX vs Manually_curated_celltype (cell level)")
    print("-" * 78)
    for scheme in ("strict", "loose"):
        print(f"\n--- clean/{scheme} vs original ---")
        tab = pd.crosstab(d["orig"], d[scheme])
        print(tab.to_string())
        print("row-normalised (%):")
        print((tab.div(tab.sum(axis=1), axis=0) * 100).round(1).to_string())
        if scheme == "loose":
            agree = float((d.orig == d.loose).mean())
            print(f"   overall agreement: {agree:.3f}")
        else:
            m = d[d.strict != "UNCLASSIFIABLE"]
            agree = float((m.orig == m.strict).mean())
            print(f"   agreement excluding UNCLASSIFIABLE "
                  f"(n={len(m)}): {agree:.3f}")

    # ---- per-organ bright share --------------------------------------
    print("\n" + "-" * 78)
    print("PER-ORGAN BRIGHT SHARE: clean minus original  (pp)")
    print("§4: reported per (tissue x assay); never merged across chemistry")
    print("-" * 78)
    rows = []
    for (tis, asy), g in d.groupby(["tissue", "assay"], observed=True):
        if len(g) < 50:
            continue
        o = float((g.orig == "bright").mean())
        s_den = (g.strict != "UNCLASSIFIABLE").sum()
        s = float((g.strict == "bright").sum() / s_den) if s_den else np.nan
        l = float((g.loose == "bright").mean())
        rows.append(dict(tissue=tis, assay=asy, n=len(g), orig=o,
                         strict=s, loose=l,
                         d_strict=100 * (s - o), d_loose=100 * (l - o),
                         unclass=float((g.strict == "UNCLASSIFIABLE").mean())))
    P = pd.DataFrame(rows).sort_values(["assay", "tissue"])
    P.to_csv(f"{OUT}/s1a_bright_share.csv", index=False)
    print(f"{'tissue':22s} {'assay':10s} {'n':>6s} {'orig':>7s} {'strict':>7s} "
          f"{'loose':>7s} {'d_str':>7s} {'d_loo':>7s} {'unclass':>8s}")
    for _, r in P.iterrows():
        print(f"{r.tissue[:22]:22s} {r.assay:10s} {r.n:6d} {r.orig:7.3f} "
              f"{r.strict:7.3f} {r.loose:7.3f} {r.d_strict:+7.1f} "
              f"{r.d_loose:+7.1f} {r.unclass:8.1%}")

    # ---- ranking check + verdict -------------------------------------
    print("\n" + "-" * 78)
    print("§S1a.2 PRE-WRITTEN VERDICT")
    print("-" * 78)
    for asy, g in P.groupby("assay"):
        if len(g) < 3:
            continue
        r_o = g.orig.rank(ascending=False)
        r_s = g.strict.rank(ascending=False)
        r_l = g.loose.rank(ascending=False)
        rho_s = g[["orig", "strict"]].corr(method="spearman").iloc[0, 1]
        rho_l = g[["orig", "loose"]].corr(method="spearman").iloc[0, 1]
        print(f"\n--- {asy} ({len(g)} organs) ---")
        print(f"   |d_strict| median {g.d_strict.abs().median():.1f} pp, "
              f"max {g.d_strict.abs().max():.1f} pp")
        print(f"   |d_loose|  median {g.d_loose.abs().median():.1f} pp, "
              f"max {g.d_loose.abs().max():.1f} pp")
        print(f"   organ-order preserved?  strict rho={rho_s:+.2f}  "
              f"loose rho={rho_l:+.2f}")
        print(f"   order changed (strict): {not (r_o.values == r_s.values).all()}"
              f" ; (loose): {not (r_o.values == r_l.values).all()}")
        for lbl, col in (("strict", "d_strict"), ("loose", "d_loose")):
            med = g[col].abs().median()
            v = ("影响有限 (<5pp)" if med < 5 else
                 "须整体重估 (>10pp)" if med > 10 else "部分可用 (5-10pp)")
            print(f"   -> {lbl}: {v}")
        agree_dir = np.sign(g.d_strict) .equals(np.sign(g.d_loose))
        print(f"   strict and loose agree in sign per organ: {agree_dir}")


if __name__ == "__main__":
    main()
