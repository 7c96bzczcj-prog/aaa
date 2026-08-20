#!/usr/bin/env python3
"""
Pre-checks demanded before ANY Kitagawa decomposition on E-MTAB-10176.

Two questions, both prior to the decomposition itself:

 (1) LABEL INDEPENDENCE. Are the NK-1/2/3 (and CD8T-1/2/3) labels independent
     of the Green/Red contrast? If the clustering was run on Green and Red
     cells jointly, the labels are partly DEFINED by the contrast, and a
     composition/state decomposition is circular. Diagnostic: the colour
     composition WITHIN each subcluster. A colour-pure cluster is the
     signature of a label that tracks the contrast.

 (2) REWEIGHTING IDENTIFIABILITY. Standardising one side to the other's
     composition up-weights small subgroups. Report, for every subcluster
     carrying >20% of the target weight, how many CELLS actually remain on
     each side, plus Kish effective sample size. Too few cells on either
     side -> UNIDENTIFIABLE, and no clean-looking attribution is produced.
     A subcluster present on one side and absent on the other makes
     standardisation in that direction UNDEFINED, not merely noisy.

Record only. This dataset already fails §4 admission.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)

H5AD = "data/raw/DW_T_NK.h5ad"
OUT = "results"
WEIGHT_FLAG = 0.20     # report any subcluster carrying >20% of target weight
MIN_CELLS = 30         # the project's own A4 threshold, reused as the tripwire


def main():
    a = ad.read_h5ad(H5AD)
    obs = a.obs

    print("=" * 78)
    print("PRE-CHECKS BEFORE ANY KITAGAWA DECOMPOSITION — E-MTAB-10176")
    print("=" * 78)

    # ---------- (1) label independence ---------------------------------
    print("\n" + "-" * 78)
    print("(1) LABEL INDEPENDENCE — colour composition WITHIN each subcluster")
    print("-" * 78)
    print("A cluster that is ~pure one colour is a label defined by the contrast.")

    rows = []
    for lin in ["NK", "CD8T"]:
        sub = obs[obs.main_celltype == lin]
        print(f"\n--- {lin} ---")
        print(f"   {'subcluster':14s} {'n':>6s} {'nGreen':>7s} {'nRed':>6s} "
              f"{'%Green':>7s}  purity")
        tab = pd.crosstab(sub["celltype_cluster"], sub["colour"])
        tab = tab.loc[tab.sum(axis=1) > 0]
        for cl in tab.index:
            ng = int(tab.loc[cl, "Green"]) if "Green" in tab.columns else 0
            nr = int(tab.loc[cl, "Red"]) if "Red" in tab.columns else 0
            n = ng + nr
            pg = ng / n if n else np.nan
            purity = max(pg, 1 - pg)
            flag = ""
            if ng == 0 or nr == 0:
                flag = "  <== COLOUR-PURE (one side has ZERO cells)"
            elif purity >= 0.90:
                flag = "  <== >=90% one colour"
            print(f"   {str(cl):14s} {n:6d} {ng:7d} {nr:6d} {100*pg:7.1f}{flag}")
            rows.append(dict(lineage=lin, subcluster=str(cl), n=n, n_green=ng,
                             n_red=nr, frac_green=pg, purity=purity,
                             colour_pure=(ng == 0 or nr == 0)))
    L = pd.DataFrame(rows)
    L.to_csv(f"{OUT}/a3_label_independence.csv", index=False)

    # how much does colour explain subcluster assignment?
    print("\n   --- summary ---")
    for lin in ["NK", "CD8T"]:
        s = L[L.lineage == lin]
        npure = int(s.colour_pure.sum())
        n90 = int((s.purity >= 0.90).sum())
        print(f"   {lin}: {len(s)} subclusters, {npure} colour-PURE, "
              f"{n90} at >=90% one colour, max purity {s.purity.max():.3f}")

    # evidence on how the clustering was built
    print("\n   --- provenance evidence from the object ---")
    print(f"   obsm keys: {list(a.obsm.keys())}  (X_pca over ALL "
          f"{a.n_obs} cells = joint embedding)")
    print(f"   obsp keys: {list(a.obsp.keys())}  (neighbour graph over all cells)")
    print(f"   leiden levels: {obs['leiden'].nunique()}; "
          f"leiden-2 levels: {obs['leiden-2'].nunique()}")
    print("   -> a single PCA + single neighbour graph + leiden over Green AND Red")
    print("      together means the labels were derived WITH the contrast in them.")

    # ---------- (2) reweighting identifiability -------------------------
    print("\n" + "-" * 78)
    print("(2) REWEIGHTING IDENTIFIABILITY — cells left under standardisation")
    print("-" * 78)

    ident = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            s = obs[(obs.main_celltype == lin) & (obs.hours == hr)]
            tab = pd.crosstab(s["celltype_cluster"], s["colour"])
            tab = tab.loc[tab.sum(axis=1) > 0]
            if "Green" not in tab or "Red" not in tab:
                continue
            ng, nr = tab["Green"], tab["Red"]
            pg, pr = ng / ng.sum(), nr / nr.sum()

            print(f"\n=== {lin} {hr}h ===")
            print(f"   {'subcluster':12s} {'nG':>5s} {'nR':>5s} {'pG':>6s} {'pR':>6s} "
                  f"{'w(R->G)':>8s} {'w(G->R)':>8s}")
            undefined_rg, undefined_gr = [], []
            for cl in tab.index:
                wrg = (pg[cl] / pr[cl]) if pr[cl] > 0 else np.inf
                wgr = (pr[cl] / pg[cl]) if pg[cl] > 0 else np.inf
                if not np.isfinite(wrg):
                    undefined_rg.append(str(cl))
                if not np.isfinite(wgr):
                    undefined_gr.append(str(cl))
                print(f"   {str(cl):12s} {ng[cl]:5d} {nr[cl]:5d} {pg[cl]:6.3f} "
                      f"{pr[cl]:6.3f} {wrg:8.2f} {wgr:8.2f}")

            # cells carrying >20% of the TARGET weight
            print(f"   -- subclusters carrying >{int(100*WEIGHT_FLAG)}% of target weight --")
            for cl in tab.index:
                if pg[cl] > WEIGHT_FLAG:
                    ok = nr[cl] >= MIN_CELLS
                    print(f"      standardise RED to GREEN: '{cl}' gets weight "
                          f"{pg[cl]:.3f} but RED has only {nr[cl]} cells"
                          f"{'' if ok else '   <== BELOW ' + str(MIN_CELLS)}")
                    ident.append(dict(lineage=lin, hours=hr, direction="red->green",
                                      subcluster=str(cl), target_weight=pg[cl],
                                      cells_on_reweighted_side=int(nr[cl]),
                                      sufficient=bool(ok)))
                if pr[cl] > WEIGHT_FLAG:
                    ok = ng[cl] >= MIN_CELLS
                    print(f"      standardise GREEN to RED: '{cl}' gets weight "
                          f"{pr[cl]:.3f} but GREEN has only {ng[cl]} cells"
                          f"{'' if ok else '   <== BELOW ' + str(MIN_CELLS)}")
                    ident.append(dict(lineage=lin, hours=hr, direction="green->red",
                                      subcluster=str(cl), target_weight=pr[cl],
                                      cells_on_reweighted_side=int(ng[cl]),
                                      sufficient=bool(ok)))

            # Kish effective n under each standardisation
            def kish(counts, w):
                ww = np.repeat(w.values, counts.values)
                ww = ww[np.isfinite(ww)]
                return (ww.sum() ** 2) / (ww ** 2).sum() if len(ww) else 0.0

            wrg = (pg / pr).replace([np.inf, -np.inf], np.nan).fillna(0)
            wgr = (pr / pg).replace([np.inf, -np.inf], np.nan).fillna(0)
            neff_r = kish(nr, wrg)
            neff_g = kish(ng, wgr)
            print(f"   -- Kish effective sample size after reweighting --")
            print(f"      RED   n={nr.sum():4d} -> n_eff={neff_r:7.1f}  "
                  f"({100*neff_r/nr.sum():.1f}% retained)")
            print(f"      GREEN n={ng.sum():4d} -> n_eff={neff_g:7.1f}  "
                  f"({100*neff_g/ng.sum():.1f}% retained)")
            if undefined_rg:
                print(f"   ** standardising RED->GREEN is UNDEFINED for {undefined_rg} "
                      f"(zero RED cells) **")
            if undefined_gr:
                print(f"   ** standardising GREEN->RED is UNDEFINED for {undefined_gr} "
                      f"(zero GREEN cells) **")

            verdict = "IDENTIFIABLE"
            if undefined_rg or undefined_gr:
                verdict = "UNIDENTIFIABLE (non-overlapping support)"
            elif min(neff_r, neff_g) < MIN_CELLS:
                verdict = f"UNIDENTIFIABLE (n_eff < {MIN_CELLS})"
            elif any(not r["sufficient"] for r in ident
                     if r["lineage"] == lin and r["hours"] == hr):
                verdict = "UNIDENTIFIABLE (a >20%-weight subcluster is under-populated)"
            print(f"   >> VERDICT: {verdict}")
            ident.append(dict(lineage=lin, hours=hr, direction="__verdict__",
                              subcluster=verdict, target_weight=np.nan,
                              cells_on_reweighted_side=np.nan, sufficient=None))
    pd.DataFrame(ident).to_csv(f"{OUT}/a3_identifiability.csv", index=False)

    print("\n" + "=" * 78)
    print("Wrote a3_label_independence.csv, a3_identifiability.csv")
    print("=" * 78)


if __name__ == "__main__":
    main()
