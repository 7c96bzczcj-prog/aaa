#!/usr/bin/env python3
"""
S-1 §3 (per-organ tables, both labellings) and §5 (donor bootstrap).

§3: one row per organ x labelling x chemistry, with every qualifying
condition written out, so a reader can see whether an organ's XCL1 is low
because its bright share is low or because its bright cells are low. Three
labellings appear: the original Manually_curated_celltype, and S1a.2's
clean strict and clean loose schemes.

§5: donor-level bootstrap, 2000 resamples, 95% intervals for detection
rates, cross-organ spread, bright-restricted spread, and the rank
correlation with bright share. Pre-written rule: if a spread interval's
lower bound falls below 1.2x, the claim that the organs differ at all does
not stand for that gene.

Reads only the cache. No matching, no stratified estimate, no
decomposition.
"""
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 240)

OUT = "results"
GENES = ["CCL3", "CCL4", "CCL5", "XCL1", "XCL2"]
MIN_DONORS, MIN_NK, FLOOR = 3, 200, 0.05
NBOOT = 2000
RNG = np.random.default_rng(20260821)


def label_sets(d):
    ncam, fcgr = d["NCAM1"].values, d["FCGR3A"].values
    orig = np.where(d.author_ct.astype(str) == "NK_CD56bright_CD16-",
                    "bright", "dim")
    strict = np.where(fcgr > 0, "dim", np.where(ncam > 0, "bright", "UNC"))
    loose = np.where(fcgr > 0, "dim", "bright")
    return {"original": orig, "clean_strict": strict, "clean_loose": loose}


def main():
    c = pd.read_parquet(f"{OUT}/organ_cells_cache.parquet")
    nk = c.author_ct.astype(str).isin(["NK_CD16+", "NK_CD56bright_CD16-"])
    d = c[nk].reset_index(drop=True)
    labs = label_sets(d)

    # ---------- §3 per-organ tables --------------------------------
    rows = []
    for lname, lab in labs.items():
        for (tis, asy, don), gi in d.groupby(["tissue", "assay", "donor"],
                                             observed=True).indices.items():
            sub = d.iloc[gi]
            l = lab[gi]
            isbr = (l == "bright")
            rec = dict(labelling=lname, tissue=tis, assay=asy, donor=don,
                       n_nk=len(gi), n_bright=int(isbr.sum()),
                       n_unclass=int((l == "UNC").sum()),
                       median_umi=float(np.median(sub.total_umi)),
                       bright_frac=float(isbr.mean()),
                       stress=float(sub.stress.mean()))
            for g in GENES:
                v = sub[g].values
                rec[g] = float((v > 0).mean())
                rec[g + "_bright"] = float((v[isbr] > 0).mean()) if isbr.any() else np.nan
            rows.append(rec)
    P = pd.DataFrame(rows)
    P.to_csv(f"{OUT}/s1a_per_donor_bylabel.csv", index=False)

    agg = []
    for (lname, tis, asy), g in P.groupby(["labelling", "tissue", "assay"],
                                          observed=True):
        rec = dict(labelling=lname, tissue=tis, assay=asy,
                   donors=g.donor.nunique(), n_nk=int(g.n_nk.sum()),
                   n_bright=int(g.n_bright.sum()),
                   n_unclass=int(g.n_unclass.sum()),
                   median_umi=float(g.median_umi.median()),
                   bright_frac=float(g.bright_frac.median()),
                   stress=float(g.stress.median()))
        for gn in GENES:
            rec[gn] = float(g[gn].median())
            rec[gn + "_bright"] = float(g[gn + "_bright"].median())
        rec["reportable"] = (rec["donors"] >= MIN_DONORS) and (rec["n_nk"] >= MIN_NK)
        agg.append(rec)
    A = pd.DataFrame(agg)
    A.to_csv(f"{OUT}/s1a_organ_table_v2.csv", index=False)
    ok = A[A.reportable]
    print(f"§3: {len(A)} strata, reportable {len(ok)} "
          f"({ok.labelling.value_counts().to_dict()})")

    # ---------- §5 bootstrap ---------------------------------------
    print("\n" + "=" * 78)
    print(f"§5 DONOR-LEVEL BOOTSTRAP  ({NBOOT} resamples)")
    print("rule: if a spread CI lower bound < 1.2x, 'the organs differ' does "
          "not stand")
    print("=" * 78)
    boot_rows = []
    for lname in labs:
        for asy in sorted(P.assay.unique()):
            s = ok[(ok.labelling == lname) & (ok.assay == asy)]
            if len(s) < 3:
                continue
            organs = sorted(s.tissue)
            pool = P[(P.labelling == lname) & (P.assay == asy) &
                     (P.tissue.isin(organs))]
            print(f"\n--- {lname} / {asy}  ({len(organs)} organs) ---")
            for gn in GENES:
                spreads, spreads_b, rhos, per_org = [], [], [], {o: [] for o in organs}
                for _ in range(NBOOT):
                    meds, meds_b, brs = [], [], []
                    for o in organs:
                        g = pool[pool.tissue == o]
                        take = g.iloc[RNG.integers(0, len(g), len(g))]
                        m = np.nanmedian(take[gn].values)
                        meds.append(m); per_org[o].append(m)
                        meds_b.append(np.nanmedian(take[gn + "_bright"].values))
                        brs.append(np.nanmedian(take.bright_frac.values))
                    meds = np.array(meds); meds_b = np.array(meds_b)
                    if np.nanmin(meds) > 0:
                        spreads.append(np.nanmax(meds) / np.nanmin(meds))
                    if np.nanmin(meds_b) > 0:
                        spreads_b.append(np.nanmax(meds_b) / np.nanmin(meds_b))
                    if np.isfinite(meds).all() and np.isfinite(brs).all():
                        rhos.append(stats.spearmanr(meds, brs).statistic)
                def ci(a):
                    a = np.asarray([x for x in a if np.isfinite(x)])
                    return (np.percentile(a, 2.5), np.percentile(a, 97.5)) if len(a) else (np.nan, np.nan)
                sl, sh = ci(spreads); bl, bh = ci(spreads_b); rl, rh = ci(rhos)
                fails = (sl < 1.2)
                print(f"   {gn:5s} spread {np.median(spreads):5.2f}x "
                      f"[{sl:.2f}, {sh:.2f}]"
                      + ("  <== LOWER BOUND < 1.2x: organs do NOT differ" if fails else "")
                      + f"   bright-only [{bl:.2f}, {bh:.2f}]"
                      + f"   rho(bright) [{rl:+.2f}, {rh:+.2f}]")
                boot_rows.append(dict(labelling=lname, assay=asy, gene=gn,
                                      n_organs=len(organs),
                                      spread_med=float(np.median(spreads)),
                                      spread_lo=sl, spread_hi=sh,
                                      spread_bright_lo=bl, spread_bright_hi=bh,
                                      rho_lo=rl, rho_hi=rh,
                                      organs_differ=not fails))
                for o in organs:
                    a = np.array(per_org[o])
                    boot_rows.append(dict(labelling=lname, assay=asy, gene=gn,
                                          tissue=o,
                                          det_med=float(np.nanmedian(a)),
                                          det_lo=float(np.nanpercentile(a, 2.5)),
                                          det_hi=float(np.nanpercentile(a, 97.5))))
    pd.DataFrame(boot_rows).to_csv(f"{OUT}/s1a_bootstrap.csv", index=False)
    print("\nwrote s1a_organ_table_v2.csv, s1a_bootstrap.csv")


if __name__ == "__main__":
    main()
