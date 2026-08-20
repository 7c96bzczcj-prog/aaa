#!/usr/bin/env python3
"""
Amendment A-1, §A1.4 — polarity consistency check.

Question (pre-registered, three branches fixed before looking):
    In E-MTAB-10176 (Li 2022, Kaede/MC38 total TIL, T and NK sorted together),
    what is the DIRECTION of the CD8 T module ratio under Green vs Red?

    Green = Kaede-Green  = "Infiltrating"  = newly entered
    Red   = Kaede-Red    = "Resident"      = retained in tumour

    branch 1  CD8T ratio also FALLS (Red < Green)  -> consistent with [M] 47%
    branch 2  CD8T ratio RISES     (Red > Green)   -> contradicts [M] 47%
    branch 3  within the A1.3c replicate-noise floor -> UNMEASURED

Modules (v1.0 §6):
    regulatory : Ccl3 Ccl4 Ccl5 Xcl1 Xcl2 (Flt3l)
    killing    : Prf1 Gzma Gzmb Gzmk Nkg7

Hard rules honoured here:
    §5   every gene gets a detection-rate band; floor (<5%) and ceiling (>85%)
         genes NEVER enter a module mean.
    §6.1 per-gene band reported for every module member.
    §6.2 leave-one-out sensitivity on every ratio.
    A1.3a NK and CD8T reported SIDE BY SIDE. No subtraction anywhere.
    A1.3c the replicate-noise floor must be measured in THIS dataset before
         any difference is read. If it cannot be measured, differences are
         not readable.
    §3.1 nothing is left blank; unmeasurable things are named.

Units: detection rate (fraction of cells with a non-zero count), which is the
unit v1.0 §7 mandates. raw.X is log-normalised but zero-preserving, so
detection rate is exact from it. adata.X is scaled/centred (has no zeros) and
is therefore UNUSABLE for detection and is never touched.
"""
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

warnings.filterwarnings("ignore")
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 60)

H5AD = "data/raw/DW_T_NK.h5ad"
OUT = "results"

REGULATORY = ["Ccl3", "Ccl4", "Ccl5", "Xcl1", "Xcl2", "Flt3l"]
KILLING = ["Prf1", "Gzma", "Gzmb", "Gzmk", "Nkg7"]
FLOOR, CEIL = 0.05, 0.85


def band(p):
    if np.isnan(p):
        return "ABSENT"
    if p < FLOOR:
        return "FLOOR"
    if p > CEIL:
        return "CEILING"
    return "MEASURABLE"


def main():
    a = ad.read_h5ad(H5AD)

    # ---- matrix choice, stated explicitly -------------------------------
    assert a.raw is not None, "raw.X absent; detection rate not computable"
    R = a.raw.X
    R = R.tocsr() if sp.issparse(R) else sp.csr_matrix(np.asarray(R))
    genes = pd.Index(a.raw.var_names)
    obs = a.obs.copy()

    print("=" * 78)
    print("A1.4 POLARITY CONSISTENCY CHECK — E-MTAB-10176 (Li 2022)")
    print("=" * 78)
    print(f"cells {R.shape[0]}   genes(raw) {R.shape[1]}   "
          f"genes(X, HVG+scaled, unusable) {a.n_vars}")

    # ---- §4(d) donor structure: the decisive admission check ------------
    print("\n" + "-" * 78)
    print("§4(d)  DONOR vs LIBRARY STRUCTURE")
    print("-" * 78)
    ct = pd.crosstab(obs["sample"], [obs["colour"], obs["hours"]])
    print(ct)
    n_lib = obs["sample"].nunique()
    conf = obs.groupby("sample", observed=True)["colour"].nunique().max()
    print(f"\nlibraries: {n_lib}")
    print(f"max distinct colours within any one library: {conf}")
    print("colour is CONFOUNDED with library" if conf == 1
          else "colour varies within a library (separable)")
    donor_cols = [c for c in obs.columns
                  if any(k in c.lower() for k in ("mouse", "donor", "animal", "replicate"))]
    print(f"donor/mouse/replicate columns present in obs: "
          f"{donor_cols if donor_cols else 'NONE'}")

    # ---- cell counts per lineage x colour x hour ------------------------
    print("\n" + "-" * 78)
    print("CELL COUNTS  lineage x colour x hour")
    print("-" * 78)
    counts = pd.crosstab(obs["main_celltype"], [obs["hours"], obs["colour"]])
    print(counts)

    # ---- per-gene detection rates ---------------------------------------
    allg = REGULATORY + KILLING
    present = [g for g in allg if g in set(genes)]
    absent = [g for g in allg if g not in set(genes)]
    gi = {g: genes.get_loc(g) for g in present}

    def det(mask, g):
        if g not in gi:
            return np.nan
        col = R[:, gi[g]]
        v = np.asarray(col[mask].todense()).ravel()
        return float((v > 0).mean()) if v.size else np.nan

    rows = []
    for lin in ["NK", "CD8T", "CD4T", "Treg"]:
        for hr in ["24", "72"]:
            for col in ["Green", "Red"]:
                m = ((obs["main_celltype"] == lin) &
                     (obs["hours"] == hr) &
                     (obs["colour"] == col)).values
                if m.sum() == 0:
                    continue
                for g in allg:
                    p = det(m, g)
                    rows.append(dict(lineage=lin, hours=hr, colour=col,
                                     n_cells=int(m.sum()), gene=g,
                                     module=("regulatory" if g in REGULATORY else "killing"),
                                     detection=p, band=band(p)))
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/a1_4_detection_rates.csv", index=False)

    print("\n" + "-" * 78)
    print("§5  DETECTION-RATE BANDS  (per gene, per lineage, per compartment)")
    print("-" * 78)
    print(f"genes absent from the deposited matrix entirely: "
          f"{absent if absent else 'none'}")
    for lin in ["NK", "CD8T"]:
        print(f"\n--- {lin} ---")
        piv = D[D.lineage == lin].pivot_table(
            index=["module", "gene"], columns=["hours", "colour"],
            values="detection", observed=True)
        bpv = D[D.lineage == lin].pivot_table(
            index=["module", "gene"], columns=["hours", "colour"],
            values="band", aggfunc="first", observed=True)
        print(piv.round(4))
        print(bpv)

    # ---- §5 admission: a gene enters only if MEASURABLE in BOTH sides ----
    print("\n" + "-" * 78)
    print("§5 ADMISSION  (gene must be MEASURABLE in BOTH Green and Red)")
    print("-" * 78)
    admit = {}
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            keep = []
            for g in allg:
                sub = D[(D.lineage == lin) & (D.hours == hr) & (D.gene == g)]
                if len(sub) == 2 and (sub["band"] == "MEASURABLE").all():
                    keep.append(g)
            admit[(lin, hr)] = keep
            excl = [g for g in allg if g not in keep]
            print(f"{lin} {hr}h  admitted={keep}")
            print(f"{'':>9} excluded={excl}")

    # ---- module ratio ----------------------------------------------------
    def ratio(lin, hr, col, drop=None, admitted=None):
        keep = admitted if admitted is not None else admit[(lin, hr)]
        if drop:
            keep = [g for g in keep if g != drop]
            reg = [g for g in keep if g in REGULATORY]
            kil = [g for g in keep if g in KILLING]
        reg = [g for g in keep if g in REGULATORY]
        kil = [g for g in keep if g in KILLING]
        if not reg or not kil:
            return np.nan, reg, kil
        s = D[(D.lineage == lin) & (D.hours == hr) & (D.colour == col)]
        rv = s[s.gene.isin(reg)]["detection"].mean()
        kv = s[s.gene.isin(kil)]["detection"].mean()
        return (rv / kv if kv > 0 else np.nan), reg, kil

    print("\n" + "-" * 78)
    print("A1.4  MODULE RATIO  regulatory / killing   (detection-rate units)")
    print("A1.3a  NK and CD8T reported SIDE BY SIDE — never subtracted")
    print("-" * 78)
    res = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            rg, reg, kil = ratio(lin, hr, "Green")
            rr, _, _ = ratio(lin, hr, "Red")
            if np.isnan(rg) or np.isnan(rr):
                direction = "UNCOMPUTABLE (a module lost all its genes to §5)"
                delta = np.nan
            else:
                delta = rr - rg
                direction = "FALLS (Red<Green)" if delta < 0 else "RISES (Red>Green)"
            ng = int(((obs.main_celltype == lin) & (obs.hours == hr) &
                      (obs.colour == "Green")).sum())
            nr = int(((obs.main_celltype == lin) & (obs.hours == hr) &
                      (obs.colour == "Red")).sum())
            res.append(dict(lineage=lin, hours=hr, n_green=ng, n_red=nr,
                            reg_genes=",".join(reg), kill_genes=",".join(kil),
                            ratio_green=rg, ratio_red=rr, delta=delta,
                            direction=direction))
            print(f"\n{lin} {hr}h   n_green={ng:4d}  n_red={nr:4d}")
            print(f"   reg genes used : {reg}")
            print(f"   kill genes used: {kil}")
            print(f"   ratio Green={rg:.4f}   ratio Red={rr:.4f}   "
                  f"delta={delta:+.4f}   {direction}")
    Rres = pd.DataFrame(res)
    Rres.to_csv(f"{OUT}/a1_4_module_ratio.csv", index=False)

    # ---- §6.2 leave-one-out ---------------------------------------------
    print("\n" + "-" * 78)
    print("§6.2  LEAVE-ONE-OUT SENSITIVITY  (does dropping any gene flip it?)")
    print("-" * 78)
    lrows = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            base = [r for r in res if r["lineage"] == lin and r["hours"] == hr][0]
            if np.isnan(base["delta"]):
                continue
            base_sign = np.sign(base["delta"])
            print(f"\n{lin} {hr}h  baseline delta={base['delta']:+.4f}")
            for g in admit[(lin, hr)]:
                keep = [x for x in admit[(lin, hr)] if x != g]
                rg, reg, kil = ratio(lin, hr, "Green", admitted=keep)
                rr, _, _ = ratio(lin, hr, "Red", admitted=keep)
                if np.isnan(rg) or np.isnan(rr):
                    print(f"   drop {g:6s} -> UNCOMPUTABLE (module emptied)")
                    lrows.append(dict(lineage=lin, hours=hr, dropped=g,
                                      delta=np.nan, flips=None))
                    continue
                d = rr - rg
                flip = np.sign(d) != base_sign
                print(f"   drop {g:6s} -> delta={d:+.4f}  {'** FLIPS **' if flip else 'stable'}")
                lrows.append(dict(lineage=lin, hours=hr, dropped=g,
                                  delta=d, flips=bool(flip)))
    pd.DataFrame(lrows).to_csv(f"{OUT}/a1_4_leave_one_out.csv", index=False)

    # ---- A1.3c replicate-noise floor ------------------------------------
    print("\n" + "-" * 78)
    print("A1.3c  REPLICATE-NOISE FLOOR (must be measured in THIS dataset)")
    print("-" * 78)
    same_cond = (obs.groupby(["colour", "hours"], observed=True)["sample"]
                 .nunique())
    print("libraries per (colour,hours) cell:")
    print(same_cond)
    if same_cond.max() < 2:
        print("\n-> NO same-condition replicate libraries exist in this dataset.")
        print("-> The replicate-noise floor CANNOT be measured here.")
        print("-> Per A1.3c, no difference computed above is readable against a floor.")

    # ---- composition diagnostic (bears on §7, not a decomposition) -------
    print("\n" + "-" * 78)
    print("COMPOSITION DIAGNOSTIC  (sub-cluster mix shifts Green->Red)")
    print("-" * 78)
    for lin in ["NK", "CD8T"]:
        sub = obs[obs.main_celltype == lin]
        tab = pd.crosstab(sub["celltype_cluster"], [sub["hours"], sub["colour"]],
                          normalize="columns")
        tab = tab.loc[(tab.sum(axis=1) > 0)]
        print(f"\n--- {lin} sub-cluster fractions ---")
        print(tab.round(3))

    print("\n" + "=" * 78)
    print("Wrote results/a1_4_detection_rates.csv, a1_4_module_ratio.csv, "
          "a1_4_leave_one_out.csv")
    print("=" * 78)


if __name__ == "__main__":
    main()
