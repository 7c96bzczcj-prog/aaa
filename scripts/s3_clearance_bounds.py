#!/usr/bin/env python3
"""
S-3 §2 — clearance-fraction bounds from the Kaede green/red contrast.

LIMIT THAT DOES NOT GO AWAY (§7): the green/red contrast cannot separate
clearance from reprogramming. Everything below answers "if it were
clearance, how large would it have to be", never "it is clearance".

Two structural facts about the dataset actually on disk force departures
from §2.4, both reported rather than papered over:

  1. Colour is perfectly confounded with library. The four libraries are
     MH1_G24, MH2_R24, MH3_G72, MH4_R72 -- one library per colour x hour.
     There is no same-colour replicate anywhere, so no technical noise
     floor can be estimated here and the §2.4 signal-to-noise gate cannot
     be built on this dataset.
  2. "Bootstrap over three libraries" is therefore impossible: resampling
     libraries would resample the colour itself. The intervals below come
     from resampling CELLS within colour x hour. They are SAMPLING
     intervals only -- a LOWER BOUND on the true uncertainty, which also
     contains library and animal variance that this design cannot see.

In place of the SNR gate, the project's own measurability bands are used
(<5% floor, 5-85% measurable, >85% ceiling), plus a per-gene cell-support
minimum. The gate is stated, not tuned.

The two timepoints (24 h, 72 h) are carried through separately. They are
the only near-replication this dataset offers: each is an independent
green-vs-red pair, and disagreement between them is itself a result.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

H5 = "data/raw/DW_T_NK.h5ad"
OUT = "results"
NBOOT = 2000
RNG = np.random.default_rng(20260821)
FLOOR, CEIL = 0.05, 0.85
MIN_POS = 10          # a gene needs >=10 positive cells on the larger side
PPT = {               # §2.3, the numbers to be re-checked (green %, red %)
    "Ccl4": (37.85, 12.28), "Ccl3": (12.59, 2.80), "Ccl5": (90.62, 36.52),
    "Itga1": (8.76, 19.34), "Cd160": (7.78, 17.24), "Ldha": (75.07, 88.93),
    "Xcl1": (76.02, 81.95),
}


def f_min(g, r):
    """clearance fraction needed if a fall were entirely clearance."""
    return (g - r) / (1.0 - r) if r < 1 else np.nan


def f_est(g, r):
    """clearance fraction implied if a rise were entirely passive."""
    return 1.0 - g / r if r > 0 else np.nan


def rates(mat, mask):
    """detection rate per gene over the cells selected by mask."""
    if mask.sum() == 0:
        return np.full(mat.shape[1], np.nan)
    return np.asarray((mat[mask] > 0).mean(axis=0)).ravel()


def main():
    a = ad.read_h5ad(H5)
    R = a.raw
    names = np.array(R.var_names)
    X = R.X.tocsr()
    obs = a.obs
    print("=" * 78)
    print("S-3 §2  CLEARANCE-FRACTION BOUNDS")
    print("=" * 78)
    print(f"dataset on disk: E-MTAB-10176 (Kaede, 24/72 h)  "
          f"{a.n_obs} cells x {len(names)} genes (raw)")
    print("libraries:", dict(obs["sample"].value_counts()))
    print("NOTE: colour is perfectly confounded with library "
          "(one library per colour x hour)")

    nk = (obs.main_celltype.astype(str) == "NK").values
    cd8 = (obs.main_celltype.astype(str) == "CD8T").values
    colour = obs.colour.astype(str).values
    hours = obs.hours.astype(str).values
    print(f"\nNK cells: {nk.sum()}   CD8T cells: {cd8.sum()}")
    print(pd.crosstab(colour[nk], hours[nk]).to_string())

    # ---------- 1. §2.3 re-check -----------------------------------
    print("\n" + "-" * 78)
    print("§2.3 RE-CHECK of the PPT-derived rough figures")
    print("-" * 78)
    gi = {g: int(np.where(names == g)[0][0]) for g in PPT if (names == g).any()}
    cols = list(gi)
    sub = X[:, [gi[g] for g in cols]]
    det = np.asarray(sub.todense()) > 0
    D = pd.DataFrame(det, columns=cols)
    D["colour"], D["hours"], D["nk"] = colour, hours, nk
    pooled = D[D.nk].groupby("colour")[cols].mean() * 100
    byhr = D[D.nk].groupby(["colour", "hours"])[cols].mean() * 100
    rows = []
    for g in cols:
        pg, pr = PPT[g]
        mg, mr = pooled.loc["Green", g], pooled.loc["Red", g]
        pf = f_min(pg / 100, pr / 100) if pg > pr else f_est(pg / 100, pr / 100)
        mf = f_min(mg / 100, mr / 100) if mg > mr else f_est(mg / 100, mr / 100)
        rows.append(dict(gene=g, ppt_green=pg, ppt_red=pr, ppt_f=pf,
                         obs_green=round(mg, 2), obs_red=round(mr, 2),
                         obs_f=mf,
                         g24=round(byhr.loc[("Green", "24"), g], 2),
                         r24=round(byhr.loc[("Red", "24"), g], 2),
                         g72=round(byhr.loc[("Green", "72"), g], 2),
                         r72=round(byhr.loc[("Red", "72"), g], 2),
                         dir_ppt="fall" if pg > pr else "rise",
                         dir_obs="fall" if mg > mr else "rise"))
    C = pd.DataFrame(rows)
    C["f_shift"] = C.obs_f - C.ppt_f
    C["dir_agrees"] = C.dir_ppt == C.dir_obs
    C.to_csv(f"{OUT}/s3_ppt_recheck.csv", index=False)
    print(C.round(3).to_string(index=False))

    # ---------- 2. gene-wide sweep ---------------------------------
    print("\n" + "-" * 78)
    print("§2.4 GENE-WIDE SWEEP")
    print("-" * 78)
    res = {}
    for lbl, sel in (("pooled", np.ones(len(hours), bool)),
                     ("24h", hours == "24"), ("72h", hours == "72")):
        gm = nk & sel & (colour == "Green")
        rm = nk & sel & (colour == "Red")
        g, r = rates(X, gm), rates(X, rm)
        gp = np.asarray((X[gm] > 0).sum(axis=0)).ravel()
        rp = np.asarray((X[rm] > 0).sum(axis=0)).ravel()
        tg = rates(X, cd8 & sel & (colour == "Green"))
        tr = rates(X, cd8 & sel & (colour == "Red"))
        res[lbl] = dict(g=g, r=r, gp=gp, rp=rp, tg=tg, tr=tr,
                        ng=int(gm.sum()), nr=int(rm.sum()))
        print(f"  {lbl}: green NK n={gm.sum()}  red NK n={rm.sum()}  "
              f"green CD8T n={(cd8&sel&(colour=='Green')).sum()}  "
              f"red CD8T n={(cd8&sel&(colour=='Red')).sum()}")

    P = res["pooled"]
    g, r = P["g"], P["r"]
    both_up = (g >= FLOOR) & (r >= FLOOR)
    support = np.maximum(P["gp"], P["rp"]) >= MIN_POS
    gate = both_up & support
    print(f"\ngate: both colours >= {FLOOR:.0%} detection AND "
          f">= {MIN_POS} positive cells on the larger side")
    print(f"      genes passing: {int(gate.sum())} / {len(names)}")
    ceil_flag = (g > CEIL) | (r > CEIL)
    print(f"      of those, touching the {CEIL:.0%} ceiling on one side: "
          f"{int((gate & ceil_flag).sum())}")

    idx = np.flatnonzero(gate)
    dNK = (g - r) * 100
    dT = (P["tg"] - P["tr"]) * 100
    tdiff = dNK - dT
    same_sign = np.sign(dNK) == np.sign(dT)
    cosyn = same_sign & (np.abs(dT) >= 0.5 * np.abs(dNK))

    recs = []
    for i in idx:
        gg, rr = g[i], r[i]
        rise = rr > gg
        f = f_est(gg, rr) if rise else f_min(gg, rr)
        rec = dict(gene=names[i], direction="rise" if rise else "fall",
                   green=gg, red=rr, delta_pp=dNK[i],
                   cd8t_green=P["tg"][i], cd8t_red=P["tr"][i],
                   cd8t_delta_pp=dT[i], T_diff_pp=tdiff[i],
                   cosynchronous_with_CD8T=bool(cosyn[i]),
                   n_pos_green=int(P["gp"][i]), n_pos_red=int(P["rp"][i]),
                   ceiling_touched=bool(ceil_flag[i]), f=f)
        for lbl in ("24h", "72h"):
            q = res[lbl]
            a2, b2 = q["g"][i], q["r"][i]
            rec[f"green_{lbl}"], rec[f"red_{lbl}"] = a2, b2
            rec[f"f_{lbl}"] = f_est(a2, b2) if b2 > a2 else f_min(a2, b2)
            rec[f"dir_{lbl}"] = "rise" if b2 > a2 else "fall"
        recs.append(rec)
    G = pd.DataFrame(recs)
    G["dir_stable"] = (G.dir_24h == G.dir_72h) & (G.dir_24h == G.direction)

    # ---------- 3. cell-level bootstrap ----------------------------
    print("\nbootstrap: resampling CELLS within colour x hour "
          f"({NBOOT} reps) -- SAMPLING uncertainty only, a LOWER BOUND")
    gsel = np.flatnonzero(nk & (colour == "Green"))
    rsel = np.flatnonzero(nk & (colour == "Red"))
    ghr, rhr = hours[gsel], hours[rsel]
    gsub = (X[gsel][:, idx] > 0).toarray()
    rsub = (X[rsel][:, idx] > 0).toarray()
    gstrat = [np.flatnonzero(ghr == h) for h in ("24", "72")]
    rstrat = [np.flatnonzero(rhr == h) for h in ("24", "72")]
    boot = np.empty((NBOOT, len(idx)), dtype=np.float32)
    for b in range(NBOOT):
        gi_ = np.concatenate([s[RNG.integers(0, len(s), len(s))] for s in gstrat])
        ri_ = np.concatenate([s[RNG.integers(0, len(s), len(s))] for s in rstrat])
        gb, rb = gsub[gi_].mean(axis=0), rsub[ri_].mean(axis=0)
        rise = rb > gb
        with np.errstate(divide="ignore", invalid="ignore"):
            boot[b] = np.where(rise, 1.0 - gb / rb, (gb - rb) / (1.0 - rb))
    G["f_lo"] = np.nanpercentile(boot, 2.5, axis=0)
    G["f_hi"] = np.nanpercentile(boot, 97.5, axis=0)
    G.to_csv(f"{OUT}/s3_clearance_per_gene.csv", index=False)

    # ---------- 4. §2.5 verdict ------------------------------------
    print("\n" + "-" * 78)
    print("§2.5 PRE-WRITTEN VERDICT")
    print("-" * 78)
    rise = G[(G.direction == "rise") & (~G.cosynchronous_with_CD8T)]
    rise_ex = G[(G.direction == "rise") & (G.cosynchronous_with_CD8T)]
    fall = G[G.direction == "fall"]
    print(f"rise-side genes usable for f_est: {len(rise)}   "
          f"excluded as co-synchronous with CD8T: {len(rise_ex)}")
    print(f"fall-side genes giving f_min:      {len(fall)}")
    q = rise.f.quantile([0.25, 0.5, 0.75])
    iqr = q[0.75] - q[0.25]
    print(f"\nf_est  median {q[0.5]:.3f}   IQR [{q[0.25]:.3f}, {q[0.75]:.3f}] "
          f"= {iqr:.3f}   range [{rise.f.min():.3f}, {rise.f.max():.3f}]")
    print(f"f_min  max {fall.f.max():.3f} ({fall.loc[fall.f.idxmax(),'gene']})   "
          f"median {fall.f.median():.3f}")
    above = int((fall.f > q[0.5]).sum())
    print(f"\nfall-side genes whose f_min EXCEEDS the f_est median: "
          f"{above} / {len(fall)}")
    if iqr > 0.25:
        v = "row 3: f_est is DISPERSED -> passive-rise explanation does not stand"
    elif q[0.5] >= fall.f.max():
        v = "row 1: a self-consistent clearance fraction exists"
    else:
        v = ("row 2: f_est is concentrated but BELOW some f_min -> the shortfall "
             "must be carried by state change")
    print(f"\nVERDICT -> {v}")
    print(f"(IQR rule: dispersed if IQR > 0.25)")

    with open(f"{OUT}/s3_verdict.txt", "w") as fh:
        fh.write(v + "\n")
    print("\nwrote s3_ppt_recheck.csv, s3_clearance_per_gene.csv, s3_verdict.txt")


if __name__ == "__main__":
    main()
