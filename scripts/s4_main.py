#!/usr/bin/env python3
"""
S-4 §2 (main), §3.2 (regression to the mean), §3.3 (expression-level
stratification) and §4 (where the chemokines sit).

Everything here is read against the §3.1 null, never against zero: the
null alone produces Spearman(p_green, drop) = -0.70 to -0.77, so the raw
sign of the correlation carries no information about targeting.

Gate is BASELINE-side only (p_green >= 0.05, >= 10 positive green cells).
Gating on the red side would be gating on the outcome and would delete
exactly the genes with the largest falls.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

H5 = "data/raw/DW_T_NK.h5ad"
OUT = "results"
TARGET_SUM = 4648.0
RNG = np.random.default_rng(20260821)
MIN_P, MIN_POS = 0.05, 10
NSIM = 25
NBOOT = 1000
NSPLIT = 200
CHEMO = ["Ccl3", "Ccl4", "Ccl5", "Xcl1", "Xcl2"]


def load():
    a = ad.read_h5ad(H5)
    R = a.raw
    names = np.array(R.var_names)
    X = R.X.tocsr()
    nc = a.obs.n_counts.values.astype(np.float64)
    obs = a.obs
    nk = (obs.main_celltype.astype(str) == "NK").values
    col = obs.colour.astype(str).values
    hrs = obs.hours.astype(str).values

    def cts(mask):
        i = np.flatnonzero(mask)
        return np.rint(np.expm1(X[i].toarray()) *
                       (nc[i][:, None] / TARGET_SUM)).astype(np.int32)

    gm, rm = nk & (col == "Green"), nk & (col == "Red")
    return names, cts(gm), cts(rm), hrs[gm], hrs[rm]


def thin(c, m, rng):
    o = np.zeros(m.shape, np.int32)
    nz = m > 0
    o[nz] = rng.binomial(m[nz], c)
    return o


def det(m):
    return (m > 0).mean(axis=0)


def null_drop(G, pg, c, nsim=NSIM, rng=RNG):
    acc = np.zeros((nsim, G.shape[1]))
    for b in range(nsim):
        acc[b] = det(thin(c, G, rng))
    return (pg - acc.mean(axis=0)) / pg


def main():
    names, G, R, hg, hr = load()
    pg, pr = det(G), det(R)
    npos = (G > 0).sum(axis=0)
    gate = (pg >= MIN_P) & (npos >= MIN_POS)
    idx = np.flatnonzero(gate)
    g, r = pg[idx], pr[idx]
    drop = (g - r) / g
    gn = names[idx]
    C_FIT = 0.87                       # level-matched, from §3.1

    print("=" * 78)
    print("S-4 §2  MAIN ANALYSIS")
    print("=" * 78)
    print(f"genes through the baseline-side gate: {len(idx)}")
    print(f"drop ratio: median {np.median(drop):+.4f}  "
          f"IQR [{np.percentile(drop,25):+.4f}, {np.percentile(drop,75):+.4f}]  "
          f"range [{drop.min():+.2f}, {drop.max():+.2f}]")
    print(f"genes that FALL (drop > 0): {int((drop>0).sum())} ({(drop>0).mean():.1%})   "
          f"RISE (drop < 0): {int((drop<0).sum())} ({(drop<0).mean():.1%})")
    rho = stats.spearmanr(g, drop)
    print(f"\nSpearman(p_green, drop ratio) = {rho.statistic:+.4f}  p = {rho.pvalue:.3g}")

    # ---- §2.3 recheck of the two named points --------------------
    print("\n--- §2.3 re-check of the two quoted points ---")
    print(f"{'gene':6s} {'§2.3 drop':>10} {'measured green':>15} {'measured red':>13} "
          f"{'measured drop':>14} {'null drop':>10} {'residual':>9}")
    dn = null_drop(G, pg, C_FIT)[idx]
    quoted = {"Ccl3": (12.59, 2.80), "Ccl5": (90.62, 36.52)}
    for gene, (a_, b_) in quoted.items():
        k = int(np.where(gn == gene)[0][0])
        qd = (a_ - b_) / a_
        print(f"{gene:6s} {qd:10.3f} {g[k]:15.4f} {r[k]:13.4f} "
              f"{drop[k]:14.4f} {dn[k]:10.4f} {drop[k]-dn[k]:9.4f}")

    # ---- bootstrap CIs -------------------------------------------
    print("\n--- bootstrap intervals on Spearman(p_green, drop) ---")
    bg = np.empty(NBOOT)
    for b in range(NBOOT):
        s = RNG.integers(0, len(idx), len(idx))
        bg[b] = stats.spearmanr(g[s], drop[s]).statistic
    print(f"  gene bootstrap  (n={NBOOT}): {np.median(bg):+.4f} "
          f"[{np.percentile(bg,2.5):+.4f}, {np.percentile(bg,97.5):+.4f}]")

    gs = [np.flatnonzero(hg == h) for h in ("24", "72")]
    rs = [np.flatnonzero(hr == h) for h in ("24", "72")]
    bc = np.empty(300)
    for b in range(300):
        gi = np.concatenate([s[RNG.integers(0, len(s), len(s))] for s in gs])
        ri = np.concatenate([s[RNG.integers(0, len(s), len(s))] for s in rs])
        pgb, prb = det(G[gi])[idx], det(R[ri])[idx]
        ok = pgb > 0
        bc[b] = stats.spearmanr(pgb[ok], ((pgb - prb) / pgb)[ok]).statistic
    print(f"  cell bootstrap  (n=300): {np.median(bc):+.4f} "
          f"[{np.percentile(bc,2.5):+.4f}, {np.percentile(bc,97.5):+.4f}]")
    print("  (cell bootstrap carries sampling noise only; colour is confounded")
    print("   with library, so library/animal variance is invisible here)")

    resid = drop - dn
    rr = stats.spearmanr(g, resid)
    print(f"\nSpearman(p_green, RESIDUAL vs null c={C_FIT}) = "
          f"{rr.statistic:+.4f}  p = {rr.pvalue:.3g}")
    br = np.empty(NBOOT)
    for b in range(NBOOT):
        s = RNG.integers(0, len(idx), len(idx))
        br[b] = stats.spearmanr(g[s], resid[s]).statistic
    print(f"  gene bootstrap: {np.median(br):+.4f} "
          f"[{np.percentile(br,2.5):+.4f}, {np.percentile(br,97.5):+.4f}]")

    # ---- §3.2 regression to the mean -----------------------------
    print("\n" + "=" * 78)
    print("S-4 §3.2  REGRESSION TO THE MEAN (split-half green)")
    print("=" * 78)
    print(f"green cells split at random within each timepoint, {NSPLIT} repeats;")
    print("baseline from half A, drop computed from half B -- the baseline")
    print("estimate and the drop numerator then share no cells.")
    raw_s, res_s, naive_s = [], [], []
    for _ in range(NSPLIT):
        A, B = [], []
        for s in gs:
            p = RNG.permutation(s)
            A.append(p[: len(p) // 2]); B.append(p[len(p) // 2:])
        A, B = np.concatenate(A), np.concatenate(B)
        pa, pb = det(G[A])[idx], det(G[B])[idx]
        ok = (pa > 0) & (pb > 0)
        dB = (pb[ok] - r[ok]) / pb[ok]
        raw_s.append(stats.spearmanr(pa[ok], dB).statistic)
        naive_s.append(stats.spearmanr(pb[ok], dB).statistic)
        dnB = null_drop(G[np.ix_(B, idx)], pb, C_FIT, nsim=5)[ok]
        res_s.append(stats.spearmanr(pa[ok], dB - dnB).statistic)
    for lab, v in (("naive (same half both sides)", naive_s),
                   ("split-half, raw drop", raw_s),
                   ("split-half, residual vs null", res_s)):
        v = np.array(v)
        print(f"  {lab:32s} {v.mean():+.4f}  "
              f"[{np.percentile(v,2.5):+.4f}, {np.percentile(v,97.5):+.4f}]")

    # ---- §3.3 stratify by expression level -----------------------
    print("\n" + "=" * 78)
    print("S-4 §3.3  STRATIFIED BY MEAN EXPRESSION")
    print("=" * 78)
    meanexp = G[:, idx].mean(axis=0)
    lm = np.log10(meanexp + 1e-6)
    q = np.quantile(lm, np.linspace(0, 1, 11))
    print(f"{'decile':>7} {'n':>6} {'mean UMI/cell':>14} {'p_green range':>18} "
          f"{'rho(raw)':>10} {'rho(resid)':>12} {'p(resid)':>10}")
    srows = []
    for k in range(10):
        m = (lm >= q[k]) & (lm <= q[k + 1] if k == 9 else lm < q[k + 1])
        if m.sum() < 30:
            continue
        a1 = stats.spearmanr(g[m], drop[m])
        a2 = stats.spearmanr(g[m], resid[m])
        print(f"{k+1:7d} {int(m.sum()):6d} {meanexp[m].mean():14.4f} "
              f"{g[m].min():7.3f}-{g[m].max():<9.3f} "
              f"{a1.statistic:10.4f} {a2.statistic:12.4f} {a2.pvalue:10.2g}")
        srows.append(dict(decile=k + 1, n=int(m.sum()),
                          mean_umi=float(meanexp[m].mean()),
                          pg_lo=float(g[m].min()), pg_hi=float(g[m].max()),
                          rho_raw=float(a1.statistic),
                          rho_resid=float(a2.statistic), p_resid=float(a2.pvalue)))
    S = pd.DataFrame(srows); S.to_csv(f"{OUT}/s4_strata.csv", index=False)
    print(f"\n  rho(resid) across deciles: median {S.rho_resid.median():+.4f}, "
          f"range [{S.rho_resid.min():+.4f}, {S.rho_resid.max():+.4f}], "
          f"negative in {int((S.rho_resid<0).sum())}/{len(S)}")

    # ---- §4 chemokine position -----------------------------------
    print("\n" + "=" * 78)
    print("S-4 §4  WHERE THE CHEMOKINES SIT")
    print("=" * 78)
    P = pd.DataFrame(dict(gene=gn, p_green=g, p_red=r, drop_ratio=drop,
                          null_drop=dn, resid=resid, mean_umi=meanexp,
                          n_pos_green=npos[idx]))
    P["resid_pct"] = P["resid"].rank(pct=True)
    P.to_csv(f"{OUT}/s4_per_gene_full.csv", index=False)
    P["drop_pct"] = P["drop_ratio"].rank(pct=True)
    P.to_csv(f"{OUT}/s4_per_gene_full.csv", index=False)
    print(f"{'gene':7s} {'p_green':>8} {'p_red':>7} {'drop':>8} {'null':>8} "
          f"{'residual':>9} {'resid pct':>10} {'drop pct':>9}")
    for gene in CHEMO:
        row = P[P.gene == gene]
        if not len(row):
            print(f"{gene:7s}  -- did not pass the baseline gate "
                  f"(p_green = {pg[int(np.where(names==gene)[0][0])]:.4f})")
            continue
        x = row.iloc[0]
        print(f"{gene:7s} {x.p_green:8.4f} {x.p_red:7.4f} {x.drop_ratio:8.4f} "
              f"{x.null_drop:8.4f} {x.resid:9.4f} {x.resid_pct:10.4f} "
              f"{x.drop_pct:9.4f}")

    # local comparison: genes with similar baseline
    print("\n  --- compared with genes of similar baseline (+/- 0.05 in p_green) ---")
    for gene in CHEMO:
        row = P[P.gene == gene]
        if not len(row):
            continue
        x = row.iloc[0]
        nb = P[(P.p_green - x.p_green).abs() <= 0.05]
        pctl = (nb.resid < x.resid).mean()
        print(f"  {gene:7s} n_neighbours={len(nb):5d}  neighbour residual median "
              f"{nb.resid.median():+.4f}  {gene} residual {x.resid:+.4f}  "
              f"-> percentile {pctl:.4f}"
              + ("   <== ABOVE 99th" if pctl > .99 else
                 "   <== BELOW 1st" if pctl < .01 else ""))
    print("\nwrote s4_per_gene_full.csv, s4_strata.csv")


if __name__ == "__main__":
    main()
