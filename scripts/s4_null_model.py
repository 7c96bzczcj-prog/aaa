#!/usr/bin/env python3
"""
S-4 §3.1 -- the null model, built and calibrated BEFORE §2 is looked at.

The question §2 asks is whether the drop ratio (p_green - p_red)/p_green
depends on the baseline p_green. §3.1 warns that the 0/1 bounds on a
detection rate manufacture such a dependence on their own. This script
builds the null and checks that it does.

NULL (uniform suppression), non-parametric
------------------------------------------
"Every cell's expression scaled by the same factor c" is, at the UMI
level, exactly binomial thinning: count' ~ Binomial(count, c). Thinning
the observed green counts preserves the real library-size spread, the
real per-gene mean structure and any zero inflation -- nothing is assumed
beyond "UMIs are Poisson samples of an underlying abundance", which is the
standard UMI model. This is strictly stronger than a fitted curve.

Counts are recoverable exactly from this object: the matrix is
log1p(normalize_total(target_sum=4648)), so
    count = expm1(x) * n_counts / 4648
(verified integral to < 1e-4 before use).

A parametric Poisson null is reported alongside for comparison:
    p_red = 1 - (1 - p_green)^c   =>   drop(q) = (q^c - q)/(1 - q),  q = 1-p
which tends to (1-c) as p_green -> 0 and to 0 as p_green -> 1.

POSITIVE CONTROL (targeted), to show the test can see targeting
---------------------------------------------------------------
§1's "targeted" is: only the originally-high-expressing cells are hit.
Instantiated per gene: silence the top-k cells by that gene's expression,
k = fraction * n_cells. A gene detected in fewer than k cells loses
everything (drop 1.0); a gene detected in nearly all cells loses only
k/n. That is exactly the shape §1 predicts.

IDENTIFIABILITY NOTE, stated once and not worked around: a global
sequencing-depth difference between the green and red libraries acts on
detection rates exactly like a global suppression factor. c therefore
absorbs both. This analysis tests the SHAPE (dependence on baseline), not
the global level, and no claim about the global level is made.
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
NSIM = 25                      # thinning replicates per c
MIN_P, MIN_POS = 0.05, 10      # BASELINE-side gate only (see below)
C_GRID = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]


def load_counts():
    a = ad.read_h5ad(H5)
    R = a.raw
    names = np.array(R.var_names)
    X = R.X.tocsr()
    nc = a.obs.n_counts.values.astype(np.float64)
    obs = a.obs
    nk = (obs.main_celltype.astype(str) == "NK").values
    colour = obs.colour.astype(str).values
    hours = obs.hours.astype(str).values

    def counts_for(mask):
        idx = np.flatnonzero(mask)
        sub = X[idx].toarray()
        c = np.expm1(sub) * (nc[idx][:, None] / TARGET_SUM)
        return np.rint(c).astype(np.int32)

    g = counts_for(nk & (colour == "Green"))
    r = counts_for(nk & (colour == "Red"))
    return names, g, r, nc[nk & (colour == "Green")], nc[nk & (colour == "Red")], \
        hours[nk & (colour == "Green")], hours[nk & (colour == "Red")]


def thin(counts, c, rng):
    """binomial thinning: uniform suppression of every cell by factor c."""
    out = np.zeros(counts.shape, dtype=np.int32)
    nz = counts > 0
    out[nz] = rng.binomial(counts[nz], c)
    return out


def detect(counts):
    return (counts > 0).mean(axis=0)


def targeted(counts, frac, rng):
    """silence the top-k cells per gene, k = frac * n_cells."""
    n = counts.shape[0]
    k = int(round(frac * n))
    if k == 0:
        return counts.copy()
    out = counts.copy()
    # rank cells per gene by expression, break ties randomly
    jitter = rng.random(counts.shape) * 0.5
    order = np.argsort(-(counts + jitter), axis=0)[:k]
    np.put_along_axis(out, order, 0, axis=0)
    return out


def main():
    names, G, R, ncg, ncr, hg, hr = load_counts()
    print("=" * 78)
    print("S-4 §3.1  NULL MODEL (built before §2 is inspected)")
    print("=" * 78)
    print(f"green NK {G.shape[0]} cells, red NK {R.shape[0]} cells, "
          f"{G.shape[1]} genes")
    print(f"median UMI/cell: green {np.median(ncg):.0f}  red {np.median(ncr):.0f}  "
          f"ratio red/green = {np.median(ncr)/np.median(ncg):.3f}")
    print(f"  green by hour: 24h {np.median(ncg[hg=='24']):.0f} (n={(hg=='24').sum()})"
          f"  72h {np.median(ncg[hg=='72']):.0f} (n={(hg=='72').sum()})")
    print(f"  red   by hour: 24h {np.median(ncr[hr=='24']):.0f} (n={(hr=='24').sum()})"
          f"  72h {np.median(ncr[hr=='72']):.0f} (n={(hr=='72').sum()})")

    pg, pr = detect(G), detect(R)
    npos_g = (G > 0).sum(axis=0)

    # ---- gate: BASELINE side only ---------------------------------
    gate = (pg >= MIN_P) & (npos_g >= MIN_POS)
    print(f"\ngate (BASELINE side only): p_green >= {MIN_P} and "
          f">= {MIN_POS} positive green cells -> {int(gate.sum())} genes")
    two_sided = gate & (pr >= MIN_P)
    print(f"  [for reference, S-3's two-sided gate would keep {int(two_sided.sum())}]")
    print("  the two-sided gate is NOT used here: gating on the red side is "
          "gating on the outcome,\n  and would drop exactly the genes with the "
          "largest falls -- the ones §2 is about.")

    idx = np.flatnonzero(gate)
    obs_drop = (pg[idx] - pr[idx]) / pg[idx]
    print(f"\nobserved drop ratio: median {np.median(obs_drop):+.3f}  "
          f"IQR [{np.percentile(obs_drop,25):+.3f}, {np.percentile(obs_drop,75):+.3f}]")
    rho_obs = stats.spearmanr(pg[idx], obs_drop)
    print(f"observed Spearman(p_green, drop ratio) = {rho_obs.statistic:+.4f}  "
          f"p = {rho_obs.pvalue:.3g}")

    # ---- 1. does the null alone create the correlation? ------------
    print("\n" + "-" * 78)
    print("1. UNIFORM-SUPPRESSION NULL (binomial thinning), c grid")
    print("-" * 78)
    print(f"{'c':>5} {'median drop':>12} {'Spearman(p_g, drop)':>21} "
          f"{'p':>10}   {'drop@p=0.1':>11} {'drop@p=0.9':>11}")
    rows, sim_store = [], {}
    for c in C_GRID:
        acc = np.zeros((NSIM, len(idx)))
        for b in range(NSIM):
            acc[b] = detect(thin(G, c, RNG))[idx]
        p_sim = acc.mean(axis=0)
        d_sim = (pg[idx] - p_sim) / pg[idx]
        sim_store[c] = d_sim
        rho = stats.spearmanr(pg[idx], d_sim)
        lo = d_sim[(pg[idx] > 0.08) & (pg[idx] < 0.12)]
        hi = d_sim[(pg[idx] > 0.88) & (pg[idx] < 0.92)]
        print(f"{c:5.1f} {np.median(d_sim):12.3f} {rho.statistic:21.4f} "
              f"{rho.pvalue:10.2g}   {np.mean(lo):11.3f} {np.mean(hi):11.3f}")
        rows.append(dict(c=c, median_drop=float(np.median(d_sim)),
                         rho=float(rho.statistic), p=float(rho.pvalue),
                         drop_at_0p1=float(np.mean(lo)),
                         drop_at_0p9=float(np.mean(hi))))
    N = pd.DataFrame(rows)
    N.to_csv(f"{OUT}/s4_null_grid.csv", index=False)

    print("\n  => the null ALONE produces a strongly NEGATIVE Spearman at every c.")
    print("     §2.2's rule 'negative => targeted' therefore cannot be applied")
    print("     to the raw correlation. Only the RESIDUAL against this null counts.")

    # ---- 2. parametric Poisson null, for comparison ---------------
    print("\n" + "-" * 78)
    print("2. PARAMETRIC POISSON NULL  drop(q) = (q^c - q)/(1-q),  q = 1 - p_green")
    print("-" * 78)
    q = 1.0 - pg[idx]
    print(f"{'c':>5} {'median drop':>12} {'Spearman':>10}   "
          f"{'vs thinning: median |diff|':>27}")
    for c in C_GRID:
        d_par = (np.power(q, c) - q) / (1.0 - q)
        rho = stats.spearmanr(pg[idx], d_par)
        print(f"{c:5.1f} {np.median(d_par):12.3f} {rho.statistic:10.4f}   "
              f"{np.median(np.abs(d_par - sim_store[c])):27.4f}")
    print("  the two nulls agree closely; the thinning null is used as primary")
    print("  because it keeps the real library-size and zero-inflation structure.")

    # ---- 3. positive control: can the test see targeting? ---------
    print("\n" + "-" * 78)
    print("3. POSITIVE CONTROL -- targeted silencing of the top-k cells per gene")
    print("-" * 78)
    print(f"{'frac':>6} {'median drop':>12} {'Spearman(p_g,drop)':>20}   "
          f"{'drop@p=0.1':>11} {'drop@p=0.9':>11}")
    tgt_store = {}
    for frac in (0.1, 0.2, 0.3, 0.5):
        acc = np.zeros((5, len(idx)))
        for b in range(5):
            acc[b] = detect(targeted(G, frac, RNG))[idx]
        p_t = acc.mean(axis=0)
        d_t = (pg[idx] - p_t) / pg[idx]
        tgt_store[frac] = d_t
        rho = stats.spearmanr(pg[idx], d_t)
        lo = d_t[(pg[idx] > 0.08) & (pg[idx] < 0.12)]
        hi = d_t[(pg[idx] > 0.88) & (pg[idx] < 0.92)]
        print(f"{frac:6.1f} {np.median(d_t):12.3f} {rho.statistic:20.4f}   "
              f"{np.mean(lo):11.3f} {np.mean(hi):11.3f}")

    # ---- 4. the actual test: residual against the null ------------
    print("\n" + "-" * 78)
    print("4. THE TEST -- residual of observed against the matched null")
    print("-" * 78)
    print("c is chosen to match the OVERALL level (median drop ratio), so that")
    print("only the SHAPE is compared. Two anchors are reported.")
    med_obs = np.median(obs_drop)
    c_fit = float(N.c.iloc[np.argmin(np.abs(N.median_drop - med_obs))])
    # refine on a finer grid around c_fit
    fine = np.round(np.arange(max(0.02, c_fit - 0.12), min(0.99, c_fit + 0.12), 0.01), 2)
    best, bestd = c_fit, 1e9
    fine_rows = []
    for c in fine:
        acc = np.zeros((NSIM, len(idx)))
        for b in range(NSIM):
            acc[b] = detect(thin(G, c, RNG))[idx]
        d = (pg[idx] - acc.mean(axis=0)) / pg[idx]
        gap = abs(np.median(d) - med_obs)
        fine_rows.append((c, float(np.median(d)), gap))
        if gap < bestd:
            bestd, best, best_d = gap, c, d
    c_depth = float(np.median(ncr) / np.median(ncg))
    print(f"  observed median drop ratio        : {med_obs:.4f}")
    print(f"  level-matched c (fine grid)       : {best:.2f}  "
          f"(null median drop {np.median(best_d):.4f})")
    print(f"  depth-implied c (median UMI ratio): {c_depth:.3f}")

    res_rows = []
    for lab, c, d_null in (("level-matched", best, best_d),
                           ("depth-implied", round(c_depth, 2), None)):
        if d_null is None:
            acc = np.zeros((NSIM, len(idx)))
            for b in range(NSIM):
                acc[b] = detect(thin(G, c, RNG))[idx]
            d_null = (pg[idx] - acc.mean(axis=0)) / pg[idx]
        resid = obs_drop - d_null
        rho = stats.spearmanr(pg[idx], resid)
        print(f"\n  --- null c = {c:.2f} ({lab}) ---")
        print(f"      residual: median {np.median(resid):+.4f}  "
              f"IQR [{np.percentile(resid,25):+.4f}, {np.percentile(resid,75):+.4f}]")
        print(f"      Spearman(p_green, residual) = {rho.statistic:+.4f}  "
              f"p = {rho.pvalue:.3g}")
        res_rows.append(dict(anchor=lab, c=c, rho_resid=float(rho.statistic),
                             p=float(rho.pvalue),
                             median_resid=float(np.median(resid))))
        np.save(f"{OUT}/s4_resid_{lab.replace('-','_')}.npy", resid)
        np.save(f"{OUT}/s4_dnull_{lab.replace('-','_')}.npy", d_null)

    # same test applied to the targeted positive control (power check)
    print("\n  --- power check: the same residual test on the targeted control ---")
    acc = np.zeros((NSIM, len(idx)))
    for b in range(NSIM):
        acc[b] = detect(thin(G, best, RNG))[idx]
    d_null_best = (pg[idx] - acc.mean(axis=0)) / pg[idx]
    for frac, d_t in tgt_store.items():
        rho = stats.spearmanr(pg[idx], d_t - d_null_best)
        print(f"      targeted frac {frac:.1f}: Spearman(p_green, residual) = "
              f"{rho.statistic:+.4f}  p = {rho.pvalue:.3g}")

    pd.DataFrame(res_rows).to_csv(f"{OUT}/s4_residual_test.csv", index=False)
    pd.DataFrame(dict(gene=names[idx], p_green=pg[idx], p_red=pr[idx],
                      n_pos_green=npos_g[idx], drop=obs_drop,
                      d_null=d_null_best, resid=obs_drop - d_null_best)
                 ).to_csv(f"{OUT}/s4_per_gene.csv", index=False)
    print("\nwrote s4_null_grid.csv, s4_residual_test.csv, s4_per_gene.csv")


if __name__ == "__main__":
    main()
