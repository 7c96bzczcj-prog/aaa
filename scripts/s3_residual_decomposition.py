#!/usr/bin/env python3
"""
S-3 §4 -- weighted-matching residual decomposition.

LIMIT THAT DOES NOT GO AWAY (§7): green/red cannot separate clearance from
reprogramming. This asks only "how much of the green->red difference can
ANY removal of green cells account for", never "it is clearance".

Question, stated exactly: is there a subset of the green cells whose
removal leaves a remainder matching red on the target genes? Removing a
subset is relaxed to weighting the green cells (continuous weights on the
simplex). The relaxation is strictly more permissive than hard subsetting,
so every residual reported here is a LOWER BOUND on the residual under
literal subset removal.

Two solvers, both on the detection (0/1) matrix:

  A. Entropy balancing (Hainmueller). Maximum-entropy weights subject to
     EXACT moment constraints, solved in the dual:
         min_lambda  log sum_i exp(-lambda . x_i) + lambda . t
     Convergence with small residual => a clearance pattern exists that
     reproduces red exactly on those genes. Divergence => infeasible.
     Used for the small chemokine/residency gene set.

  B. Best-achievable match. w = softmax(theta), minimise
         sum_j ( (X' w)_j - t_j )^2
     over the simplex. This is the most gap ANY reweighting can close.
     Used for the full gated gene set, where exact balancing on thousands
     of moments with a few hundred cells is not a well-posed problem.

Both report the Kish effective sample size ESS = (sum w)^2 / sum w^2.
Stated rule (a choice of mine, not a pre-registered number): ESS < 30 ->
UNRELIABLE.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

H5 = "data/raw/DW_T_NK.h5ad"
OUT = "results"
ESS_MIN = 30
CHEMO = ["Ccl3", "Ccl4", "Ccl5", "Xcl1", "Xcl2"]
RESID = ["Itga1", "Itgae", "Cd69", "Cxcr6", "Rgs1", "Zfp683", "Cd160"]


def kish(w):
    return float(w.sum() ** 2 / (w ** 2).sum())


def entropy_balance(Xg, t):
    """max-entropy weights on green cells with exact moments t."""
    def obj(lam):
        z = -Xg @ lam
        z -= z.max()
        e = np.exp(z)
        return np.log(e.sum()) + z.max() + lam @ t

    def grad(lam):
        z = -Xg @ lam
        z -= z.max()
        w = np.exp(z)
        w /= w.sum()
        return t - Xg.T @ w

    r = minimize(obj, np.zeros(Xg.shape[1]), jac=grad, method="L-BFGS-B",
                 options=dict(maxiter=5000, ftol=1e-14, gtol=1e-12))
    z = -Xg @ r.x
    z -= z.max()
    w = np.exp(z)
    w /= w.sum()
    return w, r


def best_match(Xg, t, iters=4000):
    """w = softmax(theta); minimise squared moment error over the simplex."""
    n = Xg.shape[0]

    def obj(th):
        th = th - th.max()
        e = np.exp(th)
        w = e / e.sum()
        r = Xg.T @ w - t
        return float(r @ r)

    def grad(th):
        th = th - th.max()
        e = np.exp(th)
        s = e.sum()
        w = e / s
        r = Xg.T @ w - t
        g = 2.0 * (Xg @ r)                    # d/dw
        return w * (g - float(g @ w))         # softmax Jacobian

    res = minimize(obj, np.zeros(n), jac=grad, method="L-BFGS-B",
                   options=dict(maxiter=iters, ftol=1e-16, gtol=1e-14))
    th = res.x - res.x.max()
    e = np.exp(th)
    w = e / e.sum()
    return w, res


def report(name, w, Xg, green, t, genes, res_obj, fh):
    ach = Xg.T @ w
    d0 = green - t                 # original gap
    d1 = ach - t                   # residual gap
    l1 = np.abs(d1).sum() / np.abs(d0).sum()
    l2 = float(np.sqrt((d1 ** 2).sum() / (d0 ** 2).sum()))
    ess = kish(w)
    ok = ess >= ESS_MIN
    lines = [
        f"\n--- {name} ---",
        f"  genes matched            : {len(genes)}",
        f"  green cells              : {Xg.shape[0]}",
        f"  Kish ESS                 : {ess:.1f}  "
        f"({'OK' if ok else f'< {ESS_MIN} -> UNRELIABLE'})",
        f"  max |w| / mean |w|       : {w.max()/w.mean():.1f}",
        f"  residual share (L1)      : {l1:.3f}   "
        f"=> clearance can close {1-l1:.1%} of the gap",
        f"  residual share (L2)      : {l2:.3f}",
        f"  optimiser converged      : {res_obj.success} ({res_obj.message})",
    ]
    for ln in lines:
        print(ln)
        fh.write(ln + "\n")
    return dict(gene_set=name, n_genes=len(genes), n_green=Xg.shape[0],
                ess=ess, reliable=ok, residual_L1=l1, residual_L2=l2,
                closed_share=1 - l1), ach


def main():
    a = ad.read_h5ad(H5)
    R = a.raw
    names = np.array(R.var_names)
    X = R.X.tocsr()
    obs = a.obs
    nk = (obs.main_celltype.astype(str) == "NK").values
    colour = obs.colour.astype(str).values
    g_idx = np.flatnonzero(nk & (colour == "Green"))
    r_idx = np.flatnonzero(nk & (colour == "Red"))

    gated = pd.read_csv(f"{OUT}/s3_clearance_per_gene.csv").gene.tolist()
    pos = {g: i for i, g in enumerate(names)}
    setA = [g for g in gated if g in pos]
    setB = [g for g in CHEMO + RESID if g in set(gated)]
    print("=" * 78)
    print("S-3 §4  WEIGHTED-MATCHING RESIDUAL DECOMPOSITION")
    print("=" * 78)
    print(f"green NK {len(g_idx)}   red NK {len(r_idx)}")
    print(f"gene set A (all gated)          : {len(setA)}")
    print(f"gene set B (chemokine+residency): {len(setB)}  {setB}")
    missing = [g for g in CHEMO + RESID if g not in set(gated)]
    print(f"  of the named markers, not gated (excluded): {missing}")

    rows, ach_store = [], {}
    with open(f"{OUT}/s3_residual.log", "w") as fh:
        fh.write(f"green NK {len(g_idx)} red NK {len(r_idx)}\n")
        for nm, gl, solver in (("B: chemokine+residency (entropy balancing)",
                                setB, "eb"),
                               ("B: chemokine+residency (best match)",
                                setB, "bm"),
                               ("A: all gated genes (best match)",
                                setA, "bm")):
            cols = [pos[g] for g in gl]
            Xg = np.asarray((X[g_idx][:, cols] > 0).todense(), dtype=float)
            Xr = np.asarray((X[r_idx][:, cols] > 0).todense(), dtype=float)
            green, t = Xg.mean(axis=0), Xr.mean(axis=0)
            if solver == "eb":
                w, res = entropy_balance(Xg, t)
            else:
                w, res = best_match(Xg, t)
            rec, ach = report(nm, w, Xg, green, t, gl, res, fh)
            rows.append(rec)
            ach_store[nm] = (gl, green, t, ach)

    pd.DataFrame(rows).to_csv(f"{OUT}/s3_residual_summary.csv", index=False)

    print("\n" + "-" * 78)
    print("PER-GENE RESIDUAL on the chemokine / residency set (best match)")
    print("-" * 78)
    gl, green, t, ach = ach_store["B: chemokine+residency (best match)"]
    per = pd.DataFrame(dict(gene=gl, green=green, red=t, matched=ach))
    per["gap_pp"] = (per.green - per.red) * 100
    per["residual_pp"] = (per.matched - per.red) * 100
    per["closed"] = 1 - np.abs(per.residual_pp) / np.abs(per.gap_pp)
    print(per.round(3).to_string(index=False))
    per.to_csv(f"{OUT}/s3_residual_per_gene.csv", index=False)

    glA, greenA, tA, achA = ach_store["A: all gated genes (best match)"]
    perA = pd.DataFrame(dict(gene=glA, green=greenA, red=tA, matched=achA))
    perA["gap_pp"] = (perA.green - perA.red) * 100
    perA["residual_pp"] = (perA.matched - perA.red) * 100
    perA.to_csv(f"{OUT}/s3_residual_per_gene_all.csv", index=False)
    key = perA[perA.gene.isin(CHEMO + RESID)]
    print("\nthe same genes inside the FULL-gene-set match:")
    print(key.round(3).to_string(index=False))
    print("\nwrote s3_residual_summary.csv, s3_residual_per_gene.csv, "
          "s3_residual_per_gene_all.csv, s3_residual.log")


if __name__ == "__main__":
    main()
