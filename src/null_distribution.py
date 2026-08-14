#!/usr/bin/env python
"""T4 null_distribution.tsv — an empirical null from expression-matched genes.

R6.5: >= 400 random genes matched to the panel on expression level, pushed
through the identical donor-level pipeline. The resulting distribution of
effect sizes is the baseline against which target effects are read.

The baseline is NOT assumed to sit at zero. Earlier work in this project
measured baselines ranging from -0.002 +/- 0.037 to -0.184 +/- 0.252 depending
on which gene set was used, which is exactly why it has to be measured.

    python src/null_distribution.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import (as_csr, cell_qc, detection_and_cpm,  # noqa: E402
                            downsample_marginal)
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402
from dnkchem.stats import (benjamini_hochberg, continuity_rate,  # noqa: E402
                           exact_wilcoxon_signed_rank, logit)

DNK_TRIO = ["dNK1", "dNK2", "dNK3"]


def match_expression(mean_expr, target_cols, n_per_target, rng, exclude, window=0):
    """Pick genes whose pooled expression brackets each target gene.

    The window is auto-sized when not given: a 10,000-gene null needs a wider
    rank neighbourhood than a 400-gene one, and silently returning fewer genes
    would cap the empirical-p resolution without saying so.
    """
    order = np.argsort(mean_expr)
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    if window <= 0:
        window = max(250, int(n_per_target * 3))
    # Null genes are SHARED between targets, not consumed. Banning a gene once
    # it is used starves whichever targets are processed last -- and starves
    # them worst exactly where the neighbourhood is thin, i.e. at the extremes
    # of expression, which is where the headline genes live. Each target's null
    # set is evaluated on its own, so overlap between sets is harmless.
    banned = set(exclude)
    per_target = {}
    for c in target_cols:
        r = rank[c]
        lo, hi = max(0, r - window), min(len(order), r + window)
        pool = [int(order[k]) for k in range(lo, hi)
                if int(order[k]) not in banned and mean_expr[int(order[k])] > 0]
        rng.shuffle(pool)
        per_target[int(c)] = pool[:n_per_target]
    return per_target


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--n-null", type=int, default=10000)
    ap.add_argument("--match-window", type=int, default=0,
                    help="rank half-width for expression matching; "
                         "0 = auto-size to satisfy --n-null")
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    meta = json.load(open(os.path.join(outdir, "detection_meta.json")))
    depth = meta["primary_depth_floor"]

    ds = load_dataset(mf)
    hits, _, _, _ = ds.panel_columns(panel)
    X = as_csr(ds.X)
    keep, _ = cell_qc(X, ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs

    sel = keep & obs["subset"].isin(DNK_TRIO).to_numpy() & (obs["compartment"] == "Decidua").to_numpy()
    rows_sel = np.arange(len(obs))[sel]
    Xs = X[rows_sel]
    obs_s = obs.loc[sel]
    print(f"[null] decidual dNK1-3 cells after QC: {Xs.shape[0]}")

    pooled = np.asarray(Xs.sum(axis=0)).ravel()
    mean_expr = pooled / max(1.0, pooled.sum()) * 1e6

    target_cols = [hits[g] for g in panel[panel.role == "primary_target"]["gene_symbol"]
                   if g in hits]
    n_per = max(1, int(np.ceil(args.n_null / max(1, len(target_cols)))))
    rng = np.random.default_rng(args.seed)
    per_target = match_expression(mean_expr, target_cols, n_per, rng,
                                  exclude=set(hits.values()),
                                  window=args.match_window)
    null_cols = sorted({c for v in per_target.values() for c in v})
    col_to_gene = {hits[g]: g for g in hits}
    # each target gets its OWN expression-matched null set: a gene at 2300 CPM
    # must not be judged against a null pool with median 5 CPM, or the null's
    # variance is understated and every z is inflated.
    gene_to_nulls = {col_to_gene[c]: v for c, v in per_target.items() if c in col_to_gene}
    print(f"[null] {len(null_cols)} expression-matched null genes "
          f"({n_per} per target, target CPM range "
          f"{mean_expr[target_cols].min():.2f}-{mean_expr[target_cols].max():.2f})")
    if len(null_cols) < args.n_null:
        print(f"[null] WARNING: only {len(null_cols)} null genes available "
              f"(requested {args.n_null}); empirical-p resolution is "
              f"1/{len(null_cols)+1} = {1/(len(null_cols)+1):.5f}")
    # expression matching quality, reported rather than assumed
    tq = np.quantile(mean_expr[target_cols], [0.1, 0.5, 0.9])
    nq = np.quantile(mean_expr[null_cols], [0.1, 0.5, 0.9])
    print(f"[null] CPM q10/q50/q90 -- targets {tq.round(3)} | nulls {nq.round(3)}")
    sizes = {col_to_gene[c]: len(v) for c, v in per_target.items() if c in col_to_gene}
    thin = {g: n for g, n in sizes.items() if n < n_per}
    print(f"[null] per-target null-set size: min {min(sizes.values())}, "
          f"median {int(np.median(list(sizes.values())))}, max {max(sizes.values())}")
    if thin:
        print(f"[null] targets with a THIN null set (resolution 1/(n+1) is worse "
              f"for these): {dict(sorted(thin.items(), key=lambda kv: kv[1])[:8])}")

    # per donor x subset detection rates for the null genes, same depth, same seed
    recs = []
    for (donor, ss), r in unit_indices(obs_s, np.ones(len(obs_s), bool), ["donor", "subset"]):
        counts, kept = downsample_marginal(Xs[r], null_cols, depth, seed=args.seed)
        n_post = int(kept.sum())
        if n_post < args.min_cells:
            continue
        n_det, det, _ = detection_and_cpm(counts, depth)
        for j in range(len(null_cols)):
            recs.append({"donor": donor, "subset": ss, "col": null_cols[j],
                         "n_cells": n_post, "n_detected": int(n_det[j]),
                         "detection_rate": float(det[j])})
    D = pd.DataFrame(recs)
    if D.empty:
        raise SystemExit("[null] no admissible units")

    T2 = pd.read_csv(os.path.join(outdir, "donor_level_tests.tsv"), sep="\t")
    out_rows, summary_rows = [], []
    for sa, sb in itertools.combinations(DNK_TRIO, 2):
        cid = f"A_decidua_{sa}_vs_{sb}"
        effects = []
        for col, g in D.groupby("col"):
            a = g[g.subset == sa].set_index("donor")
            b = g[g.subset == sb].set_index("donor")
            donors = sorted(set(a.index) & set(b.index))
            if len(donors) < 2:
                continue
            a, b = a.loc[donors], b.loc[donors]
            diff_pp = float(np.mean((b.detection_rate.to_numpy() -
                                     a.detection_rate.to_numpy()) * 100))
            la = logit(continuity_rate(a.n_detected, a.n_cells))
            lb = logit(continuity_rate(b.n_detected, b.n_cells))
            _, p, _ = exact_wilcoxon_signed_rank(lb - la)
            effects.append({"col": int(col), "mean_diff_pp": diff_pp, "p_raw": p,
                            "n_donors": len(donors)})
        E = pd.DataFrame(effects)
        if E.empty:
            continue
        mu_all = float(E.mean_diff_pp.mean())
        sd_all = float(E.mean_diff_pp.std(ddof=1))
        summary_rows.append({"dataset_id": mf.dataset_id, "comparison_id": cid,
                             "n_null_genes": len(E), "null_mean_diff_pp": mu_all,
                             "null_sd_diff_pp": sd_all,
                             "null_q025": float(E.mean_diff_pp.quantile(0.025)),
                             "null_q975": float(E.mean_diff_pp.quantile(0.975)),
                             "n_donors_median": float(E.n_donors.median())})
        print(f"[null] {cid}: pooled baseline {mu_all:+.3f} +/- {sd_all:.3f} pp over "
              f"{len(E)} genes (NOT assumed to be zero; per-gene sets used below)")

        by_col = dict(zip(E.col, E.mean_diff_pp))
        tt = T2[T2.comparison_id == cid]
        for _, row in tt.iterrows():
            if not np.isfinite(row.mean_diff_pp):
                continue
            # THIS gene's own expression-matched null set, not the pooled one
            own = [by_col[c] for c in gene_to_nulls.get(row.gene, []) if c in by_col]
            if len(own) >= 30:
                obs_eff = np.asarray(own)
                matched = True
            else:
                obs_eff = E.mean_diff_pp.to_numpy()
                matched = False
            mu = float(obs_eff.mean())
            sd = float(obs_eff.std(ddof=1))
            # rank-based empirical p with the standard +1 correction: with N
            # null genes it cannot resolve below 1/(N+1). Reporting anything
            # smaller (e.g. a normal-tail p from z) would claim precision the
            # resampling cannot support.
            n_ge = int((np.abs(obs_eff - mu) >= abs(row.mean_diff_pp - mu)).sum())
            emp_p = (n_ge + 1) / (len(obs_eff) + 1)
            z = float((row.mean_diff_pp - mu) / sd) if sd > 0 else np.nan
            out_rows.append({"dataset_id": mf.dataset_id, "comparison_id": cid,
                             "gene": row.gene, "observed_diff_pp": row.mean_diff_pp,
                             "null_mean_diff_pp": mu, "null_sd_diff_pp": sd,
                             "n_null_genes": len(obs_eff),
                             "null_set_is_expression_matched": matched,
                             "empirical_p": emp_p,
                             "empirical_p_resolution_floor": 1.0 / (len(obs_eff) + 1),
                             "empirical_p_at_resolution_floor": bool(n_ge == 0),
                             "empirical_z": z})

    T4 = pd.DataFrame(out_rows)

    # --- v1.1: BH across the primary-target family, on the EMPIRICAL p -------
    # At these donor counts the signed-rank p is floor-limited and BH over it
    # can never reject. The empirical p against the matched null is not
    # floor-limited in the same way (its resolution is 1/n_null), so this is
    # the multiplicity correction that actually has power here.
    panel_primary = set(panel[panel.role == "primary_target"]["gene_symbol"])
    T4["q_bh_empirical"] = np.nan
    T4["q_bh_empirical_per_comparison"] = np.nan
    if len(T4):
        fam = T4.gene.isin(panel_primary)
        # v1.2: the multiplicity family is EVERY analysis-A test at once
        # (27 genes x 3 pairwise contrasts). The three contrasts are not
        # independent hypothesis families -- they are three views of one set
        # of cells -- so correcting within each separately understates
        # multiplicity by up to 3x. The per-comparison version is retained
        # alongside so the difference is visible rather than assumed away.
        T4.loc[fam, "q_bh_empirical"] = benjamini_hochberg(
            T4.loc[fam, "empirical_p"].to_numpy())
        for cid, g in T4[fam].groupby("comparison_id"):
            T4.loc[g.index, "q_bh_empirical_per_comparison"] = benjamini_hochberg(
                g["empirical_p"].to_numpy())
        m_fam = int(fam.sum())
        print(f"[null] empirical-p BH, ONE family of {m_fam} tests "
              f"({len(panel_primary)} genes x {T4[fam].comparison_id.nunique()} contrasts)")
        n_sig = int((T4["q_bh_empirical"] <= 0.05).sum())
        n_sig_pc = int((T4["q_bh_empirical_per_comparison"] <= 0.05).sum())
        print(f"[null] rows at q <= 0.05: {n_sig} (pooled family) vs "
              f"{n_sig_pc} (per-comparison family)")
        for _, r in T4[T4.q_bh_empirical_per_comparison <= 0.05].sort_values(
                "empirical_p").iterrows():
            mark = "" if r.q_bh_empirical <= 0.05 else "   <- drops out under the pooled family"
            print(f"       {r.comparison_id:26s} {r.gene:7s} "
                  f"diff={r.observed_diff_pp:+7.2f}pp z={r.empirical_z:+6.2f} "
                  f"emp_p={r.empirical_p:.4f} q_pooled={r.q_bh_empirical:.4f}{mark}")
    T4.to_csv(os.path.join(outdir, "null_distribution.tsv"), sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(outdir, "null_distribution_summary.tsv"), sep="\t", index=False)
    print(f"\n[T4] {len(T4)} rows -> {outdir}/null_distribution.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
