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

from dnkchem.counts import as_csr, cell_qc, detection_and_cpm, downsample_columns  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402
from dnkchem.stats import continuity_rate, exact_wilcoxon_signed_rank, logit  # noqa: E402

DNK_TRIO = ["dNK1", "dNK2", "dNK3"]


def match_expression(mean_expr, target_cols, n_per_target, rng, exclude):
    """Pick genes whose pooled expression brackets each target gene."""
    order = np.argsort(mean_expr)
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    banned = set(exclude)
    chosen = []
    for c in target_cols:
        r = rank[c]
        lo, hi = max(0, r - 250), min(len(order), r + 250)
        pool = [int(order[k]) for k in range(lo, hi)
                if int(order[k]) not in banned and mean_expr[int(order[k])] > 0]
        rng.shuffle(pool)
        take = pool[:n_per_target]
        chosen.extend(take)
        banned.update(take)
    return sorted(set(chosen))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--n-null", type=int, default=400)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=20260813)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = os.path.join("out", mf.dataset_id)
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
    null_cols = match_expression(mean_expr, target_cols, n_per, rng,
                                 exclude=set(hits.values()))
    print(f"[null] {len(null_cols)} expression-matched null genes "
          f"({n_per} per target, target CPM range "
          f"{mean_expr[target_cols].min():.2f}-{mean_expr[target_cols].max():.2f})")
    if len(null_cols) < args.n_null:
        print(f"[null] WARNING: only {len(null_cols)} null genes available "
              f"(requested {args.n_null})")

    # per donor x subset detection rates for the null genes, same depth, same seed
    recs = []
    for (donor, ss), idx in obs_s.groupby(["donor", "subset"], observed=True).groups.items():
        r = obs_s.index.get_indexer(idx)
        counts, kept = downsample_columns(Xs[r], null_cols, depth, seed=args.seed)
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
        mu, sd = float(E.mean_diff_pp.mean()), float(E.mean_diff_pp.std(ddof=1))
        summary_rows.append({"dataset_id": mf.dataset_id, "comparison_id": cid,
                             "n_null_genes": len(E), "null_mean_diff_pp": mu,
                             "null_sd_diff_pp": sd,
                             "null_q025": float(E.mean_diff_pp.quantile(0.025)),
                             "null_q975": float(E.mean_diff_pp.quantile(0.975)),
                             "n_donors_median": float(E.n_donors.median())})
        print(f"[null] {cid}: baseline {mu:+.3f} +/- {sd:.3f} pp over {len(E)} genes "
              f"(NOT assumed to be zero)")

        obs_eff = E.mean_diff_pp.to_numpy()
        tt = T2[T2.comparison_id == cid]
        for _, row in tt.iterrows():
            if not np.isfinite(row.mean_diff_pp):
                continue
            emp_p = float((np.abs(obs_eff - mu) >= abs(row.mean_diff_pp - mu)).mean())
            z = float((row.mean_diff_pp - mu) / sd) if sd > 0 else np.nan
            out_rows.append({"dataset_id": mf.dataset_id, "comparison_id": cid,
                             "gene": row.gene, "observed_diff_pp": row.mean_diff_pp,
                             "null_mean_diff_pp": mu, "null_sd_diff_pp": sd,
                             "n_null_genes": len(E), "empirical_p": emp_p,
                             "empirical_z": z})

    T4 = pd.DataFrame(out_rows)
    T4.to_csv(os.path.join(outdir, "null_distribution.tsv"), sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(outdir, "null_distribution_summary.tsv"), sep="\t", index=False)
    print(f"\n[T4] {len(T4)} rows -> {outdir}/null_distribution.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
