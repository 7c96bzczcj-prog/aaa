#!/usr/bin/env python
"""T4 null_distribution.tsv — an empirical null from expression-matched genes.

R6.5: random genes matched to the panel, pushed through the identical
donor-level pipeline. The resulting distribution of effect sizes is the
baseline against which target effects are read.

MATCHING AXIS (v1.3). Null genes are matched on the REFERENCE ARM'S BASELINE
DETECTION RATE, not on CPM. The readout is a detection rate, which is bounded
in [0,1]: a gene detected in 2% of cells cannot move more than 2 points down,
and its sampling variance is p(1-p)/n, set by the rate and not by CPM. Matching
on CPM therefore hands a low-detection target a null set whose variance is
compressed against the zero floor for a different reason than the target's own,
and inflates its z. That is what produced a "significant" 1.5-point S1PR5
effect under CPM matching.

NO FDR IS REPORTED. See docs/PREREGISTRATION_v1.3.md: at this n the number of
rows passing BH is a function of the null-set size, not of the effects, so the
pre-registered FDR criterion is not estimable. Descriptive quantities are
reported instead: standardised effect (z) against the matched null, the effect
in percentage points, donor sign concordance, and ruler-A margin.

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


def match_on_axis(axis, target_cols, n_per_target, rng, exclude, window=0):
    """Pick null genes whose value on `axis` brackets each target gene.

    `axis` is the reference arm's baseline detection rate (v1.3), so the null
    set for a target shares the target's position on the bounded scale the
    effect is measured on.

    The window is auto-sized when not given: a larger null needs a wider rank
    neighbourhood, and silently returning fewer genes would cap the resolution
    without saying so.
    """
    mean_expr = axis
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
    ap.add_argument("--min-null-per-target", type=int, default=30,
                    help="a target with fewer matched nulls than this is not "
                         "reported: its null spread would be unestimable")
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

    target_cols = [hits[g] for g in panel[panel.role == "primary_target"]["gene_symbol"]
                   if g in hits]
    col_to_gene = {hits[g]: g for g in hits}
    n_per = max(1, int(np.ceil(args.n_null / max(1, len(target_cols)))))
    subs_arr = obs_s["subset"].astype(str).to_numpy()
    n_genes = Xs.shape[1]

    def baseline_detection(subset_name):
        """Pooled detection rate of EVERY gene in one subset, at matched depth.

        This is the axis null genes are matched on (v1.3). Computed in column
        blocks so the whole transcriptome never has to be densified at once.
        """
        r = np.flatnonzero(subs_arr == subset_name)
        if r.size == 0:
            return np.zeros(n_genes)
        det = np.zeros(n_genes)
        n_kept = None
        for start in range(0, n_genes, 4000):
            cc = list(range(start, min(start + 4000, n_genes)))
            counts, kept = downsample_marginal(Xs[r], cc, depth, seed=args.seed)
            if n_kept is None:
                n_kept = max(1, int(kept.sum()))
            if counts.shape[0]:
                det[start:start + len(cc)] = (counts >= 1).sum(axis=0) / n_kept
        return det

    print("[null] computing whole-transcriptome baseline detection per subset "
          "(the v1.3 matching axis)", flush=True)
    baseline = {ss: baseline_detection(ss) for ss in DNK_TRIO}
    for ss in DNK_TRIO:
        b = baseline[ss]
        print(f"       {ss}: median {np.median(b):.4f}, "
              f"targets span {b[target_cols].min():.4f}-{b[target_cols].max():.4f}")

    rng = np.random.default_rng(args.seed)

    T2 = pd.read_csv(os.path.join(outdir, "donor_level_tests.tsv"), sep="\t")
    T2 = pd.read_csv(os.path.join(outdir, "donor_level_tests.tsv"), sep="\t")
    out_rows, summary_rows = [], []
    for sa, sb in itertools.combinations(DNK_TRIO, 2):
        cid = f"A_decidua_{sa}_vs_{sb}"
        # match on the REFERENCE arm's baseline detection rate
        axis = baseline[sa]
        per_target = match_on_axis(axis, target_cols, n_per, rng,
                                   exclude=set(hits.values()),
                                   window=args.match_window)
        null_cols = sorted({c for v in per_target.values() for c in v})
        gene_to_nulls = {col_to_gene[c]: v for c, v in per_target.items()
                         if c in col_to_gene}
        sizes = {g: len(v) for g, v in gene_to_nulls.items()}
        print(f"\n[{cid}] reference arm {sa}: {len(null_cols)} null genes, "
              f"per-target min {min(sizes.values())} / median "
              f"{int(np.median(list(sizes.values())))}")
        tq = np.quantile(axis[target_cols], [0.1, 0.5, 0.9])
        nq = np.quantile(axis[null_cols], [0.1, 0.5, 0.9])
        print(f"[{cid}] baseline-detection q10/q50/q90 -- targets {tq.round(3)} "
              f"| nulls {nq.round(3)}")

        # detection rates for this comparison's null genes
        recs = []
        for (donor, ss), r in unit_indices(obs_s, np.ones(len(obs_s), bool),
                                           ["donor", "subset"]):
            if ss not in (sa, sb):
                continue
            counts, kept = downsample_marginal(Xs[r], null_cols, depth, seed=args.seed)
            n_post = int(kept.sum())
            if n_post < args.min_cells:
                continue
            n_det, det, _ = detection_and_cpm(counts, depth)
            for j, c in enumerate(null_cols):
                recs.append({"donor": donor, "subset": ss, "col": c,
                             "n_cells": n_post, "n_detected": int(n_det[j]),
                             "detection_rate": float(det[j])})
        D = pd.DataFrame(recs)
        if D.empty:
            continue

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
            effects.append({"col": int(col), "mean_diff_pp": diff_pp,
                            "n_donors": len(donors)})
        E = pd.DataFrame(effects)
        if E.empty:
            continue
        summary_rows.append({"dataset_id": mf.dataset_id, "comparison_id": cid,
                             "reference_arm": sa, "n_null_genes": len(E),
                             "matching_axis": "baseline_detection_rate",
                             "null_mean_diff_pp": float(E.mean_diff_pp.mean()),
                             "null_sd_diff_pp": float(E.mean_diff_pp.std(ddof=1))})
        print(f"[{cid}] pooled null baseline {E.mean_diff_pp.mean():+.3f} +/- "
              f"{E.mean_diff_pp.std(ddof=1):.3f} pp (per-gene sets used below)")

        by_col = dict(zip(E.col, E.mean_diff_pp))
        tt = T2[T2.comparison_id == cid]
        for _, row in tt.iterrows():
            if not np.isfinite(row.mean_diff_pp):
                continue
            own = [by_col[c] for c in gene_to_nulls.get(row.gene, []) if c in by_col]
            if len(own) < args.min_null_per_target:
                continue
            obs_eff = np.asarray(own)
            mu, sd = float(obs_eff.mean()), float(obs_eff.std(ddof=1))
            n_ge = int((np.abs(obs_eff - mu) >= abs(row.mean_diff_pp - mu)).sum())
            out_rows.append({
                "dataset_id": mf.dataset_id, "comparison_id": cid,
                "reference_arm": sa, "gene": row.gene,
                "observed_diff_pp": row.mean_diff_pp,
                "baseline_detection_ref_arm": float(axis[hits[row.gene]])
                    if row.gene in hits else np.nan,
                "null_mean_diff_pp": mu, "null_sd_diff_pp": sd,
                "n_null_genes": len(obs_eff),
                "null_set_is_matched": True,
                "empirical_z": float((row.mean_diff_pp - mu) / sd) if sd > 0 else np.nan,
                "n_null_at_least_as_extreme": n_ge,
                "rank_of_observed": n_ge + 1})

    T4 = pd.DataFrame(out_rows, columns=[
        "dataset_id", "comparison_id", "reference_arm", "gene",
        "observed_diff_pp", "baseline_detection_ref_arm", "null_mean_diff_pp",
        "null_sd_diff_pp", "n_null_genes", "null_set_is_matched",
        "empirical_z", "n_null_at_least_as_extreme", "rank_of_observed"])
    if T4.empty:
        print("[null] no target had enough matched null genes to report "
              f"(threshold {args.min_null_per_target}); writing an empty table")
        T4["fdr_estimable"] = pd.Series(dtype=bool)
        T4["fdr_not_estimable_reason"] = pd.Series(dtype=str)

    # --- v1.3: NO FDR. The pre-registered criterion is not estimable here. ---
    #
    # Walk the arithmetic. BH at m = 81 needs the smallest p to reach
    # 0.05 * k / 81 for k rows tied at the resolution floor 1/(N+1):
    #     N = 112 -> floor 0.00885 -> needs k >= 15
    #     N = 317 -> floor 0.00314 -> needs k >=  6
    #     N = 399 -> floor 0.00251 -> needs k >=  5
    # There are only 8 candidate rows. So the number of rows that "pass" is a
    # function of the null-set size, not of the effects. That is not a rule
    # needing a tuned parameter; it is a rule that does not hold at this n.
    #
    # Reported instead: standardised effect against the matched null, the
    # effect in percentage points, donor sign concordance, and ruler-A margin.
    # No replacement decision rule is invented -- see docs/PREREGISTRATION_v1.3.md.
    if len(T4):
        T4["fdr_estimable"] = False
        T4["fdr_not_estimable_reason"] = (
            "pass count is a function of null-set size, not of effect size, at "
            "this n; see PREREGISTRATION_v1.3.md")
        n_min = int(T4.n_null_genes.min())
        print(f"\n[null] NO FDR REPORTED. Null-set size {n_min}-"
              f"{int(T4.n_null_genes.max())} per target; the rank floor "
              f"1/{n_min + 1} = {1 / (n_min + 1):.4f} would require "
              f"{int(np.ceil(81 * (1 / (n_min + 1)) / 0.05))} of 81 rows tied at "
              f"it before any could clear q <= 0.05, and only 8 rows are "
              f"candidates. The criterion is not estimable, not merely unmet.")

        T2f = T2.set_index(["comparison_id", "gene"])
        rows = []
        for _, r in T4.iterrows():
            key = (r.comparison_id, r.gene)
            src = T2f.loc[key] if key in T2f.index else None
            rows.append({
                "comparison": r.comparison_id.replace("A_decidua_", ""),
                "gene": r.gene,
                "effect_pp": r.observed_diff_pp,
                "z_vs_matched_null": r.empirical_z,
                "n_donors": int(src.n_donors) if src is not None else np.nan,
                "all_donors_same_sign": bool(src.on_test_floor) if src is not None else False,
                "baseline_det_ref": r.baseline_detection_ref_arm,
                "n_nulls": int(r.n_null_genes),
                "null_rank": int(r.rank_of_observed)})
        R = pd.DataFrame(rows)
        R = R.reindex(R.z_vs_matched_null.abs().sort_values(ascending=False).index)
        print("\n[null] descriptive readout, primary targets, |z| >= 3 "
              "(z is a standardised effect size and is NOT converted to a p):")
        show = R[(R.z_vs_matched_null.abs() >= 3)]
        print(show.round(3).to_string(index=False))
        R.to_csv(os.path.join(outdir, "null_descriptive.tsv"), sep="\t", index=False)

    T4.to_csv(os.path.join(outdir, "null_distribution.tsv"), sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(outdir, "null_distribution_summary.tsv"), sep="\t", index=False)
    print(f"\n[T4] {len(T4)} rows -> {outdir}/null_distribution.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
