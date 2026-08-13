#!/usr/bin/env python
"""Is the dNK chemokine gradient real, or a restatement of the cluster definition?

THE PROBLEM. dNK1/dNK2/dNK3 come from Vento-Tormo's own unsupervised
whole-transcriptome clustering, and the subsets were then characterised by
differential expression. CCL5, CXCR4 and XCL1 are among the genes that
characterisation named. Re-discovering, in the same data and under the same
cluster labels, that dNK3 is CCL5/CXCR4-high therefore partly restates how the
clusters were drawn. Nobody thresholded on CCL5 -- but CCL5 helped set the
cluster boundary, because the boundary came from the whole transcriptome.

THE TEST. Throw the published labels away. Re-group the same decidual NK cells
using ONLY surface-marker genes that map onto the Huhn CyTOF gating -- markers
that are independent of the chemokines being measured -- and re-run the same
donor-level pipeline.

  gate_1 (~dNK1) : ENTPD1+                 (CD39+)
  gate_3 (~dNK3) : ENTPD1- ITGAE+          (CD39- CD103+)
  gate_2 (~dNK2) : ENTPD1- ITGAE-          (CD39- CD103-)

ITGB2 (CD18) is reported alongside but not used to gate, so that a third
marker is available as an independent check on the gate assignment.

CD160 and KLRB1 are deliberately NOT used: they are themselves dNK3
characterisation genes, so gating on them would rebuild the circularity the
test exists to break.

READING IT. Gradient survives -> the effect is not an artefact of the
clustering. Gradient collapses -> it was circular.

This does NOT re-label anything: the published annotation is untouched, the
gates are a parallel grouping written to their own table. No output of this
script feeds the main analysis.

    python src/circularity_test.py manifests/VT2018.yaml
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
                            downsample_columns)
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402
from dnkchem.stats import (bootstrap_ci, continuity_rate,  # noqa: E402
                           exact_wilcoxon_signed_rank, logit, min_achievable_p)

GATE_GENES = ["ENTPD1", "ITGAE"]          # CD39, CD103
REPORT_GENE = "ITGB2"                     # CD18, reported not gated
BARRED = ["CD160", "KLRB1"]               # dNK3 characterisation genes
TARGETS = ["XCL1", "XCL2", "CCL5", "CXCR4", "CCL3", "CCL4", "CCL4L2", "CCL3L1"]
TRIO = ["dNK1", "dNK2", "dNK3"]


def assign_gates(counts, gi):
    """Marker-gated grouping, computed on depth-matched counts."""
    cd39 = counts[:, gi["ENTPD1"]] >= 1
    cd103 = counts[:, gi["ITGAE"]] >= 1
    gate = np.full(counts.shape[0], "", dtype=object)
    gate[cd39] = "gate1_CD39pos"
    gate[~cd39 & cd103] = "gate3_CD39neg_CD103pos"
    gate[~cd39 & ~cd103] = "gate2_CD39neg_CD103neg"
    return gate


def paired_test(df, ga, gb, gene, min_cells, seed):
    a = df[(df.group == ga) & (df.gene == gene) & (df.n_cells >= min_cells)]
    b = df[(df.group == gb) & (df.gene == gene) & (df.n_cells >= min_cells)]
    a = a.set_index("donor"); b = b.set_index("donor")
    donors = sorted(set(a.index) & set(b.index))
    if not donors:
        return None
    a, b = a.loc[donors], b.loc[donors]
    diffs = (b.detection_rate.to_numpy() - a.detection_rate.to_numpy()) * 100
    la = logit(continuity_rate(a.n_detected, a.n_cells))
    lb = logit(continuity_rate(b.n_detected, b.n_cells))
    _, p, n_used = exact_wilcoxon_signed_rank(lb - la)
    lo, hi = bootstrap_ci(diffs, seed=seed)
    return {"gene": gene, "n_donors": len(donors), "mean_diff_pp": float(diffs.mean()),
            "ci_low": lo, "ci_high": hi, "p_raw": float(p),
            "min_achievable_p": float(min_achievable_p(n_used)),
            "all_donors_same_sign": bool(np.all(diffs > 0) or np.all(diffs < 0))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    depth = json.load(open(os.path.join(outdir, "detection_meta.json")))["primary_depth_floor"]

    ds = load_dataset(mf)
    hits, _, _, _ = ds.panel_columns(panel)
    X = as_csr(ds.X)
    keep, _ = cell_qc(X, ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs

    need = GATE_GENES + [REPORT_GENE] + TARGETS
    need = [g for g in need if g in hits]
    cols = [hits[g] for g in need]
    gi = {g: k for k, g in enumerate(need)}
    print(f"[gate] gating on {GATE_GENES} (reported: {REPORT_GENE}); "
          f"barred as circular: {BARRED}")
    print(f"[gate] targets never used for gating: {[t for t in TARGETS if t in gi]}")

    sel = (keep & (obs["compartment"] == "Decidua").to_numpy()
           & obs["subset"].isin(TRIO).to_numpy())
    print(f"[gate] decidual dNK1-3 cells after QC: {int(sel.sum())}")

    # depth-match FIRST, so the gate itself is not a function of depth
    rows, xtab = [], []
    for (donor,), idx in unit_indices(obs, sel, ["donor"]):
        counts, kept = downsample_columns(X[idx], cols, depth, seed=args.seed)
        if kept.sum() == 0:
            continue
        orig = obs["subset"].to_numpy()[idx][kept]
        gate = assign_gates(counts, gi)
        for g in np.unique(gate):
            m = gate == g
            n_det, det, _ = detection_and_cpm(counts[m], depth)
            for t in TARGETS:
                if t not in gi:
                    continue
                rows.append({"donor": donor, "group": g, "gene": t,
                             "n_cells": int(m.sum()), "n_detected": int(n_det[gi[t]]),
                             "detection_rate": float(det[gi[t]])})
        for o, g in zip(orig, gate):
            xtab.append({"donor": donor, "original": o, "gate": g})

    D = pd.DataFrame(rows)
    XT = pd.DataFrame(xtab)
    D.to_csv(os.path.join(outdir, "circularity_gate_rates.tsv"), sep="\t", index=False)

    print("\n=== gate vs published label (all donors pooled) ===")
    ct = pd.crosstab(XT.original, XT.gate)
    print(ct.to_string())
    print("\nrow-normalised (what fraction of each published subset lands in each gate):")
    print((ct.div(ct.sum(axis=1), axis=0) * 100).round(1).to_string())
    ct.to_csv(os.path.join(outdir, "circularity_crosstab.tsv"), sep="\t")

    # agreement: is the gate recovering the published structure at all?
    best = {"dNK1": "gate1_CD39pos", "dNK2": "gate2_CD39neg_CD103neg",
            "dNK3": "gate3_CD39neg_CD103pos"}
    agree = sum(ct.loc[o, g] for o, g in best.items() if o in ct.index and g in ct.columns)
    print(f"\nconcordance with the published labels: {agree}/{ct.values.sum()} "
          f"= {agree / ct.values.sum():.1%}")

    print(f"\n=== cells per donor x gate (admission threshold {args.min_cells}) ===")
    per = D[D.gene == TARGETS[0]].pivot_table(index="group", columns="donor",
                                              values="n_cells", aggfunc="sum")
    print(per.fillna(0).astype(int).to_string())

    print("\n=== the gradient, re-measured on marker gates ===")
    print("(gate assignment never looks at any gene below)")
    out = []
    pairs = [("gate1_CD39pos", "gate2_CD39neg_CD103neg", "gate1->gate2 (~dNK1->dNK2)"),
             ("gate1_CD39pos", "gate3_CD39neg_CD103pos", "gate1->gate3 (~dNK1->dNK3)"),
             ("gate2_CD39neg_CD103neg", "gate3_CD39neg_CD103pos", "gate2->gate3 (~dNK2->dNK3)")]
    for ga, gb, label in pairs:
        for t in TARGETS:
            if t not in gi:
                continue
            r = paired_test(D, ga, gb, t, args.min_cells, args.seed)
            if r:
                r.update({"comparison": label, "group_a": ga, "group_b": gb})
                out.append(r)
    R = pd.DataFrame(out)
    R["dataset_id"] = mf.dataset_id
    R.to_csv(os.path.join(outdir, "circularity_tests.tsv"), sep="\t", index=False)

    T2 = pd.read_csv(os.path.join(outdir, "donor_level_tests.tsv"), sep="\t")
    orig_map = {"gate1->gate2 (~dNK1->dNK2)": "A_decidua_dNK1_vs_dNK2",
                "gate1->gate3 (~dNK1->dNK3)": "A_decidua_dNK1_vs_dNK3",
                "gate2->gate3 (~dNK2->dNK3)": "A_decidua_dNK2_vs_dNK3"}
    R["original_comparison"] = R.comparison.map(orig_map)
    m = R.merge(T2[["comparison_id", "gene", "mean_diff_pp", "n_donors", "p_raw"]],
                left_on=["original_comparison", "gene"],
                right_on=["comparison_id", "gene"], suffixes=("_gate", "_label"))
    m["retained_frac"] = m.mean_diff_pp_gate / m.mean_diff_pp_label
    for label in [p[2] for p in pairs]:
        s = m[m.comparison == label]
        if not len(s):
            continue
        print(f"\n--- {label} ---")
        print(s[["gene", "n_donors_gate", "mean_diff_pp_gate", "ci_low", "ci_high",
                 "p_raw_gate", "all_donors_same_sign", "mean_diff_pp_label",
                 "retained_frac"]].round(3).to_string(index=False))
    m.to_csv(os.path.join(outdir, "circularity_comparison.tsv"), sep="\t", index=False)

    # ---------------------------------------------------------------------
    # WITHIN-LABEL test: the one that cannot be circular.
    #
    # A gate-vs-gate difference computed INSIDE a single published cluster
    # cannot be a restatement of how that cluster was drawn -- the cluster is
    # held fixed, and the gate splits it on markers the chemokines never
    # informed. It also removes the dilution confound: a shrunken effect in
    # the pooled test can be explained by impure gates (gate3 is only 35%
    # dNK3), but a surviving effect within one label cannot.
    # ---------------------------------------------------------------------
    print("\n=== WITHIN-LABEL test: gate contrast inside one published cluster ===")
    wrows = []
    for (donor, orig), idx in unit_indices(obs, sel, ["donor", "subset"]):
        counts, kept = downsample_columns(X[idx], cols, depth, seed=args.seed)
        if kept.sum() == 0:
            continue
        gate = assign_gates(counts, gi)
        for g in np.unique(gate):
            gm = gate == g
            if gm.sum() < args.min_cells:
                continue
            n_det, det, _ = detection_and_cpm(counts[gm], depth)
            for t in TARGETS:
                if t not in gi:
                    continue
                wrows.append({"donor": donor, "original": orig, "group": g,
                              "n_cells": int(gm.sum()),
                              "n_detected": int(n_det[gi[t]]),
                              "detection_rate": float(det[gi[t]]), "gene": t})
    W = pd.DataFrame(wrows)
    W.to_csv(os.path.join(outdir, "circularity_within_label_rates.tsv"),
             sep="\t", index=False)
    if len(W):
        print("cells per donor x label x gate that clear the threshold:")
        print(W[W.gene == TARGETS[0]].pivot_table(
            index=["original", "group"], columns="donor",
            values="n_cells", aggfunc="sum").fillna(0).astype(int).to_string())

    wout = []
    for orig in TRIO:
        sub = W[W.original == orig]
        if not len(sub):
            continue
        for ga, gb in itertools.combinations(sorted(sub.group.unique()), 2):
            for t in TARGETS:
                if t not in gi:
                    continue
                r = paired_test(sub, ga, gb, t, args.min_cells, args.seed)
                if r and r["n_donors"] >= 3:
                    r.update({"original_label": orig, "group_a": ga, "group_b": gb})
                    wout.append(r)
    WT = pd.DataFrame(wout)
    if len(WT):
        WT["dataset_id"] = mf.dataset_id
        WT.to_csv(os.path.join(outdir, "circularity_within_label_tests.tsv"),
                  sep="\t", index=False)
        key = WT[WT.gene.isin(["XCL1", "XCL2", "CCL5", "CXCR4"])]
        print("\ncontrasts with n >= 3 donors, inside a single published cluster:")
        print(key[["original_label", "group_a", "group_b", "gene", "n_donors",
                   "mean_diff_pp", "ci_low", "ci_high", "p_raw",
                   "min_achievable_p", "all_donors_same_sign"]].round(3).to_string(index=False))
    else:
        print("no within-label contrast reaches 3 donors at the admission threshold")

    print("\n=== verdict ===")
    for t in ["XCL1", "XCL2", "CCL5", "CXCR4"]:
        s = m[m.gene == t]
        if not len(s):
            continue
        surv = s[(np.sign(s.mean_diff_pp_gate) == np.sign(s.mean_diff_pp_label))
                 & (s.mean_diff_pp_gate.abs() > 5)]
        print(f"  {t:6s}: sign preserved and |effect| > 5pp in "
              f"{len(surv)}/{len(s)} contrasts; median retained "
              f"{s.retained_frac.median():.2f}x of the label-based effect")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
