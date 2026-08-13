#!/usr/bin/env python
"""T2 donor_level_tests.tsv — every test consumes one number per donor (R1).

Analysis A (primary): dNK1 vs dNK2 vs dNK3 inside decidua, donor-paired.
  Same library, same dissociation, same soup pool, so ambient RNA is very
  nearly a common additive offset and the between-subset contrast is close to
  immune to it.
  dNKp is NOT tested against the other three: a proliferating population has
  different RNA content and transcriptome structure. Its rates appear in T1.

Analysis B (secondary): decidual dNK vs the SAME donor's blood NK.
  Every row carries cross_compartment_soup_caveat = TRUE, unconditionally.
  CD56-bright and CD56-dim blood NK are never merged.

Every row reports min_achievable_p and on_test_floor (R2).

    python src/donor_tests.py manifests/VT2018.yaml
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

from dnkchem.dataset import load_panel  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402
from dnkchem.stats import (benjamini_hochberg, bootstrap_ci,  # noqa: E402
                           continuity_rate, exact_wilcoxon_signed_rank, logit,
                           min_achievable_p)

DNK_TRIO = ["dNK1", "dNK2", "dNK3"]
BLOOD_NK = ["pbNK_CD56bright", "pbNK_CD56dim"]
CEILING = 0.85   # must match detection_rates.CEILING
FLOOR = 0.05     # must match detection_rates.FLOOR


def paired_frame(T1, comp_a, sub_a, comp_b, sub_b, gene, min_cells):
    a = T1[(T1.compartment == comp_a) & (T1.subset == sub_a) & (T1.gene == gene)]
    b = T1[(T1.compartment == comp_b) & (T1.subset == sub_b) & (T1.gene == gene)]
    a = a[a.n_cells >= min_cells].set_index("donor_id")
    b = b[b.n_cells >= min_cells].set_index("donor_id")
    donors = sorted(set(a.index) & set(b.index))
    return a.loc[donors], b.loc[donors], donors


def run_comparison(T1, comparison_id, comp_a, sub_a, comp_b, sub_b, genes,
                   min_cells, cross_compartment, seed):
    out = []
    for gene in genes:
        a, b, donors = paired_frame(T1, comp_a, sub_a, comp_b, sub_b, gene, min_cells)
        n = len(donors)
        if n == 0:
            out.append(dict(gene=gene, n_donors=0, mean_diff_pp=np.nan,
                            ci_low=np.nan, ci_high=np.nan, p_raw=np.nan,
                            min_achievable_p=np.nan, on_test_floor=False,
                            uninterpretable_ceiling=False, uninterpretable_floor=False,
                            ceiling_arm_a=False, ceiling_arm_b=False,
                            magnitude_is_lower_bound=False))
            continue
        ra = a["detection_rate"].to_numpy()
        rb = b["detection_rate"].to_numpy()
        diffs_pp = (rb - ra) * 100.0
        # test on the logit scale; report the effect in percentage points
        la = logit(continuity_rate(a["n_detected"], a["n_cells"]))
        lb = logit(continuity_rate(b["n_detected"], b["n_cells"]))
        _, p, n_used = exact_wilcoxon_signed_rank(lb - la)
        minp = min_achievable_p(n_used)
        lo, hi = bootstrap_ci(diffs_pp, seed=seed)
        eff = float(np.mean(diffs_pp))

        # --- R5, v1.1: the ceiling rule is per ARM ---------------------------
        # v1.0 flagged a contrast uninterpretable when EITHER arm exceeded 0.85.
        # That is right when both arms saturate (a null is then unreadable, the
        # lesson from the 93.75% case). It is wrong when only the HIGH arm
        # saturates: saturation compresses the difference toward zero, so a
        # large observed difference is a LOWER BOUND, not an unreadable number.
        ceil_a = float(np.mean(ra)) > CEILING
        ceil_b = float(np.mean(rb)) > CEILING
        floor_a = float(np.mean(ra)) < FLOOR
        floor_b = float(np.mean(rb)) < FLOOR
        # one arm at the ceiling, and the effect points TOWARD that arm
        one_sided_ceiling = ((ceil_b and not ceil_a and eff > 0) or
                             (ceil_a and not ceil_b and eff < 0))
        one_sided_floor = ((floor_b and not floor_a and eff < 0) or
                           (floor_a and not floor_b and eff > 0))

        out.append(dict(
            gene=gene, n_donors=n,
            mean_diff_pp=eff, ci_low=lo, ci_high=hi,
            p_raw=float(p), min_achievable_p=float(minp),
            on_test_floor=bool(np.isfinite(p) and np.isfinite(minp)
                               and abs(p - minp) < 1e-12),
            uninterpretable_ceiling=bool(ceil_a and ceil_b),
            uninterpretable_floor=bool(floor_a and floor_b),
            ceiling_arm_a=bool(ceil_a), ceiling_arm_b=bool(ceil_b),
            magnitude_is_lower_bound=bool(one_sided_ceiling or one_sided_floor)))
    df = pd.DataFrame(out)
    df.insert(0, "comparison_id", comparison_id)
    df.insert(1, "group_a", f"{comp_a}:{sub_a}")
    df.insert(2, "group_b", f"{comp_b}:{sub_b}")
    df["cross_compartment_soup_caveat"] = bool(cross_compartment)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    meta = json.load(open(os.path.join(outdir, "detection_meta.json")))
    depth = meta["primary_depth_floor"]

    T1_all = pd.read_csv(os.path.join(outdir, "detection_rates.tsv"), sep="\t")
    T1 = T1_all[T1_all.depth_floor == depth]

    # BH is applied across the primary-target set only. XCR1 is a control, not
    # a hypothesis; markers, ambient controls and prospective genes are not
    # tested. Pre-registered in docs/PREREGISTRATION.md section 2.
    primary = panel[(panel.role == "primary_target")]["gene_symbol"].tolist()
    reported = panel[panel.category.isin(
        ["ligand", "receptor", "prospective"])]["gene_symbol"].tolist()
    genes = [g for g in reported if g in set(T1.gene)]
    print(f"[genes] {len(genes)} reported, {len(primary)} in the BH family")

    decidua = "Decidua"
    blocks = []

    # --- Analysis A -------------------------------------------------------
    for sa, sb in itertools.combinations(DNK_TRIO, 2):
        df = run_comparison(T1, f"A_decidua_{sa}_vs_{sb}", decidua, sa, decidua, sb,
                            genes, args.min_cells, False, args.seed)
        blocks.append(df)
        n = df["n_donors"].max()
        print(f"[A] {sa} vs {sb}: n_donors up to {n}, "
              f"min achievable p = {min_achievable_p(n) if n else float('nan'):.4f}")

    # --- Analysis B -------------------------------------------------------
    for sa in DNK_TRIO:
        for sb in BLOOD_NK:
            df = run_comparison(T1, f"B_{sa}_vs_blood_{sb}", decidua, sa,
                                "Blood", sb, genes, args.min_cells, True, args.seed)
            blocks.append(df)
            n = df["n_donors"].max()
            print(f"[B] {sa} vs {sb}: n_donors up to {n}, "
                  f"min achievable p = {min_achievable_p(n) if n else float('nan'):.4f}")

    T2 = pd.concat(blocks, ignore_index=True)

    # BH within each comparison, across the primary-target family only
    T2["q_bh"] = np.nan
    fam = T2.gene.isin(primary)
    for cid, g in T2[fam].groupby("comparison_id"):
        T2.loc[g.index, "q_bh"] = benjamini_hochberg(g["p_raw"].to_numpy())

    T2["dataset_id"] = mf.dataset_id
    T2 = T2[["dataset_id", "comparison_id", "group_a", "group_b", "gene", "n_donors",
             "mean_diff_pp", "ci_low", "ci_high", "p_raw", "q_bh",
             "min_achievable_p", "on_test_floor", "uninterpretable_ceiling",
             "uninterpretable_floor", "cross_compartment_soup_caveat",
             # v1.1 additions, appended after the frozen v1.0 schema
             "ceiling_arm_a", "ceiling_arm_b", "magnitude_is_lower_bound"]]
    T2.to_csv(os.path.join(outdir, "donor_level_tests.tsv"), sep="\t", index=False)
    print(f"\n[T2] {len(T2)} rows -> {outdir}/donor_level_tests.tsv")

    # --- readability gate for analysis B (pre-registered) -----------------
    xcr1 = T1[(T1.gene == "XCR1") & (T1.n_cells >= args.min_cells) &
              (T1.subset.str.startswith(("dNK", "pbNK")))]
    gate = {}
    if len(xcr1):
        by = xcr1.groupby("compartment")["detection_rate"].mean()
        if {"Decidua", "Blood"} <= set(by.index):
            delta_pp = float(abs(by["Decidua"] - by["Blood"]) * 100)
            gate = {"xcr1_decidua_pp": float(by["Decidua"] * 100),
                    "xcr1_blood_pp": float(by["Blood"] * 100),
                    "delta_pp": delta_pp, "threshold_pp": 2.0,
                    "analysis_B_readable": bool(delta_pp < 2.0)}
            print(f"[B gate] XCR1 in NK gate: decidua {by['Decidua']*100:.2f}pp vs "
                  f"blood {by['Blood']*100:.2f}pp, delta {delta_pp:.2f}pp "
                  f"-> B is {'READABLE' if delta_pp < 2.0 else 'NOT readable; report only'}")
    with open(os.path.join(outdir, "analysis_b_gate.json"), "w") as fh:
        json.dump(gate, fh, indent=2)

    # --- R6 sensitivity: repeat every comparison at all depth floors -------
    # Mandatory whenever downsampling retention differs by >10 percentage
    # points across compared groups, which it does here.
    sens = []
    for d in sorted(T1_all.depth_floor.unique()):
        Td = T1_all[T1_all.depth_floor == d]
        blocks_d = []
        for sa, sb in itertools.combinations(DNK_TRIO, 2):
            blocks_d.append(run_comparison(Td, f"A_decidua_{sa}_vs_{sb}", decidua, sa,
                                           decidua, sb, genes, args.min_cells,
                                           False, args.seed))
        for sa in DNK_TRIO:
            for sb in BLOOD_NK:
                blocks_d.append(run_comparison(Td, f"B_{sa}_vs_blood_{sb}", decidua, sa,
                                               "Blood", sb, genes, args.min_cells,
                                               True, args.seed))
        b = pd.concat(blocks_d, ignore_index=True)
        b["depth_floor"] = d
        b["is_primary_floor"] = (d == depth)
        sens.append(b)
    S = pd.concat(sens, ignore_index=True)
    S["dataset_id"] = mf.dataset_id
    S = S[["dataset_id", "depth_floor", "is_primary_floor", "comparison_id",
           "group_a", "group_b", "gene", "n_donors", "mean_diff_pp", "ci_low",
           "ci_high", "p_raw", "min_achievable_p", "on_test_floor",
           "uninterpretable_ceiling", "uninterpretable_floor",
           "cross_compartment_soup_caveat", "magnitude_is_lower_bound"]]
    S.to_csv(os.path.join(outdir, "donor_level_tests_depth_sensitivity.tsv"),
             sep="\t", index=False)
    print(f"[R6] depth sensitivity: {len(S)} rows across floors "
          f"{sorted(T1_all.depth_floor.unique())} -> "
          f"{outdir}/donor_level_tests_depth_sensitivity.tsv")
    # does any conclusion flip sign across floors?
    flip = (S[S.n_donors > 0].groupby(["comparison_id", "gene"])["mean_diff_pp"]
            .agg(lambda v: (v > 0).any() and (v < 0).any()))
    nflip = int(flip.sum())
    print(f"[R6] gene x comparison cells whose effect changes SIGN across floors: "
          f"{nflip} / {len(flip)}")

    sig = T2[(T2.q_bh <= 0.05) & T2.gene.isin(primary)]
    print(f"\n[summary] rows with q_bh <= 0.05 in the primary family: {len(sig)} "
          f"(v1.1: BH is reported, not adjudicating -- see PREREGISTRATION_v1.1.md)")
    print(f"[summary] contrasts with BOTH arms saturated (uninterpretable): "
          f"{int(T2.uninterpretable_ceiling.sum())}")
    print(f"[summary] contrasts with ONE arm saturated (direction kept, "
          f"magnitude is a lower bound): {int(T2.magnitude_is_lower_bound.sum())}")
    floored = T2[T2.on_test_floor & T2.gene.isin(primary)]
    print(f"[summary] rows sitting exactly on the test floor: {len(floored)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
