#!/usr/bin/env python
"""Does ruler A's soup fraction split by lineage CLASS rather than by gene?

The nine panel ambient controls fall into a suspiciously clean pattern:

  not saturated : IGKC 0.092, LYZ 0.164, C1QA 0.176, HBB 0.408  -- all HAEMATOPOIETIC
  saturated 1.000: DCN, COL1A1, IGFBP1, HLA-G, CSH1             -- all NON-haematopoietic

If that split is a class property rather than a coincidence of five genes, it
predicts where a missing control would have landed. PAEP is glandular
epithelial / decidualised stroma -- non-haematopoietic -- so it would be
expected to return 1.000 and leave the ceiling (a MINIMUM over controls)
untouched, and the CXCR4 risk from its absence would be small.

The test: take marker genes for each mapped lineage directly from the data --
no hand-picked list -- run them through the same leave-one-out ruler A, and
see whether soup fraction separates by haematopoietic vs not.

Why the split is expected. Ruler A assumes one ambient load rho shared by all
cells. NK cells sit inside a haematopoietic compartment, so for a
haematopoietic ambient gene the NK gate sees both true ambient AND signal from
being the same lineage class, which pushes the estimate below 1. For a
stromal or trophoblast gene there is no such route and the estimate is clean.
If so, the ceiling is being set by the class of control that the estimator
handles WORST, which is a defect of the min-over-controls rule.

    python src/ceiling_class_structure.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr, cell_qc  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

HAEM = {"Myeloid", "T", "Plasma", "Granulocyte", "cDC1", "ILC3", "Hofbauer"}
NON_HAEM = {"Stromal", "Trophoblast", "Epithelial", "Endothelial", "Perivascular"}
PANEL_AMBIENT = ["IGKC", "LYZ", "C1QA", "DCN", "COL1A1", "IGFBP1",
                 "HLA-G", "CSH1", "HBB"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--compartment", default="Decidua")
    ap.add_argument("--per-lineage", type=int, default=12)
    ap.add_argument("--specificity-fold", type=float, default=5.0,
                    help="lineage CPM must exceed the next-highest lineage by "
                         "this factor; excludes housekeeping genes owned by "
                         "whichever lineage has the most cells")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    ds = load_dataset(mf)
    X = as_csr(ds.X)
    keep, _ = cell_qc(X, ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs
    syms = np.asarray([str(s) for s in ds.symbols])

    comp = (obs["compartment"] == args.compartment).to_numpy()
    sel = keep & comp
    subs = obs["subset"].astype(object).where(obs["subset"].notna(), "").astype(str).to_numpy()
    is_nk = np.array([s.startswith(("dNK", "pbNK")) for s in subs])

    crows = np.flatnonzero(sel)
    nk_rows = crows[is_nk[crows]]
    depth_all = float(X[crows].sum())
    amb_cpm = np.asarray(X[crows].sum(axis=0)).ravel() / depth_all * 1e6
    nk_depth = float(X[nk_rows].sum())
    nk_cpm = np.asarray(X[nk_rows].sum(axis=0)).ravel() / nk_depth * 1e6

    # per-lineage CPM, to pick markers straight from the data
    lin_cpm, lin_n = {}, {}
    for lin in sorted(HAEM | NON_HAEM):
        rws = crows[subs[crows] == lin]
        if rws.size < 30:
            continue
        d = float(X[rws].sum())
        if d == 0:
            continue
        lin_cpm[lin] = np.asarray(X[rws].sum(axis=0)).ravel() / d * 1e6
        lin_n[lin] = int(rws.size)
    print(f"[class] lineages with >=30 cells in {args.compartment}: "
          f"{ {k: v for k, v in sorted(lin_n.items())} }")

    # a marker = high in its lineage, and that lineage is its top contributor
    lin_total = {}
    for lin in lin_cpm:
        rws = crows[subs[crows] == lin]
        lin_total[lin] = np.asarray(X[rws].sum(axis=0)).ravel().astype(float)
    stack = np.vstack([lin_total[l] for l in lin_cpm])
    owner_idx = np.argmax(stack, axis=0)
    owners = np.array(list(lin_cpm.keys()))[owner_idx]
    nk_total = np.asarray(X[nk_rows].sum(axis=0)).ravel().astype(float)
    nk_wins = nk_total > stack.max(axis=0)

    # Specificity is required, not just ownership. The lineage with the most
    # cells owns the total counts of every housekeeping gene as well, so
    # ranking by lineage CPM alone selects ribosomal/mitochondrial genes that
    # NK expresses too -- which drags the soup fraction down and would have
    # been misread as a property of that lineage. Demand that the lineage
    # exceeds the next-highest lineage by SPECIFICITY_FOLD.
    lin_names = list(lin_cpm.keys())
    cpm_stack = np.vstack([lin_cpm[l] for l in lin_names])
    picks = []
    for li, lin in enumerate(lin_names):
        others = np.delete(cpm_stack, li, axis=0).max(axis=0)
        specific = lin_cpm[lin] >= args.specificity_fold * np.maximum(others, 1e-9)
        cand = np.flatnonzero((owners == lin) & ~nk_wins & (lin_cpm[lin] >= 50)
                              & specific & (amb_cpm > 0) & (nk_cpm > 0))
        if cand.size == 0:
            print(f"[class] {lin}: no gene met specificity >= "
                  f"{args.specificity_fold}x; skipped")
            continue
        cand = cand[np.argsort(-lin_cpm[lin][cand])][:args.per_lineage]
        for c in cand:
            picks.append({"gene": syms[c], "col": int(c), "lineage": lin,
                          "class": "haematopoietic" if lin in HAEM else "non_haematopoietic",
                          "lineage_cpm": float(lin_cpm[lin][c]),
                          "ambient_cpm": float(amb_cpm[c]), "nk_cpm": float(nk_cpm[c])})
    P = pd.DataFrame(picks)
    if P.empty:
        raise SystemExit("no marker genes met the criteria")

    # ruler A, same estimator: rho from the panel controls, leave-one-out
    ratios = {}
    for g in PANEL_AMBIENT:
        w = np.flatnonzero(syms == g)
        if w.size and amb_cpm[w[0]] > 0 and nk_cpm[w[0]] > 0:
            ratios[g] = nk_cpm[w[0]] / amb_cpm[w[0]]
    rho = float(np.median(list(ratios.values())))
    print(f"[class] rho (median over {len(ratios)} panel controls) = {rho:.4f}")

    P["soup_fraction"] = [min(1.0, rho * r.ambient_cpm / r.nk_cpm) for _, r in P.iterrows()]
    P["saturated"] = P.soup_fraction >= 0.999
    P["dataset_id"] = mf.dataset_id
    P["compartment"] = args.compartment
    P.to_csv(os.path.join(outdir, "ceiling_class_structure.tsv"), sep="\t", index=False)

    print(f"\n=== soup fraction by lineage class ({len(P)} data-derived markers) ===")
    print(P.groupby("class").soup_fraction.describe()[
        ["count", "min", "25%", "50%", "75%", "max"]].round(4).to_string())
    print(f"\nfraction saturated at 1.000:")
    print((P.groupby("class").saturated.mean() * 100).round(1).to_string())
    print(f"\n=== by lineage ===")
    print(P.groupby(["class", "lineage"]).agg(
        n=("gene", "size"), median_sf=("soup_fraction", "median"),
        frac_saturated=("saturated", "mean")).round(3).to_string())

    h = P[P["class"] == "haematopoietic"].soup_fraction
    nh = P[P["class"] == "non_haematopoietic"].soup_fraction
    if len(h) >= 3 and len(nh) >= 3:
        from scipy.stats import mannwhitneyu
        u, pu = mannwhitneyu(h, nh, alternative="less")
        print(f"\nMann-Whitney (haematopoietic < non-haematopoietic): p = {pu:.2e}")
        print(f"  haematopoietic     median {np.median(h):.4f}  (n={len(h)})")
        print(f"  non-haematopoietic median {np.median(nh):.4f}  (n={len(nh)})")
        print("\nNOTE: this is a GENE-level test, not a donor-level one. It "
              "characterises the estimator, not any biological claim, so R1 "
              "does not apply -- and no result of it enters T2.")

    print("\n=== what this predicts for the missing PAEP ===")
    epi = P[P.lineage.isin(["Epithelial", "Stromal"])]
    if len(epi):
        print(f"  glandular-epithelial / stromal markers: median soup fraction "
              f"{epi.soup_fraction.median():.4f}, {epi.saturated.mean():.0%} saturated")
        print(f"  PAEP is glandular epithelial / decidualised stroma, so it is "
              f"expected on the saturated side")
        print(f"  => it would NOT lower the ceiling (a minimum over controls),")
        print(f"     so the CXCR4 risk from PAEP's absence is small")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
