#!/usr/bin/env python
"""Ruler B for PTN/OGN: what would carryover from decidual stroma actually look like?

GSE184719 shows that sorted decidual NK libraries carry stromal transcripts.
That establishes contamination but not its consequence: the question is whether
the *observed* PTN and OGN in those libraries is what that much contamination
predicts, or whether it exceeds it.

Ruler B in its v1.2 form: a gene's pickup in a target gate is
`f x source_lineage_CPM`, where `f` is the contaminating fraction, estimated
from a gene the target cannot transcribe.  This needs the source lineage,
which the Netskar NK-only atlas does not contain — but VT2018 does, as
dS1/dS2/dS3 decidual stromal cells from the same tissue.

So: estimate the PTN:collagen and OGN:collagen ratios in real first-trimester
decidual stroma here, and use them to predict the PTN/OGN that GSE184719's
measured collagen load implies.

**This is a cross-dataset prediction, not a merge (R9).** No count matrix is
combined.  What crosses the boundary is a within-lineage expression *ratio*
used as a prior, and it is stated as a prior with its platform stamp, not as a
measurement.  The two platforms differ (10x 3' v2 vs bulk poly-A), so the
prediction is order-of-magnitude and is read that way.

    python src/stromal_ruler.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr  # noqa: E402
from dnkchem.dataset import load_dataset  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

OUT = "out/GSE184719"

GENES = ["PTN", "OGN", "SPP1", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
         "PAEP", "VIM", "PTPRC"]

# measured in GSE184719, out/GSE184719/panel_cpm_per_library.csv
DNK_BULK = {
    "dNK1": {"PTN": 1.72, "OGN": 1.95, "COL1A1": 149.58, "COL3A1": 33.67,
             "DCN": 108.50, "LUM": 66.60},
    "dNK2": {"PTN": 15.08, "OGN": 52.09, "COL1A1": 388.39, "COL3A1": 493.94,
             "DCN": 1496.32, "LUM": 1807.66},
    "dNK3": {"PTN": 2.27, "OGN": 1.28, "COL1A1": 43.15, "COL3A1": 32.20,
             "DCN": 94.98, "LUM": 140.83},
    "dNK4": {"PTN": 0.57, "OGN": 0.71, "COL1A1": 28.74, "COL3A1": 45.45,
             "DCN": 44.39, "LUM": 57.70},
}
ANCHORS = ["COL1A1", "COL3A1", "DCN", "LUM"]


def cpm_by_group(ds, cols, groups):
    """Mean CPM per group, pooled over cells (a pseudobulk, which is what a
    sorted bulk library is)."""
    X = as_csr(ds.X)
    tot = np.asarray(X.sum(axis=1)).ravel()
    out = {}
    for name, mask in groups.items():
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            continue
        sub = X[idx][:, cols]
        gene_sum = np.asarray(sub.sum(axis=0)).ravel()
        out[name] = gene_sum / tot[idx].sum() * 1e6
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    ds = load_dataset(load_manifest(args.manifest))
    sym = np.asarray(ds.symbols)
    pos = {}
    for g in GENES:
        w = np.flatnonzero(sym == g)
        if w.size:
            pos[g] = int(w[0])
    missing = [g for g in GENES if g not in pos]
    rate = len(pos) / len(GENES)
    print(f"panel {len(pos)}/{len(GENES)} matched ({rate:.1%}); missing {missing}")
    if rate < 0.90:
        raise SystemExit(f"panel match rate {rate:.1%} < 90% (R11)")
    names = [g for g in GENES if g in pos]
    cols = np.array([pos[g] for g in names])

    obs = ds.obs
    dec = obs["compartment"].to_numpy() == "Decidua"
    # `celltype` keeps the published label; `subset` is the manifest's mapped
    # lineage, which is the one with a Stromal class.
    ct = obs["subset"].to_numpy()
    groups = {
        "Stromal": dec & (ct == "Stromal"),
        "Perivascular": dec & (ct == "Perivascular"),
        "Endothelial": dec & (ct == "Endothelial"),
        "Epithelial": dec & (ct == "Epithelial"),
        "Trophoblast": dec & (ct == "Trophoblast"),
        "dNK1": dec & (ct == "dNK1"),
        "dNK2": dec & (ct == "dNK2"),
        "dNK3": dec & (ct == "dNK3"),
    }
    print("\ncells per group:")
    for k, m in groups.items():
        print(f"  {k:14s} {int(m.sum()):6d}")

    cpm = pd.DataFrame(cpm_by_group(ds, cols, groups), index=names)
    print("\n=== VT2018 decidua, pseudobulk CPM ===")
    print(cpm.round(2).to_string())
    cpm.round(4).to_csv(f"{OUT}/vt2018_lineage_cpm.csv")

    stroma = cpm["Stromal"]
    print("\n=== ratio inside decidual stroma (the prior ruler B needs) ===")
    for t in ["PTN", "OGN"]:
        for a in ANCHORS:
            print(f"  {t}/{a} = {stroma[t] / stroma[a]:.4f}")

    # --- the platform-matched version, which needs no cross-dataset step ---
    # Ruler B proper: a gene's level in the NK gate as a fraction of its level
    # in the source lineage.  Genes NK cannot transcribe fix the carryover
    # envelope; a gene NK really transcribes must sit above it.
    print("\n=== ruler B inside VT2018 alone: NK level as % of decidual stroma ===")
    envelope = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]
    frac = pd.DataFrame(
        {arm: (cpm[arm] / cpm["Stromal"] * 100.0) for arm in ["dNK1", "dNK2", "dNK3"]}
    )
    print(frac.round(3).to_string())
    lo = frac.loc[envelope].to_numpy().min()
    hi = frac.loc[envelope].to_numpy().max()
    print(f"\ncarryover envelope from the 5 genes NK cannot transcribe: "
          f"{lo:.2f}% - {hi:.2f}% of stromal level")
    # Direction matters: only *above* the envelope is evidence of transcription.
    # Below it is carryover with room to spare, which is the same verdict as
    # inside — not a separate finding.
    for t in ["PTN", "OGN", "SPP1", "VIM", "PTPRC"]:
        vals = frac.loc[t].to_numpy()
        if (vals > hi).all():
            verdict = "ABOVE the envelope in every arm — carryover cannot account for it"
        elif (vals > hi).any():
            verdict = "above the envelope in some arms only"
        elif (vals < lo).all():
            verdict = "below the envelope — less than pure carryover would give"
        else:
            verdict = "inside the carryover envelope — no excess to explain"
        print(f"  {t:6s} {np.array2string(vals, precision=2)} %  -> {verdict}")
    frac.round(4).to_csv(f"{OUT}/vt2018_nk_over_stroma_pct.csv")

    # --- the prediction ----------------------------------------------------
    rows = []
    for lib, meas in DNK_BULK.items():
        for a in ANCHORS:
            f = meas[a] / stroma[a]          # contaminating stromal fraction
            for t in ["PTN", "OGN"]:
                pred = f * stroma[t]
                rows.append({
                    "library": lib, "anchor": a, "target": t,
                    "stromal_fraction_implied": f,
                    "predicted_cpm": pred,
                    "observed_cpm": meas[t],
                    "observed_over_predicted": meas[t] / pred if pred > 0 else np.nan,
                })
    pred = pd.DataFrame(rows)
    pred.to_csv(f"{OUT}/ptn_ogn_carryover_prediction.csv", index=False)

    print("\n=== carryover prediction, per anchor gene ===")
    piv = pred.pivot_table(index=["library", "target"], columns="anchor",
                           values="observed_over_predicted")
    print(piv.round(2).to_string())
    print("\n(observed / predicted: ~1 means the level is what this much stromal "
          "carryover already produces; >>1 means it exceeds carryover)")

    print("\nimplied stromal fraction per library (from each anchor):")
    fr = pred.pivot_table(index="library", columns="anchor",
                          values="stromal_fraction_implied")
    print((fr * 100).round(3).to_string())

    # --- rank concordance, with its floor ---------------------------------
    print("\n=== does PTN/OGN track the stromal load across the 4 libraries? ===")
    load = {k: sum(DNK_BULK[k][a] for a in ANCHORS) for k in DNK_BULK}
    order = sorted(load, key=lambda k: load[k], reverse=True)
    print("stromal load rank:", " > ".join(f"{k}({load[k]:.0f})" for k in order))
    for t in ["PTN", "OGN"]:
        o = sorted(DNK_BULK, key=lambda k: DNK_BULK[k][t], reverse=True)
        print(f"{t:4s} rank:        ", " > ".join(f"{k}({DNK_BULK[k][t]:.2f})" for k in o))
        print(f"      concordant with stromal load: {o == order}")
    n = len(DNK_BULK)
    from math import factorial
    print(f"\nn = {n} libraries. A perfect rank concordance has exact two-sided "
          f"p = 2/{factorial(n)} = {2 / factorial(n):.4f} — the floor of the test "
          f"(R2). on_test_floor = True. This ordering cannot reach p < 0.05 at "
          f"n = 4 and is reported as an ordering, not as a significant result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
