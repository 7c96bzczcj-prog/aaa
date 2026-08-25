#!/usr/bin/env python3
"""E1 -- magnitude reference for the adaptive-NK BM-PB difference, on the one
admitted same-donor paired dataset (GSE233304).

This is a MAGNITUDE READING ON PUBLIC DATA, not a result and not a prediction
about anything else.

Frozen before the run:

  gate (the T3 rule, applied to every candidate gene)
      det_NK inside 5%-85%                     (admission criterion 2)
      det_NK - floor >= 0.05                   (the gap itself clears the floor)
      det_NK >= 3 x floor                      (signal at least 3x ambient)
    floor = the gene's detection rate in cells that do not express it, measured
    in the same libraries.  B cells are the floor population for every gene:
    they express none of the eight, whereas myeloid cells express FCER1G and
    CX3CR1 and so cannot serve.  Where erythroid cells are numerous enough their
    rate is reported too and the LARGER of the two is used, which is the
    conservative choice.
    Fewer than 3 genes through the gate -> E1 is unmeasurable, full stop.

  scoring (depth-symmetric, per T1)
      score = mean log1p(CP10K) over passing HIGH-direction genes
            - mean log1p(CP10K) over passing LOW-direction genes
      adaptive if score > 0.  Both terms are measured in the same cell, so
      depth and ambient level cancel.  A presence/absence variant is reported
      as a sensitivity arm; it is NOT primary, because an absolute zero-count
      requirement is the failure mode T1 records.

  statistics
      unit = donor.  Delta = marrow fraction - blood fraction, per donor.
      No significance test: n is too small and a p value would be misread.
"""
from __future__ import annotations

import json, os, sys

import numpy as np, pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

RES = "/home/user/aaa/results"
HIGH = ["KLRC2", "B3GAT1", "ZEB2", "CX3CR1"]          # adaptive-direction up
LOW = ["FCER1G", "KLRC1", "IL7R", "GZMK"]             # adaptive-direction down
CANDIDATES = HIGH + LOW
FLOOR_POPS = ["B", "Erythroid"]
GATE_DET_LO, GATE_DET_HI = pf.DETECTION_LO, pf.DETECTION_HI
GATE_GAP = 0.05
GATE_FOLD = 3.0
MIN_GENES = 3
MIN_NK = pf.MIN_NK_CELLS


def main():
    idx = json.load(open("index/gse233304.json"))
    cache, per_lib = {}, []
    for e in idx:
        if e["labels"] not in cache:
            z = np.load(e["labels"], allow_pickle=True)
            cache[e["labels"]] = (z["barcode"], z["library"], z["lineage"])
        bc, libcol, lin = cache[e["labels"]]
        sel = libcol == e["library"]
        mat = pf.read_10x_mtx(*e["path"].split("|"), e["library"])
        g2i = {}
        for i, g in enumerate(mat.genes):
            g2i.setdefault(g, []).append(i)
        pos = {b: i for i, b in enumerate(mat.barcodes)}
        rec = {"donor": e["donor"], "compartment": e["compartment"],
               "library": e["library"], "group": e["group"]}
        cols = {}
        for pop in ["NK"] + FLOOR_POPS:
            m = sel & (lin == pop)
            ix = [pos[b] for b in bc[m] if b in pos]
            cols[pop] = np.array(sorted(ix), dtype=int)
            rec[f"n_{pop}"] = len(ix)
        for g in CANDIDATES:
            rows = g2i.get(g)
            for pop in ["NK"] + FLOOR_POPS:
                c = cols[pop]
                if rows is None or len(c) < MIN_NK:
                    rec[f"det_{pop}_{g}"] = np.nan
                    continue
                v = np.asarray(mat.mat[rows][:, c].sum(axis=0)).ravel()
                rec[f"det_{pop}_{g}"] = float((v > 0).mean())
        # per-cell normalised expression for the NK gate, kept for scoring
        c = cols["NK"]
        if len(c) >= MIN_NK:
            sub = mat.mat[:, c]
            libsize = np.asarray(sub.sum(axis=0)).ravel().astype(float)
            libsize[libsize == 0] = 1.0
            for g in CANDIDATES:
                rows = g2i.get(g)
                if rows is None:
                    continue
                x = np.asarray(sub[rows].sum(axis=0)).ravel()
                rec[f"expr_{g}"] = json.dumps(
                    np.log1p(x / libsize * 1e4).round(4).tolist())
                rec[f"cnt_{g}"] = json.dumps(x.astype(int).tolist())
        per_lib.append(rec)
        del mat
        print(f"  {e['donor']:<5}{e['compartment']:<5}{e['library']:<5} "
              f"NK={rec['n_NK']:<5} B={rec['n_B']:<5} Ery={rec['n_Erythroid']}",
              flush=True)
    df = pd.DataFrame(per_lib)

    # ---- the gate -----------------------------------------------------------
    use = df[df.n_NK >= MIN_NK]
    gate = []
    for g in CANDIDATES:
        w = use.n_NK.astype(float)
        det = float(np.average(use[f"det_NK_{g}"], weights=w)) \
            if use[f"det_NK_{g}"].notna().all() else np.nan
        floors = {}
        for pop in FLOOR_POPS:
            sub = df[(df[f"n_{pop}"] >= MIN_NK) & df[f"det_{pop}_{g}"].notna()]
            floors[pop] = (float(np.average(sub[f"det_{pop}_{g}"],
                                            weights=sub[f"n_{pop}"].astype(float)))
                           if len(sub) else np.nan)
        fl = np.nanmax(list(floors.values()))
        ok = (GATE_DET_LO <= det <= GATE_DET_HI) and (det - fl >= GATE_GAP) \
            and (det >= GATE_FOLD * fl)
        gate.append({"gene": g, "direction": "high" if g in HIGH else "low",
                     "det_NK": det, "floor_B": floors["B"],
                     "floor_Erythroid": floors["Erythroid"], "floor_used": fl,
                     "gap": det - fl, "fold": det / fl if fl > 0 else np.inf,
                     "in_5_85_window": GATE_DET_LO <= det <= GATE_DET_HI,
                     "gap_ge_0.05": det - fl >= GATE_GAP,
                     "fold_ge_3": det >= GATE_FOLD * fl,
                     "passes": bool(ok)})
    gt = pd.DataFrame(gate)
    print("\n=== T3 gate on the adaptive-NK marker panel ===")
    print(gt.to_string(index=False))
    passed = gt[gt.passes]
    ph = [g for g in passed.gene if g in HIGH]
    pl = [g for g in passed.gene if g in LOW]
    print(f"\npassing genes: {len(passed)}  high={ph}  low={pl}")

    out = {"gate": gt.to_dict("records"), "n_passing": int(len(passed)),
           "high_passing": ph, "low_passing": pl}
    if len(passed) < MIN_GENES or not ph or not pl:
        out["verdict"] = "E1_UNMEASURABLE"
        gt.to_csv(os.path.join(RES, "e1_adaptive_delta_gse233304.tsv"),
                  sep="\t", index=False)
        json.dump(out, open(os.path.join(RES, "e1_gate.json"), "w"), indent=1)
        print("\nE1 UNMEASURABLE: fewer than 3 genes clear the gate, or one "
              "direction is empty. Stopping as the frozen table requires.")
        return

    # ---- per-donor fractions ------------------------------------------------
    rows = []
    for _, r in use.iterrows():
        eh = np.mean([np.array(json.loads(r[f"expr_{g}"])) for g in ph], axis=0)
        el = np.mean([np.array(json.loads(r[f"expr_{g}"])) for g in pl], axis=0)
        score = eh - el
        ch = np.sum([np.array(json.loads(r[f"cnt_{g}"])) for g in ph], axis=0)
        cl = np.sum([np.array(json.loads(r[f"cnt_{g}"])) for g in pl], axis=0)
        rows.append({"donor": r.donor, "compartment": r.compartment,
                     "library": r.library, "group": r.group, "n_NK": r.n_NK,
                     "n_adaptive_score": int((score > 0).sum()),
                     "frac_adaptive_score": float((score > 0).mean()),
                     "n_adaptive_strict": int(((ch > 0) & (cl == 0)).sum()),
                     "frac_adaptive_strict": float(((ch > 0) & (cl == 0)).mean())})
    lib = pd.DataFrame(rows)
    dn = (lib.groupby(["donor", "compartment"])
             .apply(lambda g: pd.Series({
                 "group": g.group.iloc[0], "n_NK": int(g.n_NK.sum()),
                 "frac_adaptive_score": float(np.average(g.frac_adaptive_score,
                                                         weights=g.n_NK)),
                 "frac_adaptive_strict": float(np.average(g.frac_adaptive_strict,
                                                          weights=g.n_NK))}),
                    include_groups=False).reset_index())
    bm = dn[dn.compartment == "BM"].set_index("donor")
    pb = dn[dn.compartment == "PB"].set_index("donor")
    shared = sorted(set(bm.index) & set(pb.index))
    delta = pd.DataFrame({
        "donor": shared, "group": [bm.loc[d, "group"] for d in shared],
        "n_NK_BM": [int(bm.loc[d, "n_NK"]) for d in shared],
        "n_NK_PB": [int(pb.loc[d, "n_NK"]) for d in shared],
        "frac_BM_score": [bm.loc[d, "frac_adaptive_score"] for d in shared],
        "frac_PB_score": [pb.loc[d, "frac_adaptive_score"] for d in shared],
        "frac_BM_strict": [bm.loc[d, "frac_adaptive_strict"] for d in shared],
        "frac_PB_strict": [pb.loc[d, "frac_adaptive_strict"] for d in shared]})
    delta["delta_score"] = delta.frac_BM_score - delta.frac_PB_score
    delta["delta_strict"] = delta.frac_BM_strict - delta.frac_PB_strict

    with open(os.path.join(RES, "e1_adaptive_delta_gse233304.tsv"), "w") as fh:
        fh.write("# E1 -- magnitude reference on public data only.\n")
        fh.write("# section 1: the T3 gate on every candidate gene\n")
        gt.to_csv(fh, sep="\t", index=False)
        fh.write("\n# section 2: per-donor adaptive fraction and Delta = BM - PB\n")
        delta.to_csv(fh, sep="\t", index=False)
    d = delta.delta_score
    out.update({"verdict": "E1_MEASURED", "n_donors": len(delta),
                "delta_median": float(d.median()),
                "delta_q1": float(d.quantile(0.25)),
                "delta_q3": float(d.quantile(0.75)),
                "n_positive": int((d > 0).sum()), "n_negative": int((d < 0).sum()),
                "iqr_spans_zero": bool(d.quantile(0.25) < 0 < d.quantile(0.75))})
    json.dump(out, open(os.path.join(RES, "e1_gate.json"), "w"), indent=1)
    print("\n=== per-donor Delta (BM - PB), adaptive fraction of NK ===")
    print(delta.to_string(index=False))
    print(f"\nmedian {d.median():+.4f}  IQR [{d.quantile(0.25):+.4f}, "
          f"{d.quantile(0.75):+.4f}]  sign: {int((d>0).sum())} positive / "
          f"{int((d<0).sum())} negative of {len(d)}")


if __name__ == "__main__":
    main()
