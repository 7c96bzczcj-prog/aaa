"""Ten-minute premise check: does NK express intracellular complement C3?

The complosome hypothesis needs NK to transcribe `C3` (and to carry the
receptor `C3AR1`) in the first place. That premise is checkable on data
already on disk and has never been run.

Why this needs more than a detection call. C3 is a canonical myeloid /
stromal transcript, so any NK signal is a prime candidate for ambient
contamination -- and this project has already measured that NK carries
the second-highest soup burden of any lineage. The check therefore
reports three things per gene:

  detection    fraction of that lineage's pseudobulk samples with any count
  level        mean CPM, so NK can be read against myeloid
  soup_frac    rho * total * soup_g / observed_g, the same estimator used
               in scripts/ambient_regression.py

`soup_frac` near 1.0 means the lineage's entire observed signal for that
gene is accounted for by ambient, i.e. detection is not evidence of
expression. Two positive controls (KLRD1, NKG7 -- genuinely NK) and two
negative controls (C1QA, LYZ -- myeloid, and known to leak) calibrate the
scale, so a low soup fraction cannot be dismissed as a quirk of the
estimator.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.pseudobulk import PseudobulkSet  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")

TARGETS = ["C3", "C3AR1", "C5AR1", "CD46", "CTSL", "CFB"]
POS_CTRL = ["KLRD1", "NKG7", "GNLY"]      # genuinely NK
NEG_CTRL = ["C1QA", "LYZ", "IGKC"]        # myeloid / B, known to leak


def main():
    d = np.load(os.path.join(OUT, "pseudobulk_uncorrected_matched.npz"),
                allow_pickle=True)
    ps = PseudobulkSet(d["counts"].astype(float), d["genes"].astype(str),
                       d["patient"], d["condition"], d["lineage"],
                       np.ones(len(d["patient"]))).complete_pairs()
    sp = np.load(os.path.join(OUT, "soup_profile_mean.npz"), allow_pickle=True)
    soup_all = dict(zip(sp["genes"].astype(str), sp["soup"]))
    rho_lin = (pd.read_csv(os.path.join(OUT, "rho_variants.csv"))
               .groupby("lineage").rho_per_sample.median().to_dict())

    genes = list(ps.genes)
    panel = TARGETS + POS_CTRL + NEG_CTRL
    idx = {g: genes.index(g) for g in panel if g in genes}
    missing = [g for g in panel if g not in idx]
    if missing:
        print(f"not in matrix at all: {missing}", flush=True)

    rows = []
    for lin in LINEAGES:
        s = ps.subset_lineage(lin)
        lib = s.counts.sum(axis=0)
        cpm = s.counts / np.maximum(lib, 1) * 1e6
        total = float(lib.mean())
        rho = float(rho_lin.get(lin, 0.1))
        for g, i in idx.items():
            obs_mean = float(s.counts[i].mean())
            f = rho * total * soup_all.get(g, 0.0) / max(obs_mean, 1e-9)
            rows.append({
                "gene": g, "lineage": lin,
                "detection": float((s.counts[i] > 0).mean()),
                "mean_CPM": float(cpm[i].mean()),
                "soup_frac": float(min(f, 1.0)),
            })
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(OUT, "complosome_premise.csv"), index=False)

    def show(title, gl, col, fmt):
        print(f"\n=== {title} ===", flush=True)
        p = r[r.gene.isin(gl)].pivot(index="gene", columns="lineage", values=col)
        print(p.reindex([g for g in gl if g in idx])[list(LINEAGES)]
              .round(fmt).to_string(), flush=True)

    show("detection rate (fraction of pseudobulk samples with any count)",
         panel, "detection", 3)
    show("mean CPM", panel, "mean_CPM", 2)
    show("estimated soup fraction (1.0 = signal fully explained by ambient)",
         panel, "soup_frac", 3)

    # ---- calibration, because the raw soup fraction cannot be read at face
    # value.  Run 2 of the ambient work established that in NK this estimator
    # reads ~0.29 for genes whose TRUE soup fraction is 1.0 (immunoglobulin).
    # So "0.19" does not mean "19% ambient"; it means "below the point where
    # a fully-ambient gene reads".  Both ends of the scale are measured here
    # from genes in the same lineage rather than assumed.
    nkr = r[r.lineage == "NK"].set_index("gene")
    myr = r[r.lineage == "Myeloid"].set_index("gene")
    soup_ceiling = float(nkr.loc[[g for g in NEG_CTRL if g in idx],
                                 "soup_frac"].min())
    real_floor = float(nkr.loc[[g for g in POS_CTRL if g in idx],
                               "soup_frac"].max())
    print(f"\n=== calibration inside NK ===", flush=True)
    print(f"  genuinely-NK genes read up to soup_frac {real_floor:.3f}", flush=True)
    print(f"  known-ambient genes read down to soup_frac {soup_ceiling:.3f} "
          f"(their true value is 1.0)", flush=True)

    # A second, calibration-free discriminator: for a gene NK does not
    # express, NK's CPM is ambient pickup, so NK/Myeloid should sit at the
    # ratio seen for genes only myeloid makes.
    band = [float(nkr.loc[g, "mean_CPM"] / max(myr.loc[g, "mean_CPM"], 1e-9))
            for g in ("C1QA", "LYZ") if g in idx]
    print(f"  NK/Myeloid CPM for myeloid-only genes (pure pickup): "
          f"{'  '.join(f'{b:.3f}' for b in band)}", flush=True)

    print("\n=== VERDICT ===", flush=True)
    for g in TARGETS:
        if g not in idx:
            continue
        nk, my = nkr.loc[g], myr.loc[g]
        ratio = float(nk.mean_CPM / max(my.mean_CPM, 1e-9))
        if nk.detection < 0.5:
            call = "BELOW DETECTION"
        elif nk.soup_frac >= soup_ceiling:
            call = "INDISTINGUISHABLE FROM AMBIENT"
        elif ratio <= max(band) * 1.5:
            call = "INDISTINGUISHABLE FROM AMBIENT (ratio)"
        else:
            call = "above ambient on both scales"
        print(f"  {g:6s} det={nk.detection:.2f} CPM={nk.mean_CPM:7.2f} "
              f"soup={nk.soup_frac:.3f} NK/Mye={ratio:6.3f} -> {call}",
              flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
