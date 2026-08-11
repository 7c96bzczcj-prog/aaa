"""One pass over the libraries producing several rho estimators.

Pre-registered question for this run, fixed before it was launched:

    with rho pooled so that it is no longer a per-sample quantity,
    does the standard error return to baseline?

    SE back to ~0.12 (NK) / ~0.09 (others), markers still bad
        -> mechanism (a) is fixed, (b) is the binding constraint,
           stop tuning subtraction and move to the regression approach
    SE not back
        -> diagnosis (a) was wrong, stop this line entirely
    both good
        -> unexpected; continue

Three estimators are computed in the same pass, since the cost is
dominated by reading the tarballs:

  per_sample   rho per (library x lineage x condition), 3 markers.
               This is what failed; kept as the reference arm.
  pooled_cond  rho per (library x lineage), both conditions together.
               Removes the per-condition sampling noise. Note what it
               also removes: ambient only distorts a *paired* contrast
               when rho differs between conditions, so pooling across
               conditions discards the very term that does the damage.
               It is run because it isolates the variance question.
  wide_markers rho per (library x lineage x condition) as before, but
               estimated on ~20 lineage-foreign genes instead of 3.
               This keeps the differential term and attacks the noise
               directly.

The mean soup profile is also saved, for the regression approach that
does not touch counts at all.
"""

from __future__ import annotations

import glob
import os
import pickle
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.gating import decontaminate, estimate_contamination, read_10x_tar, soup_profile  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "raw")
CACHE = os.path.join(ROOT, "data", "labels")
OUT = os.path.join(ROOT, "results")

NARROW = {
    "NK": ["IGKC", "IGHG1", "IGHA1"], "CD8T": ["IGKC", "IGHG1", "IGHA1"],
    "CD4T": ["IGKC", "IGHG1", "IGHA1"], "Myeloid": ["IGKC", "IGHG1", "IGHA1"],
    "B": ["GZMB", "GZMA", "GZMH"],
}
# ~20 genes each lineage should not express, for a lower-variance rho
_IG = ["IGKC", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM",
       "IGLC1", "IGLC2", "IGLC3", "IGLC7", "IGKV3-20", "JCHAIN", "MZB1", "DERL3",
       "TNFRSF17", "SDC1", "XBP1", "PRDX4"]
_EPI = ["SFTPC", "SFTPB", "SFTPA1", "SCGB1A1", "SCGB3A2", "NAPSA", "KRT19",
        "KRT18", "KRT8", "EPCAM", "AGER", "CLDN18", "MUC1", "SLPI", "CAV1"]
WIDE = {
    "NK": _IG + _EPI, "CD8T": _IG + _EPI, "CD4T": _IG + _EPI,
    "Myeloid": _IG + _EPI,
    "B": ["GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "NKG7", "KLRD1", "GNLY",
          "KLRF1", "CTSW", "CD3D", "CD3E", "CD3G", "TRAC", "CD8A"] + _EPI,
}


def one_library(batch: int):
    pk = os.path.join(CACHE, f"{batch}.pkl")
    tar = os.path.join(RAW, f"GSE154826_amp_batch_ID_{batch}.tar.gz")
    if not (os.path.exists(pk) and os.path.exists(tar)):
        return None
    with open(pk, "rb") as fh:
        lab = pickle.load(fh)
    prof, prof_genes, _ = soup_profile(tar)
    lib = read_10x_tar(tar)
    if list(lib.gene_names) != list(prof_genes):
        return None
    x = lib.rna.tocsc()

    obs = {}
    for cond in ("Tumor", "Normal"):
        for lin in LINEAGES:
            m = np.flatnonzero((lab["lineage"] == lin) & (lab["condition"] == cond))
            if len(m) >= 5:
                obs[(cond, lin)] = np.asarray(x[:, m].sum(axis=1)).ravel().astype(float)

    out, rows = {}, []
    for lin in LINEAGES:
        pair = [(c, obs[(c, lin)]) for c in ("Tumor", "Normal") if (c, lin) in obs]
        if not pair:
            continue
        pooled = np.sum([v for _, v in pair], axis=0)
        rho_pool = estimate_contamination(pooled, lib.gene_names, prof, NARROW[lin])
        for cond, o in pair:
            r_ps = estimate_contamination(o, lib.gene_names, prof, NARROW[lin])
            r_wd = estimate_contamination(o, lib.gene_names, prof, WIDE[lin])
            out[("per_sample", cond, lin)] = decontaminate(o, prof, r_ps)
            out[("pooled_cond", cond, lin)] = decontaminate(o, prof, rho_pool)
            out[("wide_markers", cond, lin)] = decontaminate(o, prof, r_wd)
            rows.append({"batch": batch, "condition": cond, "lineage": lin,
                         "rho_per_sample": r_ps, "rho_pooled": rho_pool,
                         "rho_wide": r_wd})
    return {"out": out, "rows": rows, "genes": lib.gene_names,
            "patient": lab.get("patient", "NA"), "soup": prof,
            "weight": float(sum(v.sum() for v in obs.values()))}


def main():
    batches = sorted(int(os.path.basename(p)[:-4])
                     for p in glob.glob(os.path.join(CACHE, "*.pkl")))
    print(f"rho variants over {len(batches)} libraries", flush=True)
    store, rows, genes = {}, [], None
    soup_acc, soup_w = None, 0.0
    with Pool(3) as pool:
        for k, r in enumerate(pool.imap_unordered(one_library, batches), 1):
            if r is None:
                continue
            genes = r["genes"]
            rows.extend(r["rows"])
            soup_acc = r["soup"] * r["weight"] if soup_acc is None else soup_acc + r["soup"] * r["weight"]
            soup_w += r["weight"]
            for (variant, cond, lin), vec in r["out"].items():
                store.setdefault(variant, {}).setdefault(
                    (r["patient"], cond, lin), []).append(vec.astype(np.float32))
            if k % 15 == 0:
                print(f"  {k}/{len(batches)}", flush=True)

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "rho_variants.csv"), index=False)
    soup = soup_acc / soup_w
    np.savez_compressed(os.path.join(OUT, "soup_profile_mean.npz"),
                        soup=soup, genes=np.asarray(genes))

    for variant, groups in store.items():
        keys = sorted(groups)
        counts = np.column_stack([np.sum(groups[k], axis=0) for k in keys])
        np.savez_compressed(
            os.path.join(OUT, f"pseudobulk_{variant}.npz"),
            counts=counts, genes=np.asarray(genes),
            patient=np.array([k[0] for k in keys]),
            condition=np.array([k[1] for k in keys]),
            lineage=np.array([k[2] for k in keys]))
        print(f"  wrote pseudobulk_{variant}.npz {counts.shape}", flush=True)

    d = pd.DataFrame(rows)
    print("\n=== rho by estimator (median) ===", flush=True)
    print(d.groupby("lineage")[["rho_per_sample", "rho_pooled", "rho_wide"]]
          .median().reindex(LINEAGES).round(4).to_string(), flush=True)
    print("\n=== within-(library,lineage) SD of rho across the two conditions ===",
          flush=True)
    sd = (d.groupby(["batch", "lineage"])[["rho_per_sample", "rho_wide"]]
            .std().groupby("lineage").median().reindex(LINEAGES))
    print(sd.round(4).to_string(), flush=True)
    print("  (this is the per-sample noise that gets injected into every gene)",
          flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
