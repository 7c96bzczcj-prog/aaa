"""Ambient correction at scale, with a held-out acceptance test.

CellBender is the principled tool and is not affordable here: 322 s per
epoch on 4 CPU cores, i.e. ~4.5 h for one library and ~300 h for all 73.
One library is being run separately as a benchmark; this is what can
actually be applied to the whole dataset.

The method is SoupX's non-expressed-gene estimator, made concrete by
this project's own finding. The soup profile is measured directly from
empty droplets (this dataset ships the full 737,280-barcode whitelist,
so the empty droplets are right there). The contamination fraction is
then read off genes a lineage cannot express — immunoglobulin in NK,
CD8 T, CD4 T and myeloid; granzymes in B — whose entire observed signal
is therefore soup.

Correction is applied at the (library x lineage x condition) level.
Summation is linear, so subtracting `rho * total * soup` from a summed
profile equals summing the per-cell corrections, provided rho is taken
as constant within the group.

ACCEPTANCE TEST, held out
-------------------------
Estimation and validation markers are disjoint. rho for NK is fitted on
immunoglobulin; the test is whether B-cell and myeloid genes
(MS4A1, CD79A, LYZ) also fall in NK. Nothing in the fit knows about
them, so their behaviour is an independent check that the correction
removed soup rather than fitting the markers it was given.
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

from nkmine.gating import (  # noqa: E402
    SOUP_ESTIMATION_MARKERS, SOUP_VALIDATION_MARKERS,
    decontaminate, estimate_contamination, read_10x_tar, soup_profile,
)
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "raw")
CACHE = os.path.join(ROOT, "data", "labels")
OUT = os.path.join(ROOT, "results")


def marker_frac(profile: np.ndarray, genes: np.ndarray, markers: list) -> float:
    """Fraction of a profile's counts sitting in `markers`."""
    idx = {g: i for i, g in enumerate(genes)}
    rows = [idx[m] for m in markers if m in idx]
    if not rows or profile.sum() <= 0:
        return float("nan")
    return float(profile[rows].sum() / profile.sum())


def one_library(batch: int):
    pk = os.path.join(CACHE, f"{batch}.pkl")
    tar = os.path.join(RAW, f"GSE154826_amp_batch_ID_{batch}.tar.gz")
    if not (os.path.exists(pk) and os.path.exists(tar)):
        return None
    with open(pk, "rb") as fh:
        lab = pickle.load(fh)

    prof, prof_genes, n_empty = soup_profile(tar)
    lib = read_10x_tar(tar)
    if list(lib.gene_names) != list(prof_genes):
        return None
    x = lib.rna.tocsc()

    rows, groups = [], {}
    for cond in ("Tumor", "Normal"):
        for lin in LINEAGES:
            m = np.flatnonzero((lab["lineage"] == lin) & (lab["condition"] == cond))
            if len(m) < 5:
                continue
            obs = np.asarray(x[:, m].sum(axis=1)).ravel().astype(float)
            rho = estimate_contamination(obs, lib.gene_names, prof,
                                         SOUP_ESTIMATION_MARKERS[lin])
            cor = decontaminate(obs, prof, rho)
            groups[(cond, lin)] = (cor, obs)   # keep raw for a matched A/B
            rows.append({
                "batch": batch, "condition": cond, "lineage": lin,
                "n_cells": int(len(m)), "n_empty_droplets": n_empty,
                "total_counts": float(obs.sum()), "rho": rho,
                "removed_frac": float(1 - cor.sum() / max(obs.sum(), 1)),
                # held-out validation markers, before and after
                "val_before": marker_frac(obs, lib.gene_names, SOUP_VALIDATION_MARKERS[lin]),
                "val_after": marker_frac(cor, lib.gene_names, SOUP_VALIDATION_MARKERS[lin]),
                "est_before": marker_frac(obs, lib.gene_names, SOUP_ESTIMATION_MARKERS[lin]),
                "est_after": marker_frac(cor, lib.gene_names, SOUP_ESTIMATION_MARKERS[lin]),
            })
    patient = lab.get("patient", "NA")
    return {"rows": rows, "groups": groups, "genes": lib.gene_names,
            "batch": batch, "patient": patient}


def main():
    batches = sorted(int(os.path.basename(p)[:-4])
                     for p in glob.glob(os.path.join(CACHE, "*.pkl")))
    print(f"decontaminating {len(batches)} libraries", flush=True)

    diag, store, genes = [], {}, None
    with Pool(3) as pool:
        for k, res in enumerate(pool.imap_unordered(one_library, batches), 1):
            if res is None:
                continue
            diag.extend(res["rows"])
            genes = res["genes"]
            for (cond, lin), (cor, obs) in res["groups"].items():
                store[(res["patient"], cond, lin, res["batch"])] = (
                    cor.astype(np.float32), obs.astype(np.float32))
            if k % 10 == 0:
                print(f"  {k}/{len(batches)}", flush=True)

    d = pd.DataFrame(diag)
    d.to_csv(os.path.join(OUT, "decontamination_diagnostics.csv"), index=False)

    # sum the per-library corrected profiles into (patient x condition x lineage)
    keys = sorted({k[:3] for k in store})
    cor_m = np.column_stack([
        np.sum([v[0] for kk, v in store.items() if kk[:3] == k], axis=0) for k in keys])
    raw_m = np.column_stack([
        np.sum([v[1] for kk, v in store.items() if kk[:3] == k], axis=0) for k in keys])
    meta = dict(genes=np.asarray(genes),
                patient=np.array([k[0] for k in keys]),
                condition=np.array([k[1] for k in keys]),
                lineage=np.array([k[2] for k in keys]))
    np.savez_compressed(os.path.join(OUT, "pseudobulk_decontaminated.npz"),
                        counts=cor_m, **meta)
    # identical aggregation, no correction: isolates the decontamination effect
    np.savez_compressed(os.path.join(OUT, "pseudobulk_uncorrected_matched.npz"),
                        counts=raw_m, **meta)
    print(f"\nwrote corrected + matched-uncorrected: {cor_m.shape[0]} genes x "
          f"{cor_m.shape[1]} samples", flush=True)

    print("\n=== contamination fraction rho, by lineage ===", flush=True)
    s = d.groupby("lineage").agg(n=("rho", "size"), rho_med=("rho", "median"),
                                 removed_med=("removed_frac", "median")).reindex(LINEAGES)
    print(s.round(4).to_string(), flush=True)

    print("\n=== ACCEPTANCE TEST (held-out markers, never used for fitting) ===",
          flush=True)
    v = d.groupby("lineage").agg(val_before=("val_before", "median"),
                                 val_after=("val_after", "median"),
                                 est_before=("est_before", "median"),
                                 est_after=("est_after", "median")).reindex(LINEAGES)
    v["val_drop"] = 1 - v.val_after / v.val_before
    v["est_drop"] = 1 - v.est_after / v.est_before
    print((v * 1).round(5).to_string(), flush=True)
    ok = bool((v.val_drop > 0.5).all())
    print(f"\nheld-out markers fell by >50% in every lineage: {ok}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
