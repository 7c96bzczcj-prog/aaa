"""Phase 2 across every admitted GSE154826 library.

Produces the two things Phase 1 could not settle from metadata alone:
the per-(individual x condition) NK cell counts that decide admission
criterion A4, and the purity report.

Hashed libraries pool two conditions in one droplet emulsion, so their
cells are demultiplexed by HTO and split; unhashed libraries carry a
single condition taken from the GEO sample annotation.
"""

from __future__ import annotations

import collections
import glob
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.gating import (  # noqa: E402
    demux_hto,
    excluded_genes,
    gate_lineages_rna,
    purity_report,
    read_10x_tar,
)

ROOT = os.path.join(os.path.dirname(__file__), "..")
ANNOT = os.path.join(ROOT, "phase1_registry", "raw_meta", "GSE154826_sample_annots.csv")
RAW = os.path.join(ROOT, "data", "raw")

# T-lineage transcripts deliberately NOT used by the gate, so they can
# serve as an independent contamination probe inside the NK gate.
T_PROBE = ["LCK", "CD2", "THEMIS", "SKAP1", "TRAT1", "ITK"]


def hto_map(annots: pd.DataFrame, batch: int) -> dict:
    """HTO tag number -> tissue, for one library."""
    rows = annots[annots.amp_batch_ID == batch]
    out = {}
    for _, r in rows.iterrows():
        if pd.isna(r.HTO):
            continue
        for tag in str(r.HTO).split(","):
            tag = tag.strip()
            if tag:
                out[f"HTO_{tag}"] = r.tissue
    return out


def process(batch: int) -> dict:
    annots = pd.read_csv(ANNOT)
    path = os.path.join(RAW, f"GSE154826_amp_batch_ID_{batch}.tar.gz")
    rows = annots[annots.amp_batch_ID == batch]
    patients = sorted(set(rows.patient_ID.astype(str)))
    tissues = sorted(set(rows.tissue))

    lib = read_10x_tar(path)
    lineage, _ = gate_lineages_rna(lib)
    pur = purity_report(lib, lineage)

    # --- independent contamination probe inside the NK gate ----------
    g2i = {}
    for i, g in enumerate(lib.gene_names):
        g2i.setdefault(g, []).append(i)
    probe_rows = [i for g in T_PROBE for i in g2i.get(g, [])]
    probe_frac = {}
    if probe_rows:
        pm = np.asarray(lib.rna[probe_rows].todense()) > 0
        any_probe = pm.any(axis=0)
        for l in ("NK", "CD8T", "B", "Myeloid"):
            m = lineage == l
            probe_frac[l] = float(any_probe[m].mean()) if m.sum() else np.nan

    # --- assign each cell a condition --------------------------------
    hmap = hto_map(annots, batch)
    if lib.hto is not None and len(hmap) > 0 and len(tissues) > 1:
        calls = demux_hto(lib)
        cond = np.array([hmap.get(c, "ambiguous") for c in calls], dtype=object)
        mode = "hashed"
    else:
        cond = np.array([tissues[0]] * lib.rna.shape[1], dtype=object)
        mode = "single"

    patient = patients[0] if len(patients) == 1 else "|".join(patients)

    counts = []
    for c in sorted(set(cond)):
        for l in ("NK", "CD8T", "CD4T", "B", "Myeloid", "unassigned"):
            n = int(((cond == c) & (lineage == l)).sum())
            counts.append({"batch": batch, "patient": patient, "condition": c,
                           "lineage": l, "n_cells": n, "mode": mode})

    return {
        "counts": counts,
        "purity": pur.assign(batch=batch, patient=patient,
                             probe_pos_frac=[probe_frac.get(l, np.nan) for l in pur.lineage])
        if len(pur) else pd.DataFrame(),
        "meta": {"batch": batch, "patient": patient, "mode": mode,
                 "n_cells": int(lib.rna.shape[1]), "n_droplets": int(lib.n_droplets),
                 "has_adt": lib.adt is not None,
                 "med_umi": float(np.median(np.asarray(lib.rna.sum(axis=0)).ravel()))},
    }


def main():
    batches = [int(x) for x in open(
        os.path.join(ROOT, "phase1_registry", "batches_to_fetch.txt")).read().split()]
    have = {int(os.path.basename(p).split("_ID_")[1].split(".")[0])
            for p in glob.glob(os.path.join(RAW, "*.tar.gz"))}
    batches = [b for b in batches if b in have]
    print(f"processing {len(batches)} libraries", flush=True)

    with Pool(4) as pool:
        results = pool.map(process, batches)

    counts = pd.DataFrame([r for res in results for r in res["counts"]])
    purity = pd.concat([res["purity"] for res in results if len(res["purity"])],
                       ignore_index=True)
    meta = pd.DataFrame([res["meta"] for res in results])

    out = os.path.join(ROOT, "results")
    os.makedirs(out, exist_ok=True)
    counts.to_csv(os.path.join(out, "cell_counts_by_group.csv"), index=False)
    purity.to_csv(os.path.join(out, "purity_report_GSE154826.csv"), index=False)
    meta.to_csv(os.path.join(out, "library_meta.csv"), index=False)
    with open(os.path.join(ROOT, "excluded_genes.txt"), "w") as fh:
        fh.write("\n".join(excluded_genes() + T_PROBE) + "\n")
    print("WROTE results/cell_counts_by_group.csv", flush=True)


if __name__ == "__main__":
    main()
