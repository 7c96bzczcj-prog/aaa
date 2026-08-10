"""Phases 3-5 on real GSE154826 pseudobulks.

Two passes over the libraries: the first assigns lineage and condition
to every called cell and records the labels; the second, once the global
per-(patient x condition) minimum lineage size is known, re-reads each
library and sums the sampled cells into balanced pseudobulks.

Two passes rather than one because Phase 2.4 balancing needs a quantity
(the rarest lineage in a patient-condition group) that is only defined
after every library contributing to that group has been seen.
"""

from __future__ import annotations

import os
import sys
import glob
import pickle
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.gating import (  # noqa: E402
    demux_hto, excluded_genes, gate_lineages_rna, read_10x_tar,
)
from nkmine.pseudobulk import PseudobulkSet, detection_filter, se_balance_report  # noqa: E402
from nkmine.de import paired_de  # noqa: E402
from nkmine.null_calibration import calibrate, cross_lineage_noise_summary  # noqa: E402
from nkmine.quadrant import LINEAGES, classify_table  # noqa: E402
from nkmine.controls import control_check, evaluate_stop_rule_s4  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts_common import hto_map  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "raw")
CACHE = os.path.join(ROOT, "data", "labels")
OUT = os.path.join(ROOT, "results")
N_REPS = 5
MIN_CELLS = 30


def pass1(batch: int):
    """Gate one library and cache its per-cell labels."""
    dst = os.path.join(CACHE, f"{batch}.pkl")
    if os.path.exists(dst):
        return dst
    annots = pd.read_csv(os.path.join(ROOT, "phase1_registry", "raw_meta",
                                      "GSE154826_sample_annots.csv"))
    rows = annots[annots.amp_batch_ID == batch]
    tissues = sorted(set(rows.tissue))
    pats = [str(p) for p in sorted(set(rows.patient_ID.astype(str))) if str(p).isdigit()]

    lib = read_10x_tar(os.path.join(RAW, f"GSE154826_amp_batch_ID_{batch}.tar.gz"))
    lineage, _ = gate_lineages_rna(lib)
    hmap = hto_map(annots, batch)
    if lib.hto is not None and hmap and len(tissues) > 1:
        cond = np.array([hmap.get(c, "ambiguous") for c in demux_hto(lib)], dtype=object)
    else:
        cond = np.array([tissues[0]] * lib.rna.shape[1], dtype=object)

    with open(dst, "wb") as fh:
        pickle.dump({"batch": batch, "patient": pats[0] if pats else "NA",
                     "lineage": lineage, "condition": cond.astype(str),
                     "genes": lib.gene_names}, fh)
    return dst


def pass2(args):
    """Sum the selected cells of one library into pseudobulk vectors."""
    batch, picks = args
    lib = read_10x_tar(os.path.join(RAW, f"GSE154826_amp_batch_ID_{batch}.tar.gz"))
    x = lib.rna.tocsc()
    out = {}
    for key, idx in picks.items():
        if len(idx) == 0:
            continue
        out[key] = np.asarray(x[:, np.asarray(idx)].sum(axis=1)).ravel().astype(np.float32)
    return out


def main():
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    batches = sorted(int(os.path.basename(p).split("_ID_")[1].split(".")[0])
                     for p in glob.glob(os.path.join(RAW, "*.tar.gz")))

    print(f"pass 1: gating {len(batches)} libraries", flush=True)
    with Pool(4) as pool:
        pool.map(pass1, batches)

    labels = {}
    genes = None
    for b in batches:
        with open(os.path.join(CACHE, f"{b}.pkl"), "rb") as fh:
            labels[b] = pickle.load(fh)
            genes = labels[b]["genes"]

    # ---- global per-(patient x condition x lineage) inventory -------
    inv = {}
    for b, d in labels.items():
        for i, (l, c) in enumerate(zip(d["lineage"], d["condition"])):
            if c not in ("Tumor", "Normal") or l not in LINEAGES:
                continue
            inv.setdefault((d["patient"], c, l), []).append((b, i))

    groups = sorted({(p, c) for (p, c, _) in inv})
    rng = np.random.default_rng(0)
    keep_groups = []
    for (p, c) in groups:
        sizes = {l: len(inv.get((p, c, l), [])) for l in LINEAGES}
        if min(sizes.values()) >= MIN_CELLS:
            keep_groups.append(((p, c), min(sizes.values())))
    print(f"pass 2: {len(keep_groups)}/{len(groups)} groups clear "
          f"{MIN_CELLS} cells in every lineage", flush=True)

    # Balanced sampling: median over N_REPS subsamples.  All reps are
    # drawn up front and summed in a single pass over the libraries --
    # re-reading each tarball once per rep would cost 10x73 decompressions.
    picks_by_batch = {}
    for r in range(N_REPS):
        for (p, c), n_min in keep_groups:
            for l in LINEAGES:
                cells = inv[(p, c, l)]
                sel = rng.choice(len(cells), size=n_min, replace=False)
                for j in sel:
                    b, i = cells[j]
                    picks_by_batch.setdefault(b, {}).setdefault((r, p, c, l), []).append(i)

    # Stream the per-library results in and fold them into the running
    # accumulator one at a time.  Holding all 73 result dicts at once
    # (each a few hundred float64 vectors of length ~34k) is what killed
    # the first attempt.
    acc = {}
    with Pool(3) as pool:
        for done, part in enumerate(
            pool.imap_unordered(pass2, list(picks_by_batch.items())), 1
        ):
            for k, v in part.items():
                if k in acc:
                    acc[k] += v
                else:
                    acc[k] = v
            del part
            if done % 10 == 0:
                print(f"  summed {done}/{len(picks_by_batch)} libraries", flush=True)
    print(f"  {len(acc)} (rep x group x lineage) pseudobulks", flush=True)

    keys = sorted({k[1:] for k in acc})
    counts = np.column_stack([
        np.median([acc[(r,) + k] for r in range(N_REPS) if (r,) + k in acc], axis=0)
        for k in keys
    ])
    ps = PseudobulkSet(
        counts=counts, genes=np.asarray(genes),
        patient=np.array([k[0] for k in keys]),
        condition=np.array([k[1] for k in keys]),
        lineage=np.array([k[2] for k in keys]),
        n_cells=np.array([1] * len(keys)),
    ).complete_pairs()

    # ---- gating genes are barred from testing (Phase 2.1) ----------
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    keep = detection_filter(ps)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    print(f"pseudobulk: {ps.counts.shape[0]} genes x {ps.counts.shape[1]} samples "
          f"({len(set(ps.patient))} patients)", flush=True)

    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    rep = se_balance_report(de)
    print("S3:", rep["verdict"], flush=True)
    pd.DataFrame([rep["median_SE"]]).to_csv(os.path.join(OUT, "power_check.csv"), index=False)

    print("Phase 3: null calibration", flush=True)
    null = calibrate(ps, case_label="Tumor", n_perm=200, rng=np.random.default_rng(1))
    cross_lineage_noise_summary(null).to_csv(
        os.path.join(OUT, "null_distribution_summary.csv"), index=False)

    for delta in (0.25, 0.5, 1.0):
        for ek in ("p50", "p95"):
            res = classify_table(de, delta=delta, null_stats=null, equiv_key=ek)
            res.to_csv(os.path.join(OUT, f"quadrants_GSE154826_d{delta}_{ek}.csv"),
                       index=False)
            if abs(delta - 0.5) < 1e-9:
                chk = control_check(res, excluded=excl)
                chk.to_csv(os.path.join(OUT, f"control_check_{ek}.csv"), index=False)
                print(f"delta={delta} {ek}: ",
                      res.quadrant.value_counts().to_dict(), flush=True)
                print("  S4:", evaluate_stop_rule_s4(chk)["verdict"], flush=True)
            else:
                print(f"delta={delta} {ek}: ",
                      res.quadrant.value_counts().to_dict(), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
