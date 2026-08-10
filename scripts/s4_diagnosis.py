"""Why stop rule S4 fires, and what it takes to satisfy it.

The first real-data run produced Q4 candidates but failed its own
positive control. This script isolates the two causes and tests each
fix, because "S4 failed" is only useful if it comes with a reason.

Cause 1 (Phase 2.5). The detection floor requires every gene to be
expressed in all five lineages. TCR-proximal genes are T-restricted, so
the entire Q4 positive-control panel is removed before testing.

Cause 2 (Phase 5.2 vs 5.3). The Q4 rule requires >= 2 witnesses drawn
from {CD8T, B, Myeloid}. A TCR-driven gene moves only in T lineages, so
it can supply at most one witness and can never be called Q4 -- the
protocol's stated Q4 rule and its stated Q4 positive control are
mutually unsatisfiable.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.controls import control_check, evaluate_stop_rule_s4, Q4_CONTROLS  # noqa: E402
from nkmine.de import paired_de  # noqa: E402
from nkmine.null_calibration import calibrate  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
import pickle, glob  # noqa: E402
from nkmine.quadrant import LINEAGES, classify_table  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")


def load() -> PseudobulkSet:
    d = np.load(os.path.join(OUT, "pseudobulk_raw.npz"), allow_pickle=True)
    return PseudobulkSet(d["counts"], d["genes"], d["patient"], d["condition"],
                         d["lineage"], np.ones(len(d["patient"])))


def nk_saturation(genes: np.ndarray) -> np.ndarray:
    """Phase 2.6 flags, computed from the cached cell-level gating labels.

    A gene whose NK detection rate is >90% or <5% in BOTH conditions has no
    room to move, so an apparent "no change" there is uninformative and it
    must not be eligible for a Q4 call.
    """
    cache = os.path.join(ROOT, "data", "nk_detection.npz")
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)
        rate = {c: dict(zip(d["genes"], d[c])) for c in ("Tumor", "Normal")}
    else:
        from nkmine.gating import read_10x_tar
        hits = {c: None for c in ("Tumor", "Normal")}
        tot = {c: 0 for c in ("Tumor", "Normal")}
        gnames = None
        for pk in sorted(glob.glob(os.path.join(ROOT, "data", "labels", "*.pkl"))):
            with open(pk, "rb") as fh:
                lab = pickle.load(fh)
            b = lab["batch"]
            tar = os.path.join(ROOT, "data", "raw", f"GSE154826_amp_batch_ID_{b}.tar.gz")
            if not os.path.exists(tar):
                continue
            lib = read_10x_tar(tar)
            x = lib.rna.tocsc()
            gnames = lib.gene_names
            for c in ("Tumor", "Normal"):
                m = np.flatnonzero((lab["lineage"] == "NK") & (lab["condition"] == c))
                if len(m) == 0:
                    continue
                cnt = np.asarray((x[:, m] > 0).sum(axis=1)).ravel()
                hits[c] = cnt if hits[c] is None else hits[c] + cnt
                tot[c] += len(m)
        rate = {c: dict(zip(gnames, hits[c] / max(tot[c], 1))) for c in hits}
        np.savez_compressed(cache, genes=np.asarray(gnames),
                            Tumor=hits["Tumor"] / max(tot["Tumor"], 1),
                            Normal=hits["Normal"] / max(tot["Normal"], 1))
    t = np.array([rate["Tumor"].get(g, 0.0) for g in genes])
    n = np.array([rate["Normal"].get(g, 0.0) for g in genes])
    return ((t > 0.90) & (n > 0.90)) | ((t < 0.05) & (n < 0.05))


def run(ps: PseudobulkSet, keep: np.ndarray, label: str, witnesses, excl: set):
    sub = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                        ps.condition, ps.lineage, ps.n_cells)
    de = {}
    for l in LINEAGES:
        s = sub.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    null = calibrate(sub, case_label="Tumor", n_perm=100,
                     rng=np.random.default_rng(1))
    sat = nk_saturation(sub.genes)
    res = classify_table(de, delta=0.5, null_stats=null, equiv_key="p95",
                         witnesses=witnesses,
                         saturated={l: (sat if l == "NK" else np.zeros(len(sat), bool))
                                    for l in LINEAGES})
    print(f"  NK saturated genes (Phase 2.6): {int(sat.sum())}", flush=True)
    chk = control_check(res, excluded=excl)
    ev = evaluate_stop_rule_s4(chk)
    present = [g for g in Q4_CONTROLS if g in set(res.gene)]
    print(f"\n--- {label} ---", flush=True)
    print(f"  genes tested            : {len(res)}")
    print(f"  Q4 controls in panel    : {len(present)}/{len(Q4_CONTROLS)} {present}")
    print(f"  quadrant counts         : {res.quadrant.value_counts().to_dict()}")
    if present:
        got = res[res.gene.isin(present)][["gene", "quadrant"]]
        print(f"  where they landed       : {dict(zip(got.gene, got.quadrant))}")
    print(f"  S4                      : {ev['verdict'][:110]}")
    return res, chk, ev


def main():
    ps = load()
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)

    strict = detection_filter(ps)
    relaxed = detection_filter(ps, require_nk=True, min_lineages=2)
    print(f"detection filter: strict(all 5 lineages)={strict.sum()} genes, "
          f"relaxed(NK + >=2 lineages)={relaxed.sum()} genes", flush=True)

    run(ps, strict, "A. protocol as written (strict filter, witnesses CD8T/B/Myeloid)",
        ("CD8T", "B", "Myeloid"), excl)
    run(ps, relaxed, "B. fix 1 only (relaxed filter, witnesses CD8T/B/Myeloid)",
        ("CD8T", "B", "Myeloid"), excl)
    res, chk, ev = run(
        ps, relaxed, "C. fix 1 + fix 2 (relaxed filter, witnesses CD8T/CD4T/B/Myeloid)",
        ("CD8T", "CD4T", "B", "Myeloid"), excl)
    res.to_csv(os.path.join(OUT, "quadrants_relaxed_d0.5_p95.csv"), index=False)
    chk.to_csv(os.path.join(OUT, "control_check_relaxed.csv"), index=False)
    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
