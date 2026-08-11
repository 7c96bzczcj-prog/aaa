"""Standard error under each rho estimator, against the matched baseline.

This is the read-out for the pre-registered question in
`rho_variants.py`: pooling rho across conditions removes the per-sample
sampling noise, so if mechanism (a) -- "rho is a noisy per-sample
quantity and subtraction injects that noise into every gene" -- is the
binding constraint, the pooled arm's SE should come back to baseline.

The gene panel is fixed on the uncorrected data and shared by all four
arms, so the arms differ only in the counts, never in which genes are
tested.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.de import paired_de  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")

ARMS = [
    ("baseline", "pseudobulk_uncorrected_matched.npz"),
    ("rho_per_sample", "pseudobulk_per_sample.npz"),
    ("rho_pooled", "pseudobulk_pooled_cond.npz"),
    ("rho_wide_markers", "pseudobulk_wide_markers.npz"),
]
IG = ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"]
CYTO_IN_B = ["GZMB", "GZMA", "PRF1", "CTSW", "GZMH"]


def load(name: str) -> PseudobulkSet:
    d = np.load(os.path.join(OUT, name), allow_pickle=True)
    return PseudobulkSet(d["counts"].astype(float), d["genes"], d["patient"],
                         d["condition"], d["lineage"],
                         np.ones(len(d["patient"]))).complete_pairs()


def main():
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    base = load(ARMS[0][1])
    ok = ~np.isin(base.genes, list(excl))
    tmp = PseudobulkSet(base.counts[ok], base.genes[ok], base.patient,
                        base.condition, base.lineage, base.n_cells)
    keep = detection_filter(tmp, require_nk=True, min_lineages=2)
    shared = tmp.genes[keep]
    print(f"shared panel fixed on uncorrected counts: {len(shared)} genes",
          flush=True)

    rows = []
    for arm, fname in ARMS:
        ps = load(fname)
        sel = np.isin(ps.genes, shared)
        ps = PseudobulkSet(ps.counts[sel], ps.genes[sel], ps.patient,
                           ps.condition, ps.lineage, ps.n_cells)
        for lin in LINEAGES:
            s = ps.subset_lineage(lin)
            de = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
            ctrl = IG if lin != "B" else CYTO_IN_B
            m = np.isin(s.genes, ctrl)
            rows.append({
                "arm": arm, "lineage": lin,
                "median_SE": float(np.median(de.se)),
                "zero_truth_abs_log2FC":
                    float(np.median(np.abs(de.log2fc[m]))) if m.any() else np.nan,
            })
        print(f"  {arm} done", flush=True)

    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(OUT, "rho_variant_se.csv"), index=False)
    print("\n=== median SE by arm ===", flush=True)
    print(r.pivot(index="arm", columns="lineage", values="median_SE")
          .reindex([a for a, _ in ARMS])[list(LINEAGES)].round(3).to_string(), flush=True)
    print("\n=== known-zero controls, |log2FC| (Ig; granzymes in B) ===", flush=True)
    print(r.pivot(index="arm", columns="lineage", values="zero_truth_abs_log2FC")
          .reindex([a for a, _ in ARMS])[list(LINEAGES)].round(3).to_string(), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
