"""How many Q4 calls does this pipeline invent from noise?

Stop rule S4 currently blocks the Q4 list because the protocol's
positive control (TOX) turns out to be floor-saturated in NK and cannot
serve (D10). A positive control answers "can the pipeline see a real Q4
gene?". This script answers the complementary question, which is
answerable with the data in hand: "how many Q4 genes does the pipeline
report when there is nothing to find?"

The null preserves everything except the thing being tested. Condition
labels are flipped *within* each individual, so pairing, library sizes,
lineage composition, gene-level abundance and the ambient structure are
all untouched; only the tumour/normal assignment is randomised. Any Q4
call surviving that is manufactured.

This does not substitute for a positive control -- a pipeline that
reports nothing under the null and nothing under signal would look
identical here. It bounds the false-positive side only.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nkmine.de import paired_de  # noqa: E402
from nkmine.null_calibration import calibrate  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
from nkmine.quadrant import LINEAGES, classify_table  # noqa: E402
from s4_diagnosis import nk_saturation  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
N_PERM = 12
DELTA = 0.5


def classify(ps: PseudobulkSet, sat: np.ndarray, n_perm_null: int = 60) -> dict:
    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    null = calibrate(ps, case_label="Tumor", n_perm=n_perm_null,
                     rng=np.random.default_rng(1))
    res = classify_table(
        de, delta=DELTA, null_stats=null, equiv_key="p95",
        saturated={l: (sat if l == "NK" else np.zeros(len(sat), bool))
                   for l in LINEAGES},
    )
    return res.quadrant.value_counts().to_dict()


def main():
    d = np.load(os.path.join(OUT, "pseudobulk_raw.npz"), allow_pickle=True)
    ps = PseudobulkSet(d["counts"], d["genes"], d["patient"], d["condition"],
                       d["lineage"], np.ones(len(d["patient"])))
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    keep = detection_filter(ps, require_nk=True, min_lineages=2)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                       ps.condition, ps.lineage, ps.n_cells)
    sat = nk_saturation(ps.genes)

    obs = classify(ps, sat)
    print(f"OBSERVED  : {obs}", flush=True)

    rng = np.random.default_rng(7)
    rows = []
    for k in range(N_PERM):
        cond = ps.condition.copy()
        # flip the condition label within each individual, for ALL lineages
        # together, so cross-lineage structure is preserved
        for pat in np.unique(ps.patient):
            if rng.random() < 0.5:
                m = ps.patient == pat
                cond[m] = np.where(cond[m] == "Tumor", "Normal", "Tumor")
        perm = PseudobulkSet(ps.counts, ps.genes, ps.patient, cond,
                             ps.lineage, ps.n_cells)
        c = classify(perm, sat)
        rows.append({"perm": k, **{q: c.get(q, 0) for q in
                                   ("Q4", "Q4_attenuated", "Q3", "Q1",
                                    "Q3_NK_indeterminate")}})
        print(f"  perm {k+1}/{N_PERM}: Q4={rows[-1]['Q4']} Q3={rows[-1]['Q3']} "
              f"Q1={rows[-1]['Q1']}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "q4_null_rate.csv"), index=False)
    print("\n=== Q4 calls under a within-individual label flip ===", flush=True)
    for q in ("Q4", "Q3", "Q1"):
        o = obs.get(q, 0)
        n = df[q]
        print(f"  {q:3s} observed={o:4d}   null mean={n.mean():6.1f} "
              f"[{n.min()}, {n.max()}]   ratio={o/max(n.mean(),0.5):.1f}x", flush=True)
    print(f"\n  permutations reaching the observed Q4 count: "
          f"{int((df.Q4 >= obs.get('Q4', 0)).sum())}/{N_PERM}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
