"""The same semi-synthetic Q4 recall test, run in the REPLICATION cohort.

Why this is not optional. `q4_synthetic_recall.py` measured recall in the
discovery cohort (GSE154826, n = 27) and got 76%. The pre-registered rule
then reads "design has Q4 power -> the zero replication is a true
negative". That inference silently crosses two cohorts:

  * discovery (n = 27) did NOT find zero -- it found 46 Q4 candidates;
  * the zero is in GSE131907 (n = 8, relaxed 20-cell threshold).

Recall measured where the candidates were found says nothing about
whether the cohort that failed to reproduce them could have reproduced
anything at all. So the same construction is applied to the replication
cohort's own pseudobulk, using genes selected in that cohort.

  recall_kim >= 50%   the replication cohort can detect Q4 -> its zero is
                      a TRUE NEGATIVE and the 46 candidates are refuted
  recall_kim <  20%   the replication cohort cannot test Q4 -> its zero is
                      UNINFORMATIVE, and Q4 is unreplicated rather than
                      refuted
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nkmine.pseudobulk import PseudobulkSet  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402
from q4_synthetic_recall import (  # noqa: E402
    DELTA, WITNESSES, classify, pick_source, rescale_nk,
)

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
LAMBDAS = (0.0, 0.25, 0.5)


def main():
    f = os.path.join(OUT, "pseudobulk_kim.npz")
    if not os.path.exists(f):
        print("missing pseudobulk_kim.npz -- run phase7_replicate_kim.py first")
        return
    d = np.load(f, allow_pickle=True)
    ps = PseudobulkSet(d["counts"], d["genes"], d["patient"], d["condition"],
                       d["lineage"], np.ones(len(d["patient"])))
    print(f"GSE131907 panel: {ps.counts.shape[0]} genes x "
          f"{ps.counts.shape[1]} samples ({len(set(ps.patient))} individuals)",
          flush=True)

    # No saturation flags for this cohort: NK detection rates were computed
    # from GSE154826 cells and do not transfer.  Passing all-False would be
    # silently permissive, so it is stated rather than hidden.
    sat = np.zeros(ps.counts.shape[0], bool)
    print("NOTE: Phase 2.6 saturation flags unavailable for this cohort "
          "(computed from GSE154826 cells); recall here is therefore an "
          "UPPER bound.", flush=True)

    base = classify(ps, sat)
    base.to_csv(os.path.join(OUT, "synth_kim_baseline_quadrants.csv"),
                index=False)
    src = pick_source(base)
    beta = base["NK_log2FC"].values
    print(f"source genes in GSE131907 (witnesses>=2 moved concordantly AND "
          f"NK moved): {int(src.sum())}", flush=True)

    if src.sum() < 10:
        print("\n" + "=" * 66, flush=True)
        print(f"RESULT: only {int(src.sum())} genes in the replication cohort "
              "even have moving witnesses alongside a moving NK.", flush=True)
        print("A synthetic Q4 gene cannot be constructed there in useful "
              "numbers, which is itself the answer: the cohort cannot "
              "present the pipeline with a detectable Q4 case.", flush=True)
        print("VERDICT: replication cohort CANNOT TEST Q4 -> its zero is "
              "UNINFORMATIVE. Q4 is unreplicated, not refuted.", flush=True)
        print("=" * 66, flush=True)
        with open(os.path.join(OUT, "q4_synthetic_recall_kim_verdict.txt"),
                  "w") as fh:
            fh.write(f"n_source\t{int(src.sum())}\nrecall_lambda0\tNA\n"
                     "verdict\tcannot test Q4; zero replication uninformative\n")
        print("DONE", flush=True)
        return

    rows = []
    for lam in LAMBDAS:
        res = classify(rescale_nk(ps, src, beta, lam), sat)
        q = res.quadrant.values
        rec = float((q[src] == "Q4").mean())
        rows.append({"lambda": lam, "n_source": int(src.sum()),
                     "recall_Q4": rec,
                     "recall_Q4_or_attenuated":
                         float(np.isin(q[src], ["Q4", "Q4_attenuated"]).mean())})
        print(f"  lambda={lam}: recall(Q4)={rec:.1%}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "q4_synthetic_recall_kim.csv"), index=False)
    rec0 = float(df[df["lambda"] == 0.0].recall_Q4.iloc[0])
    if rec0 >= 0.50:
        verdict = ("replication cohort CAN detect Q4 -> its zero is a TRUE "
                   "NEGATIVE and the 46 candidates are refuted")
    elif rec0 >= 0.20:
        verdict = "borderline; the zero is weak evidence at best"
    else:
        verdict = ("replication cohort CANNOT test Q4 -> its zero is "
                   "UNINFORMATIVE. Q4 is unreplicated, not refuted")
    print("\n" + "=" * 66, flush=True)
    print(f"RESULT: recall at lambda=0 in GSE131907 = {rec0:.1%} "
          f"(n={int(src.sum())})", flush=True)
    print(f"VERDICT: {verdict}", flush=True)
    print("=" * 66, flush=True)
    with open(os.path.join(OUT, "q4_synthetic_recall_kim_verdict.txt"),
              "w") as fh:
        fh.write(f"recall_lambda0\t{rec0:.4f}\nn_source\t{int(src.sum())}\n"
                 f"verdict\t{verdict}\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
