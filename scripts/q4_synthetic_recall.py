"""Semi-synthetic Q4 truth set: does this design have Q4 sensitivity at all?

The protocol's Q4 positive control was self-contradictory. It asked for a
gene guaranteed to land in Q4, and `TOX` was chosen *because* NK has no
TCR -- which is the definition of mundane explanation B3, the opposite of
what Q4 claims. The correct property is "NK **could** have changed
(expresses the gene, carries the receptor, the signal arrives) **and did
not**". And Phase 0 established that no such gene is described in the
literature, so no known instance exists to borrow. A biological positive
control is therefore not merely hard to pick -- it is unavailable in
principle.

So construct one instead.

    take real genes where the witness lineages demonstrably move
      AND NK also moves
    remove NK's condition effect, leaving NK's variance structure,
      cell counts, depth and dispersion untouched
    -> genes that are Q4 by construction

Recall on that set is the number that converts "Q4 is untested" into "Q4
is tested". That distinction is the whole point: S2 and S5 are power
failures and S4 is unknown sensitivity, none of which is evidence that
Q4 is rare or absent.

HOW THE NK EFFECT IS REMOVED, AND WHY NOT BY PERMUTATION
--------------------------------------------------------
The obvious construction -- permute the tumour/normal label within each
individual in the NK pseudobulk -- is wrong here, and in a direction that
would understate recall. For a gene whose true NK effect is d, flipping
the sign of each within-pair difference leaves realised differences of
+/-d. The mean goes to zero as intended, but d^2 is absorbed into the
residual variance, so the standard error inflates and TOST fails. A
genuine Q4 gene has no such inflation, so permutation would measure the
pipeline on a harder problem than the real one.

Instead the fitted NK effect is divided out symmetrically: tumour samples
are scaled by 2^(-(1-lambda)*beta/2) and normal samples by
2^(+(1-lambda)*beta/2). The mean effect becomes lambda*beta, the residual
variance is preserved, and the gene's average expression level is
unchanged, so library sizes barely move.

Permutation is still run at lambda=0 as a secondary, and both are
reported, because the two disagree in an informative way.

DOSE CURVE
----------
lambda = 0 / 0.25 / 0.5 / 0.75 shrinks rather than abolishes the NK
effect, giving the detection floor for Q4 expressed as an effect size:
the largest residual NK effect the pipeline still calls Q4.

PRE-REGISTERED DECISION RULE (fixed before the run, all three stop)
-------------------------------------------------------------------
  recall >= 50%   design has Q4 power -> the 46 candidates' zero
                  replication is a TRUE NEGATIVE, i.e. Q4 is rare or
                  absent in this NSCLC tumour/normal contrast
  20-50%          borderline; record the number
  recall < 20%    27 individuals cannot test Q4 -> record as UNTESTED,
                  never as negative
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
DELTA = 0.5
WITNESSES = ("CD8T", "B", "Myeloid")
LAMBDAS = (0.0, 0.25, 0.5, 0.75)
N_PERM_NULL = 100


def load() -> PseudobulkSet:
    d = np.load(os.path.join(OUT, "pseudobulk_raw.npz"), allow_pickle=True)
    ps = PseudobulkSet(d["counts"], d["genes"], d["patient"], d["condition"],
                       d["lineage"], np.ones(len(d["patient"])))
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    keep = detection_filter(ps, require_nk=True, min_lineages=2)
    return PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                         ps.condition, ps.lineage, ps.n_cells)


def fit_all(ps: PseudobulkSet) -> dict:
    out = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        out[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    return out


def classify(ps: PseudobulkSet, sat: np.ndarray) -> pd.DataFrame:
    de = fit_all(ps)
    null = calibrate(ps, case_label="Tumor", n_perm=N_PERM_NULL,
                     rng=np.random.default_rng(1))
    return classify_table(
        de, delta=DELTA, null_stats=null, equiv_key="p95",
        saturated={l: (sat if l == "NK" else np.zeros(len(sat), bool))
                   for l in LINEAGES},
    )


def pick_source(res: pd.DataFrame) -> np.ndarray:
    """Genes where the witnesses move AND NK moves.

    NK must move, otherwise removing its effect is not an intervention
    and the gene would have been Q4 already.
    """
    moved = np.sum([res[f"{l}_changed"].values for l in WITNESSES], axis=0) >= 2
    signs = [np.sign(res[f"{l}_log2FC"].values) * res[f"{l}_changed"].values
             for l in WITNESSES]
    concordant = (np.abs(np.sum(signs, axis=0)) ==
                  np.sum([res[f"{l}_changed"].values for l in WITNESSES], axis=0))
    return moved & concordant & res["NK_changed"].values


def rescale_nk(ps: PseudobulkSet, src: np.ndarray, beta: np.ndarray,
               lam: float) -> PseudobulkSet:
    """Shrink the NK condition effect to `lam` x its fitted value."""
    counts = ps.counts.copy().astype(float)
    nk = ps.lineage == "NK"
    tum = nk & (ps.condition == "Tumor")
    nor = nk & (ps.condition == "Normal")
    half = (1.0 - lam) * beta / 2.0
    counts[np.ix_(src, tum)] *= (2.0 ** (-half[src]))[:, None]
    counts[np.ix_(src, nor)] *= (2.0 ** (+half[src]))[:, None]
    return PseudobulkSet(counts, ps.genes, ps.patient, ps.condition,
                         ps.lineage, ps.n_cells)


def permute_nk(ps: PseudobulkSet, src: np.ndarray, seed: int = 3) -> PseudobulkSet:
    """Secondary construction: flip the condition label within individual,
    NK only. Retained for comparison; see the module docstring for why it
    is not the primary."""
    rng = np.random.default_rng(seed)
    counts = ps.counts.copy().astype(float)
    nk = np.flatnonzero(ps.lineage == "NK")
    for pat in np.unique(ps.patient):
        idx = [i for i in nk if ps.patient[i] == pat]
        if len(idx) == 2 and rng.random() < 0.5:
            a, b = idx
            tmp = counts[np.ix_(src, [a])].copy()
            counts[np.ix_(src, [a])] = counts[np.ix_(src, [b])]
            counts[np.ix_(src, [b])] = tmp
    return PseudobulkSet(counts, ps.genes, ps.patient, ps.condition,
                         ps.lineage, ps.n_cells)


def report(res: pd.DataFrame, src: np.ndarray, label: str, rows: list):
    q = res.quadrant.values
    hit_strict = (q[src] == "Q4").mean()
    hit_any = np.isin(q[src], ["Q4", "Q4_attenuated"]).mean()
    fp = (q[~src] == "Q4").sum()
    rows.append({"construction": label, "n_source": int(src.sum()),
                 "recall_Q4": float(hit_strict),
                 "recall_Q4_or_attenuated": float(hit_any),
                 "false_Q4_outside_source": int(fp),
                 "source_Q3": int((q[src] == "Q3").sum()),
                 "source_unclassified": int((q[src] == "unclassified").sum())})
    print(f"  {label:26s} recall(Q4)={hit_strict:6.1%}  "
          f"(+attenuated)={hit_any:6.1%}  false-Q4 outside={fp}", flush=True)


def main():
    ps = load()
    sat = nk_saturation(ps.genes)
    print(f"panel: {ps.counts.shape[0]} genes x {ps.counts.shape[1]} samples "
          f"({len(set(ps.patient))} individuals)", flush=True)

    base = classify(ps, sat)
    base.to_csv(os.path.join(OUT, "synth_baseline_quadrants.csv"), index=False)
    src = pick_source(base)
    beta = base["NK_log2FC"].values
    print(f"source genes (witnesses>=2 moved concordantly AND NK moved): "
          f"{int(src.sum())}", flush=True)
    print(f"  median |NK log2FC| in source: "
          f"{np.median(np.abs(beta[src])):.3f}", flush=True)
    if src.sum() < 30:
        print("STOP: source set too small to estimate recall", flush=True)
        return

    rows = []
    print("\n=== dose curve: NK effect shrunk to lambda x its fitted value ===",
          flush=True)
    for lam in LAMBDAS:
        mod = rescale_nk(ps, src, beta, lam)
        # library-size perturbation introduced by the construction
        d_lib = np.abs(mod.counts.sum(axis=0) / ps.counts.sum(axis=0) - 1).max()
        res = classify(mod, sat)
        report(res, src, f"rescale lambda={lam}", rows)
        rows[-1]["lambda"] = lam
        rows[-1]["max_libsize_shift"] = float(d_lib)
        if lam == 0.0:
            res.to_csv(os.path.join(OUT, "synth_lambda0_quadrants.csv"),
                       index=False)

    print("\n=== secondary construction (label permutation, NK only) ===",
          flush=True)
    res = classify(permute_nk(ps, src), sat)
    report(res, src, "permute (lambda=0)", rows)
    rows[-1]["lambda"] = 0.0

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "q4_synthetic_recall.csv"), index=False)

    recall = float(df[(df.construction == "rescale lambda=0.0")]
                   .recall_Q4.iloc[0])
    print("\n" + "=" * 66, flush=True)
    print(f"PRIMARY RESULT: recall at lambda=0 = {recall:.1%} "
          f"(n={int(src.sum())} known-true Q4 genes)", flush=True)
    if recall >= 0.50:
        verdict = ("DESIGN HAS Q4 POWER -> the 46 candidates' zero replication "
                   "is a TRUE NEGATIVE: Q4 is rare or absent in this contrast")
    elif recall >= 0.20:
        verdict = "BORDERLINE -- record the number, claim nothing stronger"
    else:
        verdict = ("27 individuals CANNOT TEST Q4 -> record as UNTESTED, "
                   "not as a negative result")
    print(f"VERDICT: {verdict}", flush=True)
    print("=" * 66, flush=True)
    with open(os.path.join(OUT, "q4_synthetic_recall_verdict.txt"), "w") as fh:
        fh.write(f"recall_lambda0\t{recall:.4f}\nn_source\t{int(src.sum())}\n"
                 f"verdict\t{verdict}\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
