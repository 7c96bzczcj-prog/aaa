"""Cross-dataset replication and mundane-explanation screening
(Protocol Phases 7 and 8).

Phase 7 is a hard gate: a candidate that does not reproduce in a second,
independently-collected cohort does not advance, however artefact-proof
it looked in the first.  Resistance to artefacts and reproducibility are
separate tests, and passing one says nothing about the other.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def replicate(a: pd.DataFrame, b: pd.DataFrame, fold_tolerance: float = 2.0,
              lineage: str = "NK") -> pd.DataFrame:
    """Match quadrant calls across two datasets.

    Replication requires the same quadrant, the same sign, and effect
    sizes within `fold_tolerance` of each other -- the protocol's three
    conditions, all of which must hold.
    """
    key = f"{lineage}_log2FC"
    cols = ["gene", "quadrant", key]
    m = a[cols].merge(b[cols], on="gene", suffixes=("_a", "_b"))

    same_quadrant = m["quadrant_a"] == m["quadrant_b"]
    fa, fb = m[f"{key}_a"], m[f"{key}_b"]
    same_sign = np.sign(fa) == np.sign(fb)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.abs(fa) / np.abs(fb)
    within = (ratio <= fold_tolerance) & (ratio >= 1 / fold_tolerance)

    m["same_quadrant"] = same_quadrant
    m["same_sign"] = same_sign
    m["within_fold"] = within.fillna(False)
    m["replicated"] = same_quadrant & same_sign & m["within_fold"]
    return m


def replication_rate_by_quadrant(rep: pd.DataFrame) -> pd.DataFrame:
    """Replication rate per quadrant, with Q1 and Q3 as the baseline.

    Stop rule S6: if Q4 reproduces markedly worse than Q1 and Q3 on the
    same data and the same criteria, the Q4 list is mostly noise, and
    the comparison to those baselines is the only way to know.
    """
    out = (rep.groupby("quadrant_a")
              .agg(n=("replicated", "size"), n_replicated=("replicated", "sum"))
              .reset_index()
              .rename(columns={"quadrant_a": "quadrant"}))
    out["rate"] = out.n_replicated / out.n
    return out


def evaluate_stop_rules_s5_s6(rates: pd.DataFrame, min_q4: int = 10,
                              margin: float = 0.5) -> dict:
    """Stop rules S5 (too few replicated Q4) and S6 (Q4 worse than Q1/Q3)."""
    def get(q, col):
        r = rates[rates.quadrant == q]
        return float(r[col].iloc[0]) if len(r) else np.nan

    n_q4 = get("Q4", "n_replicated")
    rate_q4 = get("Q4", "rate")
    baseline = np.nanmean([get("Q1", "rate"), get("Q3", "rate")])

    problems = []
    if not np.isnan(n_q4) and n_q4 < min_q4:
        problems.append(f"S5: only {n_q4:.0f} Q4 genes replicated (< {min_q4})")
    if not np.isnan(rate_q4) and not np.isnan(baseline) and rate_q4 < margin * baseline:
        problems.append(
            f"S6: Q4 replication rate {rate_q4:.2f} is below {margin:g}x the "
            f"Q1/Q3 baseline {baseline:.2f} -- Q4 is behaving like noise"
        )
    return {"n_q4_replicated": n_q4, "q4_rate": rate_q4,
            "q1q3_baseline_rate": baseline, "problems": problems,
            "stop": bool(problems),
            "verdict": "PASS" if not problems else "STOP - " + "; ".join(problems)}


# --------------------------------------------------------------------------
# Phase 8
# --------------------------------------------------------------------------
MUNDANE = {
    "B1_not_exposed": "NK excluded from the tumour core; check spatial data or "
                      "vascular/marginal-zone marker proximity",
    "B2_fast_turnover": "NK population continually replaced; check MKI67, SELL, "
                        "S1PR1 higher in NK than other lineages on the tumour side",
    "B3_missing_receptor": "the driving receptor is not expressed by NK; check "
                           "pathway receptor expression in the NK gate",
    "B4_ceiling_floor": "no room to move; re-check detection rate and expression "
                        "distribution (should already be caught by Phase 2.6)",
    "B5_underpowered": "NK SE inflated relative to other lineages (should already "
                       "be caught by Phase 2.4)",
    "B6_ambient": "driven by soup; check ambient susceptibility and the hashed "
                  "shared-soup contrast",
}


def mundane_screen(candidates: pd.DataFrame, evidence: dict) -> pd.DataFrame:
    """Tick off B1-B6 for each Q4 candidate.

    `evidence` maps a B-code to a dict of gene -> bool, True meaning the
    mundane explanation *is* supported for that gene.  A gene survives
    only when every code is explicitly refuted; codes with no evidence
    supplied are recorded as `unassessed`, which blocks survival rather
    than permitting it.  "NK resistance" is the last explanation
    standing, never the first one reached for.
    """
    rows = []
    for g in candidates.gene:
        rec = {"gene": g}
        unresolved, explained = [], []
        for code in MUNDANE:
            ev = evidence.get(code)
            if ev is None or g not in ev:
                rec[code] = "unassessed"
                unresolved.append(code)
            elif ev[g]:
                rec[code] = "EXPLAINS"
                explained.append(code)
            else:
                rec[code] = "refuted"
        rec["n_explained"] = len(explained)
        rec["n_unassessed"] = len(unresolved)
        rec["survives"] = (not explained) and (not unresolved)
        rec["status"] = ("survives" if rec["survives"]
                         else "explained_by:" + ",".join(explained) if explained
                         else "incomplete")
        rows.append(rec)
    return pd.DataFrame(rows)
