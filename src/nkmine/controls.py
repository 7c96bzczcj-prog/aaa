"""Built-in positive and negative controls (Protocol Phase 5.3).

Run immediately after classification and before looking at any new
result.  Stop rule S4 fires on failure: without these there is no way to
tell an empty Q4 list from a pipeline with no power.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# TCR-proximal / exhaustion programme.  Driven by TCR signalling, which
# NK cells do not have, so T cells should move and NK should not.  This
# is the positive control for Q4 itself.
Q4_CONTROLS = ["TOX", "CD3D", "CD3E", "LAT", "ZAP70", "ITK", "CD28",
               "TNFRSF9", "PDCD1", "CTLA4", "LAG3", "TIGIT", "NR4A1", "NR4A2"]

# Dissociation stress / immediate-early genes.  Shared across every cell
# type (Marsh et al. 2022), so they must land in Q3.
Q3_CONTROLS = ["FOS", "JUN", "JUNB", "JUND", "EGR1", "DUSP1", "ZFP36",
               "HSPA1A", "HSPA1B", "IER2", "KLF6", "SOCS3"]

# Cytotoxic effector programme shared by NK and CD8 T, absent from B and
# myeloid -- the Q2 expectation.
Q2_CONTROLS = ["GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "GNLY", "KLRG1", "CTSW"]

# Housekeeping: should not land in any quadrant.
NEG_CONTROLS = ["ACTB", "GAPDH", "B2M", "TPT1", "EEF1A1", "RPL13A"]

EXPECTED = {"Q4": Q4_CONTROLS, "Q3": Q3_CONTROLS,
            "Q2": Q2_CONTROLS, "unclassified": NEG_CONTROLS}


def control_check(quadrants: pd.DataFrame, excluded: set | None = None) -> pd.DataFrame:
    """Compare control genes against their expected quadrant.

    Genes used for lineage assignment are excluded from testing by
    Phase 2.1, so some controls are simply unavailable.  This is a real
    tension in the protocol rather than an implementation detail: the
    Q2 positive control list (GZMB, PRF1, GNLY) overlaps the NK
    lineage-defining panel, and the Q4 control list overlaps the T-cell
    panel.  Unavailable controls are reported as `excluded_by_gating`,
    not silently dropped, because a control list that quietly shrinks to
    nothing is indistinguishable from one that passes.
    """
    excluded = excluded or set()
    obs = dict(zip(quadrants.gene, quadrants.quadrant))
    rows = []
    for expected_q, genes in EXPECTED.items():
        for g in genes:
            if g in excluded:
                status = "excluded_by_gating"
                got = None
            elif g not in obs:
                status = "absent_from_panel"
                got = None
            else:
                got = obs[g]
                status = "pass" if got == expected_q else "FAIL"
            rows.append({"gene": g, "expected": expected_q,
                         "observed": got, "status": status})
    return pd.DataFrame(rows)


def evaluate_stop_rule_s4(check: pd.DataFrame) -> dict:
    """Stop rule S4.

    Two specific failures are fatal rather than merely disappointing:
    an immediate-early gene landing in Q1 or Q4 means lineage-specific
    signal is being read out of a universal artefact, and TOX failing to
    reach Q4 means the pipeline cannot detect the very pattern it exists
    to find.
    """
    testable = check[check.status.isin(["pass", "FAIL"])]
    ieg = testable[testable.expected == "Q3"]
    ieg_misfiled = ieg[ieg.observed.isin(["Q1", "Q4"])]

    tox = testable[testable.gene == "TOX"]
    tox_ok = bool(len(tox) and (tox.observed == "Q4").all())

    q4_pos = testable[testable.expected == "Q4"]
    q4_rate = float((q4_pos.observed == "Q4").mean()) if len(q4_pos) else np.nan

    problems = []
    if len(ieg_misfiled):
        problems.append(
            f"immediate-early genes in Q1/Q4: {sorted(ieg_misfiled.gene)}"
        )
    if len(tox) and not tox_ok:
        problems.append(f"TOX did not reach Q4 (observed {tox.observed.iloc[0]})")
    if len(q4_pos) and q4_rate == 0.0:
        problems.append("no TCR-proximal control reached Q4 -- no demonstrated Q4 power")

    return {
        "n_testable": int(len(testable)),
        "n_excluded_by_gating": int((check.status == "excluded_by_gating").sum()),
        "n_absent": int((check.status == "absent_from_panel").sum()),
        "pass_rate": float((testable.status == "pass").mean()) if len(testable) else np.nan,
        "q4_positive_control_rate": q4_rate,
        "problems": problems,
        "stop": bool(problems),
        "verdict": ("PASS - controls behave as expected"
                    if not problems else
                    "STOP (S4) - fix the pipeline before reading new results: "
                    + "; ".join(problems)),
    }
