"""Artefact filters for Q4/Q2 candidates (Protocol Phase 6).

Three independent confounders, in the order the protocol ranks them:
ambient RNA (6.1), composition-versus-state (6.2), and depth/selection
effects (6.3).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# 6.2  Composition vs state
# --------------------------------------------------------------------------
def kitagawa_decompose(
    prop_ref: np.ndarray,
    prop_alt: np.ndarray,
    mean_ref: np.ndarray,
    mean_alt: np.ndarray,
    interaction_flag: float = 0.25,
) -> dict:
    """Three-term Kitagawa/Oaxaca decomposition of a lineage-level change.

    A lineage's mean expression is `sum_k p_k * m_k` over sub-clusters k.
    The change between conditions splits into

        within      = sum_k p_ref_k  * (m_alt_k - m_ref_k)
        between     = sum_k (p_alt_k - p_ref_k) * m_ref_k
        interaction = sum_k (p_alt_k - p_ref_k) * (m_alt_k - m_ref_k)

    and the three sum exactly to the total.  The protocol insists the
    interaction term be reported: with only `within` and `between` the
    percentages do not add up, and the split looks far more decisive
    than it is.

    When `|interaction| / |total|` exceeds `interaction_flag`, the
    attribution between "cells changed" and "cell mixture changed"
    depends on which condition is used as the reference weight, so the
    result is returned as `not_identified` rather than as a number that
    invites over-reading.
    """
    prop_ref = np.asarray(prop_ref, float)
    prop_alt = np.asarray(prop_alt, float)
    mean_ref = np.asarray(mean_ref, float)
    mean_alt = np.asarray(mean_alt, float)
    if not (len(prop_ref) == len(prop_alt) == len(mean_ref) == len(mean_alt)):
        raise ValueError("all four inputs must have one entry per sub-cluster")
    for p, nm in ((prop_ref, "prop_ref"), (prop_alt, "prop_alt")):
        if not np.isclose(p.sum(), 1.0, atol=1e-6):
            raise ValueError(f"{nm} must sum to 1 (got {p.sum():.6f})")

    dp = prop_alt - prop_ref
    dm = mean_alt - mean_ref
    within = float(np.sum(prop_ref * dm))
    between = float(np.sum(dp * mean_ref))
    interaction = float(np.sum(dp * dm))
    total = float(np.sum(prop_alt * mean_alt) - np.sum(prop_ref * mean_ref))

    assert np.isclose(within + between + interaction, total, atol=1e-8), (
        "decomposition must be exact"
    )

    denom = abs(total) if abs(total) > 1e-12 else np.nan
    frac_int = abs(interaction) / denom if denom == denom else np.nan
    identified = bool(frac_int <= interaction_flag) if frac_int == frac_int else False

    return {
        "within": within,
        "between": between,
        "interaction": interaction,
        "total": total,
        "within_pct": 100 * within / denom if denom == denom else np.nan,
        "between_pct": 100 * between / denom if denom == denom else np.nan,
        "interaction_pct": 100 * interaction / denom if denom == denom else np.nan,
        "identified": identified,
        "flag": "" if identified else "not_identified",
    }


def decompose_reference_swap(prop_ref, prop_alt, mean_ref, mean_alt) -> dict:
    """Run the decomposition both ways round.

    The index-number problem is easiest to see directly: swapping which
    condition supplies the weights moves the within/between split.  If
    the two orderings disagree materially, no single attribution is
    defensible.
    """
    a = kitagawa_decompose(prop_ref, prop_alt, mean_ref, mean_alt)
    b = kitagawa_decompose(prop_alt, prop_ref, mean_alt, mean_ref)
    swing = abs(a["within_pct"] - (-b["within_pct"]))
    return {"forward": a, "reverse": b, "within_pct_swing": swing,
            "stable": bool(swing < 20.0)}


# --------------------------------------------------------------------------
# 6.1  Ambient RNA
# --------------------------------------------------------------------------
def ambient_susceptibility(
    expr_by_celltype: pd.DataFrame,
    abundance_ref: pd.Series,
    abundance_alt: pd.Series,
) -> pd.Series:
    """Heuristic ambient risk score, for datasets with no raw droplets.

    A gene is ambient-risky when it is highly expressed by a cell type
    whose abundance changes between conditions: that combination shifts
    the soup, and the shift is added to *every* lineage's profile at
    once.  `expr_by_celltype` is genes x celltypes mean expression.
    """
    shared = [c for c in expr_by_celltype.columns
              if c in abundance_ref.index and c in abundance_alt.index]
    if not shared:
        raise ValueError("no cell types shared between expression and abundance")
    d_abund = (abundance_alt[shared] - abundance_ref[shared]).abs()
    top_expr = expr_by_celltype[shared].max(axis=1)
    top_type = expr_by_celltype[shared].idxmax(axis=1)
    return (top_expr * top_type.map(d_abund)).rename("ambient_susceptibility")


def shared_soup_contrast(shared_pool: pd.DataFrame, separate_pool: pd.DataFrame,
                         key: str = "gene") -> pd.DataFrame:
    """Protocol 6.1's cleanest control, needing no decontamination tool.

    In hashed libraries both conditions share one droplet emulsion, so
    ambient RNA contributes the *same* offset to both and largely
    cancels from the within-library contrast.  In separately-built
    libraries it does not cancel, and it is added to every lineage
    alike, which inflates apparent cross-lineage sharing.

    The prediction is therefore directional: cross-lineage concordance
    should be **higher** in the separately-built set, and the difference
    estimates ambient's contribution.  Both frames need columns
    `gene` and `concordance`.
    """
    m = shared_pool.merge(separate_pool, on=key, suffixes=("_shared", "_separate"))
    m["ambient_contribution"] = m["concordance_separate"] - m["concordance_shared"]
    return m


# --------------------------------------------------------------------------
# 6.3  Depth / selection
# --------------------------------------------------------------------------
def umi_floor_scan(run_fn, floors=(1000, 2000, 4000, 6000)) -> pd.DataFrame:
    """Re-run an effect estimate across UMI floors.

    A minimum-UMI floor applied before downsampling is itself a cell
    filter: it removes shallow cells, which are not a random subset.
    `run_fn(floor)` should return a dict of effect sizes keyed by name.
    """
    rows = []
    for f in floors:
        res = run_fn(f)
        rows.append({"umi_floor": f, **res})
    df = pd.DataFrame(rows)
    num = df.drop(columns=["umi_floor"]).select_dtypes("number")
    if len(num):
        rng = num.max() - num.min()
        df.attrs["max_range_across_floors"] = float(rng.max())
        df.attrs["sensitive_to_floor"] = bool(rng.max() > 0.3)
    return df
