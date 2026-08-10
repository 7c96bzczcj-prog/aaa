"""Four-quadrant classification with equivalence testing (Protocol Phase 5).

The load-bearing idea of the whole protocol: Q4 ("other lineages move,
NK does not") and Q2 both rest on *accepting* a null on one side.  A
non-significant p-value cannot do that job -- NK is a rare lineage, so
p > 0.05 is exactly what low power produces.  Every "did not change"
call here is therefore a TOST equivalence result (Schuirmann 1987):
the 90% CI must sit entirely inside (-delta, +delta).

The two sides of the decision are deliberately asymmetric.  Claiming
"changed" needs the 95% CI to exclude zero; claiming "unchanged" needs
the 90% CI inside the equivalence margin *and* an effect below the
lineage-and-gene-specific noise median from Phase 3.  Asserting absence
is the stronger claim and is held to the stronger standard.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

# The five lineages the protocol tracks in one droplet pool.
LINEAGES = ("NK", "CD8T", "CD4T", "B", "Myeloid")

# Lineages that must move for a Q4 call (NK's "control group").
# CD4T is excluded here on purpose: the protocol's Q4 rule names
# CD8T / B / Myeloid, so that a Q4 gene is one that moves in both
# lymphoid and myeloid compartments while NK sits still.
Q4_WITNESSES = ("CD8T", "B", "Myeloid")


@dataclass
class LineageStat:
    """One gene, one lineage: the (effect, uncertainty) vector.

    Deliberately no ratios anywhere -- Phase 4 of the protocol bans
    them, because a ratio imports the denominator's replicate noise
    into the decision.
    """

    log2fc: float
    se: float
    df: float
    ci95: tuple[float, float]
    ci90: tuple[float, float]
    null_change: float = np.inf  # Phase 3 ceiling an effect must EXCEED to be "changed"
    null_equiv: float = 0.0      # Phase 3 bound an effect must sit UNDER to be "equivalent"
    saturated: bool = False    # Phase 2.6 floor/ceiling flag
    detected: bool = True      # Phase 2.5 detection floor


def tost(log2fc: float, se: float, df: float, delta: float) -> tuple[float, bool]:
    """Two one-sided tests for equivalence to zero within +/- delta.

    Returns (p_TOST, equivalent_at_0.05).  The boolean is identical to
    "the 90% CI lies inside (-delta, delta)", which is how the protocol
    words the rule; we compute both and assert they agree.
    """
    if not np.isfinite(se) or se <= 0:
        return 1.0, False
    if np.isinf(df):
        sf = stats.norm.sf
    else:
        def sf(x):
            return stats.t.sf(x, df)
    p_lower = sf((log2fc + delta) / se)          # H0: effect <= -delta
    p_upper = sf((delta - log2fc) / se)          # H0: effect >= +delta
    p = float(max(p_lower, p_upper))
    return p, p < 0.05


def is_changed(s: LineageStat, delta: float) -> bool:
    """Protocol 5.1: CI95 excludes 0, effect exceeds delta, and effect
    exceeds this lineage's own 99th-percentile null for this gene."""
    if not s.detected:
        return False
    lo, hi = s.ci95
    excludes_zero = (lo > 0) or (hi < 0)
    return bool(
        excludes_zero
        and abs(s.log2fc) >= delta
        and abs(s.log2fc) > s.null_change
    )


def is_equivalent(s: LineageStat, delta: float) -> bool:
    """Protocol 5.1: TOST-equivalent AND below the gene's null median.

    `saturated` genes (Phase 2.6 ceiling/floor) can never be called
    equivalent, because "no room to move" is not the same as "did not
    move" -- that conflation is one of the errors this protocol was
    written to avoid.
    """
    if s.saturated or not s.detected:
        return False
    lo, hi = s.ci90
    ci_inside = (lo > -delta) and (hi < delta)
    _, tost_ok = tost(s.log2fc, s.se, s.df, delta)
    return bool(ci_inside and tost_ok and abs(s.log2fc) < s.null_equiv)


def _same_sign(values) -> bool:
    vals = [v for v in values if v != 0]
    if not vals:
        return False
    return all(v > 0 for v in vals) or all(v < 0 for v in vals)


@dataclass
class QuadrantCall:
    gene: str
    quadrant: str
    changed: dict = field(default_factory=dict)
    equivalent: dict = field(default_factory=dict)
    reason: str = ""


def classify_gene(gene: str, stats_by_lineage: dict, delta: float,
                  witnesses=Q4_WITNESSES) -> QuadrantCall:
    """Assign one gene to Q1 / Q2 / Q3 / Q4 / unclassified.

    Most genes land in `unclassified`; the protocol expects that and it
    is not a failure mode.  Order matters: Q3 (everything moves) is
    tested before the specific quadrants so that a global shift cannot
    masquerade as a lineage-specific one.
    """
    missing = [l for l in LINEAGES if l not in stats_by_lineage]
    if missing:
        raise KeyError(f"{gene}: missing lineages {missing}")

    changed = {l: is_changed(stats_by_lineage[l], delta) for l in LINEAGES}
    equiv = {l: is_equivalent(stats_by_lineage[l], delta) for l in LINEAGES}
    fc = {l: stats_by_lineage[l].log2fc for l in LINEAGES}

    def call(q, why):
        return QuadrantCall(gene, q, changed, equiv, why)

    changed_lineages = [l for l in LINEAGES if changed[l]]

    # --- Q3: universal / environmental / technical -----------------
    # The protocol's summary table defines Q3 as "ALL lineages move in the
    # same direction", but its operational rule says ">= 4 lineages moved".
    # With five lineages those disagree exactly on the case of interest: a
    # true Q4 gene (NK still, the other four moving) satisfies ">= 4" and
    # would be absorbed into Q3, which on the planted-truth simulation
    # swallowed 25/25 Q4 genes.  A gene whose NK effect is demonstrably
    # equivalent to zero is not a universal shift by definition, so Q3 is
    # required not to fire when NK is TOST-equivalent.  See docs/DEVIATIONS.md.
    if (
        len(changed_lineages) >= 4
        and _same_sign([fc[l] for l in changed_lineages])
        and not equiv["NK"]
    ):
        return call("Q3", f"{len(changed_lineages)} lineages moved concordantly")

    # --- Q4: NK resistance (the primary target) --------------------
    if equiv["NK"]:
        moved = [l for l in witnesses if changed[l]]
        if len(moved) >= 2 and _same_sign([fc[l] for l in moved]):
            return call(
                "Q4",
                f"NK equivalent to zero; {'+'.join(moved)} moved concordantly",
            )

    # --- Q1: NK-specific -------------------------------------------
    if changed["NK"] and all(equiv[l] for l in LINEAGES if l != "NK"):
        return call("Q1", "only NK moved; all other lineages equivalent")

    # --- Q2: shared cytotoxic-lymphocyte program -------------------
    if (
        changed["NK"]
        and changed["CD8T"]
        and _same_sign([fc["NK"], fc["CD8T"]])
        and equiv["B"]
        and equiv["Myeloid"]
    ):
        return call("Q2", "NK and CD8T moved concordantly; B and Myeloid equivalent")

    return call("unclassified", "no quadrant rule satisfied")


def classify_table(de_by_lineage: dict, delta: float, null_stats: dict | None = None,
                   saturated: dict | None = None, detected: dict | None = None,
                   change_key: str = "p99", equiv_key: str = "p50",
                   witnesses=Q4_WITNESSES):
    """Vectorised driver over a gene x lineage panel of DEResults.

    `de_by_lineage` maps lineage -> DEResult sharing a common gene
    order.  `null_stats` maps lineage -> dict of Phase 3 percentile
    arrays.  Returns a pandas DataFrame, one row per gene.

    `equiv_key` selects the Phase 3 percentile an effect must fall under
    to count as unchanged.  The protocol specifies "p50", but note what
    that costs: a gene with a genuinely zero effect sits below its own
    null *median* only 50% of the time, so demanding it simultaneously
    across the four non-NK lineages retains ~6% of true Q1 genes, and it
    halves Q4 sensitivity on the NK side.  TOST already controls the
    false-equivalence rate at 5%; this extra screen buys robustness to
    non-parametric noise the model's SE misses, at a steep and
    non-obvious power cost.  `scripts/sensitivity_equiv_bound.py`
    quantifies the trade-off; "p95" is the calibrated alternative.
    """
    import pandas as pd

    lineages = list(LINEAGES)
    ref = de_by_lineage[lineages[0]]
    genes = list(ref.genes)
    for l in lineages:
        if list(de_by_lineage[l].genes) != genes:
            raise ValueError(f"gene order mismatch in lineage {l}")

    if null_stats is None:
        # Phase 3 is mandatory in the protocol; running without it is only
        # meaningful for unit tests, so make the fallback explicit and loud
        # rather than silently permissive on one side and strict on the other.
        warnings.warn(
            "classify_table called without Phase 3 null calibration; "
            "the empirical noise-floor criteria are disabled. This is not a "
            "protocol-valid run.",
            RuntimeWarning,
            stacklevel=2,
        )

    n_genes = len(genes)
    permissive_p99 = np.zeros(n_genes)      # no noise floor to clear
    permissive_p50 = np.full(n_genes, np.inf)  # no noise ceiling to sit under

    rows = []
    for i, g in enumerate(genes):
        sbl = {}
        for l in lineages:
            d = de_by_lineage[l]
            ns = (null_stats or {}).get(l, {})
            sbl[l] = LineageStat(
                log2fc=float(d.log2fc[i]),
                se=float(d.se[i]),
                df=float(d.df_total),
                ci95=(float(d.ci_lo_95[i]), float(d.ci_hi_95[i])),
                ci90=(float(d.ci_lo_90[i]), float(d.ci_hi_90[i])),
                null_change=float(ns.get(change_key, permissive_p99)[i]),
                null_equiv=float(ns.get(equiv_key, permissive_p50)[i]),
                saturated=bool((saturated or {}).get(l, np.zeros(n_genes, bool))[i]),
                detected=bool((detected or {}).get(l, np.ones(n_genes, bool))[i]),
            )
        c = classify_gene(g, sbl, delta, witnesses=witnesses)
        row = {"gene": g, "quadrant": c.quadrant, "reason": c.reason, "delta": delta}
        for l in lineages:
            row[f"{l}_log2FC"] = sbl[l].log2fc
            row[f"{l}_SE"] = sbl[l].se
            row[f"{l}_CI95_lo"] = sbl[l].ci95[0]
            row[f"{l}_CI95_hi"] = sbl[l].ci95[1]
            row[f"{l}_changed"] = c.changed[l]
            row[f"{l}_equivalent"] = c.equivalent[l]
            row[f"{l}_TOST_p"] = tost(sbl[l].log2fc, sbl[l].se, sbl[l].df, delta)[0]
        rows.append(row)

    return pd.DataFrame(rows)
