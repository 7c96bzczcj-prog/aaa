"""Validation suite for the four-quadrant mining pipeline.

These are not smoke tests.  Protocol Phase 5.3 says that a pipeline
without positive controls cannot be known to have power, and stop rule
S4 says to fix the pipeline before looking at any new result.  The
tests below are the synthetic equivalent of those controls: planted
ground truth, run end to end, with the answer known in advance.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nkmine import simulate as S  # noqa: E402
from nkmine.de import _wls, fit_f_dist, paired_de  # noqa: E402
from nkmine.null_calibration import calibrate  # noqa: E402
from nkmine.pseudobulk import (  # noqa: E402
    PseudobulkSet,
    balanced_pseudobulk,
    detection_filter,
    se_balance_report,
    simple_pseudobulk,
)
from nkmine.quadrant import (  # noqa: E402
    LINEAGES,
    LineageStat,
    classify_gene,
    classify_table,
    is_equivalent,
    tost,
)

warnings.simplefilter("ignore")


# ---------------------------------------------------------------- helpers
def _fit_all(ps):
    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    return de


def _pipeline(sim, balanced=True, min_count=10):
    build = balanced_pseudobulk if balanced else simple_pseudobulk
    kw = dict(rng=np.random.default_rng(0), n_reps=5) if balanced else {}
    ps = build(sim.counts, sim.genes, sim.cell_patient, sim.cell_condition,
               sim.cell_lineage, **kw)
    keep = detection_filter(ps, min_count=min_count)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                       ps.condition, ps.lineage, ps.n_cells)
    return ps, _fit_all(ps)


@pytest.fixture(scope="module")
def planted():
    sim = S.simulate(n_patients=12, n_genes_per_class=25, n_null_genes=150,
                     rng=np.random.default_rng(7))
    ps, de = _pipeline(sim)
    null = calibrate(ps, case_label="Tumor", n_perm=100, rng=np.random.default_rng(3))
    return sim, ps, de, null


# ---------------------------------------------------------------- linear algebra
def test_vectorised_wls_matches_per_gene_reference():
    rng = np.random.default_rng(1)
    X = np.column_stack([np.ones(20), rng.normal(size=(20, 3))])
    y = rng.normal(size=(40, 20))
    w = rng.uniform(0.5, 2.0, size=(40, 20))
    beta, s2, df, su, _ = _wls(y, X, w)
    for g in range(y.shape[0]):
        sw = np.sqrt(w[g])
        b_ref, *_ = np.linalg.lstsq(X * sw[:, None], y[g] * sw, rcond=None)
        assert np.allclose(beta[g], b_ref, atol=1e-10)
    assert df == 16


def test_collinear_design_is_rejected_not_silently_fitted():
    """Admission criterion A2: condition confounded with the blocking
    factor must fail loudly, not return an uninterpretable estimate."""
    counts = np.random.default_rng(0).poisson(50, size=(20, 6)).astype(float)
    genes = np.array([f"g{i}" for i in range(20)])
    patient = np.array(["A", "A", "B", "B", "C", "C"])
    # every patient contributes only one condition -> confounded
    condition = np.array(["Tumor"] * 3 + ["Normal"] * 3)
    patient_confounded = np.array(["A", "B", "C", "D", "E", "F"])
    with pytest.raises(ValueError):
        paired_de(counts, genes, patient_confounded, condition, "Tumor")


def test_fit_f_dist_recovers_known_prior():
    rng = np.random.default_rng(4)
    d0_true, s02_true, df = 8.0, 0.4, 12.0
    s2_prior = s02_true * d0_true / rng.chisquare(d0_true, size=6000)
    sigma2 = s2_prior * rng.chisquare(df, size=6000) / df
    d0, s02 = fit_f_dist(sigma2, df)
    assert 0.5 * d0_true < d0 < 2.0 * d0_true
    assert 0.7 * s02_true < s02 < 1.4 * s02_true


# ---------------------------------------------------------------- TOST
def test_tost_agrees_with_the_90pct_ci_rule():
    """The protocol words equivalence as "90% CI inside (-d, d)"; TOST at
    alpha=0.05 is the same statement. They must never disagree."""
    from scipy import stats as st

    rng = np.random.default_rng(2)
    for _ in range(500):
        lfc = rng.normal(0, 0.5)
        se = abs(rng.normal(0.2, 0.1)) + 1e-3
        df = rng.integers(5, 60)
        d = 0.5
        q90 = st.t.ppf(0.95, df)
        ci_rule = (lfc - q90 * se > -d) and (lfc + q90 * se < d)
        assert tost(lfc, se, df, d)[1] == ci_rule


def test_wide_interval_never_counts_as_equivalent():
    """A large SE means 'we do not know', which must not be recorded as
    'no effect'. This is the failure mode the whole protocol targets."""
    s = LineageStat(log2fc=0.0, se=5.0, df=20, ci95=(-10, 10), ci90=(-8, 8),
                    null_equiv=np.inf)
    assert not is_equivalent(s, delta=0.5)


def test_saturated_gene_cannot_be_equivalent():
    """Phase 2.6: a gene at the detection ceiling has no room to move, so
    'did not move' is not evidence of resistance."""
    s = LineageStat(log2fc=0.01, se=0.02, df=30, ci95=(-0.03, 0.05),
                    ci90=(-0.02, 0.04), null_equiv=np.inf, saturated=True)
    assert not is_equivalent(s, delta=0.5)
    s.saturated = False
    assert is_equivalent(s, delta=0.5)


# ---------------------------------------------------------------- estimation
def test_effect_estimates_are_unbiased(planted):
    sim, _, de, _ = planted
    for l in LINEAGES:
        d = de[l]
        truth = np.array([sim.planted_lfc[g][l] for g in d.genes])
        err = d.log2fc - truth
        assert abs(err.mean()) < 0.10, f"{l} bias {err.mean():.3f}"
        assert np.corrcoef(d.log2fc, truth)[0, 1] > 0.95


def test_power_balancing_equalises_standard_errors(planted):
    _, _, de, _ = planted
    rep = se_balance_report(de)
    assert rep["passes"], rep["verdict"]


# ---------------------------------------------------------------- classification
def test_recovers_planted_quadrants(planted):
    sim, _, de, null = planted
    res = classify_table(de, delta=0.5, null_stats=null, equiv_key="p95")
    truth = np.array([sim.truth[g] for g in res.gene])
    for q, floor in [("Q1", 0.5), ("Q2", 0.6), ("Q3", 0.9), ("Q4", 0.8)]:
        m = truth == q
        rate = (res.quadrant[m] == q).mean()
        assert rate >= floor, f"{q} recovery {rate:.2f} below {floor}"


def test_true_null_genes_are_not_assigned_a_quadrant(planted):
    sim, _, de, null = planted
    res = classify_table(de, delta=0.5, null_stats=null, equiv_key="p95")
    truth = np.array([sim.truth[g] for g in res.gene])
    m = truth == "unclassified"
    assert (res.quadrant[m] == "unclassified").mean() > 0.98


def test_q3_does_not_swallow_q4():
    """Regression test. A Q4 gene has four lineages moving, which also
    satisfies the literal Q3 rule (">= 4 lineages moved"); before the
    precedence fix this misfiled every planted Q4 gene as Q3."""
    def stat(lfc, equivalent):
        se = 0.05
        return LineageStat(log2fc=lfc, se=se, df=40,
                           ci95=(lfc - 2 * se, lfc + 2 * se),
                           ci90=(lfc - 1.7 * se, lfc + 1.7 * se),
                           null_change=0.05 if not equivalent else np.inf,
                           null_equiv=np.inf if equivalent else 0.0)
    sbl = {"NK": stat(0.0, True)}
    for l in ("CD8T", "CD4T", "B", "Myeloid"):
        sbl[l] = stat(1.0, False)
    assert classify_gene("g", sbl, delta=0.5).quadrant == "Q4"


def test_q3_still_fires_when_nk_also_moves():
    def stat(lfc):
        se = 0.05
        return LineageStat(log2fc=lfc, se=se, df=40,
                           ci95=(lfc - 2 * se, lfc + 2 * se),
                           ci90=(lfc - 1.7 * se, lfc + 1.7 * se),
                           null_change=0.05, null_equiv=0.0)
    sbl = {l: stat(1.0) for l in LINEAGES}
    assert classify_gene("g", sbl, delta=0.5).quadrant == "Q3"


def test_opposite_directions_are_not_a_quadrant():
    """Q3/Q4 both require concordant direction; a gene that goes up in
    one lineage and down in another is not a shared programme."""
    def stat(lfc):
        se = 0.05
        return LineageStat(log2fc=lfc, se=se, df=40,
                           ci95=(lfc - 2 * se, lfc + 2 * se),
                           ci90=(lfc - 1.7 * se, lfc + 1.7 * se),
                           null_change=0.05, null_equiv=0.0)
    sbl = {"NK": stat(1.0), "CD8T": stat(-1.0), "CD4T": stat(1.0),
           "B": stat(-1.0), "Myeloid": stat(1.0)}
    assert classify_gene("g", sbl, delta=0.5).quadrant == "unclassified"


# ---------------------------------------------------------------- the core claim
@pytest.mark.parametrize("nk_cells,max_false_tost", [(3, 2), (5, 3)])
def test_low_nk_power_does_not_manufacture_q4(nk_cells, max_false_tost):
    """The reason this protocol exists.

    Every planted gene here is a true Q3: all five lineages move,
    including NK. NK is then starved of cells. Any Q4 call is false by
    construction. The naive rule ("NK p > 0.05, therefore unchanged")
    should fail catastrophically; TOST should hold.
    """
    sim = S.simulate(n_patients=12, n_genes_per_class=60, n_null_genes=0,
                     effect=0.6, dispersion=1.2, rng=np.random.default_rng(21),
                     cells_per_lineage={"NK": nk_cells, "CD8T": 400, "CD4T": 400,
                                        "B": 400, "Myeloid": 400})
    q3 = np.array([g.startswith("Q3") for g in sim.genes])
    sim.counts, sim.genes = sim.counts[q3], sim.genes[q3]

    _, de = _pipeline(sim, balanced=False, min_count=5)
    nk, delta = de["NK"], 0.5

    witnesses = np.sum(
        [(((de[l].ci_lo_95 > 0) | (de[l].ci_hi_95 < 0)) & (np.abs(de[l].log2fc) >= delta))
         for l in ("CD8T", "B", "Myeloid")], axis=0) >= 2

    naive_unchanged = nk.pvalue > 0.05
    tost_equiv = np.array([tost(nk.log2fc[i], nk.se[i], nk.df_total, delta)[1]
                           for i in range(len(nk.genes))])

    false_naive = int((naive_unchanged & witnesses).sum())
    false_tost = int((tost_equiv & witnesses).sum())

    assert false_naive > 10, (
        "the naive rule was expected to fail here; if it did not, the "
        "simulation no longer reproduces the underpowered regime"
    )
    assert false_tost <= max_false_tost, (
        f"TOST produced {false_tost} false Q4 calls at {nk_cells} NK cells"
    )
    assert false_tost < false_naive


def test_se_balance_stop_rule_fires_when_nk_is_starved():
    """Stop rule S3 must actually trigger on a dataset that cannot
    support a Q4 claim."""
    sim = S.simulate(n_patients=10, n_genes_per_class=20, n_null_genes=40,
                     dispersion=1.2, rng=np.random.default_rng(5),
                     cells_per_lineage={"NK": 4, "CD8T": 400, "CD4T": 400,
                                        "B": 400, "Myeloid": 400})
    _, de = _pipeline(sim, balanced=False, min_count=5)
    rep = se_balance_report(de)
    assert not rep["passes"]
    assert "STOP" in rep["verdict"]
