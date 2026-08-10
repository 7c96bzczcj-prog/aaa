"""Tests for Phases 6-8: artefact filters, replication, controls."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nkmine.artifacts import (  # noqa: E402
    ambient_susceptibility,
    decompose_reference_swap,
    kitagawa_decompose,
)
from nkmine.controls import control_check, evaluate_stop_rule_s4  # noqa: E402
from nkmine.replication import (  # noqa: E402
    evaluate_stop_rules_s5_s6,
    mundane_screen,
    replicate,
    replication_rate_by_quadrant,
)


# ---------------------------------------------------------------- Kitagawa
def test_decomposition_is_exact_and_reports_interaction():
    """The three terms must sum to the total. The protocol's stated bug
    was a CSV with only two columns, whose percentages did not add up."""
    res = kitagawa_decompose([0.6, 0.4], [0.3, 0.7], [1.0, 4.0], [1.5, 5.0])
    assert np.isclose(res["within"] + res["between"] + res["interaction"],
                      res["total"])
    assert "interaction" in res and res["interaction"] != 0


def test_pure_composition_change_has_zero_within_term():
    res = kitagawa_decompose([0.5, 0.5], [0.1, 0.9], [1.0, 3.0], [1.0, 3.0])
    assert np.isclose(res["within"], 0.0)
    assert np.isclose(res["between"], res["total"])
    assert np.isclose(res["interaction"], 0.0)


def test_pure_state_change_has_zero_between_term():
    res = kitagawa_decompose([0.5, 0.5], [0.5, 0.5], [1.0, 3.0], [2.0, 4.0])
    assert np.isclose(res["between"], 0.0)
    assert np.isclose(res["within"], res["total"])


def test_large_interaction_is_flagged_not_identified():
    """Protocol 6.2: when the interaction term dominates, the within /
    between split depends on the reference weights and must not be
    reported as a number."""
    res = kitagawa_decompose([0.9, 0.1], [0.1, 0.9], [0.0, 1.0], [1.0, 6.0])
    assert abs(res["interaction_pct"]) > 25
    assert not res["identified"]
    assert res["flag"] == "not_identified"


def test_reference_swap_exposes_the_index_number_problem():
    swap = decompose_reference_swap([0.9, 0.1], [0.1, 0.9], [0.0, 1.0], [1.0, 6.0])
    assert swap["within_pct_swing"] > 20
    assert not swap["stable"]


def test_proportions_must_be_normalised():
    with pytest.raises(ValueError):
        kitagawa_decompose([0.6, 0.5], [0.3, 0.7], [1.0, 4.0], [1.5, 5.0])


# ---------------------------------------------------------------- ambient
def test_ambient_susceptibility_tracks_abundance_shift():
    expr = pd.DataFrame({"Tumour": [10.0, 0.1], "Bcell": [0.1, 10.0]},
                        index=["TUMGENE", "BGENE"])
    ref = pd.Series({"Tumour": 0.1, "Bcell": 0.5})
    alt = pd.Series({"Tumour": 0.6, "Bcell": 0.5})   # only tumour shifts
    s = ambient_susceptibility(expr, ref, alt)
    assert s["TUMGENE"] > s["BGENE"]
    assert np.isclose(s["BGENE"], 0.0)


# ---------------------------------------------------------------- replication
def test_replication_requires_quadrant_sign_and_magnitude():
    a = pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"],
                      "quadrant": ["Q4", "Q4", "Q4", "Q4"],
                      "NK_log2FC": [0.02, 0.02, 0.02, 0.02]})
    b = pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"],
                      "quadrant": ["Q4", "Q3", "Q4", "Q4"],
                      "NK_log2FC": [0.03, 0.02, -0.02, 0.50]})
    rep = replicate(a, b)
    got = dict(zip(rep.gene, rep.replicated))
    assert got["g1"]            # same quadrant, sign, within 2x
    assert not got["g2"]        # different quadrant
    assert not got["g3"]        # opposite sign
    assert not got["g4"]        # magnitude out by 25x


def test_stop_rule_s6_fires_when_q4_underperforms_baseline():
    rates = pd.DataFrame({"quadrant": ["Q1", "Q3", "Q4"],
                          "n": [100, 100, 100],
                          "n_replicated": [80, 90, 5],
                          "rate": [0.80, 0.90, 0.05]})
    ev = evaluate_stop_rules_s5_s6(rates)
    assert ev["stop"]
    assert any("S5" in p for p in ev["problems"])
    assert any("S6" in p for p in ev["problems"])


def test_healthy_replication_passes():
    rates = pd.DataFrame({"quadrant": ["Q1", "Q3", "Q4"],
                          "n": [100, 100, 100],
                          "n_replicated": [70, 80, 60],
                          "rate": [0.70, 0.80, 0.60]})
    assert not evaluate_stop_rules_s5_s6(rates)["stop"]


def test_replication_rate_table_covers_every_quadrant():
    a = pd.DataFrame({"gene": list("abcd"), "quadrant": ["Q1", "Q3", "Q4", "Q4"],
                      "NK_log2FC": [1.0, 1.0, 0.01, 0.01]})
    b = pd.DataFrame({"gene": list("abcd"), "quadrant": ["Q1", "Q3", "Q4", "Q1"],
                      "NK_log2FC": [1.1, 0.9, 0.012, 0.01]})
    r = replication_rate_by_quadrant(replicate(a, b))
    assert set(r.quadrant) == {"Q1", "Q3", "Q4"}
    assert float(r[r.quadrant == "Q4"].rate.iloc[0]) == 0.5


# ---------------------------------------------------------------- controls
def test_control_check_flags_ieg_landing_in_q4():
    q = pd.DataFrame({"gene": ["FOS", "TOX"], "quadrant": ["Q4", "Q4"]})
    ev = evaluate_stop_rule_s4(control_check(q))
    assert ev["stop"]
    assert any("immediate-early" in p for p in ev["problems"])


def test_control_check_flags_tox_missing_from_q4():
    q = pd.DataFrame({"gene": ["TOX", "FOS"], "quadrant": ["unclassified", "Q3"]})
    ev = evaluate_stop_rule_s4(control_check(q))
    assert ev["stop"]
    assert any("TOX" in p for p in ev["problems"])


def test_controls_pass_when_everything_lands_correctly():
    q = pd.DataFrame({"gene": ["TOX", "FOS", "JUN", "GZMB", "ACTB"],
                      "quadrant": ["Q4", "Q3", "Q3", "Q2", "unclassified"]})
    ev = evaluate_stop_rule_s4(control_check(q))
    assert not ev["stop"], ev["problems"]
    assert ev["pass_rate"] == 1.0


def test_gating_genes_are_reported_not_silently_dropped():
    """The Q2 control list overlaps the NK gating panel; a control list
    that quietly shrinks looks the same as one that passes."""
    q = pd.DataFrame({"gene": ["TOX"], "quadrant": ["Q4"]})
    chk = control_check(q, excluded={"GNLY", "PRF1"})
    assert set(chk[chk.status == "excluded_by_gating"].gene) == {"GNLY", "PRF1"}


# ---------------------------------------------------------------- Phase 8
def test_gene_survives_only_when_every_mundane_cause_is_refuted():
    cand = pd.DataFrame({"gene": ["good", "ambient", "partial"]})
    refuted = {g: False for g in ["good", "ambient", "partial"]}
    evidence = {code: dict(refuted) for code in
                ["B1_not_exposed", "B2_fast_turnover", "B3_missing_receptor",
                 "B4_ceiling_floor", "B5_underpowered", "B6_ambient"]}
    evidence["B6_ambient"]["ambient"] = True      # explained away
    del evidence["B1_not_exposed"]["partial"]     # never assessed

    out = mundane_screen(cand, evidence).set_index("gene")
    assert out.loc["good", "survives"]
    assert not out.loc["ambient", "survives"]
    assert "B6_ambient" in out.loc["ambient", "status"]
    assert not out.loc["partial", "survives"]
    assert out.loc["partial", "status"] == "incomplete"


def test_unassessed_blocks_survival_rather_than_permitting_it():
    cand = pd.DataFrame({"gene": ["x"]})
    out = mundane_screen(cand, {}).set_index("gene")
    assert not out.loc["x", "survives"]
    assert out.loc["x", "n_unassessed"] == 6


# ---------------------------------------------------------------- D7
def test_detection_filter_can_keep_lineage_restricted_genes():
    """Regression test for D7.

    The strict all-lineage floor removes any T-restricted gene, which on
    real data deleted the entire Q4 positive-control panel and made stop
    rule S4 unsatisfiable. The relaxed form must keep a gene that is
    well measured in NK and two other lineages but absent from a third.
    """
    import numpy as np
    from nkmine.pseudobulk import PseudobulkSet, detection_filter

    lineages = ["NK", "CD8T", "CD4T", "B", "Myeloid"]
    counts, lin, pat, cond = [], [], [], []
    for l in lineages:
        for p in ("p1", "p2"):
            for c in ("Tumor", "Normal"):
                # gene 0 everywhere; gene 1 absent from B and Myeloid
                g1 = 0 if l in ("B", "Myeloid") else 100
                counts.append([100, g1])
                lin.append(l); pat.append(p); cond.append(c)
    ps = PseudobulkSet(np.array(counts).T.astype(float),
                       np.array(["ubiquitous", "T_restricted"]),
                       np.array(pat), np.array(cond), np.array(lin),
                       np.ones(len(lin)))

    strict = detection_filter(ps)
    assert list(strict) == [True, False]

    relaxed = detection_filter(ps, require_nk=True, min_lineages=3)
    assert list(relaxed) == [True, True]

    # a gene undetectable in NK must still be dropped: "NK did not
    # change" must never be a restatement of "NK does not express it"
    ps.counts[1, ps.lineage == "NK"] = 0
    assert not detection_filter(ps, require_nk=True, min_lineages=2)[1]
