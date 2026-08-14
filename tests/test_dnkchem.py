"""Tests for the dNK chemokine pipeline.

Run:  .venv/bin/python -m pytest tests/test_dnkchem.py -q
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from dnkchem.counts import (detection_and_cpm, downsample_columns,  # noqa: E402
                            triple_positive_fraction, verify_integer_counts)
from dnkchem.manifest import (PanelMatchError, match_panel,  # noqa: E402
                              split_var_names)
from dnkchem.stats import (benjamini_hochberg, continuity_rate,  # noqa: E402
                           exact_wilcoxon_signed_rank, kitagawa_decomposition,
                           logit, min_achievable_p)


# --- R2: the floor is real and the test actually reaches it ----------------
@pytest.mark.parametrize("n,expected", [(3, 0.25), (4, 0.125), (5, 0.0625),
                                        (6, 0.03125), (7, 0.015625)])
def test_min_achievable_p_is_attained(n, expected):
    assert min_achievable_p(n) == pytest.approx(expected)
    _, p, _ = exact_wilcoxon_signed_rank(np.ones(n))
    assert p == pytest.approx(expected), "all-same-sign must hit the floor exactly"


def test_p_below_05_impossible_at_n4():
    """The pre-registered claim: analysis B cannot reach p < 0.05."""
    assert min_achievable_p(4) > 0.05


def test_exact_wilcoxon_matches_scipy():
    rng = np.random.default_rng(0)
    for n in (4, 5, 6, 7, 8, 10):
        for _ in range(25):
            d = rng.normal(size=n)
            _, p, _ = exact_wilcoxon_signed_rank(d)
            assert p == pytest.approx(wilcoxon(d, method="exact").pvalue)


def test_wilcoxon_drops_zero_differences():
    _, p, n_used = exact_wilcoxon_signed_rank([0.0, 1.0, 2.0, 3.0])
    assert n_used == 3
    assert p == pytest.approx(0.25)


def test_bh_is_monotone_and_bounded():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
    q = benjamini_hochberg(p)
    assert np.all(np.diff(q) >= -1e-12)
    assert np.all(q <= 1.0)
    assert np.all(q >= p - 1e-12)


def test_bh_carries_nan_through():
    q = benjamini_hochberg(np.array([0.01, np.nan, 0.5]))
    assert np.isnan(q[1]) and np.isfinite(q[0])


# --- downsampling ----------------------------------------------------------
def test_downsample_never_exceeds_depth_or_available():
    rng = np.random.default_rng(1)
    X = sp.csr_matrix(rng.poisson(2.0, size=(200, 40)).astype(np.int64))
    cols = list(range(40))
    counts, keep = downsample_columns(X, cols, depth=20, seed=3)
    assert counts.sum(axis=1).max() <= 20
    orig = np.asarray(X[keep].todense())
    assert np.all(counts <= orig), "cannot draw more copies of a gene than exist"


def test_downsample_drops_cells_below_depth():
    X = sp.csr_matrix(np.array([[10, 0], [1, 1], [50, 50]], dtype=np.int64))
    counts, keep = downsample_columns(X, [0, 1], depth=10, seed=0)
    assert keep.tolist() == [True, False, True]
    assert counts.shape[0] == 2


def test_downsample_marginal_matches_hypergeometric_expectation():
    """A cell thinned to D reads should yield D * (c_g / T) copies of gene g."""
    c = np.array([[300, 100, 600]], dtype=np.int64)
    X = sp.csr_matrix(np.repeat(c, 4000, axis=0))
    counts, _ = downsample_columns(X, [0, 1, 2], depth=100, seed=7)
    got = counts.mean(axis=0)
    want = 100 * c[0] / c[0].sum()
    assert np.allclose(got, want, rtol=0.05), (got, want)
    assert np.all(counts.sum(axis=1) == 100), "every retained cell holds exactly D reads"


def test_downsample_is_seed_reproducible():
    rng = np.random.default_rng(2)
    X = sp.csr_matrix(rng.poisson(3.0, size=(100, 10)).astype(np.int64))
    a, _ = downsample_columns(X, list(range(10)), 15, seed=42)
    b, _ = downsample_columns(X, list(range(10)), 15, seed=42)
    assert np.array_equal(a, b)


def test_detection_rate_counts_at_least_one():
    counts = np.array([[0, 1], [0, 3], [1, 0], [0, 0]])
    n_det, det, cpm = detection_and_cpm(counts, depth=10)
    assert n_det.tolist() == [1, 2]
    assert det.tolist() == [0.25, 0.5]
    assert cpm[1] == pytest.approx((1 + 3) / 4 / 10 * 1e6)


# --- purity gate -----------------------------------------------------------
def test_triple_positive_requires_all_three():
    X = sp.csr_matrix(np.array([[1, 1, 1], [1, 1, 0], [0, 0, 0], [5, 2, 9]]))
    frac, n = triple_positive_fraction(X, [0, 1, 2])
    assert n == 2 and frac == pytest.approx(0.5)


# --- gene identity ---------------------------------------------------------
class _MF:
    raw = {"gene_id": {"var_index_type": "symbol_ensembl_concat",
                       "id_separator": "_", "symbol_field": 0, "ensembl_field": -1}}


def test_split_handles_symbols_containing_the_separator():
    names = ["CCL5_ENSG00000271503", "HLA-G_ENSG00000204632",
             "HLA_DRA_ENSG00000204287"]
    syms, ens = split_var_names(names, _MF())
    assert syms == ["CCL5", "HLA-G", "HLA_DRA"]
    assert ens == ["ENSG00000271503", "ENSG00000204632", "ENSG00000204287"]


def test_panel_match_raises_below_threshold():
    """R11: a near-total miss must abort, never print a table anyway."""
    with pytest.raises(PanelMatchError):
        match_panel(["A", "B", "C", "D"], ["", "", "", ""],
                    ["Z", "Y"], ["", ""], min_rate=0.90, label="t")


def test_panel_match_falls_back_to_ensembl():
    hits, misses, rate, how = match_panel(
        ["CCL5"], ["ENSG00000271503"], ["SOMETHINGELSE"], ["ENSG00000271503.5"],
        min_rate=0.5, label="t")
    assert hits["CCL5"] == 0 and how["CCL5"] == "ensembl" and rate == 1.0


# --- misc ------------------------------------------------------------------
def test_continuity_rate_keeps_logit_finite():
    assert np.isfinite(logit(continuity_rate(0, 50)))
    assert np.isfinite(logit(continuity_rate(50, 50)))


def test_kitagawa_three_terms_sum_exactly():
    rng = np.random.default_rng(5)
    for _ in range(50):
        ca = rng.dirichlet(np.ones(4)); cb = rng.dirichlet(np.ones(4))
        ra = rng.random(4); rb = rng.random(4)
        k = kitagawa_decomposition(ra, ca, rb, cb)
        assert abs(k["residual"]) < 1e-6
        assert abs(k["within"] + k["between"] + k["interaction"] - k["total"]) < 1e-6


def test_verify_integer_counts_rejects_normalised_matrix():
    ok, _, _ = verify_integer_counts(sp.csr_matrix(np.array([[1.5, 2.0], [0.0, 3.3]])))
    assert not ok
    ok2, _, _ = verify_integer_counts(sp.csr_matrix(np.array([[1, 2], [0, 3]])))
    assert ok2


# --- marginal downsampler (used by the 10k-gene null) ----------------------
def test_marginal_downsampler_matches_sequential_marginals():
    """The whole point: same per-gene distribution, different joint."""
    from dnkchem.counts import downsample_marginal
    c = np.array([[400, 100, 500]], dtype=np.int64)
    X = sp.csr_matrix(np.repeat(c, 6000, axis=0))
    seq, _ = downsample_columns(X, [0, 1, 2], depth=100, seed=11)
    mar, _ = downsample_marginal(X, [0, 1, 2], depth=100, seed=11)
    want = 100 * c[0] / c[0].sum()
    assert np.allclose(seq.mean(axis=0), want, rtol=0.05)
    assert np.allclose(mar.mean(axis=0), want, rtol=0.05), mar.mean(axis=0)
    # detection rates agree, which is what the null actually consumes
    assert np.allclose((seq >= 1).mean(axis=0), (mar >= 1).mean(axis=0), atol=0.02)
    # the joint differs by construction: sequential sums to exactly D
    assert np.all(seq.sum(axis=1) == 100)
    assert not np.all(mar.sum(axis=1) == 100)


def test_marginal_downsampler_block_size_changes_draws_not_distribution():
    """Block size changes the RNG call sequence, so per-gene draws differ.

    What must NOT differ is the distribution they come from. Asserting
    per-gene equality would be asserting something false; assert the pooled
    detection rate instead, which is what the null consumes.
    """
    from dnkchem.counts import downsample_marginal
    rng = np.random.default_rng(4)
    X = sp.csr_matrix(rng.poisson(2.0, size=(800, 60)).astype(np.int64))
    a, _ = downsample_marginal(X, list(range(60)), 20, seed=7, block=60)
    b, _ = downsample_marginal(X, list(range(60)), 20, seed=7, block=7)
    assert a.shape == b.shape
    assert abs((a >= 1).mean() - (b >= 1).mean()) < 0.02
    assert abs(a.mean() - b.mean()) < 0.05


def test_marginal_downsampler_is_seed_reproducible_at_fixed_block():
    from dnkchem.counts import downsample_marginal
    rng = np.random.default_rng(9)
    X = sp.csr_matrix(rng.poisson(2.0, size=(200, 30)).astype(np.int64))
    a, _ = downsample_marginal(X, list(range(30)), 15, seed=3, block=10)
    b, _ = downsample_marginal(X, list(range(30)), 15, seed=3, block=10)
    assert np.array_equal(a, b)
