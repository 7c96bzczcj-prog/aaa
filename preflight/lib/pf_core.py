"""BM-PB-NK preflight -- shared instrument.

Everything that decides a verdict is a frozen constant in this file.  The
spec (section 6) forbids changing any of them after a run; a change needs
a new SPEC_VERSION and a written reason.

The rho estimator and the empty-droplet window are inherited verbatim from
the repository's existing ambient work (`src/nkmine/gating.py`), so the
numbers produced here are on the same scale as the already-published
rho_NK = 0.151 / rho_CD8T = 0.058.
"""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

SPEC_VERSION = "BM-PB-NK-preflight-v1.0"

# ---------------------------------------------------------------- frozen ---
# empty-droplet window for the ambient profile (inherited: CellBender put the
# empty prior at 23 counts on GSE154826 library 48)
EMPTY_UMI_MIN, EMPTY_UMI_MAX = 1, 20
# cell calling (inherited from read_10x_tar defaults)
CELL_MIN_UMI, CELL_MIN_GENES = 500, 200
# admission criterion 2
DETECTION_LO, DETECTION_HI = 0.05, 0.85
# lineage assignment: winner must beat runner-up by this z-score margin
ASSIGN_MARGIN = 0.15
# minimum NK cells for a sample to contribute a donor-level value
MIN_NK_CELLS = 30
# task A judgement thresholds
DELTA_AMBIENT_EFFECT = 0.10      # |Delta_ambient| below this is "near zero"
FDR_ALPHA = 0.05
RHO_RATIO_LARGE = 2.0            # rho_BM >> rho_PB
MAJORITY_FRAC = 0.50             # "most genes"
MIN_ADMITTED_DATASETS = 2        # below this task A is unmeasurable
# task C
TASK_C_X_TARGET = 0.001          # 0.1% -- X must be at or below this to conclude

# rho markers: genes the lineage does not express, so their counts in that
# lineage are entirely ambient.  Primary set inherited verbatim.
RHO_MARKERS_NK = ["IGKC", "IGHG1", "IGHA1"]
# secondary, wider set -- reported as a sensitivity arm only, never as the
# primary number.  Bone marrow soup is rich in erythroid and granulocyte
# transcripts that blood soup lacks, so a marrow-only marker set would not be
# comparable across compartments; these are genes NK lacks in BOTH.
RHO_MARKERS_NK_WIDE = ["IGKC", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGHA2",
                       "IGHM", "IGLC1", "IGLC2", "IGLC3", "JCHAIN", "MZB1",
                       "DERL3", "TNFRSF17", "MS4A1", "CD79A", "LYZ", "CD14",
                       "S100A8", "S100A9", "ELANE", "MPO", "PRTN3", "AZU1",
                       "HBB", "HBD", "ALAS2", "AHSP", "PPBP", "PF4"]

# ------------------------------------------------------- target gene list ---
# section 3.3 of the spec, verbatim and in order
GENES_YANG2019 = ["CX3CR1", "HAVCR2", "ZEB2"]
GENES_MELSEN2018 = ["CD69", "CXCR6", "S1PR1", "SELPLG", "SELL", "TIGIT",
                    "CD96", "CD226", "GZMB", "GZMH", "GNLY"]
GENES_CLUSTERING = ["NCAM1", "FCGR3A", "KLRC1", "KLRD1"]
TARGET_GENES = GENES_YANG2019 + GENES_MELSEN2018 + GENES_CLUSTERING
GENE_ALIAS = {"CD226": "DNAM1"}

# ------------------------------------------------------------- panels -----
# The NK panel deliberately excludes every gene in TARGET_GENES so that the
# ambient readout is not measured on the genes that defined the gate.  It
# cannot exclude the *negative* constraints, which is stated as a limit.
NK_PANEL = ["KLRF1", "NKG7", "PRF1", "CTSW", "GZMA", "SPON2", "CLIC3", "XCL1"]
T_PANEL = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "IL32", "CD2"]
PANELS = {
    "NK": NK_PANEL,
    "T": T_PANEL,
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "TCL1A", "VPREB1", "DNTT", "IGLL1"],
    "Plasma": ["MZB1", "JCHAIN", "DERL3", "TNFRSF17", "XBP1"],
    "Myeloid": ["LYZ", "CD14", "FCN1", "S100A8", "S100A9", "AIF1", "CST3", "MNDA"],
    "GranProg": ["ELANE", "MPO", "PRTN3", "AZU1", "DEFA4", "CTSG", "RNASE2"],
    "Erythroid": ["HBB", "HBD", "ALAS2", "AHSP", "CA1", "GYPA", "SLC4A1"],
    "HSPC": ["CD34", "SPINK2", "PRSS57", "EGFL7", "SOX4", "CRHBP", "NPR3"],
    "Mk": ["PPBP", "PF4", "ITGA2B", "GP9"],
    "pDC": ["LILRA4", "CLEC4C", "SCT", "PTCRA", "IRF8"],
}
# an NK call additionally requires zero counts of every one of these
CD3_ZERO = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"]


# --------------------------------------------------------------- readers ---
class Matrix:
    """genes x barcodes CSC, plus gene names and barcodes."""

    def __init__(self, mat, genes, barcodes, name):
        self.mat = mat.tocsc()
        self.genes = np.asarray(genes, dtype=str)
        self.barcodes = np.asarray(barcodes, dtype=str)
        self.name = name

    @property
    def umi(self):
        if not hasattr(self, "_umi"):
            self._umi = np.asarray(self.mat.sum(axis=0)).ravel()
        return self._umi

    @property
    def n_genes_det(self):
        if not hasattr(self, "_ng"):
            self._ng = np.asarray((self.mat > 0).sum(axis=0)).ravel()
        return self._ng


def read_10x_h5(path: str) -> Matrix:
    import h5py
    with h5py.File(path, "r") as f:
        g = f["matrix"] if "matrix" in f else f[list(f.keys())[0]]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])
        bc = np.array([b.decode() for b in g["barcodes"][:]])
        if "features" in g:
            names = np.array([b.decode() for b in g["features"]["name"][:]])
            ftype = np.array([b.decode() for b in g["features"]["feature_type"][:]])
        else:
            names = np.array([b.decode() for b in g["gene_names"][:]])
            ftype = np.array(["Gene Expression"] * len(names))
    m = sp.csc_matrix((data, indices, indptr), shape=shape)
    keep = ftype == "Gene Expression"
    return Matrix(m[keep], names[keep], bc, os.path.basename(path))


def read_10x_mtx(barcodes_gz: str, features_gz: str, matrix_gz: str, name: str) -> Matrix:
    from scipy.io import mmread
    with gzip.open(barcodes_gz, "rt") as fh:
        bc = np.array([l.strip() for l in fh])
    feats = pd.read_csv(features_gz, sep="\t", header=None, compression="gzip")
    names = feats[1].astype(str).values if feats.shape[1] > 1 else feats[0].astype(str).values
    ftype = feats[2].astype(str).values if feats.shape[1] > 2 else np.array(
        ["Gene Expression"] * len(names))
    with gzip.open(matrix_gz, "rb") as fh:
        m = mmread(fh).tocsc()
    keep = ftype == "Gene Expression"
    return Matrix(m[keep], names[keep], bc, name)


# --------------------------------------------------------------- ambient ---
def soup_profile(m: Matrix, min_umi=EMPTY_UMI_MIN, max_umi=EMPTY_UMI_MAX):
    """Ambient profile from empty droplets; returns (profile summing to 1, n)."""
    umi = m.umi
    empty = (umi >= min_umi) & (umi <= max_umi)
    if empty.sum() == 0:
        raise ValueError(f"{m.name}: no droplets in the empty window")
    prof = np.asarray(m.mat[:, empty].sum(axis=1)).ravel().astype(float)
    tot = prof.sum()
    if tot <= 0:
        raise ValueError(f"{m.name}: empty droplets carry no counts")
    return prof / tot, int(empty.sum())


def estimate_rho(observed: np.ndarray, genes: np.ndarray, profile: np.ndarray,
                 markers: list) -> float:
    """Inherited verbatim: rho = observed marker counts / counts expected if
    the whole profile were soup."""
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(g, i)
    rows = [idx[g] for g in markers if g in idx]
    if not rows:
        return float("nan")
    obs = float(observed[rows].sum())
    exp = float(observed.sum()) * float(profile[rows].sum())
    if exp <= 0:
        return float("nan")
    return float(np.clip(obs / exp, 0.0, 1.0))


def ambient_fraction(observed: np.ndarray, profile: np.ndarray, rho: float,
                     rows: list) -> np.ndarray:
    """Per gene: the share of that gene's counts in the pseudobulk that the
    ambient model attributes to soup.  Clipped to [0, 1]; NaN where the gene
    has no counts at all (nothing to apportion)."""
    tot = float(observed.sum())
    out = np.full(len(rows), np.nan)
    for k, r in enumerate(rows):
        obs = float(observed[r])
        if obs <= 0:
            continue
        out[k] = min(1.0, rho * tot * float(profile[r]) / obs)
    return out


# -------------------------------------------------------------- gating ----
def _gene_rows(genes: np.ndarray) -> dict:
    d = {}
    for i, g in enumerate(genes):
        d.setdefault(g, []).append(i)
    return d


def score_panels(m: Matrix, cells: np.ndarray) -> dict:
    x = m.mat[:, cells]
    lib = np.asarray(x.sum(axis=0)).ravel().astype(float)
    lib[lib == 0] = 1.0
    g2i = _gene_rows(m.genes)
    out = {}
    for k, panel in PANELS.items():
        rows = [i for g in panel for i in g2i.get(g, [])]
        if not rows:
            out[k] = np.zeros(x.shape[1])
            continue
        sub = np.asarray(x[rows].todense())
        norm = np.log1p(sub / lib * 1e4)
        s = norm.mean(axis=0)
        sd = s.std()
        out[k] = (s - s.mean()) / sd if sd > 0 else s * 0.0
    return out


def call_cells(m: Matrix) -> np.ndarray:
    return np.flatnonzero((m.umi >= CELL_MIN_UMI) & (m.n_genes_det >= CELL_MIN_GENES))


def assign_lineage(m: Matrix, cells: np.ndarray) -> np.ndarray:
    sc = score_panels(m, cells)
    keys = list(sc)
    mat = np.stack([sc[k] for k in keys])
    order = np.argsort(mat, axis=0)
    best = mat[order[-1], np.arange(mat.shape[1])]
    second = mat[order[-2], np.arange(mat.shape[1])]
    call = np.array([keys[i] for i in order[-1]], dtype=object)
    call[(best - second) < ASSIGN_MARGIN] = "unassigned"
    # an NK call requires zero T-receptor / CD3 transcripts
    g2i = _gene_rows(m.genes)
    rows = [i for g in CD3_ZERO for i in g2i.get(g, [])]
    if rows:
        cd3 = np.asarray(m.mat[rows][:, cells].todense())
        call[(call == "NK") & (cd3 > 0).any(axis=0)] = "unassigned"
    return call.astype(str)


def detection_rate(m: Matrix, cells: np.ndarray, gene: str) -> float:
    g2i = _gene_rows(m.genes)
    rows = g2i.get(gene)
    if not rows:
        return float("nan")
    sub = m.mat[rows][:, cells]
    return float((np.asarray(sub.sum(axis=0)).ravel() > 0).mean())


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    v = p[ok]
    n = len(v)
    if n == 0:
        return q
    order = np.argsort(v)
    ranked = v[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def clopper_pearson_upper(k: int, n: int, alpha: float = 0.05) -> float:
    """One-sided 95% upper bound on a binomial proportion."""
    from scipy.stats import beta
    if n == 0:
        return float("nan")
    if k >= n:
        return 1.0
    return float(beta.ppf(1 - alpha, k + 1, n - k))


# ---------------------------------------------------------------------------
# Instrument acceptance (added before any Delta_ambient was computed; the
# reason is recorded in docs/TRANSFERABLE.md).
#
# The first acceptance rule -- "the ambient share of MS4A1/CD79A/LYZ/HBB in NK
# must exceed 0.5" -- failed on bone marrow (MS4A1 median 0.24).  The cause was
# diagnosed before any endpoint was read: on MantonBM1 lane 1 the NK pseudobulk
# carried 8 MS4A1 counts of which 5 sat in ONE cell, and 42 LYZ counts of which
# 35 sat in one cell.  Those are doublets, not a broken soup model: the per-gene
# ambient ratios r_g of the foreign panel concentrate tightly at rho
# (median 0.0123 in BM, 0.0145 in PB, against rho of 0.0130 / 0.0119).
#
# Two consequences, applied to every dataset and both compartments alike:
#   * predicted doublets are removed before any gate is formed;
#   * acceptance is judged on the SPREAD of r_g around rho over foreign genes
#     that carry enough predicted ambient counts to be more than Poisson noise,
#     instead of on four hand-picked genes.
RHO_FOREIGN_PANEL = [
    "IGKC", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGHA2", "IGHM", "IGLC1",
    "IGLC2", "IGLC3", "JCHAIN", "MZB1", "DERL3", "TNFRSF17",
    "MS4A1", "CD79A", "CD79B", "BANK1", "VPREB1", "DNTT",
    "LYZ", "CD14", "FCN1", "S100A8", "S100A9", "MNDA",
    "ELANE", "MPO", "PRTN3", "AZU1", "CTSG",
    "HBB", "HBD", "ALAS2", "AHSP", "GYPA", "CA1",
    "PF4", "PPBP", "ITGA2B", "LILRA4", "CLEC4C", "CD34", "GATA1",
]
# SPINK2 is deliberately NOT in the panel: it is not NK-foreign.  Measured
# r_g was 1.11 (BM) and 1.69 (PB), i.e. NK carry more SPINK2 than the whole
# pseudobulk could hold as soup, so treating it as foreign would inflate rho.
# A gene may only be read as an ambient measurement when the model expects a
# meaningful number of ambient counts in it.  The floor is on rho * pred (the
# counts the soup model actually predicts), NOT on pred (the counts the gene
# would carry if the pseudobulk were pure soup).  The first version used pred,
# which is wrong by a factor of rho: on GSE233304, where rho is about 0.004,
# it declared 17-40 genes testable whose predicted ambient was under 0.1 counts,
# every one of which came back at r_g = 0 and failed acceptance by construction.
MIN_EXPECTED_AMBIENT_COUNTS = 5.0
MIN_TESTABLE_FOREIGN_GENES = 3     # below this, acceptance is indeterminate
ACCEPT_FOLD = 3.0                  # r_g must sit within this factor of rho
ACCEPT_FRAC = 2.0 / 3.0            # and this share of testable genes must


def foreign_obs_pred(observed, genes, profile, panel=None):
    """Per foreign gene: (observed counts, counts if the whole pseudobulk were
    soup).  Returned raw so that libraries can be POOLED before the ratio is
    formed -- a ratio of small counts per library is far noisier than the ratio
    of their sums."""
    panel = RHO_FOREIGN_PANEL if panel is None else panel
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(g, i)
    tot = float(observed.sum())
    obs, pred = {}, {}
    for g in panel:
        i = idx.get(g)
        if i is None:
            continue
        obs[g] = float(observed[i])
        pred[g] = tot * float(profile[i])
    return obs, pred


def per_gene_ratio(observed, genes, profile, panel=None):
    """r_g = observed_g / (total_observed * soup_g).  If the soup model holds,
    every NK-foreign gene returns the same value, namely rho."""
    obs, pred = foreign_obs_pred(observed, genes, profile, panel)
    return {g: ((obs[g] / pred[g]) if pred[g] > 0 else np.nan, pred[g]) for g in obs}


def acceptance_from_pooled(obs, pred, rho):
    """Acceptance on pooled counts.  obs and pred are dicts gene -> summed
    counts; rho is the pooled soup fraction.  Returns
    (passes, n_testable, frac_within, median_ratio)."""
    if not np.isfinite(rho) or rho <= 0:
        return False, 0, np.nan, np.nan
    vals = []
    for g in obs:
        if rho * pred.get(g, 0.0) >= MIN_EXPECTED_AMBIENT_COUNTS and pred[g] > 0:
            vals.append(obs[g] / pred[g])
    if len(vals) < MIN_TESTABLE_FOREIGN_GENES:
        return False, len(vals), np.nan, (float(np.median(vals)) if vals else np.nan)
    v = np.array(vals)
    within = float(np.mean((v <= rho * ACCEPT_FOLD) & (v >= rho / ACCEPT_FOLD)))
    return bool(within >= ACCEPT_FRAC), int(len(v)), within, float(np.median(v))


def acceptance(observed, genes, profile, rho):
    obs, pred = foreign_obs_pred(observed, genes, profile)
    return acceptance_from_pooled(obs, pred, rho)


def rho_median(observed, genes, profile, rho_provisional=None):
    """Alternative rho: median per-gene ratio over the foreign genes whose
    expected ambient counts clear the floor.  The floor needs a provisional rho;
    the inherited narrow estimator supplies it."""
    obs, pred = foreign_obs_pred(observed, genes, profile)
    r = rho_provisional
    if r is None or not np.isfinite(r) or r <= 0:
        return float("nan")
    vals = [obs[g] / pred[g] for g in obs
            if pred.get(g, 0.0) > 0 and r * pred[g] >= MIN_EXPECTED_AMBIENT_COUNTS]
    return float(np.median(vals)) if len(vals) >= MIN_TESTABLE_FOREIGN_GENES else float("nan")
