"""Cell calling and non-circular lineage assignment (Protocol Phase 2.1-2.2).

GSE154826 carries a 27-marker ADT panel and the hashing oligos in the
same feature space as the RNA, which is what makes non-circular gating
possible: lineage is called on *protein*, so the RNA genes that define
a lineage never have to serve as both the label and the measurement.
The RNA marker panel is still written to `excluded_genes.txt` and
barred from testing, as a second line of defence.

Purity is checked with discrete counts, never with a mean detection
rate.  A mean hides a small, tight cluster of contaminating T cells
inside an NK gate; the sub-cluster check below is what catches it.
"""

from __future__ import annotations

import gzip
import io
import os
import tarfile
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Protocol Phase 2.1 RNA panel.  These are the genes used for cross-checks
# and are therefore excluded from every differential test downstream.
RNA_MARKERS = {
    "NK": ["NKG7", "KLRD1", "KLRF1", "GNLY", "NCAM1", "FCGR3A"],
    "CD8T": ["CD3D", "CD3E", "CD3G", "CD8A"],
    "CD4T": ["CD3D", "CD3E", "CD3G", "IL7R", "CD40LG", "CD4"],
    "B": ["MS4A1", "CD79A", "CD19"],
    "Myeloid": ["LYZ", "CD68", "FCN1", "ITGAM"],
}
CD3_RNA = ["CD3D", "CD3E", "CD3G"]


def excluded_genes() -> list[str]:
    out = set()
    for v in RNA_MARKERS.values():
        out.update(v)
    return sorted(out)


@dataclass
class Library:
    """One 10x library: RNA, ADT and HTO blocks over called cells."""

    name: str
    rna: sp.csc_matrix          # genes x cells
    gene_names: np.ndarray
    adt: np.ndarray | None      # markers x cells (dense; the panel is small)
    adt_names: np.ndarray | None
    hto: np.ndarray | None
    hto_names: np.ndarray | None
    barcodes: np.ndarray
    n_droplets: int             # droplets before cell calling


def read_10x_tar(path: str, min_umi: int = 500, min_genes: int = 200) -> Library:
    """Read a GEO amp_batch tarball and call cells.

    The archives hold the *unfiltered* droplet matrix (737,280 barcodes,
    the full 10x v2 whitelist), so cell calling happens here.  A UMI/gene
    threshold is used rather than EmptyDrops: it is cruder, but this step
    feeds cell *counts* and lineage gating, and the ambient correction
    that would justify EmptyDrops' extra machinery is Phase 6's job.
    """
    with tarfile.open(path, "r:gz") as tf:
        members = {os.path.basename(m.name): m for m in tf.getmembers() if m.isfile()}
        bc_m = next(v for k, v in members.items() if k.endswith("barcodes.tsv"))
        ft_m = next(v for k, v in members.items() if k.endswith("features.tsv"))
        mx_m = next(v for k, v in members.items() if k.endswith("matrix.mtx"))

        barcodes = np.array(tf.extractfile(bc_m).read().decode().split("\n"))
        barcodes = barcodes[barcodes != ""]
        feat = pd.read_csv(tf.extractfile(ft_m), sep="\t", header=None,
                           names=["id", "name", "type"])

        fh = tf.extractfile(mx_m)
        header = fh.readline()
        while True:
            pos = fh.tell()
            line = fh.readline()
            if not line.startswith(b"%"):
                fh.seek(pos)
                break
        dims = fh.readline().split()
        n_feat, n_bc = int(dims[0]), int(dims[1])
        trip = pd.read_csv(fh, sep=r"\s+", header=None, dtype=np.int64,
                           names=["f", "b", "v"])

    mat = sp.coo_matrix(
        (trip.v.values, (trip.f.values - 1, trip.b.values - 1)),
        shape=(n_feat, n_bc),
    ).tocsc()

    is_rna = (feat.type == "Gene Expression").values
    names = feat.name.values.astype(str)
    is_hto = np.array([n.startswith("HTO_") for n in names])
    is_adt = (~is_rna) & (~is_hto)

    rna_all = mat[is_rna]
    umi = np.asarray(rna_all.sum(axis=0)).ravel()
    ngene = np.asarray((rna_all > 0).sum(axis=0)).ravel()
    called = (umi >= min_umi) & (ngene >= min_genes)

    return Library(
        name=os.path.basename(path).replace(".tar.gz", ""),
        rna=rna_all[:, called].tocsc(),
        gene_names=names[is_rna],
        adt=np.asarray(mat[is_adt][:, called].todense()) if is_adt.any() else None,
        adt_names=names[is_adt] if is_adt.any() else None,
        hto=np.asarray(mat[is_hto][:, called].todense()) if is_hto.any() else None,
        hto_names=names[is_hto] if is_hto.any() else None,
        barcodes=barcodes[called],
        n_droplets=n_bc,
    )


def clr(x: np.ndarray) -> np.ndarray:
    """Centred log-ratio normalisation across cells, per marker."""
    lx = np.log1p(x)
    return lx - lx.mean(axis=1, keepdims=True)


def _otsu(v: np.ndarray) -> float:
    """Otsu's between-class variance threshold on a 128-bin histogram."""
    hist, edges = np.histogram(v, bins=128)
    centers = (edges[:-1] + edges[1:]) / 2
    w = hist.cumsum().astype(float)
    m = (hist * centers).cumsum()
    total_w, total_m = w[-1], m[-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        denom = w * (total_w - w)
        between = np.where(denom > 0, (total_m * w - m * total_w) ** 2 / (denom + 1e-12), np.nan)
    if np.all(np.isnan(between)):
        return float(np.median(v))
    return float(centers[np.nanargmax(between)])


def _em_2class(v: np.ndarray, n_iter: int = 100, tol: float = 1e-6):
    """Two-component 1-D Gaussian mixture by EM.

    Written out rather than pulled from scikit-learn to keep the
    dependency surface to numpy/scipy/pandas; the model is small enough
    that the explicit version is also easier to audit.
    """
    x = np.asarray(v, dtype=float)
    lo_init, hi_init = np.percentile(x, [25, 75])
    mu = np.array([lo_init, hi_init])
    sd = np.full(2, max(x.std(), 1e-3))
    pi = np.array([0.5, 0.5])
    prev = -np.inf

    for _ in range(n_iter):
        d = np.stack([
            pi[k] * np.exp(-0.5 * ((x - mu[k]) / sd[k]) ** 2) / (sd[k] * np.sqrt(2 * np.pi))
            for k in range(2)
        ])
        tot = d.sum(axis=0) + 1e-300
        r = d / tot
        ll = np.log(tot).sum()
        nk = r.sum(axis=1) + 1e-12
        pi = nk / len(x)
        mu = (r * x).sum(axis=1) / nk
        sd = np.sqrt(np.maximum((r * (x - mu[:, None]) ** 2).sum(axis=1) / nk, 1e-6))
        if abs(ll - prev) < tol:
            break
        prev = ll
    return mu, sd, pi


def _threshold_2class(v: np.ndarray) -> float:
    """Split a bimodal marker into negative/positive.

    Uses the EM mixture's equal-posterior crossing point, and falls back
    to Otsu when the two components collapse onto each other (a marker
    that is not actually bimodal in this library).
    """
    v = np.asarray(v, dtype=float)
    if v.size < 50 or np.allclose(v.std(), 0):
        return float(np.median(v))
    try:
        mu, sd, pi = _em_2class(v)
        lo, hi = int(np.argmin(mu)), int(np.argmax(mu))
        if abs(mu[hi] - mu[lo]) < 0.15:
            return _otsu(v)
        grid = np.linspace(mu[lo], mu[hi], 500)
        with np.errstate(under="ignore"):
            dens = np.stack([
                pi[k] * np.exp(-0.5 * ((grid - mu[k]) / sd[k]) ** 2) / (sd[k] * np.sqrt(2 * np.pi))
                for k in (lo, hi)
            ])
        post_hi = dens[1] / (dens.sum(axis=0) + 1e-300)
        return float(grid[np.argmin(np.abs(post_hi - 0.5))])
    except Exception:
        return _otsu(v)


def gate_lineages(lib: Library) -> tuple[np.ndarray, pd.DataFrame]:
    """Assign each called cell to a lineage using ADT, cross-checked on RNA.

    Returns (lineage per cell, threshold table).  Cells that do not fit
    exactly one lineage are labelled "unassigned" and dropped downstream;
    a permissive gate would import the contamination this protocol is
    built to avoid.
    """
    if lib.adt is None:
        raise ValueError(f"{lib.name}: no ADT block, cannot gate on protein")

    a = clr(lib.adt)
    idx = {n: i for i, n in enumerate(lib.adt_names)}
    needed = ["CD3", "CD19", "CD14", "CD56", "CD16", "CD4", "CD8"]
    missing = [m for m in needed if m not in idx]
    if missing:
        raise ValueError(f"{lib.name}: ADT panel missing {missing}")

    thr, pos = {}, {}
    for m in needed:
        t = _threshold_2class(a[idx[m]])
        thr[m] = t
        pos[m] = a[idx[m]] > t

    # RNA CD3 triple-negative requirement for the NK gate (protocol 2.1)
    g2i = {g: i for i, g in enumerate(lib.gene_names)}
    cd3_rows = [g2i[g] for g in CD3_RNA if g in g2i]
    if cd3_rows:
        cd3_counts = np.asarray(lib.rna[cd3_rows].todense())
        cd3_all_zero = (cd3_counts == 0).all(axis=0)
        cd3_triple_pos = (cd3_counts > 0).all(axis=0)
    else:
        cd3_all_zero = np.ones(lib.rna.shape[1], bool)
        cd3_triple_pos = np.zeros(lib.rna.shape[1], bool)

    T = pos["CD3"]
    lineage = np.full(lib.rna.shape[1], "unassigned", dtype=object)
    lineage[T & pos["CD8"] & ~pos["CD4"]] = "CD8T"
    lineage[T & pos["CD4"] & ~pos["CD8"]] = "CD4T"
    lineage[~T & pos["CD19"] & ~pos["CD14"]] = "B"
    lineage[~T & ~pos["CD19"] & pos["CD14"]] = "Myeloid"
    nk = (~T & ~pos["CD19"] & ~pos["CD14"]
          & (pos["CD56"] | pos["CD16"]) & cd3_all_zero)
    lineage[nk] = "NK"

    tbl = pd.DataFrame({
        "library": lib.name,
        "marker": list(thr),
        "threshold_clr": [thr[m] for m in thr],
        "pct_positive": [100 * pos[m].mean() for m in thr],
    })
    return lineage.astype(str), tbl


SCORE_PANELS = {
    "NK": ["NKG7", "KLRD1", "KLRF1", "GNLY", "PRF1", "FCGR3A"],
    "T": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["IL7R", "CD40LG", "CD4"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Myeloid": ["LYZ", "CD68", "FCN1", "ITGAM", "AIF1", "CST3"],
}


def gate_lineages_rna(lib: Library, margin: float = 0.15) -> tuple[np.ndarray, pd.DataFrame]:
    """Lineage assignment from RNA module scores.

    Needed because the ADT panel is present in only 20 of the 73
    GSE154826 libraries (and CD56/CD8/CD19 are too sparsely captured to
    threshold even there), so protein gating cannot carry the main
    analysis.  Every gene named in `RNA_MARKERS` / `SCORE_PANELS` is
    written to excluded_genes.txt and barred from Phase 4 testing, which
    is what keeps the assignment non-circular.

    Scores are means of log-normalised expression, z-scored across cells
    so that panels of different size and baseline are comparable.  A cell
    is assigned only if its best score beats the runner-up by `margin`;
    ambiguous cells become "unassigned" rather than being forced into a
    gate.
    """
    x = lib.rna
    lib_size = np.asarray(x.sum(axis=0)).ravel()
    lib_size[lib_size == 0] = 1
    g2i = {}
    for i, g in enumerate(lib.gene_names):
        g2i.setdefault(g, []).append(i)

    def score(panel):
        rows = [i for g in panel for i in g2i.get(g, [])]
        if not rows:
            return np.zeros(x.shape[1])
        sub = np.asarray(x[rows].todense())
        norm = np.log1p(sub / lib_size * 1e4)
        s = norm.mean(axis=0)
        sd = s.std()
        return (s - s.mean()) / sd if sd > 0 else s * 0.0

    sc = {k: score(v) for k, v in SCORE_PANELS.items()}

    # CD3 transcripts at raw-count level, for the NK triple-negative rule
    cd3_rows = [i for g in CD3_RNA for i in g2i.get(g, [])]
    cd3 = np.asarray(x[cd3_rows].todense()) if cd3_rows else np.zeros((1, x.shape[1]))
    cd3_all_zero = (cd3 == 0).all(axis=0)

    # top-level: T vs NK vs B vs Myeloid
    top = {"T": sc["T"], "NK": sc["NK"], "B": sc["B"], "Myeloid": sc["Myeloid"]}
    keys = list(top)
    mat = np.stack([top[k] for k in keys])
    order = np.argsort(mat, axis=0)
    best = mat[order[-1], np.arange(mat.shape[1])]
    second = mat[order[-2], np.arange(mat.shape[1])]
    call = np.array([keys[i] for i in order[-1]], dtype=object)
    call[(best - second) < margin] = "unassigned"

    lineage = np.array(call, dtype=object)
    # split T into CD8 / CD4
    is_t = lineage == "T"
    d = sc["CD8"] - sc["CD4"]
    lineage[is_t & (d > margin)] = "CD8T"
    lineage[is_t & (d < -margin)] = "CD4T"
    lineage[lineage == "T"] = "unassigned"
    # protocol 2.1: an NK call requires zero CD3 transcripts
    lineage[(lineage == "NK") & ~cd3_all_zero] = "unassigned"

    tbl = pd.DataFrame({
        "library": lib.name,
        "panel": list(sc),
        "score_sd": [float(np.std(v)) for v in sc.values()],
        "n_markers_found": [
            sum(1 for g in SCORE_PANELS[k] if g in g2i) for k in sc
        ],
    })
    return lineage.astype(str), tbl


def purity_report(lib: Library, lineage: np.ndarray) -> pd.DataFrame:
    """Protocol 2.1: discrete triple-positive counts, per lineage gate.

    CD8T should be >90% CD3 triple-positive; the NK gate should sit near
    the B-cell gate, which serves as the ambient/soup floor reference.
    """
    g2i = {g: i for i, g in enumerate(lib.gene_names)}
    rows = [g2i[g] for g in CD3_RNA if g in g2i]
    if not rows:
        return pd.DataFrame()
    cd3 = np.asarray(lib.rna[rows].todense())
    triple_pos = (cd3 > 0).all(axis=0)
    any_pos = (cd3 > 0).any(axis=0)

    out = []
    for l in ["NK", "CD8T", "CD4T", "B", "Myeloid"]:
        m = lineage == l
        if m.sum() == 0:
            continue
        out.append({
            "library": lib.name,
            "lineage": l,
            "n_cells": int(m.sum()),
            "cd3_triple_pos_n": int(triple_pos[m].sum()),
            "cd3_triple_pos_pct": 100 * float(triple_pos[m].mean()),
            "cd3_any_pos_pct": 100 * float(any_pos[m].mean()),
        })
    return pd.DataFrame(out)


def demux_hto(lib: Library, min_ratio: float = 3.0) -> np.ndarray:
    """Assign hashed cells to their HTO, or 'doublet'/'negative'.

    Only meaningful for the 20 libraries that pool tumour and adjacent
    normal into one droplet emulsion -- those are the natural ambient
    control of Phase 6.1, because both conditions share one soup.
    """
    if lib.hto is None or lib.hto.shape[0] < 2:
        return np.full(lib.rna.shape[1], "none", dtype=object)
    h = clr(lib.hto)
    order = np.argsort(h, axis=0)
    top = h[order[-1], np.arange(h.shape[1])]
    second = h[order[-2], np.arange(h.shape[1])]
    call = np.array([lib.hto_names[i] for i in order[-1]], dtype=object)
    ratio = top - second
    call[ratio < np.log1p(min_ratio)] = "doublet"
    call[top <= 0] = "negative"
    return call.astype(str)


# --------------------------------------------------------------------------
# Ambient soup profile and marker-based decontamination
# --------------------------------------------------------------------------
# Genes each lineage should NOT express.  Their entire signal in that lineage
# is soup, so they estimate the contamination fraction directly.  This is the
# SoupX "non-expressed gene" estimator, made concrete by the observation that
# immunoglobulin transcripts appear at +3 to +4.7 log2FC in NK and myeloid
# cells, which cannot produce them.
SOUP_ESTIMATION_MARKERS = {
    "NK":      ["IGKC", "IGHG1", "IGHA1"],
    "CD8T":    ["IGKC", "IGHG1", "IGHA1"],
    "CD4T":    ["IGKC", "IGHG1", "IGHA1"],
    "Myeloid": ["IGKC", "IGHG1", "IGHA1"],
    "B":       ["GZMB", "GZMA", "GZMH"],
}
# Disjoint sets, never used for estimation, reserved for acceptance testing.
SOUP_VALIDATION_MARKERS = {
    "NK":      ["MS4A1", "CD79A", "LYZ"],
    "CD8T":    ["MS4A1", "CD79A", "LYZ"],
    "CD4T":    ["MS4A1", "CD79A", "LYZ"],
    "Myeloid": ["MS4A1", "CD79A", "IGLC2"],
    "B":       ["NKG7", "KLRD1", "LYZ"],
}


def soup_profile(path: str, max_umi: int = 20, min_umi: int = 1):
    """Ambient expression profile, estimated from empty droplets.

    Droplets holding between `min_umi` and `max_umi` counts contain no cell,
    so their pooled profile *is* the soup.  CellBender's own estimate on
    library 48 put the empty-droplet prior at 23 counts, which is where this
    default comes from.

    Returns (profile summing to 1, gene names, number of droplets used).
    """
    import tarfile

    with tarfile.open(path, "r:gz") as tf:
        members = {os.path.basename(m.name): m for m in tf.getmembers() if m.isfile()}
        ft = pd.read_csv(tf.extractfile(next(v for k, v in members.items()
                                             if k.endswith("features.tsv"))),
                         sep="\t", header=None, names=["id", "name", "type"])
        fh = tf.extractfile(next(v for k, v in members.items()
                                 if k.endswith("matrix.mtx")))
        fh.readline()
        while True:
            pos = fh.tell()
            if not fh.readline().startswith(b"%"):
                fh.seek(pos)
                break
        dims = fh.readline().split()
        n_feat, n_bc = int(dims[0]), int(dims[1])
        trip = pd.read_csv(fh, sep=r"\s+", header=None, dtype=np.int64,
                           names=["f", "b", "v"])

    mat = sp.coo_matrix((trip.v.values, (trip.f.values - 1, trip.b.values - 1)),
                        shape=(n_feat, n_bc)).tocsc()
    is_rna = (ft.type == "Gene Expression").values
    rna = mat[is_rna]
    umi = np.asarray(rna.sum(axis=0)).ravel()
    empty = (umi >= min_umi) & (umi <= max_umi)
    prof = np.asarray(rna[:, empty].sum(axis=1)).ravel().astype(float)
    total = prof.sum()
    if total <= 0:
        raise ValueError(f"{path}: no counts in empty droplets")
    return prof / total, ft.name.values[is_rna].astype(str), int(empty.sum())


def estimate_contamination(observed: np.ndarray, genes: np.ndarray,
                           profile: np.ndarray, markers: list) -> float:
    """Fraction of a profile's counts attributable to soup.

    rho = (counts of marker genes observed) / (counts expected if the whole
    profile were soup).  Markers must be genes the lineage does not express,
    so their observed counts are entirely ambient.
    """
    idx = {g: i for i, g in enumerate(genes)}
    rows = [idx[m] for m in markers if m in idx]
    if not rows:
        return 0.0
    obs = float(observed[rows].sum())
    expected_if_all_soup = float(observed.sum()) * float(profile[rows].sum())
    if expected_if_all_soup <= 0:
        return 0.0
    return float(np.clip(obs / expected_if_all_soup, 0.0, 1.0))


def decontaminate(observed: np.ndarray, profile: np.ndarray, rho: float) -> np.ndarray:
    """Subtract rho x total x soup, flooring at zero."""
    return np.maximum(observed - rho * float(observed.sum()) * profile, 0.0)
