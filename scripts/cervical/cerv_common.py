"""
Shared machinery for the cervical-cancer NK-proportion recount.

Design constraints imposed by the task (see docs/CERVICAL_NK_METHODS.md):
  - author annotations are never used; the T/NK compartment is re-clustered from counts
  - NK calls are made on BOTH positive and negative markers
  - ambient RNA is corrected before any gating
  - doublets are removed before any gating
  - every number is produced per donor, never on a pooled object
"""

from __future__ import annotations

import gzip
import os
import re

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

DATA_ROOT = "/home/user/cervical_work"

# ---------------------------------------------------------------- marker sets

# NK positive markers (task step 2)
NK_POS = ["NCAM1", "KLRD1", "KLRF1", "NKG7", "GNLY", "PRF1"]

# T markers an NK cell must be negative for.
#
# This panel is the OUTPUT of scripts/cervical/validate_markers.py, not an
# assumption. The task asks to verify TRBC1 in the data before relying on it;
# applying that same test to every candidate shows the proposed panel
# {CD3D, CD3E, CD3G, TRBC2} does not survive contact with these libraries.
# Inside a CD56+ cytotoxic cluster, measured as a fraction of the T-cell level:
#
#     CD3D 0.06   CD3G 0.05   CD6 0.08     <- usable
#     CD3E 0.64   TRBC2 0.72   TRAC 0.79   TRBC1 0.90   <- do NOT go negative
#
# The failures are the highest-abundance T transcripts, which retain the largest
# ambient residue, and CD3E/CD3zeta machinery that NK cells genuinely transcribe.
# Candidates were also required to be positive in >=80% of T clusters: CD5,
# CD28, CD40LG, THEMIS and IL7R are low in CD8 effector/TEMRA cells, so they
# would let an effector-T cluster collect "negative" votes and be miscalled NK.
T_NEG = ["CD3D", "CD3G", "CD6"]
T_NEG_REJECTED = ["CD3E", "TRBC2", "TRAC", "TRBC1", "LCK", "SKAP1", "ITK", "TRAT1"]

# a cluster counts as negative for a marker below this fraction of the T level
NEG_MAX_FRAC_OF_T = 0.25
# and the T level is taken as this percentile of per-cluster means
T_REF_PCTL = 90

# compartment markers used only to split major lineages before the T/NK subset
LINEAGE = {
    "immune": ["PTPRC"],
    "T_NK": ["CD3D", "CD3E", "CD2", "TRAC", "NKG7", "KLRD1", "IL7R", "CD7"],
    "B": ["CD79A", "CD79B", "MS4A1", "BANK1"],
    "plasma": ["JCHAIN", "MZB1", "DERL3", "XBP1"],
    "myeloid": ["LYZ", "CD68", "CD14", "FCGR3A", "AIF1", "ITGAX", "C1QA"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "pDC": ["LILRA4", "IRF7", "CLEC4C", "GZMB"],
    "epithelial": ["EPCAM", "KRT5", "KRT13", "KRT14", "KRT17", "KRT18", "KRT19", "SFN"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "muscle": ["ACTA2", "MYH11", "DES", "TAGLN"],
}

# lineages counted as CD45+ immune, and as lymphocytes, for the two extra denominators
IMMUNE_LINEAGES = {"T_NK", "NK", "B", "plasma", "myeloid", "mast", "pDC"}
LYMPHOID_LINEAGES = {"T_NK", "NK", "B", "plasma"}


# ------------------------------------------------------------ sample metadata

def emtab_samples() -> pd.DataFrame:
    """
    E-MTAB-12305 (Li & Hua 2022; accession corrected by the 2024 corrigendum,
    Front Immunol 15:1386072).

    Donor assignment is read off the SDRF, not guessed: Characteristics[age]
    ties N1/T1 (45y), N2/T2 (50y), N3/T3 (51y) into three paired
    adjacent-normal/tumour donors, and T4/L1 (48y) into one tumour + its own
    metastatic node. Those pairs are the only strictly valid comparisons here.
    """
    rows = [
        ("N1", "normal_adj", "P45"), ("N2", "normal_adj", "P50"), ("N3", "normal_adj", "P51"),
        ("H1", "HSIL", "P38"), ("H2", "HSIL", "P42"),
        ("T1", "tumor", "P45"), ("T2", "tumor", "P50"), ("T3", "tumor", "P51"),
        ("T4", "tumor", "P48"), ("L1", "lymph_node", "P48"),
    ]
    return pd.DataFrame(
        [{"dataset": "E-MTAB-12305", "sample": s, "tissue": t,
          "donor": d, "histology": "SCC" if t in ("tumor", "lymph_node") else "-",
          "prefix": os.path.join(DATA_ROOT, "E-MTAB-12305", f"{s}_"),
          "mtx": "matrix.mtx.gz", "bc": "barcodes.tsv.gz", "ft": "features.tsv.gz"}
         for s, t, d in rows])


def _geo_samples(gse: str, spec: list[tuple[str, str, str, str]], sep: str) -> pd.DataFrame:
    """spec entries: (gsm_prefix, sample_label, tissue, histology)."""
    out = []
    for gsm, label, tissue, hist in spec:
        out.append({
            "dataset": gse, "sample": label, "tissue": tissue, "donor": label,
            "histology": hist,
            "prefix": os.path.join(DATA_ROOT, gse, f"{gsm}{sep}"),
            "mtx": "matrix.mtx.gz", "bc": "barcodes.tsv.gz",
            "ft": "genes.tsv.gz" if gse == "GSE173231" else "features.tsv.gz",
        })
    return pd.DataFrame(out)


def gse208653_samples() -> pd.DataFrame:
    """
    GSE208653 (PMID 36967539). The GEO sample *titles* carry HPV status that the
    supplementary file names hide: GSM6360682/683 are named "N_1/N_2" on disk but
    titled "N_HPV" in GEO, i.e. histologically normal cervix that is HPV-infected,
    not HPV-negative. Both are normal tissue; the distinction is kept in `hpv`.
    """
    spec = [
        ("GSM6360680_N_HPV_NEG_1", "N_HPVneg_1", "normal", "-"),
        ("GSM6360681_N_HPV_NEG_2", "N_HPVneg_2", "normal", "-"),
        ("GSM6360682_N_1", "N_HPVpos_1", "normal", "-"),
        ("GSM6360683_N_2", "N_HPVpos_2", "normal", "-"),
        ("GSM6360684_HSIL_1", "HSIL_1", "HSIL", "-"),
        ("GSM6360685_HSIL_2", "HSIL_2", "HSIL", "-"),
        ("GSM6360686_SCC_4", "SCC_4", "tumor", "SCC"),
        ("GSM6360687_SCC_5", "SCC_5", "tumor", "SCC"),
        ("GSM6360688_ADC_6", "ADC_6", "tumor", "ADC"),
    ]
    df = _geo_samples("GSE208653", spec, ".")
    df["hpv"] = ["neg", "neg", "pos", "pos", "pos", "pos", "pos", "pos", "pos"]
    return df


def gse197461_samples() -> pd.DataFrame:
    spec = [
        ("GSM5917937_SCC_1", "SCC_1", "tumor", "SCC"),
        ("GSM5917938_SCC_2", "SCC_2", "tumor", "SCC"),
        ("GSM5917939_SCC_3", "SCC_3", "tumor", "SCC"),
        ("GSM5917940_ADC_1", "ADC_1", "tumor", "ADC"),
        ("GSM5917941_ADC_2", "ADC_2", "tumor", "ADC"),
        ("GSM5917942_ADC_3", "ADC_3", "tumor", "ADC"),
        ("GSM5917943_ADC_4", "ADC_4", "tumor", "ADC"),
        ("GSM5917944_ADC_5", "ADC_5", "tumor", "ADC"),
    ]
    return _geo_samples("GSE197461", spec, "_")


def gse173231_samples() -> pd.DataFrame:
    """Healthy ectocervix biopsies (gene-expression libraries only; the *_2 GSMs are TCR)."""
    spec = [
        ("GSM5263982_CX1_1", "CX1", "normal_healthy", "-"),
        ("GSM5263984_CX2_1", "CX2", "normal_healthy", "-"),
        ("GSM5263986_CX3_1", "CX3", "normal_healthy", "-"),
        ("GSM5263988_CX4_1", "CX4", "normal_healthy", "-"),
        ("GSM5263990_CX5_1", "CX5", "normal_healthy", "-"),
        ("GSM5263992_CX6A_1", "CX6A", "normal_healthy", "-"),
        ("GSM5263994_CX6B_1", "CX6B", "normal_healthy", "-"),
        ("GSM5263996_CX7A_1", "CX7A", "normal_healthy", "-"),
        ("GSM5263998_CX7B_1", "CX7B", "normal_healthy", "-"),
        ("GSM5264000_CX8_1", "CX8", "normal_healthy", "-"),
    ]
    df = _geo_samples("GSE173231", spec, "_")
    # CX6A/CX6B and CX7A/CX7B are two biopsies from one subject -> one donor
    df["donor"] = df["sample"].str.replace(r"([67])[AB]$", r"\1", regex=True)
    return df


def all_samples() -> pd.DataFrame:
    df = pd.concat(
        [emtab_samples(), gse208653_samples(), gse197461_samples(), gse173231_samples()],
        ignore_index=True,
    )
    if "hpv" not in df:
        df["hpv"] = "-"
    return df.fillna({"hpv": "-"})


# ------------------------------------------------------------------- loading

def read_10x(prefix: str, mtx: str, bc: str, ft: str):
    """Read a 10x triplet. Returns (counts CSR cells x genes, barcodes, gene symbols)."""
    m = scipy.io.mmread(prefix + mtx)          # genes x cells
    m.data = m.data.astype(np.float32)          # halve before the CSR conversion
    X = sp.csr_matrix(m.T)                      # -> cells x genes
    del m
    with gzip.open(prefix + bc, "rt") as fh:
        barcodes = [ln.strip() for ln in fh if ln.strip()]
    genes = []
    with gzip.open(prefix + ft, "rt") as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    assert X.shape == (len(barcodes), len(genes)), (X.shape, len(barcodes), len(genes))
    X.data = X.data.astype(np.float32)
    return X, np.array(barcodes), np.array(genes)


def collapse_duplicate_genes(X, genes):
    """Sum columns that share a gene symbol so every dataset has unique symbols."""
    uniq, inv = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        # already unique: reordering would only cost a full sparse copy
        return X, genes
    M = sp.csr_matrix(
        (np.ones(len(genes), dtype=np.float32), (np.arange(len(genes)), inv)),
        shape=(len(genes), len(uniq)),
    )
    return sp.csr_matrix(X @ M), uniq


# --------------------------------------------------------------- ambient RNA

def decontx(X, clusters, max_iter=150, delta=(10.0, 10.0), tol=1e-4, verbose=False):
    """
    decontX-style ambient-RNA correction (Yang et al., Genome Biology 2020).

    SoupX proper needs the raw (empty-droplet) matrix to build the soup profile.
    Every cervical dataset here is distributed as cellranger *filtered* matrices
    only, so the empty droplets do not exist to be read. decontX is the standard
    equivalent that estimates the contamination distribution from the cell
    matrix itself, using the other clusters as the ambient pool.

    Model, for cell i in cluster k:
        x_ig ~ (1 - theta_i) * phi_kg  +  theta_i * eta_kg
    where phi_k is the native profile of cluster k and eta_k is the ambient
    profile seen by cluster k (all *other* clusters, size-weighted).

    Returns (theta per cell, decontaminated counts CSR float32).
    """
    X = sp.csr_matrix(X)
    n_cells, n_genes = X.shape
    clusters = np.asarray(clusters)
    ks = np.unique(clusters)
    idx_by_k = {k: np.where(clusters == k)[0] for k in ks}

    n_i = np.asarray(X.sum(axis=1)).ravel()
    n_i[n_i == 0] = 1.0

    # cluster count profiles
    def cluster_sums(mat):
        out = np.zeros((len(ks), n_genes), dtype=np.float64)
        for j, k in enumerate(ks):
            out[j] = np.asarray(mat[idx_by_k[k]].sum(axis=0)).ravel()
        return out

    csum = cluster_sums(X)
    eps = 1e-10

    def normed(a):
        s = a.sum(axis=1, keepdims=True)
        s[s == 0] = 1.0
        return a / s

    phi = normed(csum + eps)
    theta = np.full(n_cells, 0.5)

    # pre-split the matrix by cluster once (COO for per-nonzero arithmetic)
    sub_coo, sub_rows = {}, {}
    for j, k in enumerate(ks):
        rows = idx_by_k[k]
        sub_rows[k] = rows
        sub_coo[k] = sp.coo_matrix(X[rows])

    prev = theta.copy()
    for it in range(max_iter):
        # ambient profile for cluster k = size-weighted mix of the OTHER clusters
        tot = csum.sum(axis=0)
        eta = np.empty_like(phi)
        for j, k in enumerate(ks):
            eta[j] = tot - csum[j]
        eta = normed(eta + eps)

        new_csum = np.zeros_like(csum)
        contam_tot = np.zeros(n_cells)

        t_over = theta / np.clip(1.0 - theta, 1e-8, None)   # theta/(1-theta)

        for j, k in enumerate(ks):
            coo = sub_coo[k]
            rows_global = sub_rows[k][coo.row]
            s_kg = eta[j] / np.clip(phi[j], eps, None)       # ambient/native ratio
            # P(native) for each nonzero = 1 / (1 + t_i * s_kg)
            r = 1.0 / (1.0 + t_over[rows_global] * s_kg[coo.col])
            native = coo.data * r
            contam = coo.data - native
            # bincount, not np.add.at: identical arithmetic, but add.at is
            # unbuffered and dominates the runtime on 30M+ nonzero matrices
            new_csum[j] = np.bincount(coo.col, weights=native, minlength=n_genes)
            contam_tot += np.bincount(rows_global, weights=contam, minlength=n_cells)

        phi = normed(new_csum + eps)
        csum = new_csum
        # Beta prior on theta, as decontX does, to keep it off the boundary
        theta = (contam_tot + delta[0]) / (n_i + delta[0] + delta[1])

        shift = np.abs(theta - prev).max()
        if verbose and (it % 25 == 0 or shift < tol):
            print(f"      decontX iter {it:3d} max|dtheta|={shift:.2e} "
                  f"median theta={np.median(theta):.3f}", flush=True)
        if shift < tol:
            break
        prev = theta.copy()

    # final E-step -> decontaminated counts
    t_over = theta / np.clip(1.0 - theta, 1e-8, None)
    tot = csum.sum(axis=0)
    eta = normed(np.stack([tot - csum[j] for j in range(len(ks))]) + eps)
    rows_all, cols_all, vals_all = [], [], []
    for j, k in enumerate(ks):
        coo = sub_coo[k]
        rows_global = sub_rows[k][coo.row]
        s_kg = eta[j] / np.clip(phi[j], eps, None)
        r = 1.0 / (1.0 + t_over[rows_global] * s_kg[coo.col])
        rows_all.append(rows_global)
        cols_all.append(coo.col)
        vals_all.append(coo.data * r)

    Xd = sp.coo_matrix(
        (np.concatenate(vals_all).astype(np.float32),
         (np.concatenate(rows_all), np.concatenate(cols_all))),
        shape=X.shape,
    ).tocsr()
    Xd.data = np.round(Xd.data)
    Xd.eliminate_zeros()
    return theta, Xd


# ------------------------------------------------------------------- gating

def harmony_embed(sub, batch_key, n_iter=20):
    """
    Batch-correct the PCA embedding across libraries.

    scanpy's `external.pp.harmony_integrate` transposes harmonypy's output, which
    is correct for harmonypy 0.0.x but wrong for 2.x, where `Z_corr` is already
    (cells, PCs). That mismatch raises rather than corrupting silently, but it
    means the wrapper is unusable here; harmonypy is called directly and the
    orientation is checked instead of assumed.
    """
    import harmonypy

    Xp = np.asarray(sub.obsm["X_pca"])
    ho = harmonypy.run_harmony(Xp, sub.obs, [batch_key], max_iter_harmony=n_iter)
    Z = np.asarray(ho.Z_corr)
    if Z.shape != Xp.shape:
        Z = Z.T
    if Z.shape != Xp.shape:
        raise RuntimeError(f"harmony returned {Z.shape}, expected {Xp.shape}")
    sub.obsm["X_pca_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)
    return "X_pca_harmony"


def refine_nk_subclusters(tnk, nk_mask, receptor_frac=0.25, resolution=0.6,
                          min_cells=60, seed=0):
    """
    Second gating stage: re-cluster the candidate NK cells and drop subclusters
    that carry no NK receptor.

    Cluster-mean marker levels are partly a product of the clustering itself, so
    a single pass can hand back a "CD56+ cluster" that is internally mixed.
    Re-clustering the candidates exposes that: in the tumour-containing datasets
    a large subpopulation appears with NCAM1 and KLRF1 essentially absent but
    TRAC high, which is not NK under any receptor definition. It survives the
    first pass only because the parent cluster's mean is carried by the genuine
    NK cells sitting beside it.

    A subcluster is kept when NCAM1 or KLRF1 reaches `receptor_frac` of the best
    NK subcluster's level, i.e. the receptor signal has to be comparable to the
    main NK population rather than merely above the T-cell floor.

    Returns (refined mask, per-subcluster table).
    """
    import scanpy as _sc

    if int(nk_mask.sum()) < min_cells:
        return nk_mask, None

    sub = tnk[nk_mask].copy()
    sub.X = sub.layers["lognorm"].copy()
    _sc.pp.highly_variable_genes(sub, n_top_genes=min(1500, sub.n_vars - 1))
    s2 = sub[:, sub.var.highly_variable].copy()
    _sc.pp.scale(s2, max_value=10)
    _sc.tl.pca(s2, n_comps=min(20, s2.n_obs - 1, s2.n_vars - 1), svd_solver="arpack")
    _sc.pp.neighbors(s2, n_neighbors=min(15, s2.n_obs - 1))
    _sc.tl.leiden(s2, resolution=resolution, key_added="nksub",
                  flavor="igraph", n_iterations=2, directed=False, random_state=seed)
    sub.obs["nksub"] = s2.obs["nksub"].values

    ex = cluster_expression(sub, "nksub",
                            ["NCAM1", "KLRF1", "KLRD1", "NKG7", "GNLY", "PRF1",
                             "CD3D", "CD3G", "CD6", "CD3E", "TRAC", "FCGR3A", "IL7R"])
    ex["n"] = sub.obs.groupby("nksub", observed=True).size()
    keep = np.zeros(len(ex), dtype=bool)
    for g in ("NCAM1", "KLRF1"):
        m = ex[g].values.max()
        if m > 1e-9:
            keep |= ex[g].values >= receptor_frac * m
    ex["kept"] = keep

    kept_ids = set(ex.index[keep])
    sub_keys = np.asarray(sub.obs["nksub"].values).astype(str)
    refined = nk_mask.copy()
    idx = np.where(nk_mask)[0]
    refined[idx[~np.isin(sub_keys, list(kept_ids))]] = False
    return refined, ex


def detection_rates(adata, genes, layer="counts"):
    """Fraction of cells with >0 counts, per gene, as a DataFrame column set."""
    import scipy.sparse as _sp
    M = adata.layers[layer] if layer else adata.X
    out = {}
    for g in genes:
        if g not in adata.var_names:
            out[g] = np.zeros(adata.n_obs, dtype=bool)
            continue
        col = M[:, adata.var_names.get_loc(g)]
        col = col.toarray().ravel() if _sp.issparse(col) else np.asarray(col).ravel()
        out[g] = col > 0
    return pd.DataFrame(out, index=adata.obs_names)


def cluster_marker_table(adata, cluster_key, genes, layer="counts"):
    """Per-cluster detection rate for each gene."""
    det = detection_rates(adata, genes, layer=layer)
    det[cluster_key] = adata.obs[cluster_key].values
    return det.groupby(cluster_key, observed=True).mean()


def cluster_expression(adata, cluster_key, genes, layer="lognorm"):
    """Per-cluster mean expression for each gene."""
    M = adata.layers[layer] if layer else adata.X
    keys = np.asarray(adata.obs[cluster_key].values).astype(str)
    out = {}
    for g in genes:
        if g not in adata.var_names:
            out[g] = np.zeros(adata.n_obs)
            continue
        col = M[:, adata.var_names.get_loc(g)]
        out[g] = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
    df = pd.DataFrame(out)
    df["_k"] = keys
    return df.groupby("_k", observed=True).mean()


def call_nk_clusters(tnk, cluster_key="tnk_cluster",
                     neg_frac=NEG_MAX_FRAC_OF_T, min_neg=None,
                     receptor_frac=0.25):
    """
    Decide which T/NK clusters are NK.

    Negativity is judged *relative to the T-cell level in the same dataset*
    rather than against an absolute detection rate. Detection rate is strongly
    depth-dependent -- in a deeply sequenced library a single ambient molecule
    counts as "detected" -- so an absolute cutoff silently tightens as depth
    rises. The T reference is the 90th percentile of per-cluster means, i.e. the
    level in a genuine T cluster.

    Returns (per-cluster table, set of NK cluster ids).
    """
    genes = sorted(set(NK_POS + T_NEG + ["FCGR3A", "CD2", "IL7R", "TRAC", "CD3E"]))
    ex = cluster_expression(tnk, cluster_key, genes)
    min_neg = len(T_NEG) if min_neg is None else min_neg

    tab = ex.copy()
    n_neg = np.zeros(len(ex), dtype=int)
    for g in T_NEG:
        ref = float(np.percentile(ex[g].values, T_REF_PCTL))
        neg = (ex[g].values < neg_frac * ref) if ref > 1e-9 else np.ones(len(ex), bool)
        tab[f"neg_{g}"] = neg
        n_neg += neg.astype(int)
    tab["n_T_neg"] = n_neg

    # cytotoxic-granule requirement keeps ILC1 (which lack them) out of the count
    cyto = ((ex["NKG7"].values > np.median(ex["NKG7"].values))
            & (ex["GNLY"].values > np.median(ex["GNLY"].values)))
    # and at least one NK receptor clearly enriched relative to its own maximum
    rec = np.zeros(len(ex), dtype=bool)
    for g in ("NCAM1", "KLRF1", "KLRD1"):
        m = ex[g].values.max()
        if m > 1e-9:
            rec |= ex[g].values >= receptor_frac * m
    tab["cytotoxic"] = cyto
    tab["receptor"] = rec
    tab["is_NK"] = (n_neg >= min_neg) & cyto & rec
    return tab, set(tab.index[tab.is_NK])
