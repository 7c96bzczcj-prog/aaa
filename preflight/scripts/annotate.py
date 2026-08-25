#!/usr/bin/env python3
"""Per-donor clustering and cluster-level lineage annotation.

Cell-level winner-take-all on z-scored panel scores was tried first and
rejected: in PBMC it called 5.6% HSPC and 7.1% plasma cells, because a panel
with almost no expression still produces z-scores and can win by noise.
Clusters carry enough cells for the panel means to be meaningful, so the call
is made once per cluster and inherited by its cells.

Clustering is per DONOR, over all that donor's libraries at once, so a donor's
compartments (where both are present) share one cluster definition.  The
ambient profile stays strictly per LIBRARY, because a soup belongs to an
emulsion.

Writes one npz per donor: cell barcode, library, cluster, lineage.
"""
from __future__ import annotations

import argparse, json, os, sys, warnings

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import pf_core as pf  # noqa: E402

warnings.filterwarnings("ignore")

LEIDEN_RES = 1.0
N_PCS = 30
N_NEIGHBORS = 15
HVG_N = 2000
# a cluster may only take a lineage label if the winning panel is actually
# expressed: mean log1p(CP10K) over the panel must clear this floor
PANEL_FLOOR = 0.15
# and it must beat the runner-up panel by this much (mean-score units)
CLUSTER_MARGIN = 0.05


def build(libs):
    """libs: list of (library_name, Matrix). Returns AnnData of called cells."""
    import anndata as ad
    Xs, obs_lib, obs_bc = [], [], []
    genes = None
    for name, m in libs:
        cells = pf.call_cells(m)
        if genes is None:
            genes = m.genes
        elif list(genes) != list(m.genes):
            raise ValueError("gene order differs between libraries")
        Xs.append(m.mat[:, cells].T.tocsr())
        obs_lib += [name] * len(cells)
        obs_bc += list(m.barcodes[cells])
    X = sp.vstack(Xs).tocsr()
    a = ad.AnnData(X)
    a.var_names = [f"{g}#{i}" for i, g in enumerate(genes)]
    a.var["gene"] = list(genes)
    a.obs["library"] = obs_lib
    a.obs["barcode"] = obs_bc
    return a


def cluster(a):
    import scanpy as sc
    a.layers["counts"] = a.X.copy()
    # Doublet detection, per library, on raw counts.  Needed because a handful
    # of NK-B and NK-myeloid doublets dominate the low-count foreign genes that
    # the ambient acceptance test is measured on (see lib/pf_core.py).  Scrublet
    # simulates its doublets from each library's own cells, so the procedure is
    # symmetric between compartments.
    try:
        sc.pp.scrublet(a, batch_key="library")
        a.obs["is_doublet"] = a.obs["predicted_doublet"].astype(bool).values
    except Exception as e:  # never let doublet calling silently pass
        print(f"[warn] scrublet failed: {e}", flush=True)
        a.obs["is_doublet"] = False
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=HVG_N)
    a.raw = a
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=min(N_PCS, min(b.shape) - 1), svd_solver="arpack")
    sc.pp.neighbors(b, n_neighbors=N_NEIGHBORS, n_pcs=min(N_PCS, b.obsm["X_pca"].shape[1]))
    sc.tl.leiden(b, resolution=LEIDEN_RES, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    a.obs["leiden"] = b.obs["leiden"].values
    return a


def annotate(a):
    """Cluster-level lineage call from panel means of log1p(CP10K)."""
    gene = np.asarray(a.var["gene"])
    g2i = {}
    for i, g in enumerate(gene):
        g2i.setdefault(g, []).append(i)
    Xn = a.raw.X if a.raw is not None else a.X          # log1p CP10K, all genes
    clusters = np.asarray(a.obs["leiden"])
    uniq = sorted(set(clusters), key=lambda s: int(s))
    panel_mean = {}
    for k, panel in pf.PANELS.items():
        rows = [i for g in panel for i in g2i.get(g, [])]
        if not rows:
            panel_mean[k] = np.zeros(len(uniq))
            continue
        v = np.asarray(Xn[:, rows].mean(axis=1)).ravel()
        panel_mean[k] = np.array([v[clusters == c].mean() for c in uniq])
    keys = list(pf.PANELS)
    M = np.stack([panel_mean[k] for k in keys])          # panels x clusters
    # z-score each panel ACROSS CLUSTERS so panels of different baseline compare
    Z = (M - M.mean(axis=1, keepdims=True)) / np.where(M.std(axis=1, keepdims=True) > 0,
                                                       M.std(axis=1, keepdims=True), 1)
    lab = {}
    detail = []
    for j, c in enumerate(uniq):
        order = np.argsort(Z[:, j])
        top, second = keys[order[-1]], keys[order[-2]]
        gap = Z[order[-1], j] - Z[order[-2], j]
        expressed = M[order[-1], j] >= PANEL_FLOOR
        call = top if (gap >= CLUSTER_MARGIN and expressed) else "unassigned"
        # NK-specific hard rule.  The NK panel (NKG7/PRF1/CTSW/GZMA...) is
        # shared with cytotoxic CD8 T, and on BM donor 1 a plainly T cluster
        # (mean T panel 1.62 vs NK 1.19) won the NK z-argmax.  A cluster may
        # therefore only be called NK if the NK panel is also the argmax of the
        # RAW panel means.  If it is not, the cluster takes the raw argmax
        # label instead.  Applied identically in both compartments.
        raw_top = keys[int(np.argmax(M[:, j]))]
        if call == "NK" and raw_top != "NK":
            call = raw_top
        lab[c] = call
        detail.append({"cluster": c, "n": int((clusters == c).sum()), "call": call,
                       "top": top, "second": second, "gap": float(gap),
                       "panel_mean_top": float(M[order[-1], j]),
                       "raw_argmax": keys[int(np.argmax(M[:, j]))],
                       **{f"m_{k}": float(panel_mean[k][j]) for k in keys}})
    lineage = np.array([lab[c] for c in clusters], dtype=object)

    # --- purity rules inside NK clusters -----------------------------------
    # The inherited protocol-2.1 rule ("an NK call requires ZERO CD3/TCR
    # counts") turns out to be depth- and ambient-dependent: on HCA donor 1 it
    # removed 46.5% of NK-cluster cells in bone marrow but 73.7% in blood,
    # because blood soup is T-dominated and blood NK are sequenced deeper, so
    # both raise the chance of >=1 stray TCR count.  A gate whose stringency is
    # a function of the compartment's ambient composition would launder the
    # very effect task A measures into the gate.  The primary rule is therefore
    # depth-symmetric: keep a cell if its NK panel outscores its T panel on
    # log1p(CP10K).  The zero rule is kept as a labelled sensitivity arm.
    # This choice was made before any Delta_ambient was computed; the only
    # number that drove it is the drop-rate diagnostic quoted above.
    def panel_mean_cells(panel):
        rows = [i for g in panel for i in g2i.get(g, [])]
        if not rows:
            return np.zeros(Xn.shape[0])
        return np.asarray(Xn[:, rows].mean(axis=1)).ravel()

    doublet = np.asarray(a.obs.get("is_doublet", np.zeros(len(lineage), bool))).astype(bool)
    lineage[doublet] = "doublet"

    nk_s = panel_mean_cells(pf.NK_PANEL)
    t_s = panel_mean_cells(pf.T_PANEL)
    lin_ratio = lineage.copy()
    lin_ratio[(lineage == "NK") & (t_s >= nk_s)] = "unassigned"

    rows = [i for g in pf.CD3_ZERO for i in g2i.get(g, [])]
    lin_zero = lineage.copy()
    if rows:
        has = np.asarray(a.layers["counts"][:, rows].sum(axis=1)).ravel() > 0
        lin_zero[(lineage == "NK") & has] = "unassigned"

    return (lin_ratio.astype(str), lineage.astype(str), lin_zero.astype(str), detail)



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True)
    ap.add_argument("--files", required=True, help="JSON list of [name, path, kind]")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    spec = json.loads(args.files)
    libs = []
    for name, path, kind in spec:
        if kind == "h5":
            libs.append((name, pf.read_10x_h5(path)))
        else:
            libs.append((name, pf.read_10x_mtx(*path.split("|"), name)))
    a = build(libs)
    a = cluster(a)
    lineage, lin_cluster, lin_zero, detail = annotate(a)
    np.savez_compressed(args.out, barcode=np.asarray(a.obs["barcode"]),
                        library=np.asarray(a.obs["library"]),
                        leiden=np.asarray(a.obs["leiden"]).astype(str),
                        lineage=lineage, lineage_cluster=lin_cluster,
                        lineage_cd3zero=lin_zero,
                        is_doublet=np.asarray(a.obs["is_doublet"]).astype(bool))
    json.dump(detail, open(args.out.replace(".npz", "_clusters.json"), "w"), indent=1)
    import collections
    c = collections.Counter(lineage)
    nkc = int((lin_cluster == "NK").sum())
    ndbl = int(np.asarray(a.obs["is_doublet"]).astype(bool).sum())
    print(f"[{args.donor}] cells={len(lineage)} doublets={ndbl} NK_cluster={nkc} "
          f"NK_ratio={int((lineage == 'NK').sum())} "
          f"NK_cd3zero={int((lin_zero == 'NK').sum())} | " +
          " ".join(f"{k}={v}" for k, v in c.most_common()), flush=True)


if __name__ == "__main__":
    main()
