"""Source-format readers.

Branching here is on FILE FORMAT, never on dataset identity. Adding a dataset
whose format is already supported must require no change to this file.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_dense_tsv_genes_by_cells(path, chunk_rows=256, dtype=np.int32,
                                  progress=None):
    """Read a dense genes x cells TSV into a sparse cells x genes CSR matrix.

    Row labels become gene ids, the header row becomes cell ids. Read in row
    chunks so peak memory stays proportional to chunk_rows x n_cells.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    n_cells = len(cell_ids)

    gene_ids = []
    blocks = []
    reader = pd.read_csv(path, sep="\t", header=0, index_col=0,
                         chunksize=chunk_rows, dtype=dtype, engine="c")
    seen = 0
    for chunk in reader:
        gene_ids.extend(chunk.index.astype(str).tolist())
        # chunk is genes x cells -> keep as sparse, transpose at the end
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=dtype)))
        seen += chunk.shape[0]
        if progress and seen % (chunk_rows * 20) == 0:
            progress(seen)
    X_gc = sp.vstack(blocks, format="csr")
    del blocks
    if X_gc.shape[1] != n_cells:
        raise ValueError(f"header has {n_cells} cells but matrix has {X_gc.shape[1]} columns")
    X = X_gc.T.tocsr()
    return X, np.array(cell_ids, dtype=object), np.array(gene_ids, dtype=object)


def read_mtx_dir(path):
    """10x-style directory: matrix.mtx(.gz) + barcodes + features/genes."""
    from scipy.io import mmread
    files = os.listdir(path)

    def pick(*cands):
        for c in cands:
            for f in files:
                if f == c or f == c + ".gz":
                    return os.path.join(path, f)
        raise FileNotFoundError(f"none of {cands} in {path}")

    mtx = pick("matrix.mtx")
    bc = pick("barcodes.tsv")
    ft = pick("features.tsv", "genes.tsv")
    X = sp.csr_matrix(mmread(mtx))
    barcodes = pd.read_csv(bc, sep="\t", header=None)[0].astype(str).to_numpy()
    feat = pd.read_csv(ft, sep="\t", header=None)
    gene_ids = feat[0].astype(str).to_numpy()
    # 10x mtx is genes x cells
    if X.shape[0] == gene_ids.size and X.shape[1] == barcodes.size:
        X = X.T.tocsr()
    return X, barcodes, gene_ids


def read_h5ad_counts(path, counts_layer="X"):
    import anndata as ad
    a = ad.read_h5ad(path)
    X = a.X if counts_layer in ("X", None) else a.layers[counts_layer]
    return sp.csr_matrix(X), a.obs_names.to_numpy(), a.var_names.to_numpy(), a


def read_source(fmt, path, counts_layer="X", progress=None):
    if fmt == "dense_tsv_genes_by_cells":
        X, cells, genes = read_dense_tsv_genes_by_cells(path, progress=progress)
        return X, cells, genes, None
    if fmt == "mtx_dir":
        X, cells, genes = read_mtx_dir(path)
        return X, cells, genes, None
    if fmt == "h5ad":
        return read_h5ad_counts(path, counts_layer)
    raise ValueError(f"unsupported file_format {fmt!r}")


def read_obs_table(path, index_col=0, sep="\t"):
    return pd.read_csv(path, sep=sep, index_col=index_col)
