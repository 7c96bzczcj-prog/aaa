"""
Convert GSE131907's single dense UMI table into per-sample 10x triplets.

The matrix ships as one gzipped tab-separated table, genes in rows and all
208,506 cells in columns. Writing it back out per sample lets the lung arm run
through exactly the same build_dataset.py as the cervical arm, rather than a
parallel code path that could drift from it.

Only the primary-site lung samples are kept (tLung tumour, nLung normal); brain
metastases, lymph nodes and pleural effusions are different tissues and would
not belong in a tumour-vs-tumour comparison.
"""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

SRC = "/home/user/lung_work/GSE131907"
DST = "/home/user/lung_work/GSE131907_samples"
KEEP_ORIGIN = {"tLung", "nLung"}
CHUNK = 250          # gene rows per read; ~250 x 208506 int32 is about 200 MB


def main():
    os.makedirs(DST, exist_ok=True)
    ann = pd.read_csv(f"{SRC}/GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    ann = ann[ann.Sample_Origin.isin(KEEP_ORIGIN)]
    print(f"cells to keep: {len(ann)} from {ann.Sample.nunique()} samples")
    print(ann.groupby(["Sample_Origin"]).Sample.nunique().to_string())

    # A plain set, not np.isin: both arrays are object dtype, and np.isin on
    # object arrays degrades to a quadratic scan that does not finish here.
    wanted = set(ann["Index"].astype(str).tolist())

    path = f"{SRC}/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cells = np.array(header[1:], dtype=object)
    print(f"matrix columns: {len(cells)}", flush=True)

    sel = np.array([i for i, c in enumerate(cells) if c in wanted], dtype=np.int64)
    print(f"columns selected: {len(sel)}", flush=True)
    sel_cells = cells[sel]

    # usecols pushes the column filter into the C parser, so the 120k discarded
    # columns are never converted to Python objects at all
    usecols = [0] + (sel + 1).tolist()

    blocks, genes = [], []
    # no global dtype= here: it would be applied to the gene-name index column too
    reader = pd.read_csv(path, sep="\t", header=0, index_col=0, usecols=usecols,
                         chunksize=CHUNK, engine="c")
    n = 0
    for chunk in reader:
        if not blocks:
            # usecols returns columns in file order, not in the order given;
            # sel is ascending so they should agree, but assert rather than trust
            assert list(chunk.columns) == list(sel_cells), "column order mismatch"
        genes.extend(chunk.index.astype(str).tolist())
        blocks.append(sp.csr_matrix(chunk.values.astype(np.float32)))
        n += chunk.shape[0]
        if n % 2500 < CHUNK:
            nnz = sum(b.nnz for b in blocks)
            print(f"  {n} genes read, nnz so far {nnz/1e6:.0f}M", flush=True)

    M = sp.vstack(blocks, format="csr")     # genes x cells
    del blocks
    genes = np.array(genes)
    print(f"assembled {M.shape} nnz={M.nnz/1e6:.0f}M")

    X = sp.csc_matrix(M)                    # column slicing by cell
    del M

    meta = ann.set_index("Index")
    order = pd.Index(sel_cells)
    smp = meta.loc[order, "Sample"].values
    origin = meta.loc[order, "Sample_Origin"].values

    rows = []
    for s in pd.unique(smp):
        idx = np.where(smp == s)[0]
        sub = sp.csc_matrix(X[:, idx])
        pre = os.path.join(DST, f"{s}_")
        scipy.io.mmwrite(pre + "matrix.mtx", sub.tocoo(), field="integer")
        os.system(f"gzip -f {pre}matrix.mtx")
        with gzip.open(pre + "barcodes.tsv.gz", "wt") as fh:
            for b in order[idx]:
                fh.write(b + "\n")
        with gzip.open(pre + "features.tsv.gz", "wt") as fh:
            for g in genes:
                fh.write(f"{g}\t{g}\tGene Expression\n")
        rows.append({"sample": s, "origin": origin[idx][0], "n_cells": len(idx)})
        print(f"  wrote {s}: {len(idx)} cells ({origin[idx][0]})", flush=True)

    pd.DataFrame(rows).to_csv(f"{DST}/samples.csv", index=False)
    print("\ndone")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
