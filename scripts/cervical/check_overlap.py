"""
Are any of these libraries the same libraries?

E-MTAB-12305 (Li & Hua, Fudan), GSE197461 (Qiu ... Hua K, Fudan) and GSE208653
(cites the same PMID 37794698) come out of overlapping author groups. If two
"independent" datasets actually re-deposit the same library, treating them as
replication is circular. Two signals decide it:

  1. cell-barcode overlap. 10x barcodes are drawn from a ~737k whitelist, so two
     unrelated libraries of n1 and n2 cells share about n1*n2/737280 barcodes by
     chance. The same library shares nearly all of them.
  2. pseudobulk correlation on shared genes. Re-deposits correlate ~1.0.

Both must be unremarkable before the datasets are called independent.
"""
from __future__ import annotations

import gzip
import itertools
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

WHITELIST = 737280.0
OUT = "/home/user/aaa/results/cervical"


def read_barcodes(prefix, bc):
    with gzip.open(prefix + bc, "rt") as fh:
        return {ln.strip().split("-")[0] for ln in fh if ln.strip()}


def pseudobulk(prefix, mtx, bc, ft):
    X, barcodes, genes = cc.read_10x(prefix, mtx, bc, ft)
    X, genes = cc.collapse_duplicate_genes(X, genes)
    v = np.asarray(X.sum(axis=0)).ravel()
    return pd.Series(v, index=genes)


def main():
    reg = cc.all_samples()
    print("reading barcodes ...")
    bcs, sizes = {}, {}
    for _, r in reg.iterrows():
        key = f"{r['dataset']}:{r['sample']}"
        bcs[key] = read_barcodes(r.prefix, r.bc)
        sizes[key] = len(bcs[key])

    keys = list(bcs)
    rows = []
    for a, b in itertools.combinations(keys, 2):
        da, db = a.split(":")[0], b.split(":")[0]
        if da == db:
            continue                      # within-dataset pairs are not the question
        inter = len(bcs[a] & bcs[b])
        exp = sizes[a] * sizes[b] / WHITELIST
        rows.append({
            "a": a, "b": b, "n_a": sizes[a], "n_b": sizes[b],
            "observed_overlap": inter,
            "expected_by_chance": round(exp, 1),
            "ratio_obs_over_exp": round(inter / exp, 2) if exp > 0 else np.nan,
            "jaccard": round(inter / len(bcs[a] | bcs[b]), 4),
        })
    df = pd.DataFrame(rows).sort_values("ratio_obs_over_exp", ascending=False)
    df.to_csv(f"{OUT}/donor_overlap_barcodes.csv", index=False)

    print("\n=== top cross-dataset barcode overlaps ===")
    print(df.head(12).to_string(index=False))
    print("\n=== distribution of obs/expected ===")
    print(df.ratio_obs_over_exp.describe().round(3).to_string())

    suspicious = df[df.ratio_obs_over_exp > 3.0]
    print(f"\npairs with >3x chance overlap: {len(suspicious)}")

    # pseudobulk correlation among the Fudan-linked datasets
    fud = reg[reg.dataset.isin(["E-MTAB-12305", "GSE197461", "GSE208653"])]
    print(f"\ncomputing pseudobulk for {len(fud)} Fudan-linked samples ...")
    pbs = {}
    for _, r in fud.iterrows():
        key = f"{r['dataset']}:{r['sample']}"
        pbs[key] = pseudobulk(r.prefix, r.mtx, r.bc, r.ft)
        print(f"  {key}", flush=True)

    pk = list(pbs)
    cor = pd.DataFrame(index=pk, columns=pk, dtype=float)
    for a, b in itertools.combinations_with_replacement(pk, 2):
        sa, sb = pbs[a], pbs[b]
        shared = sa.index.intersection(sb.index)
        x = np.log1p(sa[shared] / max(sa.sum(), 1) * 1e6)
        y = np.log1p(sb[shared] / max(sb.sum(), 1) * 1e6)
        c = float(np.corrcoef(x, y)[0, 1])
        cor.loc[a, b] = cor.loc[b, a] = round(c, 4)
    cor.to_csv(f"{OUT}/donor_overlap_pseudobulk_corr.csv")

    print("\n=== highest cross-dataset pseudobulk correlations ===")
    pairs = []
    for a, b in itertools.combinations(pk, 2):
        if a.split(":")[0] != b.split(":")[0]:
            pairs.append({"a": a, "b": b, "pearson_log_cpm": cor.loc[a, b]})
    pdf = pd.DataFrame(pairs).sort_values("pearson_log_cpm", ascending=False)
    print(pdf.head(12).to_string(index=False))
    pdf.to_csv(f"{OUT}/donor_overlap_pairs.csv", index=False)

    print("\nVERDICT:")
    if len(suspicious) == 0 and pdf.pearson_log_cpm.max() < 0.98:
        print("  no evidence of shared libraries across datasets "
              f"(max obs/exp barcode ratio {df.ratio_obs_over_exp.max():.2f}, "
              f"max cross-dataset pseudobulk r {pdf.pearson_log_cpm.max():.3f})")
    else:
        print("  POSSIBLE SHARED LIBRARIES -- inspect the tables before "
              "treating these datasets as independent")


if __name__ == "__main__":
    main()
