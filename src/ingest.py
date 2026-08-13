#!/usr/bin/env python
"""Convert a manifest's source into the canonical h5ad at `local_path`.

Canonical form: cells x genes, sparse CSR, raw integer counts in X, unified
obs columns. Run once per dataset; every later script reads only the h5ad.

    python src/ingest.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import verify_integer_counts  # noqa: E402
from dnkchem.dataset import load_dataset  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    out = mf.local_path
    if os.path.exists(out) and not args.force:
        print(f"{out} exists; --force to rebuild")
        return 0

    t0 = time.time()

    def progress(n):
        print(f"  ... {n} genes read ({time.time() - t0:.0f}s)", flush=True)

    print(f"[{mf.dataset_id}] reading source", flush=True)
    ds = load_dataset(mf, progress=progress)
    print(f"[{mf.dataset_id}] matrix {ds.X.shape[0]} cells x {ds.X.shape[1]} genes "
          f"({ds.X.nnz} nonzero, {time.time() - t0:.0f}s)", flush=True)

    is_int, mx, checked = verify_integer_counts(ds.X)
    if not is_int:
        raise SystemExit(
            f"[{mf.dataset_id}] counts_layer is not integer counts (max {mx}). "
            f"The manifest must point at RAW counts.")
    print(f"[{mf.dataset_id}] integer counts verified on {checked} cells (max {mx:.0f})")

    import anndata as ad
    import pandas as pd

    var = pd.DataFrame({"gene_symbol": ds.symbols, "ensembl_id": ds.ensembl},
                       index=pd.Index([str(s) for s in ds.symbols], name="gene_symbol_index"))
    var.index = pd.Index(
        pd.Series([str(s) for s in ds.symbols]).where(
            ~pd.Series([str(s) for s in ds.symbols]).duplicated(),
            pd.Series([f"{s}__{e}" for s, e in zip(ds.symbols, ds.ensembl)])).to_numpy(),
        name="gene")

    obs = ds.obs.copy()
    # library surrogate: the barcode prefix identifies the sequencing run when
    # the source publishes no per-cell library column (R6 needs some grouping).
    if (obs["library"] == "unspecified").all():
        pref = [str(c).split("_")[0] for c in obs.index]
        if len(set(pref)) > 1:
            obs["library"] = pref
            print(f"[{mf.dataset_id}] library taken from barcode prefix: "
                  f"{len(set(pref))} distinct runs")

    a = ad.AnnData(X=sp.csr_matrix(ds.X.astype(np.int32)), obs=obs, var=var)
    a.uns["dataset_id"] = mf.dataset_id
    a.uns["accession"] = mf.raw.get("accession", "")
    a.uns["ingested_from"] = str(mf.get("source", "path"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    a.write_h5ad(out, compression="gzip")
    print(f"[{mf.dataset_id}] wrote {out} ({os.path.getsize(out)/1e6:.0f} MB, "
          f"{time.time() - t0:.0f}s total)")

    print(f"\nSet raw_counts_verified: true in {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
