#!/usr/bin/env python
"""Convert a manifest's source into the canonical h5ad at `local_path`.

Canonical form: cells x genes, sparse CSR, raw integer counts in X, unified
obs columns, var carrying gene_symbol and ensembl_id. Run once per dataset;
every analysis script reads only this file.

    python src/ingest.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import verify_integer_counts  # noqa: E402
from dnkchem.dataset import OPTIONAL_FIELDS, load_source_dataset  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402


def clean_obs(obs, dataset_id):
    """Make obs h5ad-writable without losing meaning.

    An all-null column cannot be written as a vlen string, and a partly-null
    object column writes inconsistently. Columns published as empty are
    dropped (load_dataset re-creates them as NaN); the rest become
    categoricals, where the empty string carries "unmapped".
    """
    obs = obs.copy()
    dropped = []
    for c in list(obs.columns):
        if obs[c].isna().all():
            if c in OPTIONAL_FIELDS:
                obs = obs.drop(columns=[c])
                dropped.append(c)
            else:
                raise SystemExit(
                    f"[{dataset_id}] required obs column {c!r} is entirely empty")
    if dropped:
        print(f"[{dataset_id}] obs columns not published by this dataset, dropped: {dropped}")
    for c in obs.columns:
        obs[c] = obs[c].astype(object).where(obs[c].notna(), "").astype(str).astype("category")
    return obs


def build_var(symbols, ensembl):
    """var index must be unique; disambiguate repeats with the Ensembl id."""
    sym = pd.Series([str(s) for s in symbols])
    ens = pd.Series([str(e) for e in ensembl])
    dup = sym.duplicated(keep=False)
    idx = sym.where(~dup, sym + "__" + ens)
    if idx.duplicated().any():
        n = int(idx.duplicated().sum())
        idx = idx.where(~idx.duplicated(), idx + "__" + pd.Series(np.arange(len(idx)), dtype=str))
        print(f"  note: {n} var names still collided after adding the Ensembl id; "
              f"suffixed by position")
    return pd.DataFrame({"gene_symbol": sym.to_numpy(), "ensembl_id": ens.to_numpy()},
                        index=pd.Index(idx.to_numpy(), name="gene"))


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
    ds = load_source_dataset(mf, progress=progress)
    print(f"[{mf.dataset_id}] matrix {ds.X.shape[0]} cells x {ds.X.shape[1]} genes "
          f"({ds.X.nnz} nonzero, {time.time() - t0:.0f}s)", flush=True)

    is_int, mx, checked = verify_integer_counts(ds.X)
    if not is_int:
        raise SystemExit(
            f"[{mf.dataset_id}] counts_layer is not integer counts (max {mx}). "
            f"The manifest must point at RAW counts.")
    print(f"[{mf.dataset_id}] integer counts verified on {checked} cells (max {mx:.0f})")

    obs = ds.obs.copy()
    # library surrogate: the barcode prefix identifies the sequencing run when
    # the source publishes no per-cell library column (R6 needs a grouping).
    if "library" not in obs.columns or (obs["library"] == "unspecified").all():
        pref = [str(c).split("_")[0] for c in obs.index]
        if len(set(pref)) > 1:
            obs["library"] = pref
            print(f"[{mf.dataset_id}] library taken from barcode prefix: "
                  f"{len(set(pref))} distinct runs")

    n_unmapped = int(obs["subset"].isna().sum())
    print(f"[{mf.dataset_id}] cells with a mapped subset: "
          f"{len(obs) - n_unmapped}/{len(obs)} ({n_unmapped} unmapped)")

    import anndata as ad
    a = ad.AnnData(X=sp.csr_matrix(ds.X.astype(np.int32)),
                   obs=clean_obs(obs, mf.dataset_id),
                   var=build_var(ds.symbols, ds.ensembl))
    a.uns["dataset_id"] = mf.dataset_id
    a.uns["accession"] = str(mf.raw.get("accession", ""))
    a.uns["panel_agnostic_ingest"] = "1"

    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".tmp"
    a.write_h5ad(tmp, compression="gzip")
    os.replace(tmp, out)
    print(f"[{mf.dataset_id}] wrote {out} ({os.path.getsize(out)/1e6:.0f} MB, "
          f"{time.time() - t0:.0f}s total)")

    # read back: the canonical file is what every later script sees
    from dnkchem.dataset import load_dataset
    chk = load_dataset(mf)
    assert chk.X.shape == ds.X.shape, (chk.X.shape, ds.X.shape)
    assert int(chk.X.sum()) == int(ds.X.sum()), "count total changed on round-trip"
    assert int(chk.obs["subset"].notna().sum()) == len(obs) - n_unmapped
    print(f"[{mf.dataset_id}] round-trip verified: {chk.X.shape[0]} cells, "
          f"{int(chk.X.sum())} counts, {int(chk.obs['subset'].notna().sum())} mapped cells")

    print(f"\nSet raw_counts_verified: true in {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
