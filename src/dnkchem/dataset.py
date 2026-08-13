"""Manifest -> in-memory dataset.

Two entry points, deliberately separate:

  load_source_dataset()  reads the ORIGINAL published files and applies the
                         manifest's column and label mapping. Used only by
                         src/ingest.py.
  load_dataset()         reads the canonical h5ad written by ingest. Used by
                         every analysis script.

Every dataset passes through ingest, so analysis code never re-parses a
source format and never sees a dataset-specific column name. The canonical
obs schema is: donor, compartment, celltype, celltype_raw, subset, library
(+ sex, gestational_week, cycle_phase where published).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

from . import loaders
from .manifest import Manifest, ManifestError, match_panel, split_var_names

CANON_FIELDS = ["donor", "compartment", "celltype", "library",
                "sex", "gestational_week", "cycle_phase"]
OPTIONAL_FIELDS = ["sex", "gestational_week", "cycle_phase"]


class Dataset:
    def __init__(self, X, obs, symbols, ensembl, manifest: Manifest):
        self.X = X
        self.obs = obs
        self.symbols = np.asarray(symbols, dtype=object)
        self.ensembl = np.asarray(ensembl, dtype=object)
        self.manifest = manifest

    @property
    def dataset_id(self):
        return self.manifest.dataset_id

    def panel_columns(self, panel, min_rate=0.90):
        """R11 gate. Raises PanelMatchError below the threshold."""
        return match_panel(panel["gene_symbol"].tolist(),
                           panel["ensembl_id"].fillna("").tolist(),
                           self.symbols, self.ensembl,
                           min_rate=min_rate, label=self.dataset_id)


def _resolve_obs_field(obs_raw, manifest: Manifest, field_name, n_cells):
    """obs column if named, else the fixed constant, else None."""
    col = manifest.obs_col(field_name)
    if col:
        if col not in obs_raw.columns:
            raise ManifestError(
                f"{manifest.dataset_id}: obs_columns.{field_name} = {col!r} "
                f"not found in the obs table (available: {list(obs_raw.columns)})")
        return obs_raw[col].astype(str).to_numpy()
    fixed = manifest.fixed_value(field_name)
    if fixed is not None:
        return np.array([str(fixed)] * n_cells, dtype=object)
    return None


def load_source_dataset(manifest: Manifest, progress=None) -> Dataset:
    """Read the original published files. Only src/ingest.py should call this."""
    src = manifest.get("source") or {}
    fmt = src.get("format", manifest.raw["file_format"])
    path = manifest.resolve_path(src.get("path", manifest.raw["local_path"]))
    if not os.path.exists(path):
        raise FileNotFoundError(f"{manifest.dataset_id}: source not found at {path}")

    X, cell_ids, gene_ids, adata = loaders.read_source(
        fmt, path, counts_layer=manifest.raw.get("counts_layer", "X"), progress=progress)

    obs_path = src.get("obs_table")
    if obs_path:
        obs_raw = loaders.read_obs_table(manifest.resolve_path(obs_path))
        if obs_raw.index.duplicated().any():
            n = int(obs_raw.index.duplicated().sum())
            raise ManifestError(
                f"{manifest.dataset_id}: obs table {obs_path} has {n} duplicate cell "
                f"ids; the matrix-to-metadata join would be ambiguous.")
        obs_raw = obs_raw.reindex(pd.Index(cell_ids))
        n_missing = int(obs_raw.isna().all(axis=1).sum())
        if n_missing:
            raise ManifestError(
                f"{manifest.dataset_id}: {n_missing} cells in the matrix have no row "
                f"in the obs table {obs_path}. Refusing to guess.")
    elif adata is not None:
        obs_raw = adata.obs.copy()
    else:
        raise ManifestError(f"{manifest.dataset_id}: no obs table available")

    n = len(cell_ids)
    obs = pd.DataFrame(index=pd.Index([str(c) for c in cell_ids], name="cell_id"))
    for f in CANON_FIELDS:
        v = _resolve_obs_field(obs_raw, manifest, f, n)
        if v is not None:
            obs[f] = v
    obs["celltype_raw"] = obs["celltype"]
    obs["subset"] = obs["celltype_raw"].map(manifest.celltype_map)
    if "library" not in obs.columns:
        obs["library"] = "unspecified"

    gid = manifest.raw["gene_id"]
    sym_col_vals = None
    if gid["var_index_type"] == "ensembl":
        col = gid["symbol_column"]
        if adata is None or col not in adata.var.columns:
            raise ManifestError(
                f"{manifest.dataset_id}: gene_id.symbol_column={col!r} not in var")
        sym_col_vals = adata.var[col].astype(str).to_numpy()
    symbols, ensembl = split_var_names(gene_ids, manifest, sym_col_vals)

    return Dataset(X, obs, symbols, ensembl, manifest)


def load_dataset(manifest: Manifest) -> Dataset:
    """Read the canonical h5ad produced by src/ingest.py."""
    import anndata as ad
    import scipy.sparse as sp

    path = manifest.local_path
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{manifest.dataset_id}: canonical h5ad missing at {path}. "
            f"Run:  python src/ingest.py {manifest.path}")
    a = ad.read_h5ad(path)
    obs = a.obs.copy()
    for c in obs.columns:
        if str(obs[c].dtype) == "category":
            obs[c] = obs[c].astype(object)
    # ingest stores unmapped cells as an empty subset; restore the NaN so that
    # `.notna()` continues to mean "has a mapped subset" everywhere downstream.
    if "subset" in obs.columns:
        obs["subset"] = obs["subset"].replace("", np.nan)
    for f in OPTIONAL_FIELDS:
        if f not in obs.columns:
            obs[f] = np.nan
    return Dataset(sp.csr_matrix(a.X), obs,
                   a.var["gene_symbol"].to_numpy(),
                   a.var["ensembl_id"].to_numpy(), manifest)


def unit_indices(obs, mask, by):
    """Yield (key_tuple, positional row indices) for each group under `mask`.

    Positional throughout: label-based lookup silently breaks when cell ids
    are not unique, and a published obs table is not guaranteed to have unique
    ids. Groups are yielded in sorted key order for reproducible output.
    """
    cols = [obs[c].astype(object).where(obs[c].notna(), "").astype(str).to_numpy()
            for c in by]
    keys = list(zip(*cols)) if cols else []
    groups = {}
    for i in np.flatnonzero(np.asarray(mask, dtype=bool)):
        groups.setdefault(keys[i], []).append(i)
    for k in sorted(groups):
        yield k, np.asarray(groups[k], dtype=int)


def load_panel(path):
    p = pd.read_csv(path, sep="\t")
    need = {"gene_symbol", "ensembl_id", "category", "role", "notes"}
    missing = need - set(p.columns)
    if missing:
        raise ValueError(f"{path}: panel missing columns {sorted(missing)}")
    p["ensembl_id"] = p["ensembl_id"].fillna("")
    p["notes"] = p["notes"].fillna("")
    return p


def eprint(*a):
    print(*a, file=sys.stderr, flush=True)
