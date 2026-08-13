"""Manifest loading and gene-ID resolution.

The analysis code reaches data ONLY through a manifest. No dataset-specific
column name, path or label may appear anywhere in src/. Dataset peculiarities
live in the manifest or in a documented exclusion decision.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import yaml

REQUIRED_TOP = ["dataset_id", "local_path", "file_format", "species",
                "gene_id", "obs_columns", "celltype_map", "counts_layer"]
REQUIRED_OBS = ["donor", "celltype"]

VAR_INDEX_TYPES = {"symbol", "ensembl", "symbol_ensembl_concat"}


class ManifestError(Exception):
    """Raised when a manifest is malformed or contradicts the data."""


class PanelMatchError(Exception):
    """R11: panel match rate below threshold. Never downgraded to a warning."""


@dataclass
class Manifest:
    path: str
    raw: dict
    dataset_id: str = ""
    root: str = ""

    def __post_init__(self):
        self.dataset_id = self.raw["dataset_id"]
        self.root = os.path.dirname(os.path.abspath(self.path))

    # --- accessors -------------------------------------------------------
    def get(self, *keys, default=None):
        node = self.raw
        for k in keys:
            if not isinstance(node, dict) or k not in node:
                return default
            node = node[k]
        return node

    def obs_col(self, field_name):
        """Column name in .obs for a logical field, or None if absent."""
        return self.get("obs_columns", field_name)

    def fixed_value(self, field_name):
        return self.get("fixed_values", field_name)

    def resolve_path(self, p):
        if p is None:
            return None
        if os.path.isabs(p):
            return p
        # manifest paths are relative to the repository root (parent of manifests/)
        return os.path.normpath(os.path.join(self.root, "..", p))

    @property
    def local_path(self):
        return self.resolve_path(self.raw["local_path"])

    @property
    def celltype_map(self):
        return dict(self.raw["celltype_map"])


def load_manifest(path) -> Manifest:
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    missing = [k for k in REQUIRED_TOP if k not in raw]
    if missing:
        raise ManifestError(f"{path}: missing required keys {missing}")

    obs = raw.get("obs_columns") or {}
    fixed = raw.get("fixed_values") or {}
    for k in REQUIRED_OBS:
        if not obs.get(k):
            raise ManifestError(f"{path}: obs_columns.{k} is required and may not be null")
    # compartment may be null in obs only if supplied as a constant
    if not obs.get("compartment") and not fixed.get("compartment"):
        raise ManifestError(
            f"{path}: compartment must be given either as obs_columns.compartment "
            f"or as fixed_values.compartment")

    gid = raw.get("gene_id") or {}
    vit = gid.get("var_index_type")
    if vit not in VAR_INDEX_TYPES:
        raise ManifestError(
            f"{path}: gene_id.var_index_type must be one of {sorted(VAR_INDEX_TYPES)}, got {vit!r}")
    if vit == "ensembl" and not gid.get("symbol_column"):
        raise ManifestError(
            f"{path}: gene_id.symbol_column is required when var_index_type='ensembl'")
    if vit == "symbol_ensembl_concat":
        for k in ("id_separator", "symbol_field", "ensembl_field"):
            if gid.get(k) is None:
                raise ManifestError(
                    f"{path}: gene_id.{k} is required when var_index_type='symbol_ensembl_concat'")
    return Manifest(path=str(path), raw=raw)


# --- gene identity ------------------------------------------------------

def split_var_names(var_names, manifest: Manifest, symbol_column_values=None):
    """Return (symbols, ensembl_ids) aligned to var_names.

    Missing side is returned as a list of empty strings.
    """
    gid = manifest.raw["gene_id"]
    kind = gid["var_index_type"]
    n = len(var_names)
    if kind == "symbol":
        return list(var_names), [""] * n
    if kind == "ensembl":
        if symbol_column_values is None:
            raise ManifestError("ensembl var_index_type requires the symbol column values")
        return list(symbol_column_values), list(var_names)
    # symbol_ensembl_concat, e.g. "CCL5_ENSG00000271503"
    sep = gid["id_separator"]
    si, ei = int(gid["symbol_field"]), int(gid["ensembl_field"])
    syms, ens = [], []
    for v in var_names:
        parts = str(v).split(sep)
        # a symbol may itself contain the separator (HLA-G with '-', TRBC2 etc.).
        # Anchor on the ensembl field position counted from the correct end.
        if ei == -1 or ei == len(parts) - 1:
            e = parts[-1] if len(parts) > 1 else ""
            s = sep.join(parts[:-1]) if len(parts) > 1 else str(v)
        else:
            e = parts[ei] if ei < len(parts) else ""
            s = parts[si] if si < len(parts) else str(v)
        syms.append(s)
        ens.append(e)
    return syms, ens


def match_panel(panel_symbols, panel_ensembl, symbols, ensembl_ids,
                min_rate=0.90, label=""):
    """Map panel genes to matrix column indices.

    R11: raises PanelMatchError below `min_rate`. Silent continuation is barred.
    Symbol match is tried first, then Ensembl (version suffixes stripped).
    """
    sym_index = {}
    for i, s in enumerate(symbols):
        if s and s not in sym_index:
            sym_index[s] = i
    ens_index = {}
    for i, e in enumerate(ensembl_ids):
        if not e:
            continue
        base = e.split(".")[0]
        if base not in ens_index:
            ens_index[base] = i

    hits, misses, how = {}, [], {}
    for sym, ens in zip(panel_symbols, panel_ensembl):
        idx = sym_index.get(sym)
        route = "symbol"
        if idx is None and ens:
            idx = ens_index.get(str(ens).split(".")[0])
            route = "ensembl"
        if idx is None:
            misses.append(sym)
        else:
            hits[sym] = idx
            how[sym] = route
    rate = len(hits) / max(1, len(panel_symbols))
    if rate < min_rate:
        raise PanelMatchError(
            f"{label}: panel match rate {rate:.1%} < {min_rate:.0%} "
            f"({len(hits)}/{len(panel_symbols)} matched). Missing: {misses[:30]}. "
            f"Refusing to continue (R11).")
    return hits, misses, rate, how
