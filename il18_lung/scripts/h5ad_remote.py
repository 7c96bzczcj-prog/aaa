"""Read selected obs columns from a CELLxGENE source h5ad over HTTPS range requests (no full download)."""
import fsspec
import h5py
import pandas as pd

URL = "https://cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com/cell-census/{ver}/h5ads/{did}.h5ad"


def _col(g):
    if isinstance(g, h5py.Group):  # categorical
        cats = g["categories"][...]
        cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
        return pd.Categorical.from_codes(g["codes"][...], categories=cats)
    v = g[...]
    return [x.decode() if isinstance(x, bytes) else x for x in v] if v.dtype.kind in "OS" else v


def obs_columns(did, ver, cols=None):
    f = fsspec.open(URL.format(ver=ver, did=did), "rb", block_size=8 * 2**20).open()
    with h5py.File(f, "r") as h:
        obs = h["obs"]
        if cols is None:
            return sorted(obs.keys())
        return pd.DataFrame({c: _col(obs[c]) for c in cols if c in obs})
