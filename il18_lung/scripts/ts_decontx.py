"""Extract IL18 from the Tabula Sapiens author decontX layer (ambient-corrected) for lung cells, via HTTP range reads."""
import sys
from pathlib import Path

import fsspec
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import CENSUS_VERSION, GENES
from h5ad_remote import URL, _col

ROOT = Path(__file__).resolve().parents[1]
DID = "53d208b0-2cfd-4366-9866-c3c6114081bc"


def main():
    f = fsspec.open(URL.format(ver=CENSUS_VERSION, did=DID), "rb", block_size=32 * 2**20).open()
    with h5py.File(f, "r") as h:
        var_ids = np.array([x.decode() for x in h["var"]["_index"][...]])
        hit = np.where(var_ids == GENES["homo_sapiens"])[0]
        if len(hit) != 1:
            raise RuntimeError("Tabula Sapiens var: Ensembl ID not matched")
        gi = int(hit[0])
        tissue = np.asarray(_col(h["obs"]["tissue_in_publication"]))
        jid = h["obs"]["observation_joinid"]
        jid = np.asarray(_col(jid))
        rows = np.where(tissue == "Lung")[0]
        L = h["layers"]["decontXcounts"]
        indptr = L["indptr"][...]
        vals = np.zeros(len(rows), dtype=np.float32)
        # read contiguous runs of lung rows
        breaks = np.where(np.diff(rows) != 1)[0] + 1
        for run in np.split(np.arange(len(rows)), breaks):
            r0, r1 = rows[run[0]], rows[run[-1]] + 1
            s, e = indptr[r0], indptr[r1]
            ind = L["indices"][s:e]
            dat = L["data"][s:e]
            rp = indptr[r0:r1 + 1] - s
            for k, ridx in enumerate(run):
                seg = slice(rp[k], rp[k + 1])
                m = ind[seg] == gi
                if m.any():
                    vals[ridx] = dat[seg][m][0]
            # decontX library size per cell is needed for CPM
        libs = np.zeros(len(rows), dtype=np.float64)
        for run in np.split(np.arange(len(rows)), breaks):
            r0, r1 = rows[run[0]], rows[run[-1]] + 1
            s, e = indptr[r0], indptr[r1]
            dat = L["data"][s:e]
            rp = indptr[r0:r1 + 1] - s
            libs[run] = np.add.reduceat(dat, rp[:-1]) if len(dat) else 0
    out = pd.DataFrame({"observation_joinid": jid[rows], "il18_decontx": vals, "lib_decontx": libs})
    out.to_parquet(ROOT / "data" / "ts_lung_decontx.parquet")
    print(len(out), "lung cells; runs", len(breaks) + 1, "; IL18>0:", int((vals > 0).sum()))


if __name__ == "__main__":
    main()
