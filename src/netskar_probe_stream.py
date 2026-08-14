"""Stream the NETSKAR2024 NK atlas counts layer and keep only what is needed.

Question (check 1): are PTN and OGN detected in NK cells across tissues, and
does that detection track ambient stromal pickup rather than the tissue?

The atlas is NK-only, so ruler B (NK / source-lineage CPM) cannot be built —
there is no stromal compartment in the file.  Ruler A survives in its
prior-selected form: COL1A1/COL1A2/COL3A1/DCN/LUM cannot be transcribed by an
NK cell, so their entire signal inside an NK barcode is pickup.  That is the
ambient index used here, and it is admissible for exactly the reason set out
in docs/TRANSFERABLE.md section 5 — the controls are chosen from prior
exclusivity, not from a statistic, and they cross a boundary the target cannot
cross.

VIM is deliberately *not* an ambient control: leukocytes transcribe vimentin.
It is carried only so that the contrast with the collagens is on the record.

The matrix is CSR over 89,216 cells x 11,866 genes with 70.5 M non-zeros, so
two gene columns cannot be sliced out; the row index must be walked.  This
walks it once (~564 MB of range requests, no disk) and writes a per-cell table
of: total counts over retained genes, plus counts for each probe gene.

No scanpy (R10).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from remote_h5 import open_remote, read_str_array  # noqa: E402

URL = "https://zenodo.org/api/records/8434224/files/all_nk_cells.h5ad/content"
OUT = Path("out/NETSKAR2024")

# --- probe panel, fixed here before any expression value is read ------------
TARGETS = ["PTN", "OGN", "SPP1"]          # SPP1 = OPN, the third Fu 2017 factor
AMBIENT_STROMAL = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]
AMBIENT_EPITHELIAL = ["EPCAM", "KRT19"]
NOT_A_CONTROL = ["VIM"]                    # leukocytes transcribe it; see above
NK_IDENTITY = ["NKG7", "KLRD1", "GNLY", "PRF1", "KLRF1"]

OBS_KEEP = ["sample", "dataset", "source", "tumor_type", "subset",
            "low_res_subset", "cell_type", "total_counts"]

ROW_BLOCK = 4000


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    f, fobj = open_remote(URL, block=8 << 20, max_cache=16)

    genes = read_str_array(f["var"]["_index"])
    gene_pos = {g: i for i, g in enumerate(genes)}

    panel = TARGETS + AMBIENT_STROMAL + AMBIENT_EPITHELIAL + NOT_A_CONTROL + NK_IDENTITY
    missing = [g for g in panel if g not in gene_pos]
    present = [g for g in panel if g in gene_pos]
    rate = len(present) / len(panel)
    print(f"panel {len(present)}/{len(panel)} matched ({rate:.1%}); missing: {missing}")
    # R11: an unexplained match failure is an exception, not a silent continue.
    if rate < 0.90:
        raise SystemExit(f"panel match rate {rate:.1%} < 90% (R11)")

    cols = np.array([gene_pos[g] for g in present], dtype=np.int64)
    col_rank = {int(c): k for k, c in enumerate(cols)}

    grp = f["layers"]["counts"]
    indptr = grp["indptr"][:].astype(np.int64)
    d_idx, d_dat = grp["indices"], grp["data"]
    n_cells = indptr.shape[0] - 1

    total = np.zeros(n_cells, dtype=np.int64)
    n_genes_det = np.zeros(n_cells, dtype=np.int32)
    probe = np.zeros((n_cells, len(present)), dtype=np.int32)

    for start in range(0, n_cells, ROW_BLOCK):
        stop = min(start + ROW_BLOCK, n_cells)
        lo, hi = int(indptr[start]), int(indptr[stop])
        idx = d_idx[lo:hi]
        dat = d_dat[lo:hi]
        if not np.allclose(dat, np.round(dat)):
            raise SystemExit("counts layer is not integral — not raw counts")
        dat = np.round(dat).astype(np.int64)

        rows = np.repeat(np.arange(start, stop),
                         np.diff(indptr[start:stop + 1]).astype(np.int64))
        np.add.at(total, rows, dat)
        np.add.at(n_genes_det, rows, 1)

        hit = np.isin(idx, cols)
        if hit.any():
            hr = rows[hit]
            hc = np.array([col_rank[int(c)] for c in idx[hit]], dtype=np.int64)
            np.add.at(probe, (hr, hc), dat[hit])

        pct = 100.0 * stop / n_cells
        print(f"  {stop:6d}/{n_cells} cells  {pct:5.1f}%  "
              f"fetched {fobj.bytes_fetched/1e6:7.1f} MB", flush=True)

    obs = {}
    for c in OBS_KEEP:
        if c not in f["obs"]:
            print(f"  obs column absent, skipped: {c}")
            continue
        node = f["obs"][c]
        if hasattr(node, "keys") and "categories" in node:
            obs[c] = read_str_array(node)
        else:
            obs[c] = node[:]

    df = pd.DataFrame(obs)
    df["retained_total"] = total
    df["retained_n_genes"] = n_genes_det
    for k, g in enumerate(present):
        df[f"n_{g}"] = probe[:, k]

    path = OUT / "cell_probe_counts.csv.gz"
    df.to_csv(path, index=False)
    print(f"\nwrote {path}  rows={len(df)}  "
          f"total fetched {fobj.bytes_fetched/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
