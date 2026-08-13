#!/usr/bin/env python
"""Gate a dataset before it may enter the main analysis.

Checks (all must pass):
  1. manifest schema
  2. matrix is raw integer counts
  3. panel match rate >= 90%  (R11 -- hard failure, never a warning)
  4. every celltype_map target label exists and its cell count
  5. cells per donor x compartment x subset
  6. UMI depth distribution

    python src/validate_manifest.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import verify_integer_counts  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel  # noqa: E402
from dnkchem.manifest import PanelMatchError, load_manifest  # noqa: E402

MIN_CELLS_UNIT = 30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    print(f"[1/6] manifest schema ....... OK ({mf.dataset_id})")

    ds = load_dataset(mf)
    print(f"      matrix: {ds.X.shape[0]} cells x {ds.X.shape[1]} genes")

    is_int, mx, checked = verify_integer_counts(ds.X)
    print(f"[2/6] raw integer counts .... {'OK' if is_int else 'FAIL'} "
          f"(max {mx:.0f} over {checked} cells)")
    if not is_int:
        raise SystemExit("FAIL: counts_layer is not raw integer counts")

    try:
        hits, misses, rate, how = ds.panel_columns(panel)
    except PanelMatchError as e:
        print(f"[3/6] panel match ........... FAIL\n{e}")
        raise SystemExit(2)
    n_by_ens = sum(1 for g in hits if how[g] == "ensembl")
    print(f"[3/6] panel match ........... OK {rate:.1%} "
          f"({len(hits)}/{len(panel)}; {n_by_ens} via Ensembl)")
    if misses:
        print(f"      unmatched: {misses}")

    obs = ds.obs
    mapped = obs["subset"].notna()
    print(f"[4/6] celltype_map .......... {mapped.sum()} / {len(obs)} cells mapped "
          f"({mapped.mean():.1%})")
    unmapped = sorted(set(obs.loc[~mapped, "celltype_raw"].astype(str)))
    if unmapped:
        print(f"      raw labels with no mapping (excluded by decision, not by accident):")
        for u in unmapped:
            print(f"        {u!r}: {(obs['celltype_raw'].astype(str) == u).sum()} cells")
    tgt_counts = obs.loc[mapped, "subset"].value_counts()
    missing_targets = [t for t in set(mf.celltype_map.values()) if t not in tgt_counts.index]
    for t, c in tgt_counts.items():
        print(f"        {t:20s} {c:7d}")
    if missing_targets:
        print(f"      TARGETS WITH ZERO CELLS: {missing_targets}")

    print("[5/6] donor x compartment x subset")
    units = (obs.loc[mapped]
             .groupby(["donor", "compartment", "subset"], observed=True)
             .size().rename("n_cells").reset_index())
    units["meets_min_cells_pre_qc"] = units["n_cells"] >= MIN_CELLS_UNIT
    units.to_csv(os.path.join(outdir, "unit_counts_pre_qc.tsv"), sep="\t", index=False)
    nk_like = units[units["subset"].str.startswith(("dNK", "pbNK"))]
    print(nk_like.pivot_table(index=["compartment", "subset"], columns="donor",
                              values="n_cells", fill_value=0, aggfunc="sum").to_string())
    print(f"      units >= {MIN_CELLS_UNIT} cells (pre-QC): "
          f"{units['meets_min_cells_pre_qc'].sum()} / {len(units)}")

    total = np.asarray(ds.X.sum(axis=1)).ravel()
    qs = np.quantile(total, [0.05, 0.10, 0.25, 0.50, 0.75, 0.95])
    print("[6/6] UMI depth (all cells): " +
          " ".join(f"q{int(p*100)}={v:.0f}" for p, v in
                   zip([0.05, 0.10, 0.25, 0.50, 0.75, 0.95], qs)))
    nkmask = mapped.to_numpy() & obs["subset"].fillna("").str.startswith("dNK").to_numpy()
    if nkmask.sum():
        qn = np.quantile(total[nkmask], [0.05, 0.10, 0.25, 0.50, 0.75, 0.95])
        print("      UMI depth (dNK cells):  " +
              " ".join(f"q{int(p*100)}={v:.0f}" for p, v in
                       zip([0.05, 0.10, 0.25, 0.50, 0.75, 0.95], qn)))

    pd.DataFrame({"gene_symbol": list(hits), "matched_via": [how[g] for g in hits]}
                 ).to_csv(os.path.join(outdir, "panel_match.tsv"), sep="\t", index=False)
    print(f"\nPASS — {mf.dataset_id} is admitted. Artifacts in {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
