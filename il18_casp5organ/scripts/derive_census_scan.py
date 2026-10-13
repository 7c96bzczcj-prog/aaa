"""Produce census_scan.tsv (spec §5 deliverable) from the per-cell data of the curated candidate strata.

Phase-1 scope note (recorded in run_log): the broad five-organ landscape was mapped via Census coverage
queries; the scan table itself is built over the confirmed candidate strata (user-approved scoping),
per dataset × tissue × suspension_type × donor × author cell type.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import GENE_ORDER, write_tsv

ROOT = Path(__file__).resolve().parents[1]


def main():
    spec = {e["id"]: e for e in yaml.safe_load(open(ROOT / "datasets.yaml"))["datasets"]}
    rows = []
    for eid, e in spec.items():
        p = ROOT / "data" / "cells" / f"{eid}.parquet"
        if not p.exists():
            continue
        c = pd.read_parquet(p)
        for (donor, ct), g in c.groupby(["donor_id", "celltype_original"], observed=True):
            lib = g.lib.sum()
            r = {"species": e["species"], "dataset_id": e["accession"].split()[-1], "tissue": e["tissue"],
                 "suspension_type": e["suspension_type"], "donor_id": donor, "cell_type_original": ct,
                 "assay": e["platform"], "n_cells": len(g)}
            for gene in GENE_ORDER:
                v = g[gene].values
                r[f"n_cells_{gene}_pos"] = int((v >= 1).sum()) if not np.isnan(v).all() else np.nan
                r[f"pct_{gene}"] = float((v >= 1).mean()) if not np.isnan(v).all() else np.nan
                r[f"logcpm_{gene}"] = float(np.log2(np.nansum(v) / lib * 1e6 + 1)) if lib > 0 and not np.isnan(v).all() else np.nan
            rows.append(r)
    out = pd.DataFrame(rows)
    cols = (["species", "dataset_id", "tissue", "suspension_type", "donor_id", "cell_type_original",
             "assay", "n_cells"] + [f"n_cells_{g}_pos" for g in GENE_ORDER]
            + [f"pct_{g}" for g in GENE_ORDER] + [f"logcpm_{g}" for g in GENE_ORDER])
    write_tsv(out[cols], ROOT / "results" / "census_scan.tsv")
    print(f"census_scan.tsv: {len(out)} rows, {out.dataset_id.nunique()} datasets, "
          f"{out.groupby('species').tissue.nunique().to_dict()}")


if __name__ == "__main__":
    main()
