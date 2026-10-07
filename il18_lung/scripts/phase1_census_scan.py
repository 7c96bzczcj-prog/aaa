"""Phase 1 (spec §3): scan CELLxGENE Census lung for IL18 raw counts, human and mouse separately.

Pulls one gene column plus obs metadata; per-cell IL18 counts are cached to data/ for Phase 2.
"""
import sys
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import CENSUS_VERSION, GENES, OBS_COLS, SPECIES_LABEL, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OBS_FILTER = "tissue_general == 'lung' and is_primary_data == True"


def pull(census, organism):
    gid = GENES[organism]
    ad = cc.get_anndata(census, organism=organism, X_name="raw",
                        obs_value_filter=OBS_FILTER,
                        var_value_filter=f"feature_id == '{gid}'",
                        obs_column_names=OBS_COLS + ["soma_joinid"])
    if ad.n_vars != 1 or ad.var["feature_id"].iloc[0] != gid:
        raise RuntimeError(f"{organism}: Ensembl ID {gid} not matched in Census var")
    obs = ad.obs.copy()
    obs["il18"] = np.asarray(ad.X.todense()).ravel()
    obs["feature_name"] = ad.var["feature_name"].iloc[0]
    return obs


def scan(obs, species):
    keys = ["dataset_id", "donor_id", "cell_type", "disease", "assay"]
    g = obs.groupby(keys, observed=True)
    out = g.agg(n_cells=("il18", "size"),
                n_cells_IL18_pos=("il18", lambda x: int((x >= 1).sum())),
                il18_sum=("il18", "sum"), lib_sum=("raw_sum", "sum")).reset_index()
    out["pct_expressing"] = out.n_cells_IL18_pos / out.n_cells
    out["mean_logcpm"] = np.log2(out.il18_sum / out.lib_sum * 1e6 + 1)
    out.insert(0, "species", species)
    out = out.rename(columns={"cell_type": "cell_type_original"})
    return out[["species", "dataset_id", "donor_id", "cell_type_original", "disease", "assay",
                "n_cells", "n_cells_IL18_pos", "pct_expressing", "mean_logcpm"]]


def main():
    (ROOT / "data").mkdir(exist_ok=True)
    scans, log = [], []
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
        datasets.to_csv(ROOT / "data" / "census_datasets.tsv", sep="\t", index=False)
        for org in GENES:
            obs = pull(census, org)
            obs.to_parquet(ROOT / "data" / f"census_lung_{org}.parquet")
            scans.append(scan(obs, SPECIES_LABEL[org]))
            log.append(f"{org}: {len(obs)} cells, {obs.dataset_id.nunique()} datasets, "
                       f"gene {GENES[org]} -> symbol {obs.feature_name.iloc[0]}")
    write_tsv(pd.concat(scans), ROOT / "results" / "census_scan.tsv")
    print("\n".join(log))


if __name__ == "__main__":
    main()
