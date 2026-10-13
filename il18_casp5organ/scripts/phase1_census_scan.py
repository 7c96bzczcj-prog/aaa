"""Phase 1 (spec §5): Census scan of IL18/CASP1/CASP3 across 5 organs, human and mouse separately,
stratified by suspension_type. Per-dataset processing to stay within memory. Normal cells only
(disease=='normal'); excluded-count recorded. Output: results/census_scan.tsv (locate only, not for
candidate judgement).

Scope (user-approved): lung/liver/pancreas scan ALL normal primary datasets; brain/spinal cord scan
the approved reference atlases only (Siletti human, Yao+cortex-snRNA mouse) — Census brain has 180+
datasets, out of scope to scan exhaustively.
"""
import sys
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import CENSUS_VERSION, GENES, GENE_ORDER, SPECIES_LABEL, TISSUES, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OBS_COLS = ["dataset_id", "donor_id", "cell_type", "tissue", "tissue_general", "assay",
            "suspension_type", "disease", "development_stage", "sex"]

# brain/spinal: only these reference datasets are scanned (approved scope)
BRAIN_SPINAL_SCOPE = {
    "homo_sapiens": {
        "brain": ["e1f595f6-ba2c-495e-9bee-7056f116b1e4"],         # Siletti MTG (nucleus)
        "spinal cord": ["090da8ea-46e8-40df-bffc-1f78e1538d27"],   # Siletti spinal cord (nucleus)
    },
    "mus_musculus": {
        "brain": ["d7291f04-fbbb-4d65-990a-f01fa44e915b",          # Yao 2023 cortex+hippocampus (cell)
                  "1229ecc2-b067-4664-91da-0251aec31574"],         # cortex WT snRNA (nucleus)
        "spinal cord": [],                                          # none in Census -> GEO/gap
    },
}
LOG = []


def datasets_in(census, org, tissue):
    """Normal primary dataset_ids for a tissue; for brain/spinal restrict to the approved scope."""
    obs = census["census_data"][org].obs
    o = obs.read(value_filter=f"tissue_general == '{tissue}' and is_primary_data == True and disease == 'normal'",
                 column_names=["dataset_id"]).concat().to_pandas()
    ids = sorted(o.dataset_id.unique())
    if tissue in ("brain", "spinal cord"):
        scope = set(BRAIN_SPINAL_SCOPE[org][tissue])
        ids = [d for d in ids if d in scope]
    return ids


def scan_dataset(census, org, tissue, did):
    gids = GENES[org]
    ad = cc.get_anndata(census, organism=org, X_name="raw",
                        obs_value_filter=f"dataset_id == '{did}' and tissue_general == '{tissue}' "
                                         f"and is_primary_data == True",
                        var_value_filter="feature_id in ['%s']" % "','".join(gids.values()),
                        obs_column_names=OBS_COLS + ["raw_sum"])
    fid_to_key = {fid: g for g, fid in gids.items()}
    present = [g for g, fid in gids.items() if fid in set(ad.var.feature_id)]
    if not present:
        LOG.append(f"{org}/{tissue}/{did}: no target gene in var — skipped")
        return None
    X = pd.DataFrame(np.asarray(ad.X.todense()), columns=[fid_to_key[c] for c in ad.var.feature_id])
    obs = ad.obs.reset_index(drop=True)
    for g in GENE_ORDER:               # attach genes while obs and X are aligned (before any filter)
        obs[g] = X[g].values if g in X.columns else np.nan
    n_total = len(obs)
    obs = obs[obs.disease == "normal"].copy()
    excl = n_total - len(obs)
    obs["lib"] = obs.pop("raw_sum").astype(float)
    keys = ["dataset_id", "tissue", "suspension_type", "donor_id", "cell_type", "assay"]
    rows = []
    for key, grp in obs.groupby(keys, observed=True):
        row = dict(zip(keys, key))
        row["species"] = SPECIES_LABEL[org]
        row["cell_type_original"] = row.pop("cell_type")
        row["n_cells"] = len(grp)
        lib = grp.lib.sum()
        for g in GENE_ORDER:
            if g in present:
                v = grp[g].values
                row[f"n_cells_{g}_pos"] = int((v >= 1).sum())
                row[f"pct_{g}"] = float((v >= 1).mean())
                row[f"logcpm_{g}"] = float(np.log2(v.sum() / lib * 1e6 + 1)) if lib > 0 else np.nan
            else:
                row[f"n_cells_{g}_pos"] = np.nan
                row[f"pct_{g}"] = np.nan
                row[f"logcpm_{g}"] = np.nan
        rows.append(row)
    LOG.append(f"{org}/{tissue}/{did}: genes {present}; {n_total} cells ({excl} non-normal excluded); "
               f"{obs.suspension_type.nunique()} susp; {len(rows)} groups")
    return pd.DataFrame(rows)


def main():
    frames = []
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        for org in GENES:
            for tissue in TISSUES:
                ids = datasets_in(census, org, tissue)
                LOG.append(f"== {org}/{tissue}: {len(ids)} datasets to scan")
                print(LOG[-1], flush=True)
                for did in ids:
                    df = scan_dataset(census, org, tissue, did)
                    if df is not None:
                        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    cols = ["species", "dataset_id", "tissue", "suspension_type", "donor_id", "cell_type_original", "assay",
            "n_cells"] + [f"n_cells_{g}_pos" for g in GENE_ORDER] + [f"pct_{g}" for g in GENE_ORDER] + \
           [f"logcpm_{g}" for g in GENE_ORDER]
    write_tsv(out[cols], ROOT / "results" / "census_scan.tsv")
    (ROOT / "results" / "phase1_scan.log").write_text("\n".join(LOG) + "\n")
    print(f"\nwrote census_scan.tsv: {len(out)} rows, {out.dataset_id.nunique()} datasets")


if __name__ == "__main__":
    main()
