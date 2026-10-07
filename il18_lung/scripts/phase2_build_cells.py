"""Phase 2 step 1: per-cell IL18 table for every entry in datasets.yaml (and nothing else).

For each entry: raw IL18 counts + per-cell library size, joined to the author's own cell-type
annotation. Output: data/cells/<entry id>.parquet with one row per cell.
"""
import gzip
import sys
from pathlib import Path

import cellxgene_census as cc
import h5py
import numpy as np
import pandas as pd
import scipy.io
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import CENSUS_VERSION, GENES
from h5ad_remote import obs_columns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "cells"
SPECIES_ORG = {"human": "homo_sapiens", "mouse": "mus_musculus"}
LOG = []


def census_cells(census, org, did):
    gid = GENES[org]
    ad = cc.get_anndata(census, organism=org, X_name="raw",
                        obs_value_filter=f"dataset_id == '{did}' and tissue_general == 'lung' and is_primary_data == True",
                        var_value_filter=f"feature_id == '{gid}'",
                        obs_column_names=["observation_joinid", "dataset_id", "donor_id", "cell_type", "disease",
                                          "assay", "suspension_type", "tissue", "raw_sum"])
    if ad.n_vars != 1 or ad.var["feature_id"].iloc[0] != gid:
        raise RuntimeError(f"{did}: Ensembl ID {gid} not matched")
    obs = ad.obs.reset_index(drop=True)
    obs["il18"] = np.asarray(ad.X.todense()).ravel()
    obs["lib"] = obs.pop("raw_sum").astype(float)
    LOG.append(f"{did}: Census var hit {gid} -> {ad.var['feature_name'].iloc[0]}; {len(obs)} primary lung cells")
    return obs


def author_obs(did, cols):
    cache = ROOT / "data" / f"srcobs_{did}.parquet"
    have = pd.read_parquet(cache) if cache.exists() else pd.DataFrame()
    need = [c for c in ["observation_joinid"] + cols if c not in have.columns]
    if need:
        add = obs_columns(did, CENSUS_VERSION, ["observation_joinid"] + [c for c in need if c != "observation_joinid"])
        have = add if have.empty else have.merge(add, on="observation_joinid", how="outer")
        have.to_parquet(cache)
    return have


def build_census(entry, census):
    org = SPECIES_ORG[entry["species"]]
    ann = entry["x_annotation"]
    extra = [c for c in [entry.get("x_split_by"), entry.get("x_condition_detail"), entry.get("x_condition_from_col")] if c]
    frames = []
    for did in entry["x_census_ids"]:
        cells = census_cells(census, org, did)
        cols = list(dict.fromkeys([ann] + extra + _query_cols(entry.get("x_subset", ""))))
        cols = [c for c in cols if c not in cells.columns or c == "cell_type"]
        src = author_obs(did, [c for c in cols if c != "cell_type"]) if [c for c in cols if c != "cell_type"] else None
        if src is not None:
            keep = ["observation_joinid"] + [c for c in cols if c in src.columns and c != "cell_type"]
            cells = cells.merge(src[keep], on="observation_joinid", how="left", validate="one_to_one")
        frames.append(cells)
    df = pd.concat(frames, ignore_index=True)
    if entry.get("x_subset"):
        df = df.query(entry["x_subset"]).copy()
    df["celltype_original"] = df[ann].astype(str)
    cmap = entry.get("x_condition_from")
    if cmap:
        key = df[entry.get("x_condition_from_col", "disease")].astype(str)
        df["condition"] = key.map(cmap)
        dropped = df["condition"].isna().sum()
        if dropped:
            LOG.append(f"{entry['id']}: {dropped} cells with unmapped condition key dropped: {sorted(key[df.condition.isna()].unique())}")
        df = df[df["condition"].notna()].copy()
    else:
        df["condition"] = entry["condition"]
    df["condition_detail"] = df[entry["x_condition_detail"]].astype(str) if entry.get("x_condition_detail") else ""
    df["sub_dataset"] = df[entry["x_split_by"]].astype(str) if entry.get("x_split_by") else ""
    df["gene_id_unavailable"] = False
    if entry["id"] == "H03_TabulaSapiens_lung":
        df = _apply_ts_decontx(df)
    return df


def _query_cols(q):
    import re
    q = re.sub(r"'[^']*'|\"[^\"]*\"", " ", q)  # drop quoted literals
    return [w for w in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", q)
            if w not in {"and", "or", "not", "in", "assay"} and not w[0].isupper()]


def _apply_ts_decontx(df):
    d = pd.read_parquet(ROOT / "data" / "ts_lung_decontx.parquet")
    n0 = len(df)
    df = df.merge(d, on="observation_joinid", how="inner", validate="one_to_one")
    LOG.append(f"H03_TabulaSapiens_lung: decontX layer joined for {len(df)}/{n0} cells; raw counts replaced by decontXcounts")
    df["il18"] = df.pop("il18_decontx")
    df["lib"] = df.pop("lib_decontx")
    return df


def _geo_features_hit(ids, names, gid, tag):
    hit = np.where(np.asarray(ids) == gid)[0]
    if len(hit) != 1:
        raise RuntimeError(f"{tag}: Ensembl ID {gid} not found in features")
    LOG.append(f"{tag}: features hit {gid} -> {names[hit[0]]}")
    return int(hit[0])


def build_geo_mtx(entry):
    gdir = ROOT / "data" / "geo" / entry["accession"]
    ann = pd.read_csv(gdir / f"{entry['accession']}_cell_annotation.csv.gz")
    gid = GENES["mus_musculus"]
    frames = []
    for tag, (cond, detail) in entry["x_samples"].items():
        mtx = next(gdir.glob(f"GSM*_{tag}_matrix.mtx.gz"))
        feat = next(gdir.glob(f"GSM*_{tag}_features.tsv.gz"))
        bc = next(gdir.glob(f"GSM*_{tag}_barcodes.tsv.gz"))
        f = pd.read_csv(feat, sep="\t", header=None)
        gi = _geo_features_hit(f[0], f[1], gid, f"{entry['id']}/{tag}")
        m = scipy.io.mmread(gzip.open(mtx)).tocsr()
        barcodes = pd.read_csv(bc, header=None)[0]
        cells = pd.DataFrame({"barcode": tag + "_" + barcodes,
                              "il18": m[gi].toarray().ravel(),
                              "lib": np.asarray(m.sum(axis=0)).ravel().astype(float)})
        cells = cells.merge(ann[["rowname", "celltype"]], left_on="barcode", right_on="rowname", how="inner")
        LOG.append(f"{entry['id']}/{tag}: {len(barcodes)} barcodes, {len(cells)} with author annotation")
        cells["condition"], cells["condition_detail"] = cond, detail
        cells["donor_id"] = tag
        frames.append(cells)
    df = pd.concat(frames, ignore_index=True)
    df["celltype_original"] = df["celltype"].astype(str)
    df["sub_dataset"] = ""
    df["gene_id_unavailable"] = False
    return df


def build_geo_h5(entry):
    gdir = ROOT / "data" / "geo" / entry["accession"]
    meta = pd.read_csv(gdir / f"{entry['accession']}_cell_metadata.csv.gz", index_col=0)
    gid = GENES["mus_musculus"]
    # The author metadata suffixes barcodes per sample: ES1 "_1_1", ES2 "_2_1", ES3 "_1", ES4 "_2"
    # (two-level Seurat merge). Mapped explicitly per GSM; match rates are logged as a check.
    suffix = {"GSM9366391": "_1_1", "GSM9366392": "_2_1", "GSM9366393": "_1", "GSM9366394": "_2"}
    files = sorted(gdir.glob("GSM*_filtered_feature_bc_matrix.h5"))
    frames = []
    for k, path in enumerate(files, start=1):
        with h5py.File(path, "r") as h:
            g = h["matrix"]
            ids = [x.decode() for x in g["features"]["id"][...]]
            names = [x.decode() for x in g["features"]["name"][...]]
            gi = _geo_features_hit(ids, names, gid, f"{entry['id']}/{path.name}")
            data, ind, ptr = g["data"][...], g["indices"][...], g["indptr"][...]
            shape = g["shape"][...]
            bcs = [x.decode() for x in g["barcodes"][...]]
        import scipy.sparse as sp
        m = sp.csc_matrix((data, ind, ptr), shape=shape)
        cells = pd.DataFrame({"barcode_raw": bcs, "il18": m[gi].toarray().ravel(),
                              "lib": np.asarray(m.sum(axis=0)).ravel().astype(float)})
        cells["key"] = cells.barcode_raw + suffix[path.name.split("_")[0]]
        frames.append(cells.assign(file=path.name))
    df = pd.concat(frames, ignore_index=True)
    df = df.merge(meta[["orig.ident", "celltype"]], left_on="key", right_index=True, how="inner")
    LOG.append(f"{entry['id']}: {len(df)} cells matched to author metadata "
               f"({df.groupby('orig.ident').size().to_dict()}); per-file origin {df.groupby('file')['orig.ident'].unique().to_dict()}")
    df = df[df["orig.ident"].isin(entry["x_samples"].keys())].copy()
    df["condition"] = df["orig.ident"].map(lambda s: entry["x_samples"][s][0])
    df["condition_detail"] = df["orig.ident"].map(lambda s: entry["x_samples"][s][1])
    df["donor_id"] = df["orig.ident"]
    df["celltype_original"] = df["celltype"].astype(str)
    df["sub_dataset"] = ""
    df["gene_id_unavailable"] = False
    return df


def main():
    spec = yaml.safe_load(open(ROOT / "datasets.yaml"))
    OUT.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        for e in spec["datasets"]:
            if only and e["id"] not in only:
                continue
            if not e.get("accession_verified"):
                LOG.append(f"{e['id']}: SKIPPED, accession_verified is false")
                continue
            if e["x_source"] == "census":
                df = build_census(e, census)
            elif e["x_source"] == "geo_10x_mtx":
                df = build_geo_mtx(e)
            else:
                df = build_geo_h5(e)
            df["donor_id"] = df["donor_id"].astype(str)
            keep = ["donor_id", "celltype_original", "condition", "condition_detail", "sub_dataset",
                    "il18", "lib", "gene_id_unavailable"]
            df[keep].to_parquet(OUT / f"{e['id']}.parquet")
            LOG.append(f"{e['id']}: wrote {len(df)} cells, {df.donor_id.nunique()} donors, "
                       f"conditions {df.condition.value_counts().to_dict()}")
            print(LOG[-1], flush=True)
    with open(ROOT / "results" / "phase2_build.log", "a") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
