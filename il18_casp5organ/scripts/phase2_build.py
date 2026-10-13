"""Phase 2 build (spec §6): per-cell IL18/CASP1/CASP3 table for every entry in datasets.yaml.

Genes taken by Ensembl ID only (§1.2). Normal cells only, filtered to the entry's suspension_type.
Author cell-type annotation kept verbatim (celltype_original). Output: data/cells/<id>.parquet.
Census source via get_anndata; the one GEO source (GSE103892) via downloaded matrix + annotation.
"""
import gzip
import sys
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import CENSUS_VERSION, GENES, GENE_ORDER

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "cells"
ORG = {"human": "homo_sapiens", "mouse": "mus_musculus"}
LOG = []


def build_census(census, e):
    org = ORG[e["species"]]
    gids = GENES[org]
    tis = e["tissue"].replace("_", " ")
    did = e["accession"].split()[-1]
    flt = (f"dataset_id == '{did}' and tissue_general == '{tis}' and is_primary_data == True "
           f"and disease == 'normal' and suspension_type == '{e['suspension_type']}'")
    ad = cc.get_anndata(census, organism=org, X_name="raw", obs_value_filter=flt,
                        var_value_filter="feature_id in ['%s']" % "','".join(gids.values()),
                        obs_column_names=["donor_id", "cell_type", "suspension_type", "assay", "raw_sum"])
    hit = {r.feature_id: r.feature_name for r in ad.var.itertuples()}
    present = [g for g, fid in gids.items() if fid in hit]
    X = pd.DataFrame(np.asarray(ad.X.todense()), columns=[hit[c] for c in ad.var.feature_id])
    obs = ad.obs.reset_index(drop=True)
    df = pd.DataFrame({"donor_id": obs.donor_id.astype(str).values,
                       "celltype_original": obs.cell_type.astype(str).values,
                       "lib": obs.raw_sum.astype(float).values})
    for g in GENE_ORDER:
        df[g] = X[g].values if g in present else np.nan
    df["gene_missing"] = len(present) < 3
    df["gene_id_unavailable"] = False
    LOG.append(f"{e['id']}: Census {did}/{tis}/{e['suspension_type']} genes={present} "
               f"{len(df)} normal cells, {df.donor_id.nunique()} donors")
    return df


def build_geo_gse103892(e):
    gdir = ROOT / "data" / "geo" / "GSE103892"
    gdir.mkdir(parents=True, exist_ok=True)
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103892/suppl"
    mtx = gdir / "GSE103892_Expression_Count_Matrix.txt.gz"
    ann = gdir / "GSE103892_Sample_Cell_Cluster_Information.txt.gz"
    import urllib.request
    for f, url in [(mtx, f"{base}/GSE103892_Expression_Count_Matrix.txt.gz"),
                   (ann, f"{base}/GSE103892_Sample_Cell_Cluster_Information.txt.gz")]:
        if not f.exists():
            urllib.request.urlretrieve(url, f)
    a = pd.read_csv(ann, sep="\t", compression="gzip")
    bc_col = a.columns[0]
    a["sample"] = a[bc_col].astype(str).str.split("_").str[0]
    keep = a[a["sample"].isin(e["x_baseline_samples"]) & a["cell.type"].notna()
            & ~a["cell.type"].astype(str).isin(["discarded", ""])].copy()
    # genes x cells matrix; first column = gene symbol (no Ensembl IDs on GEO -> gene_id_unavailable)
    m = pd.read_csv(mtx, sep="\t", compression="gzip", index_col=0)
    sym = {"IL18": "Il18", "CASP1": "Casp1", "CASP3": "Casp3"}
    present = [g for g in GENE_ORDER if sym[g] in m.index]
    cols = [c for c in keep[bc_col] if c in m.columns]
    sub = m[cols]
    lib = sub.sum(axis=0)
    df = pd.DataFrame({"cellbc": cols})
    df = df.merge(keep[[bc_col, "cell.type", "sample"]], left_on="cellbc", right_on=bc_col, how="left")
    df["donor_id"] = df["sample"]
    df["celltype_original"] = df["cell.type"].astype(str)
    df["lib"] = lib.reindex(cols).values.astype(float)
    for g in GENE_ORDER:
        df[g] = sub.loc[sym[g], cols].values.astype(float) if g in present else np.nan
    df["gene_missing"] = len(present) < 3
    df["gene_id_unavailable"] = True   # matrix carries symbols only (§1.2) -> reference only, not candidate
    LOG.append(f"{e['id']}: GEO GSE103892 symbols={present} (gene_id_unavailable=TRUE) "
               f"{len(df)} baseline cells, donors {sorted(df.donor_id.unique())}")
    return df[["donor_id", "celltype_original", "lib", *GENE_ORDER, "gene_missing", "gene_id_unavailable"]]


def main():
    spec = yaml.safe_load(open(ROOT / "datasets.yaml"))
    OUT.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        for e in spec["datasets"]:
            if only and e["id"] not in only:
                continue
            if not e.get("accession_verified"):
                LOG.append(f"{e['id']}: SKIP (accession_verified false)")
                continue
            df = build_census(census, e) if e.get("x_source") != "geo" else build_geo_gse103892(e)
            df.to_parquet(OUT / f"{e['id']}.parquet")
            print(LOG[-1], flush=True)
    with open(ROOT / "results" / "phase2_build.log", "a") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
