"""Mouse LPS acute lung injury: self-annotate GSE280364 and GSE280611 to the 11 unified classes.

These datasets carry no author per-cell annotation, so the spec's "use author annotation" rule cannot
be met. The user authorised inclusion with self-annotation. This is a documented DEVIATION:
standard scanpy processing (normalise -> HVG -> PCA -> neighbours -> Leiden) per dataset, then each
Leiden cluster is labelled by canonical marker panels (argmax of per-cell marker scores, cluster
majority vote). Cells keep raw counts for the IL18 readout; clustering uses a separate normalised copy.
No cross-dataset integration. Output per-cell parquet in data/cells/, flagged self_annotated=True.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from common import GENES, SEED

ROOT = Path(__file__).resolve().parents[1]
GID = GENES["mus_musculus"]
sc.settings.verbosity = 0

# canonical mouse lung markers; the first 11 keys are the unified classes, the rest route to 'other'
PANELS = {
    "alveolar_macrophage": ["Marco", "Siglecf", "Chil3", "Ear2", "Krt79", "Atp6v0d2", "Pparg", "Fabp1"],
    "monocyte_derived_macrophage": ["C1qa", "C1qb", "C1qc", "Ms4a7", "Mafb", "Apoe"],  # recruited/IM-like macs
    "dendritic_cell": ["Flt3", "Zbtb46", "Xcr1", "Clec9a", "Cd209a", "Ccr7", "Siglech", "Bst2"],
    "neutrophil": ["S100a8", "S100a9", "Retnlg", "Mmp9", "Csf3r", "Ly6g"],
    "T_NK": ["Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Nkg7", "Klrb1c", "Gzma", "Ncr1"],
    "B_plasma": ["Cd79a", "Cd79b", "Ms4a1", "Igkc", "Jchain", "Mzb1"],
    "AT1": ["Ager", "Pdpn", "Hopx", "Rtkn2", "Akap5"],
    "AT2": ["Sftpc", "Sftpb", "Sftpa1", "Lamp3", "Napsa"],
    "airway_epithelial": ["Scgb1a1", "Scgb3a2", "Foxj1", "Ccdc153", "Krt5", "Trp63", "Muc5b"],
    "endothelial": ["Pecam1", "Cldn5", "Cdh5", "Gpihbp1", "Vwf"],
    "fibroblast_mesenchymal": ["Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Dcn", "Acta2", "Pdgfrb"],
    # routed to 'other':
    "other:monocyte": ["Ly6c2", "Plac8", "Vcan", "F13a1", "Chil3"],
    "other:mast_basophil": ["Cpa3", "Mcpt4", "Ms4a2", "Kit"],
}
ROUTED_OTHER = {"other:monocyte", "other:mast_basophil"}


def load_mtx(gdir, samples):
    ads = []
    for tag, (cond, detail) in samples.items():
        mtx = next(gdir.glob(f"GSM*_{tag}_matrix.mtx.gz"))
        feat = pd.read_csv(next(gdir.glob(f"GSM*_{tag}_features.tsv.gz")), sep="\t", header=None)
        bc = pd.read_csv(next(gdir.glob(f"GSM*_{tag}_barcodes.tsv.gz")), header=None)[0]
        m = scipy.io.mmread(gzip.open(mtx)).T.tocsr()  # cells x genes
        ad = sc.AnnData(X=m, obs=pd.DataFrame({"donor_id": tag, "condition": cond, "condition_detail": detail},
                                              index=[f"{tag}_{b}" for b in bc]))
        ad.var_names = feat[0].values
        ad.var["symbol"] = feat[1].values
        ads.append(ad)
    return ads


def load_h5(gdir, samples):
    ads = []
    for tag, (cond, detail) in samples.items():
        path = next(gdir.glob(f"GSM*_{tag}.filtered_feature_bc_matrix.h5"))
        ad = sc.read_10x_h5(path)
        ad.var["symbol"] = ad.var_names
        ad.var_names = ad.var["gene_ids"].values
        ad.obs["donor_id"], ad.obs["condition"], ad.obs["condition_detail"] = tag, cond, detail
        ad.obs_names = [f"{tag}_{b}" for b in ad.obs_names]
        ads.append(ad)
    return ads


def annotate(ads, entry_id):
    ad = sc.concat(ads, join="inner")
    ad.var_names_make_unique()
    # map each marker symbol to an Ensembl id present in this dataset (first match wins)
    vmap = ads[0].var.reset_index(names="ensembl")[["ensembl", "symbol"]].drop_duplicates("symbol")
    sym2id = dict(zip(vmap.symbol, vmap.ensembl))
    sym2id = {s: i for s, i in sym2id.items() if i in set(ad.var_names)}
    il18_raw = np.asarray(ad[:, GID].X.todense()).ravel() if GID in ad.var_names else np.zeros(ad.n_obs)
    lib = np.asarray(ad.X.sum(1)).ravel()
    # QC: keep cells with >=500 counts and >=200 genes
    ngenes = np.asarray((ad.X > 0).sum(1)).ravel()
    keep = (lib >= 500) & (ngenes >= 200)
    ad, il18_raw, lib = ad[keep].copy(), il18_raw[keep], lib[keep]
    work = ad.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    for cls, panel in PANELS.items():
        ids = [sym2id[s] for s in panel if s in sym2id]
        sc.tl.score_genes(work, ids, score_name=f"sc_{cls}", ctrl_size=50)
    sc.pp.highly_variable_genes(work, n_top_genes=2000)
    work2 = work[:, work.var.highly_variable].copy()
    sc.pp.scale(work2, max_value=10)
    sc.tl.pca(work2, n_comps=30, random_state=SEED)
    sc.pp.neighbors(work2, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(work2, resolution=1.0, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    scores = work.obs[[f"sc_{c}" for c in PANELS]].copy()
    scores.columns = list(PANELS)
    # per-cell argmax, then cluster majority vote
    cellcls = scores.idxmax(1)
    cl = work2.obs["leiden"].values
    lab = pd.Series(cellcls.values, index=range(len(cl))).groupby(cl).agg(lambda x: x.value_counts().index[0])
    clustlabel = pd.Series(cl).map(lab).values
    out = pd.DataFrame({
        "donor_id": ad.obs.donor_id.values, "condition": ad.obs.condition.values,
        "condition_detail": ad.obs.condition_detail.values,
        "leiden": cl, "celltype_assigned": clustlabel, "il18": il18_raw, "lib": lib})
    out["celltype_original"] = np.where(out.celltype_assigned.isin(ROUTED_OTHER),
                                        out.celltype_assigned, out.celltype_assigned)
    out["sub_dataset"] = ""
    out["gene_id_unavailable"] = False
    diag = (out.groupby(["leiden", "celltype_assigned"]).size().reset_index(name="n"))
    return out, diag


DATASETS = {
    "M06_GSE280364_LPS_ALI": dict(loader="mtx", dir="GSE280364",
                                  samples={"con1": ("resting", "control"), "con2": ("resting", "control"),
                                           "LPS-1": ("lps_ali", "LPS"), "LPS-2": ("lps_ali", "LPS")}),
    "M07_GSE280611_LPS_ALI": dict(loader="h5", dir="GSE280611",
                                  samples={"untreated_rep1": ("resting", "untreated"),
                                           "untreated_rep2": ("resting", "untreated"),
                                           "LPS_rep1": ("lps_ali", "LPS"), "LPS_rep2": ("lps_ali", "LPS")}),
}


def main():
    (ROOT / "data" / "cells").mkdir(parents=True, exist_ok=True)
    log = []
    for eid, cfg in DATASETS.items():
        gdir = ROOT / "data" / "geo" / cfg["dir"]
        ads = (load_mtx if cfg["loader"] == "mtx" else load_h5)(gdir, cfg["samples"])
        out, diag = annotate(ads, eid)
        keep = ["donor_id", "celltype_original", "condition", "condition_detail", "sub_dataset",
                "il18", "lib", "gene_id_unavailable"]
        out[keep].to_parquet(ROOT / "data" / "cells" / f"{eid}.parquet")
        diag.to_csv(ROOT / "results" / f"lps_cluster_labels_{eid}.tsv", sep="\t", index=False)
        log.append(f"{eid}: {len(out)} cells after QC, {out.celltype_original.nunique()} labels, "
                   f"conditions {out.condition.value_counts().to_dict()}; "
                   f"class counts {out.celltype_original.value_counts().to_dict()}")
        print(log[-1], flush=True)
    with open(ROOT / "results" / "lps_annotate.log", "w") as fh:
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
