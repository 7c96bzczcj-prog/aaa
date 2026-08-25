"""
Per-dataset NK recount for cervical tissue.

Run:  python build_dataset.py <DATASET>
where DATASET in {E-MTAB-12305, GSE208653, GSE197461, GSE173231}

Stages, in the order the task requires:
  1. load cellranger-filtered counts, per sample
  2. QC
  3. doublet detection + removal (scrublet; scDblFinder/DoubletFinder equivalent)
  4. coarse clustering -> decontX ambient correction, per sample
  5. pooled-within-dataset clustering (harmony over samples) -> major lineages
  6. T/NK compartment re-clustered from scratch, ignoring any author labels
  7. NK gating on positive AND negative markers, at cluster level
  8. per-donor counts under three denominators + CD56bright/dim split

Nothing is pooled across donors for the reported numbers; pooling is only ever
used to define a shared clustering within one dataset.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cerv_common as cc

sc.settings.verbosity = 1
sc.settings.n_jobs = 4
OUT = "/home/user/aaa/results/cervical"
CACHE = "/home/user/cervical_work/cache"
# bump whenever QC / doublet / decontX behaviour changes, to invalidate the cache
PIPELINE_VERSION = 2
os.makedirs(OUT, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

# ---- gating thresholds (primary set; scanned for sensitivity in 02_table) ----
T_NEG_MAX_DET = 0.25     # a T marker counts as "negative" below this detection rate
MIN_T_NEG = 3            # >=3 of the 4 T markers must be negative  (task step 2)
NKG7_MIN = 0.60
GNLY_MIN = 0.40
NK_RECEPTOR_MIN = {"KLRD1": 0.40, "KLRF1": 0.20, "NCAM1": 0.10}   # any one suffices


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# --------------------------------------------------------------------------
def load_sample(row):
    X, bc, genes = cc.read_10x(row.prefix, row.mtx, row.bc, row.ft)
    X, genes = cc.collapse_duplicate_genes(X, genes)
    ad = sc.AnnData(X)
    ad.var_names = genes
    ad.obs_names = [f"{row['sample']}|{b}" for b in bc]
    ad.obs["sample"] = row["sample"]
    ad.obs["dataset"] = row["dataset"]
    ad.obs["tissue"] = row["tissue"]
    ad.obs["donor"] = row["donor"]
    ad.obs["histology"] = row["histology"]
    return ad


def qc_and_decontaminate(ad, sample):
    n0 = ad.n_obs
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    keep = ((ad.obs.n_genes_by_counts >= 200)
            & (ad.obs.total_counts >= 500)
            & (ad.obs.pct_counts_mt < 20))
    ad = ad[keep].copy()
    n_qc = ad.n_obs
    if n_qc < 100:
        log(f"  {sample}: only {n_qc} cells survive QC -> dropped")
        return None, dict(n_droplets=n0, n_qc=n_qc, n_doublet=0, n_final=0, theta_median=np.nan)

    # Deliberately NOT filtering genes here. Per-sample gene filtering followed by
    # an inner join across samples can silently delete a marker from the whole
    # dataset: NCAM1 is low and patchy, so a sample with few NK cells would drop
    # it, and the intersection would then remove CD56 from every sample. The
    # gating depends on those exact genes, so the full gene set is carried through.
    ad.layers["counts"] = ad.X.copy()

    # ---- doublets (raw counts) -------------------------------------------
    # Doublet removal is mandatory here: T-NK doublets are the single largest
    # source of false NK calls. A silent fallback to "no doublets" would quietly
    # inflate every number downstream, so a failure is fatal, not tolerated.
    exp_rate = min(0.08, 0.008 * n_qc / 1000)
    sc.pp.scrublet(ad, expected_doublet_rate=exp_rate, verbose=False)
    dbl = ad.obs["predicted_doublet"].values.astype(bool)
    if dbl.sum() == 0:
        # Scrublet picks its cutoff from the bimodality of the simulated-doublet
        # score distribution. When that distribution is unimodal it silently
        # returns zero doublets, which is not the same as there being none.
        # Fall back to calling the top expected-rate fraction by score, so the
        # mandatory doublet step still removes something rather than no-op'ing.
        k = int(round(exp_rate * n_qc))
        score = ad.obs["doublet_score"].values
        thr = np.partition(score, -k)[-k] if 0 < k < len(score) else np.inf
        dbl = score >= thr
        log(f"  {sample}: scrublet found no threshold; removing top "
            f"{int(dbl.sum())} cells by doublet score ({exp_rate:.1%} expected)")
    n_dbl = int(dbl.sum())
    ad = ad[~dbl].copy()

    # ---- coarse clusters to drive decontX --------------------------------
    tmp = ad.copy()
    tmp.X = tmp.layers["counts"].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=2000)
    tmp = tmp[:, tmp.var.highly_variable].copy()
    sc.pp.scale(tmp, max_value=10)
    sc.tl.pca(tmp, n_comps=min(30, tmp.n_obs - 1, tmp.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(tmp, n_neighbors=15)
    sc.tl.leiden(tmp, resolution=1.0, key_added="lc",
                 flavor="igraph", n_iterations=2, directed=False)
    clusters = tmp.obs["lc"].values.astype(str)
    del tmp

    # ---- decontX ----------------------------------------------------------
    theta, Xd = cc.decontx(ad.layers["counts"], clusters)
    # the pre-correction matrix is deliberately not retained: nothing downstream
    # reads it, and keeping a third copy of a 36k-gene matrix was enough to push
    # the largest dataset into the OOM killer
    ad.layers["counts"] = Xd
    ad.X = Xd
    ad.obs["decontx_theta"] = theta

    # drop cells left with too little signal after decontamination
    tot = np.asarray(Xd.sum(axis=1)).ravel()
    ad = ad[tot >= 300].copy()

    stats = dict(n_droplets=n0, n_qc=n_qc, n_doublet=n_dbl, n_final=ad.n_obs,
                 theta_median=float(np.median(theta)))
    log(f"  {sample}: {n0} droplets -> {n_qc} QC -> -{n_dbl} doublets -> "
        f"{ad.n_obs} cells (theta={stats['theta_median']:.3f})")
    return ad, stats


# --------------------------------------------------------------------------
def annotate_lineages(ad):
    """Cluster the whole dataset and label each cluster with a major lineage."""
    ad.X = ad.layers["counts"].copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers["lognorm"] = ad.X.copy()
    sc.pp.highly_variable_genes(ad, n_top_genes=2500, batch_key="sample")
    hv = ad.var.highly_variable.values
    sub = ad[:, hv].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=50, svd_solver="arpack")
    rep = cc.harmony_embed(sub, "sample") if ad.obs["sample"].nunique() > 1 else "X_pca"
    sc.pp.neighbors(sub, n_neighbors=15, use_rep=rep)
    sc.tl.leiden(sub, resolution=1.0, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    ad.obs["leiden"] = sub.obs["leiden"].values
    ad.obsm[rep] = sub.obsm[rep]
    del sub

    # score every lineage signature, then take the per-cluster argmax
    for name, gl in cc.LINEAGE.items():
        present = [g for g in gl if g in ad.var_names]
        if present:
            sc.tl.score_genes(ad, present, score_name=f"sig_{name}", use_raw=False)
        else:
            ad.obs[f"sig_{name}"] = 0.0
    sig_cols = [f"sig_{n}" for n in cc.LINEAGE if n != "immune"]
    per_cluster = ad.obs.groupby("leiden", observed=True)[sig_cols + ["sig_immune"]].mean()
    lineage = per_cluster[sig_cols].idxmax(axis=1).str.replace("sig_", "", regex=False)
    ad.obs["lineage"] = ad.obs["leiden"].map(lineage).astype(str)
    ad.uns["lineage_scores"] = per_cluster
    return ad


def gate_nk(tnk):
    """Re-cluster the T/NK compartment and call NK clusters on +/- markers."""
    tnk.X = tnk.layers["counts"].copy()
    sc.pp.normalize_total(tnk, target_sum=1e4)
    sc.pp.log1p(tnk)
    tnk.layers["lognorm"] = tnk.X.copy()
    sc.pp.highly_variable_genes(tnk, n_top_genes=2000, batch_key="sample")
    sub = tnk[:, tnk.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=min(40, sub.n_obs - 1), svd_solver="arpack")
    rep = cc.harmony_embed(sub, "sample") if tnk.obs["sample"].nunique() > 1 else "X_pca"
    sc.pp.neighbors(sub, n_neighbors=15, use_rep=rep)
    sc.tl.leiden(sub, resolution=1.5, key_added="tnk_cluster",
                 flavor="igraph", n_iterations=2, directed=False)
    tnk.obs["tnk_cluster"] = sub.obs["tnk_cluster"].values
    del sub

    # every gating gene must actually exist, or the rule silently weakens
    missing = [g for g in cc.NK_POS + cc.T_NEG if g not in tnk.var_names]
    if missing:
        raise RuntimeError(f"gating markers absent from the matrix: {missing}")

    markers = sorted(set(cc.NK_POS + cc.T_NEG + ["TRBC1", "TRAC", "CD2", "IL7R",
                                                 "FCGR3A", "GZMB", "KLRC1",
                                                 "XCL1", "SELL", "MKI67", "PTPRC"]))
    tab = cc.cluster_marker_table(tnk, "tnk_cluster", markers, layer="counts")

    n_t_neg = sum((tab[g] < T_NEG_MAX_DET).astype(int) for g in cc.T_NEG if g in tab)
    receptor_ok = np.zeros(len(tab), dtype=bool)
    for g, thr in NK_RECEPTOR_MIN.items():
        if g in tab:
            receptor_ok |= (tab[g] >= thr).values
    is_nk = ((n_t_neg >= MIN_T_NEG).values
             & (tab.get("NKG7", 0) >= NKG7_MIN).values
             & (tab.get("GNLY", 0) >= GNLY_MIN).values
             & receptor_ok)
    tab["n_T_markers_negative"] = n_t_neg
    tab["is_NK"] = is_nk
    nk_clusters = set(tab.index[is_nk])
    tnk.obs["is_NK"] = tnk.obs["tnk_cluster"].isin(nk_clusters).values
    return tnk, tab


def cd56_split(tnk):
    """CD56bright (NCAM1-high / FCGR3A-low) vs CD56dim, among gated NK cells."""
    nk = tnk[tnk.obs.is_NK].copy()
    if nk.n_obs == 0:
        return pd.Series(dtype=float)
    ln = nk.layers["lognorm"]
    def vec(g):
        if g not in nk.var_names:
            return np.zeros(nk.n_obs)
        c = ln[:, nk.var_names.get_loc(g)]
        return c.toarray().ravel() if sp.issparse(c) else np.asarray(c).ravel()
    ncam1, fcgr3a = vec("NCAM1"), vec("FCGR3A")
    bright = (ncam1 > 0) & (fcgr3a <= 0)
    out = pd.Series(bright, index=nk.obs_names)
    return out


# --------------------------------------------------------------------------
def main(dataset):
    reg = cc.registry()
    reg = reg[reg.dataset == dataset].reset_index(drop=True)
    if reg.empty:
        raise SystemExit(f"unknown dataset {dataset}")
    log(f"=== {dataset}: {len(reg)} samples ===")

    ads, qc_rows = [], []
    for _, row in reg.iterrows():
        # QC + doublets + decontX are deterministic and by far the slowest part,
        # so each processed library is cached and reused across reruns. The
        # version tag invalidates the cache whenever those steps change.
        cpath = f"{CACHE}/s{PIPELINE_VERSION}_{dataset}_{row['sample']}.h5ad"
        spath = f"{CACHE}/s{PIPELINE_VERSION}_{dataset}_{row['sample']}.json"
        if os.path.exists(cpath) and os.path.exists(spath):
            ad = sc.read_h5ad(cpath)
            st = json.load(open(spath))
            log(f"  {row['sample']}: reused cache ({ad.n_obs} cells)")
        else:
            ad = load_sample(row)
            ad, st = qc_and_decontaminate(ad, row["sample"])
            if ad is not None:
                ad.write(cpath)
            json.dump(st, open(spath, "w"))
        st.update(dataset=dataset, sample=row["sample"], tissue=row["tissue"],
                  donor=row["donor"], histology=row["histology"])
        qc_rows.append(st)
        if ad is not None:
            ads.append(ad)
    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(f"{OUT}/qc_{dataset}.csv", index=False)

    log("concatenating")
    ad = sc.concat(ads, join="outer", fill_value=0, label=None)
    del ads
    ad.obs_names_make_unique()
    log(f"pooled (for clustering only): {ad.shape}")

    log("major-lineage clustering")
    ad = annotate_lineages(ad)
    log("lineage composition:\n" + str(ad.obs.lineage.value_counts()))

    tnk_mask = ad.obs.lineage.isin(["T_NK"]).values
    log(f"T/NK compartment: {tnk_mask.sum()} cells")
    tnk = ad[tnk_mask].copy()
    tnk, gate_tab = gate_nk(tnk)
    gate_tab.to_csv(f"{OUT}/nkgate_{dataset}.csv")
    log("NK clusters: " + str(sorted(gate_tab.index[gate_tab.is_NK].tolist())))
    log(f"NK cells gated: {int(tnk.obs.is_NK.sum())} / {tnk.n_obs} T-NK")

    bright = cd56_split(tnk)

    # ---- per-donor tallies ------------------------------------------------
    ad.obs["is_NK"] = False
    ad.obs.loc[tnk.obs_names[tnk.obs.is_NK.values], "is_NK"] = True
    ad.obs["is_bright"] = False
    if len(bright):
        ad.obs.loc[bright.index[bright.values], "is_bright"] = True

    ad.obs["is_immune"] = ad.obs.lineage.isin(cc.IMMUNE_LINEAGES).values
    ad.obs["is_lymph"] = ad.obs.lineage.isin(cc.LYMPHOID_LINEAGES).values

    rows = []
    for smp, g in ad.obs.groupby("sample", observed=True):
        nk = int(g.is_NK.sum())
        rows.append({
            "dataset": dataset, "sample": smp,
            "donor": g.donor.iloc[0], "tissue": g.tissue.iloc[0],
            "histology": g.histology.iloc[0],
            "n_cells_total": len(g),
            "n_immune": int(g.is_immune.sum()),
            "n_lymphocyte": int(g.is_lymph.sum()),
            "n_TNK": int((g.lineage == "T_NK").sum()),
            "n_NK": nk,
            "pct_of_all": 100 * nk / len(g) if len(g) else np.nan,
            "pct_of_immune": 100 * nk / g.is_immune.sum() if g.is_immune.sum() else np.nan,
            "pct_of_lymphocyte": 100 * nk / g.is_lymph.sum() if g.is_lymph.sum() else np.nan,
            "n_NK_bright": int(g.is_bright.sum()),
            "pct_bright_of_NK": 100 * g.is_bright.sum() / nk if nk else np.nan,
            "decontx_theta_median": float(g.decontx_theta.median()),
            "median_UMI": float(g.total_counts.median()),
        })
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/nk_per_sample_{dataset}.csv", index=False)
    log("\n" + res[["sample", "tissue", "n_cells_total", "n_NK", "pct_of_all",
                    "pct_of_immune", "pct_of_lymphocyte", "pct_bright_of_NK"]]
        .to_string(index=False))

    # keep the T/NK object: sensitivity scans re-gate from it without recomputing
    tnk.write(f"{CACHE}/tnk_{dataset}.h5ad")
    ad.obs.to_csv(f"{OUT}/obs_{dataset}.csv.gz", compression="gzip")
    log(f"=== {dataset} done ===")


if __name__ == "__main__":
    main(sys.argv[1])
