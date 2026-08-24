"""
Validate the decontX implementation on one sample before trusting it anywhere.

Acceptance test: ambient correction must strip lineage-foreign transcripts from
cells that cannot be making them (epithelial keratins inside immune cells) while
leaving the lineage's own defining transcripts essentially intact. A correction
that removes both, or neither, is not usable.
"""
import sys
import time

import numpy as np
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, "/home/user/aaa/scripts/cervical")
import cerv_common as cc

sc.settings.verbosity = 1

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "H1"
row = cc.emtab_samples().set_index("sample").loc[SAMPLE]

t0 = time.time()
X, bc, genes = cc.read_10x(row.prefix, row.mtx, row.bc, row.ft)
X, genes = cc.collapse_duplicate_genes(X, genes)
print(f"loaded {SAMPLE}: {X.shape[0]} droplets x {X.shape[1]} genes  ({time.time()-t0:.0f}s)")

ad = sc.AnnData(X)
ad.var_names = genes
ad.obs_names = bc
ad.var_names_make_unique()

# --- QC -------------------------------------------------------------------
ad.var["mt"] = ad.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
print("median genes/cell:", np.median(ad.obs.n_genes_by_counts),
      "median UMI:", np.median(ad.obs.total_counts),
      "median pct_mt:", round(float(np.median(ad.obs.pct_counts_mt)), 2))

keep = (ad.obs.n_genes_by_counts >= 200) & (ad.obs.total_counts >= 500) & (ad.obs.pct_counts_mt < 20)
ad = ad[keep].copy()
print("after QC:", ad.shape)

# --- coarse clusters to drive decontX --------------------------------------
ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000)
adh = ad[:, ad.var.highly_variable].copy()
sc.pp.scale(adh, max_value=10)
sc.tl.pca(adh, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(adh, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adh, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
ad.obs["leiden"] = adh.obs["leiden"].values
print("clusters:", ad.obs.leiden.nunique(), "sizes:", ad.obs.leiden.value_counts().head(8).to_dict())

# --- decontX ---------------------------------------------------------------
t0 = time.time()
theta, Xd = cc.decontx(ad.layers["counts"], ad.obs["leiden"].values.astype(str), verbose=True)
print(f"decontX done in {time.time()-t0:.0f}s; theta median={np.median(theta):.3f} "
      f"IQR=({np.percentile(theta,25):.3f},{np.percentile(theta,75):.3f})")

raw = sp.csr_matrix(ad.layers["counts"])
print(f"total counts removed: {1 - Xd.sum()/raw.sum():.1%}")

# --- acceptance test -------------------------------------------------------
gi = {g: ad.var_names.get_loc(g) for g in
      ["PTPRC", "CD3D", "CD3E", "KRT5", "KRT13", "KRT14", "EPCAM", "NKG7", "GNLY", "LYZ"]
      if g in ad.var_names}


def frac_pos(M, col, rows):
    v = M[rows][:, col].toarray().ravel()
    return (v > 0).mean()


# immune cells = clusters where PTPRC detection is high and keratin is not the identity
det_ptprc = np.zeros(ad.obs.leiden.nunique())
labels = np.array(sorted(ad.obs.leiden.unique(), key=int))
imm_clusters, epi_clusters = [], []
for i, k in enumerate(labels):
    rows = np.where(ad.obs.leiden.values == k)[0]
    p = frac_pos(raw, gi["PTPRC"], rows)
    e = max(frac_pos(raw, gi[g], rows) for g in ("KRT5", "KRT13", "KRT14", "EPCAM") if g in gi)
    if p > 0.6 and e < 0.5:
        imm_clusters.append(k)
    if e > 0.8 and p < 0.3:
        epi_clusters.append(k)
print("immune clusters:", imm_clusters, "| epithelial clusters:", epi_clusters)

if imm_clusters:
    rows = np.where(np.isin(ad.obs.leiden.values, imm_clusters))[0]
    print(f"\n--- inside {len(rows)} immune cells: detection rate before -> after ---")
    for g in ["KRT5", "KRT13", "KRT14", "EPCAM", "PTPRC", "CD3D", "NKG7", "GNLY"]:
        if g not in gi:
            continue
        b, a = frac_pos(raw, gi[g], rows), frac_pos(Xd, gi[g], rows)
        tag = "FOREIGN(should drop)" if g in ("KRT5", "KRT13", "KRT14", "EPCAM") else "own (should hold)"
        print(f"  {g:8s} {b:6.3f} -> {a:6.3f}   {tag}")

# --- TRBC1 reliability check (task step 2) ---------------------------------
print("\n--- TRBC1 vs TRBC2 in CD3D+ cells (raw counts) ---")
if "CD3D" in gi:
    cd3 = raw[:, gi["CD3D"]].toarray().ravel() > 0
    rows = np.where(cd3)[0]
    for g in ("TRBC1", "TRBC2", "TRAC", "CD3E"):
        if g in ad.var_names:
            print(f"  {g:6s} detected in {frac_pos(raw, ad.var_names.get_loc(g), rows):6.3f} of CD3D+ cells")
        else:
            print(f"  {g:6s} NOT IN REFERENCE")
