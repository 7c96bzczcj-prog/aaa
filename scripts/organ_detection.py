#!/usr/bin/env python3
"""
Cross-organ NK chemokine measurability check, §3-§5.

Source: Dominguez Conde et al., Science 2022, "T & innate lymphoid cells"
subset of the Tissue Immune Cell Atlas (CELLxGENE), 216,611 cells,
12 donors, 17 tissues, all `disease == normal`. Same donors across
organs, so donor variation is controlled by design.

Computes, per (tissue x assay x donor):
  - detection rate of CCL3/CCL4/CCL5/XCL1/XCL2 in NK, under BOTH
    §3 definitions (author annotation, and a uniform marker gate)
  - the same restricted to CD56bright NK  (§5 item 2)
  - CD56bright share of NK                (§5 item 3)
  - median UMI per cell, dissociation-stress score  (§4)

NOTHING is compared across organs here: no test, no fold change, no
ranking. §0 forbids it. The output is a table.

Memory: raw/X is read in row blocks, keeping only the needed gene
columns and per-cell totals, so the full matrix is never materialised.
"""
import warnings
import numpy as np
import pandas as pd
import h5py

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

H5 = "data/raw/tissue_immune_T_ILC.h5ad"
OUT = "results"

TARGETS = ["CCL3", "CCL4", "CCL5", "XCL1", "XCL2"]
GATE = ["KLRD1", "NCAM1", "CD3E", "CD3D", "FCGR3A"]
STRESS = ["FOS", "JUN", "JUNB", "EGR1", "HSPA1A", "HSPA1B"]
NEEDED = TARGETS + GATE + STRESS

CHUNK = 20000


def read_cat(f, key):
    g = f["obs"][key]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        codes = g["codes"][:]
        return pd.Categorical.from_codes(codes, cats)
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


def main():
    f = h5py.File(H5, "r")

    # ---- obs -------------------------------------------------------
    obs = pd.DataFrame({
        "donor": read_cat(f, "donor_id"),
        "tissue": read_cat(f, "tissue"),
        "assay": read_cat(f, "assay"),
        "author_ct": read_cat(f, "Manually_curated_celltype"),
        "celltypist": read_cat(f, "Majority_voting_CellTypist"),
        "cell_type": read_cat(f, "cell_type"),
    })
    n_cells = len(obs)
    print(f"cells {n_cells}, donors {obs.donor.nunique()}, "
          f"tissues {obs.tissue.nunique()}, assays {obs.assay.nunique()}")

    # ---- var (use raw/var, which matches raw/X) --------------------
    rv = f["raw/var"]
    key = "feature_name" if "feature_name" in rv else "gene_symbols"
    g = rv[key]
    if isinstance(g, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        names = np.array(cats)[g["codes"][:]]
    else:
        names = np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])
    idx = {}
    for gene in NEEDED:
        hit = np.where(names == gene)[0]
        idx[gene] = int(hit[0]) if len(hit) else None
    missing = [k for k, v in idx.items() if v is None]
    print(f"genes not in matrix: {missing if missing else 'none'}")
    cols = {v: k for k, v in idx.items() if v is not None}
    colset = set(cols)

    # ---- chunked pass over raw/X ----------------------------------
    X = f["raw/X"]
    indptr = X["indptr"][:]
    data_ds, ind_ds = X["data"], X["indices"]
    counts = {gene: np.zeros(n_cells, dtype=np.float32) for gene in cols.values()}
    total = np.zeros(n_cells, dtype=np.float64)

    print("scanning raw/X ...")
    for start in range(0, n_cells, CHUNK):
        end = min(start + CHUNK, n_cells)
        lo, hi = int(indptr[start]), int(indptr[end])
        d = data_ds[lo:hi]
        j = ind_ds[lo:hi]
        rows = np.repeat(np.arange(start, end),
                         np.diff(indptr[start:end + 1]).astype(int))
        np.add.at(total, rows, d)
        sel = np.isin(j, list(colset))
        if sel.any():
            for r, jj, vv in zip(rows[sel], j[sel], d[sel]):
                counts[cols[int(jj)]][r] = vv
        if (start // CHUNK) % 3 == 0:
            print(f"   {end}/{n_cells}")

    integer_like = np.allclose(total[:100], np.round(total[:100]))
    print(f"raw/X integer-valued (counts): {integer_like}")

    C = pd.DataFrame({g: counts[g] for g in counts})
    obs["total_umi"] = total

    # ---- §3 two NK definitions ------------------------------------
    author_nk = obs.author_ct.astype(str).isin(["NK_CD16+", "NK_CD56bright_CD16-"])
    author_bright = obs.author_ct.astype(str).eq("NK_CD56bright_CD16-")

    # uniform gate: threshold stated = raw count > 0 positive, == 0 negative
    uni_nk = ((C["KLRD1"] > 0) & (C["NCAM1"] > 0) &
              (C["CD3E"] == 0) & (C["CD3D"] == 0)).values
    uni_bright = uni_nk & (C["FCGR3A"] == 0).values

    print(f"\nauthor NK: {int(author_nk.sum())}  "
          f"(bright {int(author_bright.sum())})")
    print(f"uniform-gate NK: {int(uni_nk.sum())}  "
          f"(bright {int(uni_bright.sum())})")
    print("uniform gate threshold: KLRD1>0 & NCAM1>0 & CD3E==0 & CD3D==0 "
          "(raw counts); bright adds FCGR3A==0")

    obs["author_nk"], obs["author_bright"] = author_nk.values, author_bright.values
    obs["uni_nk"], obs["uni_bright"] = uni_nk, uni_bright
    stress_present = [g for g in STRESS if g in C.columns]
    obs["stress"] = np.log1p(
        C[stress_present].div(obs.total_umi.replace(0, np.nan), axis=0) * 1e4
    ).mean(axis=1).values

    # cache the expensive scan so aggregation can be re-run cheaply
    cache = pd.concat([obs.reset_index(drop=True),
                       C.reset_index(drop=True)], axis=1)
    cache.to_parquet(f"{OUT}/organ_cells_cache.parquet")
    print(f"cached per-cell table: {cache.shape}")

    # ---- §5 per donor x tissue x assay ----------------------------
    rows = []
    tis_a = obs.tissue.astype(str).values
    asy_a = obs.assay.astype(str).values
    don_a = obs.donor.astype(str).values
    umi_a = obs.total_umi.values
    str_a = obs.stress.values
    gene_a = {g: C[g].values for g in TARGETS if g in C.columns}

    for defn, nkmask, brmask in (
            ("author", np.asarray(obs.author_nk), np.asarray(obs.author_bright)),
            ("uniform", np.asarray(obs.uni_nk), np.asarray(obs.uni_bright))):
        keys = pd.DataFrame({"tissue": tis_a, "assay": asy_a, "donor": don_a})
        for (tis, asy, don), gi in keys[nkmask].groupby(
                ["tissue", "assay", "donor"], observed=True).indices.items():
            pos = np.flatnonzero(nkmask)[gi]      # positions in the full arrays
            n = len(pos)
            if n == 0:
                continue
            isbr = brmask[pos]
            bpos = pos[isbr]
            rec = dict(definition=defn, tissue=tis, assay=asy, donor=don,
                       n_nk=int(n), n_bright=int(len(bpos)),
                       median_umi=float(np.median(umi_a[pos])),
                       bright_frac=float(isbr.mean()),
                       stress=float(np.mean(str_a[pos])))
            for gene in TARGETS:
                if gene not in gene_a:
                    rec[gene] = np.nan; rec[gene + "_bright"] = np.nan; continue
                v = gene_a[gene]
                rec[gene] = float((v[pos] > 0).mean())
                rec[gene + "_bright"] = (float((v[bpos] > 0).mean())
                                         if len(bpos) else np.nan)
            rows.append(rec)
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/organ_per_donor.csv", index=False)
    print(f"\nwrote per-donor rows: {len(D)}")
    print(D.groupby(["definition", "tissue"], observed=True)
          .agg(donors=("donor", "nunique"), nk=("n_nk", "sum")).to_string())


if __name__ == "__main__":
    main()
