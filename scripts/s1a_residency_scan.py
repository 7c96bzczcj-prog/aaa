#!/usr/bin/env python3
"""
S-1 §4 — residency-axis feasibility. TARGETED 3-COLUMN SCAN.

Reads raw/X once, extracting ONLY ITGA1, CD69, ITGAE (3 columns of
~36k). NCAM1 and FCGR3A come from the existing cache, not re-read.
Runtime is measured and reported, per the updated rule.

CD69 IS REPORTED SEPARATELY and is NEVER folded into a residency score:
CD69 is dissociation-inducible, and one of this round's core contrasts is
tissue versus blood, so pooling it would build the artefact into the
stratifier.

Feasibility rule (§4): if ITGA1 or ITGAE sits below the 5% floor in most
organs, the residency axis is not usable at transcript level. Reported as
such; no threshold lowering, no gene substitution (S1a.3 forbids both).
"""
import time
import warnings
import numpy as np
import pandas as pd
import h5py

warnings.filterwarnings("ignore")
pd.set_option("display.width", 220)

H5 = "data/raw/tissue_immune_T_ILC.h5ad"
OUT = "results"
NEW_GENES = ["ITGA1", "CD69", "ITGAE"]     # exactly 3 columns
FLOOR = 0.05
CHUNK = 20000


def main():
    t0 = time.time()
    cache = pd.read_parquet(f"{OUT}/organ_cells_cache.parquet")
    n_cells = len(cache)

    f = h5py.File(H5, "r")
    rv = f["raw/var"]
    key = "feature_name" if "feature_name" in rv else "gene_symbols"
    g = rv[key]
    if isinstance(g, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        names = np.array(cats)[g["codes"][:]]
    else:
        names = np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])
    idx = {}
    for gene in NEW_GENES:
        hit = np.where(names == gene)[0]
        idx[gene] = int(hit[0]) if len(hit) else None
    print(f"targeted scan: {len(NEW_GENES)} columns of {len(names)} genes")
    print(f"gene -> column index: {idx}")
    missing = [k for k, v in idx.items() if v is None]
    if missing:
        print(f"MISSING FROM MATRIX: {missing}")
    cols = {v: k for k, v in idx.items() if v is not None}
    colset = set(cols)

    X = f["raw/X"]
    indptr = X["indptr"][:]
    data_ds, ind_ds = X["data"], X["indices"]
    counts = {gene: np.zeros(n_cells, dtype=np.float32) for gene in cols.values()}
    for start in range(0, n_cells, CHUNK):
        end = min(start + CHUNK, n_cells)
        lo, hi = int(indptr[start]), int(indptr[end])
        d = data_ds[lo:hi]
        j = ind_ds[lo:hi]
        sel = np.isin(j, list(colset))
        if sel.any():
            rows = np.repeat(np.arange(start, end),
                             np.diff(indptr[start:end + 1]).astype(int))
            for r, jj, vv in zip(rows[sel], j[sel], d[sel]):
                counts[cols[int(jj)]][r] = vv
    scan_s = time.time() - t0
    print(f"scan complete in {scan_s:.0f} s")

    for gene in NEW_GENES:
        cache[gene] = counts.get(gene, np.zeros(n_cells))
    cache.to_parquet(f"{OUT}/organ_cells_cache.parquet")

    # ---- author-NK only, per donor x tissue x assay -------------------
    nk = cache.author_ct.astype(str).isin(["NK_CD16+", "NK_CD56bright_CD16-"])
    d = cache[nk]
    REPORT = ["ITGA1", "ITGAE", "CD69", "NCAM1", "FCGR3A"]
    rows = []
    for (tis, asy, don), gi in d.groupby(["tissue", "assay", "donor"],
                                         observed=True).indices.items():
        rec = dict(tissue=tis, assay=asy, donor=don, n_nk=len(gi))
        for gene in REPORT:
            rec[gene] = float((d[gene].values[gi] > 0).mean())
        rows.append(rec)
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/s1a_residency_per_donor.csv", index=False)

    agg = (D.groupby(["tissue", "assay"], observed=True)
           .agg(donors=("donor", "nunique"), n_nk=("n_nk", "sum"),
                **{g: (g, "median") for g in REPORT}).reset_index())
    agg["reportable"] = (agg.donors >= 3) & (agg.n_nk >= 200)
    agg.to_csv(f"{OUT}/s1a_residency_by_organ.csv", index=False)

    print("\n" + "=" * 78)
    print("§4 RESIDENCY MARKERS — donor-median detection rate in author-NK")
    print("CD69 SHOWN SEPARATELY (dissociation-inducible; never pooled)")
    print("=" * 78)
    ok = agg[agg.reportable]
    print(f"{'tissue':22s} {'assay':10s} {'don':>4s} {'nNK':>6s} "
          f"{'ITGA1':>8s} {'ITGAE':>8s} | {'CD69':>8s} | {'NCAM1':>8s} {'FCGR3A':>8s}")
    for _, r in ok.sort_values(["assay", "tissue"]).iterrows():
        def mark(v):
            return f"{v:.3f}" + ("*" if v < FLOOR else " ")
        print(f"{r.tissue[:22]:22s} {r.assay:10s} {r.donors:4d} {r.n_nk:6d} "
              f"{mark(r.ITGA1):>8s} {mark(r.ITGAE):>8s} | {mark(r.CD69):>8s} | "
              f"{mark(r.NCAM1):>8s} {mark(r.FCGR3A):>8s}")
    print("   * = below the 5% floor band")

    print("\n" + "-" * 78)
    print("§4 FEASIBILITY VERDICT")
    print("-" * 78)
    for gene in ["ITGA1", "ITGAE", "CD69", "NCAM1", "FCGR3A"]:
        v = ok[gene]
        below = int((v < FLOOR).sum())
        print(f"   {gene:7s} median across organs {v.median():.3f}  "
              f"range {v.min():.3f}-{v.max():.3f}  "
              f"organs below 5% floor: {below}/{len(v)}"
              + ("   <== MOSTLY ON THE FLOOR" if below > len(v) / 2 else ""))
    a_bad = int((ok.ITGA1 < FLOOR).sum()) > len(ok) / 2
    e_bad = int((ok.ITGAE < FLOOR).sum()) > len(ok) / 2
    print(f"\n   ITGA1 mostly on floor: {a_bad};  ITGAE mostly on floor: {e_bad}")
    print("   -> residency axis at transcript level: "
          + ("NOT FEASIBLE" if (a_bad or e_bad) else "possibly feasible"))
    print(f"\nRUNTIME: targeted scan {scan_s:.0f} s, "
          f"{len(NEW_GENES)} columns extracted")

    # ---- author labels that are residency-related ---------------------
    print("\n" + "-" * 78)
    print("§4.1 residency-related labels already in the annotation?")
    print("-" * 78)
    labs = sorted(cache.author_ct.astype(str).unique())
    res = [l for l in labs
           if any(k in l.lower() for k in ("trnk", "ilc", "trm", "resident"))]
    print(f"   all labels: {labs}")
    print(f"   residency-related: {res if res else 'NONE (no trNK / tissue-resident NK label)'}")


if __name__ == "__main__":
    main()
