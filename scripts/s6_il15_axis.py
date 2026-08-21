#!/usr/bin/env python3
"""
S-6 -- IL-15 supply, receptor level and STAT5 activity across organs.

TARGETED READ: 8 columns of ~35k (IL15, IL15RA, IL2RB, IL2RG, CISH, BCL2,
SOCS1, PIM1). Everything already cached (CCL3/4/5, XCL1/2, covariates) is
read from the parquet, not re-extracted. Runtime and column count are
recorded, per the standing rule.

SCOPE LIMIT, established before any number is computed and carried into
the output: the object on disk is the atlas's "T & innate lymphoid cells"
subset. Its 18 author labels are all T, NK or ILC3 -- there are NO
myeloid, dendritic, stromal or epithelial cells in it. IL-15 is chiefly
made and trans-presented by exactly those absent populations, so §1.1's
"non-NK cells of each organ" can only mean "the T cells and ILC3 of each
organ" here. That is a different quantity from an organ's IL-15 supply,
and is labelled SCOPE_LIMITED throughout. §6 forbids a new download, so
the myeloid subset is not fetched.

§1.2 and §1.3 are NK-intrinsic and are unaffected by this limit.
"""
import time
import warnings
import numpy as np
import pandas as pd
import h5py
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

H5 = "data/raw/tissue_immune_T_ILC.h5ad"
OUT = "results"
NEW = ["IL15", "IL15RA", "IL2RB", "IL2RG", "CISH", "BCL2", "SOCS1", "PIM1"]
NKG = ["IL2RB", "IL2RG", "CISH", "BCL2", "SOCS1", "PIM1"]
CHEMO = ["CCL3", "CCL4", "CCL5", "XCL1"]
NKLAB = ["NK_CD16+", "NK_CD56bright_CD16-"]
MIN_DON, MIN_NK = 3, 200
FLOOR, CEIL = 0.05, 0.85
CHUNK = 20000
NBOOT = 2000
RNG = np.random.default_rng(20260821)


def band(v):
    if not np.isfinite(v):
        return "NA"
    return "floor" if v < FLOOR else ("ceiling" if v > CEIL else "measurable")


def main():
    t0 = time.time()
    cache = pd.read_parquet(f"{OUT}/organ_cells_cache.parquet")
    n = len(cache)
    have = [g for g in NEW if g in cache.columns]
    need = [g for g in NEW if g not in cache.columns]
    print(f"cached cells {n}; already cached {have}; to extract {need}")

    if need:
        f = h5py.File(H5, "r")
        rv = f["raw/var"]
        key = "feature_name" if "feature_name" in rv else "gene_symbols"
        g = rv[key]
        if isinstance(g, h5py.Group):
            cats = [x.decode() if isinstance(x, bytes) else x
                    for x in g["categories"][:]]
            names = np.array(cats)[g["codes"][:]]
        else:
            names = np.array([x.decode() if isinstance(x, bytes) else x
                              for x in g[:]])
        idx = {}
        for gene in need:
            hit = np.where(names == gene)[0]
            idx[gene] = int(hit[0]) if len(hit) else None
        print(f"TARGETED SCAN: {len(need)} columns of {len(names)} genes")
        print(f"  gene -> column: {idx}")
        missing = [k for k, v in idx.items() if v is None]
        if missing:
            print(f"  NOT_IN_MATRIX: {missing}")
        cols = {v: k for k, v in idx.items() if v is not None}
        colset = set(cols)
        X = f["raw/X"]
        indptr = X["indptr"][:]
        data_ds, ind_ds = X["data"], X["indices"]
        counts = {gg: np.zeros(n, dtype=np.float32) for gg in cols.values()}
        for start in range(0, n, CHUNK):
            end = min(start + CHUNK, n)
            lo, hi = int(indptr[start]), int(indptr[end])
            d = data_ds[lo:hi]
            j = ind_ds[lo:hi]
            sel = np.isin(j, list(colset))
            if sel.any():
                rows = np.repeat(np.arange(start, end),
                                 np.diff(indptr[start:end + 1]).astype(int))
                for r, jj, vv in zip(rows[sel], j[sel], d[sel]):
                    counts[cols[int(jj)]][r] = vv
        for gene in need:
            cache[gene] = counts.get(gene, np.zeros(n))
        cache.to_parquet(f"{OUT}/organ_cells_cache.parquet")
    scan_s = time.time() - t0
    print(f"SCAN COMPLETE in {scan_s:.0f} s "
          f"({len(need)} columns extracted, rest from cache)")

    ct = cache.author_ct.astype(str)
    isnk = ct.isin(NKLAB).values
    print(f"\nlabels present: {sorted(ct.unique())}")
    print(f"  => all labels are T / NK / ILC3. NO myeloid, DC, stromal or "
          f"epithelial cells in this object.\n     §1.1 is therefore "
          f"SCOPE_LIMITED (see the module docstring).")

    # ---------------- donor-level tables --------------------------
    rows = []
    for (tis, asy, don), gi in cache.groupby(["tissue", "assay", "donor"],
                                             observed=True).indices.items():
        sub = cache.iloc[gi]
        nk = isnk[gi]
        rec = dict(tissue=tis, assay=asy, donor=don,
                   n_nk=int(nk.sum()), n_nonnk=int((~nk).sum()))
        for gg in ("IL15", "IL15RA"):
            v = sub[gg].values[~nk]
            rec[gg + "_nonNK"] = float((v > 0).mean()) if len(v) else np.nan
        for gg in NKG + CHEMO:
            v = sub[gg].values[nk]
            rec[gg + "_NK"] = float((v > 0).mean()) if len(v) else np.nan
        rows.append(rec)
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/s6_per_donor.csv", index=False)

    aggc = {"donors": ("donor", "nunique"), "n_nk": ("n_nk", "sum"),
            "n_nonnk": ("n_nonnk", "sum")}
    val = [c for c in D.columns if c.endswith("_nonNK") or c.endswith("_NK")]
    A = D.groupby(["tissue", "assay"], observed=True).agg(
        **aggc, **{c: (c, "median") for c in val}).reset_index()
    Q1 = D.groupby(["tissue", "assay"], observed=True)[val].quantile(.25)
    Q3 = D.groupby(["tissue", "assay"], observed=True)[val].quantile(.75)
    A["reportable"] = (A.donors >= MIN_DON) & (A.n_nk >= MIN_NK)
    A.to_csv(f"{OUT}/s6_by_organ.csv", index=False)
    Q1.to_csv(f"{OUT}/s6_q1.csv"); Q3.to_csv(f"{OUT}/s6_q3.csv")

    print("\n" + "=" * 78)
    print("MAIN TABLE (donor medians; INSUFFICIENT = donors < 3 or NK < 200)")
    print("=" * 78)
    show = ["tissue", "assay", "donors", "n_nk", "n_nonnk",
            "IL15_nonNK", "IL15RA_nonNK"] + [g + "_NK" for g in NKG + CHEMO]
    print(A[show].round(4).to_string(index=False))
    print(f"\nreportable strata: {int(A.reportable.sum())} / {len(A)}")

    # ---------------- who provides IL-15 --------------------------
    print("\n" + "=" * 78)
    print("§1.1 cell-type composition of IL15+ / IL15RA+ NON-NK cells")
    print("  (SCOPE_LIMITED: T and ILC3 only -- see above)")
    print("=" * 78)
    comp_rows = []
    ok = A[A.reportable]
    for _, r in ok.iterrows():
        m = ((cache.tissue.astype(str) == r.tissue) &
             (cache.assay.astype(str) == r.assay) & (~isnk))
        sub = cache[m]
        for gg in ("IL15", "IL15RA"):
            pos = sub[sub[gg].values > 0]
            top = pos.author_ct.astype(str).value_counts(normalize=True).head(5)
            comp_rows.append(dict(tissue=r.tissue, assay=r.assay, gene=gg,
                                  n_pos=len(pos), n_nonnk=len(sub),
                                  top5="; ".join(f"{k} {v:.1%}"
                                                 for k, v in top.items())))
    C = pd.DataFrame(comp_rows)
    C.to_csv(f"{OUT}/s6_provider_composition.csv", index=False)
    for _, r in C.iterrows():
        print(f"  {r.tissue[:22]:22s} {r.assay:10s} {r.gene:7s} "
              f"n+={r.n_pos:5d}/{r.n_nonnk:6d}  {r.top5}")

    # ---------------- rank correlations ---------------------------
    print("\n" + "=" * 78)
    print("§1.3 cross-organ rank correlations (10x 5' v1, author annotation)")
    print("  donor bootstrap within organ, 2000 resamples")
    print("=" * 78)
    ASY = "10x 5' v1"
    s = ok[ok.assay == ASY]
    organs = sorted(s.tissue.astype(str))
    print(f"  organs in the reportable 5' v1 stratum ({len(organs)}): {organs}")
    pool = D[(D.assay == ASY) & (D.tissue.astype(str).isin(organs))]
    pairs = [("CISH_NK", c + "_NK") for c in CHEMO] + \
            [("IL15_nonNK", c + "_NK") for c in CHEMO] + \
            [("IL15RA_nonNK", c + "_NK") for c in CHEMO] + \
            [("IL2RB_NK", c + "_NK") for c in CHEMO] + \
            [(t + "_NK", c + "_NK") for t in ("BCL2", "SOCS1", "PIM1")
             for c in CHEMO]
    out = []
    for a, b in pairs:
        obs = stats.spearmanr(s[a], s[b]).statistic
        bs = []
        for _ in range(NBOOT):
            xa, xb = [], []
            for o in organs:
                gg = pool[pool.tissue.astype(str) == o]
                take = gg.iloc[RNG.integers(0, len(gg), len(gg))]
                xa.append(np.nanmedian(take[a].values))
                xb.append(np.nanmedian(take[b].values))
            if np.isfinite(xa).all() and np.isfinite(xb).all():
                bs.append(stats.spearmanr(xa, xb).statistic)
        bs = np.array([v for v in bs if np.isfinite(v)])
        lo, hi = (np.percentile(bs, 2.5), np.percentile(bs, 97.5)) if len(bs) else (np.nan, np.nan)
        sig = "excludes 0" if (lo > 0 or hi < 0) else "includes 0"
        print(f"  rho({a:14s}, {b:9s}) = {obs:+.3f}   "
              f"95% [{lo:+.3f}, {hi:+.3f}]   {sig}")
        out.append(dict(x=a, y=b, rho=obs, lo=lo, hi=hi, excludes_zero=(lo > 0 or hi < 0)))
    R = pd.DataFrame(out)
    R.to_csv(f"{OUT}/s6_correlations.csv", index=False)

    # ---------------- bands ---------------------------------------
    print("\n" + "=" * 78)
    print("MEASURABILITY BANDS (reportable strata only)")
    print("=" * 78)
    for gg in ["IL15_nonNK", "IL15RA_nonNK"] + [x + "_NK" for x in NKG + CHEMO]:
        v = ok[gg].values
        b = [band(x) for x in v]
        print(f"  {gg:14s} median {np.nanmedian(v):.3f}  "
              f"range {np.nanmin(v):.3f}-{np.nanmax(v):.3f}  "
              f"floor {b.count('floor')}/{len(b)}  "
              f"measurable {b.count('measurable')}/{len(b)}  "
              f"ceiling {b.count('ceiling')}/{len(b)}")
    print(f"\nRUNTIME: targeted scan {scan_s:.0f} s, {len(need)} columns")
    print("wrote s6_per_donor.csv, s6_by_organ.csv, s6_provider_composition.csv, "
          "s6_correlations.csv")


if __name__ == "__main__":
    main()
