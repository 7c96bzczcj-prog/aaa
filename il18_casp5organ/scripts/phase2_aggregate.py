"""Phase 2 aggregate + candidate logic (spec §8, §9, §10).

Stratum = dataset_id × tissue × suspension_type = one datasets.yaml entry (each entry is already a single
stratum). Statistical unit is the donor. Candidate judgement runs at BOTH cell-type levels: level-1 lineage
(§7.1) and level-2 author label (§7.2). No cross-stratum, no cross-suspension, no cross-species comparison.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import GENE_ORDER, MIN_CELLS, SEED, write_tsv

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
LINEAGES = ["parenchymal_epithelial", "neuron", "glia", "resident_macrophage", "monocyte_derived_myeloid",
            "dendritic_cell", "lymphoid", "granulocyte", "endothelial", "mesenchymal", "other"]
IL18_TOP = 1 / 3            # §8.2 IL18 top tercile
CASP1_BOT = 1 / 3          # §8.2 CASP1 bottom tercile
CASP3_FLOOR = 0.20         # §8.2 CASP3 not in lowest 20%
MIN_PCT_IL18 = 0.10        # §8.4
MIN_CELLTYPES_FOR_FLOOR = 8  # §8.3


def load_map():
    m = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#", dtype=str)
    bad = set(m.lineage) - set(LINEAGES)
    if bad:
        raise ValueError(f"celltype_map lineages not in the 11 classes: {bad}")
    return m


def pseudobulk(cells, key_cols):
    """Per (stratum, donor, celltype) pseudobulk; NA when n_cells < MIN_CELLS."""
    g = cells.groupby(key_cols, observed=True)
    agg = {"n_cells": ("lib", "size"), "lib": ("lib", "sum")}
    for gene in GENE_ORDER:
        agg[f"sum_{gene}"] = (gene, "sum")
        agg[f"pos_{gene}"] = (gene, lambda x: float((x >= 1).sum()))
    out = g.agg(**agg).reset_index()
    for gene in GENE_ORDER:
        out[f"pct_{gene}"] = out[f"pos_{gene}"] / out.n_cells
        out[f"logcpm_{gene}"] = np.log2(out[f"sum_{gene}"] / out.lib * 1e6 + 1)
    small = out.n_cells < MIN_CELLS
    for gene in GENE_ORDER:
        out.loc[small, [f"pct_{gene}", f"logcpm_{gene}"]] = np.nan
    return out.drop(columns=["lib"] + [f"sum_{g}" for g in GENE_ORDER] + [f"pos_{g}" for g in GENE_ORDER])


def donor_summary(pb, celltype_col, meta):
    """Donor-median per (stratum, celltype); keep only celltypes with >= 3 donors (§10.3)."""
    rows = []
    for (eid, ct), g in pb.groupby(["entry_id", celltype_col], observed=True):
        valid = g.dropna(subset=["logcpm_IL18"])
        r = {"entry_id": eid, "celltype": ct, "n_donors": valid.donor_id.nunique(),
             "n_cells_total": int(g.n_cells.sum())}
        for gene in GENE_ORDER:
            r[f"pct_{gene}"] = valid[f"pct_{gene}"].median()
            r[f"logcpm_{gene}"] = valid[f"logcpm_{gene}"].median()
            r[f"iqr_logcpm_{gene}"] = valid[f"logcpm_{gene}"].quantile(.75) - valid[f"logcpm_{gene}"].quantile(.25)
        rows.append(r)
    s = pd.DataFrame(rows).merge(meta, on="entry_id", how="left")
    return s


def judge(s, il18_q=IL18_TOP, casp1_q=CASP1_BOT, casp3_floor=CASP3_FLOOR):
    """Add candidate flags per stratum. Terciles/quantiles computed within stratum over celltypes with
    n_donors>=3 and non-NA medians."""
    s = s.copy()
    for col in ["rank_IL18", "rank_CASP1", "rank_CASP3", "IL18_top", "CASP1_bot", "CASP3_ok",
                "detect_ok", "n_ct_stratum", "casp3_floor_skipped", "is_candidate"]:
        s[col] = np.nan
    for eid, g in s.groupby("entry_id"):
        v = g[(g.n_donors >= 3) & g.logcpm_IL18.notna()].copy()
        n = len(v)
        if n == 0:
            continue
        idx = v.index
        s.loc[idx, "n_ct_stratum"] = n
        for gene in GENE_ORDER:
            s.loc[idx, f"rank_{gene}"] = v[f"logcpm_{gene}"].rank(ascending=False, method="min")
        il18_cut = v.logcpm_IL18.quantile(1 - il18_q)
        casp1_cut = v.logcpm_CASP1.quantile(casp1_q)
        floor_skip = n < MIN_CELLTYPES_FOR_FLOOR
        casp3_cut = v.logcpm_CASP3.quantile(casp3_floor)
        s.loc[idx, "IL18_top"] = v.logcpm_IL18 >= il18_cut
        s.loc[idx, "CASP1_bot"] = v.logcpm_CASP1 <= casp1_cut
        s.loc[idx, "CASP3_ok"] = True if floor_skip else (v.logcpm_CASP3 > casp3_cut)
        s.loc[idx, "detect_ok"] = v.pct_IL18 >= MIN_PCT_IL18
        s.loc[idx, "casp3_floor_skipped"] = floor_skip
    s["is_candidate"] = (s.IL18_top.fillna(False) & s.CASP1_bot.fillna(False)
                         & s.CASP3_ok.fillna(False) & s.detect_ok.fillna(False)
                         & (s.n_donors >= 3) & (~s.gene_id_unavailable) & (~s.gene_missing)
                         & (~s.ambient_uncorrected_excluded))
    return s


def main():
    np.random.seed(SEED)
    spec = {e["id"]: e for e in yaml.safe_load(open(ROOT / "datasets.yaml"))["datasets"]}
    cmap = load_map()
    meta_rows, cells_all = [], []
    for eid, e in spec.items():
        p = ROOT / "data" / "cells" / f"{eid}.parquet"
        if not p.exists():
            continue
        c = pd.read_parquet(p)
        m = cmap[cmap.entry_id == eid].set_index("celltype_original")["lineage"]
        miss = sorted(set(c.celltype_original) - set(m.index))
        if miss:
            raise ValueError(f"{eid}: labels missing from celltype_map.tsv: {miss}")
        c["lineage"] = c.celltype_original.map(m)
        c["entry_id"] = eid
        cells_all.append(c)
        # provisional: all Census/GEO here are ambient-uncorrected (decision ①B) -> flagged, not excluded
        meta_rows.append({"entry_id": eid, "accession": e["accession"], "species": e["species"],
                          "tissue": e["tissue"], "suspension_type": e["suspension_type"],
                          "platform": e["platform"], "gene_missing": bool(c.gene_missing.any()),
                          "gene_id_unavailable": bool(c.gene_id_unavailable.any()),
                          "ambient_uncorrected": True, "ambient_uncorrected_excluded": False})
    meta = pd.DataFrame(meta_rows)
    cells = pd.concat(cells_all, ignore_index=True)

    # ---- main table: per (stratum, donor, author label), with lineage column and within-stratum ranks ----
    pb_orig = pseudobulk(cells, ["entry_id", "donor_id", "celltype_original", "lineage"])
    pb_orig = pb_orig.merge(meta, on="entry_id")
    # ranks within stratum at the author-label level
    for gene in GENE_ORDER:
        pb_orig[f"rank_{gene}"] = pb_orig.groupby("entry_id")[f"logcpm_{gene}"].rank(ascending=False, method="min")
    main_cols = (["entry_id", "accession", "species", "tissue", "suspension_type", "platform", "donor_id",
                  "lineage", "celltype_original", "n_cells"] + [f"pct_{g}" for g in GENE_ORDER]
                 + [f"logcpm_{g}" for g in GENE_ORDER] + [f"rank_{g}" for g in GENE_ORDER]
                 + ["gene_missing", "gene_id_unavailable", "ambient_uncorrected"])
    write_tsv(pb_orig.rename(columns={"entry_id": "dataset_id"})[
        [c.replace("entry_id", "dataset_id") for c in main_cols]], RES / "il18_casp_celltype_table.tsv")

    # ---- summary + candidate judgement at BOTH levels ----
    results = {}
    for level, col in [("lineage", "lineage"), ("author", "celltype_original")]:
        pb = pseudobulk(cells, ["entry_id", "donor_id", col])
        s = donor_summary(pb.rename(columns={col: "celltype"}), "celltype", meta)
        s = judge(s)
        s.insert(1, "level", level)
        results[level] = s
    summ = pd.concat(results.values(), ignore_index=True)
    scols = (["entry_id", "level", "accession", "species", "tissue", "suspension_type", "platform",
              "celltype", "n_donors", "n_cells_total"] + [f"pct_{g}" for g in GENE_ORDER]
             + [f"logcpm_{g}" for g in GENE_ORDER] + [f"iqr_logcpm_{g}" for g in GENE_ORDER]
             + [f"rank_{g}" for g in GENE_ORDER] + ["n_ct_stratum", "IL18_top", "CASP1_bot", "CASP3_ok",
             "detect_ok", "casp3_floor_skipped", "is_candidate", "gene_missing", "gene_id_unavailable",
             "ambient_uncorrected"])
    summ_out = summ.rename(columns={"entry_id": "dataset_id"})
    write_tsv(summ_out[[c.replace("entry_id", "dataset_id") for c in scols]], RES / "il18_casp_summary.tsv")

    # ---- candidates.tsv (§9.3): §8.2 + §8.4 passers, ordered by CASP3 rank high->low within stratum (§8.5)
    cand = summ[summ.is_candidate].copy()
    cand["priority_casp3_rank"] = cand["rank_CASP3"]
    cand = cand.sort_values(["species", "tissue", "suspension_type", "level", "priority_casp3_rank"])
    cand["passed"] = "IL18_top & CASP1_bot & CASP3_ok & detect_ok"
    ccols = ["dataset_id", "level", "species", "tissue", "suspension_type", "celltype", "n_donors",
             "pct_IL18", "pct_CASP1", "pct_CASP3", "logcpm_IL18", "logcpm_CASP1", "logcpm_CASP3",
             "rank_IL18", "rank_CASP1", "rank_CASP3", "casp3_floor_skipped", "priority_casp3_rank", "passed"]
    write_tsv(cand.rename(columns={"entry_id": "dataset_id"})[ccols], RES / "candidates.tsv")

    # ---- rank_pass_detection_fail.tsv (§9.4): IL18_top & CASP1_bot & CASP3_ok but pct_IL18 < 0.10 ----
    rpf = summ[(summ.IL18_top.fillna(False)) & (summ.CASP1_bot.fillna(False)) & (summ.CASP3_ok.fillna(False))
               & (~summ.detect_ok.fillna(False)) & (summ.n_donors >= 3)].copy()
    write_tsv(rpf.rename(columns={"entry_id": "dataset_id"})[
        ["dataset_id", "level", "species", "tissue", "suspension_type", "celltype", "n_donors",
         "pct_IL18", "logcpm_IL18", "logcpm_CASP1", "logcpm_CASP3"]], RES / "rank_pass_detection_fail.tsv")

    # ---- threshold sensitivity (§8.2.2, §9.6): tercile->quartile; 20%->10%,30% ----
    variants = {"baseline(1/3,1/3,20%)": (1/3, 1/3, .20), "quartile(1/4,1/4,20%)": (1/4, 1/4, .20),
                "casp3_10%": (1/3, 1/3, .10), "casp3_30%": (1/3, 1/3, .30)}
    srows = []
    for level, col in [("lineage", "lineage"), ("author", "celltype_original")]:
        pb = pseudobulk(cells, ["entry_id", "donor_id", col])
        base = donor_summary(pb.rename(columns={col: "celltype"}), "celltype", meta)
        sets = {}
        for name, (a, b, cc3) in variants.items():
            j = judge(base, a, b, cc3)
            sets[name] = set(map(tuple, j[j.is_candidate][["entry_id", "celltype"]].values))
        allc = set().union(*sets.values())
        for eid, ct in sorted(allc):
            srows.append({"level": level, "dataset_id": eid, "celltype": ct,
                          **{n: ((eid, ct) in sets[n]) for n in variants}})
    write_tsv(pd.DataFrame(srows), RES / "threshold_sensitivity.tsv")

    print(f"strata={meta.entry_id.nunique()}  main_rows={len(pb_orig)}  "
          f"candidates={len(cand)}  rank_pass_detect_fail={len(rpf)}")


if __name__ == "__main__":
    main()
