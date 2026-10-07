"""Phase 2 step 2: main table, donor-level summary, within-study comparisons (spec §6, §7).

Statistical unit is the donor. Cells -> (dataset, condition, donor, unified celltype) pseudobulk.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from common import MIN_CELLS, SEED, write_tsv

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
UNIFIED = ["alveolar_macrophage", "monocyte_derived_macrophage", "dendritic_cell", "neutrophil", "T_NK",
           "B_plasma", "AT1", "AT2", "airway_epithelial", "endothelial", "fibroblast_mesenchymal"]
EPITHELIAL = {"AT1", "AT2", "airway_epithelial"}
PLATE_BASED = {"smartseq2"}


def load_map():
    m = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#", dtype=str)
    bad = set(m.celltype) - set(UNIFIED) - {"other"}
    if bad:
        raise ValueError(f"celltype_map.tsv has non-unified labels: {bad}")
    if m.duplicated(["entry_id", "celltype_original"]).any():
        raise ValueError("celltype_map.tsv has duplicate (entry_id, celltype_original) rows")
    return m


def pseudobulk(df):
    g = df.groupby(["dataset_id", "condition", "condition_detail", "donor_id", "celltype", "celltype_key"],
                   observed=True)
    out = g.agg(n_cells=("il18", "size"), n_cells_IL18_pos=("il18", lambda x: int((x >= 1).sum())),
                il18_sum=("il18", "sum"), lib_sum=("lib", "sum"),
                celltype_original=("celltype_original", lambda x: "|".join(sorted(set(x))))).reset_index()
    out["pct_expressing"] = out.n_cells_IL18_pos / out.n_cells
    out["mean_logcpm"] = np.log2(out.il18_sum / out.lib_sum * 1e6 + 1)
    small = out.n_cells < MIN_CELLS
    out.loc[small, ["pct_expressing", "mean_logcpm"]] = np.nan
    return out.drop(columns=["il18_sum", "lib_sum"])


def main():
    np.random.seed(SEED)
    spec = {e["id"]: e for e in yaml.safe_load(open(ROOT / "datasets.yaml"))["datasets"]}
    cmap = load_map()
    rows, log = [], []
    for eid, e in spec.items():
        path = ROOT / "data" / "cells" / f"{eid}.parquet"
        if not path.exists():
            log.append(f"{eid}: no per-cell table, not aggregated")
            continue
        c = pd.read_parquet(path)
        m = cmap[cmap.entry_id == eid].set_index("celltype_original")["celltype"]
        missing = sorted(set(c.celltype_original) - set(m.index))
        if missing:
            raise ValueError(f"{eid}: labels missing from celltype_map.tsv: {missing}")
        c["celltype"] = c.celltype_original.map(m)
        # 'other' keeps its original label as its own row; unified types pool their author sub-labels
        c["celltype_key"] = np.where(c.celltype == "other", "other:" + c.celltype_original, c.celltype)
        c["dataset_id"] = eid + np.where(c.sub_dataset != "", ":" + c.sub_dataset, "")
        pb = pseudobulk(c)
        pb["accession"], pb["source_layer"], pb["species"] = e["accession"], e["source_layer"], e["species"]
        pb["tissue"], pb["platform"] = e["tissue"], e["platform"]
        pb["ambient_uncorrected"] = not (eid == "H03_TabulaSapiens_lung" or e["platform"] in PLATE_BASED)
        pb["normalized_only"] = False
        pb["gene_id_unavailable"] = bool(c.gene_id_unavailable.any())
        n_small = int((pb.n_cells < MIN_CELLS).sum())
        log.append(f"{eid}: {len(c)} cells, {len(pb)} donor x celltype cells, {n_small} set to NA (n_cells < {MIN_CELLS})")
        rows.append(pb)
    main_tab = pd.concat(rows, ignore_index=True)

    # rank among the 11 unified types, within dataset x condition x time point x donor (1 = highest)
    uni = main_tab.celltype != "other"
    main_tab["rank_within_dataset"] = np.nan
    main_tab.loc[uni, "rank_within_dataset"] = (
        main_tab[uni].groupby(["dataset_id", "condition", "condition_detail", "donor_id"])["mean_logcpm"]
        .rank(ascending=False, method="min"))

    cols = ["dataset_id", "accession", "source_layer", "species", "condition", "condition_detail", "tissue",
            "platform", "donor_id", "celltype", "celltype_original", "n_cells", "n_cells_IL18_pos",
            "pct_expressing", "mean_logcpm", "rank_within_dataset", "ambient_uncorrected", "normalized_only",
            "gene_id_unavailable"]
    main_tab = main_tab.sort_values(["species", "dataset_id", "condition", "condition_detail", "donor_id", "celltype"])
    write_tsv(main_tab[cols + ["celltype_key"]].rename(columns={"celltype_key": "celltype_row_key"}),
              RES / "il18_celltype_table.tsv")

    # summary: donor as unit; NA cells excluded
    v = main_tab.dropna(subset=["mean_logcpm"])
    q = lambda p: (lambda x: np.nanpercentile(x, p))
    summ = v.groupby(["species", "dataset_id", "tissue", "platform", "condition", "condition_detail", "celltype",
                      "celltype_key", "ambient_uncorrected"], observed=True).agg(
        n_donors=("donor_id", "nunique"), n_cells_total=("n_cells", "sum"),
        pct_median=("pct_expressing", "median"), pct_q1=("pct_expressing", q(25)), pct_q3=("pct_expressing", q(75)),
        logcpm_median=("mean_logcpm", "median"), logcpm_q1=("mean_logcpm", q(25)), logcpm_q3=("mean_logcpm", q(75)),
        rank_median=("rank_within_dataset", "median"),
        celltype_original=("celltype_original", lambda x: "|".join(sorted(set("|".join(x).split("|")))))).reset_index()
    summ["epithelial_excluded_from_conclusions"] = summ.celltype.isin(EPITHELIAL) & summ.ambient_uncorrected
    summ = summ.rename(columns={"celltype_key": "celltype_row_key"})
    write_tsv(summ.sort_values(["species", "dataset_id", "condition", "condition_detail", "rank_median"]),
              RES / "il18_summary.tsv")

    # within-study resting vs condition; only where both arms have >= 3 donors (spec §7.3, §7.5)
    tests = []
    for (ds, ct), g in v[uni.loc[v.index]].groupby(["dataset_id", "celltype"]):
        conds = g.condition.unique()
        if "resting" not in conds:
            continue
        a = g[g.condition == "resting"]
        for cond in [x for x in conds if x != "resting"]:
            b = g[g.condition == cond]
            for arm in (a, b):  # the donor is the unit: one value per donor per arm
                if arm.donor_id.duplicated().any():
                    raise ValueError(f"{ds}/{ct}: donor repeated within one arm")
            if a.donor_id.nunique() < 3 or b.donor_id.nunique() < 3:
                continue
            for metric in ["mean_logcpm", "pct_expressing"]:
                p = mannwhitneyu(a[metric], b[metric], alternative="two-sided").pvalue
                tests.append(dict(dataset_id=ds, celltype=ct, contrast=f"{cond} vs resting", metric=metric,
                                  n_resting=a.donor_id.nunique(), n_condition=b.donor_id.nunique(),
                                  median_resting=a[metric].median(), median_condition=b[metric].median(),
                                  p_mannwhitney=p))
    t = pd.DataFrame(tests)
    if len(t):
        t["q_bh"] = np.nan
        for _, idx in t.groupby(["dataset_id", "metric"]).groups.items():
            p = t.loc[idx, "p_mannwhitney"].values
            o = np.argsort(p)
            r = np.empty(len(p))
            r[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
            t.loc[idx, "q_bh"] = np.minimum(r, 1)
        t["epithelial_excluded_from_conclusions"] = t.celltype.isin(EPITHELIAL) & ~t.dataset_id.str.startswith("H03")
    write_tsv(t, RES / "il18_within_study_tests.tsv")
    with open(RES / "phase2_aggregate.log", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
