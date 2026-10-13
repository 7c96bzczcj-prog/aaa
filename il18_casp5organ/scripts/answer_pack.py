"""Build the final answer table + answer figure from results/il18_casp1_evidence.tsv.

Why this exists: the deck/report previously annotated candidates with the DONOR-MEDIAN detection rate.
With 25-180 cells per donor, the donor median censors any true detection rate below ~1/n to exactly 0,
so populations that differ 100-fold (human lung mast cell 138/5299 positive vs mouse brain
oligodendrocyte 2/8987) both printed as "0.0%". Everything user-facing is therefore switched to the
POOLED rate = positive cells / analysed cells across all qualifying donors, with the raw counts and the
donor support shown next to it; the donor median is kept in a separate column, labelled as censoring.

No candidate is added or removed here (spec §9.5 — side-by-side presentation only). The AT1 rows come
from screen_missed_but_reference_supported.tsv and are marked as missed by the relative CASP1 rule.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import write_tsv
from evidence_integrate import pooled_stats

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"

# explicit display names (no regex, no inference)
CN = {
    ("mouse", "brain", "oligodendrocyte"): ("小鼠", "脑", "少突胶质"),
    ("mouse", "brain", "astrocyte"): ("小鼠", "脑", "星形胶质"),
    ("mouse", "brain", "central nervous system macrophage"): ("小鼠", "脑", "小胶质/CNS巨噬"),
    ("mouse", "brain", "Cajal-Retzius cell"): ("小鼠", "脑", "Cajal-Retzius 神经元"),
    ("mouse", "lung", "alveolar macrophage"): ("小鼠", "肺", "肺泡巨噬"),
    ("mouse", "pancreas", "pancreatic D cell"): ("小鼠", "胰", "胰岛 D(δ) 细胞"),
    ("mouse", "pancreas", "pancreatic ductal cell"): ("小鼠", "胰", "胰导管细胞"),
    ("human", "lung", "mast cell"): ("人", "肺", "肥大细胞"),
    ("human", "lung", "dendritic cell"): ("人", "肺", "树突状细胞"),
    ("human", "lung", "ciliated cell"): ("人", "肺", "纤毛上皮"),
    ("human", "lung", "multiciliated columnar cell of tracheobronchial tree"): ("人", "肺", "气道多纤毛上皮"),
    ("human", "lung", "pulmonary alveolar type 1 cell"): ("人", "肺", "AT1 肺泡上皮"),
    ("human", "liver", "basophil"): ("人", "肝", "嗜碱粒细胞"),
}
SUSP = {"cell": "单细胞", "nucleus": "单核"}
TIER_CN = {
    "A_supported": "A 多源一致支持",
    "B_potential": "B 潜在可能",
    "B_potential_screen_missed": "B 潜在（相对规则漏掉）",
    "C_casp1_present": "C CASP1 确实存在",
    "D_refuted": "D 方向相反，排除",
    "E_no_matched_reference": "E 仍无匹配参考",
}
TIER_COL = {"A_supported": "#1B7A43", "B_potential": "#1F5FA8", "B_potential_screen_missed": "#4C87C4",
            "C_casp1_present": "#B07A00", "D_refuted": "#9A2B2B", "E_no_matched_reference": "#8A8A8A"}
TIER_ORDER = list(TIER_CN)


def name_of(sp, ti, ct, dataset_id):
    sp_cn, ti_cn, ct_cn = CN[(sp, ti, ct)]
    if (sp, ti, ct) == ("mouse", "pancreas", "pancreatic D cell"):
        ct_cn += "（全胰）" if dataset_id == "M_pancreas_cell" else "（分选胰岛）"
    return sp_cn, ti_cn, ct_cn


def main():
    ev = pd.read_csv(ROOT / "results" / "il18_casp1_evidence.tsv", sep="\t", comment="#")
    ev = ev[ev.duplicate_of.isna() | (ev.duplicate_of == "")].copy()
    cmap = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#")

    rows = []
    for _, r in ev.iterrows():
        sp_cn, ti_cn, ct_cn = name_of(r.species, r.tissue, r.celltype, r.dataset_id)
        rows.append({"tier": r.evidence_tier, "species": r.species, "tissue": r.tissue,
                     "suspension_type": r.suspension_type, "celltype": r.celltype,
                     "sp_cn": sp_cn, "tissue_cn": ti_cn, "celltype_cn": ct_cn,
                     "dataset_id": r.dataset_id, "level": r.level,
                     "n_cells_analysed": r.n_cells_analysed, "n_donors_analysed": r.n_donors_analysed,
                     "median_cells_per_donor": r.median_cells_per_donor,
                     "n_cells_IL18_pos": r.n_cells_IL18_pos, "pooled_pct_IL18": r.pooled_pct_IL18,
                     "donors_with_IL18_pos": r.donors_with_IL18_pos, "donor_median_pct_IL18": r.sc_pct_IL18,
                     "n_cells_CASP1_pos": r.n_cells_CASP1_pos, "pooled_pct_CASP1": r.pooled_pct_CASP1,
                     "donors_with_CASP1_pos": r.donors_with_CASP1_pos, "donor_median_pct_CASP1": r.sc_pct_CASP1,
                     "ref_ratio_IL18_over_CASP1": r.best_ratio_IL18_over_CASP1,
                     "n_ref_sources": r.n_ref_sources, "tier_reason": r.tier_reason,
                     "screen_status": "进入候选"})

    # AT1: passed IL18 but was excluded by the within-stratum relative CASP1 rule (§8), reference-supported
    miss = pd.read_csv(ROOT / "results" / "screen_missed_but_reference_supported.tsv", sep="\t", comment="#")
    for _, r in miss.iterrows():
        ps = pooled_stats(r.dataset_id, "author_label", r.celltype, cmap)
        sp_cn, ti_cn, ct_cn = name_of(r.species, r.tissue, r.celltype, r.dataset_id)
        rows.append({"tier": "B_potential_screen_missed", "species": r.species, "tissue": r.tissue,
                     "suspension_type": r.suspension_type, "celltype": r.celltype,
                     "sp_cn": sp_cn, "tissue_cn": ti_cn, "celltype_cn": ct_cn,
                     "dataset_id": r.dataset_id, "level": "author_label",
                     **{k: ps[k] for k in ("n_cells_analysed", "n_donors_analysed", "median_cells_per_donor",
                                           "n_cells_IL18_pos", "pooled_pct_IL18", "donors_with_IL18_pos",
                                           "n_cells_CASP1_pos", "pooled_pct_CASP1", "donors_with_CASP1_pos")},
                     "donor_median_pct_IL18": r.pct_IL18, "donor_median_pct_CASP1": r.pct_CASP1,
                     "ref_ratio_IL18_over_CASP1": r.best_ref.split("ratio=")[-1],
                     "n_ref_sources": r.n_ref_sources_supporting,
                     "tier_reason": f"IL18 分层内排名 {int(r.rank_IL18)}/{int(r.n_ct_stratum)}，"
                                    f"参考层支持 CASP1 低；{r.excluded_because}",
                     "screen_status": "未进候选（CASP1 未进最低三分位）"})

    out = pd.DataFrame(rows)
    out["tier_rank"] = out.tier.map({t: i for i, t in enumerate(TIER_ORDER)})
    out = out.sort_values(["tier_rank", "pooled_pct_IL18"], ascending=[True, False]).drop(columns="tier_rank")
    out["tier_cn"] = out.tier.map(TIER_CN)
    write_tsv(out, ROOT / "results" / "answer_table.tsv")

    # ---- figure: two aligned panels, pooled rates, counts on the bars -------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "WenQuanYi Zen Hei", "axes.unicode_minus": False,
                         "axes.spines.top": False, "axes.spines.right": False})

    df = out.iloc[::-1].reset_index(drop=True)      # bottom-up so A ends up on top
    n = len(df)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.2, 0.40 * n + 1.9), sharey=True,
                                   gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.34})
    ypos = np.arange(n)
    cols = [TIER_COL[t] for t in df.tier]
    labels = [f"{r.sp_cn}·{r.tissue_cn} {r.celltype_cn}\n{SUSP[r.suspension_type]}" for r in df.itertuples()]

    ax1.barh(ypos, df.pooled_pct_IL18 * 100, color=cols, height=0.58)
    for i, r in enumerate(df.itertuples()):
        ax1.text(r.pooled_pct_IL18 * 100 + 0.8, i,
                 f"{r.pooled_pct_IL18*100:.1f}%  ({r.n_cells_IL18_pos}/{r.n_cells_analysed} 细胞, "
                 f"{r.donors_with_IL18_pos}/{r.n_donors_analysed} 供体)",
                 va="center", fontsize=9.2, color="#333333")
    ax1.set_yticks(ypos, labels, fontsize=9.6)
    ax1.set_xlim(0, 160)
    ax1.set_xticks([0, 20, 40, 60, 80, 100], ["0", "20", "40", "60", "80", "100%"], fontsize=9)
    ax1.set_xlabel("IL18 检出率 = 原始计数≥1 的细胞数 ÷ 分析细胞数（合并全部供体）", fontsize=10.2)
    ax1.set_title("IL18（pro-IL-18 转录本）", fontsize=13, loc="left")

    pc = df.pooled_pct_CASP1 * 100
    ax2.barh(ypos, pc.clip(lower=1e-3), color=cols, height=0.58)
    ax2.set_xscale("symlog", linthresh=0.01, linscale=0.4)
    ax2.set_xlim(0, 30000)
    ax2.set_xticks([0, 0.01, 0.1, 1, 10], ["0", "0.01%", "0.1%", "1%", "10%"], fontsize=9)
    for i, r in enumerate(df.itertuples()):
        v = r.pooled_pct_CASP1 * 100
        txt = (f"{v:.2f}%  ({r.n_cells_CASP1_pos}/{r.n_cells_analysed} 细胞, "
               f"{r.donors_with_CASP1_pos}/{r.n_donors_analysed} 供体)" if r.n_cells_CASP1_pos
               else f"0  (0/{r.n_cells_analysed} 细胞)")
        if r.donor_median_pct_CASP1 == 0 and r.n_cells_CASP1_pos:
            txt += "，供体中位数=0"
        ax2.text(max(v, 1e-3) * 1.7 if v > 0 else 0.0017, i, txt, va="center", fontsize=9.2, color="#333333")
    ax2.set_xlabel("CASP1 检出率（同一分母；对数轴，故 0.02% 与 2.6% 可区分）", fontsize=10.2)
    ax2.set_title("CASP1", fontsize=13, loc="left")

    for ax in (ax1, ax2):
        ax.set_ylim(-0.7, n - 0.3)
        ax.grid(axis="x", color="#E6E6E6", lw=0.7)
        ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=TIER_COL[t]) for t in TIER_ORDER if t in set(df.tier)]
    names = [TIER_CN[t] for t in TIER_ORDER if t in set(df.tier)]
    ax1.legend(handles, names, loc="lower right", frameon=False, fontsize=9.6)
    fig.suptitle("IL18 高 / CASP1 低 候选：合并检出率与证据级别", fontsize=15, x=0.055, ha="left",
                 y=0.995, weight="bold")
    fig.text(0.055, 0.012,
             "分母 = 该细胞类型在该分层（数据集×器官×单细胞/单核）内、每位供体 ≥20 个细胞后进入分析的细胞总数；"
             "分子 = 其中该基因原始计数 ≥1 的细胞数。不是表达量占比、不是相对最高值的百分比、不是读段占比。",
             fontsize=8.8, color="#666666")
    fig.subplots_adjust(left=0.158, right=0.997, top=0.885, bottom=0.125)
    fig.savefig(FIG / "answer_v2.png", dpi=170)
    print("wrote figures/answer_v2.png", fig.get_size_inches() * 170)

    # ---- json for deck/report ----------------------------------------------------------------------
    js = [{"tier": r.tier, "tier_cn": TIER_CN[r.tier], "sp": r.sp_cn, "organ": r.tissue_cn,
           "susp": SUSP[r.suspension_type], "name": r.celltype_cn,
           "il18": round(r.pooled_pct_IL18 * 100, 2), "il18_pos": int(r.n_cells_IL18_pos),
           "casp1": round(r.pooled_pct_CASP1 * 100, 3), "casp1_pos": int(r.n_cells_CASP1_pos),
           "n_cells": int(r.n_cells_analysed), "n_donors": int(r.n_donors_analysed),
           "casp1_donors": int(r.donors_with_CASP1_pos),
           "med_casp1": round(r.donor_median_pct_CASP1 * 100, 2),
           "ratio": ("" if pd.isna(r.ref_ratio_IL18_over_CASP1) else str(r.ref_ratio_IL18_over_CASP1)),
           "src": int(r.n_ref_sources),
           "why": r.tier_reason, "status": r.screen_status} for r in out.itertuples()]
    (ROOT / "results" / "answer_table.json").write_text(json.dumps(js, ensure_ascii=False, indent=1))
    print(f"wrote results/answer_table.tsv + .json ({len(js)} rows)")


if __name__ == "__main__":
    main()
