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


ZERO_SLOT = 0.0018   # where a true zero (0 positive cells) is drawn on the log axis
XMIN, XMAX = 0.0011, 150.0


def draw_answer_figure(out):
    """One figure, one shared log detection-rate axis, both genes on the SAME row.

    Why this shape: the question is "IL18 high AND CASP1 low", so the two genes have to be readable
    against each other per cell type, not in two separate panels. A filled dot (IL18) and a hollow dot
    (CASP1) sit on one log axis; the segment between them IS the gap. The log axis is required — CASP1
    spans 0.02%-13.3%, and on a linear axis everything under ~1% collapses onto the spine, which is
    exactly the "it all reads as 0" failure this figure exists to undo.

    Species go in separate stacked panels, never interleaved: human/mouse and scRNA/snRNA are not
    quantitatively comparable, so the figure must not invite a cross-species ranking.

    The grey tick on each row is 1/n — the smallest non-zero rate that row could possibly show. A CASP1
    dot sitting on it, or a true zero in a row whose 1/n is high, is weakly powered, and the reader can
    see that without being told.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.family": "WenQuanYi Zen Hei", "axes.unicode_minus": False,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False})

    groups = [("人", out[out.species == "human"]), ("小鼠", out[out.species == "mouse"])]
    counts = [len(g) for _, g in groups]
    n_all = sum(counts)
    fig, axes = plt.subplots(2, 1, figsize=(17.0, 0.40 * n_all + 2.2),
                             gridspec_kw={"height_ratios": counts, "hspace": 0.26})

    for ax, (sp, g) in zip(axes, groups):
        d = g.iloc[::-1].reset_index(drop=True)      # bottom-up so the A tier lands on top
        n = len(d)
        for i, r in enumerate(d.itertuples()):
            col = TIER_COL[r.tier]
            il18 = r.pooled_pct_IL18 * 100
            casp1 = r.pooled_pct_CASP1 * 100
            x_c = casp1 if r.n_cells_CASP1_pos else ZERO_SLOT
            ax.plot([x_c, il18], [i, i], color=col, lw=3.0, alpha=.22, solid_capstyle="round", zorder=1)
            ax.vlines(100.0 / r.n_cells_analysed, i - .26, i + .26, color="#B5B5B5", lw=1.0, zorder=2)
            ax.scatter([x_c], [i], s=66, facecolors="white", edgecolors=col, linewidths=2.1, zorder=4)
            ax.scatter([il18], [i], s=78, color=col, zorder=4)
            ax.text(il18 * 1.22, i, f"{il18:.1f}%", va="center", fontsize=9.8, color=col, weight="bold")
            ax.text(x_c / 1.28, i, "0" if not r.n_cells_CASP1_pos else f"{casp1:.2f}%",
                    va="center", ha="right", fontsize=9.3, color=col)
            ax.text(1.006, i, f"{r.n_cells_CASP1_pos} · {r.n_cells_IL18_pos} / {r.n_cells_analysed}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left",
                    fontsize=8.9, color="#555555")
        ax.set_xscale("log")
        ax.set_xlim(XMIN, XMAX)
        ax.set_yticks(np.arange(n), [f"{r.tissue_cn} {r.celltype_cn}（{SUSP[r.suspension_type]}）"
                                     for r in d.itertuples()], fontsize=9.9)
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.72, n - 0.28)
        ax.grid(axis="x", color="#ECECEC", lw=.7)
        ax.set_axisbelow(True)
        ax.text(-0.235, 1.0, sp, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=15, weight="bold", color="#1A1A1A")
        ax.text(1.006, 1.0, "CASP1 · IL18 阳性细胞 / 分析细胞", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=8.9, color="#555555")
        ax.set_xticks([0.01, 0.1, 1, 10, 100], ["0.01%", "0.1%", "1%", "10%", "100%"], fontsize=9.6)

    axes[0].tick_params(axis="x", labelbottom=False, length=0)
    axes[1].set_xlabel("检出率（对数轴）= 该基因原始计数 ≥1 的细胞数 ÷ 分析细胞数（全部供体合并）", fontsize=10.6)

    marks = [Line2D([], [], marker="o", ls="", ms=8.6, color="#444"),
             Line2D([], [], marker="o", ls="", ms=8.6, mfc="white", mec="#444", mew=2.1),
             Line2D([], [], marker="|", ls="", ms=9, mew=1.4, color="#B5B5B5")]
    tiers = [t for t in TIER_ORDER if t in set(out.tier)]
    handles = marks + [Line2D([], [], marker="s", ls="", ms=8, color=TIER_COL[t]) for t in tiers]
    names = ["IL18（pro-IL-18 转录本）", "CASP1", "该行可检出的最低非零值 1/n"] + [TIER_CN[t] for t in tiers]
    axes[0].legend(handles, names, loc="lower left", bbox_to_anchor=(-0.235, 1.07), ncol=3,
                   frameon=False, fontsize=9.5, handletextpad=.5, columnspacing=1.8)
    fig.text(0.012, 0.012,
             "分母 = 该细胞类型在该分层（数据集 × 器官 × 单细胞/单核）内、每位供体 ≥20 个细胞后进入分析的细胞总数；"
             "分子 = 其中该基因原始计数 ≥1 的细胞数。不是表达量占比、不是相对最高值的百分比、不是读段占比。"
             "人与小鼠、单细胞与单核之间不作定量比较。",
             fontsize=8.9, color="#666666")
    fig.subplots_adjust(left=0.185, right=0.845, top=0.862, bottom=0.105)
    fig.savefig(FIG / "answer_v2.png", dpi=170)
    print("wrote figures/answer_v2.png", fig.get_size_inches() * 170)


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

    draw_answer_figure(out)

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
