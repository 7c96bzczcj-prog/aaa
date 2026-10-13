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
from evidence_integrate import MIN_CELLS_PER_DONOR, pooled_stats

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

# Everything the reader sees collapses the six internal tiers into three verdicts. The six-way tier
# stays in the TSVs as the audit trail (it records WHY a row landed where it did), but no figure, slide
# or report table uses more than these three.
VERDICT = {
    "A_supported": "确定",                   # >= 2 independent references agree CASP1 is at/near zero
    "B_potential": "可能",                   # exactly one reference supports it
    "B_potential_screen_missed": "可能",     # reference-supported, excluded only by the relative rule
    "E_no_matched_reference": "可能",        # own data is consistent, no external reference exists yet
    "C_casp1_present": "相反",               # CASP1 is demonstrably present alongside IL18
    "D_refuted": "相反",                     # CASP1 >= IL18 in the reference layer
}
VERDICT_ORDER = ["确定", "可能", "相反"]
VERDICT_COL = {"确定": "#1B7A43", "可能": "#1F5FA8", "相反": "#9A2B2B"}
VERDICT_NOTE = {
    "确定": "≥2 个独立参考来源一致：CASP1 处于零底线或比 IL18 低 ≥20 倍",
    "可能": "单一参考来源支持，或本数据一致但尚无匹配参考",
    "相反": "参考层显示 CASP1 确实存在，或方向与前提相反",
}


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
        d = g.iloc[::-1].reset_index(drop=True)      # bottom-up so 确定 lands on top
        n = len(d)
        for i, r in enumerate(d.itertuples()):
            col = VERDICT_COL[r.verdict]
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
    vs = [v for v in VERDICT_ORDER if v in set(out.verdict)]
    handles = marks + [Line2D([], [], marker="s", ls="", ms=8, color=VERDICT_COL[v]) for v in vs]
    names = ["IL18（pro-IL-18 转录本）", "CASP1", "该行可检出的最低非零值 1/n"] + \
            [f"{v}：{VERDICT_NOTE[v]}" for v in vs]
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


BRAIN_ROWS = [   # mouse brain, scRNA, author labels — explicit, no regex, no re-annotation
    ("astrocyte", "星形胶质"),
    ("oligodendrocyte", "少突胶质"),
    ("central nervous system macrophage", "小胶质/CNS巨噬"),
    ("Cajal-Retzius cell", "Cajal-Retzius 神经元"),
    ("VIP GABAergic cortical interneuron", "VIP 中间神经元"),
    ("hippocampal neuron", "海马神经元"),
    ("L2/3-6 intratelencephalic projecting glutamatergic neuron", "L2/3-6 谷氨酸能神经元"),
    ("sst GABAergic cortical interneuron", "SST 中间神经元"),
    ("pvalb GABAergic cortical interneuron", "PV 中间神经元"),
    ("endothelial cell", "内皮细胞"),
]


def draw_brain_dotplot():
    """Mouse brain (scRNA) dot plot: rows = cell type, columns = IL18 / CASP1 / CASP3.

    Dot area = pooled detection rate (same definition as everywhere else), colour = pseudobulk
    log2(CPM+1) with a separate scale per gene, because the three genes differ by orders of magnitude
    and one shared colour scale would render CASP1 uniformly blank. Every dot carries its own number,
    so the figure is readable without decoding dot area.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update({"font.family": "WenQuanYi Zen Hei", "axes.unicode_minus": False})
    cmap = LinearSegmentedColormap.from_list(
        "seq", ["#DCE9F8", "#9EC5F4", "#5E9AE8", "#2E6FC6", "#17457F"])

    c = pd.read_parquet(ROOT / "data" / "cells" / "M_brain_cell.parquet")
    stats = {}
    for label, _ in BRAIN_ROWS:
        sub = c[c.celltype_original == label]
        keep = [d for d, x in sub.groupby("donor_id") if len(x) >= MIN_CELLS_PER_DONOR]
        s = sub[sub.donor_id.isin(keep)]
        if not len(s):
            continue
        lib = s.lib.sum()
        stats[label] = {g: ((s[g] >= 1).mean(),
                            float(np.log2(s[g].sum() / lib * 1e6 + 1)) if lib > 0 else 0.0,
                            int((s[g] >= 1).sum()), len(s))
                        for g in ("IL18", "CASP1", "CASP3")}
    rows = [(lab, cn) for lab, cn in BRAIN_ROWS if lab in stats]
    rows = sorted(rows, key=lambda r: stats[r[0]]["IL18"][0], reverse=True)
    n = len(rows)

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 0.46 * n + 2.0), sharey=True,
                             gridspec_kw={"wspace": 0.12})
    for ax, gene in zip(axes, ("IL18", "CASP1", "CASP3")):
        vals = [stats[lab][gene] for lab, _ in rows]
        hi = max(v[1] for v in vals) or 1.0
        for i, (pct_, lcpm, pos, ntot) in enumerate(vals):
            ax.scatter([0.30], [i], s=40 + 1100 * pct_, c=[lcpm], cmap=cmap, vmin=0, vmax=hi,
                       edgecolors="white", linewidths=1.0, zorder=3)
            ax.text(0.50, i, f"{pct_*100:.2f}%" if pct_ < 0.1 else f"{pct_*100:.1f}%",
                    va="center", fontsize=9.4, color="#1A1A1A")
            ax.text(1.06, i, f"{pos}/{ntot}", va="center", ha="right", fontsize=8.6, color="#777777")
        ax.set_xlim(0, 1.09)
        ax.set_xticks([])
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(n - 0.5, -0.6)
        ax.set_title(f"{gene}（pro-IL-18 转录本）" if gene == "IL18" else gene,
                     fontsize=13, loc="left", pad=10)
        ax.grid(axis="y", color="#F0F0F0", lw=.7)
        ax.set_axisbelow(True)
        for s_ in ax.spines.values():
            s_.set_visible(False)
    axes[0].set_yticks(range(n), [cn for _, cn in rows], fontsize=10.2)
    axes[0].tick_params(axis="y", length=0)
    for p in (0.05, 0.25, 0.55):
        axes[2].scatter([], [], s=40 + 1100 * p, color="#8A8A8A", label=f"{int(p*100)}%")
    axes[2].legend(title="点面积 = 检出率", fontsize=9, title_fontsize=9, frameon=False,
                   loc="center left", bbox_to_anchor=(1.0, 0.5), labelspacing=1.5, borderpad=1.0)
    fig.text(0.012, 0.025,
             "颜色深浅 = log2(CPM+1)，三个基因各自独立色标（量级相差太大，不可共用一套色标，也不可跨基因比颜色）；"
             "每个点右侧为该基因的检出率与「阳性细胞/分析细胞」。",
             fontsize=9, color="#666666")
    fig.subplots_adjust(left=0.165, right=0.885, top=0.895, bottom=0.085)
    fig.savefig(FIG / "mouse" / "brain_cell_dots.png", dpi=170)
    print("wrote figures/mouse/brain_cell_dots.png", fig.get_size_inches() * 170)


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
    out["tier_cn"] = out.tier.map(TIER_CN)
    out["verdict"] = out.tier.map(VERDICT)
    out["verdict_rank"] = out.verdict.map({v: i for i, v in enumerate(VERDICT_ORDER)})
    out = out.sort_values(["verdict_rank", "pooled_pct_IL18"],
                          ascending=[True, False]).drop(columns="verdict_rank")
    cols = ["verdict", "tier", "tier_cn"] + [c for c in out.columns if c not in ("verdict", "tier", "tier_cn")]
    out = out[cols]
    write_tsv(out, ROOT / "results" / "answer_table.tsv")

    draw_answer_figure(out)
    draw_brain_dotplot()

    # ---- json for deck/report ----------------------------------------------------------------------
    js = [{"verdict": r.verdict, "verdict_note": VERDICT_NOTE[r.verdict],
           "tier": r.tier, "tier_cn": TIER_CN[r.tier], "sp": r.sp_cn, "organ": r.tissue_cn,
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
