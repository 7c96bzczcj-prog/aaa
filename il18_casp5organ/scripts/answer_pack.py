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

# explicit display names (no regex, no inference). Keyed by (species, tissue, author label).
SP_CN = {"human": "人", "mouse": "小鼠"}
TISSUE_CN = {"lung": "肺", "liver": "肝", "pancreas": "胰", "brain": "脑", "spinal_cord": "脊髓"}
CT_CN = {
    # mouse
    ("mouse", "brain", "oligodendrocyte"): "少突胶质",
    ("mouse", "brain", "astrocyte"): "星形胶质",
    ("mouse", "brain", "central nervous system macrophage"): "小胶质/CNS巨噬",
    ("mouse", "brain", "Cajal-Retzius cell"): "Cajal-Retzius 神经元",
    ("mouse", "lung", "alveolar macrophage"): "肺泡巨噬",
    ("mouse", "pancreas", "pancreatic D cell"): "胰岛 D(δ) 细胞",
    ("mouse", "pancreas", "pancreatic ductal cell"): "胰导管细胞",
    # human — brain (new in expansion batch v0.2)
    ("human", "brain", "pvalb GABAergic cortical interneuron"): "PV 中间神经元",
    ("human", "brain", "chandelier pvalb GABAergic cortical interneuron"): "枝形 PV 中间神经元",
    ("human", "brain", "L5 extratelencephalic projecting glutamatergic cortical neuron"): "L5 ET 投射神经元",
    ("human", "brain", "caudal ganglionic eminence derived cortical interneuron"): "CGE 源中间神经元",
    ("human", "brain", "vascular leptomeningeal cell"): "血管软脑膜细胞",
    ("human", "brain", "astrocyte of the cerebral cortex"): "星形胶质",
    ("human", "brain", "oligodendrocyte"): "少突胶质",
    ("human", "brain", "microglial cell"): "小胶质",
    # human — lung
    ("human", "lung", "mast cell"): "肥大细胞",
    ("human", "lung", "dendritic cell"): "树突状细胞",
    ("human", "lung", "ciliated cell"): "纤毛上皮",
    ("human", "lung", "multiciliated columnar cell of tracheobronchial tree"): "气道多纤毛上皮",
    ("human", "lung", "lung multiciliated epithelial cell"): "气道多纤毛上皮",
    ("human", "lung", "pulmonary alveolar type 1 cell"): "AT1 肺泡上皮",
    ("human", "lung", "alveolar adventitial fibroblast"): "肺泡外膜成纤维",
    # human — liver
    ("human", "liver", "basophil"): "嗜碱粒细胞",
    ("human", "liver", "cholangiocyte"): "胆管细胞",
    ("human", "liver", "CD16-negative, CD56-bright natural killer cell, human"): "CD56bright NK",
    ("human", "liver", "megakaryocyte-erythroid progenitor cell"): "巨核红系祖细胞",
    # human — pancreas
    ("human", "pancreas", "pancreatic acinar cell"): "胰腺腺泡细胞",
    ("human", "pancreas", "type B pancreatic cell"): "胰岛 β 细胞",
    ("human", "pancreas", "pancreatic D cell"): "胰岛 D(δ) 细胞",
    ("human", "pancreas", "mast cell"): "肥大细胞",
    # ---- remaining labels of the ten per-organ three-gene panels -----------------------------------
    # Only straight renderings of the ontology term; nothing is merged, split or re-interpreted here.
    ("human", "lung", "pulmonary alveolar type 2 cell"): "AT2 肺泡上皮",
    ("human", "lung", "epithelial cell of lower respiratory tract"): "下呼吸道上皮",
    ("human", "lung", "respiratory basal cell"): "气道基底细胞",
    ("human", "lung", "respiratory tract hillock cell"): "气道 hillock 细胞",
    ("human", "lung", "nasal mucosa goblet cell"): "鼻黏膜杯状细胞",
    ("human", "lung", "multiciliated epithelial cell"): "多纤毛上皮",
    ("human", "lung", "mesothelial cell"): "间皮细胞",
    ("human", "lung", "fibroblast"): "成纤维细胞",
    ("human", "lung", "alveolar macrophage"): "肺泡巨噬",
    ("human", "lung", "lung macrophage"): "肺巨噬细胞",
    ("human", "lung", "elicited macrophage"): "募集型巨噬细胞",
    ("human", "lung", "CD1c-positive myeloid dendritic cell"): "CD1c+ 髓系树突状",
    ("human", "lung", "plasmacytoid dendritic cell"): "浆细胞样树突状",
    ("human", "lung", "classical monocyte"): "经典单核细胞",
    ("human", "lung", "non-classical monocyte"): "非经典单核细胞",
    ("human", "lung", "unknown"): "未注释细胞群",
    ("human", "liver", "macrophage"): "巨噬细胞",
    ("human", "liver", "monocyte"): "单核细胞",
    ("human", "liver", "conventional dendritic cell"): "经典树突状细胞",
    ("human", "liver", "plasmacytoid dendritic cell"): "浆细胞样树突状",
    ("human", "liver", "neutrophil"): "中性粒细胞",
    ("human", "liver", "natural killer cell"): "NK 细胞",
    ("human", "liver", "T cell"): "T 细胞",
    ("human", "liver", "B cell"): "B 细胞",
    ("human", "liver", "plasma cell"): "浆细胞",
    ("human", "liver", "endothelial cell"): "内皮细胞",
    ("human", "pancreas", "pancreatic ductal cell"): "胰导管细胞",
    ("human", "pancreas", "pancreatic stellate cell"): "胰星状细胞",
    ("human", "pancreas", "pancreatic A cell"): "胰岛 α 细胞",
    ("human", "pancreas", "pancreatic PP cell"): "胰岛 PP(γ) 细胞",
    ("human", "pancreas", "endothelial cell"): "内皮细胞",
    ("human", "pancreas", "leukocyte"): "白细胞",
    ("human", "spinal_cord", "astrocyte"): "星形胶质",
    ("human", "spinal_cord", "oligodendrocyte"): "少突胶质",
    ("human", "spinal_cord", "oligodendrocyte precursor cell"): "少突胶质前体 OPC",
    ("human", "spinal_cord", "microglial cell"): "小胶质",
    ("human", "spinal_cord", "GABAergic neuron"): "GABA 能神经元",
    ("human", "spinal_cord", "glutamatergic neuron"): "谷氨酸能神经元",
    ("human", "spinal_cord", "cerebellar granule cell"): "颗粒细胞（原注释如此）",
    ("human", "spinal_cord", "capillary endothelial cell"): "毛细血管内皮",
    ("human", "spinal_cord", "endothelial cell of artery"): "动脉内皮",
    ("human", "spinal_cord", "mural cell"): "血管壁细胞",
    ("human", "spinal_cord", "leukocyte"): "白细胞",
    ("mouse", "lung", "fibroblast of lung"): "肺成纤维细胞",
    ("mouse", "lung", "adventitial cell"): "外膜细胞",
    ("mouse", "lung", "bronchial smooth muscle cell"): "支气管平滑肌",
    ("mouse", "lung", "vein endothelial cell"): "静脉内皮",
    ("mouse", "lung", "myeloid dendritic cell"): "髓系树突状细胞",
    ("mouse", "lung", "classical monocyte"): "经典单核细胞",
    ("mouse", "lung", "intermediate monocyte"): "中间型单核细胞",
    ("mouse", "lung", "non-classical monocyte"): "非经典单核细胞",
    ("mouse", "lung", "neutrophil"): "中性粒细胞",
    ("mouse", "lung", "basophil"): "嗜碱粒细胞",
    ("mouse", "lung", "natural killer cell"): "NK 细胞",
    ("mouse", "lung", "mature NK T cell"): "NKT 细胞",
    ("mouse", "lung", "T cell"): "T 细胞",
    ("mouse", "lung", "CD4-positive, alpha-beta T cell"): "CD4+ T 细胞",
    ("mouse", "lung", "CD8-positive, alpha-beta T cell"): "CD8+ T 细胞",
    ("mouse", "lung", "B cell"): "B 细胞",
    ("mouse", "liver", "Kupffer cell"): "Kupffer 细胞",
    ("mouse", "liver", "myeloid leukocyte"): "髓系白细胞",
    ("mouse", "liver", "hepatocyte"): "肝细胞",
    ("mouse", "liver", "endothelial cell of hepatic sinusoid"): "肝窦内皮",
    ("mouse", "liver", "natural killer cell"): "NK 细胞",
    ("mouse", "liver", "B cell"): "B 细胞",
    ("mouse", "pancreas", "type B pancreatic cell"): "胰岛 β 细胞",
    ("mouse", "pancreas", "pancreatic A cell"): "胰岛 α 细胞",
    ("mouse", "pancreas", "pancreatic PP cell"): "胰岛 PP(γ) 细胞",
    ("mouse", "pancreas", "pancreatic endocrine cell"): "胰岛内分泌细胞",
    ("mouse", "pancreas", "pancreatic stellate cell"): "胰星状细胞",
    ("mouse", "pancreas", "hematopoietic cell"): "造血来源细胞",
    ("mouse", "pancreas", "Schwann cell"): "施万细胞",
    ("mouse", "pancreas", "endothelial cell"): "内皮细胞",
    ("mouse", "pancreas", "unknown"): "未注释细胞群",
    ("mouse", "spinal_cord", "astrocytes"): "星形胶质",
    ("mouse", "spinal_cord", "Oligos"): "少突胶质",
    ("mouse", "spinal_cord", "neurons"): "神经元",
    ("mouse", "spinal_cord", "meninges/schwann"): "脑膜/施万细胞",
    ("mouse", "spinal_cord", "endo/vascular"): "内皮/血管",
}
# short dataset tags — the same cell type now appears in several strata, so the row label must say which
DS_CN = {
    "H_lung_cell": "HLCA", "H_lung_cell_cellref": "CellRef", "H_lung_cell_distal": "Distal",
    "H_lung_nucleus": "肺单核", "H_liver_cell": "Guilliams", "H_liver_cell_full": "肝 full",
    "H_liver_cell_global": "肝 global", "H_liver_nucleus": "肝单核",
    "H_pancreas_cell": "TabulaSapiens", "H_pancreas_nucleus": "胰单核",
    "H_panc_islet_mo": "胰岛 multiomics", "H_panc_islet_48": "胰岛 48-donor", "H_panc_islet_xs": "胰岛 x-species",
    "H_brain_nucleus": "Siletti MTG", "H_brain_cell_hvs": "HVS 皮层", "H_brain_nucleus_seaad_mtg": "SEA-AD MTG",
    "H_brain_nucleus_seaad_dlpfc": "SEA-AD DLPFC", "H_brain_nucleus_psychad": "PsychAD",
    "H_spinal_nucleus": "Siletti 脊髓", "H_spinal_nucleus_wm": "脊髓白质",
    "M_lung_cell": "TMS 10x", "M_lung_cell_ss2": "TMS SS2", "M_liver_cell": "TMS 10x",
    "M_liver_cell_ss2": "TMS SS2", "M_pancreas_cell": "TMS 全胰", "M_pancreas_islet_cell": "小鼠胰岛图谱",
    "M_brain_cell": "Yao 皮层+海马", "M_brain_nucleus": "皮层单核", "M_spinal_nucleus": "GSE103892",
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

# What the reader sees is TWO axes, not one adjective.
#
# The earlier 确定 / 可能 / 相反 scheme mixed them: 确定-vs-可能 was about how many independent sources
# agree, while 相反 was about which way the evidence points. Collapsing both into one three-way label
# meant "strongly shown not to fit" and "weakly shown not to fit" became the same word, and "one source
# supports it" shared a label with "nobody has ever measured this cell type". It also read as a claim
# about certainty that transcript-level screening data cannot carry.
#
# So每行 now gets:
#   direction  符合 / 不符合      — does this cell type fit the premise (IL18 high AND CASP1 low)?
#   evidence   多源一致 / 单一来源 / 无参考可查   — how much external evidence is behind that direction
#   label      a sentence stating what was actually observed, with the source count in it
DIRECTION = {
    "A_supported": "符合", "B_potential": "符合", "B_potential_screen_missed": "符合",
    "C_casp1_present": "不符合", "D_refuted": "不符合",
    "E_no_matched_reference": None,      # resolved from the row's own data below
}
# 证据 is not just "how many" — "four sources, one of them dissenting" is a different state from
# "four sources, all agreeing", and the earlier wording collapsed them into 单一来源, which was wrong.
# colour carries direction + strength: green = fits with several sources, blue = fits but thin evidence,
# red = does not fit. Same three colours as before, but the label now says which is which.
LABEL_COL = {"符合·多源一致": "#1B7A43", "符合·多源有分歧": "#2E8B57", "符合·单一来源": "#1F5FA8",
             "符合·无参考可查": "#4C87C4", "不符合·多源一致": "#9A2B2B", "不符合·无参考可查": "#B5564E"}
OWN_DATA_RATIO_WEAK = 3.0   # for rows with no matched reference, judge direction from their own data


def classify(tier, n_sup, n_dis, pooled_il18, pooled_casp1, screen_missed=False):
    """Return (direction, evidence, label). The label states the observation, not a verdict adjective."""
    if tier == "E_no_matched_reference":
        own = (pooled_il18 / pooled_casp1) if pooled_casp1 > 0 else float("inf")
        if own < OWN_DATA_RATIO_WEAK:
            return ("不符合", "无参考可查",
                    f"无参考可查；本数据 CASP1 与 IL18 同量级（{pooled_casp1*100:.2f}% vs {pooled_il18*100:.1f}%）")
        return ("符合", "无参考可查",
                f"无参考可查；仅本数据方向相符（CASP1 {pooled_casp1*100:.2f}%）")
    if tier == "A_supported":
        return "符合", "多源一致", f"CASP1 多源确认近零（{n_sup} 个独立来源，无分歧）"
    if tier in ("B_potential", "B_potential_screen_missed"):
        tail = "；被分层内相对规则漏掉" if screen_missed else ""
        if n_dis:
            return ("符合", "多源有分歧",
                    f"{n_sup} 个来源支持 CASP1 近零，另有 {n_dis} 个来源显示 CASP1 可测{tail}")
        return "符合", "单一来源", f"CASP1 近零，仅 {n_sup} 个来源支持{tail}"
    if tier == "C_casp1_present":
        return "不符合", "多源一致", "参考层 CASP1 可测到，仅比 IL18 低 3–20 倍"
    return "不符合", "多源一致", "参考层 CASP1 ≥ IL18（方向相反）"


def name_of(sp, ti, ct, dataset_id):
    """Display name for one screen row. KeyError here is deliberate: a new author label must be given a
    name by hand rather than falling back to the raw ontology string."""
    return SP_CN[sp], TISSUE_CN[ti], CT_CN[(sp, ti, str(ct))], DS_CN.get(dataset_id, dataset_id)


DIR_ORDER = ["符合", "不符合"]
EVI_ORDER = ["多源一致", "多源有分歧", "单一来源", "无参考可查"]


ZERO_SLOT = 0.0018   # where a true zero (0 positive cells) is drawn on the log axis
XMIN, XMAX = 0.0011, 150.0


def draw_answer_figure(out, species=("human", "mouse"), out_name="answer_v2.png", row_h=0.34,
                       collapse=False):
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

    if collapse:
        # One row per cell type instead of per stratum: keep the stratum with the most analysed donors as
        # the representative and let the "xN 数据集" tag carry the replication. Nothing is pooled across
        # strata — the numbers shown are still one stratum's own numbers.
        out = (out.sort_values("n_donors_analysed", ascending=False)
                  .drop_duplicates(["species", "tissue", "celltype"])
                  .sort_values(["verdict_rank_fig", "pooled_pct_IL18"], ascending=[True, False]))
    groups = [(SP_CN[sp], out[out.species == sp]) for sp in species]
    groups = [(k, g) for k, g in groups if len(g)]
    counts = [len(g) for _, g in groups]
    n_all = sum(counts)
    fig, axes = plt.subplots(len(groups), 1, figsize=(17.0, row_h * n_all + 2.4), squeeze=False,
                             gridspec_kw={"height_ratios": counts, "hspace": 0.26})
    axes = axes[:, 0]

    for ax, (sp, g) in zip(axes, groups):
        d = g.iloc[::-1].reset_index(drop=True)      # bottom-up so 确定 lands on top
        n = len(d)
        for i, r in enumerate(d.itertuples()):
            col = LABEL_COL[r.verdict]
            il18 = r.pooled_pct_IL18 * 100
            casp1 = r.pooled_pct_CASP1 * 100
            x_c = casp1 if r.n_cells_CASP1_pos else ZERO_SLOT
            ax.plot([x_c, il18], [i, i], color=col, lw=3.0, alpha=.22, solid_capstyle="round", zorder=1)
            ax.vlines(100.0 / r.n_cells_analysed, i - .26, i + .26, color="#B5B5B5", lw=1.0, zorder=2)
            ax.scatter([x_c], [i], s=66, facecolors="white", edgecolors=col, linewidths=2.1, zorder=4)
            ax.scatter([il18], [i], s=78, color=col, zorder=4)
            c_txt = "0" if not r.n_cells_CASP1_pos else f"{casp1:.2f}%"
            if x_c <= il18:
                ax.text(il18 * 1.22, i, f"{il18:.1f}%", va="center", fontsize=9.8, color=col, weight="bold")
                ax.text(x_c / 1.28, i, c_txt, va="center", ha="right", fontsize=9.3, color=col)
            else:                                   # CASP1 above IL18 — swap the sides so labels never collide
                ax.text(x_c * 1.22, i, c_txt, va="center", fontsize=9.3, color=col)
                ax.text(il18 / 1.28, i, f"{il18:.1f}%", va="center", ha="right", fontsize=9.8,
                        color=col, weight="bold")
            rep = f"  ×{r.n_strata} 数据集" if r.n_strata > 1 else ""
            ax.text(1.006, i, f"{r.n_cells_CASP1_pos} · {r.n_cells_IL18_pos} / {r.n_cells_analysed}{rep}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left",
                    fontsize=8.9, color="#555555")
        ax.set_xscale("log")
        ax.set_xlim(XMIN, XMAX)
        ax.set_yticks(np.arange(n), [f"{r.tissue_cn} {r.celltype_cn}（{SUSP[r.suspension_type]}·{r.dataset_cn}）"
                                     for r in d.itertuples()], fontsize=9.4)
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.72, n - 0.28)
        ax.grid(axis="x", color="#ECECEC", lw=.7)
        ax.set_axisbelow(True)
        ax.text(-0.235, 1.0, sp, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=15, weight="bold", color="#1A1A1A")
        ax.text(1.006, 1.0, "CASP1 · IL18 阳性细胞 / 分析细胞\n×n = 在几个独立数据集里重复", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=8.9, color="#555555")
        ax.set_xticks([0.01, 0.1, 1, 10, 100], ["0.01%", "0.1%", "1%", "10%", "100%"], fontsize=9.6)

    for ax in axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False, length=0)
    axes[-1].set_xlabel("检出率（对数轴）= 该基因原始计数 ≥1 的细胞数 ÷ 分析细胞数（全部供体合并）", fontsize=10.6)

    marks = [Line2D([], [], marker="o", ls="", ms=8.6, color="#444"),
             Line2D([], [], marker="o", ls="", ms=8.6, mfc="white", mec="#444", mew=2.1),
             Line2D([], [], marker="|", ls="", ms=9, mew=1.4, color="#B5B5B5")]
    seen = list(dict.fromkeys(out.sort_values("verdict_rank_fig").verdict))
    handles = marks + [Line2D([], [], marker="s", ls="", ms=8, color=LABEL_COL[v]) for v in seen]
    names = ["IL18（pro-IL-18 转录本）", "CASP1", "该行可检出的最低非零值 1/n"] + \
            [v.replace("·", " / ") for v in seen]
    axes[0].legend(handles, names, loc="lower left", bbox_to_anchor=(-0.235, 1.07), ncol=3,
                   frameon=False, fontsize=9.5, handletextpad=.5, columnspacing=1.8)
    fig.text(0.012, 0.012,
             "分母 = 该细胞类型在该分层（数据集 × 器官 × 单细胞/单核）内、每位供体 ≥20 个细胞后进入分析的细胞总数；"
             "分子 = 其中该基因原始计数 ≥1 的细胞数。不是表达量占比、不是相对最高值的百分比、不是读段占比。"
             "人与小鼠、单细胞与单核之间不作定量比较。",
             fontsize=8.9, color="#666666")
    fig.subplots_adjust(left=0.205, right=0.825, top=0.885, bottom=0.085)
    fig.savefig(FIG / out_name, dpi=170)
    print(f"wrote figures/{out_name}", fig.get_size_inches() * 170)
    plt.close(fig)


# Human cortex scRNA (Allen HVS, 78 donors) — ALL 18 author labels in that stratum, candidates and
# non-candidates alike, so the panel is read as a whole stratum rather than a shortlist.
HUMAN_BRAIN_ROWS = [
    ("microglial cell", "小胶质"),
    ("L5 extratelencephalic projecting glutamatergic cortical neuron", "L5 ET 投射神经元"),
    ("pvalb GABAergic cortical interneuron", "PV 中间神经元"),
    ("chandelier pvalb GABAergic cortical interneuron", "枝形 PV 中间神经元"),
    ("corticothalamic-projecting glutamatergic cortical neuron", "L6 CT 投射神经元"),
    ("near-projecting glutamatergic cortical neuron", "L5/6 NP 近程投射神经元"),
    ("L6b glutamatergic cortical neuron", "L6b 神经元"),
    ("L2/3-6 intratelencephalic projecting glutamatergic neuron", "L2/3-6 IT 神经元"),
    ("sst GABAergic cortical interneuron", "SST 中间神经元"),
    ("caudal ganglionic eminence derived cortical interneuron", "CGE 源中间神经元"),
    ("lamp5 GABAergic cortical interneuron", "LAMP5 中间神经元"),
    ("sncg GABAergic cortical interneuron", "SNCG 中间神经元"),
    ("VIP GABAergic cortical interneuron", "VIP 中间神经元"),
    ("vascular leptomeningeal cell", "血管软脑膜细胞"),
    ("astrocyte of the cerebral cortex", "星形胶质"),
    ("oligodendrocyte precursor cell", "少突胶质前体 OPC"),
    ("cerebral cortex endothelial cell", "内皮细胞"),
    ("oligodendrocyte", "少突胶质"),
]

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


# one three-gene panel per species × organ, in the order the deck shows them. Entry chosen per organ
# as the stratum with the most analysable cell types at the author level (see run_log.md).
PANEL_STRATA = [
    ("H_lung_cell", "human", "lung", "human/lung_dots.png"),
    ("H_liver_cell", "human", "liver", "human/liver_dots.png"),
    ("H_panc_islet_48", "human", "pancreas", "human/pancreas_dots.png"),
    ("H_spinal_nucleus_wm", "human", "spinal_cord", "human/spinal_dots.png"),
    ("M_lung_cell", "mouse", "lung", "mouse/lung_dots.png"),
    ("M_liver_cell", "mouse", "liver", "mouse/liver_dots.png"),
    ("M_pancreas_islet_cell", "mouse", "pancreas", "mouse/pancreas_dots.png"),
    ("M_spinal_nucleus", "mouse", "spinal_cord", "mouse/spinal_dots.png"),
]
# panels whose underlying stratum carries a caveat that must travel with the figure
PANEL_NOTE = {
    "mouse/spinal_dots.png": "GSE103892 仅带基因符号（gene_id_unavailable），不进候选判定。",
    "human/spinal_dots.png": "颈段白质，不含灰质与运动神经元。",
}


def wrap_cjk(s, width_in, fontsize=9.0, margin=0.06):
    """Wrap a Chinese caption to a figure width. textwrap can't: CJK has no spaces to break on.

    Widths are counted in em (1 for CJK/full-width, ~0.55 for ASCII) and lines are broken after a
    punctuation mark where one is within reach, so a break never lands inside a number or a gene name.
    """
    budget = (width_in * (1 - 2 * margin)) * 72 / fontsize
    brk = "，。；：）」、,;)"
    lines, cur, w = [], "", 0.0
    for ch in s:
        cur += ch
        w += 1.0 if ord(ch) > 0x2E80 else 0.55
        if w >= budget:
            cut = max((cur.rfind(b) for b in brk), default=-1)
            if cut >= len(cur) * 0.55:                     # only honour a break point near the end
                lines.append(cur[:cut + 1])
                cur = cur[cut + 1:]
                w = sum(1.0 if ord(c) > 0x2E80 else 0.55 for c in cur)
            else:
                lines.append(cur)
                cur, w = "", 0.0
    if cur:
        lines.append(cur)
    return lines


def rows_for_stratum(entry_id, species, tissue, max_rows=20):
    """Every author label in one stratum that cleared the donor floor, ranked by IL18 detection.

    Used so an organ panel shows the whole stratum rather than a hand-picked shortlist. Chinese names
    are used where one has been assigned by hand; otherwise the ontology label is kept verbatim, which
    is precise and traceable — inventing a name for a label nobody has reviewed would be worse.
    """
    s = pd.read_csv(ROOT / "results" / "il18_casp_summary.tsv", sep="\t", comment="#")
    d = s[(s.dataset_id == entry_id) & (s.level == "author") & (s.n_donors >= 3)]
    d = d.sort_values("pct_IL18", ascending=False)
    labs = list(d.celltype)[:max_rows]
    return [(l, CT_CN.get((species, tissue, l), l)) for l in labs], len(d)


PANEL_IN = (11.93, 5.30)   # the image box left on a 13.33 × 7.5 slide once title and caption are placed
PANEL_FS = 8.0             # every point size below is what the reader actually sees on that slide


def draw_brain_dotplot(entry_id="M_brain_cell", spec_rows=None, out_name="mouse/brain_cell_dots.png",
                       verdicts=None, species="mouse", tissue="brain",
                       n_total=None, extra_note=""):
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

    spec_rows = spec_rows if spec_rows is not None else BRAIN_ROWS
    verdicts = verdicts or {}
    c = pd.read_parquet(ROOT / "data" / "cells" / f"{entry_id}.parquet")
    stats = {}
    for label, _ in spec_rows:
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
    rows = [(lab, cn) for lab, cn in spec_rows if lab in stats]
    rows = sorted(rows, key=lambda r: stats[r[0]]["IL18"][0], reverse=True)
    n = len(rows)

    # The canvas is exactly the image box of the slide, so a point size here is the point size the
    # reader sees. Earlier these were drawn large and let the deck shrink them, which put 20-row panels
    # on screen at under 5 pt. Row pitch and dot size follow from the row count instead of the reverse.
    fw, fh = PANEL_IN
    scope = "" if n_total is None or n_total <= n else f"本分层 {n_total} 类，按 IL18 取前 {n}。"
    foot = wrap_cjk("点面积 = 检出率，右侧为阳性细胞/分析细胞；颜色 = log2(CPM+1)，三基因各自独立色标。"
                    "行名颜色 = 判定方向（绿/蓝 = 符合，红 = 不符合，黑 = 未进候选）。" + scope + extra_note,
                    fw, PANEL_FS - 1.2, margin=0.012)
    top_in, line_in = 0.42, 0.145
    foot_in = line_in * len(foot) + 0.10
    row_pt = (fh - top_in - foot_in) / n * 72                 # vertical pitch of one row, in points
    dmax = min(0.68 * row_pt, 26.0)                           # keep the largest dot inside its row
    smax, smin = (dmax / 2) ** 2 * np.pi, 3.0

    fig, axes = plt.subplots(1, 3, figsize=(fw, fh), sharey=True, gridspec_kw={"wspace": 0.12})
    for ax, gene in zip(axes, ("IL18", "CASP1", "CASP3")):
        vals = [stats[lab][gene] for lab, _ in rows]
        hi = max(v[1] for v in vals) or 1.0
        for i, (pct_, lcpm, pos, ntot) in enumerate(vals):
            ax.scatter([0.30], [i], s=smin + (smax - smin) * pct_, c=[lcpm], cmap=cmap, vmin=0,
                       vmax=hi, edgecolors="white", linewidths=0.7, zorder=3)
            ax.text(0.50, i, f"{pct_*100:.2f}%" if pct_ < 0.1 else f"{pct_*100:.1f}%",
                    va="center", fontsize=PANEL_FS, color="#1A1A1A")
            ax.text(1.06, i, f"{pos}/{ntot}", va="center", ha="right", fontsize=PANEL_FS - 1.1,
                    color="#777777")
        ax.set_xlim(0, 1.09)
        ax.set_xticks([])
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(n - 0.5, -0.6)
        ax.set_title(f"{gene}（pro-IL-18 转录本）" if gene == "IL18" else gene,
                     fontsize=PANEL_FS + 1.8, loc="left", pad=7)
        ax.grid(axis="y", color="#F0F0F0", lw=.7)
        ax.set_axisbelow(True)
        for s_ in ax.spines.values():
            s_.set_visible(False)
    axes[0].set_yticks(range(n), [cn for _, cn in rows], fontsize=PANEL_FS)
    axes[0].tick_params(axis="y", length=0)
    # colour the row label by this cell type's verdict, so all three classes are visible at once
    for tick, (lab, _) in zip(axes[0].get_yticklabels(), rows):
        tick.set_color(LABEL_COL.get(verdicts.get((species, tissue, lab)), "#1A1A1A"))
    for p in (0.05, 0.25, 0.55):
        axes[2].scatter([], [], s=smin + (smax - smin) * p, color="#8A8A8A", label=f"{int(p*100)}%")
    axes[2].legend(title="点面积 = 检出率", fontsize=PANEL_FS - 0.6, title_fontsize=PANEL_FS - 0.6,
                   frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5),
                   labelspacing=1.4, borderpad=0.8)
    for i, line in enumerate(foot):
        fig.text(0.012, (line_in * (len(foot) - i) - 0.045) / fh, line,
                 fontsize=PANEL_FS - 1.2, color="#666666")
    fig.subplots_adjust(left=1.62 / fw, right=1 - 1.30 / fw, top=1 - top_in / fh, bottom=foot_in / fh)
    (FIG / out_name).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / out_name, dpi=240)
    print(f"wrote figures/{out_name}", (fig.get_size_inches() * 240).astype(int))
    plt.close(fig)


def main():
    ev = pd.read_csv(ROOT / "results" / "il18_casp1_evidence.tsv", sep="\t", comment="#")
    ev = ev[ev.duplicate_of.isna() | (ev.duplicate_of == "")].copy()
    cmap = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#")

    rows = []
    for _, r in ev.iterrows():
        sp_cn, ti_cn, ct_cn, ds_cn = name_of(r.species, r.tissue, r.celltype, r.dataset_id)
        rows.append({"tier": r.evidence_tier, "species": r.species, "tissue": r.tissue,
                     "suspension_type": r.suspension_type, "celltype": r.celltype,
                     "sp_cn": sp_cn, "tissue_cn": ti_cn, "celltype_cn": ct_cn, "dataset_cn": ds_cn,
                     "dataset_id": r.dataset_id, "level": r.level,
                     "n_cells_analysed": r.n_cells_analysed, "n_donors_analysed": r.n_donors_analysed,
                     "median_cells_per_donor": r.median_cells_per_donor,
                     "n_cells_IL18_pos": r.n_cells_IL18_pos, "pooled_pct_IL18": r.pooled_pct_IL18,
                     "donors_with_IL18_pos": r.donors_with_IL18_pos, "donor_median_pct_IL18": r.sc_pct_IL18,
                     "n_cells_CASP1_pos": r.n_cells_CASP1_pos, "pooled_pct_CASP1": r.pooled_pct_CASP1,
                     "donors_with_CASP1_pos": r.donors_with_CASP1_pos, "donor_median_pct_CASP1": r.sc_pct_CASP1,
                     "ref_ratio_IL18_over_CASP1": r.best_ratio_IL18_over_CASP1,
                     "n_ref_sources": r.n_ref_sources,
                     "n_sources_supporting": r.n_sources_supporting,
                     "n_sources_dissenting": r.n_sources_dissenting, "tier_reason": r.tier_reason,
                     "screen_status": "进入候选"})

    # AT1: passed IL18 but was excluded by the within-stratum relative CASP1 rule (§8), reference-supported
    miss = pd.read_csv(ROOT / "results" / "screen_missed_but_reference_supported.tsv", sep="\t", comment="#")
    for _, r in miss.iterrows():
        ps = pooled_stats(r.dataset_id, "author_label", r.celltype, cmap)
        sp_cn, ti_cn, ct_cn, ds_cn = name_of(r.species, r.tissue, r.celltype, r.dataset_id)
        rows.append({"tier": "B_potential_screen_missed", "species": r.species, "tissue": r.tissue,
                     "suspension_type": r.suspension_type, "celltype": r.celltype,
                     "sp_cn": sp_cn, "tissue_cn": ti_cn, "celltype_cn": ct_cn, "dataset_cn": ds_cn,
                     "dataset_id": r.dataset_id, "level": "author_label",
                     **{k: ps[k] for k in ("n_cells_analysed", "n_donors_analysed", "median_cells_per_donor",
                                           "n_cells_IL18_pos", "pooled_pct_IL18", "donors_with_IL18_pos",
                                           "n_cells_CASP1_pos", "pooled_pct_CASP1", "donors_with_CASP1_pos")},
                     "donor_median_pct_IL18": r.pct_IL18, "donor_median_pct_CASP1": r.pct_CASP1,
                     "ref_ratio_IL18_over_CASP1": r.best_ref.split("ratio=")[-1],
                     "n_ref_sources": r.n_ref_sources_supporting,
                     "n_sources_supporting": r.n_ref_sources_supporting, "n_sources_dissenting": 0,
                     "tier_reason": f"IL18 分层内排名 {int(r.rank_IL18)}/{int(r.n_ct_stratum)}，"
                                    f"参考层支持 CASP1 低；{r.excluded_because}",
                     "screen_status": "未进候选（CASP1 未进最低三分位）"})

    out = pd.DataFrame(rows)
    # Replication is its own evidence axis, separate from the reference layer: the same cell type found
    # in several INDEPENDENT strata (different datasets, different donors) is harder to explain away as
    # one atlas's annotation or batch quirk. Counted over distinct rows only, never pooled across strata.
    # `out` already excludes the lineage re-aggregation rows (ev is filtered on duplicate_of above).
    rep = out.groupby(["species", "tissue", "celltype"]).dataset_id.nunique()
    out["n_strata"] = [int(rep.get((r.species, r.tissue, r.celltype), 1)) for r in out.itertuples()]
    out["tier_cn"] = out.tier.map(TIER_CN)
    cls = [classify(r.tier, int(r.n_sources_supporting or 0), int(r.n_sources_dissenting or 0),
                    r.pooled_pct_IL18, r.pooled_pct_CASP1,
                    r.tier == "B_potential_screen_missed") for r in out.itertuples()]
    out["direction"] = [c[0] for c in cls]
    out["evidence"] = [c[1] for c in cls]
    out["label"] = [c[2] for c in cls]
    out["verdict"] = out.direction + "·" + out.evidence     # the colour key, and the sort key
    out.loc[out.tier == "E_no_matched_reference", "tier_reason"] = out.loc[out.tier == "E_no_matched_reference", "label"]
    rank = {f"{d}·{e}": i * 10 + j for i, d in enumerate(DIR_ORDER) for j, e in enumerate(EVI_ORDER)}
    out["verdict_rank"] = out.verdict.map(rank)
    out = out.sort_values(["verdict_rank", "pooled_pct_IL18"], ascending=[True, False])
    out = out.rename(columns={"verdict_rank": "verdict_rank_fig"})
    tsv = out.drop(columns=["verdict_rank_fig", "verdict"])
    front = ["direction", "evidence", "label", "n_sources_supporting", "n_sources_dissenting",
             "n_strata", "tier", "tier_cn"]
    cols = front + [c for c in tsv.columns if c not in front]
    write_tsv(tsv[cols], ROOT / "results" / "answer_table.tsv")

    draw_answer_figure(out)                                              # combined, for the report
    draw_answer_figure(out, ("human",), "answer_human.png", row_h=0.42, collapse=True)   # one slide each
    draw_answer_figure(out, ("mouse",), "answer_mouse.png", row_h=0.42, collapse=True)
    # keyed by (species, tissue, author label): both species share label strings such as
    # "oligodendrocyte", and so do organs ("macrophage", "endothelial cell", "mast cell") — a map
    # missing either key paints one row's verdict onto an unrelated row
    # when one cell type appears in several strata of the same organ, keep the stratum with the most
    # donors — the same tie-break the collapsed overview figure uses, so the two agree
    vmap = {(r.species, r.tissue, r.celltype): r.verdict
            for r in out.sort_values("n_donors_analysed").itertuples()}
    draw_brain_dotplot("H_brain_cell_hvs", HUMAN_BRAIN_ROWS, "human/brain_cell_dots.png", vmap,
                       "human", "brain")
    draw_brain_dotplot("M_brain_cell", BRAIN_ROWS, "mouse/brain_cell_dots.png", vmap, "mouse", "brain")
    # the other four organs, same treatment — the brain-only pair was an accident of how the page was
    # first built, not a decision; the screen covers five organs in both species
    for eid, sp, tis, fname in PANEL_STRATA:
        rws, total = rows_for_stratum(eid, sp, tis)
        if not rws:
            print(f"skip {fname}: no cell type clears the donor floor")
            continue
        draw_brain_dotplot(eid, rws, fname, vmap, sp, tis, n_total=total,
                           extra_note=PANEL_NOTE.get(fname, ""))
        if len(rws) < total:
            print(f"   note: {fname} shows the top {len(rws)} of {total} cell types by IL18 detection")

    # ---- json for deck/report ----------------------------------------------------------------------
    js = [{"direction": r.direction, "evidence": r.evidence, "label": r.label,
           "n_sup": int(r.n_sources_supporting or 0), "n_dis": int(r.n_sources_dissenting or 0), "verdict": r.verdict,
           "tier": r.tier, "tier_cn": TIER_CN[r.tier], "sp": r.sp_cn, "organ": r.tissue_cn,
           "susp": SUSP[r.suspension_type], "name": r.celltype_cn,
           "il18": round(r.pooled_pct_IL18 * 100, 2), "il18_pos": int(r.n_cells_IL18_pos),
           "casp1": round(r.pooled_pct_CASP1 * 100, 3), "casp1_pos": int(r.n_cells_CASP1_pos),
           "n_cells": int(r.n_cells_analysed), "n_donors": int(r.n_donors_analysed),
           "n_strata": int(r.n_strata), "dataset_cn": r.dataset_cn,
           "casp1_donors": int(r.donors_with_CASP1_pos),
           "med_casp1": round(r.donor_median_pct_CASP1 * 100, 2),
           "ratio": ("" if pd.isna(r.ref_ratio_IL18_over_CASP1) else str(r.ref_ratio_IL18_over_CASP1)),
           "src": int(r.n_ref_sources),
           "why": r.tier_reason, "status": r.screen_status} for r in out.itertuples()]
    (ROOT / "results" / "answer_table.json").write_text(json.dumps(js, ensure_ascii=False, indent=1))
    print(f"wrote results/answer_table.tsv + .json ({len(js)} rows)")


if __name__ == "__main__":
    main()
