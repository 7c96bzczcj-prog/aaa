"""
The two figures.

Fig 1  cervical NK effector score by tissue type, violin + one point per donor.
Fig 2  cervical vs lung tumour NK effector score, split by dataset, with a
       depth-matched panel beside the CP10K panel.

Fig 2 is deliberately drawn per dataset rather than pooling each cancer type
into one violin. The datasets differ in dissociation, enrichment and sequencing
depth, so between-dataset spread within one cancer type is the yardstick any
between-cancer difference has to clear. Pooling would hide exactly that.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import effector_common as ec

OUT = "/home/user/aaa/results/xcancer"
FIGS = "/home/user/aaa/results/xcancer/figures"
os.makedirs(FIGS, exist_ok=True)

# ---- palette (validated: see dataviz validator run) -------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
# ordinal disease-stage ramp, one hue, monotone lightness (PASS, --ordinal)
STAGE_COLORS = {
    "healthy": "#86b6ef",
    "adjacent": "#3987e5",
    "HSIL": "#1c5cab",
    "cancer": "#0d366b",
}
STAGE_ORDER = ["healthy", "adjacent", "HSIL", "cancer"]
STAGE_LABEL = {
    "healthy": "healthy cervix",
    "adjacent": "tumour-adjacent",
    "HSIL": "HSIL / CIN",
    "cancer": "cervical cancer",
}
# cancer-type categorical slots 1 and 2 (PASS, --pairs all)
CANCER_COLORS = {"cervical": "#2a78d6", "lung": "#eb6834"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
})


# below this many cells a kernel density is a drawing, not a distribution
MIN_CELLS_FOR_VIOLIN = 50


def half_violin(ax, x, vals, color, width=0.34, rng=None):
    """
    Thin violin; no extrema whiskers, no box - the donor points carry detail.

    Under MIN_CELLS_FOR_VIOLIN the individual cells are drawn instead. A KDE
    over a few dozen points invents shape it has no support for, and the
    tumour-adjacent group here has 39 cells across 3 donors.
    """
    if len(vals) == 0:
        return
    if len(vals) < MIN_CELLS_FOR_VIOLIN:
        rng = rng or np.random.default_rng(0)
        xs = x + rng.uniform(-0.16, 0.16, size=len(vals))
        ax.scatter(xs, vals, s=7, facecolor=color, alpha=0.55,
                   edgecolor="none", zorder=2)
        return
    parts = ax.violinplot([vals], positions=[x], widths=width * 2,
                          showextrema=False, showmedians=False)
    for b in parts["bodies"]:
        b.set_facecolor(color)
        b.set_alpha(0.28)
        b.set_edgecolor(color)
        b.set_linewidth(1.0)


def donor_points(ax, x, donor_med, color, rng, jitter=0.10):
    xs = x + rng.uniform(-jitter, jitter, size=len(donor_med))
    ax.scatter(xs, donor_med, s=34, facecolor=color, edgecolor=SURFACE,
               linewidth=1.2, zorder=4, clip_on=False)


def fig1(cells, donors, score="score_cp10k"):
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for i, g in enumerate(STAGE_ORDER):
        sub = cells[cells.tissue_group == g]
        dsub = donors[donors.tissue_group == g]
        if not len(sub):
            continue
        vals = sub[score].dropna().values
        half_violin(ax, i, vals, STAGE_COLORS[g], rng=rng)
        donor_points(ax, i, dsub["median"].values, STAGE_COLORS[g], rng)
        med = float(np.median(vals))
        ax.plot([i - 0.36, i + 0.36], [med, med], color=STAGE_COLORS[g],
                lw=2.0, zorder=5, solid_capstyle="round")
        note = f"{len(sub):,} cells\n{len(dsub)} donors"
        if len(vals) < MIN_CELLS_FOR_VIOLIN:
            note += "\n(cells shown\nindividually)"
        ax.annotate(note, xy=(i, -0.115), xycoords=("data", "axes fraction"),
                    ha="center", va="top", fontsize=8.5,
                    color=MUTED, annotation_clip=False)

    ax.set_xticks(range(len(STAGE_ORDER)))
    ax.set_xticklabels([STAGE_LABEL[g] for g in STAGE_ORDER], fontsize=10, color=INK2)
    ax.set_ylabel("NK effector score\nmean log1p(CP10K) over 6 genes", fontsize=10)
    ax.set_title("Cervical NK effector-gene expression by tissue type",
                 fontsize=13, color=INK, pad=14, loc="left")
    ax.annotate("GZMB · PRF1 · GZMA · IFNG · NKG7 · KLRD1   ·   violin = cells, "
                "dot = one donor, bar = group median",
                xy=(0, 1.015), xycoords="axes fraction", fontsize=9, color=MUTED)
    ax.set_xlim(-0.6, len(STAGE_ORDER) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", visible=False)
    fig.subplots_adjust(bottom=0.30, top=0.86, left=0.13, right=0.97)
    p = f"{FIGS}/fig1_cervical_effector_by_tissue.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def fig2(cells, donors):
    rng = np.random.default_rng(1)
    tum_c = cells[(cells.tissue_group == "cancer")]
    tum_d = donors[(donors.tissue_group == "cancer")]
    order, colors = [], []
    for canc in ("cervical", "lung"):
        for ds in sorted(tum_c[tum_c.cancer == canc].dataset.unique()):
            order.append((canc, ds))
            colors.append(CANCER_COLORS[canc])

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), sharey=False)
    panels = [("score_cp10k", "median", "CP10K normalisation",
               "mean log1p(CP10K) over 6 genes"),
              ("score_ds", "median_ds", f"depth-matched to {ec.DOWNSAMPLE_TO} UMI",
               f"same score after downsampling every cell to {ec.DOWNSAMPLE_TO} UMI")]

    for ax, (scol, dcol, sub_t, ylab) in zip(axes, panels):
        for i, ((canc, ds), col) in enumerate(zip(order, colors)):
            sub = tum_c[(tum_c.cancer == canc) & (tum_c.dataset == ds)]
            dsub = tum_d[(tum_d.cancer == canc) & (tum_d.dataset == ds)]
            v = sub[scol].dropna().values
            if not len(v):
                continue
            half_violin(ax, i, v, col, rng=rng)
            if dcol in dsub:
                donor_points(ax, i, dsub[dcol].dropna().values, col, rng)
            med = float(np.median(v))
            ax.plot([i - 0.36, i + 0.36], [med, med], color=col, lw=2.0, zorder=5,
                    solid_capstyle="round")
            ax.annotate(f"{len(v):,} cells\n{len(dsub)} donors", xy=(i, -0.16),
                        xycoords=("data", "axes fraction"), ha="center", va="top",
                        fontsize=8, color=MUTED, annotation_clip=False)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([ds for _, ds in order], fontsize=9,
                           color=INK2, rotation=20, ha="right")
        ax.set_ylabel(ylab, fontsize=9.5)
        ax.set_title(sub_t, fontsize=10.5, color=INK2, loc="left", pad=8)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", visible=False)

    handles = [plt.Line2D([], [], marker="o", ls="", markersize=8,
                          markerfacecolor=CANCER_COLORS[c], markeredgecolor=SURFACE,
                          label=f"{c} cancer") for c in ("cervical", "lung")]
    axes[0].legend(handles=handles, frameon=False, fontsize=9.5,
                   loc="upper left", labelcolor=INK2)
    fig.suptitle("Tumour NK effector-gene expression: cervical vs lung",
                 fontsize=13.5, color=INK, x=0.007, ha="left", y=0.985)
    fig.text(0.007, 0.925,
             "one violin per dataset, not per cancer type — absolute values are "
             "not comparable across datasets; only distribution shape is. "
             "counts below each violin: cells / donors.",
             fontsize=9, color=MUTED, ha="left")
    fig.subplots_adjust(bottom=0.30, top=0.82, left=0.075, right=0.985, wspace=0.22)
    p = f"{FIGS}/fig2_cervical_vs_lung_tumour.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def main():
    parts = []
    for f in ("cervical_nk_effector_cells.parquet", "lung_nk_effector_cells.parquet"):
        p = f"{OUT}/{f}"
        if os.path.exists(p):
            parts.append(pd.read_parquet(p))
        else:
            print(f"[missing] {f}")
    cells = pd.concat(parts, ignore_index=True)
    print(f"cells: {len(cells)}")
    print(cells.groupby(["cancer", "tissue_group"], observed=True).size().to_string())

    donors = ec.summarise(cells)
    ds = ec.summarise(cells, score="score_ds").rename(columns={"median": "median_ds"})
    donors = donors.merge(ds[["dataset", "cancer", "tissue_group", "donor", "median_ds"]],
                          on=["dataset", "cancer", "tissue_group", "donor"], how="left")
    donors.to_csv(f"{OUT}/all_nk_effector_by_donor.csv", index=False)

    cerv = cells[cells.cancer == "cervical"]
    cerv_d = donors[donors.cancer == "cervical"]
    p1 = fig1(cerv, cerv_d)
    p2 = fig2(cells, donors)
    print("wrote", p1)
    print("wrote", p2)


if __name__ == "__main__":
    main()
