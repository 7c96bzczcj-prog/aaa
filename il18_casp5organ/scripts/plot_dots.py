"""Spec §11.12: one figure per (species × tissue × suspension_type); three panels (IL18, CASP1, CASP3)
side by side. x = lineage cell type, dot size = detection rate (pct), colour = pseudobulk log2(CPM+1).
Species / tissue / suspension never combined. No decoration.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from phase2_aggregate import LINEAGES
from common import GENE_ORDER

ROOT = Path(__file__).resolve().parents[1]
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq", BLUE)
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#ececec"
LAB = [l for l in LINEAGES if l != "other"] + ["other"]


def main():
    s = pd.read_csv(ROOT / "results" / "il18_casp_summary.tsv", sep="\t", comment="#")
    s = s[s.level == "lineage"].copy()
    (ROOT / "figures").mkdir(exist_ok=True)
    for (sp, tis, susp), g in s.groupby(["species", "tissue", "suspension_type"]):
        present = [l for l in LAB if l in set(g.celltype)]
        if not present:
            continue
        yi = {c: i for i, c in enumerate(present)}
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 0.7 + 0.32 * len(present)), sharey=True)
        for ax, gene in zip(axes, GENE_ORDER):
            d = g[g.celltype.isin(present)]
            vmin, vmax = d[f"logcpm_{gene}"].min(), d[f"logcpm_{gene}"].max()
            sc = ax.scatter(d[f"logcpm_{gene}"] * 0 + 0, d.celltype.map(yi),
                            s=20 + 500 * d[f"pct_{gene}"].fillna(0), c=d[f"logcpm_{gene}"],
                            cmap=CMAP, vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1,
                            edgecolors="white", linewidths=0.8)
            ax.set_title(gene, fontsize=11, color=INK)
            ax.set_xticks([])
            ax.set_xlim(-0.5, 0.5)
            ax.grid(axis="y", color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            for sp_ in ax.spines.values():
                sp_.set_visible(False)
            cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.04)
            cb.set_label("log2(CPM+1)", fontsize=7, color=MUTED)
            cb.ax.tick_params(labelsize=6, length=0)
            cb.outline.set_visible(False)
        axes[0].set_yticks(range(len(present)), present, fontsize=8, color=INK)
        axes[0].set_ylim(len(present) - 0.5, -0.5)
        for p in [0.1, 0.3, 0.6]:
            axes[2].scatter([], [], s=20 + 500 * p, color=MUTED, label=f"{int(p*100)}%")
        axes[2].legend(title="% cells detected", fontsize=7, title_fontsize=7, frameon=False,
                       loc="center left", bbox_to_anchor=(1.25, 0.5), labelspacing=1.1)
        fig.suptitle(f"{sp} · {tis} · {susp}  —  IL18 / CASP1 / CASP3 (pseudobulk, lineage level)",
                     fontsize=10, color=INK, x=0.02, ha="left")
        out = ROOT / "figures" / sp
        out.mkdir(exist_ok=True)
        fig.savefig(out / f"{tis}_{susp}.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main()
