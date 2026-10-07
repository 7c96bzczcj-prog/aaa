"""Spec §8.8: one dot plot per dataset; x = unified cell type, size = pct_expressing, colour = mean_logcpm.

Human and mouse go to separate folders and are never drawn together. Values are donor medians from
il18_summary.tsv. The colour scale is per dataset, since logCPM is not comparable across datasets (§7.4).
Hollow dots = epithelial classes in ambient-uncorrected data (excluded from conclusions, §4.3).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from phase2_aggregate import UNIFIED

ROOT = Path(__file__).resolve().parents[1]
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]  # sequential ramp
CMAP = LinearSegmentedColormap.from_list("seq_blue", BLUE)
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e6e6e6"


def main():
    s = pd.read_csv(ROOT / "results" / "il18_summary.tsv", sep="\t", comment="#", keep_default_na=False,
                    na_values=[""])
    s = s[s.celltype != "other"].copy()
    s["condition_detail"] = s.condition_detail.astype(str).replace("nan", "")
    s["row"] = s.condition + s.condition_detail.map(lambda d: f" | {d}" if d else "")
    for (sp, ds), g in s.groupby(["species", "dataset_id"]):
        rows = sorted(g.row.unique(), key=lambda r: (not r.startswith("resting"), _num(r)))
        fig, ax = plt.subplots(figsize=(8.5, 1.2 + 0.42 * len(rows)))
        xi = {c: i for i, c in enumerate(UNIFIED)}
        yi = {r: i for i, r in enumerate(rows)}
        vmin, vmax = g.logcpm_median.min(), g.logcpm_median.max()
        filled = g[~g.epithelial_excluded_from_conclusions]
        hollow = g[g.epithelial_excluded_from_conclusions]
        sc = ax.scatter(filled.celltype.map(xi), filled.row.map(yi), s=20 + 600 * filled.pct_median,
                        c=filled.logcpm_median, cmap=CMAP, vmin=vmin, vmax=vmax, edgecolors="white", linewidths=1)
        ax.scatter(hollow.celltype.map(xi), hollow.row.map(yi), s=20 + 600 * hollow.pct_median,
                   facecolors="none", edgecolors=MUTED, linewidths=1)
        ax.set_xticks(range(len(UNIFIED)), UNIFIED, rotation=45, ha="right", fontsize=8, color=INK)
        ax.set_yticks(range(len(rows)), rows, fontsize=8, color=INK)
        ax.set_xlim(-0.6, len(UNIFIED) - 0.4)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for sp_ in ax.spines.values():
            sp_.set_visible(False)
        ax.set_title(ds, fontsize=9, color=INK, loc="left")
        cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
        cb.set_label("IL18 transcript (pro-IL-18)\ndonor-median log2(CPM+1)", fontsize=7, color=MUTED)
        cb.ax.tick_params(labelsize=7)
        for p in [0.1, 0.3, 0.6]:
            ax.scatter([], [], s=20 + 600 * p, color=MUTED, label=f"{int(p * 100)}%")
        ax.legend(title="% cells\nIL18 transcript+", fontsize=7, title_fontsize=7, frameon=False,
                  loc="upper left", bbox_to_anchor=(1.12, 1), labelspacing=1.2)
        out = ROOT / "figures" / sp
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{ds.replace(':', '__')}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def _num(r):
    try:
        return float(r.split("|")[-1])
    except ValueError:
        return 0.0


if __name__ == "__main__":
    main()
