"""Two master dot plots (human, mouse): every dataset x condition row, 11 cell-type columns.

size = donor-median % cells IL18+ ; colour = within-dataset rank (1 = highest, cross-dataset safe, §7.4).
Decoration kept to a minimum. Human and mouse are separate figures (§4.5).
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
LABELS = ["Alv.Mac", "Mono-der.Mac", "DC", "Neutrophil", "T/NK", "B/Plasma", "AT1", "AT2",
          "Airway epi.", "Endothelial", "Fibro/Mesen."]
# rank 1 (highest) -> dark; rank 10/11 -> light
BLUE = ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"]
CMAP = LinearSegmentedColormap.from_list("rank", BLUE)
INK, MUTED, GRID = "#1a1a1a", "#777777", "#ececec"

COND = {"resting": "Resting", "influenza": "Influenza", "sars_cov_2": "SARS-CoV-2",
        "bacterial": "Bacterial", "lps_ali": "LPS-ALI", "pneumonia_unspecified": "Pneumonia (unspec.)"}
ROWLABEL = {
    "H01_HLCA_core:Banovich_Kropski_2020": "HLCA Banovich/Kropski", "H01_HLCA_core:Lafyatis_Rojas_2019": "HLCA Lafyatis/Rojas",
    "H01_HLCA_core:Misharin_2021": "HLCA Misharin", "H01_HLCA_core:Misharin_Budinger_2018": "HLCA Misharin/Budinger",
    "H01_HLCA_core:Teichmann_Meyer_2019": "HLCA Teichmann/Meyer", "H02_IPF_atlas_controls": "IPF-atlas controls",
    "H03_TabulaSapiens_lung": "Tabula Sapiens (decontX)", "H04_Travaglini_10x": "Travaglini (Krasnow)",
    "H05_CrossTissueImmune_lung": "Cross-tissue immune", "H06_LungMAP_CellRef": "LungMAP CellRef",
    "H07_Melms_COVID_autopsy": "Melms autopsy", "H08_Delorey_COVID_autopsy": "Delorey autopsy",
    "H09_Bharat_COVID_autopsy": "Bharat autopsy", "H10_Liao_BAL": "Liao BAL", "H11_Wauters_BAL": "Wauters BAL",
    "H12_Grant_BAL": "Grant BAL", "H13_Wu_exvivo_SARS2": "Wu ex-vivo lung",
    "M01_TMS_10x_lung": "Tabula Muris Senis 10x", "M02_TMS_SS2_lung": "Tabula Muris Senis SS2",
    "M03_PostFlu_timeseries": "Post-flu time series", "M04_GSE236342_Spn": "S. pneumoniae (GSE236342)",
    "M05_GSE313334_Abaumannii": "A. baumannii (GSE313334)", "M06_GSE280364_LPS_ALI": "LPS-ALI (GSE280364)*",
    "M07_GSE280611_LPS_ALI": "LPS-ALI (GSE280611)*"}
CONDORDER = {"resting": 0, "influenza": 1, "bacterial": 2, "lps_ali": 3, "sars_cov_2": 4, "pneumonia_unspecified": 5}


def rowkey(r):
    det = "" if str(r.condition_detail) in ("nan", "") else f" {r.condition_detail}"
    base = ROWLABEL.get(r.dataset_id, r.dataset_id)
    return f"{base} — {COND.get(r.condition, r.condition)}{det}"


def build(species, path):
    s = pd.read_csv(ROOT / "results" / "il18_summary.tsv", sep="\t", comment="#",
                    keep_default_na=False, na_values=[""])
    s = s[(s.species == species) & (s.celltype != "other")].copy()
    s["condition_detail"] = s.condition_detail.astype(str).replace("nan", "")
    # collapse flu time points to one representative (day 42, most donors) to keep rows readable
    if species == "mouse":
        s = s[~((s.dataset_id == "M03_PostFlu_timeseries") & (s.condition == "influenza") & (s.condition_detail != "42"))]
        s.loc[s.condition_detail == "42", "condition_detail"] = "d42"
    s["rk"] = (s.dataset_id + "||" + s.condition + "||" + s.condition_detail).map(
        lambda k: rowkey(s[(s.dataset_id + "||" + s.condition + "||" + s.condition_detail) == k].iloc[0]))
    s["ord"] = s.condition.map(CONDORDER) * 1000
    rows = (s[["rk", "ord", "dataset_id"]].drop_duplicates().sort_values(["dataset_id", "ord"]).rk.tolist())
    yi = {r: i for i, r in enumerate(rows)}
    xi = {c: i for i, c in enumerate(UNIFIED)}
    fig, ax = plt.subplots(figsize=(9.2, 0.9 + 0.34 * len(rows)))
    filled = s[~s.epithelial_excluded_from_conclusions]
    hollow = s[s.epithelial_excluded_from_conclusions]
    sc = ax.scatter(filled.celltype.map(xi), filled.rk.map(yi), s=10 + 430 * filled.pct_median.fillna(0),
                    c=filled.rank_median, cmap=CMAP, vmin=1, vmax=11, edgecolors="white", linewidths=0.8)
    ax.scatter(hollow.celltype.map(xi), hollow.rk.map(yi), s=10 + 430 * hollow.pct_median.fillna(0),
               facecolors="none", edgecolors="#bbbbbb", linewidths=0.9)
    ax.set_xticks(range(len(UNIFIED)), LABELS, rotation=40, ha="right", fontsize=8.5, color=INK)
    ax.set_yticks(range(len(rows)), rows, fontsize=8, color=INK)
    ax.set_xlim(-0.6, len(UNIFIED) - 0.4); ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.tick_params(length=0)
    ax.grid(color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.015)
    cb.set_label("IL18 transcript (pro-IL-18)\nrank within dataset  (1 = highest)", fontsize=7.5, color=MUTED)
    cb.ax.tick_params(labelsize=7, length=0); cb.outline.set_visible(False)
    for p in [0.1, 0.3, 0.6]:
        ax.scatter([], [], s=10 + 430 * p, color="#9ba7b4", label=f"{int(p*100)}%")
    ax.legend(title="% cells\nIL18 transcript+", fontsize=7.5, title_fontsize=7.5, frameon=False,
              loc="upper left", bbox_to_anchor=(1.09, 1.0), labelspacing=1.3, borderpad=0)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(species, len(rows), "rows ->", path)


def main():
    (ROOT / "figures").mkdir(exist_ok=True)
    build("human", ROOT / "figures" / "master_human.png")
    build("mouse", ROOT / "figures" / "master_mouse.png")


if __name__ == "__main__":
    main()
