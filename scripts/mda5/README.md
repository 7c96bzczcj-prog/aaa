# MDA5 / IFIH1 analysis scripts

Run order (all expect the GEO/Census downloads described in
`docs/MDA5_IFIH1_CNS_AND_MS.md` under Methods):

| script | what it does |
|---|---|
| `libsize.py` | one full pass over the GSE180759 dense matrix to get per-nucleus total counts (needed for CP10k; the matrix ships no library sizes) |
| `final180759.py` | GSE180759: depth check, baseline per cell type, donor-level CP10k, depth-matched re-test, lesion-region gradient, RLR/ISG specificity, P2RY12+ composition control, ADAR panel |
| `lerma_analyze.py` | GSE279180: chunked extraction of the gene panel plus per-nucleus totals from the nine per-cell-type h5ad files |
| `census_ifih1.py` | CELLxGENE Census pull of IFIH1 across all human brain cells (normal + MS) |

Gene panel extraction from the GSE180759 matrix is a one-line awk filter on the
first CSV field; see the Methods section. Outputs land in `results/mda5/`.

## 图与 PPT 生成

| 脚本 | 作用 |
|---|---|
| `figstyle.py` | matplotlib 中文字体与统一配色（文泉驿正黑；病灶分期用有序蓝色梯度） |
| `mkfig12.py` | 图1 基线三技术对照；图2 MS vs 对照（小提琴+散点 与 供体水平森林图） |
| `mkfig345.py` | 图3 灰质/白质；图4 病灶分期与小胶质状态（含 UMAP 特征图）；图5 三项对照 |
| `../results/mda5/build_ppt.js` | 用 pptxgenjs 组装 9 页中文 PPT |

图输出到 `results/mda5/figs/`，PPT 为 `results/mda5/MDA5_IFIH1_人脑与MS.pptx`。
