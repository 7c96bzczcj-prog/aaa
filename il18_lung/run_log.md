# run_log.md — IL18 肺细胞类型分析运行记录

Census LTS `2025-11-08`；`cellxgene_census` 1.18.0 / `tiledbsoma` 2.3.0；随机种子 20261007。
基因严格以 Ensembl ID 取：人 `ENSG00000150782`、鼠 `ENSMUSG00000039217`。两物种在每个数据集的 var/features 中均命中（`gene_id_unavailable=False`），无静默跳过。

## 阶段一 — Census 全库扫描

- 筛选 `tissue_general=='lung' and is_primary_data==True`，取 `raw` 层。
- 人：6,167,731 primary 肺细胞，45 数据集；鼠：220,706 primary 肺细胞，4 数据集。
- Census 疾病标签覆盖：人有 normal / COVID-19 / pneumonia，**无 influenza**；鼠有 normal / influenza，**无细菌、无 LPS**。
- 输出 `results/census_scan.tsv`。此阶段仅用于定位候选，数值不进统计。

## 环境可用性（规格 §9）

- 网络：CELLxGENE（S3 + API）、EMBL-EBI（GXA、ArrayExpress）、GEO（网页/E-utilities/FTP）、Ensembl 均返回 200。
- 磁盘：会话配额约 30 GB 可写；内存 15 GB；4 核；无 GPU。
- 因此 CellBender（需 GPU）不可行；SoupX 需空液滴矩阵，多数数据集不提供——故除作者已校正者外，10x 数据标 `ambient_uncorrected=TRUE`，其上皮结果不进结论（§4.3）。

## 阶段二 — 各数据集处理

| 数据集 | accession | 物种 | 条件 | 平台 | 组织 | 供体 | 细胞 | 格子 | NA格子 | 环境RNA |
|---|---|---|---|---|---|---|---|---|---|---|
| H01_HLCA_core:Banovich_Kropski_2020 | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 38 | 121894 | 602 | 247 | 未校正 |
| H01_HLCA_core:Lafyatis_Rojas_2019 | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 6 | 24181 | 102 | 32 | 未校正 |
| H01_HLCA_core:Misharin_2021 | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 2 | 55426 | 32 | 3 | 未校正 |
| H01_HLCA_core:Misharin_Budinger_2018 | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 8 | 41220 | 124 | 39 | 未校正 |
| H01_HLCA_core:Teichmann_Meyer_2019 | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 6 | 12231 | 91 | 46 | 未校正 |
| H02_IPF_atlas_controls | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting | 10x_3p | lung_parenchyma | 28 | 95303 | 411 | 179 | 未校正 |
| H03_TabulaSapiens_lung | CELLxGENE 53d208b0-2cfd-4366-9866-c3c6114081bc | human | resting | 10x_3p | lung_parenchyma | 4 | 62954 | 67 | 15 | decontX(作者) |
| H04_Travaglini_10x | CELLxGENE 8c42cfd0-0b0a-46d5-910c-fc833d83c45e | human | resting | 10x_3p | lung_parenchyma | 3 | 60993 | 48 | 3 | 未校正 |
| H05_CrossTissueImmune_lung | CELLxGENE 1b9d8702-5af8-4142-85ed-020eb06ec4f6 | human | resting | 10x_5p | lung_parenchyma | 10 | 35419 | 115 | 59 | 未校正 |
| H06_LungMAP_CellRef | CELLxGENE eb499fd8-7000-419f-8854-926b9b61b11f | human | resting | 10x_3p | lung_parenchyma | 21 | 122445 | 299 | 113 | 未校正 |
| H07_Melms_COVID_autopsy | CELLxGENE d8da613f-e681-4c69-b463-e94f5e66847f | human | resting,sars_cov_2 | 10x_3p | lung_parenchyma | 27 | 116313 | 442 | 85 | 未校正 |
| H08_Delorey_COVID_autopsy | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | sars_cov_2 | 10x_3p | lung_parenchyma | 16 | 96060 | 266 | 47 | 未校正 |
| H09_Bharat_COVID_autopsy | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | sars_cov_2 | 10x_3p | lung_parenchyma | 3 | 91980 | 47 | 2 | 未校正 |
| H10_Liao_BAL | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | resting,sars_cov_2 | 10x_5p | bal | 12 | 62469 | 78 | 23 | 未校正 |
| H11_Wauters_BAL | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | pneumonia_unspecified,sars_cov_2 | 10x_5p | bal | 35 | 65166 | 347 | 134 | 未校正 |
| H12_Grant_BAL | CELLxGENE 9f222629-9e39-47d0-b83f-e08d610c7479 | human | sars_cov_2 | 10x_5p | bal | 10 | 77146 | 104 | 36 | 未校正 |
| H13_Wu_exvivo_SARS2 | CELLxGENE 2ac76f1b / 7f7faf6b / 7b55fe5c / 7c4 | human | resting,sars_cov_2 | 10x_3p | lung_parenchyma | 4 | 271035 | 175 | 23 | 未校正 |
| M01_TMS_10x_lung | CELLxGENE 48b37086-25f7-4ecd-be66-f5bb378e3aea | mouse | resting | 10x_3p | lung_parenchyma | 16 | 24540 | 211 | 91 | 未校正 |
| M02_TMS_SS2_lung | CELLxGENE 98e5ea9f-16d6-47ec-a529-686e76515e39 | mouse | resting | smartseq2 | lung_parenchyma | 14 | 5218 | 166 | 113 | 不适用(plate) |
| M03_PostFlu_timeseries | CELLxGENE 89840cec-5734-455d-b68e-34d27650e148 | mouse | influenza,resting | 10x_3p | lung_parenchyma | 25 | 123154 | 446 | 96 | 未校正 |
| M04_GSE236342_Spn | GSE236342 | mouse | bacterial,resting | 10x_3p | lung_parenchyma | 3 | 16070 | 42 | 8 | 未校正 |
| M05_GSE313334_Abaumannii | GSE313334 | mouse | bacterial,resting | 10x_5p | lung_parenchyma | 2 | 4461 | 18 | 5 | 未校正 |

## Ensembl ID 命中记录（逐数据集，来自 phase2_build.log）

Census 数据集：`get_anndata` 以 `var_value_filter=feature_id==<ENSG/ENSMUSG>` 取基因，返回 var 的 feature_id 与请求一致、feature_name 为 IL18/Il18 时才通过。
GEO 数据集：在 features.tsv / h5 的 features/id 中精确匹配 Ensembl ID，各样本命中 1 行。

## 数据集核对（规格 §2.1）

- **GSE236342**（肺炎链球菌）：本会话拉取 GEO series + 8 个 GSM SOFT 文本确认——Mus musculus、lung、S.pn. vs PBS、10x；附 `GSE236342_cell_annotation.csv.gz`。仅用 WT 臂（PBS/12h/36h）；60 hpi 无注释，未用。
- **GSE313334**（鲍曼不动杆菌）：GEO series + 4 GSM SOFT 确认——Mus musculus、lung、infected vs mock、10x 5'；附 `GSE313334_cell_metadata.csv.gz`。barcode 后缀按样本显式映射（ES1`_1_1`/ES2`_2_1`/ES3`_1`/ES4`_2`），8541 细胞全部匹配到作者注释。悬液富集 CD45−，注释以免疫细胞为主，故仅用 no-antibiotics 的 mock 与感染臂。
- **GSE107947**（Steuerman 2018 流感双物种）：已确认存在，但其设计文件只有 per-well 表、无细胞类型注释，按「用作者注释、不重新聚类」本轮未纳入；流感由 Census 中有注释的 Post-Flu Time Series 覆盖。

## 排除的数据集及原因

- 人甲型流感肺组织：几乎不存在；唯一候选 GSE268542（单 COPD 供体离体感染，注释在未读 RDS）未纳入主表 → 记为数据缺口。
- 小鼠 LPS 急性肺损伤单细胞（GSE280364 / GSE280611 / GSE303253 等）：全肺、有对照，但 GEO 无作者逐细胞注释文件 → 按 §4.2 本轮排除；以 ImmGen 分选 LPS 数据作旁证。
- 仅体外细胞系、仅 BALF 上清、无细胞类型注释的研究：按 §2.1 排除。

## Census 内按作者来源研究拆分（HLCA）

HLCA（dataset 9f222629）被拆为作者 `study` 子研究以避免跨研究混池：
- 静息主力取 `core` 且 `Healthy` 且 `donor_lung`/`surgical_resection` 的 5 个子研究（H01）。
- COVID 尸检单独取 Regev_2021（Delorey，H08）、Budinger_2020（Bharat，H09）。
- COVID/肺炎 BALF 取 Zhang_2021（Liao，H10）、Lambrechts_2021（Wauters，H11）、Wunderink_2021（Grant，H12）。
- IPF Cell Atlas 对照取 Kaminski_2020 的 Healthy（H02）。

## 追加：小鼠 LPS 急性肺损伤（用户授权纳入，自行注释——对规格的偏离）

- M06 GSE280364、M07 GSE280611：全肺 10x，各 2 对照 vs 2 LPS（GSE280611 的 LPS+DMH1 药物臂已排除）。Il18 以 Ensembl ID 命中。
- 因 GEO 无作者逐细胞注释，用 scanpy 标准流程（normalize→HVG→PCA→neighbors→Leiden res=1.0）聚类，再以犬齿清晰的 marker panel 对每个簇评分、按细胞 argmax 的簇内多数票归入 11 类（`scripts/lps_annotate.py`）。QC：counts≥500 且 genes≥200。随机种子 20261007。
- 这是对规格 §4.2「用作者注释、不重新聚类」的偏离；相关行在主表/汇总表标 `self_annotated=TRUE`，簇→类诊断见 `results/lps_cluster_labels_*.tsv`。
- 结果与其余数据一致：肺泡巨噬细胞 IL18 阳性 63–95% 居首；LPS 后肺泡巨噬细胞被耗竭（细胞数锐减）但存活者仍高表达。

## 幻灯片交付

`IL18_lung_celltype.pptx`（10 页，无装饰）：结论、数据与方法、人/小鼠总览散点图、上皮真伪、流感时间序列、零值判读、数据缺口与偏离。两张总览散点图（`figures/master_human.png`、`figures/master_mouse.png`）将每个单细胞数据集作为一行，点大小=阳性比例、颜色=数据集内排序。
