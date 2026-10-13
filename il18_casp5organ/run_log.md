# run_log.md — 跨器官 IL18/CASP1/CASP3 筛选运行记录

Census LTS 2025-11-08；cellxgene_census 1.18.0 / tiledbsoma 2.3.0；随机种子 20261013。检索与核对时间：2026-10-13。

## 基因 ID 解析与比对（§1.1）
六个 Ensembl ID 经 Ensembl REST (`xrefs/symbol`) 解析并与 Census feature 表交叉核对，**全部一致，无不一致项**。详见 `results/gene_id_resolution.md`。取基因一律用 Ensembl ID。

## §12 开工前四项确认
1. 网络：CELLxGENE / EMBL-EBI / GEO(E-utilities) / Ensembl / brainrnaseq.org 均可访问。
2. cellxgene_census 1.18.0；阶段一用 get_anndata(organism, X_name="raw", obs_value_filter, var_value_filter="feature_id in [...]", obs_column_names)，含 suspension_type 字段。
3. 六个基因 ID 全部解析一致（见上）。
4. 磁盘可写约 28 GB；人脑 2900 万细胞不一次性入内存，逐数据集处理；脑按用户确认仅取参考图谱。

## 环境 RNA 校正
本机无 GPU（CellBender 不可行）；Census 源矩阵已过滤、无空液滴（SoupX 不可行）。经用户确认（①B），全部数据标 ambient_uncorrected=TRUE，候选为临时/待验证结果。

## 每分层细胞数、供体数、基因命中（normal 细胞；disease≠normal 已在取数时排除）
| 分层 id | accession | 物种 | 器官 | susp | 平台 | 供体 | 细胞 | 基因命中 | gene_id缺失 |
|---|---|---|---|---|---|---|---|---|---|
| H_brain_nucleus | e1f595f6-ba2c-495e-9bee-7056f116b1e4 | human | brain | nucleus | snRNA_10x | 3 | 107301 | IL18/CASP1/CASP3 | False |
| H_liver_cell | e84f2780-51e8-4cfa-8aa0-13bbfef677c7 | human | liver | cell | 10x_3p | 14 | 152535 | IL18/CASP1/CASP3 | False |
| H_liver_nucleus | ddb22b3d-a75c-4dd1-9730-dff7fc8ca530 | human | liver | nucleus | snRNA_10x | 4 | 43863 | IL18/CASP1/CASP3 | False |
| H_lung_cell | 9f222629-9e39-47d0-b83f-e08d610c7479 | human | lung | cell | 10x_3p | 163 | 575893 | IL18/CASP1/CASP3 | False |
| H_lung_nucleus | 3de0ad6d-4378-4f62-b37b-ec0b75a50d94 | human | lung | nucleus | snRNA_10x | 9 | 46500 | IL18/CASP1/CASP3 | False |
| H_pancreas_cell | 53d208b0-2cfd-4366-9866-c3c6114081bc | human | pancreas | cell | 10x_3p | 3 | 14140 | IL18/CASP1/CASP3 | False |
| H_pancreas_nucleus | 31f657dc-1875-4c4b-a5ca-ce63b3ef3a82 | human | pancreas | nucleus | snRNA_10x | 5 | 79367 | IL18/CASP1/CASP3 | False |
| H_spinal_nucleus | 090da8ea-46e8-40df-bffc-1f78e1538d27 | human | spinal_cord | nucleus | snRNA_10x | 1 | 24190 | IL18/CASP1/CASP3 | False |
| M_brain_cell | d7291f04-fbbb-4d65-990a-f01fa44e915b | mouse | brain | cell | 10x_3p | 55 | 1169213 | IL18/CASP1/CASP3 | False |
| M_brain_nucleus | 1229ecc2-b067-4664-91da-0251aec31574 | mouse | brain | nucleus | snRNA_10x | 10 | 205722 | IL18/CASP1/CASP3 | False |
| M_liver_cell | 48b37086-25f7-4ecd-be66-f5bb378e3aea | mouse | liver | cell | 10x_3p | 12 | 7294 | IL18/CASP1/CASP3 | False |
| M_lung_cell | 48b37086-25f7-4ecd-be66-f5bb378e3aea | mouse | lung | cell | 10x_3p | 16 | 24540 | IL18/CASP1/CASP3 | False |
| M_pancreas_cell | 48b37086-25f7-4ecd-be66-f5bb378e3aea | mouse | pancreas | cell | 10x_3p | 8 | 6201 | IL18/CASP1/CASP3 | False |
| M_pancreas_islet_cell | 49e4ffcc-5444-406d-bdee-577127404ba8 | mouse | pancreas | cell | 10x_3p | 26 | 134590 | IL18/CASP1/CASP3 | False |
| M_spinal_nucleus | GSE103892 | mouse | spinal_cord | nucleus | snRNA_dropseq | 6 | 3554 | IL18/CASP1/CASP3 | True |

## 分层选择与范围（用户确认）
- 肺/肝/胰：取 Census 中成年正常、作者注释的参考图谱（HLCA、Guilliams、Tabula Sapiens、LungMAP、健康肝整合、胰腺协调参考、小鼠 Tabula Muris Senis、小鼠胰岛图谱）。
- 脑：人用 Siletti Human Brain Cell Atlas v1.0（MTG dissection, nucleus）；小鼠用 Yao 2023（皮层+海马, cell）与皮层 WT 单核（nucleus）。Census 人脑 180+ 数据集，经用户确认不全扫，仅取参考图谱。
- 脊髓：人用 Siletti 脊髓 dissection（nucleus，仅 1 供体 → 仅描述统计）；小鼠 Census 无数据，GEO 核对后取 GSE103892。

## 数据集核对（§4.1）
- 全部 Census 分层：以 Census obs 查询确认 dataset_id 存在、物种正确、器官正确、suspension_type 可定、normal 细胞存在（2026-10-13）。
- GSE103892（小鼠脊髓）：GEO SOFT + FTP suppl 核对——Mus musculus、腰段脊髓、单核；原始计数矩阵（24 MB）+ 作者 cell.type 注释（112 KB）；取 naive 基线样本 f1/m1/M4/M5/F3/F4（排除 *cx 皮层对照、formalin、rotarod）。矩阵仅基因符号、无 Ensembl ID → gene_id_unavailable=TRUE，仅参考不进候选。

## 排除与缺口
- 人成年正常「脑 cell」：Census 无（大样本为胎儿/类器官），记为缺口；人脑仅 nucleus。
- 人脊髓 nucleus：仅 1 供体，n_donors<3，不产候选（§10.3），仅描述统计。
- 小鼠脊髓：Census 无 → GSE103892（符号，gene_id_unavailable）。另核对的 GSE158380（Russ，单细胞）原始矩阵 >50 MB 且 GEO 上无逐细胞注释，未用。拒绝：GSE172167（损伤）、GSE99590（芯片）等。
- GSE103892 因 gene_id_unavailable=TRUE 不进候选集；其结果仅作小鼠脊髓的参考呈现。

## 过滤
n_cells<20 的「供体×细胞类型」格子置 NA（§10.2）；n_donors<3 的细胞类型不进候选（§10.3）；候选还需 pct_IL18≥0.10（§8.4）。
