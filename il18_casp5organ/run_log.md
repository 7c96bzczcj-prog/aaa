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

## 修订：判读重写与检出率口径更正（2026-10-13）

**（1）废弃 `zero_check.py`。** 该脚本从未打开参考文件、也从未使用候选自身的 CASP1 值，判读结果是细胞谱系标签的纯函数：星形胶质（CASP1 合并检出 2/3766）被判「可信」，肺泡巨噬（178/1341 = 13.3%，全部候选中最高）被判「dropout 伪影」。替代脚本 `evidence_integrate.py` 读取 `reference/ext/unified_reference.tsv` 与候选自身数值，IL18:CASP1 比值只在「同一来源 × 同一单位 × 同一状态」内计算，按 §9.5 只并列呈现、不做统计检验、不删除候选。

**（2）参考层扩展。** 3 → 约 10 个权威来源：Human Protein Atlas（含 IHC 蛋白层）、GTEx、Monaco 2019、Haemopedia、ImmGen 扩展面板、Allen Brain Map、brainrnaseq/Zhang、CELLxGENE WMG，以及新检索到的 FACS 分选 GEO 数据集（GSE80673 分选胰岛 α/β/δ、GSE156665、GSE264633、GSE222990、GSE226103、GSE254235、GSE99074、GSE67835）。统一参考表 1624 行，只有 exact/close 参与判读。扩展后，前一版因「无分选对照」而不可判读的 5 个实质细胞群全部变为可判读。

**（3）检出率口径更正。** 对外展示一律改用合并检出率（阳性细胞数 ÷ 分析细胞数），并标注分子、分母与有阳性的供体数；供体中位数单列保留。原因：每位供体仅数十到上百个该类型细胞，真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，中位数被截断为 0。例（CASP1）：

| 群体 | 阳性/分析细胞 | 合并 | 供体中位数 | 有阳性供体 |
|---|---|---|---|---|
| 人肺 肥大细胞 | 138 / 5299 | 2.60% | 0.00% | 21 / 44 |
| 小鼠脑 小胶质/CNS 巨噬 | 15 / 818 | 1.83% | 0.00% | 4 / 20 |
| 人肺 纤毛上皮（核） | 7 / 660 | 1.06% | 0.00% | 3 / 8 |
| 小鼠胰 导管细胞 | 2 / 1947 | 0.10% | 0.00% | 2 / 9 |
| 小鼠脑 星形胶质 | 2 / 3766 | 0.05% | 0.00% | 2 / 40 |
| 小鼠脑 少突胶质 | 2 / 8987 | 0.02% | 0.00% | 2 / 55 |

六行在旧图中全部显示 0.0%，真实值相差一百倍以上。

**（4）结论变化。** A 级（多源一致）小鼠脑星形胶质、少突胶质；B 级（潜在）人肺肥大细胞、小鼠胰岛 D(δ) 细胞、小鼠胰导管细胞，以及被分层内相对规则漏掉的人肺 AT1；C 级（CASP1 确实存在，非测不准）小鼠肺泡巨噬、人肺气道多纤毛上皮、小鼠脑小胶质、人肺纤毛上皮；D 级（方向相反）人肝嗜碱粒细胞、人肺树突状细胞；E 级（仍无匹配参考）小鼠脑 Cajal-Retzius 神经元。没有任何候选再以「你的零不可信」为由被否定。

**（5）物种差异。** A 级为小鼠特异：HPA 人脑单核数据中 IL18 压倒性集中于小胶质（338.5 nCPM），星形胶质 3.7、少突胶质 2.6；人小胶质 CASP1 41.2 nCPM（比值 8.2），不满足前提。人脑无对应 A 级候选。

## 判定收敛为三类 + 第五页重做（2026-10-13）

**三类判定。** 对外一律只分 **确定 / 可能 / 相反**：确定 = ≥2 个独立参考来源一致显示 CASP1 处于零底线或比 IL18 低 ≥20 倍；可能 = 只有单一参考来源支持，或本数据一致但尚无匹配参考；相反 = 参考层显示 CASP1 确实存在，或方向与前提相反。原六档内部判读（A/B/B′/C/D/E）保留在 `answer_table.tsv` 的 `tier` / `tier_cn` 列与 `il18_casp1_evidence.tsv` 中作为留痕，不再出现在图或结论里。映射：A→确定；B、B′（相对规则漏掉）、E（无匹配参考）→可能；C（CASP1 确实存在）、D（方向相反）→相反。结果为 确定 2 个、可能 7 个、相反 6 个。

**第五页重做。** 原图 `figures/mouse/brain_cell.png` 的总标题与第一个子图标题相互重叠，且为谱系级英文标签（只有 glia，看不到星形胶质与少突胶质）。新增 `scripts/answer_pack.py::draw_brain_dotplot()`，输出 `figures/mouse/brain_cell_dots.png`：小鼠脑单细胞、作者标签级、中文名、10 个细胞类型 × 3 个基因，点面积 = 合并检出率，颜色 = log2(CPM+1)（三基因各自独立色标），每个点旁直接标出检出率与「阳性细胞/分析细胞」。无标题重叠、无装饰。

## 扩充批次 v0.2 结果（2026-10-13）

分层 15 → 29，候选 17 → 32（去重后 32 行，确定 2 / 可能 21 / 相反 9）。人侧候选由 5 个增至 23 个。

**两处判定修正（扩充后新数据暴露出来的）**
1. `screen_missed.py` 原先只要求 IL18 进入分层内最高三分位，未要求 `detect_ok`（检出率 ≥10%）。人胰岛 β 细胞因此被列为「被相对规则漏掉」，而它的 IL18 检出率只有 1.0%——它不是被相对规则漏掉，是根本不表达。已加上检出率门槛，并按实际情况区分「仅 CASP1 相对规则」与「CASP1 相对规则 + CASP3 下限」两种落选原因。
2. 无匹配参考（E 档）原先一律映射为「可能」。人肝 CD56bright NK 细胞 CASP1 检出 25.2% 反而高于 IL18 的 19.1%，自身数据就与前提相反，却被判为「可能」。已改为：E 档行按自身数据判，若 IL18:CASP1 < 3 则判「相反」，并在 tier_reason 中写明依据的是自身数据而非参考层。

**四个此前无法回答的问题**
1. 人脑里 IL18 高的是神经元，不是胶质。PV 中间神经元（HVS 单细胞 78 供体 + PsychAD 单核 234 供体）与 L5 ET 投射神经元（HVS + SEA-AD MTG + SEA-AD DLPFC）IL18 检出 21–43%、CASP1 合并检出 0.07–0.22%，跨单细胞与单核两类技术、四套独立数据集重复。
2. 小鼠的星形胶质/少突胶质结论是真实物种差异，不是单核测序灵敏度问题。同为真单细胞 scRNA：人皮层星形胶质 IL18 检出 4.1%（分层内排 16/18，78 供体）、少突胶质 1.4%（17/18）；小鼠 55.4%（1/24）与 27.4%（2/24）。差 13–20 倍。
3. 人胰岛 δ 细胞不支持小鼠的发现。三套人胰岛数据（19/17/5 供体，δ 细胞 1217/4677/2889 个）中 δ 细胞 IL18 检出率均为 0.0–0.5%。
4. 人脊髓从「无结论」变为「明确阴性」。19 供体 15 类白质数据中 IL18 几乎只在小胶质（52.7%），而小胶质 CASP1 检出 10.5%、分层内排第 2；其余细胞类型 IL18 ≤3.3%。

**图**：新增 `n_strata`（同一细胞类型在几个独立分层里成为候选）作为与参考层相互独立的第二条证据轴；总览图按物种拆成 `answer_human.png` / `answer_mouse.png` 供 PPT 用（每个细胞类型取供体最多的一个分层作代表，不跨分层合并数值），合并版 `answer_v2.png` 供报告用。CASP1 高于 IL18 时标签自动换边，避免重叠。
