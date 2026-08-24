# NPC 单细胞占位审计 — 占位问题答案

**run date**: 2026-08-24  
**检索来源**: Europe PMC REST (Q1–Q5, Q7–Q8) + NCBI E-utilities db=gds (Q6)  
**included 记录数**: 38  |  **secondary_analysis**: 58  |  **excluded**: 2342

> Q7/Q8 是 spec §5 强制的 query 修复：anchor #15 (PMID 37607536) 在 Europe PMC 中非开放获取
> (isOpenAccess=N, 无 PMCID)，全文未被索引，任何以 TITLE_ABS:"nasopharyngeal" 为轴的 query
> 都不可能召回它。Q7/Q8 改以 pan-cancer 单细胞图谱语言为轴，召回后再按 NPC 样本证据筛选。

---

## Q-A 有几篇做了 NPC 外周血的单细胞，且 NK_depth ≥ subclustered？

**答：4 篇。**

| pmid / doi | year | modality | n_patients | NK_depth | NK_claim |
|---|---|---|---|---|---|
| 41992060 | 2026 | scRNA | NOT_FOUND | NK_is_theme | Higher levels of CD16+CD57+ NK cells in blood correlated with better patient outcomes |
| 39438442 | 2024 | CyTOF | 24 | subclustered | CD56+ NK cells resolved into phenotypic subpopulations by CyTOF |
| 37607536 | 2023 | scRNA | 10 | NK_is_theme | tumor-associated NK cells enriched in tumors show impaired anti-tumor functions |
| 33531485 | 2021 | scRNA;scTCR | 10 | subclustered | NK cells were decreased in the tumours compared to the PBMC |

PMID 列表：41992060, 39438442, 37607536, 33531485

**关键限定**（决定这个答案是「没人做过」还是「图谱里已经有了」）：

- 血 + NK 深度分析的记录里，只有 **PMID 41992060**(Br J Cancer 2026) 是以 NPC 外周血 NK 亚群
  为主题的原始研究；其摘要明确 `CD16+CD57+ NK cells in blood correlated with better patient outcomes`。
- **PMID 37607536**(anchor #15, Cell 2023 泛癌 NK 图谱) 确实含 NPC 循环 NK：从 Zenodo record
  8275845 的 `comb_CD56_CD16_NK_blood.h5ad` 直接读出 —— `meta_histology` 取值
  `Nasopharyngeal Carcinoma(NPC)` 共 **9512 个细胞 / 10 位病人 / 10 个样本**，
  全部来自单一数据集 **GSE162025**，即 anchor #1 (PMID 33531485) 的 10 对瘤-血配对样本。
  这些细胞在图谱里被细分到 13 个 NK 亚群（最大三个：CD56dimCD16hi-c1-CX3CR1 5177、CD56dimCD16hi-c7-LAG3 2087、CD56dimCD16hi-c0-IL32 637）。
- 因此：NPC 循环 NK 的**单细胞数据**在公开图谱里已经存在并已被亚群化，但其**原始样本**
  仍只有 GSE162025 这一份 10 人队列；没有第二份独立的 NPC 外周血 NK 单细胞队列。

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

## Q-B 有几篇纳入了远处转移灶（M1）样本？其中几篇同时有血？

**答：5 篇纳入远处转移灶样本；其中 1 篇同时有血。**

| pmid / doi | year | tissue_distant_met | M_stage | tissue_blood | n_patients |
|---|---|---|---|---|---|
| 42230534 | 2026 | Y (site: not specified - metastatic vs non-metastatic NPC co | M1_included | N | NOT_FOUND |
| 42456326 | 2026 | Y (site: liver) | M1_included | N | NOT_FOUND |
| 41619722 | 2026 | Y (site: bone/spine) | M1_included | N | 3 |
| 40265092 | 2025 | Y (site: cervical lymph node metastasis) | unstated | N | 1 |
| 36739462 | 2023 | Y (site: liver/lung/bone/other - de novo metastatic) | M1_included | Y | 2 |

PMID 列表（M1）：42230534, 42456326, 41619722, 40265092, 36739462
PMID 列表（M1 且有血）：36739462

注：PMID 36739462 的 163 份配对血+组织样本是**基因组**层面（WES/WGS）；其单细胞部分只有
2 位病人的 11 个样本，血未进入单细胞。这是唯一一篇 M1 + 血同时出现的记录，但两者不在同一层。

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

## Q-C 有几篇对同一病人做了纵向重复取血的单细胞？时间点分别是什么？

**答：2 篇。**（另有 2 篇是纵向的但取的是组织而非血）

| pmid / doi | year | modality | 时间点（原文） |
|---|---|---|---|
| doi:10.1101/2025.11.12.687945 | 2025 | scRNA | Y (before and after 18 Gy RT; leukocyte counts at week 1-2, 3-4, 5-6) |
| 39438442 | 2024 | CyTOF | Y (baseline and during treatment) |

PMID 列表：doi:10.1101/2025.11.12.687945, 39438442

纵向但非血的记录：39415331, 37280275

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

## Q-D 有几篇把免疫单细胞状态对连续血浆 EBV DNA 载量做了回归（而非阳性/阴性二分）？

**答：1 篇。**

| pmid / doi | year | EBV_DNA_link | EBV 变量类型 | 原文依据 |
|---|---|---|---|---|
| 32686767 | 2020 | continuous | plasma_EBV_DNA | We further explored their correlations with plasma EBV DNA concentrations. |

PMID 列表：32686767

其余 included 记录的 EBV_DNA_link 分布：categorical 6 篇、none 31 篇。

**`categorical` 的两种子类型必须分开看** —— spec 的字段名是 `EBV_DNA_link`，但被判 categorical
的 6 篇里只有一部分真的用了**血浆 EBV DNA**，其余用的是**肿瘤组织 EBV 状态**(EBER ISH /
EBV-high vs EBV-low)，两者不是同一个变量：

| pmid | EBV 变量类型 | 原文依据 |
|---|---|---|
| 42245668 | tumour_EBV_status | EBER in situ hybridization (ISH) was performed on formalin-fixed, paraffin-embedded (FFPE) tissue specimens |
| 40059116 | tumour_EBV_status | NPC scRNA-seq data was retrieved from a published study, which classified NPC cells into EBV-high and EBV-low |
| 39231979 | tumour_EBV_status | All NPC tumours were EBV positive, as confirmed using in situ hybridisation of EBV encoded small RNAs (EBERs) in tumour biopsy. |
| 39415331 | plasma_EBV_DNA | one of the patients with negative plasma EBV DNA at the baseline (P5) showed decreased CD8T_IFN cells on the treatment with induct |
| 33531485 | tumour_EBV_status | malignant cells with different Epstein-Barr virus infection status |
| 32901110 | tumour_EBV_status | we performed single-cell RNA sequencing on ~104,000 cells from 19 EBV + NPCs and 7 nonmalignant nasopharyngeal biopsies |

即：**真正把单细胞结果连到血浆 EBV DNA 的只有 2 篇**（39415331, 32686767），其中只有 1 篇是连续量。

注：被判为 `continuous` 的唯一一篇是 PMID 32686767（Cell Res 2020），其全文写明
`We further explored their correlations with plasma EBV DNA concentrations.`——把免疫细胞
signature 对**连续**血浆 EBV DNA 浓度作相关，而非 Sero+/Sero− 二分。作为对照，anchor #12
(PMID 36860875) 全篇以 EBV DNA **血清阳性/阴性**分组，属 `categorical`；但该记录在本轮
被判为 secondary_analysis（其 scRNA-seq 数据取自 GSE150430），故不在 included 表内。

**没有任何一篇**把免疫单细胞状态对连续血浆 EBV DNA 拷贝数做**回归建模**；32686767 做的是
相关性探索。若把标准严格限定为「回归」，答案是 0 篇。

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

## Q-E 有几篇明确用 metastatic cascade 分阶段的框架来组织免疫学结论？

**答：严格口径 1 篇。**（下面同时给出较宽的结构性口径）

**严格口径** —— 检索文本里出现明确的 cascade / 分阶段路径措辞：

| pmid / doi | year | 原文依据 |
|---|---|---|
| 36739462 | 2023 | To track the evolutionary route of metastasis, here we perform an integrative genomic analysis of 163 matched blood and primary, regional lymph node m |

PMID 列表：36739462

**结构性口径** —— 不要求用 cascade 这个词，只要求同一研究在单细胞层面同时纳入
解剖学上 ≥2 个阶段（原发灶 / 区域淋巴结转移 / 远处转移），共 5 篇：

| pmid / doi | year | primary | LN_met | distant_met | blood |
|---|---|---|---|---|---|
| 42230534 | 2026 | Y | N | Y (site: not specified - met | N |
| 42456326 | 2026 | Y | N | Y (site: liver) | N |
| 42625209 | 2026 | Y | Y | N | N |
| 40265092 | 2025 | N | Y | Y (site: cervical lymph node | N |
| 36739462 | 2023 | Y | Y | Y (site: liver/lung/bone/oth | Y |

PMID 列表：42230534, 42456326, 42625209, 40265092, 36739462

两个口径的差距本身就是答案：**采样上跨越 ≥2 个转移阶段的有 5 篇，但只有 1 篇真的把结论
按 cascade 阶段组织**；其余都是两点对比（原发 vs 转移），没有分阶段框架。

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

## Q-F 公开可下载、且含外周血的 NPC 单细胞数据集有哪些？

**答：2 个。**

| accession | repository | n_samples | 组织类型（据 GEO sample title） | 关联 PMID |
|---|---|---|---|---|
| GSE162025 | GEO | 40 | primary tumour + PBMC | originating: 33531485 (anchor #1, Nat Commun 2021); re-used by: 37607536,40818459 |
| GSE202604 | GEO | 8 | Single-cell RNA sequencing of right atrial blood mononuclear | originating: 无关联论文（Europe PMC 检索 'GSE202604' 命中 0 条，属仅存 GEO 的数据集）; re-used by: - |

- **GSE162025** — Tumour Heterogeneity and Intercellular Networks of Nasopharyngeal Carcinoma at Single Cell Resolution  
  sample titles: `NPC_SC_1807_Tumor_TCR NPC_SC_1807_Tumor_cDNA NPC_SC_1811_Tumor_TCR NPC_SC_1811_Tumor_cDNA NPC_SC_1810_PBMC_TCR NPC_SC_1810_PBMC_cDNA NPC_SC_1816_Tumor_cDNA NPC_SC_1815_PBMC_cDNA NPC_SC_1815_PBMC_TCR NPC_SC_1816_Tumor_TCR NPC_SC_1806_Tumor_cDNA NPC_SC_1808_PBMC`
- **GSE202604** — Single-cell RNA sequencing of right atrial blood mononuclear cells offers a novel approach to studying nasopharyngeal carcinoma.  
  sample titles: `NPC_TM437 NPC_TM419 NPC_TM440 NPC_TM438 NPC_TM420 NPC_TM439 NPC_TM421 NPC_TM418`

**注意 GSE202604 的血不是外周血**：其 series title 写明取的是 *right atrial blood*
（右心房血，术中/导管取样），不是肘静脉外周血。若把 Q-F 严格限定为**外周血**，
公开可下载的 NPC 单细胞数据集只有 **GSE162025 一个**。该数据集同时也是 anchor #15
泛癌 NK 图谱中全部 NPC 循环 NK 细胞的唯一来源 —— 也就是说，
**目前全世界公开可下载的 NPC 外周血单细胞数据，追溯到底只有这一份 10 人队列。**

**同样含血、但不是自由下载的**（受控访问，需申请）：

| accession | repository | 关联 PMID | 说明 |
|---|---|---|---|
| HRA006885 | GSA-Human | 39231979 | anchor #11 的 77 份活检+血 scRNA-seq 原始数据；GSA-Human 为受控访问 |
| HRA008590 | GSA-Human | 40315843 | TIL-ACT 试验的 scRNA/scTCR；含 PBMC 分离流程 |

**重要更正**：GSE206245（anchor #11 的 GEO 条目）的 series summary 复述了论文摘要中的
「77 biopsy and blood samples」，但其**实际存放的 sample 全部是空间转录组**（NPC_ST1–ST19）
加一条 BCR 记录，**不含血**。本审计的 `includes_blood` 一律以 sample title / series title
判定，不采信 series summary，正是为避免这一类误判。

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。

---

## anchor #15 的额外核实（spec §5 第 15 项）

Zenodo record **8275845**（`10.5281/zenodo.8275845`）:

- 文件 `comb_CD56_CD16_NK_blood.h5ad`：`meta_tissue` 唯一取值 `Blood`，共 **84904** 个循环 NK 细胞。
- 其中 **Nasopharyngeal Carcinoma(NPC) = 9512 个细胞 / 10 位病人 / 10 个样本**，
  patient ID 为 `P1805`, `P1806`, `P1807`, `P1808`, `P1810`, `P1811`, `P1813`, `P1815`, `P1816`, `Ppbmc`。
- NPC 细胞的来源数据集：**GSE162025**（该文 24 个瘤种确含 nasopharyngeal cancer，此为直接证据）。
- NPC 循环 NK 在图谱中的亚群分布：

| NK subset | cells |
|---|---|
| CD56dimCD16hi-c1-CX3CR1 | 5177 |
| CD56dimCD16hi-c7-LAG3 | 2087 |
| CD56dimCD16hi-c0-IL32 | 637 |
| CD56dimCD16hi-c3-NFKBIA | 552 |
| CD56brightCD16lo-c0-FGFBP2-RGS1lo | 416 |
| CD56dimCD16hi-c2-KLRF1 | 357 |
| CD56brightCD16lo-c4-SPTSSB | 86 |
| CD56dimCD16hi-c4-MKI67 | 59 |
| CD56brightCD16lo-c2-CCL3 | 50 |
| CD56brightCD16lo-c1-SPTSSB-RGS1lo | 50 |
| CD56dimCD16hi-c6-CREM | 25 |
| CD56brightCD16hi | 15 |
| CD56brightCD16lo-c5-CREM | 1 |

全部瘤种在该 blood 对象中的占比：

| histology | cells | patients |
|---|---|---|
| Breast Cancer(BRCA) | 26656 | 23 |
| Healthy donor | 19143 | 70 |
| Melanoma(MELA) | 14501 | 5 |
| Nasopharyngeal Carcinoma(NPC) | 9512 | 10 |
| Head and Neck Squamous Cell Carcinoma(HNSCC) | 6078 | 41 |
| Colorectal Cancer(CRC) | 3233 | 12 |
| Pancreatic Cancer(PACA) | 2038 | 15 |
| Renal Carcinoma(RC) | 1022 | 7 |
| Lung Cancer(LC) | 893 | 6 |
| Hepatocellular Carcinoma(HCC) | 768 | 4 |
| Gastric Cancer(GC) | 613 | 2 |
| Intrahepatic cholangiocarcinoma(ICC) | 447 | 1 |

同一 record 的另一个文件 `comb_CD56_CD16_NK.h5ad`（全部组织，142304 个细胞）中，NPC 共 **10205 个细胞 / 23 位病人**，按组织拆开：

| meta_tissue | cells | patients |
|---|---|---|
| Blood | 9513 | 10 |
| Tumor | 692 | 22 |

即：这份泛癌 NK 图谱里 NPC 的 NK 细胞**绝大多数来自血**（9513/10205），
瘤内 NK 只有 692 个、分散在 22 位病人 —— 这与 scRNA-seq 对瘤内 NK 捕获率低的已知现象一致，
也说明「图谱里已经有了」这句话对**血**成立、对**瘤内 NK** 只在极稀疏的意义上成立。

（两个文件的 NPC 血细胞数相差 1 个细胞：blood 专用对象 9512，全组织对象 9513。
此处如实记录，未做调和。）

本结论基于 2026-08-24 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。
