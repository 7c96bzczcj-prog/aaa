# NONTUMOR_KAEDE_SEARCH.md —— S-5 §1：非肿瘤 Kaede 对照的检索

> ## 结论
> **存在可用的非肿瘤 Kaede 对照，且其设计与 `E-MTAB-10176` 同型。**
> 主选：**`GSE228629`**（Kaede，IL-23 诱导的银屑病皮肤，光转换于 **−48 h 与 −24 h**，第 7 天取皮肤与关节，10x + 哈希标签）。
> 该数据集在**三个方面优于**肿瘤数据：颜色由哈希标签编码而**非与文库共线**；有**两个品系文库**可作复制；绿的时间窗（<24 h）与 `E-MTAB-10176` 同量级。
> **§0 的第三种结果（无可用数据）不成立**，本任务继续执行 §2、§3。

---

## 1. 检索路径与结果

### 1.1 已执行的检索

| # | 库 | 检索式 | 命中 | 状态 |
|---|---|---|---|---|
| 1 | Europe PMC | `Kaede photoconversion lymphocyte single-cell RNA` | 52 | 已阅前 25 |
| 2 | Europe PMC | `Kaede mice photoconversion migration transcriptome` | 60 | 已阅前 25 |
| 3 | Europe PMC | `AUTH:"Withers DR" AND Kaede` | **15** | 全阅（`E-MTAB-10176`／`GSE221064` 作者群） |
| 4 | GEO DataSets | `Kaede[All Fields] AND Mus musculus[Organism]` | 256 | 已阅前 60 条记录 |
| 5 | GEO DataSets | `photoconvertible[All Fields] AND Mus musculus[Organism]` | 33 | 全阅 |
| 6 | GEO DataSets | `Kaede photoconversion` | 36 | 全阅 |
| 7 | GEO DataSets | `KikGR[All Fields]` | 109 | 全阅（哺乳类相关者仅 4 条） |
| 8 | GEO DataSets | `31152090[PMID]`、`brachial lymph node Kaede`、`innate lymphoid cell photoconver* lymph node` | **0** | 见 §3 的 `NOT_AVAILABLE` |
| 9 | BioStudies/ArrayExpress | `Kaede photoconversion` | 441 | 前 20 条全为蛋白质光物理与斑马鱼／文昌鱼胚胎，与本任务无关 |

### 1.2 未执行的检索（`NOT_SEARCHED`，如实记录）

| 项 | 状态 | 说明 |
|---|---|---|
| SRA／ENA 的原始 FASTQ 层检索 | **`NOT_SEARCHED`** —— 本轮只在处理后矩阵层检索；若主选数据集不可用才需下沉到该层 |
| 中文文献库、预印本服务器（bioRxiv） | **`NOT_SEARCHED`** —— 未检索；预印本的数据可得性通常更差，且本轮已找到可用数据集 |
| 单细胞数据门户（CELLxGENE、Human Cell Atlas、Single Cell Portal） | **`NOT_SEARCHED`** —— 这些门户以人类数据为主，Kaede 为小鼠工具 |
| 光转换以外的时间分辨设计（活体标记后分选、TdT／Confetti 谱系追踪） | **`NOT_SEARCHED`** —— §1 允许作次选，但主选已找到，未再扩大范围 |

---

## 2. 找到的非肿瘤 Kaede／等价设计数据集（全部列出）

| 数据集 | 组织／模型 | 是否肿瘤 | 设计 | 平台 | 与 `E-MTAB-10176` 同型？ | 判定 |
|---|---|---|---|---|---|---|
| **`GSE228629`** | **银屑病皮肤 + 踝关节，IL-23 过表达** | **否** | 光转换**皮肤**（−48 h、−24 h），第 7 天取**皮肤**与关节；哈希标签分 6 组 | 10x，CellRanger 5.0.0 | **是**（皮肤室内绿 vs 红＝新到 vs 滞留） | **主选，已下载并分析** |
| `GSE235702` | 肠系膜淋巴结 | 否 | Dendra2 光转换，**7 天后**分选 D-Red⁺（滞留）vs D-Red⁻（新到）CD4⁺CD62L^lo | 10x | 部分——同型对比，但**间隔 7 天**且非 Kaede；颜色与文库共线（2 个样本） | **次选，未分析**（间隔差异过大，§3 会降级） |
| `GSE241820` | 痘苗病毒感染皮肤 → 引流淋巴结 | 否 | 光转换耳皮肤，24 h 后取**引流淋巴结**的 CD8 T | 10x | **否**——比较的是「自皮肤迁出者 vs 淋巴结本地细胞」，**不同室、极性相反** | 不用于 §2 主判 |
| `GSE216928` | 结肠 → 脾（Sci Immunol 2024 迁移图谱） | 否 | 光转换**降结肠**，24/48/72 h 后取**脾**的 CD45⁺ 绿与红 | 10x 5′ v3 | **否**——红＝自结肠迁出者，绿＝脾本地细胞，**不同室、极性相反** | 不用于 §2 主判 |
| `GSE195817` | 炎症关节／耳（关节炎） | 否 | 光转换后分选「已迁移／未迁移」白细胞 | **bulk RNA-seq** | **否**——**bulk 无法计算逐细胞检出率**，S-4 的整条流程不适用 | **技术上不可用**，如实记录 |
| `GSE221064` | MC38 肿瘤，NK | **是** | 本项目既有记录：仅 filtered matrix，判据 (a) 不通过 | 10x | —— | 仍 **`NOT_RETRIEVED`**（`CLAIMS_S3.md` A1） |
| `GSE203557`、`GSE193654`、`GSE221513`、`GSE298881`、`GSE322572`、`GSE296291` | 各类小鼠肿瘤 | **是** | 光转换肿瘤 | 各异 | —— | **均为肿瘤，不能作非肿瘤对照** |

---

## 3. 概念上最贴切、但**无数据**的一项

**Withers 实验室 2019 年 *Science Immunology* 的 `PMID 31152090`**：*Peripheral lymph nodes contain migratory and resident innate lymphoid cell populations*。

- 设计：**光转换整个臂丛淋巴结**，追踪新迁入的 ILC 与经 S1P 受体依赖迁出的 ILC。
- 关键契合点：**「大多数迁移性 ILC 是 ILC1，直接自循环进入，行为如同常规 NK 细胞」**——这是**非肿瘤组织中、最接近 NK 谱系、同型（同一组织内新到 vs 滞留）的对照**。
- **数据可得性：`NOT_AVAILABLE`。** GEO 以 `31152090[PMID]`、`brachial lymph node Kaede`、`innate lymphoid cell photoconver* lymph node` 检索均为 **0 命中**。该文非开放获取，其数据可得性声明未能核实，记 **`[待核实]`，不进入结论**。

**这一条须留在记录里**：它说明「非肿瘤组织中 NK／ILC1 的同型 Kaede 对照」在**实验层面已经做过**，只是**转录组数据未公开或不存在**。若日后需要最贴切的对照，向该作者群索取数据是成本最低的路径。

---

## 4. 主选数据集 `GSE228629` 的取数与验证

| 项 | 值 |
|---|---|
| 来源 | `ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228629/suppl/GSE228629_RAW.tar`（180 MB） |
| 论文 | *Skin-derived myeloid precursors and joint-resident fibroblasts spread psoriatic disease from skin to joints*，**Nat Immunol 2026**，PMID 41482544 / PMC12764428（开放获取，全文已读） |
| 文件 | 两个文库：`GSM7134137_B6`（C57BL/6）与 `GSM7134139_C`（BALB/c），各含 barcodes／features／matrix |
| 矩阵 | 32,291 行 = **32,285 基因 + 6 条哈希标签**（`Antibody Capture`）；B6 14,923 个细胞，BALB/c 13,072 个 |
| 哈希标签 | `skin_CD45neg_GREEN`、**`skin_CD45pos_RED`**、**`skin_CD45pos_GREEN`**、`joint_CD45neg_GREEN`、`joint_CD45pos_RED`、`joint_CD45pos_GREEN` |
| **计数类型** | **CellRanger 原始整数计数**——二项稀释零模型可直接使用，无需还原 |

### 4.1 解复用规则（声明，非调出）

细胞的哈希标签总计数 ≥ 20，且最高标签 ≥ 3 × 次高标签。B6 命中 **9,866 / 14,923（66.1%）**，BALB/c **7,967 / 13,072（60.9%）**。

**皮肤 CD45⁺ 的可用细胞数**（§2 主比较）：

| 文库 | 绿（<24 h 新到） | 红（滞留） |
|---|---|---|
| C57BL/6 | **1,068** | **3,062** |
| BALB/c | **1,199** | **2,012** |
| 合计 | **2,267** | **5,074** |

作为对照，`E-MTAB-10176` 的 NK 为绿 385 / 红 650。

---

## 5. 未做

- 未使用 bulk RNA-seq 数据集（`GSE195817`）——逐细胞检出率不可计算
- 未把「不同室、极性相反」的数据集（`GSE216928`、`GSE241820`）用于 §2 主判
- 未向作者索取 `PMID 31152090` 的数据
- 未下沉到 SRA／ENA 原始 FASTQ 层
- 未因主选已找到而停止记录其余候选——全部列于 §2
