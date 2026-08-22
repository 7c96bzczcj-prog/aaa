# 分选 NK bulk 数据集检索结果 —— 任务 P-2，§1

## §0.1 必须置于开头的限制

> 离子通道为催化型元件，单个通道每次开放可通过约 10⁶ 个离子，故其功能重要性与表达丰度
> 不成比例。更重要的是，`PIEZO1` 的活性取决于质膜张力——同样的表达量，在悬浮球形细胞与
> 被基质挤压的细胞上，通道活性可相差数个数量级，而转录组对此完全不可见。
>
> **因此本任务的产出只能支持「值得进一步测定」，不能支持「该通道重要」。**

---

## 一、四组对比的检索结果（一句话版）

| 序 | 对比 | 预期可得性 | **实际结果** | 采用的数据集 |
|---|---|---|---|---|
| **一** | CD56^bright 对 CD56^dim | 较高 | **找到 2 个** | `GSE236394`（6 供者，配对）、`GSE133383`（4 供者 × 5 组织，配对） |
| **二** | 静息对细胞因子活化 | 较高 | **找到 2 个（含小鼠臂共 3 臂）** | `GSE140035` 人（6 供者 × 4 条件，配对）与小鼠（3 重复 × 8 条件）、`GSE242941` bulk 臂（3 vs 3，弱） |
| **三** | 外周血对组织驻留 | 较低 | **找到 2 个** | `GSE133383`（血 vs 骨髓／脾／肺／肺门淋巴结，供者内）、`GSE200319`（血 vs 肝灌注液，非配对） |
| **四** | 肿瘤浸润对外周血 | 最低 | **找到 1 个，但它通不过自己的状态对照** | `GSE205492`（4 例软组织肉瘤，同一患者配对） |

**第一、二组均有可用数据**，故按交接语第 2 条继续执行，未停下。

第四组须特别说明：`GSE205492` 是唯一找到的「分选 NK 的配对肿瘤／外周血 bulk」，
但其 9 个状态对照基因**无一在两臂间显著不同**（`IFNG` −0.12、`CD69` −0.19、`ITGAE` −0.14，
均 q > 0.9，n = 4）。按 §3.2 的精神，**这一组的两臂未被证明处在不同状态**，
因此该组的任何阴性结果都不能读作「肿瘤中 NK 的 PIEZO1 不变」，只能读作「这一组没有信息」。

---

## 二、采用的六个数据集（§3.2 要求的元信息）

| 数据集 | 分选方式 | 材料 | 样本 / 供者 | 平台与建库 | 配对 | 测序深度 | 单位 |
|---|---|---|---|---|---|---|---|
| **GSE236394** | FACS：先按 CD56/CD16 分 bright／dim，再按 CD8a±；dump gate 排除 CD3/CD19/CD14 | 外周血 NK，**分选前在 1 ng/mL IL-15 中过夜静置** | 23 / 6 | GPL24676 NovaSeq 6000 | 是（供者内） | 未提交计数，深度不可得 `[待核实]` | 已提交的 moderated log2 CPM |
| **GSE133383** | FACS：CD56^bright CD16⁻ 与 CD56^dim CD16⁺ | 血、骨髓、脾、肺、肺门淋巴结，4 位器官捐献者 | 30 / 4 | GPL16791 HiSeq 2500 | 是（供者内，各供者组织覆盖不全） | 中位 21,098,240（9.9 M–40.7 M）；提交方已过滤至 11,943 基因 | log2(CPM+1)（由提交的计数算得） |
| **GSE140035**（人） | 「ex vivo 人 NK」，GEO 记录未写分选标记 `[待核实]` | 健康供者 NK，UNSTIM／IFNA／IL12+IL18／IL2+IL15 | 24 / 6 | GPL16791 HiSeq 2500 | 是（供者 4–9 在每个条件都在） | 中位 24,327,660（20.5 M–26.8 M） | log2(CPM+1) |
| **GSE140035**（小鼠） | 同上 `[待核实]` | WT 小鼠 NK，同一细胞因子面板 8 种组合 | 24 / 每条件 3 重复 | GPL21103 HiSeq 4000 | 否 | 中位 26,053,348（22.1 M–31.1 M） | log2(CPM+1) |
| **GSE242941**（bulk 臂） | 「从 PBMC 分选的 NK」，门控未写 `[待核实]` | 血 NK，未处理 vs 细胞因子处理；三个重复是否为三位供者未写 `[待核实]` | 5 / 3 vs 3（**其一被剔除，见下**） | GPL24676 NovaSeq 6000 | 否 | 中位约 1.2 M，**比本表其余数据集浅一个数量级** | log2(CPM+1) |
| **GSE200319** | FACS：CD56^hi CD16⁻、CD56^hi CD16⁺、CD56^lo CD16⁺ | 肝灌注液（5 供者）vs 健康外周血（5 供者） | 23 / 5+5，**不同个体** | GPL18573 NextSeq 500，改良 SMART-Seq2 低起始量 | 否 | 未提交 `[待核实]` | 已提交的 log2 归一化值（**非 CPM**） |
| **GSE205492** | FACS：CD3⁻CD56⁺ NK，同一患者的血与肿瘤 | 软组织肉瘤（脂肪肉瘤、黏液纤维肉瘤、UPS） | 8 / 4 | GPL20301 HiSeq 4000 | 是（同一患者） | 未提交（仅归一化计数）`[待核实]` | log2(归一化计数+1)（**非 CPM**） |

**分选方式的差异必须被记住**（§3.2）：`GSE205492` 用 CD3⁻CD56⁺，`GSE236394` 用 CD56/CD16
加 CD3/CD19/CD14 dump gate，`GSE133383` 与 `GSE200319` 用 CD56/CD16 双参数门。
这些群体不完全相同，**这正是 §3.1 禁止跨数据集比较绝对值的理由之一**。

### 一个自补门槛，如实报告

`GSE242941` 的 `control_1` 文库总计数 **32,276**，而同批其余文库为 1.1 M–2.2 M；
它与其余样本的相关系数为 0.14–0.40，而其余样本彼此为 0.63–0.93。**这是一条失败文库。**
执行处自补了一条 QC：**总计数 < 1×10⁵ 的样本剔除**（该阈值比次低文库低一个数量级，
只剔除这一个样本）。

门槛变动的影响：保留该样本时，`PIEZO1` 的区间为 **[−5.64, +8.36]**（纯粹由该文库的
噪声撑开）；剔除后为 **[−0.88, +0.30]**。两者的**判读相同（无显著差异）**，但前者的
区间宽到不携带任何信息。该数据集在剔除后仅剩 2 vs 3，是本任务中最弱的一个。

---

## 三、§1 点名的候选数据源，逐个核实

任务要求「逐个核实，不得预设可用」。

| 数据源 | 核实结果 | 能否用于四组对比 |
|---|---|---|
| **DICE** | **已核实。** 其 bulk RNA-seq 细胞类型清单中 NK 只有**一个**群体：`NK cell, CD56dim CD16+`。无 bright／dim 拆分 | **否**——只有单一 NK 群体，四组对比一组也构不成 |
| **Monaco 等 2019**（`GSE107011`） | **已核实。** 29 种免疫细胞、4 位供者，NK 样本共 4 个（`925L_NK`、`9JD4_NK`、`DZQV_NK`、`G4YW_NK`），**未拆分亚群、无活化臂、无组织臂** | **否**——同上 |
| **BLUEPRINT** | **未能核实，记 `[待核实]`。** 门户存在（BLUEPRINT Data Analysis Portal），但（a）未能确认其 RNA-seq 是否含 bright／dim 拆分的 NK；（b）BLUEPRINT 的个体层数据多经 EGA 受控访问，本会话无法开放下载 | 未定 |
| **ImmGen**（`GSE109125`，ULI RNA-seq） | **已核实为存在。** 含小鼠 NK 成熟亚群 `NK.27+11b−`、`NK.27+11b+`、`NK.27−11b+`，来源骨髓与脾，**每组 n = 2** | **可得但未采用**——每组 n = 2 低于任何可作推断的下限。记 `AVAILABLE_NOT_USED`，理由为样本量，而非未检索 |
| **Human Protein Atlas** | **已核实**（任务 P-1 即已取得）。分选血液 NK 的 `PIEZO1` = 1.0 nTPM，血液细胞特异性 "Low immune cell specificity" | **否**——只有单一 NK 群体，无对比轴；只能作绝对值参照，不能作差异检验 |

---

## 四、检索到但未采用的数据集，及理由

「没找到」记 `NOT_FOUND`，永不留空；此处记录的是**找到了但不适用**的，理由逐条写明。

| 数据集 | 为何不采用 |
|---|---|
| `GSE122324` | CD94⁻CD56dim／CD94⁺CD56dim／CD94⁺CD56hi 三个群体，但 **n = 3（每群一个样本）**，无重复，不能作检验 |
| `GSE300685` | 题名为「肿瘤浸润 NK 分型」，但 8 个样本的元数据只有 `tissue: Lung`，**无肿瘤／癌旁标注、无分选 NK 标注**，实为组织 bulk |
| `GSE227664` | **NK-92 细胞系** ＋ HDAC 抑制剂，不是原代 NK 的细胞因子活化 |
| `GSE255486` | 卵巢癌 PD-1⁺／PD-1⁻ NK，但肿瘤侧只有 **2 位患者**，且对比轴是 PD-1 而非肿瘤／外周血 |
| `GSE242941` scRNA 臂（24 个样本） | 单细胞，P-1 已证明 `PIEZO1` 在液滴数据的 NK 中落地板带；只取其 bulk 臂 |
| `GSE159624`、`GSE130430`、`GSE246994`、`GSE212890` 等 | 单细胞数据集，同上 |
| `GSE165461`、`GSE168212`、`GSE116178` | 题名相关，但未逐一取回样本级元数据核实分选与对比轴，记 `NOT_VERIFIED`——它们不进入任何结论 |

---

## 五、检索式（全部于 2026-08 执行，NCBI E-utilities `db=gds`）

1. `natural killer[Title] AND (CD56bright[All Fields] OR CD56dim[All Fields])` → 16 条
2. `CD56bright AND CD56dim AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter]` → 15 条
3. `natural killer[Title] AND (IL-15 OR IL-2 OR IL-12 OR IL-18 OR cytokine[Title] OR activation[Title] OR stimulated) AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter]` → 51 条
4. `NK cell[Title] AND (liver OR decidua OR uterine OR tissue-resident) AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter]` → 21 条
5. `natural killer[Title] AND (tumor-infiltrating OR intratumoral) AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter]` → 9 条
6. `(NK cells[Title] OR natural killer[Title]) AND (sorted OR FACS OR purified) AND (tumor OR ascites OR carcinoma) AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter] NOT single-cell` → 35 条

另：Ensembl Compara REST（直系同源逐个核实，`results/p2_orthologs.csv`）；
DICE 网站细胞类型清单；BLUEPRINT 门户（未取得数据）。

**未检索的方向**（记 `NOT_SEARCHED`，非「没找到」）：
- 蜕膜 NK 的分选 bulk（第三组的一个子方向）——检索式 4 命中了几个蜕膜相关系列
  （`GSE330176`、`GSE319493`、`GSE330178`），但均未取回样本级元数据核实是否为分选 bulk；
- ArrayExpress／ENA 独立检索（只走了 GEO）；
- 非英文文献与未提交 GEO 的数据。

---

## 六、复现

```bash
python3 scripts/p2_fetch.py        # 下载六个数据集的表达矩阵（约 15 MB）
python3 scripts/p2_orthologs.py    # 人／小鼠直系同源逐个核实
python3 scripts/p2_bulk_analysis.py
python3 scripts/verify_claims_p2.py
```
