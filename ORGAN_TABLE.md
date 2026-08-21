# ORGAN_TABLE.md —— 跨器官 NK 趋化因子检出率

> ## ⚠ 本文件已被 `ORGAN_TABLE_v2.md` 取代（S-1a）
> 表中的**数字仍然有效**（全 NK 检出率与标签无关，v2 逐位复现），但**基于本表的 §6 判读已改**：
> - **§6 第二行（XCL1 器官差异主要为亚群构成）已降级为「不可判定」**，归因规范设计缺陷（`CLAIMS_ORGAN.md` H1）
> - **`CCL3`／`CCL4`／`CCL5` 的跨器官离散度经自助法判为不成立**（区间下界 < 1.2×，H7）
> - 本表未含 bright 细胞数，读者无法分辨「限 bright 的高检出率」是生物学还是 9 个细胞（H17）
>
> **请以 `ORGAN_TABLE_v2.md` 与 `RESIDENCY_FEASIBILITY.md` 为准。**

**唯一数据源**：Domínguez Conde 等，*Science* 2022，Tissue Immune Cell Atlas 的 `T & innate lymphoid cells` 子集（CELLxGENE，216,611 细胞，12 供者，17 组织，`disease == normal`）。**同一供者跨器官取样**，供者变异由设计控制。

**未过滤 droplet 矩阵**：不可得（CELLxGENE 发布的是 cell-called 矩阵）。按本规范 §1，本轮可用，但**检出率含未校正的环境 RNA 成分**。

**统计单元 = 供者**（§5）：每供者一个检出率，表中为供者层面**中位数**。细胞未跨供者合并。

**§4 硬性**：不同化学版本**不合并、不并列排序**。每行是 `组织 × 化学版本`。

**§0 硬性**：本表不含任何跨器官比较、检验或倍数变化。

---

## 通过 §7 的层

| 器官 | 平台/化学 | NK 界定 | 供者数 | NK 细胞数 | 中位 UMI | CD56^bright 占比 | 解离应激 | CCL3 | CCL4 | CCL5 | XCL1 | XCL2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bone marrow | 10x 5' v1 | author | 5 | 601 | 4302 | 0.113 | 1.21 | **0.817** | **0.902** | **0.838** | **0.316** | **0.377** |
| liver | 10x 5' v1 | author | 5 | 1515 | 4436 | 0.594 | 2.15 | **0.781** | **0.875** | **0.790** | **0.723** | **0.507** |
| lung | 10x 5' v1 | author | 6 | 324 | 3282 | 0.119 | 1.74 | **0.667** | **0.882** | **0.847** | **0.556** | **0.446** |
| spleen | 10x 5' v1 | author | 6 | 1361 | 5133 | 0.558 | 2.38 | **0.726** | **0.884** | **0.800** | **0.692** | **0.678** |
| thoracic lymph node | 10x 5' v1 | author | 6 | 257 | 4260 | 0.982 | 2.59 | **0.661** | **0.764** | **0.773** | **0.879** | **0.773** |
| spleen | 10x 5' v2 | author | 4 | 1751 | 4040 | 0.520 | 2.07 | **0.666** | **0.781** | **0.790** | **0.424** | **0.473** |
| spleen | 10x 5' v2 | uniform | 4 | 410 | 3976 | 0.493 | 1.95 | **0.727** | **0.813** | **0.764** | **0.515** | **0.537** |

### 同一层的 CD56^bright 限定值（§5 第 2 项）

| 器官 | 化学 | NK 界定 | CD56^bright 占比 | CCL3 (bright) | CCL4 (bright) | CCL5 (bright) | XCL1 (bright) | XCL2 (bright) |
|---|---|---|---|---|---|---|---|---|
| bone marrow | 10x 5' v1 | author | 0.113 | 0.714 | 0.905 | 0.889 | 0.800 | 0.500 |
| liver | 10x 5' v1 | author | 0.594 | 0.866 | 0.889 | 0.753 | 0.901 | 0.584 |
| lung | 10x 5' v1 | author | 0.119 | 0.600 | 0.889 | 0.778 | 1.000 | 0.750 |
| spleen | 10x 5' v1 | author | 0.558 | 0.760 | 0.873 | 0.762 | 0.854 | 0.805 |
| thoracic lymph node | 10x 5' v1 | author | 0.982 | 0.653 | 0.762 | 0.769 | 0.912 | 0.783 |
| spleen | 10x 5' v2 | author | 0.520 | 0.706 | 0.738 | 0.741 | 0.773 | 0.785 |
| spleen | 10x 5' v2 | uniform | 0.493 | 0.724 | 0.763 | 0.723 | 0.750 | 0.712 |

---

## 判为 `INSUFFICIENT` 的层（§7，如实列出，不外推）

| 器官 | 化学 | NK 界定 | 供者数 | NK 细胞数 | 卡在哪一条 |
|---|---|---|---|---|---|
| blood | 10x 3' v3 | author | 2 | 3511 | donors 2<3; one donor 73% of cells |
| blood | 10x 5' v1 | author | 2 | 206 | donors 2<3; one donor 61% of cells |
| blood | 10x 5' v2 | author | 2 | 267 | donors 2<3; one donor 99% of cells |
| bone marrow | 10x 3' v3 | author | 2 | 6872 | donors 2<3; one donor 73% of cells |
| bone marrow | 10x 5' v2 | author | 3 | 753 | one donor 51% of cells |
| caecum | 10x 5' v1 | author | 1 | 2 | donors 1<3; NK 2<200; one donor 100% of cells |
| duodenum | 10x 5' v1 | author | 2 | 5 | donors 2<3; NK 5<200; one donor 80% of cells |
| ileum | 10x 5' v1 | author | 1 | 1 | donors 1<3; NK 1<200; one donor 100% of cells |
| jejunal epithelium | 10x 3' v3 | author | 2 | 60 | donors 2<3; NK 60<200; one donor 93% of cells |
| jejunal epithelium | 10x 5' v2 | author | 3 | 66 | NK 66<200; one donor 53% of cells |
| jejunum lamina propria | 10x 3' v3 | author | 2 | 72 | donors 2<3; NK 72<200; one donor 69% of cells |
| jejunum lamina propria | 10x 5' v2 | author | 3 | 71 | NK 71<200; one donor 77% of cells |
| liver | 10x 5' v2 | author | 2 | 2072 | donors 2<3; one donor 54% of cells |
| lung | 10x 3' v3 | author | 2 | 2880 | donors 2<3; one donor 72% of cells |
| lung | 10x 5' v2 | author | 2 | 29 | donors 2<3; NK 29<200; one donor 66% of cells |
| mesenteric lymph node | 10x 5' v1 | author | 5 | 41 | NK 41<200 |
| mesenteric lymph node | 10x 5' v2 | author | 3 | 950 | one donor 84% of cells |
| omentum | 10x 5' v1 | author | 4 | 33 | NK 33<200; one donor 55% of cells |
| sigmoid colon | 10x 5' v1 | author | 1 | 1 | donors 1<3; NK 1<200; one donor 100% of cells |
| skeletal muscle tissue | 10x 5' v1 | author | 3 | 134 | NK 134<200; one donor 95% of cells |
| spleen | 10x 3' v3 | author | 2 | 5098 | donors 2<3; one donor 71% of cells |
| thoracic lymph node | 10x 3' v3 | author | 2 | 377 | donors 2<3; one donor 90% of cells |
| thoracic lymph node | 10x 5' v2 | author | 4 | 177 | NK 177<200 |
| thymus | 10x 5' v1 | author | 1 | 2 | donors 1<3; NK 2<200; one donor 100% of cells |
| transverse colon | 10x 5' v1 | author | 1 | 4 | donors 1<3; NK 4<200; one donor 100% of cells |
| blood | 10x 3' v3 | uniform | 2 | 330 | donors 2<3; one donor 75% of cells |
| blood | 10x 5' v1 | uniform | 2 | 48 | donors 2<3; NK 48<200; one donor 71% of cells |
| blood | 10x 5' v2 | uniform | 1 | 57 | donors 1<3; NK 57<200; one donor 100% of cells |
| bone marrow | 10x 3' v3 | uniform | 2 | 1318 | donors 2<3; one donor 88% of cells |
| bone marrow | 10x 5' v1 | uniform | 4 | 99 | NK 99<200 |
| bone marrow | 10x 5' v2 | uniform | 3 | 126 | NK 126<200; one donor 54% of cells |
| duodenum | 10x 5' v1 | uniform | 1 | 1 | donors 1<3; NK 1<200; one donor 100% of cells |
| ileum | 10x 5' v1 | uniform | 1 | 1 | donors 1<3; NK 1<200; one donor 100% of cells |
| jejunal epithelium | 10x 3' v3 | uniform | 2 | 12 | donors 2<3; NK 12<200; one donor 83% of cells |
| jejunal epithelium | 10x 5' v2 | uniform | 2 | 19 | donors 2<3; NK 19<200; one donor 58% of cells |
| jejunum lamina propria | 10x 3' v3 | uniform | 2 | 17 | donors 2<3; NK 17<200; one donor 88% of cells |
| jejunum lamina propria | 10x 5' v2 | uniform | 3 | 15 | NK 15<200; one donor 80% of cells |
| liver | 10x 5' v1 | uniform | 5 | 186 | NK 186<200; one donor 55% of cells |
| liver | 10x 5' v2 | uniform | 2 | 457 | donors 2<3; one donor 56% of cells |
| lung | 10x 3' v3 | uniform | 2 | 486 | donors 2<3; one donor 80% of cells |
| lung | 10x 5' v1 | uniform | 6 | 45 | NK 45<200 |
| lung | 10x 5' v2 | uniform | 2 | 5 | donors 2<3; NK 5<200; one donor 60% of cells |
| mesenteric lymph node | 10x 5' v1 | uniform | 3 | 14 | NK 14<200; one donor 71% of cells |
| mesenteric lymph node | 10x 5' v2 | uniform | 3 | 241 | one donor 82% of cells |
| omentum | 10x 5' v1 | uniform | 2 | 9 | donors 2<3; NK 9<200; one donor 67% of cells |
| skeletal muscle tissue | 10x 5' v1 | uniform | 2 | 31 | donors 2<3; NK 31<200; one donor 97% of cells |
| spleen | 10x 3' v3 | uniform | 2 | 1134 | donors 2<3; one donor 82% of cells |
| spleen | 10x 5' v1 | uniform | 6 | 268 | one donor 53% of cells |
| thoracic lymph node | 10x 3' v3 | uniform | 2 | 89 | donors 2<3; NK 89<200; one donor 92% of cells |
| thoracic lymph node | 10x 5' v1 | uniform | 4 | 41 | NK 41<200 |
| thoracic lymph node | 10x 5' v2 | uniform | 4 | 39 | NK 39<200 |
| transverse colon | 10x 5' v1 | uniform | 1 | 1 | donors 1<3; NK 1<200; one donor 100% of cells |

共 59 层，通过 §7 者 **7** 层。

---

## §6 判读（只能落在四行之一）

**可评估的层只有一个**：`10x 5' v1` × author 界定，含 **5 个器官**（bone marrow / liver / lung / spleen / thoracic lymph node），供者数 5–6。其余层因器官数 < 3 无法评估 §6（§4 禁止跨化学版本合并）。

| 基因 | 全 NK 跨器官离散（max/min） | 限定 CD56^bright 后 | ρ(深度) | ρ(bright 占比) | **判读** |
|---|---|---|---|---|---|
| CCL3 | 1.24× | 1.44 | +0.50 | −0.60 | **不可归因**（第四行） |
| CCL4 | 1.18× | 1.19 | +0.30 | −0.90 | **不可归因**（第四行） |
| CCL5 | 1.10× | 1.18 | −0.30 | −0.90 | **不可归因**（第四行） |
| **XCL1** | **2.78×** | **1.25** | +0.10 | **+1.00** | **主要为亚群构成**（第二行） |
| **XCL2** | **2.05×** | **1.61** | +0.20 | **+0.90** | **主要为亚群构成**（第二行） |

### XCL1／XCL2 —— 第二行：器官差异**主要为亚群构成**

XCL1 的跨器官离散 **2.78×**，限定 CD56^bright 后塌到 **1.25×**；与 CD56^bright 占比的秩相关为 **+1.00**（五个器官完全同序）。XCL2 同向（2.05× → 1.61×，ρ = +0.90）。**判据要求的"限定 CD56^bright 后差异大幅缩小"成立。**

### CCL3／CCL4／CCL5 —— 第四行：**不可归因**

三者的跨器官离散本就很小（1.10–1.24×）：不随深度单调（ρ = −0.30 至 +0.50，方向不一致），限定 CD56^bright 后也不缩小（反而略增），与 bright 占比呈**负**相关（−0.60 至 −0.90，即 bright 比例低的器官反而更高，与 XCL1 相反）。三条判据互不一致 → 记 **不可归因**，不作选择性解读。

---

## 必须随判读一同陈述的四条限制

1. **蜕膜不在本图谱内。** 17 个组织无生殖道样本——与本项目既有记录一致。**因此"蜕膜与肺之间 XCL1 相差 8 倍"这一具体问题，本轮并未回答。** 本轮回答的是：在单一化学版本下的五个器官之间，XCL1 的器官离散由 CD56^bright 构成支配。
2. **血液在每一种化学版本下都判 `INSUFFICIENT`**（各仅 2 位供者）。故 §4 特别要求的**组织 vs 血液解离应激对照无法做**。已通过 §7 的五个器官**全部经解离**，其应激分数 1.21–2.59，彼此无未解离的参照。
3. **统一标记集几乎不可用。** `KLRD1>0 & NCAM1>0 & CD3E==0 & CD3D==0` 只召回 **5,499 / 29,493 = 18.6%** 的作者注释 NK（`NCAM1` 在 10x 上大量 dropout）。**§3 要求的成对呈现只在 spleen `10x 5' v2` 一层成立**：该层两套界定的 XCL1 为 author 0.424 / uniform 0.515，方向一致，未出现方向性冲突。
4. **未过滤矩阵不可得**，检出率含未校正的环境 RNA 成分（§1 已声明）。

