# CIRCULARITY_CHECK.md —— S-1 §2

**问题**：Tissue Immune Cell Atlas 在界定 `NK_CD56bright_CD16-` 时，是否使用了 `XCL1` 或 `XCL2`？

**结论：**

> ## 命中。`XCL1` 实质参与了 bright／dim 的界定。
> 按 §2 预写判读，**上一轮的 ρ = +1.00 不得作为结论使用，降级为「部分由定义产生」。**

---

## 四条查证路径，逐条记录

### 路径 1 —— 注释模型的特征权重（**命中**）

`Domínguez Conde 等 2022` **确为 CellTypist 原始论文**，原文：「we developed **CellTypist**, a machine learning tool for rapid and precise cell type annotation」`[V]`。模型为逻辑回归，权重可直接读取。

下载 `Immune_All_Low.pkl` 并读取系数矩阵（6,639 个特征）：

| 类别 | XCL1 | XCL2 | NCAM1 | FCGR3A |
|---|---|---|---|---|
| **CD16- NK cells**（bright 侧） | **+0.1557，排名 15 / 6,639** | +0.0745，排名 424 | +0.1135，排名 96 | −0.1928，排名 6,625 |
| **CD16+ NK cells**（dim 侧） | **−0.2499，排名 6,635 / 6,639** | −0.0306，排名 5,667 | −0.1067，排名 6,576 | +0.2322，排名 8 |

**`XCL1` 在 bright 侧位列前 15、在 dim 侧位列倒数第 5。** 它不是边缘特征，而是该分类器用来分开两个 NK 亚群的主力特征之一。

**附带发现（另一个目标基因）**：`CCL3` 在 `CD16- NK cells` 的权重为 **+0.292，排名 3 / 6,639**——即本轮五个目标基因中，**有两个**出现在 bright 类的高权重特征里。

`[M]`，`scripts/` 外的一次性读取，权重值见本文件。

### 路径 2 —— 论文正文与补充材料（**未命中**）

取 Europe PMC 全文 XML（`PMC7612735`，CC BY，72,065 字符）：

- **`XCL1` 出现 0 次，`XCL2` 出现 0 次。**
- NK 亚群的界定表述为：「Lastly, **NK cells in our data were represented by two clusters with high expression of either `FCGR3A` (encoding CD16) or `NCAM1` (encoding CD56)**.」`[V]`
- 正文中 `CD56bright` 字样出现 0 次（该文用 CD16±，不用 bright／dim 措辞）。

**故就论文的书面理据而言，NK 亚群是按 `FCGR3A` / `NCAM1` 划的，与 XCL1／XCL2 无关。**

### 路径 3 —— h5ad 的 `.uns`（**不可得**）

`.uns` 仅含 CELLxGENE schema 字段（`citation`、`schema_version`、颜色映射等），**无 `rank_genes_groups` 或任何标记基因信息**——已被发布流程剥除。该路径记 `NOT_DETERMINED`。

### 路径 4 —— 人工注释所用标记（**部分可得，且与路径 1 相连**）

论文原文：manual curation 是**建立在 CellTypist 注释之上**的（"we built on CellTypist's annotation by performing additional manual curation"）。因此关键在于**两套标签的重合度**——若人工标签基本沿用 CellTypist 标签，路径 1 的循环论证即传导至我实际使用的 `Manually_curated_celltype`。

实测（本地缓存，29,493 个作者注释 NK 细胞）：

| `Manually_curated_celltype` | → CellTypist `CD16+ NK cells` | → `CD16- NK cells` | → `NK cells` | → `Cytotoxic T cells` |
|---|---|---|---|---|
| `NK_CD16+`（n = 20,587） | **95.1%** | 0.0% | 0.1% | 4.8% |
| `NK_CD56bright_CD16-`（n = 8,843） | 1.2% | **78.0%** | 19.9% | 0.8% |

**二分一致率（人工 bright ↔ CellTypist CD16−）= 0.932。** `[M]`

---

## 判读（按 §2 预写表）

| 判据 | 命中情况 |
|---|---|
| XCL1／XCL2 参与了 bright 的界定 | **是** —— 路径 1 命中（XCL1 排名 15/6,639 正、6,635/6,639 负），并经路径 4 的 **93.2% 一致率**传导至本轮实际使用的标签 |

→ **上一轮 `XCL1` 与 CD56^bright 占比 ρ = +1.00 的结果，降级为「部分由定义产生」，不得作为结论使用。**

### 降级的精确范围（不要过度外推）

1. **受影响最重的是 `XCL1`**（排名 15）。`XCL2` 排名 424，受影响明显较轻，但同向，一并降级。
2. **`CCL3` 也在高权重之列（排名 3）**，但上一轮实测其与 bright 占比的相关是 **−0.60（负）**——**与循环论证会产生的方向相反**。故循环传导并非对所有基因一致，这一点是实测支持的，不是推测。
3. **论文的书面理据（`FCGR3A`／`NCAM1`）本身是干净的。** 循环来自注释**流程**（CellTypist 模型权重），不是来自论文所述的划分依据。二者须分开陈述。
4. **本条不推翻「XCL1 的器官差异与 bright 构成高度共变」这一观察本身**——它推翻的是把该共变**当作独立证据**的资格。真实的生物学共变与定义产生的共变，在本数据集内无法分离。

---

## 对上一轮结论的处置

| 上一轮陈述 | 现状态 |
|---|---|
| XCL1 落 §6 第二行（器官差异主要为亚群构成） | **降级** —— 判读方向不变，但**证据等级下降**：ρ = +1.00 部分由定义产生 |
| XCL1 跨器官离散 2.78× → 限定 bright 后 1.25× | **仍可用** —— 这是**离散度塌缩**，不是相关系数；但须注意 bright 子集本身由含 XCL1 的分类器划定 |
| CCL3／CCL4／CCL5 落第四行（不可归因） | **不受影响** —— 见上第 2 点，CCL3 的相关方向与循环论证相反 |
| 「未匹配 bright 占比的跨区室 XCL1 比较，测的主要是 bright 比例」 | **须改写** —— 现应为：「**该数据集内，XCL1 与 bright 构成不可分离**，原因兼含生物学与注释定义」 |

---

## 未做

| 项 | 状态 |
|---|---|
| 论文补充材料中的逐细胞类型标记基因表（table S / fig. S 的点图） | `NOT_SEARCHED` —— 全文 XML 不含补充表；正文已给出 NK 的划分依据，且路径 1 已足以判定 |
| 其他 CellTypist 模型（如 `Immune_All_High.pkl`） | `NOT_SEARCHED` —— 本图谱用的是 low（细粒度）层级 |
| 用不含 XCL1 的自建 bright 定义重算 | **未做，且不应在本轮做** —— 属 §7 禁止的分层操作，需先由人裁决 |
