# CLAIMS.md

一行一条断言。列：断言 / 证据类型 / 等级（M·V·L·X）/ 出处 / **是否推翻初判**（允许为"是"）。

范围：修正案 A-1 交接语第 2、3 项——Kersten 2022 极性来源的核实，与 A1.4 极性一致性检验。
未覆盖：§6 主计算（未跑）、G-KILL-2（未跑）。

等级：`[M]` 本项目实测 / `[V]` 原文已核实（Claude 本人打开） / `[L]` 模型知识未核实 / `[X]` 已被检验推翻。

---

## A. Kersten 2022 极性来源的核实（交接语第 2 项）

| # | 断言 | 证据类型 | 等级 | 出处 | 是否推翻初判 |
|---|---|---|---|---|---|
| K1 | Kersten 2022 的 CCL3/CCL4/CCL5 上调，其**火山图对比是"瘤内 Tex d14 **vs** 脾脏 naïve"**。图注原文："Volcano plot showing differential gene expression in tumor-infiltrating CD44+ OT-I CD8+ **Tex d14** cells (red) compared to **splenic CD44− OT-I CD8+ Tnaïve** cells (grey)" | primary，全文（Claude 本人两次定向打开 PMC 全文） | **[V]** | Kersten et al. 2022, *Cancer Cell* 40(6):624-638（PMC9197962），Fig 2A 图注 | **是** |
| K2 | **Fig 2C 的三组全部以 Tnaïve 归一**，且其中的 Teff 是**体外激活**产物、不取自肿瘤。图注："Average gene expression in Teff (black), Tex d4 (blue) and Tex d14 (red) **normalized to Tnaïve** as determined by RNA-seq"；Methods："OT-I T cells were **activated in vitro** … OT-I lymph node cells were stimulated with B6 splenocytes pulsed with SL8 peptide" | primary，全文 | **[V]** | 同上，Fig 2C 图注 + Methods | **是** |
| K3 | 该文**没有**任何"同一效应/记忆亚群在瘤内 vs 瘤外"的趋化因子对比 | primary，全文（定向检索，返回 NOT_STATED） | **[V]** | 同上 | **是** |
| K4 | **因此 Kersten 的极性不属于 A1.0 所列三个对比中的任何一个。** 它是第四种、且更弱的一种："瘤内耗竭 vs 脾脏**未致敏**"——与**分化阶段**（naïve→effector）完全混杂。CCL3/4/5 是经典效应/记忆基因，naïve CD8 T 近乎不表达，故**任何**部位的效应 CD8 T 对 naïve 都会显示 CCL3/4/5 上调，该对比从构造上无法分离肿瘤效应 | 推论，基于 K1–K3 的 [V] | **[V]**（前提）/ 推论 | 本报告 | **是** |
| K5 | 但该文**确实包含一个结构上可比的子层**：**Tex d4 → Tex d14**，两组均取自肿瘤，是瘤内**停留时长**轴，与本项目的绿/红同构。原文："most of these **increased with prolonged residence in the TME**"，且"the majority of these changes was **not** observed in effector CD8+ T cells (Teff)" | primary，全文 | **[V]** | 同上，Results | 否（部分支持 A1.1） |
| K6 | K5 的这一子层里，Ccl3/Ccl4/Ccl5 的方向是**随停留时间上升** | primary，全文 | **[V]** | 同上 | 否 |
| **K7** | **判定（按 A1.0 预写规则）**：Kersten 的极性论据**降级为"待定"**。A1.0 写明"若 Kersten 的极性只在'耗竭 vs 非耗竭'这一层成立，则 A1.1 的极性论据降级为待定"；实测结果比该条件**更弱**（连"耗竭 vs 非耗竭"都不是，是"耗竭 vs naïve"）。故**降级成立**。**A1.3–A1.6 照常执行**，其成立不依赖极性 | 推论 | **[V]** | 本报告 | **是**（对 A1.1 理由 1） |

---

## B. A1.4 极性一致性检验（交接语第 3 项）

数据集：`E-MTAB-10176`（Li 2022，Kaede/MC38 total TIL）。矩阵：`DW_T_NK.h5ad`，4,637 细胞。
Green = Kaede-Green = `Infiltrating` = 新进入；Red = Kaede-Red = `Resident` = 滞留。
单位：检出率（非零计数细胞占比），取自 `raw.X`（16,227 基因，log-normalized，零保留）。
**`adata.X` 已被 scale/center（min = −3.85，零占比 0.0000），检出率从构造上不可从它计算，全程未使用。**

| # | 断言 | 证据类型 | 等级 | 出处 | 是否推翻初判 |
|---|---|---|---|---|---|
| P1 | **`colour` 与文库完全混杂。** 4 个文库（MH1_G24 / MH2_R24 / MH3_G72 / MH4_R72），每个文库只含一种 colour；SDRF 显示 4 个 Source Name = 4 个 ENA sample = 4 个 BioSample，`Comment[technical replicate group]` 为 group 1–4 各一。obs 中**无任何** mouse/donor/animal/replicate 列 | 实测 + 归档元数据 | **[M]** | `scripts/a1_4_polarity_check.py`；`data/raw/E-MTAB-10176.sdrf.txt` | **是** |
| P2 | **绿/红每侧 n = 1 个文库。** 这正是 §4(d) 设立来抓的 library-as-donor 结构 | 实测 | **[M]** | 同上 | **是** |
| P3 | **A1.3c 的复孔噪声底在本数据集内无法测量**——不存在同条件复孔文库（每个 colour×hours 格子恰好 1 个文库）。A1.3c 明写"实测值不迁移"，故本项目 Dean 数据的 [M] 锚（XCL1 ±15.5、CCL4 ±11.3、CCL5 ±6.6 pp）**不能**充当本数据集的底 | 实测 | **[M]** | 同上 | 否（执行 A1.3c） |
| P4 | **NK 的模块比值在本数据集 UNCOMPUTABLE。** 杀伤模块五个基因全部落在 §5 禁区：Gzma / Gzmb / Nkg7 / Prf1 天花板（>85%），Gzmk 地板（<5%）。§5 后杀伤模块**零基因存活** | 实测 | **[M]** | `results/a1_4_detection_rates.csv` | **是** |
| P5 | 因此 **A1.5 的第一条预测（NK 比值下降）在本数据集内不可检验** | 推论，基于 P4 | **[M]** | 本报告 | **是** |
| P6 | **CD8T 的两个时间点方向相反。** 24h：Δ = **+0.0197**（Red>Green，上升）；72h：Δ = **−0.0512**（Red<Green，下降）。同谱系、同基因集、同数据集 | 实测 | **[M]** | `results/a1_4_module_ratio.csv` | **是** |
| P7 | **留一法（§6.2）在两个时间点都翻转。** 24h：去掉 Ccl4 / Flt3l / Gzma 任一即翻转（8 基因中 3 个）；72h：去掉 Ccl5 或 Prf1 即翻转。去掉 Prf1 使 72h 的 Δ 从 −0.0512 变为 **+1.3900** | 实测 | **[M]** | `results/a1_4_leave_one_out.csv` | 否（§6.2 按设计生效） |
| **P8** | **A1.4 判定 = 第三分支：`UNMEASURED`。** 且是**超定的**——P3（无噪声底）、P6（两时点反向）、P7（留一法翻转）、P1–P2（n=1 且与文库混杂）四条各自独立即可判 UNMEASURED。按 A1.4 预写规则，**"CD8T 下降"与"CD8T 上升"两条都不得引用** | 实测 | **[M]** | 本报告 | 否（按预写分支执行） |

---

## C. 检验过程中产生的、推翻规范自身初判的实测发现

| # | 断言 | 证据类型 | 等级 | 出处 | 是否推翻初判 |
|---|---|---|---|---|---|
| **N1** | **v1.0 §6 的"调节模块"在 NK 内部不是一个同向块。** 绿→红逐基因（pp）：24h Ccl3 **−34.4**、Ccl4 **−44.8**、Ccl5 **−62.0**，而 Xcl1 **+7.9**、Flt3l **+11.9**；72h Ccl3 −29.3、Ccl4 −44.8、Ccl5 −54.4，而 Xcl1 **+8.9**、Flt3l −3.3。**Xcl1 与 Ccl3/4/5 系统性反向。** 把它们平均成一个"模块均值"，是在把 −45pp 的块和 +9pp 的块对冲 | 实测 | **[M]** | `results/a1_4_detection_rates.csv` | **是** —— 直接质疑 §6 的模块构造本身 |
| **N2** | **A1.1 理由 2（基因集不可比）获得比原论证更强的实测支持。** 原论证是"XCL1 在 CD8 T 落地板"；实测是**两个谱系在 §5 后根本不共享基因集**：NK admitted = {Ccl3, Ccl4, Xcl1, Flt3l} + 杀伤模块**空集**；CD8T admitted = {Ccl3, Ccl4, Ccl5, Xcl1, Flt3l} + {Prf1, Gzma, Gzmk}。**"在 CD8T 上跑同一计算"在本数据集内字面意义上不存在** | 实测 | **[M]** | 同上 | 否 —— **加强** A1.1（且不依赖极性） |
| N3 | **NK 侧最大的那个降幅从天花板起跳，§5 明令不可作阴性证据。** NK Ccl5 Green = **1.0000**（24h）/ 0.8715（72h），均在天花板带；−62.0 pp 正是 §5"天花板压缩…不可作为阴性证据"要挡的数 | 实测 | **[M]** | 同上 | 否（§5 按设计生效） |
| N4 | **绿→红几乎完全是亚群构成的更换，不是同一群细胞换状态。** NK：24h Green 84.5% 是 NK-1，Red 仅 3.5% 是 NK-1、89.5% 是 NK-2。CD8T：72h Green 26.3/28.4/45.3%（CD8T-1/2/3），Red 0.5/5.2/**94.2%**。按 §7，between-composition 项将占主导，结论须改写为构成陈述 | 实测 | **[M]** | `results/a1_4_polarity_check.log` 构成诊断段 | 否 —— 预警 §7 |
| N5 | **`Xcl2` 不在该归档矩阵中**（16,227 基因里没有），本数据集内记 `ABSENT`，非 `NOT_FOUND` 亦非地板 | 实测 | **[M]** | 同上 | 否 |
| N6 | 描述性、且**不得引用**：CD8T 的 Ccl3/Ccl4 随停留上升（24h +3.9/+10.9 pp，72h +7.2/+9.7 pp），方向与 K5/K6 的 Kersten d4→d14 子层一致；CD8T 的 Ccl5 在 72h 下降 −23.0 pp，与 NK 的 Ccl5 同向。这**与 [M] 47% 趋同度定性相容**，但按 P8 全部记 `UNMEASURED` | 实测 | **[M]**（数值）；**推论不成立**（无噪声底） | 同上 | 否 —— 但**不得**进入结论 |

---

## D. 未做 / 未搜（§3.1，不留空）

| 项 | 状态 |
|---|---|
| §6 主计算（NK 与 CD8T 的模块比值正式对比） | `NOT_RUN` —— 交接语明令停下 |
| G-KILL-2（Q4 是否已被占位） | `NOT_RUN` —— 交接语明令停下 |
| Kitagawa 三项分解（§7） | `NOT_RUN` —— 仅出具了 N4 的构成诊断，未做分解 |
| Kersten 2022 补充材料 / 逐基因 DE 表 | `NOT_SEARCHED` —— 仅主文 Results/Methods/图注 |
| 从 `PRJNA912695` / E-MTAB-10176 FASTQ 重跑 CellRanger 以取得未过滤 droplet 矩阵 | `NOT_ATTEMPTED` —— 见 `DATASETS.md`；即使重跑也不能造出不存在的生物学重复 |
| Dean 2024 的 `GSE221064` 上做同一 A1.4 检验 | `NOT_RUN` —— 仅 filtered matrix，且同属 n 问题，待确认后再定 |
