# 跨器官（肺·肝·胰·脑·脊髓）单细胞 IL18 高 / CASP1 低细胞类型筛选

规格 v3.0。数据集清单 v0.1（15 个分层，经用户确认推进）。CELLxGENE Census LTS 2025-11-08。随机种子 20261013。

## 本分析的定位（规格 §0.2，逐字声明）

> 一、IL18 基因只有一个蛋白产物，即 24 kDa 的前体 pro-IL-18。成熟 IL-18 由蛋白酶切割前体产生，不存在对应的独立转录本。故本分析测得的 IL18 转录本丰度即 pro-IL-18 的转录信号。但 pro-IL-18 属储备型蛋白，半衰期长，转录水平与细胞内前体蛋白储量不成线性关系。
>
> 二、实际起筛选作用的是 IL18 与 CASP1 两个基因。CASP3 在绝大多数细胞类型中组成型表达，不具备细胞类型间的区分力，在本分析中仅用于排除最低档和给候选集内部排序，不参与候选集的定义。本分析不是三重阳性门。
>
> 三、CASP3 的 mRNA 水平不反映 caspase-3 的酶活性。caspase-3 以无活性酶原形式存在，需被起始型 caspase 切割活化，且活化后受 XIAP 等凋亡抑制蛋白直接抑制。这两层调控均在转录后，本数据看不见。因此本分析无法判断任何细胞类型中 caspase-3 是否真的作用于 pro-IL-18。
>
> 四、本分析为静息基线筛选。一个 CASP1 低的细胞类型在受到刺激后能否诱导出 CASP1，本数据无法回答。CASP1 为可诱导基因，此限制不能通过增加基线转录组基因来弥补，只能通过带刺激条件的数据集或实验解决。

**本分析是底物充足且经典加工酶缺失的交集筛选，用途是为后续蛋白与功能层面的验证挑选器官和细胞类型的优先级。它不是任何机制结论的证据。**

## 单核 vs 单细胞不可混比（规格 §3，逐字）

> 由于脑与脊髓以 nucleus 为主、其余三器官以 cell 为主，本分析在严格意义上无法做五器官之间的直接定量比较。

跨器官只能比较「候选集成员名单」这一层级的定性结论，不能比较表达量或排序值。所有排序、分位数、候选判定均在「dataset × tissue × suspension_type」分层内部进行，不跨 suspension_type、不跨数据集、不跨物种比较表达量。

## 基因 ID（规格 §1.1，全部一致）

六个 Ensembl ID 经 Ensembl REST 解析并与 Census feature 表交叉核对，全部与规格参考表一致，无不一致项（详见 `results/gene_id_resolution.md`）。取基因一律用 Ensembl ID，不用符号匹配。唯一例外：小鼠脊髓 GSE103892 的矩阵只带基因符号、无 Ensembl ID，故标 `gene_id_unavailable=TRUE`，仅作参考、不进候选判定（§1.2）。

## 环境 RNA 校正（关键限制）

本机无 GPU（CellBender 不可行），Census 源矩阵为已过滤矩阵、无空液滴（SoupX 不可行）。经用户确认（决策 ①B），**候选集在未做环境 RNA 校正的数据上计算，全部标 `ambient_uncorrected=TRUE`，属临时/待验证结果**，不作为机制证据。环境污染是本分析首要假阳性来源（巨噬细胞游离 RNA 污染上皮/神经元簇等），相关风险在候选表与零值判读中标注。

## 覆盖与缺口

- 人：肺、肝、胰各有 cell + nucleus 两套；脑、脊髓为 nucleus（人成年正常脑在 Census 中无 cell 数据集，记为缺口）。
- 小鼠：肺、肝、胰为 cell（Census 无 nucleus）；脑有 cell + nucleus；脊髓 Census 无数据，改用 GEO 的 GSE103892（Sathyamurthy 2018，单核腰段脊髓，naive 基线），不以脑/背根神经节顶替。
- 人脊髓 nucleus 仅 1 供体（Siletti），按 §10.3 只作描述统计，不产候选。

## 列定义（检出率的百分比到底是什么的百分比）

| 列 | 定义 |
|---|---|
| `pct_IL18` / `pct_CASP1` / `pct_CASP3` | **检出率（供体中位数口径）**。分母 = 某分层（数据集 × 器官 × suspension）内、**某一位供体、某一个细胞类型所含的细胞个数**；分子 = 这些细胞中**该基因原始计数 ≥ 1 的细胞个数**。汇总表再对该分层内各供体取**中位数**。即「有多少比例的细胞里测得到这个基因」——**不是**表达量占比，**不是**相对最高值的百分比，**不是**读段占比。 |
| `pooled_pct_IL18` / `pooled_pct_CASP1` | **检出率（合并口径，对外展示一律用这个）**。分母 = 该细胞类型在该分层内、**每位供体 ≥20 细胞后进入分析的全部细胞**；分子 = 其中原始计数 ≥1 的细胞数。与上一行同一个意思，只是不先按供体取中位数。<br>**为什么必须并列给出**：每位供体只有数十到上百个该类型细胞，当真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，**供体中位数被截断为 0**。例：人肺肥大细胞 CASP1 合并 138/5299 = 2.60%（21/44 供体有阳性），小鼠脑少突胶质 CASP1 合并 2/8987 = 0.02%，两者的供体中位数都写作 0.00%，真实值却差一百倍以上。 |
| `n_cells_IL18_pos` / `n_cells_CASP1_pos`、`n_cells_analysed`、`donors_with_*_pos`、`n_donors_analysed` | 合并检出率的分子与分母本身，以及有阳性细胞的供体数 / 进入分析的供体数。图表与报告中的每个百分比都带这两个原始计数，便于判断「0%」是真低还是样本不足。 |
| `logcpm_IL18` 等 | **表达量**（伪bulk）。= log2( 该格子内该基因计数之和 ÷ 该格子内文库大小之和 × 10⁶ + 1 )。 |
| `rank_IL18` 等 | 该细胞类型在**同一分层内**按 logcpm 的排名，1 = 最高。排名只在分层内有意义，不跨分层/物种比较。 |
| `n_cells` | 该「供体 × 细胞类型」格子的细胞个数（<20 时 pct 与 logcpm 置 NA）。 |
| `n_donors` | 该细胞类型在该分层内、格子有效（n_cells ≥ 20）的供体数。<3 不进候选。 |

## 筛选逻辑（§8）

分层内，对第一级谱系（11 类）与第二级作者标签各算一遍：

```
候选 = IL18 分层内排序处于最高三分位
       AND CASP1 分层内排序处于最低三分位
       AND NOT (CASP3 处于分层内最低 20%)      # 分层内细胞类型 <8 时跳过此下限，打 casp3_floor_skipped
       AND pct_IL18 ≥ 0.10                      # §8.4 检出率支撑
       AND n_donors ≥ 3                         # §10.3
```

三分位/20% 为工程取值，无生物学依据。阈值敏感性（§9.6）：候选数在 基线17 → 四分位12 → CASP3-10%:17 → CASP3-30%:16。候选集内部按 CASP3 分层内排序从高到低排优先级（§8.5，仅工程排队，非生物学强度）。

## 判定：两个轴，不用形容词

早先用「确定 / 可能 / 相反」三个词下判断，问题是它把两件不同的事压进了同一个三选一：确定与可能说的是**有多少独立来源支持**，相反说的是**证据指向哪个方向**。结果「证据很强地指向不符合」与「证据很弱地指向不符合」同名；「有一个来源支持」与「从没人测过这个细胞类型」也同名。而「确定」对转录本层面的筛选结果本就是 overclaim，「相反」在中文里又容易被误读成「IL18 低」（它实际说的是 CASP1 存在）。

现在分两列：

| 轴 | 取值 | 含义 |
|---|---|---|
| **方向** | 符合 / 不符合 | 这个细胞类型是否符合「IL18 高 + CASP1 低」这一前提 |
| **证据** | 多源一致 / 多源有分歧 / 单一来源 / 无参考可查 | 支持该方向的独立参考来源有多少、彼此是否一致 |

再加一列 **`label`（观察到的事实）**，直接写出观察结果与来源数，不用形容词。例如：

- `CASP1 多源确认近零（4 个独立来源，无分歧）`
- `4 个来源支持 CASP1 近零，另有 1 个来源显示 CASP1 可测`
- `CASP1 近零，仅 1 个来源支持；被分层内相对规则漏掉`
- `无参考可查；仅本数据方向相符（CASP1 0.75%）`
- `参考层 CASP1 可测到，仅比 IL18 低 3–20 倍`
- `参考层 CASP1 ≥ IL18（方向相反）`

判定阈值（工程取值，无生物学依据，完整比值在 `answer_table.tsv` 里可自行重切）：参考层同源同单位的 IL18:CASP1 比值 ≥20 倍或 CASP1 处于该来源零底线 → 支持「CASP1 近零」；3–20 倍 → CASP1 确实可测；<3 倍 → 方向相反。无匹配参考时按本行自身数据的 IL18:CASP1 判方向，阈值同为 3。

`answer_table.tsv` 的 `tier` / `tier_cn` 与 `il18_casp1_evidence.tsv` 的 `evidence_tier` 保留六档内部判读（A/B/B′/C/D/E）作为留痕，说明每一行为何落入上述方向与证据的组合。

## 交付物

`census_scan.tsv`（阶段一扫描，定位用，不进判定）· `datasets.yaml`（15 分层）· `celltype_map.tsv`（319 行逐条映射）· `il18_casp_celltype_table.tsv`（主表）· `il18_casp_summary.tsv`（供体为单位汇总，谱系级+作者级）· `candidates.tsv`（候选）· `rank_pass_detection_fail.tsv` · `threshold_sensitivity.tsv` · `reference/ext/unified_reference.tsv`（统一参考层，1624 行，带 exact/close/tissue_level/extrapolated 匹配质量）· `il18_casp1_evidence.tsv`（参考层整合判读 + 合并检出率）· `screen_missed_but_reference_supported.tsv`（被分层内相对规则漏掉但有参考支持者，如人肺 AT1）· `answer_table.tsv` / `answer_table.json`（对外答案表，两套检出率口径并列）· `run_log.md` · `figures/`（`answer_v2.png` 候选总览，`answer_human.png` / `answer_mouse.png` 分物种总览，以及 **五器官 × 两物种共 10 张三基因点图** `{human,mouse}/{lung,liver,pancreas,brain_cell,spinal}_dots.png`——每张列出同一个分层内达到 ≥3 供体的全部细胞类型，不论判定如何一律列出，行名颜色即方向/证据；人肺该分层有 39 类，图中按 IL18 检出率取前 20，其余见 `answer_table.tsv`）· `IL18_CASP_五器官筛选.pptx` / `…报告.docx`。

`casp1_zero_check.tsv` 为**已废弃**的旧判读产物，保留仅作留痕：生成它的 `scripts/zero_check.py` 从未读取过参考文件、也从未使用候选自身的 CASP1 值，判读结果是细胞谱系标签的纯函数。请改用 `il18_casp1_evidence.tsv`。

代码在 `scripts/`，可重跑：`phase1/derive_census_scan → phase2_build → make_celltype_map → phase2_aggregate → evidence_integrate → screen_missed → answer_pack → plot_dots`，文档由 `make_deck.js` / `make_report.js` 生成（读 `results/answer_table.json`）。
