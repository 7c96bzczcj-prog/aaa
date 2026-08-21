# PIEZO1 蛋白层与功能读数的可行性核查 —— 任务 P-1，Q7

## §0 必须置于开头的限制

> PIEZO1 是机械敏感离子通道，其活性取决于膜张力与局部力学环境，**不仅取决于表达量**。
> 同样的转录本水平，在不同环境下的通道活性可以相差很大。
>
> **本任务给出的是「表达底数」，不是「机械感受能力」。**

本文件是这条限制的直接后果：既然活性不只取决于表达量，而表达量在 NK 中又落地板带
（[`PIEZO1_BASELINE.md`](PIEZO1_BASELINE.md) Q1），那么蛋白层与功能读数是否可及，
就不是一个次要问题。

**本文件是可行性核查，不是数据分析，也不提任何方向。**

---

## 一、标注约定（自补，按 §3 报告）

任务规定「没找到」记 `NOT_SEARCHED`。执行处把它拆成两个标签，因为两者的信息量差别很大：

| 标签 | 含义 |
|---|---|
| `NOT_FOUND` | **已检索、未找到**。检索式列于该条目下 |
| `NOT_SEARCHED` | **未检索**。理由列于该条目下 |
| `[待核实]` | 找到了出处，但无法读到原文核实其内容 |
| `VENDOR_CLAIM` | 仅有厂商说明书为据，无独立验证 |

这是对规定标签的细化而非替换；`NOT_FOUND` 与 `NOT_SEARCHED` 合起来对应规定中的
`NOT_SEARCHED`。列表中永不留空。

---

## 二、蛋白层检测手段

### 2.1 Western blot / 免疫沉淀：有经同行评议的成熟方案，但对样品制备苛刻

*Bio-protocol* 2025，"Western Blotting and Immunoprecipitation of Native Human PIEZO1
Channels"（PMC12304459）。要点：

| 项 | 内容 |
|---|---|
| 抗体 | 克隆 **2-10**，Novus Biologicals **NBP2-75617** |
| 验证 | 在 `HEK293T-Piezo1⁻/⁻` 敲除细胞上作阴性对照 |
| 物种 | **不识别小鼠 PIEZO1**（原文明示）——若要做小鼠 NK，须另选抗体 |
| 已验证细胞 | HEK293T、MCF-7、HeLa、LNCaP、BJ-5ta、原代人心脏成纤维细胞、HUVEC、猪主动脉内皮细胞 |
| **未涉及** | **任何免疫细胞，包括 NK** |
| 上样量 | WB 每道 10–30 µg 总蛋白；IP 需一个 60/100 mm 皿约 80% 汇合度；接质谱需 ≥ 2 个 150 mm 皿 |
| 关键约束 | 裂解全程冰上，**不可 95 °C 煮沸**（煮沸导致聚集、信号大幅损失）；TCEP 10 mM + NEM 20 mM；改良 RIPA，SDS 仅 0.1% |

**该方案存在的原因本身就是一条结论**：作者明言「常规样品制备方法似乎不适用于 PIEZO1」——
它是一个巨大的多次跨膜蛋白，疏水面积大、翻译后修饰多，容易聚集并产生非特异相互作用。
**所以 PIEZO1 的 WB 阴性结果在没有用该方案做过之前不能当作阴性**。

对 NK 的换算：10–30 µg 总蛋白对应约 **(0.5–3)×10⁵ 个淋巴细胞**，
按 1×10⁶ PBMC ≈ 100–200 µg 总蛋白这一常用经验值折算——**该换算系数 `[待核实]`**，
未在本次检索中找到一手出处。以原代人 NK 的常规产量（一份血沉棕黄层可得 10⁶–10⁷）衡量，
上样量本身不是瓶颈；瓶颈在抗体与制备条件。

### 2.2 流式细胞术：有胞外表位抗体，但 NK 上无已验证方案

- 存在针对 PIEZO1 **胞外表位**的抗体（Alomone `APC-087` 系列；Proteintech `28511-1-AP`
  标注为 extracellular domain），原理上可用于**活细胞、不固定**的表面染色与流式。
  以上均为 `VENDOR_CLAIM`——本次未取得任何独立验证数据。
- **物种坑**：检索结果显示 Alomone 的 FITC 偶联变体 `APC-087-F` 不识别人 PIEZO1
  （针对大鼠／小鼠）。同一产品线不同偶联型号的物种反应性可以不同，**逐个型号核实**。
- **NK 细胞上经验证的 PIEZO1 流式方案：`NOT_FOUND`。**
  检索式："PIEZO1 antibody flow cytometry human NK cells validated surface staining"。
  找到的是通用抗体页与 PIEZO1-NK 功能论文，没有一篇给出 NK 上的染色对照
  （同型对照、敲除／敲低对照、荧光减一）。
- 一个必须预先说明的判读困难：PIEZO1 在 NK 中的转录本落地板带，
  若流式染出一个整体右移的单峰，**无法区分弱阳性与背景**——
  该实验必须自带 PIEZO1 敲除或敲低的阴性对照，否则不产出可用信息。

### 2.3 质谱：路径存在，本次未查证 NK 中是否检出 PIEZO1

- ImmProt（immprot.org，Rieckmann 等，*Nat Immunol* 2017）对 28 种原代人造血细胞群
  给出蛋白拷贝数，含 NK。**其中是否含 PIEZO1：`NOT_SEARCHED`**——
  该资源需交互式查询，本次未查。这是一个成本极低的后续动作，建议先做。
- 原代人 T/B/NK 的 DIA 谱图库（*Sci Data* 2024, s41597-024-03721-2）存在，
  但本次因出版方跳转鉴权未能读到正文，**其覆盖是否含 PIEZO1：`[待核实]`**。
- 一般性预期（不作为证据）：多次跨膜、低丰度的大蛋白在自下而上蛋白组中通常覆盖差。

### 2.4 报告系统

| 系统 | 内容 | 对 NK 的适用性 |
|---|---|---|
| **GenEPi** | Piezo1 融合低亲和力 GCaMP 的基因编码报告器，把 GCaMP 定位到 Piezo1 通道口的 Ca²⁺ 微区，从而只报告 **Piezo1 依赖**的钙信号；已建斑马鱼转基因系（Yaganoglu 等，*Nat Commun* 2023，PMC10356793） | 需转基因导入。原代 NK 转染／转导效率低；NK-92 可行性较高。**NK 上的既有应用：`NOT_FOUND`** |
| **Piezo1-tdTomato 敲入小鼠** | JAX **029214**，Piezo1-tdTomato 融合报告等位基因 | 可直接用于小鼠 NK 的流式与成像，绕开抗体问题。**NK 上的既有应用：`NOT_FOUND`** |
| **Piezo1-flox** | JAX **029213**（`Piezo1^tm2.1Apat`），loxP 侧翼 exon 20–23 | 与 NK 限制性 Cre 组合即得 NK 特异敲除 |
| **NK 限制性 Cre** | `Ncr1`(NKp46)-Cre / -iCre 品系已建立并广泛使用（Narni-Mancinelli 等；Merzoug 等 *Eur J Immunol* 2014） | 两个组件都现成 |
| **`Ncr1-Cre × Piezo1^fl/fl` 的已发表研究** | **`NOT_FOUND`** —— 检索式："Piezo1 fl/fl Ncr1-Cre NK cell specific conditional knockout mouse"。两个组件均存在，未检到已发表的组合 | — |

---

## 三、功能读数

### 3.1 Yoda1 + 钙成像

- **Yoda1** 是 PIEZO1 的小分子激动剂，对**人与小鼠 PIEZO1 均有效**，
  机制上是能量学上的门控修饰剂而非孔道开放剂（Botello-Smith 等 *Nat Commun* 2019；
  Lacroix 等 *Nat Commun* 2018；Nosyreva 等 *PNAS* 2022）。
- 标准读数为 Ca²⁺ 敏感染料（Fluo-4/Fluo-8 AM）负载后灌流 Yoda1 的活细胞成像，
  是成熟且低门槛的方案。
- **在 NK 上**：Yanamandra 等 2024 使用 Yoda1 激活 NK 的 PIEZO1 并观察到杀伤与浸润增强
  （综述转述：*Front Immunol* 2026, PMC13287921）。
  **所用浓度、染料、成像参数：`NOT_FOUND`**——EJI 正文（Wiley 403）、bioRxiv 预印本
  （Cloudflare 限流）、机构库副本（bot 拦截）三条路径均未取得全文；两篇综述转述时均未给出方法细节。
- **在原代人 NK 或 NK-92 上、独立于该文的 Yoda1 钙成像方案：`NOT_FOUND`。**
  检索式："Yoda1 calcium imaging NK-92 primary human NK cells PIEZO1 agonist Fluo-4"。
- Yoda1 自身的两个已知缺陷，必须在设计时计入：溶解度差；剂量-反应呈钟形，
  高浓度反而抑制。已有 4-苯甲酸修饰的改良类似物（*Br J Pharmacol* 2024，PMC10952572）。
- **判读上的关键点**：NK 是钙信号极其活跃的细胞（脱颗粒依赖钙）。
  「Yoda1 引起钙内流」本身**不足以**证明 PIEZO1 参与——须有 PIEZO1 敲除／敲低的
  阴性对照，或至少一条独立的 PIEZO1 扰动臂。Yoda1 的脱靶在高浓度下有报道。

### 3.2 膜片钳

- poking 法诱发的 PIEZO1／PIEZO2 电流是标准做法，但**其激活阈值与失活动力学
  对机械刺激参数的细微差别高度敏感**（*J Gen Physiol* 2024，PMC11734767），
  跨实验室可比性差。
- 该方法为贴壁细胞设计。NK 是小体积悬浮细胞，**poking 构型的适配性存在困难**——
  这是可行性上的实质障碍，不是措辞上的谨慎。
- **原代人 NK 或 NK-92 上的 PIEZO1 机械电流记录：`NOT_FOUND`。**
  检索式："PIEZO1 patch clamp primary human NK cells mechanically activated current NK-92"。
  检到的最接近者是白血病细胞系上内源 Piezo1 的单通道记录（PMC8346046）。

---

## 四、GsMTx4 的选择性问题（任务要求写明）

**GsMTx4 不是 PIEZO1 特异的抑制剂。用它得到的抑制结果不能归因于 PIEZO1。**

四条已核实的事实：

1. **抑制谱是一类通道，不是一个通道。** GsMTx4 是 34 氨基酸的狼蛛毒肽，
   抑制阳离子通透的机械敏感通道，**已知涵盖 PIEZO1、PIEZO2、TRPC1、TRPC6**
   （Tocris/Bio-Techne 产品说明；Bae 等 *Biochemistry* 2011, PMID 21696149）。
2. **机制是脂质介导的门控修饰，不是锁钥式结合。** 其 **D-型对映体 enGsMTx4 活性相当**
   （Suchyna 等 *Nature* 2004, 430:235），效应随膜厚与通道长度的疏水错配程度增强——
   它作用于通道周边的边界脂，改变局部曲率并把通道推向关闭态。
   **这一点决定了「脱靶清单」原则上无法穷举**：任何张力门控的阳离子通道都可能受影响，
   而不只是已被点名的那几个。
3. **它不阻断 TREK-1**（K2P 类力学门控钾通道）——所以它也不是「机械敏感通道」的全谱抑制剂。
   已知谱与未知谱之间没有清晰边界。
4. **在体内／离体的贡献可以由脱靶主导**：一项 2025 年的研究报告 GsMTx-4 抑制运动升压反射
   **主要经由 TRPC 而非 Piezo**（PMC12487608）。

**本仓库数据能对这一条补一句，但补得有限。** 在 `GSE154826` 的 NK 中，
GsMTx4 已被点名的两个 TRPC 脱靶几乎不表达：

| 基因 | 肿瘤 CPM / 检出率 / 零计数供者 | 癌旁 CPM / 检出率 / 零计数供者 |
|---|---|---|
| `TRPC1` | 0.0 / 0.07% / 23 of 27 | 0.0 / 0.11% / 22 of 27 |
| `TRPC6` | 2.0 / 0.34% / 21 of 27 | 1.1 / 0.18% / 19 of 27 |
| `PIEZO1`（对照） | 15.6 / 2.62% / 1 of 27 | 16.3 / 2.73% / 0 of 27 |

这**收窄**了 NK 上 GsMTx4 的歧义，但**不消除**它，理由有三：
（a）这是转录本，不是蛋白；（b）这三个读数全在地板带内，本身不可靠；
（c）按第 2 条，GsMTx4 的机制是脂质介导的，其影响不限于已点名的通道。

**结论**：GsMTx4 可以作为「机械敏感阳离子通道整体是否参与」的探针，
**不能**作为「PIEZO1 是否参与」的探针。要归因到 PIEZO1，须有基因层面的扰动
（siRNA／shRNA／CRISPR／条件敲除），并配 Yoda1 激动臂作方向相反的验证。

---

## 五、原始论断的核实状态

Yanamandra 等，*Eur J Immunol* 2024;54:e2350693，PMID 38279603，DOI 10.1002/eji.202350693。

| 项 | 状态 |
|---|---|
| 论文存在、题名、作者、期刊卷期 | **已核实**（Europe PMC 记录） |
| 摘要原文措辞 | **已核实**："PIEZO1 as the predominantly expressed mechanosensitive ion channel **among the examined candidates** in NK cells" |
| 所检候选通道清单 | **`[待核实]`** —— 三条取全文路径均失败（Wiley 403；bioRxiv Cloudflare 1015；机构库 bot 拦截） |
| 判定「主要表达」所用的测量方法（qPCR？RNA-seq？蛋白？电流？） | **`[待核实]`**，同上 |
| PIEZO1 蛋白在 NK 中是否被该文检出、用何抗体 | **`[待核实]`**，同上 |
| 该文的 Yoda1 / GsMTx4 浓度 | **`[待核实]`**，同上 |
| 两篇综述是否复述了上述细节 | **已核实：都没有。** *Cancers* 2024（PMC11048000）与 *Front Immunol* 2026（PMC13287921）均只转述结论，不给候选清单、不给方法、不给浓度 |

如实记录一个观察：**该论断在综述中被逐层去掉限定语**——
原文的「在所检候选之中」，到综述里变成「NK 中主要表达的机械敏感离子通道」。
[`MECHANOSENSOR_PANEL.md`](MECHANOSENSOR_PANEL.md) 的实测直接反对的是后一种表述。

---

## 六、可行性小结

| 读数 | 是否可及 | 主要障碍 |
|---|---|---|
| WB / IP 检测人 PIEZO1 蛋白 | **可及** | 须用 KO 验证过的克隆 2-10 与不煮沸的专用制备流程；NK 上无先例 |
| WB 检测小鼠 PIEZO1 | 受限 | 上述抗体不识别小鼠，须另选并另行验证 |
| 流式（胞外表位） | 原理可及，**NK 上无已验证方案** | 逐型号的物种反应性；低丰度下必须自带 KO/KD 阴性对照 |
| 质谱（NK 中是否检出 PIEZO1） | **未查证** | ImmProt 交互查询即可，成本最低，建议先做 |
| Yoda1 + 钙成像 | **可及且门槛低** | NK 钙信号本底强，须有 PIEZO1 扰动臂才能归因；Yoda1 剂量呈钟形 |
| 膜片钳（NK 上的机械电流） | **困难** | poking 构型面向贴壁细胞；小体积悬浮细胞适配困难；无 NK 先例 |
| GenEPi 报告器 | 可及性中等 | 需转基因；原代 NK 导入效率低，NK-92 更现实 |
| Piezo1-tdTomato 报告小鼠 | **可及**（小鼠） | 现成品系，绕开抗体 |
| `Ncr1-Cre × Piezo1^fl/fl` | **组件齐备，未见已发表组合** | 需自建 |
| GsMTx4 作为 PIEZO1 探针 | **不可用于归因** | 见第四节 |

**一句话**：若目标是「PIEZO1 在 NK 中是否有功能」，蛋白与功能读数都够得着，
但每一条都需要基因层面的扰动作为归因臂；若目标是「PIEZO1 在 NK 中表达多少」，
液滴式 scRNA-seq 已经给到它的极限（地板带），再往下要靠 bulk RNA、qPCR 或蛋白。

---

## 七、检索式与来源

本次实际执行的检索式（全部为 2026-08 检索）：

1. `Yanamandra PIEZO1 natural killer cells mechanosensitive European Journal of Immunology 2024`
2. `PIEZO1 antibody flow cytometry human NK cells validated surface staining`
3. `Yoda1 calcium imaging NK-92 primary human NK cells PIEZO1 agonist Fluo-4`
4. `GsMTx4 selectivity not specific Piezo1 TRPC1 TRPC6 mechanosensitive channel inhibitor peptide`
5. `Piezo1-tdTomato reporter mouse GenEPi genetically encoded PIEZO1 indicator fluorescent`
6. `PIEZO1 patch clamp primary human NK cells mechanically activated current NK-92`
7. `PIEZO1 antibody specificity problems non-specific western blot knockout validation`
8. `Piezo1 fl/fl Ncr1-Cre NK cell specific conditional knockout mouse natural killer`
9. `Rieckmann immune cell proteome ImmProt NK cell PIEZO1 copy number mass spectrometry`
10. `GsMTx4 D-enantiomer equally active gating modifier lipid bilayer Suchyna`
11. Europe PMC REST：`DOI:"10.1002/eji.202350693"`、题名检索（预印本 `PPR637801`）
12. Ensembl REST `/homology/symbol/human/{gene}`（直系同源核实，19 基因 + 2 对照）
13. Human Protein Atlas `search_download` API（`blood_RNA_NK-cell` 等，19 基因）

主要来源：

- Yanamandra AK 等. *Eur J Immunol* 2024;54:e2350693. PMID 38279603. DOI 10.1002/eji.202350693 —— <https://onlinelibrary.wiley.com/doi/10.1002/eji.202350693>
- Western Blotting and Immunoprecipitation of Native Human PIEZO1 Channels. *Bio-protocol* 2025 —— <https://pmc.ncbi.nlm.nih.gov/articles/PMC12304459/>
- Suchyna TM 等. Bilayer-dependent inhibition of mechanosensitive channels by neuroactive peptide enantiomers. *Nature* 2004;430:235 —— <https://www.nature.com/articles/nature02743>
- Bae C, Sachs F, Gottlieb PA. The mechanosensitive ion channel Piezo1 is inhibited by the peptide GsMTx4. *Biochemistry* 2011. PMID 21696149 —— <https://pubmed.ncbi.nlm.nih.gov/21696149/>
- GsMTx4 产品说明（抑制谱：Piezo、TRPC1、TRPC6）—— <https://www.tocris.com/products/gsmtx4_4912>
- GsMTx-4 抑制运动升压反射主要经由 TRPC —— <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12487608/>
- Yaganoglu S 等. Highly specific and non-invasive imaging of Piezo1-dependent activity using GenEPi. *Nat Commun* 2023 —— <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10356793/>
- Piezo1-tdTomato 报告小鼠 JAX 029214 —— <https://www.jax.org/strain/029214>；Piezo1-flox JAX 029213 —— <https://www.jax.org/strain/029213>
- Merzoug LB 等. Conditional ablation of NKp46⁺ cells using Ncr1^greenCre. *Eur J Immunol* 2014 —— <https://pubmed.ncbi.nlm.nih.gov/25142413/>
- Yoda1 作用机制 —— <https://www.nature.com/articles/s41467-019-12501-1>、<https://www.nature.com/articles/s41467-018-04405-3>、<https://www.pnas.org/doi/10.1073/pnas.2202269119>
- 4-苯甲酸修饰的改良 PIEZO1 激动剂 —— <https://pmc.ncbi.nlm.nih.gov/articles/PMC10952572/>
- poking 诱发 PIEZO1/PIEZO2 电流对刺激参数敏感 —— <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11734767/>
- NK 细胞力学感受综述 *Cancers* 2024 —— <https://pmc.ncbi.nlm.nih.gov/articles/PMC11048000/>
- Piezo1 in immune cells 综述 *Front Immunol* 2026 —— <https://pmc.ncbi.nlm.nih.gov/articles/PMC13287921/>
- Rieckmann JC 等. *Nat Immunol* 2017（ImmProt）—— <https://www.nature.com/articles/ni.3693>
- Human Protein Atlas，PIEZO1 —— <https://www.proteinatlas.org/ENSG00000103335-PIEZO1/immune+cell>
