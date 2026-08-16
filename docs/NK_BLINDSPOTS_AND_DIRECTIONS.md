# NK 研究盲区穷举 → 方向推导

日期：2026-08-16。方法：先列盲区，再从盲区反推方向，全程对抗性核验新颖性与可测性。

本文档承接本仓库前序项目（四格挖掘协议）的终局结论：**该项目证伪了自己的地基**——
NK 存在三个独立的、NK 特异的技术效应（ambient 负担、differential ambient、
肿瘤侧 RNA content 丢失 p = 6×10⁻⁸），三者全部是 Q4 形状，会把结果推向想要的答案。
因此本次推导有一条硬筛选器：

> **任何以解离 scRNA-seq 转录组为主要读出的 NK 方向，默认判死，除非它自带针对上述三效应的具体辩护。**

---

## 0. 关于本次工作的范围与诚实边界

原计划是一条六段并行工作流（穷举 → 聚类 → 验真 → 生成 → 击杀 → 综合，约 90 个 agent）。
**该工作流失败了，且失败原因与科学内容无关**，记录在案：

- 子 agent 的结构化输出通道在本环境损坏：提交合法的最小 JSON 对象仍被判"缺少必填字段"，
  重试 5 次全部失败；同批子 agent 还报告 Bash/Glob/Grep "permission handler returned malformed input"。
- 更关键的是，子 agent 的 WebSearch **返回结果中没有任何 URL**——即 6 个扫描 agent 实际检索量为零。
  692k tokens 全部消耗在重试上。
- 主进程的 WebSearch/WebFetch 正常。故改为由主进程亲自执行。

**由此产生的范围限制，明确声明：**本文档基于 **26 次定向检索 + 1 篇全文精读**，覆盖 6 个域。
这是一次**有结构的深扫，不是穷举**。下表中每条盲区都标注了证据等级：

- `[V]` = 已检索到文献直接支持该缺口存在或该共识存在
- `[R]` = 由已检索事实推理得出，未单独验证
- `[U]` = 未验证，仅作为待查项列出

未覆盖的域（诚实列出）：NK 代谢抑制、NKG2A/monalizumab 临床、KIR-HLA 教育的定量模型细节、
NK 冻存损伤、NK 组织蛋白质组、adaptive/NKG2C⁺ memory NK、NK engager。

---

## 1. 盲区清单

### A 域：蜕膜 NK（dNK）

| # | 盲区 | 场上的共识 | 为什么还开着 | 证据 |
|---|---|---|---|---|
| A1 | dNK 的个体发生未定：局部 CD34⁺ 祖细胞分化 / 局部 NK 前体 / 外周招募，三说并存 | 综述普遍并列三说而不判 | 人体无法做谱系追踪；小鼠 uNK 不等价 | [V] |
| A2 | dNK1/2/3 分型被广泛沿用，但"验证"主要来自**复用**而非独立挑战 | 该分型已被 CyTOF 与后续 scRNA-seq"证实" | 没有独立团队以证伪为目的重测；亚群的功能归属大多是推断 | [V] |
| A3 | dNK 对螺旋动脉重塑的**因果**贡献有争议，小鼠数据甚至提示 uNK 是**限制**滋养层侵袭而非驱动重塑 | 教科书写作 dNK 驱动重塑 | 人体不可干预；小鼠模型方向相反 | [V] |
| A4 | PTdNK（妊娠训练记忆 dNK）缺少独立团队复制 | Immunity 2018 提出，同组 PLOS Pathog 2024 扩展 | 需要多次妊娠配对样本，样本获取难 | [V] |
| A5 | uNK 临床检测缺乏正常值与预测效度，指南不推荐常规使用，但商业检测在做 | "NK 高导致反复种植失败" | 取样方式、RIF/RPL 定义、外周 vs 子宫 NK 混淆三重不统一 | [V] |
| **A6** | **dNK 是 CD56bright；而 CD56bright 刚被证明是"暂时驻留 + 经淋巴回流"的。dNK 是否真正驻留，成了 2025 之后才存在的新问题** | dNK 是典型组织驻留 NK | 该证据 2025 年才出现，尚未回灌到生殖免疫领域 | [V] |

### B 域：组织驻留 NK 与 ILC1 边界

| # | 盲区 | 场上的共识 | 为什么还开着 | 证据 |
|---|---|---|---|---|
| B1 | 人类 trNK 与 ILC1 无公认可操作定义；肝"trNK"与肝"ILC1"是同一群细胞的两个名字 | 两者是不同群体 | 命名先于机制；跨研究不可比 | [V] |
| B2 | 驻留标志物同时是活化标志物。CD69 在**所有受检器官**的 CD56bright 上都表达；CD49a 随组织变化 | CD69/CD49a/CD103 = 驻留 | 缺少与驻留正交的独立指标 | [V] |
| **B3** | **NK 的组织滞留机制根本未知。** TRM 的 CD69–S1PR1 轴不适用于 NK：NK 用 S1PR5 出境，而**CD69 不拮抗 S1PR5** | 默认沿用 TRM 的滞留模型 | 该不适用性在文献中被陈述过，但没人接着问"那 NK 靠什么留下" | [V] |
| B4 | 人类组织内 NK 的**驻留半衰期**从未被测量 | 肝内驻留 NK "可存活 13 年" | 只有移植嵌合、血管内标记、FTY720 三种代理指标，各有硬限制 | [V] |
| B5 | 若 CD56bright 驻留是暂时的，则大量横断面"trNK"研究拍到的是一股**流**而非一个**池** | 横断面表型 = 驻留群体 | 2025 年才出现的反证 | [R] |

### C 域：肿瘤 NK

| # | 盲区 | 场上的共识 | 为什么还开着 | 证据 |
|---|---|---|---|---|
| **C1** | **"NK 耗竭"框架在小鼠中已被证伪**：功能障碍在入瘤后 48–72 h 内快速发生、**可逆**、与免疫检查点表达**解耦**，且 ICP⁺ NK 反而**更活跃** | NK 耗竭≈T 细胞耗竭，可用检查点阻断逆转 | 术语与临床策略惯性；人体验证缺失 | [V] |
| C2 | 人肿瘤中 NK 的**入胞 / 滞留 / 转化 / 出境 / 死亡**从未被分开测量 | "NK 浸润不足" 被当作单一现象 | 小鼠可用光标记（photolabeling）拆开，人体无等价手段 | [V] |
| C3 | 泛癌 NK 图谱把 DNAJB1/HSPA1A/FOS/JUN 读作肿瘤 NK 生物学，而这基本就是解离应激基因集 | 该signature 是肿瘤浸润 NK 的特征 | 没人跑分离 artefact 与生物学的跨谱系对照 | [V] 本仓库 Phase 0 |
| C4 | 解离式 scRNA-seq 对 NK 的表型读出被三个 NK 特异技术效应污染 | 单细胞数据可直接读 NK 状态 | 需要共乳液设计或非解离读出才能拆开 | [V] 本仓库实测 |
| C5 | NK 位置排斥（多见于间质/侵袭边缘）被描述，但其成因未与"入胞失败"分离 | NK 进不去肿瘤 | 同 C2 | [V] |
| C6 | **不存在**蜕膜 NK 与肿瘤 NK 的同方案配对比较 | "肿瘤 NK 呈蜕膜样"已是流行说法 | 两个领域样本流程与团队互不重叠 | [V] 检索确认无此研究 |

### D 域：趋化因子

| # | 盲区 | 场上的共识 | 为什么还开着 | 证据 |
|---|---|---|---|---|
| D1 | 人组织内的趋化因子**梯度**无法测量。小鼠耳皮肤的 CCL21 单体梯度做到过，人组织无等价 | 梯度是既定事实 | 需活体成像 + 固定化梯度探针 | [V] |
| **D2** | **CD26/DPP4 的 N 端截短把 CXCL10 变成 CXCR3 拮抗剂**；能区分完整型与截短型的免疫亲和-MS 方法（ISTAMPA 类）已存在并验证于**血浆**和**滑液**——**从未用于人肿瘤组织或蜕膜** | CXCL9/10/11 转录本丰度 ≈ NK 招募驱动力 | 需要 MS 平台，转录组学界与蛋白型学界不交叉 | [V] |
| D3 | 配体与受体的 mRNA 都不能代表蛋白：NK 上 CCR1 有 mRNA 无表面蛋白；CXCR3 受强翻译后调控并内化 | 单细胞配受体推断可信 | 配受体推断工具主要以共表达为基准，而非蛋白真值 | [V] |
| D4 | ACKR 清道夫塑造梯度：滋养层 ACKR2 构成"胎盘趋化因子防火墙"；肿瘤中 ACKR2 由 HIF-1α 诱导并**约束 NK 迁移**（小鼠） | ACKR 是背景变量 | 人肿瘤中 ACKR2 与 NK 位置的因果未测 | [V] |
| D5 | **出境轴（S1PR5）在肿瘤中基本无人研究** | NK 浸润问题=招募问题 | 全领域注意力集中在入口 | [R] |

### E 域：模型与方法

| # | 盲区 | 场上的共识 | 为什么还开着 | 证据 |
|---|---|---|---|---|
| E1 | Ncr1/NKp46 不是 NK 特异的（ILC1/ILC3 亦表达），且 **NKp46 本身是 ILC1 发育所必需**——所以 Ncr1-Cre 在扰动 NK 的同时扰动了对照组 | Ncr1-Cre = NK 特异工具 | 没有真正 NK 特异的 Cre | [V] |
| E2 | 37 ℃ 酶解会诱导**恰好被读作肿瘤 NK 生物学的那组应激基因**，并同时耗竭敏感细胞；冷蛋白酶可多回收 NK | 解离方案是技术细节 | 已有证据但未回灌到 NK 领域的结论重估 | [V] |
| E3 | 人体内 NK 杀伤从未被直接观察；活体成像证据几乎全来自小鼠 | NK 在人肿瘤内杀伤 | 无人活体成像手段 | [V] |
| E4 | 流式门控循环论证：用 CD56/CD16 定义 NK，再"发现" CD56/CD16 生物学 | — | 定义与发现同源 | [R] |

---

## 2. 验真：哪些"新方向"其实已经死了

以下均为本次检索中**确认被占位**的想法，列出以防重复投入：

| 想法 | 判决 | 占位工作 |
|---|---|---|
| 肿瘤 NK 呈"蜕膜样"（CD9/CD49a/促血管） | **已做** | Cancer Discovery 2021 观点文；iScience 2022 移植后 CD9⁺ decidual-like NK；HCC CD49a⁺ NK 预后差 |
| ACKR2 约束肿瘤内 NK 迁移 | **已做（小鼠）** | ACKR2 缺失小鼠转移减少、KLRG1⁺ NK 招募增加；2024 年 HIF-1α 诱导 ACKR2 |
| DPP4 截短 CXCL10 抑制 NK 招募 | **已做** | 小鼠 HCC 中 DPP4i 依赖 NK 与 CXCR3；人体血浆药理学已验证 |
| TGF-β 驱动 NK→ILC1 转化 | **已做** | SMAD4 非经典通路；人膀胱癌 2024 |
| NK 入瘤后快速功能丧失 | **已做（小鼠）** | Nat Commun 2024；Sci Adv 2024 |
| 给 CAR-NK 装 CXCR2/CXCR4 以改善实体瘤归巢 | **重度占位** | CAR-T/CAR-NK 均已大量发表并进临床 |
| 再挖一轮肿瘤 NK scRNA-seq | **可测性判死** | 本仓库实测的三个 NK 特异技术效应 |

---

## 3. 存活方向

### 方向 1（最强）：把"NK 被抑制"的读出从转录组换成**突触极化能力**，并用蜕膜做阳性对照

**来自盲区**：C1 + C4 + C6 + A6

**观察到的巧合**（这是本次推导的核心）：

- dNK **有**穿孔素与颗粒酶，能与靶细胞形成结合并极化 CD2/LFA-1/actin，
  却**无法把 MTOC 与颗粒极化到突触**——裂解活性仅约外周 NK 的 15%。机器完好，投递被卡。
- 肿瘤 NK 的功能障碍**快速、可逆、与检查点解耦**，且 IL-15 可恢复。同样是"机器完好、开关被压"。

**假设（可证伪）**：肿瘤 NK 的"功能障碍"与 dNK 的生理性非裂解**是同一个病灶**——
一个位于 MTOC/颗粒极化步骤的、转录后的、可逆的投递阻断；而不是转录程序的重编，也不是检查点耗竭。

**读出**：成像流式（imaging flow cytometry）或共聚焦，对同一标准靶细胞测量
①结合率 ②MTOC–突触距离 ③颗粒酶 B 突触内占比 ④脱颗粒。
**这正是本方向的战略价值**：这是**单细胞影像读出，对 ambient RNA 和 RNA content 丢失完全免疫**，
绕开了杀死上一个项目的三个 NK 特异技术效应。

**证伪条件**：肿瘤 NK 极化正常而卡在更早（结合）或更晚（脱颗粒）的步骤；
或极化缺陷只跟 CD56bright 身份走、与组织无关。

**最小可行实验**：蜕膜（择期终止妊娠，dNK 占蜕膜白细胞 50–70%，细胞数不是限制）
+ 外周血对照，跑极化 assay，先确认该 assay 在新鲜组织 NK 上的动态范围。不需要肿瘤样本即可启动。

**新颖性**：dNK 突触缺陷 2005 年已知；肿瘤侧突触研究亦存在（PNAS 2025：癌细胞用 actin 把抑制性配体极化到突触）。
**未被占的是**：把两者作为同一病灶提出，并把"极化能力"确立为**取代转录组的 NK 抑制主指标**。

**对前序项目的意义**：上一个项目的停止规则 S4 判定"Q4 阳性对照在原理上不可得"——
即拿不到一个"确知被抑制但活着"的真值样本。**蜕膜就是那个真值样本。**

---

### 方向 2：出口，而不是入口

**来自盲区**：B3 + C2 + D5 + A6

**假设**：蜕膜中 NK 的富集与肿瘤中 NK 的稀少，主要由**出境速率**决定，而非招募速率；
蜕膜强制滞留，肿瘤放任出境。

**为什么现在可做**：三个事实刚刚拼齐——CD56bright 会经淋巴回流（2025）；
NK 的出境受体是 S1PR5；**CD69 不拮抗 S1PR5**，所以 TRM 的滞留模型对 NK 不成立，
NK 靠什么留下是一个**公开的空位**。

**读出**：配对样本（蜕膜 / 肿瘤 / 癌旁 / 血）的 S1PR5 表面蛋白与**功能性 S1P 趋化反应**；
体外 FTY720 类似物阻断；有条件时采集人输出淋巴。

**证伪条件**：各腔室间 S1PR5 水平与 S1P 反应性无差异。

**已知风险（一次打击）**：S1PR5 抗体质量差是公认问题 → 补救是以**功能性 S1P transwell 反应**为主指标，
表面染色仅作辅助。

---

### 方向 3：把 CXCL9/10/11 的"激动剂 : 拮抗剂"比例在人组织里真正测出来

**来自盲区**：D2 + D3

**假设**：人肿瘤中 NK 的浸润程度与 CXCL10/CXCL12 的**完整型 : 截短型比值**相关，
而与趋化因子的总转录本或总蛋白量不相关；蜕膜维持一套不同的蛋白型组合。

**读出**：免疫亲和-LC-MS/MS 蛋白型分析（ISTAMPA 类方法）用于肿瘤组织间质液与蜕膜组织。

**新颖性核验**：该方法已验证于人血浆（CXCL12α 五种蛋白型）与类风湿滑液（CXCL10），
**未见用于人肿瘤组织或蜕膜**。属于"成熟方法 + 未占应用"。

**为什么重要**：它给出了一个直接解释——为什么 CXCR3 轴与 NK 浸润的相关性一直偏弱，
以及为什么 DPP4 抑制剂在小鼠里有效。同时它是一个**蛋白层读出**，同样绕开 scRNA-seq 的三个陷阱。

---

## 4. 贯穿性结论

三个存活方向里有两个**换的是仪器，不是问题**。

这不是巧合。本仓库前序项目的每一次死亡——四格框架、over-dispersion 轴、Q2、
Q4 阳性对照——都不是死于问题问得不好，而是**死于读出手段测不到那个前提**。
NK 恰恰是最不适合用解离转录组去读的谱系：它稀少、脆弱、ambient 负担第二高、
且在肿瘤侧特异性丢 RNA。

因此本次推导的净输出是一句判断：

> **NK 领域的下一步不在于找新基因，而在于换掉那个把 NK 生物学与 NK 假象混在一起的仪器。
> 而蜕膜提供了这个领域一直缺的东西——一个"被抑制但活着"的 NK 真值样本。**

---

## 5. 第一步动作

1. **本周可做，零新样本**：把 A6 这个问题写成一页——"dNK 是真驻留还是被滞留的过路 CD56bright"，
   并对已公开的蜕膜 scRNA-seq 检查 S1PR5 / CX3CR1 / SELL 与滞留相关基因的表达，
   注意这只能作为**线索**，不能作为结论（三个技术效应同样适用于蜕膜数据，且蜕膜的
   滋养层/基质 ambient 负担可能更重）。
2. **建立极化 assay**：先在外周血 NK + K562 上跑通成像流式极化读出，确定动态范围与所需细胞数。
3. **接触样本来源**：择期终止妊娠蜕膜（方向 1、3），配对肿瘤/癌旁（方向 1、2、3）。
4. **找 MS 合作方**：方向 3 的全部难度在平台，不在设计。

---

## 附：主要文献

- Transient tissue residency and lymphatic egress define human CD56bright NK cell homeostasis — *Nat Immunol* 2025 — https://www.nature.com/articles/s41590-025-02290-9 ／ PMC12571907
- Human decidual NK cells form immature activating synapses and are not cytotoxic — *PNAS* 2005 — https://pnas.org/content/102/43/15563.short
- Tumor-induced NK cell dysfunction is a rapid and reversible process uncoupled from the expression of immune checkpoints — *Sci Adv* 2024 — https://www.science.org/doi/10.1126/sciadv.adn0164 ／ PMC11352832
- Rapid functional impairment of natural killer cells following tumor entry limits anti-tumor immunity — *Nat Commun* 2024 — https://www.nature.com/articles/s41467-024-44789-z
- Sphingosine 1-phosphate receptor 5 (S1PR5) regulates the peripheral retention of tissue-resident lymphocytes — *JEM* 2022 — https://rupress.org/jem/article/219/1/e20210116/212717
- Placental chemokine compartmentalisation (ACKR2 防火墙) — *PLOS Biology* 2019 — https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000287
- The Atypical Chemokine Receptor Ackr2 Constrains NK Cell Migratory Activity and Promotes Metastasis — PMC6176105
- Inhibition of DPP4 activity in humans establishes its in vivo role in CXCL10 post-translational modification — PMC4888857
- From ELISA to Immunosorbent Tandem Mass Spectrometry Proteoform Analysis (ISTAMPA) — PMC7991300
- Quantification of increased biologically active CXCL12α plasma concentrations after ACKR3 antagonist treatment in humans — *Clin Transl Sci* 2024 — PMC10818162
- Cancer cells suppress NK cell activity by actin-driven polarization of inhibitory ligands to the immunological synapse — *PNAS* 2025 — https://www.pnas.org/doi/10.1073/pnas.2503259122
- Number and function of uterine NK cells in recurrent miscarriage and implantation failure: systematic review and meta-analysis — *Hum Reprod Update* 2022 — https://academic.oup.com/humupd/article/28/4/548/6545823
- Systematic assessment of tissue dissociation and storage biases in scRNA/snRNA-seq — *Genome Biol* 2020 — https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6
- Dependence of innate lymphoid cell 1 development on NKp46 — *PLOS Biology* — https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2004867
- Decidual-Like NK Cell Polarization: From Cancer Killing to Cancer Nurturing — *Cancer Discov* 2021 — https://aacrjournals.org/cancerdiscovery/article/11/1/28/2874
- Trained Memory of Human Uterine NK Cells Enhances Their Function in Subsequent Pregnancies — *Immunity* 2018 — https://pubmed.ncbi.nlm.nih.gov/29768178/
