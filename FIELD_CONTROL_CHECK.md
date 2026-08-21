# FIELD_CONTROL_CHECK.md —— S-5 §4：既有肿瘤 Kaede 工作做过非肿瘤对照吗

> ## 结论
> **没有。两篇肿瘤 Kaede 工作都未做过时间分辨的非肿瘤对照，记 `NOT_DONE`。**
> Dean 2024 做过一项**静态**的跨组织表型比较（流式，CD49a/CD11b，对比脾／肝／小肠／结肠），但那不是同型对照——它比的是「瘤内 NK 与其它组织的 NK/ILC1 长什么样」，**不是**「进入非肿瘤组织的细胞随时间发生什么」。
> 因此 §4 的那句判断成立：**「NK 入瘤后功能丧失」这一解读，至今未与「NK 进入任何组织后的驻留化」区分开。**

---

## 1. 核查对象与方法

| 论文 | 标识 | 全文 |
|---|---|---|
| **Dean 2024** —— *Rapid functional impairment of natural killer cells following tumor entry limits anti-tumor immunity* | Nat Commun 2024，PMID 38267402，PMC10808449 | 开放获取，全文 XML 实读（103,883 字符） |
| **Li 2022** —— *In vivo labeling reveals continuous trafficking of TCF-1⁺ T cells between tumor and lymphoid tissue* | J Exp Med 2022，PMID 35472220，PMC9048291 | 开放获取，全文 XML 实读（118,559 字符） |

方法：全文关键词计数 + 逐个上下文阅读。

---

## 2. Dean 2024（NK，MC38）

### 2.1 做过什么

**唯一与非肿瘤组织有关的一项**（原文，**[V]**）：

> 「To extend this analysis and compare intratumoral NK cells with the NK/ILC1 populations found in healthy tissue, we assessed the **CD3⁻ NK1.1⁺ cells in MC38 and B16-F10 tumors, alongside spleen, liver, small intestine and colon** (Supplementary Fig. 4). These data indicated that the tumor NK cell compartment was **phenotypically distinct** to the ILC1 populations identified in non-tumor tissues」

**性质**：**流式细胞术的静态表型比较**（CD49a vs CD49b、CD49a vs CD11b）。

### 2.2 为什么这不是 §4 要问的那个对照

| 项 | Dean 2024 的跨组织比较 | §4 要求的对照 |
|---|---|---|
| 是否光转换非肿瘤组织 | **否** | 是 |
| 是否有时间轴 | **否**——单一时点的表型快照 | 是（新到 vs 滞留） |
| 读数 | 流式，2–3 个蛋白标记 | 转录组 |
| 所答的问题 | 瘤内 NK 与其它组织的 NK/ILC1 **长得像不像** | 细胞**进入某组织之后**发生什么 |

**「表型不同」恰恰不能回答本问题**：瘤内 NK 与肝 ILC1 稳态群体表型不同，与「NK 进入肝之后 48 小时内是否发生同样的驻留化程序切换」是两件事。**后者从未被测。**

### 2.3 一处对 S-4 有直接意义的旁证

Dean 2024 自己对瘤内 NK 随时间变化的描述（**[V]**）：

> 「there was an upregulation of inhibitory receptors including *Pdcd1*, *Lag3*, *Havcr2*, *Cd96*, **a change in migration-associated transcripts with a marked reduction in *Cxcr4* and *Itgam* expression, but an increase in *Itga1***, as well as **loss of expression of *Ccl3*, *Ccl4*, *Ccl5***」

以及：

> 「Other down-regulated transcription factors included *Irf8*, **_Klf2_**, *Myc*」

> **这与 S-4 在另一个数据集（`E-MTAB-10176`）上算出的残差顶端逐个吻合**：`Itgam`（99.94 百分位）、`Cxcr4`（99.75）、`Klf2`（98.59）、`Ccl3`（99.86）、`Ccl4`（99.59）、`Ccl5`（99.56）下降，`Itga1` 上升。
>
> **两点意义**：
> 1. **S-4 的计算得到了独立的外部印证**——同一组基因，不同数据集，不同方法，原作者以自己的流程也得到同一张名单。
> 2. **原作者已经把 `Cxcr4`/`Itgam`/`Klf2` 与 `Ccl3/4/5` 放在同一句话里描述**，且把整体状态称为「with features associated with **tissue residency**」（摘要原文）。**S-4 §4 提出的「这可能是一套驻留化程序」这一读法，并非本项目独有的怀疑——原论文的措辞已经指向同一处，只是没有把它当作对肿瘤特异性的威胁来检验。**

---

## 3. Li 2022（CD8 T，MC38）

### 3.1 做过什么

光转换**仅限肿瘤**；取样为**肿瘤**与**引流淋巴结**。全文出现 35 次的 `ear` 全部指向方法学讨论——引用 Torcellan et al. 2017 曾在耳廓内的**小肿瘤**上做光转换（原文：「Previous efforts to photoconvert very small tumors within the ear pinnae of Kaede mice failed to label all cells」）。**那同样是肿瘤。**

`naive` 的 16 次出现全部指 **naive CD8 T 细胞表型**（CD44⁻CD62L⁺），**不是** naive 小鼠或非肿瘤对照。

### 3.2 判定

**`NOT_DONE`。** Li 2022 无任何非肿瘤组织的光转换对照。

---

## 4. 合并判定

| 论文 | 非肿瘤组织的**时间分辨**（光转换）对照 | 非肿瘤组织的**静态**比较 |
|---|---|---|
| Dean 2024（NK） | **`NOT_DONE`** | 有——流式表型，脾／肝／小肠／结肠（Supplementary Fig. 4） |
| Li 2022（CD8 T） | **`NOT_DONE`** | **`NOT_DONE`** |

> ### §4 的那句判断成立，且比预想的更具体
> 整个领域在这两篇代表作里**都未做过**该对照。**「NK 入瘤后功能丧失」这一解读因此始终未与「NK 进入任何组织后的驻留化」区分开。**
>
> 更值得注意的是：**Dean 2024 的摘要自己就写着瘤内 NK「adopt a distinct phenotypic state with features associated with tissue residency」。** 也就是说，「驻留化」这一描述**原文已经给出**，但它被当作**肿瘤内发生的事**来陈述，而不是被当作**任何组织都会发生、因而需要对照排除**的替代解释。
>
> **这是一个关于文献的结论，独立于我们的数据能否支持它**（§4 原话）。它成立与否不取决于 S-5 §2 的结果。

---

## 5. 这条结论**不**主张什么

- **不**主张 Dean 2024 或 Li 2022 有错误。两篇论文都没有声称做过该对照；本节记录的是一处**领域层面的空缺**，不是对具体论文的指控。
- **不**主张「驻留化」就是全部解释。那需要 §2 的数据来判，两者独立。
- **不**主张 Withers 实验室没有非肿瘤 Kaede 工作——**他们有**（`PMID 31152090`，2019，臂丛淋巴结的迁移性与驻留性 ILC；`PMID 25575242`、`PMID 28295233`、`PMID 33414524` 等）。缺的是**把非肿瘤组织的同型对照放进肿瘤论文里做对照**这一步，以及那些非肿瘤工作的**公开转录组数据**（`NONTUMOR_KAEDE_SEARCH.md` §3 记 `NOT_AVAILABLE`）。

---

## 6. 证据等级

| 断言 | 等级 |
|---|---|
| 两篇论文的全文内容与引文 | **[V]** —— 全文 XML 实读 |
| 「未做过时间分辨的非肿瘤对照」 | **[V]** —— 由全文的穷尽关键词检索与逐个上下文阅读得出；限制：**补充材料的图形内容未逐一目视核查**，仅依据正文与图注 |
| 「Dean 2024 的基因名单与 S-4 的残差顶端吻合」 | **[V]**（原文措辞）+ **[M]**（S-4 实测） |
| 「这是领域层面的空缺」 | 推论 —— 基于两篇代表作，**未穷尽全部肿瘤 Kaede 文献**（Withers 群 15 篇中另有 DC、巨噬细胞、γδ T 等，未逐篇核查），故记为**指向性结论**，不作为完备性断言 |
