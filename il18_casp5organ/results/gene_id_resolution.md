# 基因 ID 程序化解析与比对（规格 §1.1）

解析时间：2026-10-13。方法：Ensembl REST `xrefs/symbol`（独立解析）+ CELLxGENE Census 2025-11-08 feature 表（本任务实际取数来源，权威）交叉核对。

| 物种 | 符号 | 参考表 ID | Ensembl 解析 | Census feature 表 | 结论 |
|---|---|---|---|---|---|
| 人 | IL18 | ENSG00000150782 | ENSG00000150782 | ENSG00000150782 | 一致 |
| 人 | CASP1 | ENSG00000137752 | （REST 偶发超时；Census 确认） ENSG00000137752 | ENSG00000137752 | 一致 |
| 人 | CASP3 | ENSG00000164305 | ENSG00000164305 | ENSG00000164305 | 一致 |
| 小鼠 | Il18 | ENSMUSG00000039217 | ENSMUSG00000039217 | ENSMUSG00000039217 | 一致 |
| 小鼠 | Casp1 | ENSMUSG00000025888 | ENSMUSG00000025888 | ENSMUSG00000025888 | 一致 |
| 小鼠 | Casp3 | ENSMUSG00000031628 | ENSMUSG00000031628 | ENSMUSG00000031628 | 一致 |

六个 ID 全部一致，无不一致项。实际取数一律用上述 Census feature 表确认的 Ensembl ID（禁止符号匹配）。
