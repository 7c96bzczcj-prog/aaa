"""Shared constants for the five-organ IL18/CASP1/CASP3 screen (spec v3.0).

Genes are taken by Ensembl ID only (§1.2). IDs below were resolved via Ensembl REST and
cross-checked against the Census feature table — all six consistent (see results/gene_id_resolution.md).
"""
CENSUS_VERSION = "2025-11-08"   # LTS, pinned for reproducibility
SEED = 20261013

GENES = {
    "homo_sapiens": {"IL18": "ENSG00000150782", "CASP1": "ENSG00000137752", "CASP3": "ENSG00000164305"},
    "mus_musculus": {"IL18": "ENSMUSG00000039217", "CASP1": "ENSMUSG00000025888", "CASP3": "ENSMUSG00000031628"},
}
GENE_ORDER = ["IL18", "CASP1", "CASP3"]
SPECIES_LABEL = {"homo_sapiens": "human", "mus_musculus": "mouse"}
TISSUES = ["lung", "liver", "pancreas", "brain", "spinal cord"]
MIN_CELLS = 20

# §0.2 — must appear verbatim in README.md; short tag for file headers / figure captions.
POSITIONING = (
    "# 本分析是底物充足且经典加工酶缺失的交集筛选（IL18 高 + CASP1 低），用途是为后续蛋白与功能验证挑选"
    "器官和细胞类型的优先级，不是任何机制结论的证据。CASP3 仅用于排除最低档与排序，不参与候选集定义。\n"
    "# 转录本丰度 ≠ 前体蛋白储量；CASP3 mRNA ≠ caspase-3 酶活性；本分析为静息基线，不反映刺激后 CASP1 诱导。\n"
    "# 脑/脊髓以单核测序为主、其余器官以单细胞为主，snRNA 与 scRNA 不可混比；跨器官只比候选名单，不比表达量。\n"
)


def write_tsv(df, path):
    with open(path, "w") as fh:
        fh.write(POSITIONING)
        df.to_csv(fh, sep="\t", index=False)
