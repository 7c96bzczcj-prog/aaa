"""Shared constants for the IL18 lung analysis. Gene access is by Ensembl ID only (spec §1.1)."""
CENSUS_VERSION = "2025-11-08"  # LTS release, pinned for reproducibility
SEED = 20261007
GENES = {"homo_sapiens": "ENSG00000150782", "mus_musculus": "ENSMUSG00000039217"}
SPECIES_LABEL = {"homo_sapiens": "human", "mus_musculus": "mouse"}
OBS_COLS = ["dataset_id", "donor_id", "cell_type", "assay", "disease", "tissue",
            "suspension_type", "raw_sum"]
MIN_CELLS = 20

HEADER = (
    "# IL18 基因只有一个蛋白产物，即 24 kDa 的前体 pro-IL-18；成熟 IL-18 由半胱天冬酶（caspase）切割前体产生，"
    "不存在对应的独立转录本或剪接变体。因此本分析得到的 IL18 转录本丰度即 pro-IL-18 的转录信号，"
    "无需也无法在转录层面区分前体与成熟形式。\n"
    "# 但 pro-IL-18 属于储备型蛋白，半衰期长，转录水平与细胞内实际前体蛋白储量不成线性关系。"
    "本分析的结论限于「哪些细胞转录 IL18」，不可外推为「哪些细胞储备 pro-IL-18 蛋白」。\n"
)


def write_tsv(df, path):
    with open(path, "w") as fh:
        fh.write(HEADER)
        df.to_csv(fh, sep="\t", index=False)
