"""Spec §9.5: casp1_zero_check.tsv. Pair each candidate's CASP1 detection with sorted-bulk CASP1
(DICE human immune / ImmGen mouse immune / brainrnaseq CNS) to judge true-low vs dropout. Presentation
only, no statistical test. Hepatocyte/acinar/islet etc. have no sorted reference -> no_sorted_reference.
"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
from common import write_tsv
ROOT = Path(__file__).resolve().parents[1]

IMMUNE = {"lymphoid", "dendritic_cell", "monocyte_derived_myeloid", "granulocyte"}
# representative sorted-bulk CASP1 (within-source units; see reference/*.tsv for all populations)
REF = {
 ("human", "immune"): ("DICE", "blood immune (mono 301 / NK 112 / B 38 / T 37-48 TPM)", ">1 TPM 全供体", "clearly expressed"),
 ("mouse", "immune"): ("ImmGen", "sorted immune (mono 776 / neu 623 / NK 457 / B 317; 1.00=zero)", "数百计数", "clearly expressed"),
 ("human", "resident_macrophage"): ("brainrnaseq/DICE", "human microglia 16 FPKM / blood mono 301 TPM", "-", "clearly expressed"),
 ("mouse", "resident_macrophage_lung"): ("ImmGen", "alveolar MF 511 (1.00=zero)", "-", "clearly expressed"),
 ("mouse", "resident_macrophage_cns"): ("ImmGen/brainrnaseq", "microglia 216-304 / 4.2 FPKM", "-", "clearly expressed"),
 ("human", "neuron"): ("brainrnaseq", "human neuron 0.6 FPKM", "-", "truly very low"),
 ("mouse", "neuron"): ("brainrnaseq", "mouse neuron 0.64 FPKM", "-", "truly very low"),
 ("human", "glia"): ("brainrnaseq", "human astro 0.53 / oligo 1.6 FPKM", "-", "truly very low"),
 ("mouse", "glia"): ("brainrnaseq", "mouse astro 1.13 / oligo 0.1-0.12 / OPC 0.66 FPKM", "-", "truly very low"),
 ("human", "endothelial_cns"): ("brainrnaseq", "human CNS endothelial 2.25 FPKM", "-", "low"),
 ("mouse", "endothelial_cns"): ("brainrnaseq", "mouse CNS endothelial 0.1 FPKM", "-", "truly very low"),
}


def classify(species, tissue, lineage, label):
    cns = tissue in ("brain", "spinal_cord")
    micro = ("microglia" in label.lower() or "central nervous system macrophage" in label.lower()
             or (lineage == "resident_macrophage" and cns))
    if micro:
        key = (species, "resident_macrophage_cns")
        return (*REF.get(key, REF[(species, "immune")]), "dropout_likely",
                "小胶质/CNS 巨噬 CASP1 确有表达 → 单核中的零多为 dropout，CASP1 低不可信")
    if lineage == "resident_macrophage" and tissue == "lung":
        return (*REF[(species, "resident_macrophage_lung")], "dropout_likely",
                "肺泡巨噬 CASP1 bulk 强表达 → 其“CASP1 低”为相对/技术性，不可作真低")
    if lineage in IMMUNE or label in ("basophil", "mast cell", "dendritic cell"):
        return (*REF[(species, "immune")], "dropout_likely", "免疫细胞 CASP1 bulk 普遍强表达 → 单细胞零多为 dropout")
    if lineage == "neuron":
        return (*REF[(species, "neuron")], "true_low_likely", "神经元 CASP1 bulk 真低 → 单核零多为真零，CASP1 低可信")
    if lineage == "glia":
        return (*REF[(species, "glia")], "true_low_likely", "星胶/少突/OPC CASP1 bulk 真低 → CASP1 低可信（注意小胶质不在此列）")
    if lineage == "endothelial" and cns:
        return (*REF[(species, "endothelial_cns")], "true_low_likely", "CNS 内皮 CASP1 bulk 低")
    return ("none", "no_sorted_reference", "-", "-", "no_sorted_reference",
            "实质上皮/间质无同等分选对照（§9.5），该候选 CASP1 低的真伪无法判定，结论降级")


def main():
    cand = pd.read_csv(ROOT / "results" / "candidates.tsv", sep="\t", comment="#")
    cm = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#")
    rows = []
    for _, r in cand.iterrows():
        if r.level == "lineage":
            lineage = r.celltype
        else:
            mm = cm[(cm.entry_id == r.dataset_id) & (cm.celltype_original == r.celltype)]
            lineage = mm.lineage.iloc[0] if len(mm) else "other"
        src, pop, frac, refexpr, verdict, note = classify(r.species, r.tissue, lineage, str(r.celltype))
        rows.append({"dataset_id": r.dataset_id, "level": r.level, "species": r.species, "tissue": r.tissue,
                     "suspension_type": r.suspension_type, "celltype": r.celltype, "lineage": lineage,
                     "sc_pct_CASP1": r.pct_CASP1, "sc_logcpm_CASP1": r.logcpm_CASP1,
                     "reference_source": src, "ref_population_CASP1": pop, "ref_CASP1_level": refexpr,
                     "verdict": verdict, "note": note})
    out = pd.DataFrame(rows).sort_values(["species", "tissue", "suspension_type", "level"])
    write_tsv(out, ROOT / "results" / "casp1_zero_check.tsv")
    print(out.verdict.value_counts().to_dict())


if __name__ == "__main__":
    main()
