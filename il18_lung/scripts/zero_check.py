"""Spec §6.3: put near-zero single-cell cell types next to sorted-bulk IL18 (DICE human, ImmGen mouse).

Presentation only, no tests. 'Near zero' = dataset-level donor-median pct_expressing < ZERO_PCT.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import write_tsv

ROOT = Path(__file__).resolve().parents[1]
ZERO_PCT = 0.05

# explicit correspondence, unified class -> sorted reference populations (no pattern matching)
DICE = {
    "T_NK": ["T cell, CD4, naive", "T cell, CD8, naive", "T cell, CD4, TH1", "T cell, CD4, naive TREG",
             "NK cell, CD56dim CD16+"],
    "B_plasma": ["B cell, naive"],
    "other:monocyte": ["Monocyte, classical", "Monocyte, non-classical"],
}
IMMGEN = {
    "alveolar_macrophage": ["MF.Alv.Lu", "MF.alv.11cp64pSiglecFp.Lu", "MF.11cpSigFp.BAL"],
    "monocyte_derived_macrophage": ["Mo.6Cp.Lu", "Mo.64p6CpIIp.LPS.d3.Lu"],
    "dendritic_cell": ["DC.103p11bnsiglecFn.Lu", "DC.103n11bpsiglecFn.Lu", "DC.pDC.Sp"],
    "neutrophil": ["GN.BM", "GN.Sp"],
    "T_NK": ["T.4.Nve.Sp", "T.8.Nve.Sp", "Treg.4.25hi.Sp", "MAIT.Lu", "NK.27+11b+.Sp", "NK.27+11b-.Sp"],
    "B_plasma": ["B.Fo.Sp", "B.MZ.Sp"],
    "other:monocyte": ["Mo.6Cp.Lu", "Mo.6Cn.Lu"],
}
MONOCYTE_NOTES = {"monocyte"}


def main():
    s = pd.read_csv(ROOT / "results" / "il18_summary.tsv", sep="\t", comment="#")
    cmap = pd.read_csv(ROOT / "celltype_map.tsv", sep="\t", comment="#", dtype=str).fillna("")
    mono = set(cmap.loc[cmap.note.isin(MONOCYTE_NOTES), "celltype_original"])
    dice = pd.read_csv(ROOT / "reference" / "dice_il18.tsv", sep="\t").set_index("population")
    imm = pd.read_csv(ROOT / "reference" / "immgen_il18.tsv", sep="\t").drop_duplicates("population").set_index("population")
    z = s[s.pct_median < ZERO_PCT].copy()
    rows = []
    for _, r in z.iterrows():
        key = r.celltype if r.celltype != "other" else (
            "other:monocyte" if r.celltype_original in mono else None)
        ref, refmap, unit = (dice, DICE, "TPM, donor median") if r.species == "human" else (
            imm, IMMGEN, "ImmGen normalized count (1.00 = 0 reads), replicate mean")
        pops = refmap.get(key, []) if key else []
        base = dict(species=r.species, dataset_id=r.dataset_id, condition=r.condition,
                    condition_detail=r.condition_detail, celltype=r.celltype, celltype_original=r.celltype_original,
                    n_donors=r.n_donors, sc_pct_median=r.pct_median, sc_logcpm_median=r.logcpm_median)
        if not pops:
            rows.append({**base, "reference": "DICE" if r.species == "human" else "ImmGen",
                         "ref_population": "none available (no sorted bulk for this class)",
                         "ref_value": None, "ref_unit": "", "ref_n": None})
        for p in pops:
            rows.append({**base, "reference": "DICE" if r.species == "human" else "ImmGen", "ref_population": p,
                         "ref_value": ref.loc[p, "summary_value"], "ref_unit": unit,
                         "ref_n": ref.loc[p, "n_reps"]})
    out = pd.DataFrame(rows)
    write_tsv(out, ROOT / "results" / "il18_zero_check.tsv")
    print(len(z), "near-zero summary rows;", len(out), "zero-check rows")


if __name__ == "__main__":
    main()
