"""Append the expansion-batch reference sources to reference/ext/unified_reference.tsv.

Every value below was computed in-session from the named GEO supplementary file (see the run log for the
per-source verification: identity markers, contamination markers, and how the group was assigned). Nothing
here is written from memory. Sources judged unusable as a reference are deliberately NOT added, and the
reason is recorded at the bottom of this file.

Idempotent: rows whose (source, cell_type_screen, gene, value) already exist are not appended again.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ext" / "unified_reference.tsv"
G = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

# (species, cell_type_screen, matched_population, source, url, unit, state, match_quality, IL18, CASP1)
NEW = [
    # --- human liver cholangiocyte: primary biliary epithelial cells, 5 donors, untreated control.
    # PTPRC 0.01 / AIF1 0.03 CPM, so the signal cannot be macrophage carryover. Cultured, hence the
    # state string; identity KRT7 5604, KRT19 2206, CFTR 31 CPM.
    ("human", "liver cholangiocyte",
     "primary human biliary epithelial cell (untreated scrambled-siRNA control) || n=5 donors",
     "GEO-GSE317449", G + "GSE317nnn/GSE317449/suppl/GSE317449_JK_T6_RNA_counts_names.txt.gz",
     "CPM (computed here from raw counts)",
     "CULTURED primary cells (10% human serum + hydrocortisone/EGF/HGF/insulin), untreated control",
     "exact", 109.05, 20.52),

    # --- human islet delta cell: LY6H+ FACS-sorted, 2 non-diabetic donors. SST 2580/55825 CPM confirms
    # delta enrichment. PTPRC 30.0/27.3 CPM = leukocyte carryover, and CASP1 is myeloid-high, so part of
    # the CASP1 signal is probably not delta-cell. Recorded as close, not exact, for that reason.
    ("human", "islet delta cell",
     "LY6H+ FACS-sorted delta-enriched islet cells || n=2 donors (SST 2580/55825 CPM; PTPRC 30/27 CPM carryover)",
     "GEO-GSE264633", G + "GSE264nnn/GSE264633/suppl/GSE264633_LY6H_combined_raw_count_matrix.txt.gz",
     "CPM (computed here from raw counts)",
     "steady-state (FACS-sorted ex vivo, non-diabetic donor islets); leukocyte carryover present",
     "close", 5.05, 13.40),

    # --- human cortical astrocyte, acutely purified (HepaCAM immunopanning), two independent series.
    ("human", "CNS astrocyte",
     "acutely purified human cortical astrocyte (HepaCAM immunopanning) || n=4 (SLC1A3 19.0, GFAP 28.1 vs PTPRC 0.03 RPKM)",
     "GEO-GSE147870", G + "GSE147nnn/GSE147870/suppl/GSE147870_RPKM.xlsx",
     "RPKM", "steady-state (acutely purified, no culture)", "exact", 0.029, 1.112),
    ("human", "CNS astrocyte",
     "acutely purified human cortical astrocyte, Control cohort || n=12 (SLC1A3 1577, AQP4 610 vs CD68 0.05 CPM)",
     "GEO-GSE168375", G + "GSE168nnn/GSE168375/suppl/GSE168375_RAW.tar",
     "CPM (computed here from raw counts)",
     "steady-state (acutely purified from surgical resection, Control samples only)",
     "exact", 0.32, 0.10),

    # --- human lung mast cell: CD117/MACS-enriched primary HLMC, 4 donors, unstimulated.
    # TPSAB1 5469-6629, CPA3 3655-4316 CPM vs AIF1 4-35 CPM. Cultured in SCF, hence the state string.
    ("human", "lung mast cell",
     "primary human lung mast cell (CD117 MACS-enriched), unstimulated || n=4 donors",
     "GEO-GSE314641", G + "GSE314nnn/GSE314641/suppl/GSE314641_salmon.merged.gene_counts_length_scaled.tsv.gz",
     "CPM (computed here from length-scaled gene counts)",
     "CULTURED (StemPro-34 + 100 ng/ml SCF), unstimulated control",
     "exact", 85.5, 19.1),

    # --- human lung AT1: there is no sorted NATIVE human AT1 bulk anywhere I could find. This is primary
    # AT2 transdifferentiated to an AT1-like state over 6 days in culture (AGER 12.7->100.3, HOPX
    # 1.25->4.39 confirm the transition). Units are gene-body mean coverage, not TPM, so only the
    # within-source ratio is meaningful. Marked close + CULTURED so it is read as a weak proxy.
    ("human", "lung AT1",
     "AT1-like cells from primary human AT2 transdifferentiated 6 d in vitro || n=2 donors (AGER 12.7->100.3)",
     "GEO-GSE140071", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4153nnn/GSM4153673/suppl/GSM4153673_Donor1_D0_RNA.bw",
     "gene-body mean coverage (normalised per 1e9 bigWig sumData); NOT TPM — within-source ratio only",
     "CULTURED (in-vitro AT2->AT1 transdifferentiation, day 6) — not native AT1",
     "close", 2.622, 2.345),
]

# Deliberately NOT added, and why:
#   GSE247159 (human A2B5+ oligodendrocytes) — 2 days in culture, GFAP 701-4868 CPM astrocyte
#     contamination, AIF1/CD68/PTPRC high, and the A2B5- fraction has HIGHER CASP1 than A2B5+, i.e. the
#     CASP1 signal tracks non-oligodendrocyte cells. Unusable as a reference for oligodendrocyte CASP1.
#   GSE36634 — microarray, CEL files only, not quantified in-session.
#   HPA IHC "IL18 not detected in alveolar cells type I" — an ordinal protein call, not a number on the
#     same scale as the RNA rows, so it is reported as a caveat in the report rather than as a ratio row.


def main():
    ref = pd.read_csv(REF, sep="\t")
    def num(v):
        try:
            return round(float(v), 6)
        except (TypeError, ValueError):
            return None          # ordinal rows ("Not detected", "Medium", ...) never collide with ours
    have = {(r.source, r.cell_type_screen, r.gene, num(r.value)) for r in ref.itertuples()}
    rows = []
    for sp, screen, pop, src, url, unit, state, qual, il18, casp1 in NEW:
        for gene, val in (("IL18", il18), ("CASP1", casp1)):
            if (src, screen, gene, round(float(val), 6)) in have:
                continue
            rows.append({"species": sp, "cell_type_screen": screen, "matched_population": pop,
                         "source": src, "url": url, "gene": gene, "value": val, "unit": unit,
                         "state": state, "match_quality": qual})
    if not rows:
        print("nothing to add (already present)")
        return
    out = pd.concat([ref, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(REF, sep="\t", index=False)
    print(f"appended {len(rows)} rows -> {len(out)} total")
    for r in rows:
        print(f"  {r['species']:6} {r['cell_type_screen']:22} {r['source']:16} {r['gene']:6} {r['value']}")


if __name__ == "__main__":
    main()
