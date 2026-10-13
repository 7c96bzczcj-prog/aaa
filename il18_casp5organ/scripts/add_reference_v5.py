"""Reference layer v5: mouse expansion, including a real boundary on the project's own "确定" result.

Three things happened in this batch.

(1) Mouse OLIGODENDROCYTE got much stronger. Three methods that share nothing agree: Tabula Muris Senis
    Smart-seq2 (2094 sorted cells, Casp1 in 18 of them), and a Plp-Cre RiboTag series whose logic is the
    cleanest available — the oligodendrocyte IP has LOWER Casp1 than its own paired whole-tissue input
    (0.82 vs 2.36 CPM), i.e. Casp1 is depleted along with the myeloid markers while Il18 is enriched 2.6x.
    Contamination cannot produce depletion.

(2) Mouse ASTROCYTE turned out to be REGION-SPECIFIC, and the screen's claim has to be narrowed. Cortex
    and hippocampus hold up (8 RiboTag libraries, 0.13-0.99 CPM). Cerebellum (Bergmann glia) and spinal
    cord do not: 8.57-18.89 CPM across three independent sources. Those are added under their OWN screen
    labels, not under "brain astrocyte", because the screen's mouse astrocytes come from a cortex +
    hippocampus dataset and folding a cerebellar measurement into that row would refute a claim nobody
    made. The decisive argument that this is not contamination: the spinal-cord astrocyte IP is 6.8x
    ENRICHED for Casp1 over the whole tissue that contains all of its microglia — carryover can only
    pull a sorted population toward the tissue value, never past it.

(3) Allen trimmed-mean rows were again offered and are again NOT taken, for the reason recorded in
    reference batch v4: in the mouse 10x file Micro-PVM reads Casp1 = 0.000, so the metric has no
    sensitivity for this gene in the one cell type that certainly expresses it.

A calibration worth keeping in view: resting mouse microglia carry only ~50-61 CPM Casp1 (17% of cells).
The dynamic range between "myeloid-high" and "floor" in mouse CNS is about 50x, not 1000x, so any claim
of a Casp1 zero has to be read against a microglia value measured by the same assay. The anchor rows
below exist for that purpose.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ext" / "unified_reference.tsv"

TMS = "TabulaMurisSenis-FACS-SmartSeq2"
TMS_URL = "https://czb-tabula-muris-senis.s3.amazonaws.com/Data-objects/"
G = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

# (cell_type_screen, matched_population, source, url, unit, state, quality, Il18, Casp1)
NEW = [
    # ---- mouse astrocyte: cortex / hippocampus — region-matched to the screen's dataset -----------
    ("brain astrocyte",
     "astrocyte, Brain_Non-Myeloid, FACS-sorted Smart-seq2 || 592 cells / 12 mice "
     "(Il18 detected in 59.3% of cells, Casp1 in 1.35%; still 0% at 18m and 24m)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice, no stimulation)",
     "exact", 184.3, 0.3),
    ("brain astrocyte",
     "astrocyte RiboTag (anti-HA) IP, CEREBRAL CORTEX || n=4 healthy mice "
     "(Slc1a3 1986-3505 CPM; per-animal Casp1 0.657/0.789/0.474/0.133)",
     "GEO-GSE100329-RiboTag", G + "GSE100nnn/GSE100329/suppl/GSE100329_rsem_count3.txt.gz",
     "CPM (recomputed here from the authors' RSEM counts)", "steady-state (healthy, uninjured)",
     "exact", 47.69, 0.51),
    ("brain astrocyte",
     "astrocyte RiboTag (anti-HA) IP, HIPPOCAMPUS || n=4 healthy mice",
     "GEO-GSE100329-RiboTag", G + "GSE100nnn/GSE100329/suppl/GSE100329_rsem_count3.txt.gz",
     "CPM (recomputed here from the authors' RSEM counts)", "steady-state (healthy, uninjured)",
     "exact", 24.28, 0.67),

    # ---- mouse astrocyte: the regions where the claim does NOT hold. Separate labels on purpose. ---
    ("cerebellum Bergmann glia (astrocyte) — REGIONAL COUNTEREXAMPLE",
     "Bergmann glia, FACS-sorted Smart-seq2 || 54 cells / 9 mice "
     "(Casp1 detected in 22.2% of cells — higher than sorted microglia at 16.7%; Csf1r 0%)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice)",
     "exact", 62.9, 22.4),
    ("cerebellum astrocyte — REGIONAL COUNTEREXAMPLE",
     "astrocyte RiboTag IP, CEREBELLUM || n=4 healthy mice (per-animal Casp1 9.26/18.22/18.03/8.57; "
     "myeloid markers the lowest in the study: Ptprc 0.01, C1qa 0.69)",
     "GEO-GSE100329-RiboTag", G + "GSE100nnn/GSE100329/suppl/GSE100329_rsem_count3.txt.gz",
     "CPM (recomputed here from the authors' RSEM counts)", "steady-state (healthy, uninjured)",
     "exact", 13.34, 13.52),
    ("spinal cord astrocyte — REGIONAL COUNTEREXAMPLE",
     "astrocyte RiboTag IP, SPINAL CORD || n=4 healthy mice. Casp1 12.45 CPM in the IP vs only 1.82 CPM "
     "in the paired whole tissue — 6.8x ENRICHED, which carryover cannot produce",
     "GEO-GSE100329-RiboTag", G + "GSE100nnn/GSE100329/suppl/GSE100329_rsem_count3.txt.gz",
     "CPM (recomputed here from the authors' RSEM counts)", "steady-state (healthy, uninjured)",
     "exact", 32.15, 12.45),
    ("spinal cord astrocyte — REGIONAL COUNTEREXAMPLE",
     "Aldh1l1-Cre;RiboTag cervical spinal cord || n=3 uninjured controls (Gfap 9350, Aqp4 4130 CPM; "
     "myeloid fraction ~0.16% so pure carryover would give 0.08 CPM — observed is 100-210x that)",
     "GEO-GSE277064-RiboTag", G + "GSE277nnn/GSE277064/suppl/GSE277064_Aldh1l1.txt.gz",
     "CPM (recomputed here from the authors' raw counts)", "steady-state (uninjured control)",
     "exact", 77.25, 16.67),

    # ---- mouse oligodendrocyte --------------------------------------------------------------------
    ("brain oligodendrocyte",
     "oligodendrocyte, Brain_Non-Myeloid, FACS-sorted Smart-seq2 || 2094 cells / 12 mice "
     "(Casp1 detected in 18 cells = 0.86%, against 16.7% in sorted microglia by the same assay)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice)",
     "exact", 73.7, 0.03),
    ("brain oligodendrocyte",
     "Plp-Cre;Rpl22-HA RiboTag IP, spinal cord, naive || n=3. Casp1 LOWER in the IP (0.82) than in its "
     "own paired whole-tissue input (2.36 CPM) while Il18 is enriched 2.6x — depletion, not carryover",
     "GEO-GSE225308-RiboTag", G + "GSE225nnn/GSE225308/suppl/GSE225308_GeneMatrix_RawCounts.txt.gz",
     "CPM (recomputed here from the authors' raw counts)", "steady-state (naive)",
     "close", 39.26, 0.82),

    # ---- mouse islet delta / ductal ---------------------------------------------------------------
    ("islet delta cell",
     "FACS-sorted mouse islet delta cells || n=6 libraries (Sst 163268 FPKM; Casp1 >1 in 0 of 6)",
     "GEO-GSE76017", G + "GSE76nnn/GSE76017/suppl/GSE76017_data_geo.txt.gz",
     "authors' per-gene expression table (FPKM)", "steady-state (healthy mice)", "exact", 4.71, 0.0),
    ("islet delta cell",
     "purified mouse delta cells, young/adult/old || n=9 libraries, 5.7-18.0 M counts each — "
     "ZERO Casp1 raw reads in all nine",
     "GEO-GSE269616", G + "GSE269nnn/GSE269616/suppl/GSE269616_raw_counts.txt.gz",
     "CPM (recomputed here from raw counts)", "steady-state (healthy mice)", "exact", 6.77, 0.0),
    ("islet delta cell",
     "delta cell, Pancreas, FACS-sorted Smart-seq2 || 191 cells / 8 mice (Sst in 99.5% of cells)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Pancreas.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice)",
     "exact", 12.1, 0.1),
    ("pancreas ductal cell",
     "pancreatic ductal cell, FACS-sorted Smart-seq2 || 313 cells / 8 mice (Krt19 in 93% of cells)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Pancreas.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice)",
     "exact", 12.7, 0.1),

    # ---- mouse alveolar macrophage: two more platforms confirming the "opposite" call --------------
    ("lung alveolar macrophage",
     "MF.Lu, double-sorted alveolar macrophage, Affymetrix GPL6246 || n=3 arrays "
     "(Siglecf 2442, Marco 590, Car4 2267; array floor 11.42)",
     "ImmGen-Microarray-GSE15907", G + "GSE15nnn/GSE15907/suppl/GSE15907_series_matrix.txt.gz",
     "linear RMA signal", "steady-state (unstimulated)", "exact", 1607.82, 489.81),
    ("lung alveolar macrophage",
     "FACS-sorted alveolar macrophage (CD45+CD11c+F4/80+), WT arm || n=3 mice",
     "GEO-GSE155965", G + "GSE155nnn/GSE155965/suppl/", "authors' FPKM table",
     "steady-state (wild-type control arm)", "exact", 256.00, 40.41),

    # ---- anchors: what Casp1 actually looks like in the myeloid cells of the same assays -----------
    ("brain microglia",
     "[ANCHOR] microglia, Brain_Myeloid, FACS-sorted Smart-seq2 || 13130 cells / 14 mice "
     "(Casp1 in 16.7% of cells; C1qa 8452 CPM)",
     TMS, TMS_URL + "tabula-muris-senis-facs-processed-official-annotations-Brain_Myeloid.h5ad",
     "mean CPM (recomputed here from raw integer counts)", "steady-state (healthy mice)",
     "exact", 42.0, 61.0),
    ("brain microglia",
     "[ANCHOR] Cx3cr1-CreER;RiboTag microglia || n=3 uninjured controls "
     "(Csf1r 8228, P2ry12 2375 CPM; carries the full inflammasome module: Nlrp3 19.6, Casp4 9.8, Il1b 5.0)",
     "GEO-GSE277064-RiboTag", G + "GSE277nnn/GSE277064/suppl/",
     "CPM (recomputed here from the authors' raw counts)", "steady-state (uninjured control)",
     "exact", 24.61, 49.97),
]

# NOT added, deliberately:
#   Allen mouse CTX+HPF trimmed means (SMART-seq and 10x) — same zero-inflated metric retired in v4;
#     the 10x file reports Casp1 = 0.000 for Micro-PVM, which settles it.
#   GSE312761 cortical astrocyte TRAP — five myeloid markers all score ABOVE Casp1, so its Casp1 cannot
#     be attributed to astrocytes. Recorded as unjudgeable rather than as support.
#   Mouse cortical PV / SST / VIP / L5 ET subclass rows — they change no verdict (no mouse neuron
#     subclass is a candidate), and the only usable values came from the retired trimmed-mean metric.
#     The finding they carry is a conclusion, not a reference row, and is written into the run log:
#     mouse cortical PV and L5 ET do sit at the Casp1 floor, but their Il18 is NOT high (median 0.000
#     and 0.442 against 8.302 for astrocytes), so the human PV / L5 ET result has no mouse counterpart.


def main():
    ref = pd.read_csv(REF, sep="\t")

    def num(v):
        try:
            return round(float(v), 6)
        except (TypeError, ValueError):
            return None
    have = {(r.source, r.cell_type_screen, r.gene, num(r.value)) for r in ref.itertuples()}
    rows = []
    for screen, pop, src, url, unit, state, qual, il18, casp1 in NEW:
        for gene, val in (("IL18", il18), ("CASP1", casp1)):
            if (src, screen, gene, num(val)) in have:
                continue
            rows.append({"species": "mouse", "cell_type_screen": screen, "matched_population": pop,
                         "source": src, "url": url, "gene": gene, "value": val, "unit": unit,
                         "state": state, "match_quality": qual})
    out = pd.concat([ref, pd.DataFrame(rows)], ignore_index=True) if rows else ref
    out.to_csv(REF, sep="\t", index=False)
    print(f"appended {len(rows)} mouse rows -> {len(out)} total")
    for s in sorted({r["cell_type_screen"] for r in rows}):
        srcs = sorted({r["source"] for r in rows if r["cell_type_screen"] == s})
        print(f"  {s[:52]:52} {', '.join(srcs)}")


if __name__ == "__main__":
    main()
