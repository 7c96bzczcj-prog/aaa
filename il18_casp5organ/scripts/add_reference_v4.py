"""Reference layer v4: human cortical neuron SUBCLASS evidence, plus retirement of an unusable source.

Two changes, both consequential.

(1) RETIRE the Allen "trimmed means" rows. Allen's trimmed mean is defined as the mean of the middle 50%
    of log2(CPM) after dropping the top and bottom 25% of cells, so any gene detected in fewer than about
    half the cells of a cluster is reported as exactly 0. The proof that this breaks the present question
    is inside our own table: those rows give CASP1 = 0.000 for human CNS MICROGLIA, the one cell type that
    unambiguously expresses CASP1 (every other source here puts microglia among the highest). A metric
    that cannot see CASP1 in microglia cannot be used as evidence that CASP1 is absent anywhere else.
    These rows are re-flagged so the ratio code (which uses exact/close only) skips them. They are kept in
    the file rather than deleted, so the record of what was used and why it was dropped stays visible.
    Impact check: mouse brain astrocyte/oligodendrocyte keep two independent supporting sources
    (CELLxGENE WMG and brainrnaseq/Zhang) after this removal, so that verdict does not depend on it.

(2) ADD subclass-level references for human cortical PV interneurons and L5 ET neurons. The screen found
    these in four independent datasets, but the reference layer only resolved human CNS to "neuron", a
    class that turns out to be heterogeneous for CASP1 — which is exactly why it looked internally
    inconsistent. Values below were computed in-session from the named files.

Numbers are per-source as noted in `unit`. The two Pressl GSE accessions share donors and are one study,
so they carry one source name and count once.
"""
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ext" / "unified_reference.tsv"

PRESSL = "GEO-GSE233386_GSE233408-Pressl2024-sFANS"
PRESSL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233386/suppl/"
ALLEN_SS = "AllenBrainMap-HumanCortex-SMARTseq-pseudobulk"
ALLEN_URL = "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv"
LCM = "GEO-GSE149154-LCM"
LCM_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4491nnn/"
NEUN = "GEO-GSE287651-NeuN-sorted-nuclei"
NEUN_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287651/suppl/GSE287651_rpkm_NeuN.csv.gz"
HPA_SC = "HumanProteinAtlas-SingleCellConsensus-v25"
HPA_URL = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"

# (cell_type_screen, matched_population, source, url, unit, state, quality, IL18, CASP1)
NEW = [
    # ---- PV interneurons -------------------------------------------------------------------------
    ("CNS PV interneuron",
     "Interneuron-PVALB, sFANS antibody-sorted nuclei, bulk RNA-seq || n=3 libraries / 2 donors "
     "(PVALB 165.8 CPM, GAD1 1105; microglia carryover ~0.61% so pure-contamination CASP1 would be 0.68)",
     PRESSL, PRESSL_URL, "CPM (computed here from per-sample raw counts)",
     "steady-state (postmortem CTRL donors, no stimulation)", "exact", 26.02, 0.046),
    ("CNS PV interneuron",
     "PVALB subclass pseudobulk, SMART-seq single nuclei || n=2800 nuclei / 3 donors "
     "(CASP1 detected in 10/2800 nuclei; 24021 vs 67 summed reads)",
     ALLEN_SS, ALLEN_URL, "mean raw reads per nucleus (exon+intron)",
     "steady-state (normal donors, no stimulation)", "exact", 8.579, 0.024),
    ("CNS PV interneuron",
     "aggrecan+ (PV) interneurons, laser-capture microdissection, 200 cells pooled || n=4 donors "
     "(PVALB 20.9x enriched; PTPRC/AIF1/CD68/CSF1R all 0 raw counts)",
     LCM, LCM_URL, "FPKM", "steady-state (healthy postmortem BA9)", "close", 0.0, 0.0),
    ("CNS PV interneuron", "brain inhibitory neurons || consensus over single-cell datasets",
     HPA_SC, HPA_URL, "nCPM", "steady-state (healthy donor tissue)", "close", 5.6, 0.0),

    # ---- L5 ET / corticospinal projection neurons ------------------------------------------------
    ("CNS L5 ET neuron",
     "Betz cells (BA4/BA9), sFANS antibody-sorted nuclei, bulk RNA-seq || n=5 CTRL libraries / 4 donors "
     "(FEZF2 132.7, BCL11B 104.8, NEFH 223.4 CPM; microglia carryover ~1.03%)",
     PRESSL, PRESSL_URL, "CPM (computed here from per-sample raw counts)",
     "steady-state (postmortem CTRL donors, no stimulation)", "exact", 17.90, 0.065),
    ("CNS L5 ET neuron",
     "von Economo neurons, sFANS antibody-sorted nuclei || n=9 CTRL libraries",
     PRESSL, PRESSL_URL, "CPM (computed here from per-sample raw counts)",
     "steady-state (postmortem CTRL donors, no stimulation)", "close", 18.04, 0.047),
    ("CNS L5 ET neuron",
     "L5 ET subclass pseudobulk, SMART-seq single nuclei || n=158 nuclei / 3 donors "
     "(BCL11B detected in 100% of nuclei; CASP1 detected in 1/158; 851 vs 4 summed reads)",
     ALLEN_SS, ALLEN_URL, "mean raw reads per nucleus (exon+intron)",
     "steady-state (normal donors, no stimulation)", "exact", 5.386, 0.025),

    # ---- positive controls: cortical neuron subclasses that DO express CASP1. Not mapped to any
    # screen row; they are here because they are the evidence that the assay can see CASP1 in neurons,
    # which is what makes the PV / L5 ET near-zero a real subclass difference rather than a detection
    # limit. They also show "CNS neuron" is not a homogeneous reference class.
    ("CNS L5a IT neuron (CASP1-positive control)",
     "L5a pyramidal, sFANS sorted nuclei || n=25 CTRL libraries, 24/25 non-zero; 7.2x its contamination expectation",
     PRESSL, PRESSL_URL, "CPM (computed here from per-sample raw counts)",
     "steady-state (postmortem CTRL donors)", "exact", 17.21, 8.34),
    ("CNS L5/6 NP neuron (CASP1-positive control)",
     "L5/6 NP subclass pseudobulk, SMART-seq || n=816 nuclei, CASP1 detected in 16.4% of nuclei, "
     "reproduced in all 3 donors",
     ALLEN_SS, ALLEN_URL, "mean raw reads per nucleus (exon+intron)",
     "steady-state (normal donors)", "exact", 4.846, 23.488),

    # ---- pan-neuron baseline ---------------------------------------------------------------------
    ("CNS neuron",
     "NeuN+ sorted nuclei, entorhinal cortex, bulk RNA-seq || n=4 donors without AD pathology",
     NEUN, NEUN_URL, "RPKM", "steady-state (asymptomatic controls)", "close", 1.821, 0.152),
]


def main():
    ref = pd.read_csv(REF, sep="\t")

    # (1) retire the trimmed-mean rows
    bad = ref.unit.astype(str).str.contains("trimmed mean", case=False, na=False)
    n_bad = int(bad.sum())
    ref.loc[bad, "match_quality"] = "unusable_trimmed_mean_zero_inflated"
    ref.loc[bad, "state"] = ref.loc[bad, "state"].astype(str) + \
        " [RETIRED: Allen trimmed mean zeroes any gene detected in <~50% of cells; it reports CASP1=0 " \
        "even for microglia, so it cannot evidence CASP1 absence]"

    # (2) append the subclass rows
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
            rows.append({"species": "human", "cell_type_screen": screen, "matched_population": pop,
                         "source": src, "url": url, "gene": gene, "value": val, "unit": unit,
                         "state": state, "match_quality": qual})
    out = pd.concat([ref, pd.DataFrame(rows)], ignore_index=True) if rows else ref
    out.to_csv(REF, sep="\t", index=False)
    print(f"retired {n_bad} trimmed-mean rows; appended {len(rows)} subclass rows -> {len(out)} total")
    for s in sorted({r["cell_type_screen"] for r in rows}):
        srcs = sorted({r["source"] for r in rows if r["cell_type_screen"] == s})
        print(f"  {s:46} sources: {', '.join(srcs)}")


if __name__ == "__main__":
    main()
