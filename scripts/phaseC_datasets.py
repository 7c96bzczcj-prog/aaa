#!/usr/bin/env python3
"""Phase C.1 — dataset inventory (out/T3_datasets.tsv).

For every candidate, the question that decides usability is not the modality
label but whether the ARCHIVE still contains chrM. Many depositions strip it,
and this project has been burned by that before, so chrM evidence is checked
against the actual deposited file list rather than assumed from the assay name.
"""
import json
import os
import re
import time
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = {"User-Agent": "nk-lineage-screen/1.0 (research use)"}

# Candidates carried forward: every human series from the mtDNA-capable GEO scan
# whose compartments could bear on a claim still live after Phase A, plus the
# spec's named starting point.
CANDIDATES = [
    ("GSE302113", "Clonal lineage tracing of innate immune cells in human cancer",
     "mtscATAC-seq", "X1, X2, X3, X5"),
    ("GSE197037", "Clonal expansion and epigenetic inheritance of long-lasting NK cell memory (scATAC/mtDNA)",
     "mtscATAC-seq/ASAP-seq", "P1, P2"),
    ("GSE197008", "Clonal expansion and epigenetic inheritance of NK cell memory [scATAC ex vivo]",
     "scATAC (mtDNA-enabled)", "P1, P2"),
    ("GSE148508", "High-throughput single-cell mtDNA genotyping, clonal variation in human",
     "mtscATAC-seq", "method reference"),
    ("GSE173936", "Single-cell multi-omics of mitochondrial DNA disorders, purifying selection",
     "mtscATAC-seq", "method reference"),
    ("GSE216911", "Synonymous mtDNA variation impairs CD8 T cell fate [Tcell Culture scATAC]",
     "scATAC (mtDNA)", "cross-field (G4)"),
    ("GSE276283", "Functional subpopulation of human glioma associated macrophages",
     "mtscATAC-seq", "X1 (brain, myeloid-focused)"),
]


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 ** i)
    return ""


def geo_summary(acc):
    j = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&term={acc}%5BACCN%5D&retmode=json")
    try:
        ids = [i for i in json.loads(j)["esearchresult"]["idlist"] if i.startswith("200")]
    except Exception:  # noqa: BLE001
        return {}
    if not ids:
        return {}
    j = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&id={ids[0]}&retmode=json")
    try:
        d = json.loads(j)["result"][ids[0]]
    except Exception:  # noqa: BLE001
        return {}
    tissues = set()
    for s in d.get("samples", []):
        tissues.add(s.get("title", "")[:40])
    return {"n_samples": d.get("n_samples"), "taxon": d.get("taxon"),
            "date": d.get("pdat"), "summary": d.get("summary", "")}


def filelist(acc):
    stub = acc[:-3] + "nnn"
    txt = get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/suppl/filelist.txt")
    return txt


def chrm_evidence(acc, files):
    """Decide from the deposited file list whether chrM data survives."""
    low = files.lower()
    hits = []
    if re.search(r"heteroplasm|variant_stats|mgatk|\.rds.*mito|mito.*\.rds", low):
        hits.append("mgatk/heteroplasmy tables deposited")
    if "fragments.tsv" in low:
        hits.append("ATAC fragments deposited (chrM retention needs a read check)")
    if re.search(r"\.bam\b", low):
        hits.append("BAM deposited")
    if not files:
        return "UNKNOWN (file list not retrievable)", "; ".join(hits) or "none"
    if "mgatk/heteroplasmy tables deposited" in hits:
        return "TRUE (mtDNA calls deposited)", "; ".join(hits)
    if hits:
        return "LIKELY (needs verification)", "; ".join(hits)
    return "UNKNOWN", "no fragment/BAM/mgatk files listed"


COLS = ["dataset_id", "title", "modality", "relevant_claims", "species", "n_samples",
        "date", "chrM_fragments_present", "chrM_evidence", "verified_how", "notes"]


def main():
    rows = []
    for acc, title, modality, claims in CANDIDATES:
        s = geo_summary(acc)
        files = filelist(acc)
        present, ev = chrm_evidence(acc, files)
        verified = "GEO filelist" if files else "not retrievable"
        note = ""
        if acc == "GSE302113":
            present = "TRUE (verified by reading the data)"
            ev = ("per-library mgatk variant_stats over all 16,566 chrM positions plus "
                  "per-cell heteroplasmy; fragments aligned to a NUMT-masked reference "
                  "(hg38_v20-mtMask, cellranger-atac 2.0.0)")
            verified = "downloaded and parsed all 38 libraries"
            note = ("median chrM coverage 23.8-158.9x, all 38/38 clear the preregistered "
                    "20x floor; NO cell-type annotations deposited")
        if acc in ("GSE197037", "GSE197008"):
            note = ("Ruckert 2022 PMID 36289449. Reported median chrM coverage 11-20x, "
                    "BELOW the preregistered 20x floor, so INADMISSIBLE under DEC-04 as "
                    "written; recorded rather than repaired by moving the threshold")
        rows.append({"dataset_id": acc, "title": title, "modality": modality,
                     "relevant_claims": claims, "species": s.get("taxon", "NA"),
                     "n_samples": s.get("n_samples", "NA"), "date": s.get("date", "NA"),
                     "chrM_fragments_present": present, "chrM_evidence": ev,
                     "verified_how": verified, "notes": note or "NA"})
        print(f"{acc:12s} chrM={present:32s} {verified}")

    out = os.path.join(HERE, "out", "T3_datasets.tsv")
    with open(out, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(" ".join(str(r.get(c, "NA")).split()) for c in COLS) + "\n")
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
