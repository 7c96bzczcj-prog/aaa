#!/usr/bin/env python3
"""Route D — re-audit the T3 dataset candidates at the VARIANT level, not the file level.

T3 classified candidates by what files exist (fragments / mgatk output present →
"chrM present"). GSE302113 then showed that is the wrong question: chrM was fully
retained, every library cleared the coverage floor, and the dataset still could not
carry the statistic, because the deposited heteroplasmy matrices are per-library
and the variant union is not deposited.

So the candidates' status is UNKNOWN, not "second best". This re-audit asks, per
dataset: is per-cell heteroplasmy deposited at all; is it per-library or per-donor
union; and does the design even have multiple compartments per donor.
"""
import json
import os
import re
import time
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = {"User-Agent": "nk-lineage-screen/1.0 (research use)"}
CANDIDATES = ["GSE302113", "GSE197037", "GSE197008", "GSE148508",
              "GSE173936", "GSE216911", "GSE276283"]

# design, from the series records — decides whether a cross-compartment question is
# even askable, independently of what is deposited
COMPARTMENTS = {
    "GSE302113": "tumour + non-involved organ + blood, 5 NSCLC donors (+5 ovarian)",
    "GSE197037": "PBMC only (CMV+/- donors) — no tissue compartment",
    "GSE197008": "PBMC only (CMV+/- donors) — no tissue compartment",
    "GSE148508": "CD34 + PBMC — no tumour, no matched normal organ",
    "GSE173936": "PBMC / mixed mtDNA-disorder samples — not a tumour design",
    "GSE216911": "cultured T cells — not a tissue design",
    "GSE276283": "glioma-associated macrophages — unretrieved",
}


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 ** i)
    return ""


def filelist(acc):
    return get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/filelist.txt")


def classify(acc, files):
    low = files.lower()
    lines = [l.split("\t")[1] for l in files.strip().split("\n")[1:] if "\t" in l]
    # mgatk output is commonly deposited as an .rds SummarizedExperiment; missing that
    # pattern produced a false negative on GSE197037 on the first run, which is the same
    # error class this whole re-audit exists to catch.
    per_cell = [f for f in lines
                if re.search(r"heteroplasm|allele|af\.tsv|cell_.*variant|mgatk", f, re.I)]
    mgatk_obj = [f for f in lines if re.search(r"mgatk.*\.(rds|h5|zip|tar)", f, re.I)]
    varstats = [f for f in lines if "variant_stats" in f.lower()]
    frags = [f for f in lines if "fragments" in f.lower()]
    bams = [f for f in lines if re.search(r"\.bam\b", f, re.I)]
    rds = [f for f in lines if re.search(r"\.rds|\.h5|\.h5ad|seurat", f, re.I)]
    ann = [f for f in lines if re.search(r"annot|metadata|celltype|cell_type|clone", f, re.I)]

    # per-library vs per-donor: count distinct GSM prefixes on the heteroplasmy files
    gsm_pref = {re.match(r"(GSM\d+)", f).group(1) for f in per_cell if re.match(r"GSM\d+", f)}
    if per_cell and len(gsm_pref) > 1:
        granularity = f"per-library ({len(gsm_pref)} GSMs carry one each)"
    elif per_cell:
        granularity = "single file (may be pooled/union — needs opening)"
    else:
        granularity = "no per-cell heteroplasmy deposited"

    return {
        "dataset_id": acc,
        "mgatk_objects": len(mgatk_obj),
        "n_files": len(lines),
        "per_cell_heteroplasmy": len(per_cell),
        "variant_stats": len(varstats),
        "fragments": len(frags),
        "bams": len(bams),
        "object_files_rds_h5": len(rds),
        "annotation_like": len(ann),
        "heteroplasmy_granularity": granularity,
        "raw_allele_counts": ("YES — mgatk object deposited, any variant evaluable in any cell"
                              if mgatk_obj else
                              "NO — only derived/selected summaries or nothing"),
    }


COLS = ["dataset_id", "n_files", "per_cell_heteroplasmy", "mgatk_objects", "variant_stats",
        "fragments", "bams", "object_files_rds_h5", "annotation_like",
        "heteroplasmy_granularity", "raw_allele_counts", "compartments", "verdict"]


def main():
    rows = []
    for acc in CANDIDATES:
        files = filelist(acc)
        if not files:
            rows.append({"dataset_id": acc, "verdict": "FILELIST UNRETRIEVABLE",
                         "union_matrix_deposited": "UNKNOWN"})
            print(f"{acc:12s} filelist unretrievable")
            continue
        r = classify(acc, files)
        r["compartments"] = COMPARTMENTS.get(acc, "unknown")
        if r["mgatk_objects"]:
            r["verdict"] = ("NO ascertainment confound — raw mgatk allele counts deposited, so any "
                            "variant can be evaluated in any cell. Usable for a cross-compartment "
                            "statistic IF the design has multiple compartments per donor")
        elif r["per_cell_heteroplasmy"]:
            r["verdict"] = ("derived/selected summaries only — SAME FAILURE MODE AS GSE302113 unless "
                            "opening shows otherwise")
        elif r["fragments"] and not r["bams"]:
            r["verdict"] = ("fragments only. Fragments carry intervals, NOT alleles, so mtDNA "
                            "variants cannot be recovered from them at all — reprocessing needs raw reads")
        else:
            r["verdict"] = "needs opening at variant level"
        rows.append(r)
        print(f"{acc:12s} het={r['per_cell_heteroplasmy']:3d} varstats={r['variant_stats']:3d} "
              f"frag={r['fragments']:3d} bam={r['bams']:3d} | {r['heteroplasmy_granularity']}")

    out = os.path.join(HERE, "out", "T3b_datasets_variant_level.tsv")
    with open(out, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(" ".join(str(r.get(c, "NA")).split()) for c in COLS) + "\n")
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
