"""Task P-2: fetch the six sorted-NK bulk expression matrices from GEO.

Only processed matrices are pulled -- counts or normalised values as the
depositors left them.  No FASTQ, no realignment: §3.1 forbids comparing
absolute values across datasets anyway, so re-quantifying them onto a common
pipeline would buy nothing the analysis is allowed to use.

Total download is about 15 MB.  Files land in data/p2/, which is git-ignored;
the derived results in results/p2_*.csv are what gets committed.

Run: python3 scripts/p2_fetch.py
"""

from __future__ import annotations

import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "data", "p2")

SERIES = [
    # group 1
    ("GSE236nnn/GSE236394", "GSE236394_all.gene_moderated_log2cpm.tsv.gz"),
    # groups 1 and 3
    ("GSE133nnn/GSE133383", "GSE133383_filtered_NKcounts_table.csv.gz"),
    # group 2
    ("GSE140nnn/GSE140035", "GSE140035_RawCounts_RNA_CytoStims_human.txt.gz"),
    ("GSE140nnn/GSE140035", "GSE140035_RawCounts_RNA_CytoStims_mouse.txt.gz"),
    # group 3
    ("GSE200nnn/GSE200319", "GSE200319_ProcessedGeneList.csv.gz"),
    # group 4
    ("GSE205nnn/GSE205492", "GSE205492_rnaseq_workshop_normalized_counts.txt.gz"),
]

# GSE242941's bulk arm is deposited per sample; the series tarball is 1 GB of
# single-cell matrices, so only the six bulk samples are pulled individually.
SAMPLES = [
    ("GSM7775nnn/GSM7775364", "GSM7775364_NK_control_rep_1.txt.gz"),
    ("GSM7775nnn/GSM7775365", "GSM7775365_NK_control_rep_2.txt.gz"),
    ("GSM7775nnn/GSM7775366", "GSM7775366_NK_control_rep_3.txt.gz"),
    ("GSM7775nnn/GSM7775370", "GSM7775370_NK_cytokine_control_1.txt.gz"),
    ("GSM7775nnn/GSM7775371", "GSM7775371_NK_cytokine_control_2.txt.gz"),
    ("GSM7775nnn/GSM7775372", "GSM7775372_NK_cytokine_control_3.txt.gz"),
]

BASE = "https://ftp.ncbi.nlm.nih.gov/geo"


def fetch(url: str, dest: str):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  have {os.path.basename(dest)}")
        return
    for attempt in range(4):
        r = subprocess.run(["curl", "-sS", "--max-time", "300", "-o", dest, url])
        if r.returncode == 0 and os.path.getsize(dest) > 0:
            print(f"  got  {os.path.basename(dest)} "
                  f"({os.path.getsize(dest):,} bytes)")
            return
        subprocess.run(["sleep", str(2 ** (attempt + 1))])
    raise RuntimeError(f"failed to fetch {url}")


def main():
    os.makedirs(DATA, exist_ok=True)
    print("series supplementary files:")
    for path, name in SERIES:
        fetch(f"{BASE}/series/{path}/suppl/{name}", os.path.join(DATA, name))
    print("GSE242941 bulk samples:")
    for path, name in SAMPLES:
        fetch(f"{BASE}/samples/{path}/suppl/{name}", os.path.join(DATA, name))
    print("\nDONE")


if __name__ == "__main__":
    main()
