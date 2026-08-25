#!/usr/bin/env python3
"""Download the bone-marrow and blood samples of GSE233304 (cranial bone
marrow + matched PBMC, CD45+, unfiltered droplet matrices).  Tumour samples
are not downloaded -- the preflight does not touch tumour NK."""
import os, re, sys, time, urllib.request

DEST = "/home/user/aaa/data/gse233304"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples"
# GSM -> (sample id, compartment, patient, group)
SAMPLES = {
    "GSM7421860": ("S1", "BM", "P1", "GBM"),   "GSM7421862": ("S3", "PB", "P1", "GBM"),
    "GSM7421863": ("S4", "BM", "P2", "GBM"),   "GSM7421865": ("S6", "PB", "P2", "GBM"),
    "GSM7421866": ("S7", "BM", "P7", "GBM"),
    "GSM7421867": ("S8", "BM", "P8", "GBM"),
    "GSM7421868": ("S9", "BM", "P9", "GBM"),   "GSM7421870": ("S11", "PB", "P9", "GBM"),
    "GSM7421871": ("S12", "BM", "c4", "control"), "GSM7421872": ("S13", "PB", "c4", "control"),
    "GSM7421873": ("S14", "BM", "c5", "control"), "GSM7421874": ("S15", "PB", "c5", "control"),
    "GSM8244370": ("S16", "BM", "c8", "control"), "GSM8244371": ("S17", "PB", "c8", "control"),
    "GSM8244372": ("S18", "BM", "c9", "control"), "GSM8244373": ("S19", "PB", "c9", "control"),
    "GSM8244374": ("S20", "BM", "c10", "control"), "GSM8244375": ("S21", "PB", "c10", "control"),
    "GSM8244376": ("S22", "BM", "P4", "GBM"),
    "GSM8244380": ("S26", "BM", "P15", "GBM"), "GSM8244384": ("S30", "PB", "P15", "GBM"),
    "GSM8244386": ("S32", "BM", "P16", "GBM"), "GSM8244390": ("S36", "PB", "P16", "GBM"),
    "GSM8244392": ("S38", "BM", "P21", "GBM"), "GSM8244396": ("S42", "BMMC", "P21", "GBM"),
    "GSM8244398": ("S44", "BM", "P22", "GBM"), "GSM8244402": ("S48", "BMMC", "P22", "GBM"),
    "GSM8244404": ("S50", "BM", "P24", "GBM"), "GSM8244408": ("S54", "BMMC", "P24", "GBM"),
}


def fetch(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "cached"
    for a in range(5):
        try:
            urllib.request.urlretrieve(url, out + ".part")
            os.rename(out + ".part", out)
            return "ok"
        except Exception as e:
            err = e
            time.sleep(2 ** a)
    return f"FAIL {os.path.basename(out)} {err}"


if __name__ == "__main__":
    os.makedirs(DEST, exist_ok=True)
    for gsm, (sid, comp, pat, grp) in SAMPLES.items():
        pre = "GSM" + gsm[3:-3] + "nnn"
        for kind in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
            fn = f"{gsm}_{sid}_{kind}"
            r = fetch(f"{FTP}/{pre}/{gsm}/suppl/{fn}", os.path.join(DEST, fn))
            if r.startswith("FAIL"):
                print(r, flush=True)
        print(f"{gsm} {sid} {comp} {pat} {grp} done", flush=True)
    print("ALL_DONE")
