#!/usr/bin/env python3
"""Emit the per-library index files the task drivers consume."""
import glob, json, os, re

D_HCA = "/home/user/aaa/data/hca_ica"
D_233 = "/home/user/aaa/data/gse233304"
D_181 = "/home/user/aaa/data/gse181543"
D_120 = "/home/user/aaa/data/gse120221"
LAB = "/home/user/aaa/preflight/labels"
OUT = "/home/user/aaa/preflight/index"
os.makedirs(OUT, exist_ok=True)

# ---- HCA -------------------------------------------------------------------
hca = []
for f in sorted(glob.glob(f"{D_HCA}/Manton*_raw_feature_bc_matrix.h5")):
    lib = os.path.basename(f).replace("_raw_feature_bc_matrix.h5", "")
    donor = lib.split("_")[0]
    if not re.match(r"^Manton(BM|BL)[1-8]$", donor):
        continue
    hca.append({"dataset": "HCA_ICA", "donor": donor,
                "compartment": "BM" if "BM" in donor else "PB",
                "library": lib, "path": f, "kind": "h5",
                "labels": f"{LAB}/{donor}.npz", "group": "healthy"})
json.dump(hca, open(f"{OUT}/hca.json", "w"), indent=1)

# ---- GSE233304 -------------------------------------------------------------
META = {"S1": ("P1", "BM", "GBM"), "S3": ("P1", "PB", "GBM"),
        "S4": ("P2", "BM", "GBM"), "S6": ("P2", "PB", "GBM"),
        "S7": ("P7", "BM", "GBM"), "S8": ("P8", "BM", "GBM"),
        "S9": ("P9", "BM", "GBM"), "S11": ("P9", "PB", "GBM"),
        "S12": ("c4", "BM", "control"), "S13": ("c4", "PB", "control"),
        "S14": ("c5", "BM", "control"), "S15": ("c5", "PB", "control"),
        "S16": ("c8", "BM", "control"), "S17": ("c8", "PB", "control"),
        "S18": ("c9", "BM", "control"), "S19": ("c9", "PB", "control"),
        "S20": ("c10", "BM", "control"), "S21": ("c10", "PB", "control"),
        "S22": ("P4", "BM", "GBM"),
        "S26": ("P15", "BM", "GBM"), "S30": ("P15", "PB", "GBM"),
        "S32": ("P16", "BM", "GBM"), "S36": ("P16", "PB", "GBM"),
        "S38": ("P21", "BM", "GBM"), "S42": ("P21", "BMMC", "GBM"),
        "S44": ("P22", "BM", "GBM"), "S48": ("P22", "BMMC", "GBM"),
        "S50": ("P24", "BM", "GBM"), "S54": ("P24", "BMMC", "GBM")}
g233 = []
for f in sorted(glob.glob(f"{D_233}/*_barcodes.tsv.gz")):
    gsm, sid = os.path.basename(f).split("_")[:2]
    donor, comp, grp = META[sid]
    stem = f"{D_233}/{gsm}_{sid}"
    g233.append({"dataset": "GSE233304", "donor": donor, "compartment": comp,
                 "library": sid,
                 "path": f"{stem}_barcodes.tsv.gz|{stem}_features.tsv.gz|{stem}_matrix.mtx.gz",
                 "kind": "mtx", "labels": f"{LAB}/GSE233304_{donor}.npz", "group": grp})
json.dump(g233, open(f"{OUT}/gse233304.json", "w"), indent=1)

# ---- GSE181543 -------------------------------------------------------------
g181 = []
for f in sorted(glob.glob(f"{D_181}/*_raw_feature_bc_matrix.h5")):
    nm = os.path.basename(f).split("_")[1]
    g181.append({"dataset": "GSE181543", "donor": nm,
                 "compartment": "BM" if nm.startswith("BM") else "PB",
                 "library": nm, "path": f, "kind": "h5",
                 "labels": f"{LAB}/GSE181543_{nm}.npz",
                 "group": "healthy;BM is B-cell enriched"})
json.dump(g181, open(f"{OUT}/gse181543.json", "w"), indent=1)

# ---- GSE120221 (filtered; task B / task C only) -----------------------------
AGE = {"A": 59, "B": 47, "C1": 60, "Ck": 59, "C2": 60, "E": 30, "F": 41, "G": 58,
       "H": 50, "J": 43, "K": 84, "L": 57, "M": 60, "N": 67, "O": 50, "P": 58,
       "Q": 66, "R": 31, "S1": 56, "Sk1": 55, "S2": 56, "Sk2": 55, "T": 24,
       "U": 46, "W": 28}
g120 = []
for f in sorted(glob.glob(f"{D_120}/*_barcodes_*.tsv.gz")):
    b = os.path.basename(f)
    gsm = b.split("_")[0]
    sid = b.replace(f"{gsm}_barcodes_", "").replace(".tsv.gz", "")
    donor = re.sub(r"(k?\d*)$", "", sid) or sid      # C1/C2/Ck -> C ; S1/Sk2 -> S
    stem = f"{D_120}/{gsm}"
    g120.append({"dataset": "GSE120221", "donor": donor, "sample": sid,
                 "compartment": "BM", "library": sid, "age": AGE[sid],
                 "path": f"{stem}_barcodes_{sid}.tsv.gz|{stem}_genes_{sid}.tsv.gz|"
                         f"{stem}_matrix_{sid}.mtx.gz",
                 "kind": "mtx", "labels": f"{LAB}/GSE120221_{donor}.npz",
                 "group": "healthy"})
json.dump(g120, open(f"{OUT}/gse120221.json", "w"), indent=1)

for n, v in [("hca", hca), ("gse233304", g233), ("gse181543", g181), ("gse120221", g120)]:
    donors = sorted({e["donor"] for e in v})
    print(f"{n:<12} libraries={len(v):<4} donors={len(donors)} {donors}")
