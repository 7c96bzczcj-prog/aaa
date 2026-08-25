#!/usr/bin/env python3
"""E2 -- donor age structure of the admitted human datasets.

Ages come only from the metadata already verified in manifest.yaml: HCA Azul
donorOrganisms.organismAge, and the GEO GSM characteristics field.  Where no
age is deposited it is recorded MISSING; nothing is inferred.
"""
import json, os
import numpy as np, pandas as pd

RES = "/home/user/aaa/results"

# HCA Census of Immune Cells -- Azul donorOrganisms.organismAge, one per donor
HCA = {"MantonBM1": (52, "BM"), "MantonBM2": (50, "BM"), "MantonBM3": (39, "BM"),
       "MantonBM4": (29, "BM"), "MantonBM5": (29, "BM"), "MantonBM6": (26, "BM"),
       "MantonBM7": (36, "BM"), "MantonBM8": (32, "BM"),
       "MantonBL1": (28, "PB"), "MantonBL2": (50, "PB"), "MantonBL3": (32, "PB"),
       "MantonBL4": (27, "PB"), "MantonBL5": (23, "PB"), "MantonBL6": (36, "PB"),
       "MantonBL7": (24, "PB"), "MantonBL8": (25, "PB")}
# GSE120221 -- GSM characteristics "age: N years old"; 25 libraries / 20 donors,
# C1/C2/Ck are one donor and S1/S2/Sk1/Sk2 another (series design: 20 donors)
G120_LIB = {"A": 59, "B": 47, "C1": 60, "Ck": 59, "C2": 60, "E": 30, "F": 41,
            "G": 58, "H": 50, "J": 43, "K": 84, "L": 57, "M": 60, "N": 67,
            "O": 50, "P": 58, "Q": 66, "R": 31, "S1": 56, "Sk1": 55, "S2": 56,
            "Sk2": 55, "T": 24, "U": 46, "W": 28}

rows = []
for d, (age, comp) in HCA.items():
    rows.append({"dataset": "HCA_ICA", "donor": d, "compartment": comp,
                 "age": age, "age_source": "HCA Azul donorOrganisms.organismAge"})
byd = {}
for lib, age in G120_LIB.items():
    donor = "".join(ch for ch in lib if ch.isalpha() and ch.isupper()) or lib[0]
    byd.setdefault(donor, []).append(age)
for donor, ages in byd.items():
    rows.append({"dataset": "GSE120221", "donor": donor, "compartment": "BM",
                 "age": float(np.median(ages)),
                 "age_source": f"GEO GSM characteristics; {len(ages)} library(ies), "
                               f"ages {sorted(ages)}"})
# datasets admitted on C1 that deposit no age at all
for ds, donors, comp in [("GSE233304", ["P1", "P2", "P4", "P7", "P8", "P9", "P15",
                                        "P16", "P21", "P22", "P24", "c4", "c5",
                                        "c8", "c9", "c10"], "BM"),
                         ("GSE181543", ["BM1", "BM2"], "BM"),
                         ("GSE181543", ["PBMC1", "PBMC2", "PBMC3"], "PB")]:
    for d in donors:
        rows.append({"dataset": ds, "donor": d, "compartment": comp,
                     "age": np.nan, "age_source": "MISSING -- no age deposited"})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES, "e2_age_distribution.tsv"), sep="\t", index=False)

DEC = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 70),
       (70, 80), (80, 90)]
summ = []
for (ds, comp), g in df.groupby(["dataset", "compartment"]):
    a = g.age.dropna()
    r = {"dataset": ds, "compartment": comp, "n_donors": len(g),
         "n_with_age": len(a),
         "age_min": a.min() if len(a) else np.nan,
         "age_max": a.max() if len(a) else np.nan}
    for lo, hi in DEC:
        r[f"n_{lo}_{hi}"] = int(((a >= lo) & (a < hi)).sum())
    summ.append(r)
s = pd.DataFrame(summ)
s.to_csv(os.path.join(RES, "e2_age_summary.tsv"), sep="\t", index=False)
print(s.to_string(index=False))

bm = df[(df.compartment == "BM")]
bm_age = bm.age.dropna()
print(f"\nBONE MARROW, admitted datasets only")
print(f"  donors total                 {len(bm)}")
print(f"  donors with an age deposited {len(bm_age)}  "
      f"({100*len(bm_age)/len(bm):.0f}%)  -> missing {len(bm)-len(bm_age)}")
print(f"  age range                    {bm_age.min():.0f} - {bm_age.max():.0f}")
print(f"  donors aged 25-55 inclusive  {int(((bm_age >= 25) & (bm_age <= 55)).sum())}")
print(f"  donors under 18              {int((bm_age < 18).sum())}")
print(f"  donors 55+                   {int((bm_age > 55).sum())}")
print(f"  donors under 25              {int((bm_age < 25).sum())}")
print(f"\n  datasets with NO age at all: "
      f"{sorted(set(bm[bm.age.isna()].dataset))}")
