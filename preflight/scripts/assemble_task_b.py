#!/usr/bin/env python3
"""Assemble results/ltnk_age.tsv and results/ltnk_detection_gate.tsv."""
import json, os, shutil
import numpy as np, pandas as pd

RES = "/home/user/aaa/results"
AGE_HCA = {"MantonBM1": 52, "MantonBM2": 50, "MantonBM3": 39, "MantonBM4": 29,
           "MantonBM5": 29, "MantonBM6": 26, "MantonBM7": 36, "MantonBM8": 32,
           "MantonBL1": 28, "MantonBL2": 50, "MantonBL3": 32, "MantonBL4": 27,
           "MantonBL5": 23, "MantonBL6": 36, "MantonBL7": 24, "MantonBL8": 25}

dn = pd.read_csv("out/ltnk_age.tsv", sep="\t")
dn["age"] = [AGE_HCA.get(d, a) for d, a in zip(dn.donor, dn.age)]
front = ["dataset", "donor", "compartment", "age", "n_libraries", "n_nk",
         "det_CD69", "det_CXCR6", "floorB_CD69", "floorErythroid_CD69",
         "floorB_CXCR6", "floorErythroid_CXCR6",
         "n_ltNK", "frac_ltNK", "n_bright", "frac_bright", "n_dim", "frac_dim"]
cols = [c for c in front if c in dn.columns] + \
       [c for c in dn.columns if c not in front]
dn[cols].to_csv(os.path.join(RES, "ltnk_age.tsv"), sep="\t", index=False)
for f in ["ltnk_detection_gate.tsv", "ltnk_age_regression.tsv"]:
    if os.path.exists(os.path.join("out", f)):
        shutil.copy(os.path.join("out", f), os.path.join(RES, f))
print("wrote results/ltnk_age.tsv", dn.shape)
g = pd.read_csv(os.path.join(RES, "ltnk_detection_gate.tsv"), sep="\t")
print("\ngate outcome per dataset x compartment:")
print(g[["dataset", "compartment", "n_donors", "n_nk", "det_CD69",
         "floor_Erythroid_CD69", "det_CXCR6", "floor_Erythroid_CXCR6",
         "gate_CD69_pass", "gate_CXCR6_pass", "gate_pass"]].to_string(index=False))
print(f"\ngate passes anywhere: {bool(g.gate_pass.any())}")
print(f"CXCR6 detection range across all six: "
      f"{100*g.det_CXCR6.min():.2f}% - {100*g.det_CXCR6.max():.2f}%  "
      f"(admission floor 5%)")
