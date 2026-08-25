#!/usr/bin/env python3
"""Assemble results/soup_by_compartment.tsv and its supporting files."""
import os, shutil
import pandas as pd

RES = "/home/user/aaa/results"
os.makedirs(RES, exist_ok=True)
for src, dst in [("out/soup_by_compartment.tsv", "soup_by_compartment.tsv"),
                 ("out/soup_by_donor.tsv", "soup_by_donor.tsv"),
                 ("out/soup_by_library.tsv", "soup_by_library.tsv"),
                 ("out/soup_model_acceptance.tsv", "soup_model_acceptance.tsv"),
                 ("out/target_gene_detection.tsv", "target_gene_detection.tsv"),
                 ("out/task_a_verdict.tsv", "task_a_verdict.tsv")]:
    shutil.copy(src, os.path.join(RES, dst))

# cross-lineage rho, both datasets in one file
frames = []
for f, ds in [("out/rho_lineage_hca.tsv", "HCA_ICA"),
              ("out/rho_lineage_gse233304.tsv", "GSE233304")]:
    if os.path.exists(f):
        frames.append(pd.read_csv(f, sep="\t"))
if frames:
    pd.concat(frames, ignore_index=True).to_csv(
        os.path.join(RES, "rho_by_lineage.tsv"), sep="\t", index=False)

d = pd.read_csv(os.path.join(RES, "soup_by_compartment.tsv"), sep="\t")
t = pd.read_csv(os.path.join(RES, "target_gene_detection.tsv"), sep="\t")
m = d[d.gene != ""].merge(t[["dataset", "gene", "testable"]], on=["dataset", "gene"])
print("=== Delta_ambient on testable genes only ===")
for ds, g in m[m.testable].groupby("dataset"):
    print(f"\n{ds}  (paired={bool(g.paired.iloc[0])})  n_testable={len(g)}")
    print(f"  |Delta| max = {g.delta.abs().max():.4f}   "
          f"pre-registered effect threshold = 0.10")
    print(f"  q_bh < 0.05 : {(g.q_bh < 0.05).sum()} genes "
          f"({', '.join(g[g.q_bh < 0.05].gene)})")
    print(f"  q_bh < 0.05 AND |Delta| >= 0.10 : "
          f"{((g.q_bh < 0.05) & (g.delta.abs() >= 0.10)).sum()}")
    print(f"  sign of Delta: {int((g.delta > 0).sum())} BM-higher, "
          f"{int((g.delta < 0).sum())} PB-higher")
