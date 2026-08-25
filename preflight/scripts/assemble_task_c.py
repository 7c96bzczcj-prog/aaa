#!/usr/bin/env python3
"""Assemble results/early_nk_stages.tsv (spec section 7)."""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
import pf_core as pf  # noqa: E402

RES = "/home/user/aaa/results"
dn = pd.read_csv("out/early_nk_stages.tsv", sep="\t")
supp = pd.read_csv("out/task_c_supp.tsv", sep="\t")

# per-donor supplement roll-up (HCA only; the supplement was run on HCA BM)
s = supp.groupby("donor").agg(
    n_postQC_supp=("n_postQC", "sum"),
    n_postQC_ge1000genes=("n_postQC_ge1000genes", "sum"),
    prog_anchor_any_n=("prog_anchor_any_n", "sum"),
    prog_anchor_noCyto_n=("prog_anchor_noCyto_n", "sum"),
    earlyNK_permissive_n_ge1000genes=("earlyNK_permissive_n_ge1000genes", "sum"),
    CLP_like_n_ge1000genes=("CLP_like_n_ge1000genes", "sum"),
    det_NK_IL2RB=("det_NK_IL2RB", "mean"),
    det_PROG_CD34=("det_PROG_CD34", "mean"),
    det_PROG_KIT=("det_PROG_KIT", "mean"),
).reset_index()
out = dn.merge(s, on="donor", how="left")

# sensitivity-corrected bound: divide the permissive bound by the detection
# rate of its NK-priming term, measured in cells that certainly carry it
out["earlyNK_permissive_postQC_X_upper95_sensCorrected"] = (
    out["earlyNK_permissive_postQC_X_upper95"] / out["det_NK_IL2RB"])
cols = ["dataset", "donor", "compartment", "n_libraries", "n_preQC", "n_postQC",
        "n_lost_to_QC", "median_genes_postQC", "frac_postQC_under_1000_genes",
        "n_postQC_ge1000genes",
        "earlyNK_permissive_preQC_n", "earlyNK_permissive_postQC_n",
        "earlyNK_permissive_postQC_X_upper95",
        "earlyNK_permissive_postQC_X_upper95_sensCorrected",
        "earlyNK_permissive_n_ge1000genes",
        "earlyNK_strict_preQC_n", "earlyNK_strict_postQC_n",
        "earlyNK_TFvariant_preQC_n", "earlyNK_TFvariant_postQC_n",
        "CLP_like_preQC_n", "CLP_like_postQC_n", "CLP_like_postQC_X_upper95",
        "prog_anchor_any_n", "CD34_any_postQC_n",
        "proB_control_preQC_n", "proB_control_postQC_n",
        "eryProg_control_preQC_n", "eryProg_control_postQC_n",
        "matureNK_marker_postQC_n", "det_NK_IL2RB", "det_PROG_CD34", "det_PROG_KIT"]
cols = [c for c in cols if c in out.columns]
os.makedirs(RES, exist_ok=True)
out[cols].to_csv(os.path.join(RES, "early_nk_stages.tsv"), sep="\t", index=False)

# pooled summary appended as a second file the verdict quotes from
pool = pd.read_csv("out/early_nk_pooled.tsv", sep="\t")
pool.to_csv(os.path.join(RES, "early_nk_stages_pooled.tsv"), sep="\t", index=False)
print("wrote results/early_nk_stages.tsv", out.shape)
h = out[out.dataset == "HCA_ICA"]
N = int(h.n_postQC.sum()); k = int(h.earlyNK_permissive_postQC_n.sum())
kp = int(h.earlyNK_permissive_preQC_n.sum()); Np = int(h.n_preQC.sum())
print(f"HCA_ICA BM  pre-QC droplets={Np:,}  post-QC cells={N:,}")
print(f"  earlyNK permissive: preQC={kp}  postQC={k}  "
      f"X95={100*pf.clopper_pearson_upper(k, N):.4f}%  "
      f"sens-corrected={100*pf.clopper_pearson_upper(k, N)/h.det_NK_IL2RB.mean():.4f}%")
print(f"  CLP-like superset : postQC={int(h.CLP_like_postQC_n.sum())}  "
      f"X95={100*pf.clopper_pearson_upper(int(h.CLP_like_postQC_n.sum()), N):.4f}%")
print(f"  controls: proB={int(h.proB_control_postQC_n.sum())}  "
      f"eryProg={int(h.eryProg_control_postQC_n.sum())}  "
      f"matureNK={int(h.matureNK_marker_postQC_n.sum())}")
