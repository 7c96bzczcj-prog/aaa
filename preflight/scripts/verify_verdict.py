#!/usr/bin/env python3
"""Re-check every number quoted in VERDICT.md against results/.

Same role as the repository's existing verify_claims.py: the verdict is
hand-written, so it needs a machine that can contradict it.
"""
import os, sys
import numpy as np, pandas as pd

RES = "/home/user/aaa/results"
fails, checks = [], 0


def eq(name, got, want, tol=None):
    global checks
    checks += 1
    ok = (abs(got - want) <= (tol if tol is not None else abs(want) * 0.02 + 1e-9))
    print(f"{'ok ' if ok else 'FAIL'} {name:<56} got={got:<12.5g} verdict={want:<12.5g}")
    if not ok:
        fails.append(name)


con = pd.read_csv(f"{RES}/soup_by_compartment.tsv", sep="\t")
tst = pd.read_csv(f"{RES}/target_gene_detection.tsv", sep="\t")
vd = pd.read_csv(f"{RES}/task_a_verdict.tsv", sep="\t")
acc = pd.read_csv(f"{RES}/soup_model_acceptance.tsv", sep="\t")

print("== task A ==")
for ds, rho_bm, rho_pb, ratio, p, dmax, nsig, ntest in [
        ("GSE233304", 0.0090, 0.0028, 3.23, 0.32, 0.0081, 0, 11),
        ("HCA_ICA", 0.0175, 0.0154, 1.13, 0.39, 0.0076, 0, 14)]:
    r = con[(con.dataset == ds) & (con.quantity == "rho_narrow")].iloc[0]
    eq(f"{ds} rho_BM", r.mean_bm, rho_bm, 0.0002)
    eq(f"{ds} rho_PB", r.mean_pb, rho_pb, 0.0002)
    eq(f"{ds} rho ratio", r.mean_bm / r.mean_pb, ratio, 0.02)
    eq(f"{ds} rho p", r.p, p, 0.01)
    g = con[(con.dataset == ds) & (con.gene != "")].merge(
        tst[(tst.dataset == ds)][["gene", "testable"]], on="gene")
    g = g[g.testable]
    eq(f"{ds} n testable genes", len(g), ntest, 0)
    eq(f"{ds} max |Delta_ambient|", g.delta.abs().max(), dmax, 0.0002)
    eq(f"{ds} genes over 0.10 effect at FDR<0.05",
       int(((g.q_bh < 0.05) & (g.delta.abs() >= 0.10)).sum()), nsig, 0)
    v = vd[vd.dataset == ds].iloc[0]
    print(f"{'ok ' if v.verdict_key == 'A_PASS' else 'FAIL'} {ds} verdict key = {v.verdict_key}")
    if v.verdict_key != "A_PASS":
        fails.append(f"{ds} verdict")
eq("HCA genes clearing FDR", int((con[(con.dataset == "HCA_ICA") & (con.gene != "")]
                                  .merge(tst[tst.dataset == "HCA_ICA"][["gene", "testable"]], on="gene")
                                  .query("testable and q_bh < 0.05")).shape[0]), 3, 0)
eq("acceptance: compartments passing", int(acc.accept_pass.sum()), len(acc), 0)
eq("acceptance: min frac within 3x", acc.frac_within_3x.min(), 0.70, 0.01)
eq("acceptance: max frac within 3x", acc.frac_within_3x.max(), 1.00, 0.01)

print("\n== task A, cross-lineage rho ==")
rl = pd.read_csv(f"{RES}/rho_by_lineage.tsv", sep="\t")
ratios = []
for (ds, comp), g in rl[rl.lineage.isin(["NK", "T"])].groupby(["dataset", "compartment"]):
    r = {l: np.average(x.rho, weights=x.n_cells) for l, x in g.groupby("lineage")}
    if "NK" in r and "T" in r:
        ratios.append(r["NK"] / r["T"])
eq("min rho_NK/rho_T", min(ratios), 0.44, 0.02)
eq("max rho_NK/rho_T", max(ratios), 1.86, 0.02)
nk = rl[rl.lineage == "NK"]
rho_nk = [np.average(g.rho, weights=g.n_cells)
          for _, g in nk.groupby(["dataset", "compartment"])]
eq("min rho_NK across dataset x compartment", min(rho_nk), 0.004, 0.0005)
eq("max rho_NK across dataset x compartment", max(rho_nk), 0.019, 0.0005)

print("\n== task B ==")
gate = pd.read_csv(f"{RES}/ltnk_detection_gate.tsv", sep="\t")
eq("task B: gate passes anywhere", int(gate.gate_pass.sum()), 0, 0)
eq("CXCR6 min detection %", 100 * gate.det_CXCR6.min(), 0.29, 0.01)
eq("CXCR6 max detection %", 100 * gate.det_CXCR6.max(), 1.32, 0.01)
for ds, comp, cx, cd, fl in [("GSE120221", "BM", 0.62, 60.5, 30.7),
                             ("HCA_ICA", "BM", 0.44, 32.9, 4.24)]:
    r = gate[(gate.dataset == ds) & (gate.compartment == comp)].iloc[0]
    eq(f"{ds} {comp} CXCR6 %", 100 * r.det_CXCR6, cx, 0.02)
    eq(f"{ds} {comp} CD69 %", 100 * r.det_CD69, cd, 0.2)
    eq(f"{ds} {comp} CD69 erythroid floor %", 100 * r.floor_Erythroid_CD69, fl, 0.2)

print("\n== task C ==")
c = pd.read_csv(f"{RES}/early_nk_stages.tsv", sep="\t")
h = c[c.dataset == "HCA_ICA"]
eq("HCA BM post-QC cells", h.n_postQC.sum(), 374305, 0)
eq("HCA BM pre-QC droplets", h.n_preQC.sum(), 441171, 0)
eq("early NK pre-QC", h.earlyNK_permissive_preQC_n.sum(), 25, 0)
eq("early NK post-QC", h.earlyNK_permissive_postQC_n.sum(), 25, 0)
eq("early NK surviving >=1000 genes", h.earlyNK_permissive_n_ge1000genes.sum(), 20, 0)
eq("all cells surviving >=1000 genes %",
   100 * h.n_postQC_ge1000genes.sum() / h.n_postQC.sum(), 18.1, 0.1)
eq("pro-B control", h.proB_control_postQC_n.sum(), 1284, 0)
eq("erythroid progenitor control", h.eryProg_control_postQC_n.sum(), 2981, 0)
eq("mature NK", h.matureNK_marker_postQC_n.sum(), 13437, 0)
eq("CLP-like superset", h.CLP_like_postQC_n.sum(), 331, 0)
sys.path.insert(0, "/home/user/aaa/preflight/lib")
import pf_core as pf
X = pf.clopper_pearson_upper(int(h.earlyNK_permissive_postQC_n.sum()), int(h.n_postQC.sum()))
eq("X upper 95 %", 100 * X, 0.0093, 0.0002)
eq("X sensitivity-corrected %", 100 * X / h.det_NK_IL2RB.mean(), 0.031, 0.001)
eq("CLP-like X upper 95 %",
   100 * pf.clopper_pearson_upper(int(h.CLP_like_postQC_n.sum()), int(h.n_postQC.sum())),
   0.0968, 0.001)
eq("IL2RB detection in mature NK", h.det_NK_IL2RB.mean(), 0.299, 0.001)
g2 = c[(c.dataset == "GSE233304") & (c.compartment == "BM")]
eq("GSE233304 marrow X upper 95 %",
   100 * pf.clopper_pearson_upper(int(g2.earlyNK_permissive_postQC_n.sum()),
                                  int(g2.n_postQC.sum())), 0.112, 0.002)

print("\n== E1 (spec v1.2) ==")
import json as _json
e1 = _json.load(open(f"{RES}/e1_gate.json"))
eq("E1 genes passing the T3 gate", e1["n_passing"], 6, 0)
eq("E1 donors", e1["n_donors"], 4, 0)
eq("E1 delta median", e1["delta_median"], -0.046, 0.001)
eq("E1 delta Q1", e1["delta_q1"], -0.088, 0.001)
eq("E1 delta Q3", e1["delta_q3"], -0.005, 0.001)
eq("E1 donors with positive delta", e1["n_positive"], 1, 0)
gt = pd.DataFrame(e1["gate"]).set_index("gene")
eq("E1 ZEB2 floor", gt.loc["ZEB2", "floor_used"], 0.547, 0.002)
eq("E1 ZEB2 fold", gt.loc["ZEB2", "fold"], 1.55, 0.02)
eq("E1 IL7R detection", gt.loc["IL7R", "det_NK"], 0.076, 0.001)
eq("E1 IL7R floor", gt.loc["IL7R", "floor_used"], 0.063, 0.001)

print("\n== E2 (spec v1.2) ==")
a = pd.read_csv(f"{RES}/e2_age_distribution.tsv", sep="\t")
bm = a[a.compartment == "BM"]
age = bm.age.dropna()
eq("E2 marrow donors", len(bm), 46, 0)
eq("E2 marrow donors with an age", len(age), 28, 0)
eq("E2 marrow donors without an age", len(bm) - len(age), 18, 0)
eq("E2 pct with age", 100 * len(age) / len(bm), 61, 1)
eq("E2 age min", age.min(), 24, 0)
eq("E2 age max", age.max(), 84, 0)
eq("E2 donors 25-55", int(((age >= 25) & (age <= 55)).sum()), 17, 0)
eq("E2 donors under 18", int((age < 18).sum()), 0, 0)
eq("E2 donors under 25", int((age < 25).sum()), 1, 0)

print(f"\n{checks} checks, {len(fails)} failed")
if fails:
    print("FAILED:", "; ".join(fails))
    sys.exit(1)
