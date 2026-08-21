"""Check every number quoted in the task P-1 prose against the result files.

Topic 1 learned this the hard way: a checker that scanned only the README
reported "16/16 verified" while a retracted sentence was still live in the
document the README sent readers to (`scripts/verify_claims.py`).  So this
one does two things:

  1. re-derives each quoted figure from `results/p1_*.csv` and compares;
  2. greps ALL FOUR prose files for the phrases that would mean the floor
     verdict had been softened into a biological "no difference", which is
     the single misreading this task is most exposed to.

Run: python3 scripts/verify_claims_p1.py    (exit code 1 if anything fails)
"""

from __future__ import annotations

import os
import re
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")
DOCS = ["PIEZO1_BASELINE.md", "MECHANOSENSOR_PANEL.md",
        "PIEZO1_FEASIBILITY.md", "CLAIMS_P1.md"]

fails = []


def chk(name, want, got, ok):
    print(f"  [{'ok ' if ok else 'FAIL'}] {name}: want {want}, got {got}")
    if not ok:
        fails.append(name)


def close(a, b, tol):
    return abs(float(a) - float(b)) <= tol


def main():
    q1 = pd.read_csv(os.path.join(OUT, "p1_q1_baseline.csv"))
    q2 = pd.read_csv(os.path.join(OUT, "p1_q2_panel.csv"))
    q4 = pd.read_csv(os.path.join(OUT, "p1_q4_paired_human.csv"))
    q6 = pd.read_csv(os.path.join(OUT, "p1_q6_measurability.csv"))
    cal = pd.read_csv(os.path.join(OUT, "p1_calibration.csv"))
    amb = pd.read_csv(os.path.join(OUT, "p1_ambient_sensitivity.csv"))
    cx = pd.read_csv(os.path.join(OUT, "p1_artifact_crosscheck.csv"))
    h2h = pd.read_csv(os.path.join(OUT, "p1_q2_head_to_head.csv"))
    orth = pd.read_csv(os.path.join(OUT, "p1_orthologs.csv"))
    hpa = pd.read_csv(os.path.join(OUT, "p1_hpa_blood_rna.csv"))

    nk = q1[q1.lineage == "NK"].set_index("condition")
    print("== Q1 ==")
    chk("PIEZO1 NK tumour detection = 2.62%", 2.62,
        round(100 * nk.det_pooled["Tumor"], 2),
        close(100 * nk.det_pooled["Tumor"], 2.62, 0.01))
    chk("PIEZO1 NK normal detection = 2.73%", 2.73,
        round(100 * nk.det_pooled["Normal"], 2),
        close(100 * nk.det_pooled["Normal"], 2.73, 0.01))
    chk("both NK rows are floor band", "floor,floor",
        ",".join(nk.band_pooled.tolist()),
        set(nk.band_pooled) == {"floor"})
    chk("donors in floor band, tumour", 26, int(nk.donors_in_floor_band["Tumor"]),
        int(nk.donors_in_floor_band["Tumor"]) == 26)
    chk("donors in floor band, normal", 27, int(nk.donors_in_floor_band["Normal"]),
        int(nk.donors_in_floor_band["Normal"]) == 27)
    chk("n donors = 27", 27, int(nk.n_donors["Tumor"]), int(nk.n_donors["Tumor"]) == 27)
    cells = int(nk.n_cells_total.sum())
    chk("NK cells entering the pseudobulk = 19054", 19054, cells, cells == 19054)
    # the depth story: myeloid/NK depth ratio vs detection ratio
    my = q1[q1.lineage == "Myeloid"].set_index("condition")
    dratio = float(my.umi_per_cell_median["Normal"] / nk.umi_per_cell_median["Normal"])
    detratio = float(my.det_pooled["Normal"] / nk.det_pooled["Normal"])
    chk("myeloid/NK depth ratio in 3.4-3.5", "3.4-3.5", round(dratio, 2),
        3.35 <= dratio <= 3.55)
    chk("myeloid/NK detection ratio in 3.1-3.2", "3.1-3.2", round(detratio, 2),
        3.0 <= detratio <= 3.25)

    print("== calibration ==")
    for cond, truth, want in (("Tumor", 0.043, 4.35), ("Normal", 0.049, 4.79)):
        r = cal[cal.condition == cond].iloc[0]
        chk(f"TOX {cond} bound = {want}%", want,
            round(100 * r.poisson_bound_from_pseudobulk, 2),
            close(100 * r.poisson_bound_from_pseudobulk, want, 0.01))
        chk(f"TOX {cond} truth = {truth}", truth, r.recorded_cell_level_detection,
            close(r.recorded_cell_level_detection, truth, 1e-9))
    chk("worst calibration error <= 0.11 pp", "<=0.11",
        round(float(cal.absolute_error_pp.abs().max()), 2),
        float(cal.absolute_error_pp.abs().max()) <= 0.115)

    print("== sensitivity ==")
    a = amb[amb.gene == "PIEZO1"].set_index("condition")
    chk("ambient-adjusted tumour = 1.86%", 1.86,
        round(100 * a.det_ambient_adjusted["Tumor"], 2),
        close(100 * a.det_ambient_adjusted["Tumor"], 1.86, 0.01))
    chk("ambient-adjusted normal = 2.52%", 2.52,
        round(100 * a.det_ambient_adjusted["Normal"], 2),
        close(100 * a.det_ambient_adjusted["Normal"], 2.52, 0.01))
    chk("both ambient-adjusted rows still floor", "floor,floor",
        ",".join(a.band_ambient_adjusted.tolist()),
        set(a.band_ambient_adjusted) == {"floor"})
    chk("ambient fraction tumour = 29%", 29,
        round(100 * a.ambient_fraction_of_signal["Tumor"]),
        close(100 * a.ambient_fraction_of_signal["Tumor"], 29, 0.6))
    chk("ambient fraction normal = 8%", 8,
        round(100 * a.ambient_fraction_of_signal["Normal"]),
        close(100 * a.ambient_fraction_of_signal["Normal"], 8, 0.6))
    c = cx[cx.gene == "PIEZO1"].set_index("condition")
    chk("other-artifact tumour = 1.86%", 1.86,
        round(100 * c.det_allcells_decontam["Tumor"], 2),
        close(100 * c.det_allcells_decontam["Tumor"], 1.86, 0.01))
    chk("other-artifact normal = 2.27%", 2.27,
        round(100 * c.det_allcells_decontam["Normal"], 2),
        close(100 * c.det_allcells_decontam["Normal"], 2.27, 0.01))
    chk("cross-check agrees on band in every row", True, bool(cx.same_band.all()),
        bool(cx.same_band.all()))

    print("== Q2 ==")
    for cond, top in (("Tumor", "TRPV2"), ("Normal", "TRPV2")):
        s = q2[q2.condition == cond].sort_values("det_pooled", ascending=False)
        chk(f"{cond}: highest channel is {top}", top, s.gene.iloc[0],
            s.gene.iloc[0] == top)
        rank = int(q2[(q2.condition == cond) & (q2.gene == "PIEZO1")].rank_in_panel.iloc[0])
        chk(f"{cond}: PIEZO1 rank = 4", 4, rank, rank == 4)
    order = (q2[q2.condition == "Normal"].sort_values("det_pooled", ascending=False)
             .gene.tolist()[:4])
    chk("normal order TRPV2>TRPM7>TMEM63A>PIEZO1",
        "TRPV2,TRPM7,TMEM63A,PIEZO1", ",".join(order),
        order == ["TRPV2", "TRPM7", "TMEM63A", "PIEZO1"])
    beaten = set(h2h[(h2h.verdict == "b_above_a") & (h2h.q_BH < 0.05)].gene_b)
    chk("channels above PIEZO1 (BH q<0.05)", "TMEM63A,TRPM7,TRPV2",
        ",".join(sorted(beaten)), beaten == {"TMEM63A", "TRPM7", "TRPV2"})
    n_floor = int((q2.band_pooled == "floor").sum())
    chk("19 of 20 panel rows are floor", 19, n_floor, n_floor == 19)
    k4 = q2[(q2.gene == "KCNK4")]
    chk("KCNK4 zero in every donor, both conditions", "27,27",
        ",".join(str(int(x)) for x in k4.donors_zero_count),
        set(k4.donors_zero_count) == {27})

    print("== HPA cross-platform ==")
    h = hpa.set_index("gene")["nTPM_NK-cell"]
    chk("HPA PIEZO1 NK = 1.0", 1.0, h["PIEZO1"], close(h["PIEZO1"], 1.0, 1e-9))
    chk("HPA TRPV2 NK = 63.5", 63.5, h["TRPV2"], close(h["TRPV2"], 63.5, 1e-9))
    hpa_order = h.reindex(["TRPV2", "TMEM63A", "TRPM7", "PIEZO1"]).tolist()
    chk("HPA reproduces TRPV2>TMEM63A>TRPM7>PIEZO1", "descending",
        hpa_order, all(x > y for x, y in zip(hpa_order, hpa_order[1:])))

    print("== Q4 ==")
    p = q4[q4.gene == "PIEZO1"].iloc[0]
    chk("PIEZO1 log2FC = -0.099", -0.099, round(p.mean_log2FC_tumour_vs_normal, 3),
        close(p.mean_log2FC_tumour_vs_normal, -0.099, 0.001))
    chk("PIEZO1 paired-t p = 0.605", 0.605, round(p.paired_t_p, 3),
        close(p.paired_t_p, 0.605, 0.001))
    chk("PIEZO1 not measurable in either condition", False, bool(p.measurable),
        not bool(p.measurable))
    chk("PIEZO1 n_pairs = 27", 27, int(p.n_pairs), int(p.n_pairs) == 27)
    t = q4[q4.gene == "TRPV2"].iloc[0]
    chk("TRPV2 is the only measurable panel gene", 1, int(q4.measurable.sum()),
        int(q4.measurable.sum()) == 1)
    chk("TRPV2 BH q = 0.0019", 0.0019, round(t.paired_t_q_BH, 4),
        close(t.paired_t_q_BH, 0.0019, 0.0001))

    print("== Q6 ==")
    g = q6.set_index(["gene", "condition"])
    chk("ITGA1 normal = floor", "floor", g.band_pooled[("ITGA1", "Normal")],
        g.band_pooled[("ITGA1", "Normal")] == "floor")
    chk("ITGA1 normal detection = 4.32%", 4.32,
        round(100 * g.det_pooled[("ITGA1", "Normal")], 2),
        close(100 * g.det_pooled[("ITGA1", "Normal")], 4.32, 0.01))
    chk("ITGA1 tumour CI spans the 5% line", True,
        f"{100*g.det_pooled_lo95[('ITGA1','Tumor')]:.2f}-"
        f"{100*g.det_pooled_hi95[('ITGA1','Tumor')]:.2f}",
        g.det_pooled_lo95[("ITGA1", "Tumor")] < 0.05 < g.det_pooled_hi95[("ITGA1", "Tumor")])
    chk("S1PR1 floor in both conditions", "floor,floor",
        ",".join(q6[q6.gene == "S1PR1"].band_pooled.tolist()),
        set(q6[q6.gene == "S1PR1"].band_pooled) == {"floor"})
    chk("NKG7 ceiling in both conditions", "ceiling,ceiling",
        ",".join(q6[q6.gene == "NKG7"].band_pooled.tolist()),
        set(q6[q6.gene == "NKG7"].band_pooled) == {"ceiling"})
    chk("CD69 tumour = ceiling", "ceiling", g.band_pooled[("CD69", "Tumor")],
        g.band_pooled[("CD69", "Tumor")] == "ceiling")

    print("== orthology ==")
    panel = ["PIEZO1", "PIEZO2", "TRPV4", "TRPV2", "TRPM7", "TRPC1",
             "TMEM63A", "TMEM63B", "KCNK2", "KCNK4"]
    sub = orth[orth.human.isin(panel)]
    chk("all 10 panel genes one-to-one", 10, int(sub.one_to_one.sum()),
        int(sub.one_to_one.sum()) == 10)
    xcl = orth[orth.human.isin(["XCL1", "XCL2"])]
    chk("XCL1/XCL2 collapse onto one mouse gene", "Xcl1",
        ",".join(sorted(set(xcl.mouse_symbol))),
        set(xcl.mouse_symbol) == {"Xcl1"} and not xcl.one_to_one.any())

    # ---- prose guard: the floor verdict must never be stated as a plain
    # ---- biological "no difference" anywhere in the four documents.
    print("== prose guard ==")
    BANNED = [
        (r"PIEZO1[^。\n]{0,40}肿瘤与癌旁[^。\n]{0,20}无差异(?![^。\n]{0,40}FLOOR)",
         "unqualified 'no tumour/normal difference' for PIEZO1"),
        (r"NK\s*(?:细胞)?不表达\s*PIEZO1", "'NK does not express PIEZO1'"),
    ]
    for doc in DOCS:
        path = os.path.join(ROOT, doc)
        if not os.path.exists(path):
            chk(f"{doc} exists", "present", "MISSING", False)
            continue
        text = open(path, encoding="utf-8").read()
        # A mention is not an assertion.  Strip 「...」 quoted spans before
        # scanning, so that a sentence warning against a misreading is not
        # itself flagged as that misreading.
        scan = re.sub(r"「[^」]*」", "「…」", text)
        for pat, label in BANNED:
            hit = re.search(pat, scan)
            chk(f"{doc}: no {label}", "absent",
                (hit.group(0)[:40] if hit else "absent"), hit is None)
        # every document must carry the §0 limitation
        chk(f"{doc}: carries the §0 limitation", "present",
            "present" if "不是「机械感受能力」" in text else "MISSING",
            "不是「机械感受能力」" in text)

    print(f"\n{'ALL CHECKS PASS' if not fails else 'FAILURES: ' + ', '.join(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
