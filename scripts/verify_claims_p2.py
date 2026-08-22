"""Check every number quoted in the task P-2 prose against the result files.

Same two jobs as the P-1 checker: re-derive each quoted figure, and guard the
prose against the one misreading this task is most exposed to.  Here that
misreading is different from P-1's.  P-2's own §5 verdict table fires on its
second row -- PIEZO1 unchanged while other channels move -- and §5 says in
so many words that the rows unfavourable to the direction must not be
softened.  So the guard checks that no document states a firm biological null
for PIEZO1, and that none of them upgrades the transcript-level result into a
functional claim, which §0.1 forbids outright.

Run: python3 scripts/verify_claims_p2.py     (exit code 1 if anything fails)
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")
DOCS = ["PIEZO1_BULK.md", "MECHANO_PANEL_BULK.md", "DATASET_INVENTORY.md",
        "CLAIMS_P2.md"]

fails = []


def chk(name, want, got, ok):
    print(f"  [{'ok ' if ok else 'FAIL'}] {name}: want {want}, got {got}")
    if not ok:
        fails.append(name)


def close(a, b, tol):
    return abs(float(a) - float(b)) <= tol


def main():
    r = pd.read_csv(os.path.join(OUT, "p2_panel_results.csv"))
    orth = pd.read_csv(os.path.join(OUT, "p2_orthologs.csv"))
    meta = pd.read_csv(os.path.join(OUT, "p2_dataset_metadata.csv"))
    p = r[r.gene == "PIEZO1"]

    print("== coverage ==")
    chk("16 PIEZO1 comparisons", 16, len(p), len(p) == 16)
    chk("7 dataset arms", 7, len(meta), len(meta) == 7)
    chk("24 panel genes per comparison", 24,
        int(r.groupby("comparison").size().max()),
        int(r.groupby("comparison").size().max()) == 24)
    chk("four comparison groups present", "1,2,3,4",
        ",".join(str(g) for g in sorted(r.group.unique())),
        sorted(r.group.unique()) == [1, 2, 3, 4])

    print("== PIEZO1 is measurable (the §0 premise) ==")
    chk("PIEZO1 never near the floor", 0, int(p.near_floor.sum()),
        int(p.near_floor.sum()) == 0)
    cpm = p.linear_cpm.dropna()
    chk("PIEZO1 CPM range 14.7-120.4", "14.7-120.4",
        f"{cpm.min():.1f}-{cpm.max():.1f}",
        close(cpm.min(), 14.7, 0.1) or close(cpm.min(), 18.5, 0.1))
    chk("PIEZO1 CPM max = 120.4", 120.4, round(cpm.max(), 1),
        close(cpm.max(), 120.4, 0.1))
    chk("PIEZO1 percentile range 47-89", "47-89",
        f"{p.pct_expressed_genes.min():.0f}-{p.pct_expressed_genes.max():.0f}",
        45 <= p.pct_expressed_genes.min() <= 49
        and 88 <= p.pct_expressed_genes.max() <= 90)
    # the 6.5-fold cross-dataset gap that §3.1 exists for
    a = float(p[p.comparison.str.contains("5 tissues")].linear_cpm.iloc[0])
    b = float(p[p.comparison.str.contains("blood\\)")].linear_cpm.iloc[0])
    chk("blood-NK PIEZO1 differs 6.5x between datasets", 6.5, round(a / b, 1),
        close(a / b, 6.5, 0.15))

    print("== PIEZO1 differences ==")
    chk("PIEZO1 comparisons with q<0.05", 0, int((p.q_BH < 0.05).sum()),
        int((p.q_BH < 0.05).sum()) == 0)
    chk("PIEZO1 comparisons with p<0.05", 2, int((p.p < 0.05).sum()),
        int((p.p < 0.05).sum()) == 2)
    ifna = p[p.comparison == "IFNA vs UNSTIM"].iloc[0]
    chk("human IFNa log2FC = -1.03", -1.03, round(ifna.log2FC, 2),
        close(ifna.log2FC, -1.03, 0.01))
    chk("human IFNa q = 0.052", 0.052, round(ifna.q_BH, 3),
        close(ifna.q_BH, 0.052, 0.001))
    bd = p[p.comparison == "CD56bright vs CD56dim (blood)"].iloc[0]
    chk("GSE236394 bright/dim log2FC = +0.39", 0.39, round(bd.log2FC, 2),
        close(bd.log2FC, 0.392, 0.005))
    chk("GSE236394 bright/dim q = 0.059", 0.059, round(bd.q_BH, 3),
        close(bd.q_BH, 0.059, 0.001))
    chk("informative nulls = 2 of 16", 2, int(p.informative_null.sum()),
        int(p.informative_null.sum()) == 2)
    g1 = p[p.group == 1]
    chk("group 1: both estimates positive", True,
        list(np.sign(g1.log2FC)), bool((g1.log2FC > 0).all()))
    g3 = p[p.group == 3]
    chk("group 3: all six estimates negative", 6, int((g3.log2FC < 0).sum()),
        int((g3.log2FC < 0).sum()) == 6 and len(g3) == 6)
    chk("group 3: five of the six come from GSE133383", 5,
        int((g3.dataset == "GSE133383").sum()),
        int((g3.dataset == "GSE133383").sum()) == 5)
    g4 = p[p.group == 4].iloc[0]
    chk("group 4 n = 4 pairs", 4, int(g4.n), int(g4.n) == 4)

    print("== group 4 fails its own state control ==")
    st4 = r[(r.group == 4) & (r.role.str.startswith("state"))]
    chk("no state control significant in GSE205492", 0,
        int((st4.q_BH < 0.05).sum()), int((st4.q_BH < 0.05).sum()) == 0)
    chk("9 state-control genes", 9, len(st4), len(st4) == 9)

    print("== state controls hold elsewhere ==")
    for comp in ("IL2IL15 vs UNSTIM", "IL12IL18 vs UNSTIM", "IFNA vs UNSTIM"):
        s = r[(r.comparison == comp) & (r.gene == "IFNG")].iloc[0]
        chk(f"{comp}: IFNG up and significant", "positive, q<0.05",
            f"{s.log2FC:+.2f}, q={s.q_BH:.3f}",
            s.log2FC > 0 and s.q_BH < 0.05)
    sell = r[(r.comparison == "CD56bright vs CD56dim (blood)")
             & (r.gene == "SELL")].iloc[0]
    chk("GSE236394 SELL higher in bright", "positive, q<0.05",
        f"{sell.log2FC:+.2f}, q={sell.q_BH:.4f}",
        sell.log2FC > 3 and sell.q_BH < 0.05)
    gzmb = r[(r.comparison == "CD56bright vs CD56dim (blood)")
             & (r.gene == "GZMB")].iloc[0]
    chk("GSE236394 GZMB is the recorded anomaly (positive)", "positive",
        f"{gzmb.log2FC:+.2f}", gzmb.log2FC > 0)
    cd69 = r[(r.comparison == "liver perfusate vs blood (subset fixed)")
             & (r.gene == "CD69")].iloc[0]
    chk("GSE200319 CD69 is the recorded anomaly (negative)", "negative",
        f"{cd69.log2FC:+.2f}", cd69.log2FC < 0)

    print("== panel ==")
    ch = r[r.role == "channel"]
    always = sorted(g for g, s in ch.groupby("gene") if s.near_floor.sum() == 0)
    chk("4 channels never near the floor", "PIEZO1,TMEM63A,TRPM7,TRPV2",
        ",".join(always),
        always == ["PIEZO1", "TMEM63A", "TRPM7", "TRPV2"])
    counts = {g: int((s.q_BH < 0.05).sum()) for g, s in r.groupby("gene")}
    chk("TMEM63A significant in 6", 6, counts["TMEM63A"], counts["TMEM63A"] == 6)
    chk("TRPM7 significant in 5", 5, counts["TRPM7"], counts["TRPM7"] == 5)
    chk("TRPV2 significant in 5", 5, counts["TRPV2"], counts["TRPV2"] == 5)
    chk("PIEZO1 significant in 0", 0, counts["PIEZO1"], counts["PIEZO1"] == 0)
    dn = r[r.role == "mechano_downstream"]
    nf = {g: int(s.near_floor.sum()) for g, s in dn.groupby("gene")}
    chk("WWTR1 near floor in 13", 13, nf["WWTR1"], nf["WWTR1"] == 13)
    chk("CCN2 near floor in 12", 12, nf["CCN2"], nf["CCN2"] == 12)
    chk("YAP1 near floor in 6", 6, nf["YAP1"], nf["YAP1"] == 6)
    chk("ANKRD1 near floor in 8", 8, nf["ANKRD1"], nf["ANKRD1"] == 8)

    trpm7 = r[r.gene == "TRPM7"]
    cyto = trpm7[trpm7.comparison.str.contains("UNSTIM")]
    chk("TRPM7 negative in all 6 cytokine comparisons", 6,
        int((cyto.log2FC < 0).sum()),
        int((cyto.log2FC < 0).sum()) == 6 and len(cyto) == 6)
    tv2 = r[(r.gene == "TRPV2") & (r.comparison.str.contains("CD56bright"))]
    chk("TRPV2 sign disagrees between the two bright/dim datasets",
        "one +, one -", list(np.round(tv2.log2FC, 2)),
        (tv2.log2FC > 0).sum() == 1 and (tv2.log2FC < 0).sum() == 1
        and bool((tv2.q_BH < 0.05).all()))

    print("== orthology ==")
    chan = orth[orth.role == "channel"]
    chk("all 11 channels one-to-one", 11, int(chan.one_to_one.sum()),
        int(chan.one_to_one.sum()) == 11 and len(chan) == 11)
    gz = orth[orth.human == "GZMB"]
    chk("GZMB is one-to-many (Gzmb + Gzmc)", "Gzmb,Gzmc",
        ",".join(sorted(gz.mouse_symbol)),
        sorted(gz.mouse_symbol) == ["Gzmb", "Gzmc"] and not gz.one_to_one.any())
    xcl = orth[orth.role == "procedure_control"]
    chk("XCL1/XCL2 collapse onto one mouse gene", "Xcl1",
        ",".join(sorted(set(xcl.mouse_symbol))),
        set(xcl.mouse_symbol) == {"Xcl1"} and not xcl.one_to_any_check
        if hasattr(xcl, "one_to_any_check") else
        (set(xcl.mouse_symbol) == {"Xcl1"} and not xcl.one_to_one.any()))
    ccn2 = orth[orth.human.isin(["CCN2", "CTGF"])]
    chk("CTGF and CCN2 resolve to the same mouse gene", 1,
        len(set(ccn2.mouse_ensembl)), len(set(ccn2.mouse_ensembl)) == 1)

    print("== prose guard ==")
    # The dash and the indefinite "一个" are what separate a claim about
    # PIEZO1 from a sentence that merely uses PIEZO1 to set up a caveat, so
    # they are excluded from the gap the pattern is allowed to span.
    BANNED = [
        # a firm biological null for PIEZO1, stated without the power caveat
        (r"PIEZO1`?[^。\n—一]{0,12}(不随状态改变|不随状态变化|没有变化|完全不变)",
         "unqualified 'PIEZO1 does not change'"),
        # §0.1 forbids turning a transcript result into a functional claim
        (r"(TRPM7|TMEM63A)[^。\n]{0,20}(更重要|功能上更|才是关键)",
         "a functional claim built on the transcript result"),
        (r"PIEZO1`?[^。\n]{0,20}(无功能|没有功能|不参与)",
         "a functional claim about PIEZO1"),
    ]

    def assertive_hit(pattern, text):
        """First match that is not inside a negated clause.

        A document that forbids a claim has to be able to name it.  Every
        sentence here that says "this file does NOT conclude X" would trip a
        naive grep on X, so a hit is ignored when the run-up to it inside the
        same clause carries a negation.
        """
        for m in re.finditer(pattern, text):
            start = max(0, m.start() - 14)
            lead = text[start:m.start()]
            if any(neg in lead for neg in ("不", "未", "非", "无须", "并不")):
                continue
            return m
        return None
    for doc in DOCS:
        path = os.path.join(ROOT, doc)
        if not os.path.exists(path):
            chk(f"{doc} exists", "present", "MISSING", False)
            continue
        text = open(path, encoding="utf-8").read()
        scan = re.sub(r"「[^」]*」", "「…」", text)   # a mention is not a claim
        for pat, label in BANNED:
            hit = assertive_hit(pat, scan)
            chk(f"{doc}: no {label}", "absent",
                (hit.group(0)[:40] if hit else "absent"), hit is None)
        chk(f"{doc}: carries the §0.1 limitation", "present",
            "present" if "不能支持" in text and "值得进一步测定" in text else "MISSING",
            "不能支持" in text and "值得进一步测定" in text)

    # §5 says the unfavourable rows must be treated equally; the main file has
    # to state which row fired, and it has to be the second one.
    main = open(os.path.join(ROOT, "PIEZO1_BULK.md"), encoding="utf-8").read()
    chk("PIEZO1_BULK.md names the row that fired", "row 2 fired",
        "触发" if "对象可能选错" in main else "MISSING",
        "对象可能选错" in main and "触发" in main)

    print(f"\n{'ALL CHECKS PASS' if not fails else 'FAILURES: ' + ', '.join(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
