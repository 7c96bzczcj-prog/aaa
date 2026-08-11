"""Check every quantitative claim in the prose docs against the result files.

Written after an adversarial audit found a headline claim resting on a
guard that was never armed. Numbers in prose drift from numbers in files
as an analysis is corrected; this makes the drift a test failure rather
than something a reader has to catch.

Run it after any re-analysis: `python3 scripts/verify_claims.py`
"""
from __future__ import annotations

import os
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")


def main() -> int:
    checks = []

    def chk(label, claimed, actual, ok):
        checks.append((bool(ok), label, str(claimed), str(actual)))

    q = pd.read_csv(os.path.join(OUT, "quadrants_relaxed_d0.5_p95.csv"))
    vc = q.quadrant.value_counts().to_dict()
    for name, want in (("Q4", 46), ("Q4_attenuated", 7), ("Q3", 43),
                       ("Q3_NK_indeterminate", 10), ("Q1", 10)):
        chk(f"quadrant {name}", want, vc.get(name), vc.get(name) == want)
    chk("genes tested", 4525, len(q), len(q) == 4525)

    a4 = pd.read_csv(os.path.join(OUT, "A4_cross_dataset_scan.csv"))
    for acc, thr, col, want in (("154826", 30, "all_5_lineages", 27),
                                ("178341", 30, "NK_only", 5),
                                ("131907", 20, "all_5_lineages", 8)):
        row = a4[(a4.dataset.str.contains(acc)) & (a4.threshold == thr)].iloc[0]
        chk(f"A4 {acc} @{thr} {col}", want, row[col], row[col] == want)

    p6 = pd.read_csv(os.path.join(OUT, "phase6_shared_soup_contrast.csv"))
    sh = float(p6[p6.group == "shared_soup"].mean_pairwise_r.iloc[0])
    mt = float(p6[p6.group == "separate_matched"].mean_pairwise_r.mean())
    chk("phase6 shared r", 0.291, round(sh, 3), abs(sh - 0.291) < 0.002)
    chk("phase6 separate r", 0.490, round(mt, 3), abs(mt - 0.490) < 0.002)

    rr = pd.read_csv(os.path.join(OUT, "replication_rates_by_quadrant.csv"))
    got = dict(zip(rr.quadrant, rr.n_replicated))
    tot = dict(zip(rr.quadrant, rr.n))
    chk("phase7 Q4", "0/46", f"{got.get('Q4')}/{tot.get('Q4')}",
        got.get("Q4") == 0 and tot.get("Q4") == 46)
    chk("phase7 Q3", "0/43", f"{got.get('Q3')}/{tot.get('Q3')}",
        got.get("Q3") == 0 and tot.get("Q3") == 43)

    nul = pd.read_csv(os.path.join(OUT, "q4_null_rate.csv"))
    chk("null Q4 max", 0, int(nul.Q4.max()), nul.Q4.max() == 0)

    rc = pd.read_csv(os.path.join(OUT, "q4_synthetic_recall.csv"))
    r0 = float(rc[rc.construction == "rescale lambda=0.0"].recall_Q4.iloc[0])
    chk("synthetic recall lambda=0", 0.760, round(r0, 3), abs(r0 - 0.760) < 0.005)
    perm = float(rc[rc.construction == "permute (lambda=0)"].recall_Q4.iloc[0])
    chk("permutation construction", 0.0, round(perm, 3), perm == 0.0)
    kimv = open(os.path.join(OUT, "q4_synthetic_recall_kim_verdict.txt")).read()
    chk("kim source genes", 6, kimv.split("n_source\t")[1].split("\n")[0],
        "n_source\t6" in kimv)

    pw = pd.read_csv(os.path.join(OUT, "power_check.csv"))
    chk("NK median SE", 0.113, round(float(pw.NK.iloc[0]), 3),
        abs(float(pw.NK.iloc[0]) - 0.113) < 0.002)

    # --- the two authorised ambient follow-up runs -----------------------
    # Run 1: the pre-registered branch turns on whether pooling rho brings
    # NK's SE back to baseline.  It does not, and that is the claim.
    se = pd.read_csv(os.path.join(OUT, "rho_variant_se.csv"))
    sev = se.set_index(["arm", "lineage"]).median_SE
    for arm, lin, want in [("baseline", "NK", 0.120), ("rho_per_sample", "NK", 0.199),
                           ("rho_pooled", "NK", 0.193), ("rho_pooled", "B", 0.270),
                           ("rho_wide_markers", "B", 0.516)]:
        got = float(sev[(arm, lin)])
        chk(f"SE {arm}/{lin}", want, round(got, 3), abs(got - want) < 0.002)
    chk("pooling rho does not restore baseline SE", "True",
        str(float(sev[("rho_pooled", "NK")]) > 1.5 * float(sev[("baseline", "NK")])),
        float(sev[("rho_pooled", "NK")]) > 1.5 * float(sev[("baseline", "NK")]))

    rv = pd.read_csv(os.path.join(OUT, "rho_variants.csv"))
    sd = (rv.groupby(["batch", "lineage"]).rho_per_sample.std()
            .groupby("lineage").median())
    chk("NK between-condition SD of rho", 0.120, round(float(sd["NK"]), 3),
        abs(float(sd["NK"]) - 0.120) < 0.002)
    chk("CD4T between-condition SD of rho", 0.005, round(float(sd["CD4T"]), 3),
        abs(float(sd["CD4T"]) - 0.005) < 0.002)

    # Run 2: the whole claim is that the correction tracks the estimated
    # soup fraction -- it works where f is right and not where f is wrong.
    ar = pd.read_csv(os.path.join(OUT, "ambient_regression_summary.csv")
                     ).set_index("lineage")
    for lin, f_want, drop_want in [("B", 1.000, 0.77), ("NK", 0.291, 0.09),
                                   ("Myeloid", 0.767, 0.12)]:
        f_got = float(ar.soup_frac_zero_truth[lin])
        drop = 1 - float(ar.after_median_abs[lin]) / float(ar.before_median_abs[lin])
        chk(f"soup fraction {lin}", f_want, round(f_got, 3),
            abs(f_got - f_want) < 0.002)
        chk(f"control reduction {lin}", drop_want, round(drop, 2),
            abs(drop - drop_want) < 0.01)
    for lin in ("CD8T", "CD4T"):
        drop = 1 - float(ar.after_median_abs[lin]) / float(ar.before_median_abs[lin])
        chk(f"control got worse in {lin}", "<0", f"{drop:+.2f}", drop < 0)

    # --- the complosome premise check (docs/LEADS_CLOSED.md) -------------
    # The whole verdict turns on C3 clearing an ambient calibration that is
    # measured inside NK rather than assumed, so both ends of that scale are
    # checked, not just C3's number.
    cp = pd.read_csv(os.path.join(OUT, "complosome_premise.csv"))
    nk = cp[cp.lineage == "NK"].set_index("gene")
    my = cp[cp.lineage == "Myeloid"].set_index("gene")
    for g, want in [("C3", 0.191), ("C3AR1", 0.388), ("IGKC", 0.402)]:
        got = float(nk.soup_frac[g])
        chk(f"NK soup_frac {g}", want, round(got, 3), abs(got - want) < 0.002)
    ceiling = float(nk.loc[["C1QA", "LYZ", "IGKC"], "soup_frac"].min())
    floor = float(nk.loc[["KLRD1", "NKG7", "GNLY"], "soup_frac"].max())
    chk("ambient calibration ceiling in NK", 0.402, round(ceiling, 3),
        abs(ceiling - 0.402) < 0.002)
    chk("genuine-NK calibration floor", 0.013, round(floor, 3),
        abs(floor - 0.013) < 0.002)
    chk("C3 below the fully-ambient reading", "True",
        str(float(nk.soup_frac["C3"]) < ceiling),
        float(nk.soup_frac["C3"]) < ceiling)
    ratio = float(nk.mean_CPM["C3"] / my.mean_CPM["C3"])
    band = max(float(nk.mean_CPM[g] / my.mean_CPM[g]) for g in ("C1QA", "LYZ"))
    chk("C3 NK/Myeloid CPM", 0.437, round(ratio, 3), abs(ratio - 0.437) < 0.002)
    chk("pure-pickup band max", 0.118, round(band, 3), abs(band - 0.118) < 0.002)
    chk("C3 above the pickup band", "True", str(ratio > band), ratio > band)
    chk("C3 detected in NK", 0.690, round(float(nk.detection["C3"]), 3),
        abs(float(nk.detection["C3"]) - 0.690) < 0.002)

    # over-dispersion death certificate: NK's rho vs CD8T's
    rr = float(rv.groupby("lineage").rho_per_sample.median()["NK"] /
               rv.groupby("lineage").rho_per_sample.median()["CD8T"])
    chk("NK rho / CD8T rho", 2.6, round(rr, 1), abs(rr - 2.6) < 0.05)

    # A retracted claim must stay retracted EVERYWHERE, not just in the
    # README.  The first version of this check scanned README.md alone and
    # reported "16/16 verified" while the retracted sentence was still live,
    # unqualified and in bold, in docs/DEVIATIONS.md D7 -- the very document
    # the README sends readers to for the argument.  An audit caught it, not
    # this checker.  Scan every prose file, and require any surviving mention
    # to sit inside an explicit retraction/supersession marker.
    RETRACTED = ("demonstrated Q4 power", "S3 and S4 both passing",
                 "cleared to produce candidates")
    MARKERS = ("retract", "SUPERSEDED", "superseded", "earlier revision",
               "earlier version")
    prose = [os.path.join(ROOT, "README.md")]
    docs = os.path.join(ROOT, "docs")
    if os.path.isdir(docs):
        prose += [os.path.join(docs, f) for f in sorted(os.listdir(docs))
                  if f.endswith(".md")]

    for path in prose:
        text = open(path).read()
        name = os.path.relpath(path, ROOT)
        for phrase in RETRACTED:
            bad = []
            for para in text.split("\n\n"):
                if phrase in para and not any(m in para for m in MARKERS):
                    bad.append(para.strip().split("\n")[0][:60])
            chk(f"retracted phrase in {name}",
                f"absent or marked", "OK" if not bad else f"UNMARKED: {bad[0]}",
                not bad)

    bad = [c for c in checks if not c[0]]
    for ok, label, claimed, actual in checks:
        print(f"  {'OK      ' if ok else 'MISMATCH'} {label:28s} "
              f"README={claimed:8s} file={actual}")
    print(f"\n{len(checks) - len(bad)}/{len(checks)} claims verified")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
