"""Check every quantitative claim in README.md against the result files.

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

    pw = pd.read_csv(os.path.join(OUT, "power_check.csv"))
    chk("NK median SE", 0.113, round(float(pw.NK.iloc[0]), 3),
        abs(float(pw.NK.iloc[0]) - 0.113) < 0.002)

    # the retracted claim must stay retracted
    readme = open(os.path.join(ROOT, "README.md")).read()
    chk("TOX-power claim absent", "absent",
        "absent" if "demonstrated Q4 power" not in readme else "PRESENT",
        "demonstrated Q4 power" not in readme)

    bad = [c for c in checks if not c[0]]
    for ok, label, claimed, actual in checks:
        print(f"  {'OK      ' if ok else 'MISMATCH'} {label:28s} "
              f"README={claimed:8s} file={actual}")
    print(f"\n{len(checks) - len(bad)}/{len(checks)} claims verified")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
