#!/usr/bin/env python3
"""Check every quantitative claim in docs/RESULTS.md against the result files.

This repository has been bitten before by a README number that no longer matched
the file it came from, so the numbers are re-derived here rather than trusted.
Exits non-zero if any check fails.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
PA = os.path.join(HERE, "phaseA")
fails = []


def tsv(name):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(OUT, name))]
    head = rows[0]
    return [dict(zip(head, r)) for r in rows[1:]]


def check(label, got, want):
    ok = got == want
    print(f"  [{'OK ' if ok else 'FAIL'}] {label}: got {got!r}, expect {want!r}")
    if not ok:
        fails.append(label)


text = open(os.path.join(HERE, "docs", "RESULTS.md")).read()

print("== table shapes ==")
t1, t2, t5 = tsv("T1_evidence_audit.tsv"), tsv("T2_decidability.tsv"), tsv("T5_shortlist.tsv")
check("T1 rows (25 claims)", len(t1), 25)
check("T2 rows (INF_* claims)", len(t2), 17)
check("T5 rows", len(t5), 25)

print("\n== shortlist ==")
short = sorted(r["claim_id"] for r in t5 if r["on_shortlist"] == "TRUE")
check("shortlist membership", short, ["S2", "S3"])
check("RESULTS says shortlist is non-empty", "shortlist is non-empty" in text.lower(), True)
for c in short:
    check(f"{c} claimed in RESULTS §1", bool(re.search(rf"\*\*{c}\*\*", text)), True)

print("\n== shortlist entries satisfy all four preregistered criteria ==")
INF = {"INF_INVITRO", "INF_MARKER", "INF_TRAJECTORY", "INF_KINETIC"}
for r in t5:
    if r["on_shortlist"] != "TRUE":
        continue
    c = r["claim_id"]
    check(f"{c} criterion1 load_bearing", r["load_bearing"], "HIGH")
    check(f"{c} criterion2 evidence in INF_*", r["evidence_type_strongest"] in INF, True)
    check(f"{c} criterion3 decidable", r["lineage_decidable"], "TRUE")
    check(f"{c} criterion3 both actionable", r["both_answers_actionable"], "TRUE")
    check(f"{c} criterion4 not traced", r["already_traced"], "FALSE")

print("\n== X1 ==")
x1 = next(r for r in t5 if r["claim_id"] == "X1")
check("X1 not on shortlist", x1["on_shortlist"], "FALSE")
check("X1 fails on load-bearing only", x1["failed_criteria"], "criterion1 load_bearing=MEDIUM")
check("X1 evidence", x1["evidence_type_strongest"], "INF_MARKER")
check("X1 decidable", x1["lineage_decidable"], "TRUE")
check("X1 both actionable", x1["both_answers_actionable"], "TRUE")
check("X1 already_traced FALSE", x1["already_traced"], "FALSE")
check("X1 mouse_only", x1["mouse_only"], "TRUE")

print("\n== Phase C power (GSE302113) ==")
lib = tsv("T4_power_GSE302113_libraries.tsv")
check("libraries profiled", len(lib), 38)
covs = [float(r["median_chrM_coverage"]) for r in lib]
check("all libraries >= 20x chrM", all(c >= 20 for c in covs), True)
check("min median chrM coverage", round(min(covs), 1), 23.8)
check("max median chrM coverage", round(max(covs), 1), 158.9)
check("chrM present in all libraries", all(r["chrM_data_present"] == "TRUE" for r in lib), True)

pw = tsv("T4_power.tsv")
lung = sorted({r["donor_id"] for r in pw if r["donor_id"].startswith("SU-L")})
x1_ok = sorted(d for d in lung
               if all(next((x for x in pw if x["donor_id"] == d and x["compartment"] == comp),
                           {"nk_pass": "FALSE"})["nk_pass"] == "TRUE"
                      for comp in ("adjacent_normal", "tumour")))
check("X1 donors passing >=30 NK in both compartments", x1_ok,
      ["SU-L-001", "SU-L-002", "SU-L-004", "SU-L-005"])
check("X1 power verdict PASS (>=3 donors)", len(x1_ok) >= 3, True)

print("\n== classifier acceptance control ==")
s = json.load(open(os.path.join(HERE, "phaseC", "calls", "summary.json")))
pos = s["GSM9096509"]["nk"] / s["GSM9096509"]["n_typed"]
neg = s["GSM9096510"]["nk"] / s["GSM9096510"]["n_typed"]
check("CD56+ NK fraction (%)", round(pos * 100, 1), 74.6)
check("CD56- NK fraction (%)", round(neg * 100, 1), 5.5)
check("ratio", round(pos / neg, 1), 13.4)
check("acceptance test passes", pos > 0.5 and pos > 3 * neg, True)

print("\n== Phase A completeness ==")
fin = json.load(open(os.path.join(PA, "final_classes.json")))
check("final classes for all 25", len(fin), 25)
pb = json.load(open(os.path.join(PA, "phaseB.json")))
check("Phase B decidability records", sum(1 for v in pb.values() if v.get("b")), 17)
check("Phase B occupancy records", sum(1 for v in pb.values() if v.get("o")), 17)
inf_claims = sorted(c for c, f in fin.items() if f["evidence_type_strongest"] in INF)
check("Phase B covers exactly the INF_* claims", sorted(pb) == inf_claims, True)
n_changed = sum(1 for f in fin.values()
                if f["evidence_type_strongest"] != f["audit_evidence_type"])
check("adversarial changed 9 of 25 classifications", n_changed, 9)

print("\n== full-text availability (early-termination gate) ==")
gate = [r for r in t1 if r["claim_id"] in ("D1", "D2", "D3", "D4", "P1", "P2")]
check("gate claims audited", len(gate), 6)
check("gate FALSE count (rule: halt if >3)",
      sum(1 for r in gate if r["fulltext_available"] == "FALSE"), 0)
false_rows = sorted(r["claim_id"] for r in t1 if r["fulltext_available"] == "FALSE")
check("claims with fulltext FALSE overall", false_rows, ["P3"])
check("claims with fulltext TRUE", sum(1 for r in t1 if r["fulltext_available"] == "TRUE"), 10)
check("claims with fulltext PARTIAL", sum(1 for r in t1 if r["fulltext_available"] == "PARTIAL"), 14)

print()
if fails:
    print(f"FAILED {len(fails)} check(s): {fails}")
    sys.exit(1)
print("ALL CHECKS PASS")


