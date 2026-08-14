#!/usr/bin/env python3
"""Assemble out/T1_evidence_audit.tsv and out/T2_decidability.tsv from the
per-claim audit records in phaseA/*.json.

Each phaseA record is {"audit": {...}, "adv": {...}} — the initial audit and the
independent adversarial pass. Where the adversarial pass overturns the audit,
the adversarial verdict wins and the change is recorded in
`adversarial_revision`, because the adversarial agent re-read the same methods
with the explicit task of refuting the classification.
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(HERE, "phaseA")
OUT = os.path.join(HERE, "out")

T1_COLS = ["claim_id", "claim_text", "primary_source_pmid", "fulltext_available", "species",
           "in_vivo", "what_was_measured", "n", "n_unit", "controls",
           "evidence_type_strongest", "evidence_types_all",
           "strongest_alternative_explanation", "counter_evidence", "load_bearing",
           "downstream_dependency", "claude_initial_call", "audit_overturns_claude",
           "adversarial_revision", "crossfield_g4_findings", "primary_source_citation"]

T2_COLS = ["claim_id", "lineage_decidable", "why", "action_if_yes", "action_if_no",
           "both_answers_actionable", "compartments_needed", "compartments_obtainable",
           "already_traced", "occupancy_unverified", "search_terms_used"]

INF = {"INF_INVITRO", "INF_MARKER", "INF_TRAJECTORY", "INF_KINETIC"}


def clean(v):
    if v is None:
        return "NA"
    if isinstance(v, (list, tuple)):
        v = "; ".join(str(x) for x in v)
    return str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ").strip() or "NA"


def load():
    recs = {}
    for p in sorted(glob.glob(os.path.join(PA, "*.json"))):
        d = json.load(open(p))
        cid = (d.get("audit") or {}).get("claim_id") or os.path.basename(p)[:-5]
        recs[cid] = d
    return recs


def merged_evidence(r):
    """Adversarial verdict wins when it overturns."""
    a = r.get("audit") or {}
    v = r.get("adv") or {}
    ev = a.get("evidence_type_strongest", "NA")
    ch = (v.get("classification_challenged") or "UPHELD").upper()
    rev = v.get("revised_evidence_type_strongest")
    if ch in ("DOWNGRADE", "UPGRADE") and rev and rev not in ("", "NA", "NONE"):
        return rev, f"{ch}: {ev} -> {rev}. {clean(v.get('challenge_reason'))}"
    return ev, f"UPHELD ({clean(v.get('challenge_reason'))[:300]})"


def merged_load_bearing(r):
    a = r.get("audit") or {}
    v = r.get("adv") or {}
    ch = (v.get("load_bearing_challenged") or "UPHELD").upper()
    if ch in ("DOWNGRADE", "UPGRADE") and v.get("load_bearing_revised"):
        return v["load_bearing_revised"]
    return a.get("load_bearing", "NA")


def merged_counter(r):
    a = r.get("audit") or {}
    v = r.get("adv") or {}
    parts = []
    if a.get("counter_evidence_status") == "NOT_SEARCHED" and v.get("counter_evidence_status") == "NOT_SEARCHED":
        return "NOT_SEARCHED"
    for src, d in (("audit", a), ("adversarial", v)):
        t = d.get("counter_evidence") or d.get("new_counter_evidence")
        if t and t.upper() not in ("NOT_SEARCHED", "NONE", "NA", ""):
            parts.append(f"[{src}] {clean(t)}")
    if not parts:
        return "NONE (searched, nothing found)"
    return " || ".join(parts)


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = load()
    claims = json.load(open(os.path.join(HERE, "phaseA", "claims.json")))

    with open(os.path.join(OUT, "T1_evidence_audit.tsv"), "w") as f:
        f.write("\t".join(T1_COLS) + "\n")
        for cid in sorted(recs, key=lambda c: (c[0], int(c[1:]) if c[1:].isdigit() else 0)):
            r = recs[cid]
            a = r.get("audit") or {}
            v = r.get("adv") or {}
            ev, rev = merged_evidence(r)
            row = {
                "claim_id": cid,
                "claim_text": claims.get(cid, {}).get("text", "NA"),
                "claude_initial_call": claims.get(cid, {}).get("initial", "NA"),
                "evidence_type_strongest": ev,
                "adversarial_revision": rev,
                "load_bearing": merged_load_bearing(r),
                "counter_evidence": merged_counter(r),
                "crossfield_g4_findings": v.get("crossfield_g4_findings", "NOT_SEARCHED"),
            }
            for c in T1_COLS:
                row.setdefault(c, a.get(c, "NA"))
            f.write("\t".join(clean(row[c]) for c in T1_COLS) + "\n")

    # T2 only exists for claims whose strongest evidence is still INF_*
    dec = {}
    p2 = os.path.join(PA, "decidability.json")
    if os.path.exists(p2):
        dec = json.load(p2 and open(p2))
    with open(os.path.join(OUT, "T2_decidability.tsv"), "w") as f:
        f.write("\t".join(T2_COLS) + "\n")
        for cid in sorted(dec, key=lambda c: (c[0], int(c[1:]) if c[1:].isdigit() else 0)):
            d = dec[cid]
            d.setdefault("claim_id", cid)
            d["occupancy_unverified"] = "TRUE" if d.get("already_traced") == "NOT_SEARCHED" else "FALSE"
            f.write("\t".join(clean(d.get(c)) for c in T2_COLS) + "\n")

    print(f"T1: {len(recs)} rows   T2: {len(dec)} rows")
    inf = [c for c in recs if merged_evidence(recs[c])[0] in INF]
    print(f"claims whose strongest evidence is INF_* (eligible for Phase B): {len(inf)} {sorted(inf)}")


if __name__ == "__main__":
    main()
