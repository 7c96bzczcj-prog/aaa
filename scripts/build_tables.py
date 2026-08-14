#!/usr/bin/env python3
"""Assemble the spec's output tables from the collected phase records.

  out/T1_evidence_audit.tsv  — Phase A, one row per claim (all 25)
  out/T2_decidability.tsv    — Phase B, one row per claim still INF_* after Phase A
  out/T5_shortlist.tsv       — the preregistered section 5 rule applied

Where the adversarial pass overturned the audit, the adversarial verdict wins
(it re-read the same methods with the explicit task of refuting the class), and
the change is recorded in `adversarial_revision` so nothing is silently replaced.
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(HERE, "phaseA")
OUT = os.path.join(HERE, "out")
INF = {"INF_INVITRO", "INF_MARKER", "INF_TRAJECTORY", "INF_KINETIC"}

T1_COLS = ["claim_id", "claim_text", "primary_source_pmid", "fulltext_available", "species",
           "in_vivo", "what_was_measured", "n", "n_unit", "controls",
           "evidence_type_strongest", "evidence_types_all",
           "strongest_alternative_explanation", "counter_evidence", "load_bearing",
           "downstream_dependency", "claude_initial_call", "audit_overturns_claude",
           "adversarial_revision", "crossfield_g4_findings", "primary_source_citation"]

T2_COLS = ["claim_id", "lineage_decidable", "why", "action_if_yes", "action_if_no",
           "both_answers_actionable", "compartments_needed", "compartments_obtainable",
           "already_traced", "occupancy_unverified", "mouse_only", "nearest_occupant",
           "crossfield_occupancy", "search_confidence", "search_terms_used"]

T5_COLS = ["claim_id", "claim_text", "load_bearing", "evidence_type_strongest",
           "lineage_decidable", "both_answers_actionable", "already_traced",
           "occupancy_unverified", "on_shortlist", "failed_criteria"]


def clean(v):
    if v is None:
        return "NA"
    if isinstance(v, (list, tuple)):
        v = "; ".join(str(x) for x in v)
    return " ".join(str(v).split()) or "NA"


def load_claims():
    return json.load(open(os.path.join(PA, "claims.json")))


def load_records():
    recs = {}
    for p in sorted(glob.glob(os.path.join(PA, "[A-Z][0-9].json"))):
        d = json.load(open(p))
        recs[d["claim_id"]] = d
    return recs


def order(cids):
    return sorted(cids, key=lambda c: (c[0], int(c[1:])))


def build_t1(recs, claims):
    with open(os.path.join(OUT, "T1_evidence_audit.tsv"), "w") as f:
        f.write("\t".join(T1_COLS) + "\n")
        for cid in order(recs):
            r = recs[cid]
            a = r.get("audit") or {}
            v = r.get("adv") or {}
            fin = r.get("final") or {}
            src = fin.get("resolution_source", "")
            if src.startswith("adversarial") or src == "adjudicated":
                rev = (f"{src.upper()}: {fin.get('audit_evidence_type')} -> "
                       f"{fin.get('evidence_type_strongest')}. {fin.get('resolution_note', '')}")
            else:
                rev = f"UPHELD. {v.get('challenge_reason', '')}"
            counter = []
            for tag, d in (("audit", a), ("adversarial", v)):
                t = d.get("counter_evidence") or d.get("new_counter_evidence")
                if t and str(t).upper() not in ("NOT_SEARCHED", "NONE", "NA", ""):
                    counter.append(f"[{tag}] {t}")
            if not counter:
                st = {a.get("counter_evidence_status"), v.get("counter_evidence_status")}
                counter = ["NOT_SEARCHED"] if "NOT_SEARCHED" in st else ["NONE (searched, none found)"]
            row = dict(a)
            row.update({
                "claim_id": cid,
                "claim_text": claims.get(cid, {}).get("text", "NA"),
                "claude_initial_call": claims.get(cid, {}).get("initial", "NA"),
                "evidence_type_strongest": fin.get("evidence_type_strongest", "NA"),
                "load_bearing": fin.get("load_bearing", "NA"),
                "adversarial_revision": rev,
                "counter_evidence": " || ".join(counter),
                "crossfield_g4_findings": v.get("crossfield_g4_findings", "NOT_SEARCHED"),
            })
            f.write("\t".join(clean(row.get(c)) for c in T1_COLS) + "\n")
    print(f"T1: {len(recs)} rows")


def build_t2(phaseb):
    with open(os.path.join(OUT, "T2_decidability.tsv"), "w") as f:
        f.write("\t".join(T2_COLS) + "\n")
        for cid in order(phaseb):
            d = phaseb[cid]
            b, o = d.get("b") or {}, d.get("o") or {}
            row = dict(b)
            row.update(o)
            row["claim_id"] = cid
            row["occupancy_unverified"] = ("TRUE" if o.get("already_traced") == "NOT_SEARCHED"
                                           else "FALSE")
            f.write("\t".join(clean(row.get(c)) for c in T2_COLS) + "\n")
    print(f"T2: {len(phaseb)} rows")


def build_t5(recs, claims, phaseb):
    """Apply the preregistered section 5 rule. No other rule is applied here."""
    rows = []
    for cid in order(recs):
        fin = recs[cid].get("final") or {}
        ev = fin.get("evidence_type_strongest")
        lb = fin.get("load_bearing")
        d = phaseb.get(cid) or {}
        b, o = d.get("b") or {}, d.get("o") or {}
        ld = b.get("lineage_decidable")
        ba = b.get("both_answers_actionable")
        at = o.get("already_traced")
        fails = []
        if lb != "HIGH":
            fails.append(f"criterion1 load_bearing={lb}")
        if ev not in INF:
            fails.append(f"criterion2 evidence={ev}")
        if cid in phaseb:
            if ld != "TRUE":
                fails.append(f"criterion3 lineage_decidable={ld}")
            if ba != "TRUE":
                fails.append(f"criterion3 both_answers_actionable={ba}")
            if at == "TRUE":
                fails.append("criterion4 already_traced=TRUE")
        elif ev in INF:
            fails.append("criterion3/4 Phase B not returned")
        rows.append({
            "claim_id": cid, "claim_text": claims.get(cid, {}).get("text"),
            "load_bearing": lb, "evidence_type_strongest": ev,
            "lineage_decidable": ld or "n/a", "both_answers_actionable": ba or "n/a",
            "already_traced": at or "n/a",
            "occupancy_unverified": "TRUE" if at == "NOT_SEARCHED" else "FALSE",
            "on_shortlist": "TRUE" if not fails else "FALSE",
            "failed_criteria": "; ".join(fails) or "none",
        })
    with open(os.path.join(OUT, "T5_shortlist.tsv"), "w") as f:
        f.write("\t".join(T5_COLS) + "\n")
        for r in rows:
            f.write("\t".join(clean(r.get(c)) for c in T5_COLS) + "\n")
    short = [r["claim_id"] for r in rows if r["on_shortlist"] == "TRUE"]
    unver = [r["claim_id"] for r in rows
             if r["on_shortlist"] == "TRUE" and r["occupancy_unverified"] == "TRUE"]
    print(f"T5: shortlist = {len(short)} {short}")
    if unver:
        print(f"    of which occupancy_unverified (already_traced=NOT_SEARCHED): {unver}")
    return short


def main():
    os.makedirs(OUT, exist_ok=True)
    claims, recs = load_claims(), load_records()
    p = os.path.join(PA, "phaseB.json")
    phaseb = json.load(open(p)) if os.path.exists(p) else {}
    build_t1(recs, claims)
    build_t2(phaseb)
    build_t5(recs, claims, phaseb)


if __name__ == "__main__":
    main()
