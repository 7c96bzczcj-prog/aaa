#!/usr/bin/env python3
"""Normalise the adversarial pass's evidence-class verdicts to bare enum codes.

The adversarial schema left `revised_evidence_type_strongest` as a free string,
so several reviewers returned a reasoned paragraph rather than a token. The
content is usable; the field is not machine-readable. This module extracts the
leading class token where the reviewer supplied one, and applies an explicit,
recorded adjudication where the reviewer's answer was conditional.

Every non-mechanical call is listed in ADJUDICATED below with its reason and is
mirrored in docs/DECISIONS.md. Nothing here is decided by string-matching alone.
"""
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(HERE, "phaseA")
CLASSES = ["OBS_DIRECT", "OBS_TRANSFER", "OBS_GENETIC",
           "INF_INVITRO", "INF_MARKER", "INF_TRAJECTORY", "INF_KINETIC"]

# Claims where the reviewer's verdict was conditional on a reading of the claim,
# so a token could not simply be lifted. Resolved under DEC-05: the claim text in
# spec section 1 is authoritative for what is being audited.
ADJUDICATED = {
    "P1": ("INF_KINETIC",
           "The reviewer gave two answers conditional on scope: INF_KINETIC for the human "
           "NKG2C/HCMV claim, OBS_DIRECT for the mouse Ly49H/MCMV homolog (Flommersfeld 2021 "
           "single-cell colour-barcode fate mapping). The claim as written in spec section 1 is "
           "'cNK -> adaptive NKG2C+ NK, HCMV-driven'. NKG2C and HCMV are human-specific; the "
           "mouse homolog uses a different receptor (Ly49H), a different virus (MCMV) and a "
           "different species. Under DEC-05 that is a neighbouring but different transition, so "
           "it is recorded as adjacent evidence and does not set this row's class. "
           "Human strongest = INF_KINETIC."),
}


def leading_class(text):
    """Return the class token a reviewer led with, if they led with one."""
    if not text:
        return None
    t = text.strip().strip('*"—- ')
    for c in CLASSES:
        if t.upper().startswith(c):
            return c
    return None


def resolve(rec):
    """Return (final_class, source, note)."""
    a = rec.get("audit") or {}
    v = rec.get("adv") or {}
    cid = rec.get("claim_id")
    audited = a.get("evidence_type_strongest")
    if not v:
        return audited, "audit_only", "adversarial pass not returned"
    ch = (v.get("classification_challenged") or "UPHELD").upper()
    if ch == "UPHELD":
        return audited, "upheld", ""
    if cid in ADJUDICATED:
        cls, why = ADJUDICATED[cid]
        return cls, "adjudicated", why
    tok = leading_class(v.get("revised_evidence_type_strongest"))
    if tok:
        return tok, f"adversarial_{ch.lower()}", (v.get("revised_evidence_type_strongest") or "")[:400]
    return audited, "unparsed_kept_audit", (
        "reviewer challenged the class but supplied no leading token; audit class kept "
        "rather than guessed")


def main():
    out = {}
    for p in sorted(glob.glob(os.path.join(PA, "[A-Z][0-9].json"))):
        rec = json.load(open(p))
        cls, src, note = resolve(rec)
        a = rec.get("audit") or {}
        v = rec.get("adv") or {}
        lbc = (v.get("load_bearing_challenged") or "UPHELD").upper()
        lb = (v.get("load_bearing_revised") if lbc in ("DOWNGRADE", "UPGRADE")
              and v.get("load_bearing_revised") else a.get("load_bearing"))
        rec["final"] = {"evidence_type_strongest": cls, "resolution_source": src,
                        "resolution_note": note, "load_bearing": lb,
                        "audit_evidence_type": a.get("evidence_type_strongest")}
        json.dump(rec, open(p, "w"), indent=1)
        out[rec["claim_id"]] = rec["final"]

    INF = {"INF_INVITRO", "INF_MARKER", "INF_TRAJECTORY", "INF_KINETIC"}
    print(f"{'id':4s} {'audit':15s} {'final':15s} {'LB':7s} {'source':22s} shortlist-eligible")
    elig = []
    for cid in sorted(out, key=lambda c: (c[0], int(c[1:]))):
        f = out[cid]
        ok = f["evidence_type_strongest"] in INF and f["load_bearing"] == "HIGH"
        if ok:
            elig.append(cid)
        print(f"{cid:4s} {str(f['audit_evidence_type']):15s} {str(f['evidence_type_strongest']):15s} "
              f"{str(f['load_bearing']):7s} {f['resolution_source']:22s} {'YES' if ok else ''}")
    print(f"\nEligible for Phase B (INF_* AND load_bearing=HIGH): {len(elig)} -> {elig}")
    json.dump(out, open(os.path.join(PA, "final_classes.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
