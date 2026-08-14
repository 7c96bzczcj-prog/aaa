#!/usr/bin/env python3
"""Collect Phase A audit + adversarial records from the workflow journals into
phaseA/<claim>.json, one file per claim.

Reads the journals rather than the workflow return values so that partial runs
are recoverable and so every row is traceable to the agent that produced it.
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(HERE, "phaseA")
BASE = ("/root/.claude/projects/-home-user-aaa/"
        "5f46e874-56b9-5c6b-9d27-d08dbfb1d5d8/subagents/workflows")
WFS = {
    "wf_48c4ecfe-cc6": "gate",
    "wf_0004b694-565": "tissue",
    "wf_4b37b92d-b2b": "tumour",
    "wf_041a111c-558": "reversal",
}


def main():
    os.makedirs(PA, exist_ok=True)
    audits, advs = {}, {}
    for wf, name in WFS.items():
        p = os.path.join(BASE, wf, "journal.jsonl")
        if not os.path.exists(p):
            continue
        for line in open(p):
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if d.get("type") != "result":
                continue
            r = d.get("result")
            if not isinstance(r, dict):
                continue
            cid = r.get("claim_id")
            if not cid:
                continue
            r["_batch"] = name
            if "classification_challenged" in r:
                advs[cid] = r
            elif "evidence_type_strongest" in r:
                audits[cid] = r

    claims = json.load(open(os.path.join(PA, "claims.json")))
    for cid in sorted(set(audits) | set(advs)):
        rec = {"claim_id": cid,
               "claim_text": claims.get(cid, {}).get("text"),
               "claude_initial_call": claims.get(cid, {}).get("initial"),
               "audit": audits.get(cid),
               "adv": advs.get(cid)}
        json.dump(rec, open(os.path.join(PA, f"{cid}.json"), "w"), indent=1)

    missing_a = sorted(set(claims) - set(audits))
    missing_v = sorted(set(claims) - set(advs))
    print(f"audits={len(audits)}/25  adversarial={len(advs)}/25")
    if missing_a:
        print(f"MISSING AUDIT: {missing_a}")
    if missing_v:
        print(f"MISSING ADVERSARIAL: {missing_v}")
    return 0 if not missing_a else 1


if __name__ == "__main__":
    sys.exit(main())
