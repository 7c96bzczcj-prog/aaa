#!/usr/bin/env python3
"""Collect Phase B decidability + occupancy records into phaseA/phaseB.json."""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(HERE, "phaseA")
BASE = ("/root/.claude/projects/-home-user-aaa/"
        "5f46e874-56b9-5c6b-9d27-d08dbfb1d5d8/subagents/workflows")
WFS = ["wf_db7eeea7-41f", "wf_77de77c1-0fc"]
EXPECTED = ["D3", "P1", "P2", "R1", "R3", "R4", "S1", "S2", "S3",
            "T1", "T3", "T5", "T6", "X1", "X3", "X4", "X5"]


def main():
    dec, occ = {}, {}
    for wf in WFS:
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
            if "already_traced" in r:
                occ[cid] = r
            elif "lineage_decidable" in r:
                dec[cid] = r

    out = {cid: {"b": dec.get(cid), "o": occ.get(cid)} for cid in sorted(set(dec) | set(occ))}
    json.dump(out, open(os.path.join(PA, "phaseB.json"), "w"), indent=1)
    md = [c for c in EXPECTED if c not in dec]
    mo = [c for c in EXPECTED if c not in occ]
    print(f"decidability={len(dec)}/{len(EXPECTED)}  occupancy={len(occ)}/{len(EXPECTED)}")
    if md:
        print(f"MISSING decidability: {md}")
    if mo:
        print(f"MISSING occupancy: {mo}")
    return 0 if not (md or mo) else 1


if __name__ == "__main__":
    sys.exit(main())
