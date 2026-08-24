#!/usr/bin/env python3
"""Emit the stage-2 review pool (group A: NPC term + modality term in title/abstract)."""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")


def pubtypes(r):
    pl = (r.get("pubTypeList") or {}).get("pubType") or []
    return pl if isinstance(pl, list) else [pl]


def main():
    s = json.load(open(os.path.join(RAW, "survivors.json")))
    pool = [r for r in s if r["_flags"]["has_nasopharyngeal_term"]]
    pool.sort(key=lambda r: (r.get("pubYear") or "0", r.get("pmid") or ""))
    json.dump([r.get("pmid") or r.get("id") for r in pool], open(os.path.join(RAW, "groupA_ids.json"), "w"))
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else len(pool)
    for i, r in enumerate(pool[lo:hi], lo + 1):
        ab = re.sub(r"\s+", " ", (r.get("abstractText") or "")).strip()
        ab = re.sub(r"<[^>]+>", "", ab)
        pt = ",".join(pubtypes(r))
        ji = (r.get("journalInfo") or {}).get("journal") or {}
        jn = r.get("journalTitle") or ji.get("medlineAbbreviation") or ji.get("title") or "?"
        print("### [%d] pmid=%s doi=%s pmcid=%s src=%s" % (i, r.get("pmid"), r.get("doi"), r.get("pmcid"), r.get("source")))
        print("    %s %s | %s | q=%s" % (r.get("pubYear"), jn[:45], pt[:70], ",".join(r.get("_queries"))))
        print("    T: %s" % re.sub(r"<[^>]+>", "", (r.get("title") or ""))[:250])
        print("    A: %s" % (ab[:1150] if ab else "<NO ABSTRACT>"))
        print()


if __name__ == "__main__":
    main()
