#!/usr/bin/env python3
"""Pull the sentences that state where each candidate's SINGLE-CELL data came from.

Include vs secondary_analysis turns on exactly one question -- did this study
generate single-cell data from its own NPC samples, or only re-analyse deposited
data? That is a Methods/Data-availability statement, so quote it rather than infer.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftlib  # noqa: E402

SC = r"(single[\s\-]?cell|single[\s\-]?nucle|scRNA|snRNA|CyTOF|mass cytometry|CITE[\s\-]?seq|spatial transcriptom|10x|10×)"
SRC = (r"(downloaded|obtained|retrieved|acquired|publicly|public\s+dat|GEO\s+database|"
       r"GSE\d+|HRA\d+|EGA|OEP\d+|CRA\d+|we\s+(?:collected|performed|enrolled|recruited|generated|profiled)|"
       r"were\s+collected|were\s+enrolled|fresh|informed\s+consent|this\s+study|our\s+(?:cohort|study|centre|center)|"
       r"deposited|available\s+(?:at|in|from)|accession)")


def main():
    ids = json.load(open(os.path.join(RAW, "candidate_ids.json")))
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else len(ids)
    for rid in ids[lo:hi]:
        r = recs[rid]
        print("=" * 108)
        print("%s | pmid=%s | doi=%s | %s" % (rid, r.get("pmid"), r.get("doi"), r.get("pubYear")))
        print("T: %s" % re.sub(r"<[^>]+>", "", r.get("title") or "")[:150])
        xml = ftlib.load_xml(r["pmcid"]) if r.get("pmcid") else None
        if not xml:
            ab = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", r.get("abstractText") or ""))
            print("  [NO FULL TEXT] abstract:")
            print("   ", ab[:900])
            print()
            continue
        p = ftlib.parse(xml)
        secs = ftlib.pick_sections(p["sections"])
        hits = []
        for s in ftlib.sentences_with(p["plain"], SC, maxn=200):
            if re.search(SRC, s, re.I):
                hits.append(s)
        seen = set()
        n = 0
        for s in hits:
            k = s[:70]
            if k in seen:
                continue
            seen.add(k)
            print("   >", re.sub(r"\s+", " ", s)[:330])
            n += 1
            if n >= 7:
                break
        da = secs.get("data_availability")
        if da:
            print("   [DATA AVAIL]", re.sub(r"\s+", " ", da)[:420])
        print()


if __name__ == "__main__":
    main()
