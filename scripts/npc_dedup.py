#!/usr/bin/env python3
"""Dedupe the union of Q1..Q5 by PMID -> DOI -> lowercased title (spec v1.0 s.1)."""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")


def norm_title(t):
    t = (t or "").lower()
    t = re.sub(r"<[^>]+>", " ", t)
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return re.sub(r"\s+", " ", t).strip()


def load():
    recs = []
    for qid in ["Q1", "Q2", "Q3", "Q4", "Q5", "Q7", "Q8"]:
        p = os.path.join(RAW, "%s_results.json" % qid)
        d = json.load(open(p))
        for r in d["results"]:
            r = dict(r)
            r["_query"] = qid
            recs.append(r)
    return recs


def main():
    recs = load()
    retrieved = len(recs)
    by_key = {}
    order = []
    dup = 0
    for r in recs:
        pmid = (r.get("pmid") or "").strip()
        doi = (r.get("doi") or "").strip().lower()
        ti = norm_title(r.get("title"))
        key = None
        for cand in (("pmid", pmid), ("doi", doi), ("title", ti)):
            if cand[1]:
                k = "%s:%s" % cand
                if k in by_key:
                    key = k
                    break
        if key is None:
            # new record: register under every identifier it has
            primary = ("pmid:" + pmid) if pmid else (("doi:" + doi) if doi else ("title:" + ti))
            r["_queries"] = [r.pop("_query")]
            by_key_entry = r
            for cand in (("pmid", pmid), ("doi", doi), ("title", ti)):
                if cand[1]:
                    by_key["%s:%s" % cand] = by_key_entry
            order.append(primary)
        else:
            dup += 1
            tgt = by_key[key]
            q = r.pop("_query")
            if q not in tgt["_queries"]:
                tgt["_queries"].append(q)
            # backfill identifiers seen on the duplicate
            for cand in (("pmid", pmid), ("doi", doi), ("title", ti)):
                if cand[1] and ("%s:%s" % cand) not in by_key:
                    by_key["%s:%s" % cand] = tgt

    deduped = [by_key[k] for k in order]
    assert retrieved - dup == len(deduped), (retrieved, dup, len(deduped))
    with open(os.path.join(RAW, "deduped.json"), "w") as f:
        json.dump(deduped, f)
    counts = {"retrieved": retrieved, "duplicates": dup, "deduped": len(deduped)}
    json.dump(counts, open(os.path.join(RAW, "dedup_counts.json"), "w"), indent=2)
    print(json.dumps(counts, indent=2))
    # per-query contribution
    from collections import Counter
    c = Counter()
    for r in deduped:
        c["+".join(sorted(r["_queries"]))] += 1
    print(json.dumps(dict(c), indent=2))


if __name__ == "__main__":
    main()
