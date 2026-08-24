#!/usr/bin/env python3
"""Recover accession numbers via the Europe PMC annotations API.

Records with no OA full text still have accession numbers text-mined from whatever
text Europe PMC holds. This is API-retrieved evidence, so it can legitimately fill
the `accession` field where the abstract alone says nothing. Note the mined set can
include accessions the article merely *cites*, so every hit is kept with its
surrounding sentence for adjudication rather than trusted blindly.
"""
import json, os, sys, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_fetch import http_get  # noqa: E402

API = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=%s&type=Accession%%20Numbers&format=JSON"


def main():
    ids = json.load(open(os.path.join(RAW, "candidate_ids.json")))
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    out = {}
    for rid in ids:
        r = recs[rid]
        src = r.get("source")
        ext = r.get("pmid") or r.get("id")
        if not src or not ext:
            continue
        aid = "%s:%s" % (src, ext)
        body = http_get(API % urllib.parse.quote(aid), retries=2, timeout=60)
        if not body:
            out[rid] = {"ok": False}
            continue
        try:
            d = json.loads(body)
        except Exception:  # noqa: BLE001
            out[rid] = {"ok": False}
            continue
        hits, seen = [], set()
        for art in d:
            for a in art.get("annotations", []):
                t = (a.get("exact") or "").strip()
                if not t or t in seen or t.startswith("10."):
                    continue
                seen.add(t)
                ctx = ((a.get("prefix") or "") + t + (a.get("postfix") or "")).strip()
                hits.append({"accession": t, "context": ctx[:200],
                             "section": a.get("section")})
        out[rid] = {"ok": True, "accessions": hits}
        sys.stderr.write("%s: %d\n" % (rid, len(hits)))
    json.dump(out, open(os.path.join(RAW, "accessions_api.json"), "w"), indent=1)
    print("records with >=1 mined accession: %d/%d" %
          (sum(1 for v in out.values() if v.get("accessions")), len(out)))


if __name__ == "__main__":
    main()
