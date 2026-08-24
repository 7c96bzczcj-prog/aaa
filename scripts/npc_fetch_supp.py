#!/usr/bin/env python3
"""Supplementary queries Q7/Q8 -- mandated by spec v1.0 section 5 (anchor repair).

Anchor #15 (Cell 2023;186:4235-4251, PMID 37607536) is closed-access in Europe PMC
(isOpenAccess=N, no PMCID, inEPMC=N), so no full text is indexed and the string
"nasopharyngeal" is absent from every field Europe PMC can search. No query pinned
on TITLE_ABS:"nasopharyngeal" (Q1-Q5) can retrieve it. The repair is a query keyed
on what IS indexed -- pan-cancer single-cell atlas language:

  Q7  pan-cancer single-cell atlas/NK/landscape studies (title/abstract), screened
      afterwards for whether the retrieved text evidences NPC samples.
  Q8  same atlas framing but with an unfielded "nasopharyngeal" term, which Europe
      PMC resolves against full text for OA records -- catches atlases that name NPC
      only in Methods.
"""
import json, os, sys, time, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_fetch import http_get, FAILURES, EPMC  # noqa: E402

QUERIES = {
 "Q7": '(TITLE_ABS:"pan-cancer" OR TITLE_ABS:"pancancer") AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"scRNA-seq") AND (TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell*" OR TITLE_ABS:"atlas" OR TITLE_ABS:"panorama" OR TITLE_ABS:"landscape")',
 "Q8": '(TITLE_ABS:"pan-cancer" OR TITLE_ABS:"pancancer" OR TITLE_ABS:"atlas") AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"scRNA-seq") AND ("nasopharyngeal")',
}


def run_query(qid, q):
    cursor, hits, page, hitcount = "*", [], 0, None
    while True:
        url = (EPMC + "search?query=" + urllib.parse.quote(q, safe="")
               + "&format=json&pageSize=100&resultType=core&cursorMark=" + urllib.parse.quote(cursor, safe=""))
        body = http_get(url)
        if body is None:
            break
        d = json.loads(body)
        if hitcount is None:
            hitcount = d.get("hitCount")
        res = d.get("resultList", {}).get("result", [])
        hits.extend(res)
        nxt = d.get("nextCursorMark")
        page += 1
        sys.stderr.write("%s page %d: +%d (total %d / %s)\n" % (qid, page, len(res), len(hits), hitcount))
        if not nxt or nxt == cursor or not res:
            break
        cursor = nxt
    json.dump({"query": q, "hitCount": hitcount, "retrieved": len(hits), "results": hits},
              open(os.path.join(RAW, "%s_results.json" % qid), "w"))
    return hits, hitcount


if __name__ == "__main__":
    summ = {}
    for qid, q in QUERIES.items():
        h, hc = run_query(qid, q)
        summ[qid] = {"query": q, "hitCount": hc, "retrieved": len(h)}
    json.dump(summ, open(os.path.join(RAW, "query_summary_supp.json"), "w"), indent=2)
    with open(os.path.join(OUT, "fetch_failures.tsv"), "a") as f:
        for u, e in FAILURES:
            f.write("search_supp\t%s\t%s\n" % (u, e))
    print(json.dumps(summ, indent=2)[:2000])
