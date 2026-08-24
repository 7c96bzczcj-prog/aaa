#!/usr/bin/env python3
"""NPC single-cell occupancy audit -- retrieval layer (spec v1.0 section 1).

Runs Q1..Q5 against Europe PMC REST (cursorMark pagination, resultType=core so
abstractText is present), Q6 against NCBI E-utilities (db=gds).
Rate limit <=3 req/s, 3 retries, failures logged to out/fetch_failures.tsv.
Nothing is inferred here -- raw API payloads are stored verbatim under raw/.
"""
import json, os, sys, time, urllib.parse, urllib.request, urllib.error, hashlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
OUT = os.path.join(ROOT, "out")
os.makedirs(RAW, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
UA = "npc-sc-audit/1.0 (research audit; mailto:wubide666006@gmail.com)"

FAILURES = []
_last_req = [0.0]
MIN_INTERVAL = 0.34  # <=3 req/s


def http_get(url, retries=3, timeout=90):
    """GET with rate limiting + retries. Returns bytes or None (logged)."""
    for attempt in range(1, retries + 1):
        wait = MIN_INTERVAL - (time.time() - _last_req[0])
        if wait > 0:
            time.sleep(wait)
        _last_req[0] = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            sys.stderr.write("  attempt %d/%d failed: %s | %s\n" % (attempt, retries, type(e).__name__, str(e)[:160]))
            if attempt == retries:
                FAILURES.append((url, "%s: %s" % (type(e).__name__, str(e)[:300])))
                return None
            time.sleep(2 ** attempt)
    return None


QUERIES = {
 "Q1": '(TITLE_ABS:"nasopharyngeal") AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"scRNA-seq" OR TITLE_ABS:"snRNA-seq" OR TITLE_ABS:"single-nucleus")',
 "Q2": '(TITLE_ABS:"nasopharyngeal") AND (TITLE_ABS:"CyTOF" OR TITLE_ABS:"mass cytometry" OR TITLE_ABS:"CITE-seq" OR TITLE_ABS:"spatial transcriptom*")',
 "Q3": '(TITLE_ABS:"nasopharyngeal") AND (TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell*")',
 "Q4": '(TITLE_ABS:"nasopharyngeal") AND (TITLE_ABS:"peripheral blood" OR TITLE_ABS:"PBMC" OR TITLE_ABS:"circulating")',
 "Q5": '(TITLE_ABS:"nasopharyngeal") AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"scRNA-seq" OR TITLE_ABS:"snRNA-seq" OR TITLE_ABS:"single-nucleus") AND (SRC:"PPR")',
}


def run_query(qid, q):
    cursor = "*"
    hits = []
    page = 0
    hitcount = None
    while True:
        url = (EPMC + "search?query=" + urllib.parse.quote(q, safe="")
               + "&format=json&pageSize=100&resultType=core&cursorMark=" + urllib.parse.quote(cursor, safe=""))
        body = http_get(url)
        if body is None:
            sys.stderr.write("%s: page %d permanently failed\n" % (qid, page))
            break
        d = json.loads(body)
        if hitcount is None:
            hitcount = d.get("hitCount")
        res = d.get("resultList", {}).get("result", [])
        hits.extend(res)
        nxt = d.get("nextCursorMark")
        page += 1
        sys.stderr.write("%s page %d: +%d (total %d / hitCount %s)\n" % (qid, page, len(res), len(hits), hitcount))
        if not nxt or nxt == cursor or not res:
            break
        cursor = nxt
    with open(os.path.join(RAW, "%s_results.json" % qid), "w") as f:
        json.dump({"query": q, "hitCount": hitcount, "retrieved": len(hits), "results": hits}, f)
    return hits, hitcount


def main():
    summary = {}
    for qid, q in QUERIES.items():
        sys.stderr.write("=== %s ===\n" % qid)
        hits, hc = run_query(qid, q)
        summary[qid] = {"query": q, "hitCount": hc, "retrieved": len(hits)}

    # --- Q6: GEO DataSets via E-utilities ---
    term = "nasopharyngeal carcinoma AND single cell"
    esearch = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
               + urllib.parse.quote(term, safe="") + "&retmax=500&retmode=json&usehistory=y")
    body = http_get(esearch)
    ids = []
    if body:
        j = json.loads(body)
        ids = j["esearchresult"].get("idlist", [])
        with open(os.path.join(RAW, "Q6_esearch.json"), "w") as f:
            json.dump(j, f)
    sys.stderr.write("Q6 esearch ids: %d\n" % len(ids))
    docs = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        esum = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id="
                + ",".join(chunk) + "&retmode=json")
        b = http_get(esum)
        if not b:
            continue
        j = json.loads(b)
        r = j.get("result", {})
        for uid in r.get("uids", []):
            docs.append(r[uid])
    with open(os.path.join(RAW, "Q6_gds.json"), "w") as f:
        json.dump(docs, f)
    summary["Q6"] = {"term": term, "esearch_ids": len(ids), "esummary_docs": len(docs)}

    with open(os.path.join(RAW, "query_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    with open(os.path.join(OUT, "fetch_failures.tsv"), "w") as f:
        f.write("stage\turl\terror\n")
        for u, e in FAILURES:
            f.write("search\t%s\t%s\n" % (u, e))
    print(json.dumps(summary, indent=2)[:4000])


if __name__ == "__main__":
    main()
