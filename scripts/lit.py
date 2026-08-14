#!/usr/bin/env python3
"""Literature retrieval helper for the NK lineage-assertion screen.

Exists to make rule G1 ("read the methods, not the abstract") mechanically
checkable: every full text pulled is archived under refs/ so a row's
classification can be traced back to the text it was derived from.

Usage
  lit.py search "<query>" [--n 15] [--src both|pubmed|epmc]
  lit.py fetch <PMID|PMCID|DOI> [--tag claim_id]
  lit.py methods <PMID|PMCID> [--tag claim_id]   # full text, methods-first
"""
import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

REFS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "refs")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UA = {"User-Agent": "nk-lineage-screen/1.0 (research use)"}


def get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"GET failed after {tries}: {url}: {last}")


def search(query, n, src):
    out = {"query": query, "pubmed": [], "epmc": []}
    if src in ("both", "pubmed"):
        q = urllib.parse.quote(query)
        js = json.loads(get(f"{EUTILS}/esearch.fcgi?db=pubmed&term={q}&retmax={n}&retmode=json"))
        ids = js["esearchresult"].get("idlist", [])
        if ids:
            sm = json.loads(get(f"{EUTILS}/esummary.fcgi?db=pubmed&id={','.join(ids)}&retmode=json"))
            for i in ids:
                d = sm["result"].get(i, {})
                out["pubmed"].append({
                    "pmid": i,
                    "title": d.get("title", ""),
                    "journal": d.get("source", ""),
                    "year": (d.get("pubdate", "") or "")[:4],
                    "doi": next((x["value"] for x in d.get("articleids", [])
                                 if x.get("idtype") == "doi"), ""),
                })
    if src in ("both", "epmc"):
        q = urllib.parse.quote(query)
        js = json.loads(get(f"{EPMC}/search?query={q}&format=json&pageSize={n}&resultType=core"))
        for r in js.get("resultList", {}).get("result", []):
            out["epmc"].append({
                "pmid": r.get("pmid", ""),
                "pmcid": r.get("pmcid", ""),
                "title": r.get("title", ""),
                "journal": r.get("journalTitle", ""),
                "year": r.get("pubYear", ""),
                "doi": r.get("doi", ""),
                "isOpenAccess": r.get("isOpenAccess", ""),
                "citedByCount": r.get("citedByCount", ""),
                "abstract": (r.get("abstractText", "") or "")[:1200],
            })
    return out


def resolve(ident):
    """Return (pmid, pmcid, title) for a PMID / PMCID / DOI."""
    ident = ident.strip()
    if ident.upper().startswith("PMC"):
        q = f"PMCID:{ident.upper()}"
    elif re.match(r"^10\.", ident):
        q = f'DOI:"{ident}"'
    else:
        q = f"EXT_ID:{ident} AND SRC:MED"
    js = json.loads(get(f"{EPMC}/search?query={urllib.parse.quote(q)}&format=json&pageSize=1&resultType=core"))
    res = js.get("resultList", {}).get("result", [])
    if not res:
        return ("", "", "")
    r = res[0]
    return (r.get("pmid", ""), r.get("pmcid", ""), r.get("title", ""))


def strip_xml(x):
    x = re.sub(r"<(ref-list|back|table-wrap|fig|xref|table)\b.*?</\1>", " ", x, flags=re.S)
    x = re.sub(r"<title[^>]*>", "\n\n## ", x)
    x = re.sub(r"</title>", "\n", x)
    x = re.sub(r"</(p|sec|abstract|article-title)>", "\n", x)
    x = re.sub(r"<[^>]+>", " ", x)
    x = re.sub(r"&lt;", "<", x)
    x = re.sub(r"&gt;", ">", x)
    x = re.sub(r"&amp;", "&", x)
    x = re.sub(r"[ \t]+", " ", x)
    x = re.sub(r"\n\s*\n\s*\n+", "\n\n", x)
    return x.strip()


def fetch_fulltext(ident):
    """Return (status, text, pmid, pmcid). status in OA_FULLTEXT / ABSTRACT_ONLY / NOT_FOUND."""
    pmid, pmcid, _ = resolve(ident)
    if not pmcid and ident.upper().startswith("PMC"):
        pmcid = ident.upper()
    if pmcid:
        try:
            xml = get(f"{EPMC}/{pmcid}/fullTextXML")
            if len(xml) > 3000 and "<article" in xml:
                return ("OA_FULLTEXT", strip_xml(xml), pmid, pmcid)
        except Exception:  # noqa: BLE001
            pass
    if pmid:
        try:
            xml = get(f"{EUTILS}/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml")
            return ("ABSTRACT_ONLY", strip_xml(xml), pmid, pmcid)
        except Exception:  # noqa: BLE001
            pass
    return ("NOT_FOUND", "", pmid, pmcid)


METHOD_HEADS = re.compile(
    r"^##\s*(.*(method|material|experimental procedure|star.?methods|"
    r"online method|supplementary method|procedure).*)$", re.I | re.M)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search"); s.add_argument("query"); s.add_argument("--n", type=int, default=15)
    s.add_argument("--src", default="both", choices=["both", "pubmed", "epmc"])
    for name in ("fetch", "methods"):
        f = sub.add_parser(name); f.add_argument("ident"); f.add_argument("--tag", default="")
        f.add_argument("--chars", type=int, default=60000)
    a = ap.parse_args()

    if a.cmd == "search":
        print(json.dumps(search(a.query, a.n, a.src), indent=1, ensure_ascii=False))
        return

    status, text, pmid, pmcid = fetch_fulltext(a.ident)
    os.makedirs(REFS, exist_ok=True)
    if text:
        tag = (a.tag + "_") if a.tag else ""
        path = os.path.join(REFS, f"{tag}{pmcid or pmid or 'unknown'}.txt")
        with open(path, "w") as fh:
            fh.write(f"# status={status} pmid={pmid} pmcid={pmcid} ident={a.ident}\n\n{text}")
        sys.stderr.write(f"[archived {path}]\n")
    print(f"STATUS={status} PMID={pmid} PMCID={pmcid}")
    if a.cmd == "methods" and status == "OA_FULLTEXT":
        hits = list(METHOD_HEADS.finditer(text))
        if hits:
            print("METHODS SECTIONS FOUND:", [h.group(1)[:70] for h in hits])
            chunks = []
            for h in hits:
                nxt = text.find("\n## ", h.end())
                # a methods section that runs to the next top-level heading
                chunks.append(text[h.start(): nxt if nxt > 0 else min(len(text), h.start() + 30000)])
            body = "\n\n".join(chunks)
            print(body[: a.chars])
            print(f"\n--- [methods {len(body)} chars; abstract+body available via `fetch`] ---")
            return
        print("NO EXPLICIT METHODS HEADING — printing full text")
    print(text[: a.chars])


if __name__ == "__main__":
    main()
