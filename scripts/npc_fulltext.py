#!/usr/bin/env python3
"""Full-text retrieval for survivors with a PMCID (spec v1.0 section 1).

GET https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML
Stores the XML gzipped under raw/ft/ so every extracted snippet stays re-verifiable,
and builds a compact index recording whether the article's full text evidences NPC.
Failures land in out/fetch_failures.tsv.
"""
import gzip, json, os, sys, re, time, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = os.path.join(ROOT, "raw"), os.path.join(ROOT, "out")
FT = os.path.join(RAW, "ft")
os.makedirs(FT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_fetch import http_get, FAILURES, EPMC  # noqa: E402

NPC_RE = re.compile(r"nasophary", re.I)
NPC_ABBR = re.compile(r"\bNPC\b")


def main():
    survivors = json.load(open(os.path.join(RAW, "survivors.json")))
    targets = [r for r in survivors if r.get("pmcid")]
    sys.stderr.write("full text targets: %d of %d survivors\n" % (len(targets), len(survivors)))
    index = {}
    for i, r in enumerate(targets, 1):
        pmcid = r["pmcid"]
        path = os.path.join(FT, pmcid + ".xml.gz")
        if os.path.exists(path):
            xml = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
        else:
            body = http_get(EPMC + urllib.parse.quote(pmcid) + "/fullTextXML")
            if body is None:
                index[pmcid] = {"pmid": r.get("pmid"), "ok": False}
                continue
            xml = body.decode("utf-8", errors="replace")
            if "<" not in xml[:200] or "not found" in xml[:400].lower():
                index[pmcid] = {"pmid": r.get("pmid"), "ok": False, "note": "no_xml_body"}
                FAILURES.append((pmcid, "no_xml_body"))
                continue
            with gzip.open(path, "wt", encoding="utf-8") as f:
                f.write(xml)
        index[pmcid] = {
            "pmid": r.get("pmid"), "ok": True, "chars": len(xml),
            "npc_full": bool(NPC_RE.search(xml)),
            "npc_abbr": bool(NPC_ABBR.search(xml)),
            "npc_count": len(NPC_RE.findall(xml)),
        }
        if i % 50 == 0:
            sys.stderr.write("  %d/%d\n" % (i, len(targets)))
    json.dump(index, open(os.path.join(RAW, "fulltext_index.json"), "w"), indent=1)
    with open(os.path.join(OUT, "fetch_failures.tsv"), "a") as f:
        for u, e in FAILURES:
            f.write("fulltext\t%s\t%s\n" % (u, e))
    ok = sum(1 for v in index.values() if v.get("ok"))
    print("fetched ok: %d / %d ; failed: %d" % (ok, len(targets), len(targets) - ok))
    print("full text mentioning 'nasophary': %d" % sum(1 for v in index.values() if v.get("npc_full")))


if __name__ == "__main__":
    main()
