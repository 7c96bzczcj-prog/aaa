#!/usr/bin/env python3
"""Per-record evidence sheet for manifest field curation (spec v1.0 section 3).

Prints, for one included record, the verbatim sentences bearing on each schema
field, so each manifest value can be written with a quotable <=200-char snippet
and the field it came from. Regex-guessed numbers are deliberately NOT printed as
answers -- only the sentences that contain them.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftlib  # noqa: E402

QUERIES = [
 ("COHORT", r"\b\d[\d,]{2,}\s+(?:high[\s\-]quality\s+)?(?:single\s+)?cells\b|\b\d{1,4}\s+(?:NPC\s+)?patients\b|"
            r"\b\d{1,4}\s+samples\b|\bn\s*=\s*\d{1,4}\b|\b\d{1,3}\s+(?:tumou?rs?|biopsies|specimens)\b"),
 ("TISSUE", r"primary\s+(?:tumou?r|NPC|lesion|site)|lymph\s+node|metasta|peripheral\s+blood|\bPBMC|"
            r"normal\s+nasophary|non[\s\-]?(?:tumou?r|malignant)|hyperplasia|\bNLH\b|nasopharyngitis|"
            r"liver|lung|bone|leptomening|cerebrospinal|\bCSF\b"),
 ("TREATMENT", r"treatment[\s\-]?naive|untreated|chemoradio|induction\s+chemo|anti[\s\-]?PD[\s\-]?[1L]|"
               r"immunotherap|recurrent|radiotherapy|gemcitabine|cisplatin|neoadjuvant|post[\s\-]?treatment|"
               r"camrelizumab|toripalimab|sintilimab|avelumab|nivolumab|pembrolizumab"),
 ("LONGITUDINAL", r"longitudinal|before\s+and\s+after|pre[\s\-]?\s*and\s+post|on[\s\-]treatment|"
                  r"time\s?points?|baseline|week\s+\d|serial|paired\s+pre|matched.{0,20}post"),
 ("MSTAGE", r"\bM[01]\b|distant\s+metasta|de\s+novo\s+metastatic|stage\s+I{1,3}V?|locoregionally\s+advanced|"
            r"metastatic\s+NPC|TNM|T1N0M0"),
 ("NK", r"natural\s+killer|\bNK\s?-?\s?cell|\bNK\d|CD56|non[\s\-]?kerat|NK[\s\-]NPC"),
 ("EBV", r"EBV\s+DNA|plasma\s+EBV|viral\s+load|copies\s*/\s*mL|sero(?:positive|negative)|EBV[\s\-]?(?:positive|negative)|EBER"),
 ("ACCESSION", r"GSE\d{4,7}|HRA\d{5,6}|CRA\d{5,6}|OEP\d{6}|OMIX\d{6}|PRJNA\d{4,}|EGAS\d{5,}|CNP\d{7}|"
               r"zenodo|SRP\d{5,}|deposited|accession"),
]


def main():
    rid = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    r = recs[rid]
    ti = re.sub(r"<[^>]+>", "", r.get("title") or "")
    ab = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", r.get("abstractText") or "")).strip()
    ji = (r.get("journalInfo") or {}).get("journal") or {}
    print("#" * 100)
    print("id=%s pmid=%s doi=%s pmcid=%s year=%s journal=%s" % (
        rid, r.get("pmid"), r.get("doi"), r.get("pmcid"), r.get("pubYear"),
        r.get("journalTitle") or ji.get("medlineAbbreviation") or ji.get("title")))
    print("AUTHORS: %s" % (r.get("authorString") or "")[:200])
    print("AFFIL  : %s" % (r.get("affiliation") or "NOT_FOUND")[:300])
    print("TITLE  : %s" % ti)
    print("ABSTRACT: %s" % ab[:1500])
    xml = ftlib.load_xml(r["pmcid"]) if r.get("pmcid") else None
    if not xml:
        print("\n[NO FULL TEXT -> abstract_only]")
        return
    p = ftlib.parse(xml)
    plain = p["plain"]
    for label, pat in QUERIES:
        print("\n--- %s ---" % label)
        got = ftlib.sentences_with(plain, pat, maxn=200)
        seen, k = set(), 0
        for s in got:
            key = s[:60]
            if key in seen:
                continue
            seen.add(key)
            print("  >", re.sub(r"\s+", " ", s)[:300])
            k += 1
            if k >= n:
                break
        if not k:
            print("  (none)")


if __name__ == "__main__":
    main()
