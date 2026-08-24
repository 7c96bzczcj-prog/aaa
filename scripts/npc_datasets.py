#!/usr/bin/env python3
"""Accession-level dataset inventory -- raw material for Q-F (spec v1.0 section 6).

Collects every accession mentioned in the retrieved text of the included and
secondary records, plus the GEO DataSets hits from Q6, then queries NCBI
E-utilities for each GSE's real sample composition so `includes_blood` is decided
from repository metadata rather than from the paper's prose.
"""
import json, os, re, sys, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = os.path.join(ROOT, "raw"), os.path.join(ROOT, "out")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_fetch import http_get  # noqa: E402

ACC_RE = re.compile(r"\b(GSE\d{4,7}|HRA\d{6}|CRA\d{6}|OEP\d{6}|OMIX\d{6}|CNP\d{7}|PRJNA\d{5,}|EGAS\d{11})\b")
# NB: no \b before PBMC -- GEO sample titles look like "NPC_SC_1810_PBMC_cDNA",
# where the preceding underscore is a word character and \b would never match.
BLOOD_RE = re.compile(r"PBMC|peripheral[\s_-]*blood|blood|leuko|buffy", re.I)
NPC_RE = re.compile(r"nasophary|\bNPC\b", re.I)
SC_RE = re.compile(r"single[\s\-]?cell|single[\s\-]?nucle|scRNA|snRNA|10x|CyTOF|mass cytometry|"
                   r"CITE[\s\-]?seq|spatial", re.I)


def geo_summary(gses):
    """Resolve GSE accessions to GEO DataSets records via E-utilities."""
    out = {}
    for i in range(0, len(gses), 20):
        chunk = gses[i:i + 20]
        term = " OR ".join("%s[ACCN]" % g for g in chunk)
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
               + urllib.parse.quote(term, safe="") + "&retmax=200&retmode=json")
        b = http_get(url, retries=2, timeout=60)
        if not b:
            continue
        ids = json.loads(b)["esearchresult"].get("idlist", [])
        if not ids:
            continue
        url2 = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id="
                + ",".join(ids) + "&retmode=json")
        b2 = http_get(url2, retries=2, timeout=90)
        if not b2:
            continue
        res = json.loads(b2).get("result", {})
        for uid in res.get("uids", []):
            d = res[uid]
            acc = d.get("accession")
            if acc and acc.startswith("GSE"):
                out[acc] = d
    return out


def main():
    finals = json.load(open(os.path.join(RAW, "final_sets.json")))
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from npc_manifest_data import M
    from npc_adjudication import VERDICT

    # where each accession was seen
    seen = {}
    for rid in finals["included"]:
        for a in ACC_RE.findall(M[rid]["acc"] + " " + M[rid].get("note", "")):
            seen.setdefault(a, {"included_in": [], "secondary_in": []})["included_in"].append(
                recs[rid].get("pmid") or rid)
    for rid in finals["secondary"]:
        for a in ACC_RE.findall(VERDICT[rid][2]):
            seen.setdefault(a, {"included_in": [], "secondary_in": []})["secondary_in"].append(
                recs[rid].get("pmid") or rid)
    # full-text mined accessions from the included records
    ext = json.load(open(os.path.join(RAW, "extracted.json")))
    for rid in finals["included"]:
        for a in ext.get(rid, {}).get("accessions", []):
            for x in ACC_RE.findall(a):
                seen.setdefault(x, {"included_in": [], "secondary_in": []})

    gses = sorted(a for a in seen if a.startswith("GSE"))
    sys.stderr.write("resolving %d GSE accessions against GEO...\n" % len(gses))
    geo = geo_summary(gses)

    # Q6 GEO DataSets sweep
    q6 = json.load(open(os.path.join(RAW, "Q6_gds.json")))
    for d in q6:
        acc = d.get("accession")
        if acc and acc.startswith("GSE"):
            geo.setdefault(acc, d)
            seen.setdefault(acc, {"included_in": [], "secondary_in": [], "from_Q6": True})
            seen[acc]["from_Q6"] = True

    rows = []
    for acc in sorted(seen):
        d = geo.get(acc)
        if d:
            blob = " ".join(str(d.get(k) or "") for k in
                            ("title", "summary", "seriestitle", "gdstype", "taxon"))
            samples = d.get("samples") or []
            sblob = " ".join((s.get("title") or "") for s in samples) if isinstance(samples, list) else ""
            n = d.get("n_samples") or (len(samples) if isinstance(samples, list) else "NOT_FOUND")
            # Decide includes_blood from SAMPLE titles (what is actually deposited),
            # never from the series summary -- GEO summaries repeat the paper abstract,
            # so a study whose blood arm sits in a different (often controlled-access)
            # repository would otherwise be scored as blood-bearing. GSE206245 is exactly
            # that case: summary says "77 biopsy and blood samples", deposited samples are
            # spatial only.
            # The series TITLE does describe the deposited material, so it counts;
            # the summary does not.
            tblob = str(d.get("title") or "") + " " + str(d.get("seriestitle") or "")
            if BLOOD_RE.search(sblob) or BLOOD_RE.search(tblob):
                blood = "Y"
            elif sblob:
                blood = "N"
            elif BLOOD_RE.search(blob):
                blood = "UNKNOWN_summary_only"
            else:
                blood = "UNKNOWN"
            rows.append({
                "accession": acc, "repository": "GEO",
                "n_samples": str(n),
                "is_NPC": "Y" if NPC_RE.search(blob) else "N",
                "is_single_cell": "Y" if SC_RE.search(blob) else "N",
                "includes_blood": blood,
                "taxon": d.get("taxon") or "NOT_FOUND",
                "title": re.sub(r"\s+", " ", (d.get("title") or ""))[:180],
                "sample_titles": re.sub(r"\s+", " ", sblob)[:400] or "NOT_FOUND",
                "cited_by_included": ",".join(seen[acc]["included_in"]) or "-",
                "cited_by_secondary": ",".join(seen[acc]["secondary_in"]) or "-",
            })
        else:
            repo = ("GSA-Human" if acc.startswith(("HRA", "CRA")) else
                    "NODE" if acc.startswith("OEP") else
                    "OMIX" if acc.startswith("OMIX") else
                    "CNGB" if acc.startswith("CNP") else
                    "SRA" if acc.startswith("PRJNA") else
                    "EGA" if acc.startswith("EGAS") else "other")
            rows.append({
                "accession": acc, "repository": repo, "n_samples": "NOT_FOUND",
                "is_NPC": "NOT_FOUND", "is_single_cell": "NOT_FOUND",
                "includes_blood": "UNKNOWN", "taxon": "NOT_FOUND",
                "title": "NOT_FOUND (not indexed in NCBI GEO; controlled-access or non-NCBI repository)",
                "sample_titles": "NOT_FOUND",
                "cited_by_included": ",".join(seen[acc]["included_in"]) or "-",
                "cited_by_secondary": ",".join(seen[acc]["secondary_in"]) or "-",
            })

    hdr = ["accession", "repository", "n_samples", "is_NPC", "is_single_cell",
           "includes_blood", "taxon", "title", "sample_titles",
           "cited_by_included", "cited_by_secondary"]
    with open(os.path.join(OUT, "npc_sc_datasets.tsv"), "w") as f:
        f.write("# Accession-level inventory (raw material for Q-F)\n")
        f.write("# run date: 2026-08-24 | GEO fields from NCBI E-utilities esummary db=gds\n")
        f.write("# includes_blood is decided from repository metadata (series/sample titles), not from paper prose\n")
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in hdr) + "\n")
    json.dump(rows, open(os.path.join(RAW, "datasets_rows.json"), "w"), indent=1)
    print("accessions inventoried:", len(rows))
    print("  GEO-resolved:", sum(1 for r in rows if r["repository"] == "GEO"))
    print("  NPC + single-cell + blood:",
          sum(1 for r in rows if r["is_NPC"] == "Y" and r["is_single_cell"] == "Y"
              and r["includes_blood"] == "Y"))


if __name__ == "__main__":
    main()
