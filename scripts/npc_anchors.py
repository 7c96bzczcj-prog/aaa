#!/usr/bin/env python3
"""Anchor recall self-check (spec v1.0 section 5)."""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")

ANCHORS = [
 (1,  "Nat Commun 2021;12:741 -- 10 NPC tumour-blood pairs, 176,447 cells, scRNA+TCR",
      {"jvp": ("nat commun", "12", "741"), "year": "2021"}),
 (2,  "Cell Res 2020 doi 10.1038/s41422-020-0374-x -- GSE150430",
      {"doi": "10.1038/s41422-020-0374-x"}),
 (3,  "PMID 32901110 -- Cell Res 2020, ~104,000 cells, 19 EBV+ NPC", {"pmid": "32901110"}),
 (4,  "PMC7943808 -- Nat Commun 2021, 66,627 cells, 14 pts, stromal dynamics", {"pmcid": "PMC7943808"}),
 (5,  "PMID 32061950 -- 2020 NPC tumour + infiltrating immune landscape", {"pmid": "32061950"}),
 (6,  "PMC8794254 -- primary vs recurrent NPC, ~60,000 cells", {"pmcid": "PMC8794254"}),
 (7,  "doi 10.1186/s12967-023-04112-8 -- J Transl Med 2023 NK exhaustion, n=3",
      {"doi": "10.1186/s12967-023-04112-8"}),
 (8,  "PMID 37349991 -- Clin Transl Med 2023 CD8+ NK cluster", {"pmid": "37349991"}),
 (9,  "doi 10.1111/jcmm.70137 -- J Cell Mol Med 2024, primary vs LN met, 47,618 cells",
      {"doi": "10.1111/jcmm.70137"}),
 (10, "doi 10.1038/s41467-023-35995-2 -- Nat Commun 2023 metastatic evolution",
      {"doi": "10.1038/s41467-023-35995-2"}),
 (11, "doi 10.1038/s41467-024-52153-4 -- Nat Commun 2024 TLS, 343,829 cells",
      {"doi": "10.1038/s41467-024-52153-4"}),
 (12, "doi 10.3389/fimmu.2023.1124066 -- Front Immunol 2023 EBV DNA sero+/-",
      {"doi": "10.3389/fimmu.2023.1124066"}),
 (13, "doi 10.1038/s41392-024-01988-w -- STTT 2024 CONTINUUM PBMC CyTOF",
      {"doi": "10.1038/s41392-024-01988-w"}),
 (14, "PMID 39415331 -- paired pre/post immunochemotherapy, 11 pts, 87,191 cells",
      {"pmid": "39415331"}),
 (15, "Cell 2023;186:4235-4251 -- Tang F et al pan-cancer NK atlas",
      {"jvp": ("cell", "186", "4235-4251"), "year": "2023"}),
]


def match(rec, spec):
    if "pmid" in spec and (rec.get("pmid") or "") == spec["pmid"]:
        return True
    if "doi" in spec and (rec.get("doi") or "").lower() == spec["doi"].lower():
        return True
    if "pmcid" in spec and (rec.get("pmcid") or "").upper() == spec["pmcid"].upper():
        return True
    if "jvp" in spec:
        # Europe PMC leaves the flat journalTitle/journalVolume null on many core
        # records; the populated copy lives under journalInfo. Check both.
        j, v, p = spec["jvp"]
        ji = rec.get("journalInfo") or {}
        jj = (ji.get("journal") or {})
        names = {(rec.get("journalTitle") or "").lower().rstrip("."),
                 (jj.get("medlineAbbreviation") or "").lower().rstrip("."),
                 (jj.get("isoabbreviation") or "").lower().rstrip("."),
                 (jj.get("title") or "").lower().rstrip(".")}
        vols = {str(rec.get("journalVolume") or ""), str(ji.get("volume") or "")}
        rp = (rec.get("pageInfo") or "")
        if (j.replace(".", "") in {n.replace(".", "") for n in names}
                and v in vols and (rp == p or rp.startswith(p.split("-")[0]))):
            return True
    return False


def main():
    recs = json.load(open(os.path.join(RAW, "deduped.json")))
    rows = []
    for num, label, spec in ANCHORS:
        hit = [r for r in recs if match(r, spec)]
        rows.append({
            "anchor": num, "label": label, "spec": {k: (list(v) if isinstance(v, tuple) else v) for k, v in spec.items()},
            "hit": bool(hit),
            "pmid": hit[0].get("pmid") if hit else None,
            "doi": hit[0].get("doi") if hit else None,
            "pmcid": hit[0].get("pmcid") if hit else None,
            "title": (hit[0].get("title") or "")[:150] if hit else None,
            "queries": hit[0].get("_queries") if hit else None,
        })
    json.dump(rows, open(os.path.join(RAW, "anchor_check.json"), "w"), indent=2)
    miss = 0
    for r in rows:
        st = "HIT " if r["hit"] else "MISS"
        if not r["hit"]:
            miss += 1
        print("%s  #%-2d %s" % (st, r["anchor"], r["label"]))
        if r["hit"]:
            print("       -> pmid=%s doi=%s q=%s" % (r["pmid"], r["doi"], r["queries"]))
    print("\nHIT %d / MISS %d of %d" % (len(rows) - miss, miss, len(rows)))


if __name__ == "__main__":
    main()
