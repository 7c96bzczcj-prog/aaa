#!/usr/bin/env python3
"""Stage-2 review dump: one compact block per candidate, with lexical signals."""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_modality import modalities  # noqa: E402

SIG = {
 "BLOOD": r"peripheral\s+blood|\bPBMC|circulating\s+(immune|lymph|NK|T\s?cell)|blood\s+samples?|paired\s+blood",
 "LNMET": r"lymph\s+node\s+met|LN\s+met|metastatic\s+lymph\s+node|nodal\s+met|regional\s+lymph",
 "DISTMET": r"distant\s+met|liver\s+met|lung\s+met|bone\s+met|metastatic\s+(lesion|site|tumou?r)|M1\b|metastas[ei]s\s+(sample|tissue)",
 "NORMAL": r"normal\s+nasophary|non[\s\-]?(cancerous|tumou?r)|nasopharyngitis|healthy\s+(donor|control)|inflam",
 "NK": r"natural\s+killer|\bNK\s?-?\s?cell",
 "NONKERAT": r"non[\s\-]?kerat",
 "EBV_DNA": r"EBV\s+DNA|Epstein[\s\-]Barr\s+virus\s+DNA|plasma\s+EBV|viral\s+load|copy\s+number.{0,20}EBV",
 "TREAT": r"treatment[\s\-]naive|chemoradio|induction\s+chemo|anti[\s\-]?PD[\s\-]?1|immunotherap|"
          r"post[\s\-]?treatment|neoadjuvant|recurrent|radiotherapy",
 "LONGIT": r"longitudinal|before\s+and\s+after|pre[\s\-]?\s?and\s+post|paired\s+pre|time\s?points?|serial",
 "PUBLIC": r"public(ly)?\s+(available|data)|\bGEO\b|GSE\d+|TCGA|re[\s\-]?analy|integrat(ed|ing)\s+.{0,25}datasets?|"
           r"downloaded|retrieved\s+from",
 "CELLLINE": r"cell\s+lines?|CNE[\s\-]?[12]|HONE|HK1|C666|SUNE|5[\s\-]?8F|6[\s\-]?10B|xenograft|PDX|nude\s+mice",
}
SIG_RE = {k: re.compile(v, re.I) for k, v in SIG.items()}


def main():
    s = json.load(open(os.path.join(RAW, "survivors.json")))
    A = [r for r in s if r["_flags"]["has_nasopharyngeal_term"]]
    A = [r for r in A if modalities((r.get("title") or "") + " " + (r.get("abstractText") or ""))]
    A.sort(key=lambda r: (-(int(r.get("pubYear") or 0)), r.get("pmid") or ""))
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else len(A)
    json.dump([r.get("pmid") or r.get("doi") or r.get("id") for r in A],
              open(os.path.join(RAW, "reviewpool_ids.json"), "w"))
    for i, r in enumerate(A[lo:hi], lo + 1):
        t = (r.get("title") or "")
        ab = re.sub(r"<[^>]+>", "", re.sub(r"\s+", " ", (r.get("abstractText") or ""))).strip()
        blob = t + " " + ab
        sig = [k for k, rx in SIG_RE.items() if rx.search(blob)]
        pl = (r.get("pubTypeList") or {}).get("pubType") or []
        ji = (r.get("journalInfo") or {}).get("journal") or {}
        jn = r.get("journalTitle") or ji.get("medlineAbbreviation") or ji.get("title") or "?"
        print("[%d] %s | pmid=%s | pmcid=%s | %s | %s" % (
            i, r.get("pubYear"), r.get("pmid"), r.get("pmcid") or "-", jn[:32], ",".join(pl)[:52]))
        print("  MOD=%s  SIG=%s" % (",".join(modalities(blob)), ",".join(sig)))
        print("  T: %s" % re.sub(r"<[^>]+>", "", t)[:200])
        print("  A: %s" % (ab[:820] if ab else "<NO ABSTRACT>"))
        print()


if __name__ == "__main__":
    main()
