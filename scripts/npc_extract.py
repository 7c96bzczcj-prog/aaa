#!/usr/bin/env python3
"""Evidence extraction from retrieved text (spec v1.0 section 3).

For each candidate record this pulls verbatim spans out of the API-retrieved
title/abstract and, where available, the Europe PMC fullTextXML -- cohort counts,
tissue sources, sample provenance, treatment, NK depth, EBV-DNA handling and
accessions -- each paired with the field it came from. Nothing here infers a value
that is not present in the retrieved text; absence is reported as absence.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftlib  # noqa: E402
from npc_modality import modalities, TRUE_RE  # noqa: E402

NUM = r"[\d][\d,\. ]*"

PATTERNS = {
 "n_cells": [
   r"([\d][\d,]{2,})\s+(?:high[\s\-]quality\s+)?(?:single\s+)?cells\b",
   r"transcriptomes?\s+of\s+([\d][\d,]{2,})\s+cells",
   r"profil\w+\s+([\d][\d,]{2,})\s+cells",
   r"a\s+total\s+of\s+([\d][\d,]{2,})\s+cells",
   r"([\d][\d,]{2,})\s+cells\s+(?:were|from|across)",
   r"~\s?([\d][\d,]{2,})\s+cells",
 ],
 "n_patients": [
   r"(?:from|of|in)\s+(\d{1,4})\s+(?:NPC\s+)?patients\b",
   r"(\d{1,4})\s+patients\s+(?:with|diagnosed|were)",
   r"n\s*=\s*(\d{1,4})\s+patients",
   r"(\d{1,3})\s+(?:treatment[\s\-]naive\s+)?(?:NPC|nasopharyngeal carcinoma)\s+(?:patients|cases)",
 ],
 "n_samples": [
   r"(\d{1,4})\s+(?:tumou?r\s+)?samples\b",
   r"(\d{1,4})\s+(?:tumou?r\s+)?specimens\b",
   r"(\d{1,4})\s+biopsy\s+(?:and\s+blood\s+)?samples",
 ],
}

TISSUE = {
 "tissue_primary": r"primary\s+(?:NPC|nasopharyngeal|tumou?r|carcinoma|lesion|site)|treatment[\s\-]naive\s+(?:NPC|tumou?r)|primary\s+lesion",
 "tissue_LN_met": r"lymph\s+node\s+metasta|metastatic\s+lymph\s+node|LN\s+metasta|nodal\s+metasta|regional\s+lymph\s+node|cervical\s+lymph\s+node",
 "tissue_distant_met": r"distant\s+metasta|liver\s+metasta|lung\s+metasta|bone\s+metasta|hepatic\s+metasta|"
                       r"pulmonary\s+metasta|metastatic\s+(?:lesion|site|tumou?r|carcinoma)|leptomeningeal",
 "tissue_normal_NP": r"normal\s+nasophary|non[\s\-]?(?:tumou?r|malignant|cancerous)\s+nasophary|nasopharyngeal\s+lymphatic\s+hyperplasia|"
                     r"\bNLH\b|chronic\s+nasopharyngitis|adjacent\s+normal|healthy\s+nasophary|non[\s\-]?malignant\s+nasopharyngeal",
 "tissue_blood": r"peripheral\s+blood|\bPBMC(?:s)?\b|blood\s+samples?|tumou?r[\s\-]blood\s+pair|paired\s+blood|"
                 r"peripheral\s+blood\s+mononuclear",
}

DIST_SITE = r"(liver|hepatic|lung|pulmonary|bone|leptomening|brain|adrenal|spleen|pleural)"

TREATMENT = {
 "naive": r"treatment[\s\-]?naive|therapy[\s\-]?naive|chemo[\s\-]?naive|untreated|before\s+(?:any\s+)?treatment|"
          r"prior\s+to\s+(?:any\s+)?(?:treatment|therapy)",
 "post-CRT": r"post[\s\-]?(?:chemo)?radio|after\s+(?:concurrent\s+)?chemoradio|post[\s\-]?CRT|following\s+radiotherapy|"
             r"recurrent\s+(?:NPC|nasopharyngeal)|after\s+induction\s+chemo|post[\s\-]?GP|after\s+gemcitabine",
 "post-ICI": r"anti[\s\-]?PD[\s\-]?1|anti[\s\-]?PD[\s\-]?L1|immune\s+checkpoint\s+(?:blockade|inhibitor)|"
             r"pembrolizumab|nivolumab|camrelizumab|toripalimab|sintilimab|tislelizumab|avelumab|immunotherap",
}

LONGIT = r"longitudinal|before\s+and\s+after|pre[\s\-]?\s*and\s+post[\s\-]?treatment|paired\s+pre[\s\-]?\s*and|" \
         r"on[\s\-]treatment|serial\s+(?:sampl|blood)|time\s?points?|week\s+\d|baseline\s+and"

M_STAGE = {
 "M1": r"\bM1\b|distant\s+metasta|metastatic\s+disease|de\s+novo\s+metastatic|stage\s+IVb|synchronous\s+metasta",
 "M0": r"\bM0\b|non[\s\-]?metastatic|locoregionally\s+advanced|without\s+distant\s+metasta|M0\s+only",
}

NK_NATURAL = r"natural\s+killer|\bNK\s?-?\s?cells?\b|\bNK\b(?!\s*[/-]\s*T)"
NK_NONKERAT = r"non[\s\-]?kerat|\bNK[\s\-]NPC\b"
NK_SUBCLUSTER = r"NK\s+(?:cell\s+)?(?:sub)?(?:cluster|subset|subpopulation|subtype)s?|" \
                r"(?:sub)?clustering\s+of\s+NK|NK\s+cells?\s+were\s+(?:further\s+)?(?:sub)?clustered|" \
                r"CD56(?:bright|dim)|NK1|NK_c\d|NK\s+cell\s+(?:heterogeneity|diversity|trajector|exhaust|dysfunction)"

EBV_CONT = r"EBV\s+DNA\s+(?:copy|load|titer|titre|level|concentration)s?|copies\/mL|plasma\s+EBV\s+DNA\s+(?:level|load|titer)|" \
           r"correlat\w+\s+with\s+(?:plasma\s+)?EBV\s+DNA|regress\w+.{0,40}EBV\s+DNA|EBV\s+DNA\s+load"
EBV_CAT = r"EBV\s+DNA\s+sero(?:positive|negative|\+|-)|EBV[\s\-]positive|EBV[\s\-]negative|EBV\+|EBV−|EBV\s+DNA\s+(?:positive|negative)"

ACCESSION = r"\b(GSE\d{4,7}|EGAS\d{5,}|EGAD\d{5,}|PRJNA\d{4,}|PRJEB\d{3,}|PRJCA\d{3,}|" \
            r"HRA\d{5,6}|CRA\d{5,6}|OEP\d{6}|SRP\d{5,}|E-MTAB-\d{3,}|phs\d{6}|" \
            r"syn\d{7,}|zenodo\.\d{6,})\b"

PROVENANCE_PUBLIC = r"downloaded\s+from|obtained\s+from\s+the\s+(?:GEO|Gene\s+Expression)|publicly\s+available|" \
                    r"public\s+(?:data|dataset|scRNA)|were\s+retrieved\s+from|re[\s\-]?analy[sz]ed|" \
                    r"from\s+the\s+GEO\s+database|acquired\s+from\s+GEO|data\s+were\s+collected\s+from\s+(?:GEO|the\s+public)"
PROVENANCE_NEW = r"we\s+(?:collected|enrolled|recruited|obtained|prospectively)|were\s+(?:collected|obtained|enrolled|recruited)\s+from\s+" \
                 r"(?:patients|the\s+\w+\s+hospital|\w+\s+Cancer\s+(?:Center|Hospital))|written\s+informed\s+consent|" \
                 r"institutional\s+review\s+board|ethic(?:s|al)\s+(?:committee|approval)|fresh\s+(?:tumou?r\s+)?(?:tissue|sample|biops)"


def first_match(text, patterns, maxlen=200):
    for p in patterns:
        m = re.search(p, text, re.I)
        if m:
            return m.group(1).replace(" ", ""), snippet_around(text, m, maxlen)
    return None, None


def snippet_around(text, m, maxlen=200):
    lo = max(0, m.start() - 90)
    hi = min(len(text), m.end() + 110)
    s = re.sub(r"\s+", " ", text[lo:hi]).strip()
    return s[:maxlen]


def flag(text, pattern, maxlen=200):
    m = re.search(pattern, text, re.I)
    if not m:
        return False, None
    return True, snippet_around(text, m, maxlen)


def extract(rec, ft_plain, ft_sections):
    ab = re.sub(r"<[^>]+>", " ", (rec.get("abstractText") or ""))
    ti = re.sub(r"<[^>]+>", " ", (rec.get("title") or ""))
    abstract_blob = re.sub(r"\s+", " ", ti + ". " + ab)
    sources = [("abstract", abstract_blob)]
    if ft_plain:
        sources.append(("fullText", ft_plain))
    out = {"has_fulltext": bool(ft_plain)}

    # counts: prefer abstract (author-headline numbers), fall back to full text
    for key, pats in PATTERNS.items():
        val = sn = src = None
        for sname, blob in sources:
            v, s = first_match(blob, pats)
            if v:
                val, sn, src = v, s, sname
                break
        out[key] = {"value": val or "NOT_FOUND", "snippet": sn, "source_field": src}

    # tissues
    for key, pat in TISSUE.items():
        got = False
        for sname, blob in sources:
            ok, sn = flag(blob, pat)
            if ok:
                out[key] = {"value": "Y", "snippet": sn, "source_field": sname}
                got = True
                break
        if not got:
            out[key] = {"value": "N", "snippet": None, "source_field": None}
    if out["tissue_distant_met"]["value"] == "Y":
        sites = set()
        for _, blob in sources:
            sites |= {s.lower() for s in re.findall(DIST_SITE, blob, re.I)}
        out["tissue_distant_met"]["sites"] = sorted(sites)

    # treatment status
    tstat = []
    tsnip = {}
    for key, pat in TREATMENT.items():
        for sname, blob in sources:
            ok, sn = flag(blob, pat)
            if ok:
                tstat.append(key)
                tsnip[key] = (sn, sname)
                break
    out["treatment_flags"] = {"flags": tstat, "snippets": tsnip}

    # longitudinal
    for sname, blob in sources:
        ok, sn = flag(blob, LONGIT)
        if ok:
            out["longitudinal"] = {"value": "Y", "snippet": sn, "source_field": sname}
            break
    else:
        out["longitudinal"] = {"value": "N", "snippet": None, "source_field": None}

    # M stage
    mflags = {}
    for key, pat in M_STAGE.items():
        for sname, blob in sources:
            ok, sn = flag(blob, pat)
            if ok:
                mflags[key] = (sn, sname)
                break
    out["M_flags"] = mflags

    # NK: natural killer vs non-keratinizing -- the spec's known trap
    nk = {}
    for sname, blob in sources:
        ok, sn = flag(blob, NK_NATURAL, 240)
        if ok and "natural" not in nk:
            nk["natural"] = (sn, sname)
        ok2, sn2 = flag(blob, NK_NONKERAT, 240)
        if ok2 and "nonkeratinizing" not in nk:
            nk["nonkeratinizing"] = (sn2, sname)
        ok3, sn3 = flag(blob, NK_SUBCLUSTER, 240)
        if ok3 and "subcluster" not in nk:
            nk["subcluster"] = (sn3, sname)
    nk["in_title"] = bool(re.search(r"natural\s+killer|\bNK\s?-?\s?cell", ti, re.I))
    nk["nonkerat_in_title"] = bool(re.search(NK_NONKERAT, ti, re.I))
    nk["nk_mentions_fulltext"] = len(re.findall(NK_NATURAL, ft_plain or "", re.I))
    out["NK"] = nk

    # EBV DNA linkage
    ebv = {}
    for key, pat in (("continuous", EBV_CONT), ("categorical", EBV_CAT)):
        for sname, blob in sources:
            ok, sn = flag(blob, pat, 240)
            if ok:
                ebv[key] = (sn, sname)
                break
    out["EBV"] = ebv

    # accessions
    acc = []
    for sname, blob in sources:
        for m in re.finditer(ACCESSION, blob):
            acc.append(m.group(1))
    out["accessions"] = sorted(set(acc))
    da = (ft_sections or {}).get("data_availability")
    out["data_availability"] = re.sub(r"\s+", " ", da)[:1200] if da else None

    # provenance: new samples vs public reuse
    prov = {}
    for key, pat in (("public", PROVENANCE_PUBLIC), ("new", PROVENANCE_NEW)):
        for sname, blob in sources:
            ok, sn = flag(blob, pat, 240)
            if ok:
                prov[key] = (sn, sname)
                break
    out["provenance"] = prov

    out["modalities"] = modalities(abstract_blob + " " + (ft_plain or "")[:200000])
    out["modalities_abstract"] = modalities(abstract_blob)
    return out


def main():
    ids = json.load(open(os.path.join(RAW, "candidate_ids.json")))
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    res = {}
    for rid in ids:
        r = recs.get(rid)
        if r is None:
            continue
        xml = ftlib.load_xml(r["pmcid"]) if r.get("pmcid") else None
        plain, secs = None, None
        if xml:
            p = ftlib.parse(xml)
            plain, secs = p["plain"], ftlib.pick_sections(p["sections"])
        res[rid] = extract(r, plain, secs)
        res[rid]["_meta"] = {
            "pmid": r.get("pmid"), "doi": r.get("doi"), "pmcid": r.get("pmcid"),
            "year": r.get("pubYear"), "title": re.sub(r"<[^>]+>", "", r.get("title") or ""),
            "source": r.get("source"),
        }
    json.dump(res, open(os.path.join(RAW, "extracted.json"), "w"), indent=1)
    print("extracted %d records (%d with full text)" %
          (len(res), sum(1 for v in res.values() if v["has_fulltext"])))


if __name__ == "__main__":
    main()
