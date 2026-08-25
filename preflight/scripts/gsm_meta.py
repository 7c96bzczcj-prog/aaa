#!/usr/bin/env python3
"""Fetch GSM-level metadata for a list of GSE and classify compartment/age."""
import json, os, re, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

BM_RE = re.compile(r"\b(bone[\s_-]?marrow|bmmc|\bBM\b|marrow|aspirate|iliac)\b", re.I)
PB_RE = re.compile(r"\b(peripheral[\s_-]?blood|pbmc|\bPB\b|whole[\s_-]?blood|apheresis|leukapheresis|buffy)\b", re.I)
CB_RE = re.compile(r"\b(cord[\s_-]?blood|umbilical|\bCB\b)\b", re.I)
AGE_RE = re.compile(r"age|yrs|years|\byo\b", re.I)


def get(url, tries=4):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if a == tries - 1:
                return ""
            time.sleep(2 ** a)
    return ""


def parse(gse):
    txt = get(f"{ACC}?acc={gse}&targ=gsm&form=text&view=brief")
    samples, cur = [], None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=")[-1].strip(), "char": [], "suppl": []}
        elif cur is None:
            continue
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[-1].strip()
        elif line.startswith("!Sample_source_name"):
            cur["source"] = line.split("=", 1)[-1].strip()
        elif line.startswith("!Sample_characteristics"):
            cur["char"].append(line.split("=", 1)[-1].strip())
        elif line.startswith("!Sample_supplementary_file"):
            cur["suppl"].append(line.split("=", 1)[-1].strip().rsplit("/", 1)[-1])
        elif line.startswith("!Sample_library_strategy"):
            cur["strategy"] = line.split("=", 1)[-1].strip()
    if cur:
        samples.append(cur)
    for s in samples:
        blob = " ".join([s.get("title", ""), s.get("source", "")] + s["char"])
        s["blob"] = blob
        s["is_cb"] = bool(CB_RE.search(blob))
        s["is_bm"] = bool(BM_RE.search(blob)) and not s["is_cb"]
        s["is_pb"] = bool(PB_RE.search(blob)) and not s["is_cb"]
        s["has_age"] = bool(AGE_RE.search(" ".join(s["char"])))
    return {"accession": gse, "n_gsm": len(samples),
            "n_bm": sum(s["is_bm"] for s in samples),
            "n_pb": sum(s["is_pb"] for s in samples),
            "n_cb": sum(s["is_cb"] for s in samples),
            "n_age": sum(s["has_age"] for s in samples),
            "samples": samples}


if __name__ == "__main__":
    gses = sys.argv[1].split(",") if "," in sys.argv[1] else json.load(open(sys.argv[1]))
    out = sys.argv[2]
    with ThreadPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(parse, gses))
    json.dump(res, open(out, "w"), indent=1)
    for r in sorted(res, key=lambda x: -(min(x["n_bm"], x["n_pb"]))):
        print(f"{r['accession']:<12} gsm={r['n_gsm']:<4} BM={r['n_bm']:<4} PB={r['n_pb']:<4} CB={r['n_cb']:<3} age_fields={r['n_age']}")
