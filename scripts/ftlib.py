#!/usr/bin/env python3
"""Full-text XML helpers: section extraction + sentence context windows.

Every string this module returns is verbatim text from the retrieved JATS XML,
so anything quoted downstream as source_snippet stays re-verifiable.
"""
import gzip, os, re, xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT = os.path.join(ROOT, "raw", "ft")


def load_xml(pmcid):
    p = os.path.join(FT, pmcid + ".xml.gz")
    if not os.path.exists(p):
        return None
    return gzip.open(p, "rt", encoding="utf-8", errors="replace").read()


def _text(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def parse(xml):
    """Return dict with plain body text and a {section_title: text} map."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        # salvage: strip the DOCTYPE/entities and retry
        cleaned = re.sub(r"<!DOCTYPE[^>]*>", "", xml)
        cleaned = re.sub(r"&(?!(amp|lt|gt|quot|apos);)[a-zA-Z0-9#]+;", " ", cleaned)
        try:
            root = ET.fromstring(cleaned)
        except ET.ParseError:
            return {"plain": re.sub(r"<[^>]+>", " ", xml), "sections": {}, "parsed": False}
    sections = {}
    for sec in root.iter("sec"):
        t = sec.find("title")
        if t is None:
            continue
        title = _text(t)
        if title:
            sections.setdefault(title, []).append(_text(sec))
    sections = {k: " ".join(v) for k, v in sections.items()}
    body = root.find(".//body")
    back = root.find(".//back")
    parts = []
    for node in (body, back):
        if node is not None:
            parts.append(_text(node))
    return {"plain": " ".join(parts), "sections": sections, "parsed": True}


SEC_PATTERNS = {
    "methods": r"^(materials?\s+and\s+methods?|methods?|experimental\s+(procedures?|methods?)|"
               r"star\s*\*?\s*methods?|patients?\s+and\s+methods?|methodology)\b",
    "data_availability": r"(data\s+availability|availability\s+of\s+data|accession\s+(codes?|numbers?)|"
                         r"data\s+and\s+code\s+availability|data\s+access)",
    "results": r"^results?\b",
    "discussion": r"^discussion\b",
}


def pick_sections(sections):
    out = {}
    for key, pat in SEC_PATTERNS.items():
        rx = re.compile(pat, re.I)
        hits = [v for k, v in sections.items() if rx.search(k.strip())]
        if hits:
            out[key] = " ".join(hits)
    return out


SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9(\[])")


def sentences_with(text, pattern, window=0, maxn=40):
    """Verbatim sentences matching `pattern` (optionally +/- `window` neighbours)."""
    rx = re.compile(pattern, re.I)
    sents = SENT_SPLIT.split(text)
    out, seen = [], set()
    for i, s in enumerate(sents):
        if rx.search(s):
            lo, hi = max(0, i - window), min(len(sents), i + window + 1)
            key = (lo, hi)
            if key in seen:
                continue
            seen.add(key)
            out.append(" ".join(sents[lo:hi]).strip())
            if len(out) >= maxn:
                break
    return out
