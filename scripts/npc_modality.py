#!/usr/bin/env python3
"""Distinguish genuine single-cell-genomics modality from lexical look-alikes.

"single cell" in an abstract very often means a comet assay, a cell suspension, a
clonal subline, or an FNA smear pattern -- none of which is single-cell resolution
profiling. This module scores each record on both, so stage-2 reading can be aimed
at real candidates while every rejection still carries a quotable reason.
"""
import re

# genuine single-cell / single-nucleus resolution profiling
TRUE_MOD = {
 "scRNA": r"scrna[\s\-]?seq|\bscrna\b|single[\s\-]?cell\s+rna|single[\s\-]?cell\s+transcriptom|"
          r"single[\s\-]?cell\s+transcriptional|single[\s\-]?cell\s+(gene\s+)?expression\s+profil|"
          r"single[\s\-]?cell\s+sequencing|single[\s\-]?cell\s+atlas|single[\s\-]?cell\s+landscape|"
          r"10[x×]\s+genomics|single[\s\-]?cell\s+resolution|single[\s\-]?cell\s+multi[\s\-]?omic|"
          r"single[\s\-]?cell\s+profiling|droplet[\s\-]?based\s+single|single[\s\-]?cell\s+analy(sis|ses)|"
          r"single[\s\-]?cell\s+data|single[\s\-]?cell\s+dataset|single[\s\-]?cell\s+omics|"
          r"single[\s\-]?cell\s+map(ping)?\b|single[\s\-]?cell\s+dissect|single[\s\-]?cell\s+level\s+"
          r"(transcriptom|expression)|single[\s\-]?cell\s+and\s+(bulk|spatial)|"
          r"single[\s\-]?cell\s+(immune\s+)?(profil|landscap)|sc[\s\-]?rna[\s\-]?seq",
 "snRNA": r"snrna[\s\-]?seq|single[\s\-]?nucle(us|i)\s+rna|single[\s\-]?nucleus\s+(sequenc|transcriptom)",
 "CITE": r"cite[\s\-]?seq|cellular\s+indexing\s+of\s+transcriptomes",
 "CyTOF": r"\bcytof\b|mass\s+cytometry|time[\s\-]of[\s\-]flight\s+cytometry",
 "scTCR": r"sc[\s\-]?tcr|single[\s\-]?cell\s+tcr|tcr[\s\-]?seq|\bvdj\b|v\(d\)j\s+sequenc|"
          r"single[\s\-]?cell\s+t\s?cell\s+receptor",
 "scBCR": r"sc[\s\-]?bcr|single[\s\-]?cell\s+bcr|single[\s\-]?cell\s+b\s?cell\s+receptor|bcr[\s\-]?seq",
 "spatial": r"spatial\s+transcriptom|10[x×]\s+visium|\bvisium\b|geomx|cosmx|merfish|slide[\s\-]?seq|"
            r"spatially\s+resolved\s+transcriptom",
 "scATAC": r"scatac|single[\s\-]?cell\s+atac",
}
TRUE_RE = {k: re.compile(v, re.I) for k, v in TRUE_MOD.items()}

# lexical look-alikes that are NOT single-cell resolution profiling
LOOKALIKE = re.compile(
 r"single[\s\-]?cell\s+gel\s+electrophoresis|single[\s\-]?cell\s+microgel|comet\s+assay|"
 r"single[\s\-]?cell\s+suspension|single[\s\-]?cell[\s\-]?derived|single[\s\-]?cell\s+clone|"
 r"single[\s\-]?cell\s+cloning|single\s+cell\s+pattern|dissociated\s+\(single\s+cell\)|"
 r"single[\s\-]?cell\s+microgel\s+electrophoresis|scge\b", re.I)


def modalities(text):
    return sorted([k for k, rx in TRUE_RE.items() if rx.search(text)])


def lookalike_only(text):
    return (not modalities(text)) and bool(LOOKALIKE.search(text))


def evidence(text, mod):
    """Return the first verbatim matching span's sentence for `mod`."""
    m = TRUE_RE[mod].search(text)
    if not m:
        return None
    lo = text.rfind(".", 0, max(0, m.start() - 1))
    hi = text.find(".", m.end())
    lo = 0 if lo < 0 else lo + 1
    hi = len(text) if hi < 0 else hi + 1
    return re.sub(r"\s+", " ", text[lo:hi]).strip()
