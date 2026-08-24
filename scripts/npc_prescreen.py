#!/usr/bin/env python3
"""Stage-1 mechanical pre-screen (spec v1.0 section 2).

Purely lexical over API-retrieved title+abstract. Its ONLY job is to split the
2438 deduped records into (a) records that cannot possibly qualify -- no
single-cell modality term anywhere in the retrieved text -- and (b) a survivor
pool that gets read individually. Nothing is included by this script; inclusion
is decided at stage 2 from the retrieved text.
"""
import json, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "raw")

# --- single-cell resolution modality evidence -------------------------------
MODALITY = [
 r"single[\s\-]?cell\s+rna", r"scrna[\s\-]?seq", r"sc[\s\-]?rna", r"single[\s\-]?cell\s+transcriptom",
 r"snrna[\s\-]?seq", r"single[\s\-]?nucle(us|ar|i)", r"single[\s\-]?cell\s+sequenc",
 r"single[\s\-]?cell\s+analys", r"single[\s\-]?cell\s+profil", r"single[\s\-]?cell\s+resolution",
 r"single[\s\-]?cell\s+atlas", r"single[\s\-]?cell\s+landscape", r"single[\s\-]?cell\s+multi",
 r"cite[\s\-]?seq", r"cytof", r"mass\s+cytometry", r"sc[\s\-]?tcr", r"sc[\s\-]?bcr",
 r"single[\s\-]?cell\s+tcr", r"single[\s\-]?cell\s+bcr", r"tcr[\s\-]?seq", r"vdj\s+sequenc",
 r"spatial\s+transcriptom", r"10x\s+genomics", r"10×\s+genomics", r"single[\s\-]?cell\s+immune",
 r"single[\s\-]?cell\b", r"single\s+cells\b", r"scatac", r"single[\s\-]?cell\s+atac",
 r"smart[\s\-]?seq", r"drop[\s\-]?seq", r"seq[\s\-]?well", r"microwell[\s\-]?seq",
 r"single[\s\-]?cell\s+level", r"single[\s\-]?cell\s+map", r"single[\s\-]?cell\s+dissect",
]
MOD_RE = re.compile("|".join(MODALITY), re.I)

NPC = re.compile(r"nasophary", re.I)

# --- publication types that are never primary research ----------------------
REVIEW_PT = re.compile(r"\b(review|meta[\s\-]?analysis|editorial|comment|letter|"
                       r"systematic review|corrigend|erratum|retract|preface|news)\b", re.I)
REVIEW_TITLE = re.compile(r"\b(a review|review of|systematic review|meta[\s\-]?analysis|"
                          r"narrative review|scoping review|mini[\s\-]?review|"
                          r"bibliometric|current (advances|progress|status)|"
                          r"recent (advances|progress)|perspectives?\b|overview of)\b", re.I)

# --- the NK trap: non-keratinizing vs natural killer -------------------------
NK_NONKERAT = re.compile(r"non[\s\-]?kerat", re.I)
NK_NATURAL = re.compile(r"natural\s+killer|\bNK\s+cell|\bNK\s?-?cells", re.I)
NKT_LYMPHOMA = re.compile(r"NK\s*/\s*T[\s\-]?cell\s+lymphoma|NKTL|natural\s+killer/?T[\s\-]?cell\s+lymphoma|"
                          r"extranodal\s+NK", re.I)


def text_of(r):
    return " ".join([(r.get("title") or ""), (r.get("abstractText") or "")])


def main():
    recs = json.load(open(os.path.join(RAW, "deduped.json")))
    survivors, dropped = [], []
    for r in recs:
        t = text_of(r)
        pt = (r.get("pubType") or "")
        has_mod = bool(MOD_RE.search(t))
        has_npc = bool(NPC.search(t))
        no_abstract = not (r.get("abstractText") or "").strip()
        r["_flags"] = {
            "has_modality_term": has_mod,
            "has_nasopharyngeal_term": has_npc,
            "no_abstract": no_abstract,
            "pubtype_review": bool(REVIEW_PT.search(pt)),
            "title_review": bool(REVIEW_TITLE.search(r.get("title") or "")),
            "nkt_lymphoma": bool(NKT_LYMPHOMA.search(t)),
            "nk_nonkeratinizing_term": bool(NK_NONKERAT.search(t)),
            "nk_natural_killer_term": bool(NK_NATURAL.search(t)),
            "queries": r.get("_queries"),
        }
        if not has_mod and not no_abstract:
            r["_drop_reason"] = "no_single_cell_modality_term_in_title_or_abstract"
            dropped.append(r)
        else:
            survivors.append(r)

    json.dump(survivors, open(os.path.join(RAW, "survivors.json"), "w"))
    json.dump(dropped, open(os.path.join(RAW, "stage1_dropped.json"), "w"))
    print("deduped   :", len(recs))
    print("survivors :", len(survivors))
    print("dropped   :", len(dropped))
    print("closure   :", len(survivors) + len(dropped) == len(recs))
    from collections import Counter
    c = Counter()
    for r in survivors:
        f = r["_flags"]
        c["npc_term"] += f["has_nasopharyngeal_term"]
        c["no_abstract"] += f["no_abstract"]
        c["pubtype_review"] += f["pubtype_review"]
        c["nkt_lymphoma"] += f["nkt_lymphoma"]
        c["nonkeratinizing"] += f["nk_nonkeratinizing_term"]
    print(dict(c))


if __name__ == "__main__":
    main()
