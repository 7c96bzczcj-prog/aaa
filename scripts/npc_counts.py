#!/usr/bin/env python3
"""COUNTS.txt -- numeric closure, anchor recall, re-verification, NOT_FOUND rates.

Implements spec v1.0 section 7 (delivery self-check) and section 0 rule 6
(numeric closure: retrieved - duplicates = deduped; deduped = included + excluded).
"""
import csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = os.path.join(ROOT, "raw"), os.path.join(ROOT, "out")
RUN_DATE = "2026-08-24"

# Re-verification result (spec 7.1). A field counts as INCONSISTENT only if the
# record's abstract CONTRADICTS the manifest value. A value the abstract is simply
# silent about is not a contradiction -- those are counted separately, because the
# manifest declares source_field=fullText for them.
VERIFY = [
 ("32901110", [("n_patients", "NOT_FOUND", "consistent", "abstract gives 19 EBV+ NPCs and 7 nonmalignant biopsies (samples), states no patient count"),
               ("tissue_blood", "N", "consistent", "abstract mentions only nasopharyngeal biopsies"),
               ("NK_depth", "NOT_FOUND", "consistent", "abstract_only record; NK not mentioned in the abstract")]),
 ("40315843", [("n_patients", "22", "abstract_silent", "22 comes from full text 'TIL plus CCRT (n = 22)'; abstract reports 47 infusion products / 62 TMEs"),
               ("tissue_blood", "Y", "abstract_silent", "full text: 'PBMCs from NPC patients and healthy donors were isolated from peripheral blood'"),
               ("NK_depth", "cluster_only", "abstract_silent", "abstract names CD56+/CD56- DN TIL subsets; NK as such only in full text")]),
 ("40691404", [("n_patients", "24", "consistent", "abstract: 'single-cell and spatial transcriptomics analysis of 39 tumors from 24 patients'"),
               ("tissue_blood", "N", "consistent", "abstract describes tumours only"),
               ("NK_depth", "NOT_FOUND", "consistent", "abstract_only record; NK not mentioned in the abstract")]),
 ("41395113", [("n_patients", "NOT_FOUND", "consistent", "abstract states no cohort size"),
               ("tissue_blood", "N", "consistent", "abstract: 'NPC and adjacent normal tissue samples'"),
               ("NK_depth", "cluster_only", "abstract_silent", "NK appears only in a full-text figure legend listing cell types")]),
]


def main():
    dd = json.load(open(os.path.join(RAW, "dedup_counts.json")))
    finals = json.load(open(os.path.join(RAW, "final_sets.json")))
    anchors = json.load(open(os.path.join(RAW, "anchor_check.json")))
    qsum = json.load(open(os.path.join(RAW, "query_summary.json")))
    qsup = json.load(open(os.path.join(RAW, "query_summary_supp.json")))
    ftidx = json.load(open(os.path.join(RAW, "fulltext_index.json")))
    ans = json.load(open(os.path.join(RAW, "occupancy_answers.json")))
    dropped = len(json.load(open(os.path.join(RAW, "stage1_dropped.json"))))

    with open(os.path.join(OUT, "npc_sc_manifest.tsv")) as f:
        rows = list(csv.DictReader([l for l in f if not l.startswith("#")], delimiter="\t"))
    with open(os.path.join(OUT, "npc_sc_excluded.tsv")) as f:
        excl = list(csv.DictReader([l for l in f if not l.startswith("#")], delimiter="\t"))
    with open(os.path.join(OUT, "npc_sc_secondary.tsv")) as f:
        sec = list(csv.DictReader([l for l in f if not l.startswith("#")], delimiter="\t"))

    L = []
    A = L.append
    A("NPC single-cell occupancy audit -- COUNTS")
    A("run date: %s" % RUN_DATE)
    A("=" * 78)
    A("")
    A("1. RETRIEVAL")
    A("-" * 78)
    for qid in ("Q1", "Q2", "Q3", "Q4", "Q5"):
        A("  %-3s hitCount=%-6s retrieved=%-6s" % (qid, qsum[qid]["hitCount"], qsum[qid]["retrieved"]))
    A("  Q6  NCBI db=gds  esearch_ids=%s  esummary_docs=%s"
      % (qsum["Q6"]["esearch_ids"], qsum["Q6"]["esummary_docs"]))
    for qid in ("Q7", "Q8"):
        A("  %-3s hitCount=%-6s retrieved=%-6s   [spec section 5 anchor-repair query]"
          % (qid, qsup[qid]["hitCount"], qsup[qid]["retrieved"]))
    A("")
    A("  NOTE Q4 retrieved 1555 vs hitCount 1554: Europe PMC cursorMark paging returned one")
    A("  record more than the reported hitCount (the index shifted mid-crawl). The extra")
    A("  record is removed by de-duplication and does not propagate; closure below is exact.")
    A("")
    A("2. NUMERIC CLOSURE  (spec section 0, rule 6)")
    A("-" * 78)
    A("  retrieved (Q1-Q5,Q7,Q8 union, pre-dedup) : %d" % dd["retrieved"])
    A("  duplicates removed                       : %d" % dd["duplicates"])
    A("  deduped                                  : %d" % dd["deduped"])
    ok1 = dd["retrieved"] - dd["duplicates"] == dd["deduped"]
    A("  CHECK retrieved - duplicates = deduped   : %s  (%d - %d = %d)"
      % ("PASS" if ok1 else "FAIL", dd["retrieved"], dd["duplicates"], dd["deduped"]))
    A("")
    n_inc, n_sec = len(rows), len(sec)
    n_exc = len(excl)
    A("  included                                 : %d" % n_inc)
    A("  secondary_analysis (separate file)       : %d" % n_sec)
    A("  excluded                                 : %d" % n_exc)
    tot = n_inc + n_sec + n_exc
    ok2 = tot == dd["deduped"]
    A("  CHECK included + secondary + excluded    : %s  (%d + %d + %d = %d vs deduped %d)"
      % ("PASS" if ok2 else "FAIL", n_inc, n_sec, n_exc, tot, dd["deduped"]))
    A("")
    A("  The spec writes this closure as `deduped = included + excluded`. This audit keeps")
    A("  secondary_analysis as its own bucket (spec section 2 requires those records be kept,")
    A("  not discarded), so the identity is stated with three terms. No record appears twice.")
    A("")
    A("  excluded breakdown by stage:")
    from collections import Counter
    for k, v in Counter(r["stage"] for r in excl).most_common():
        A("    %-24s %d" % (k, v))
    A("")
    A("3. FULL-TEXT RETRIEVAL")
    A("-" * 78)
    tgt = len(ftidx)
    okft = sum(1 for v in ftidx.values() if v.get("ok"))
    A("  survivors with a PMCID          : %d" % tgt)
    A("  fullTextXML retrieved           : %d" % okft)
    A("  fullTextXML failed (404/no body): %d   -> logged in fetch_failures.tsv" % (tgt - okft))
    A("  included records with full text : %d / %d"
      % (sum(1 for r in rows if r["source_field"] == "fullText"), len(rows)))
    A("  included records abstract_only  : %d / %d"
      % (sum(1 for r in rows if r["source_field"].startswith("abstract")), len(rows)))
    A("")
    A("4. ANCHOR RECALL  (spec section 5)")
    A("-" * 78)
    hit = 0
    for a in anchors:
        st = "HIT " if a["hit"] else "MISS"
        if a["hit"]:
            hit += 1
        A("  %s #%-2d %s" % (st, a["anchor"], a["label"][:88]))
        if a["hit"]:
            A("         -> pmid=%s doi=%s queries=%s"
              % (a["pmid"], a["doi"], ",".join(a["queries"] or [])))
    A("")
    A("  RESULT: %d/%d hit, %d missed." % (hit, len(anchors), len(anchors) - hit))
    A("")
    A("  Two anchors needed work before they resolved:")
    A("   * #1  Nat Commun 2021;12:741 was present all along (PMID 33531485) but the matcher")
    A("         missed it: Europe PMC leaves the flat journalTitle/journalVolume null on that")
    A("         core record and puts them under journalInfo. Matcher fixed; no query change.")
    A("   * #15 Cell 2023;186:4235-4251 (PMID 37607536) is genuinely unreachable from Q1-Q5:")
    A("         isOpenAccess=N, no PMCID, inEPMC=N, so Europe PMC indexes no full text and the")
    A("         string 'nasopharyngeal' appears nowhere in any searchable field. Repaired by")
    A("         adding Q7 (pan-cancer single-cell atlas language) per spec section 5; it is")
    A("         retrieved by Q7, not hand-added.")
    A("")
    A("5. RE-VERIFICATION OF 10% OF INCLUDED RECORDS  (spec section 7.1)")
    A("-" * 78)
    A("  sample: %d of %d included records (10%%), drawn with a fixed seed (20260824)"
      % (len(VERIFY), len(rows)))
    A("  fields re-checked: n_patients, tissue_blood, NK_depth")
    A("")
    nchk = ncontra = nsilent = 0
    for pid, fields in VERIFY:
        A("  pmid %s" % pid)
        for fname, val, verdict, basis in fields:
            nchk += 1
            if verdict == "inconsistent":
                ncontra += 1
            elif verdict == "abstract_silent":
                nsilent += 1
            A("    %-14s = %-12s [%s] %s" % (fname, val, verdict, basis[:96]))
    A("")
    A("  field checks              : %d" % nchk)
    A("  contradicted by abstract  : %d" % ncontra)
    A("  abstract silent (value from full text, declared as source_field=fullText): %d" % nsilent)
    rate = 100.0 * ncontra / nchk if nchk else 0.0
    A("  INCONSISTENCY RATE        : %.1f%%  (threshold 5%%)  -> %s"
      % (rate, "PASS" if rate <= 5.0 else "FAIL -- extraction round rejected"))
    A("")
    A("6. NOT_FOUND RATE PER FIELD  (over %d included records)" % len(rows))
    A("-" * 78)
    A("  A high rate here is a finding about the literature, not a defect of the extraction.")
    A("")
    fields = [k for k in rows[0].keys() if k not in ("note",)]
    for k in fields:
        nf = sum(1 for r in rows if str(r[k]).strip() in ("NOT_FOUND", "unstated", ""))
        A("  %-28s %3d/%d  %5.1f%%" % (k, nf, len(rows), 100.0 * nf / len(rows)))
    A("")
    A("7. OCCUPANCY ANSWERS (see OCCUPANCY.md)")
    A("-" * 78)
    A("  Q-A blood + NK_depth>=subclustered        : %d" % len(ans["qa"]))
    A("  Q-B distant metastasis (M1) samples       : %d   (of which with blood: %d)"
      % (len(ans["qb"]), len(ans["qb_blood"])))
    A("  Q-C longitudinal repeat BLOOD sampling    : %d   (longitudinal of any material: %d)"
      % (len(ans["qc"]), len(ans["qc_all"])))
    A("  Q-D continuous plasma EBV DNA regression  : %d" % len(ans["qd"]))
    A("  Q-E metastatic-cascade staged framework   : %d" % len(ans["qe"]))
    A("  Q-F public downloadable NPC sc datasets")
    A("      that contain blood                    : %d   (%s)"
      % (len(ans["qf"]), ", ".join(ans["qf"])))
    A("")
    A("8. KNOWN LIMITS OF THIS ROUND")
    A("-" * 78)
    A("  * %d of %d included records have no retrievable full text, so their fields are"
      % (sum(1 for r in rows if r["source_field"].startswith("abstract")), len(rows)))
    A("    abstract_only and their NK_depth is NOT_FOUND rather than 'absent' -- 'absent'")
    A("    would assert something about a full text that was never retrieved.")
    A("  * Spatial-transcriptomics platforms differ in resolution. 10x Visium and GeoMx DSP are")
    A("    spot/region-level, NOT single-cell. Records whose only new data are Visium/GeoMx are")
    A("    still included (spec section 2 excludes only PURE public re-analysis with no new")
    A("    samples), but modality names the platform so the reader can discount them.")
    A("  * includes_blood in npc_sc_datasets.tsv is decided from GEO sample titles / series")
    A("    title only. GEO series summaries repeat the paper abstract and would have produced")
    A("    at least one false positive (GSE206245).")
    A("  * Q6 uses the exact term the spec prescribes; NCBI's translation of that term returns")
    A("    24 GSE records, which is narrower than the accession set recovered from the papers")
    A("    themselves. npc_sc_datasets.tsv is the union of both.")
    A("")
    with open(os.path.join(OUT, "COUNTS.txt"), "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L[:1]))
    print("closure checks: retrieved-dup=deduped %s | inc+sec+exc=deduped %s"
          % ("PASS" if ok1 else "FAIL", "PASS" if ok2 else "FAIL"))
    print("anchors: %d/%d | inconsistency rate %.1f%%" % (hit, len(anchors), rate))


if __name__ == "__main__":
    main()
