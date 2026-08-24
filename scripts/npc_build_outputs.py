#!/usr/bin/env python3
"""Build the deliverables in ./out (spec v1.0 section 6)."""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = os.path.join(ROOT, "raw"), os.path.join(ROOT, "out")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from npc_adjudication import VERDICT, INCLUDE_BASIS  # noqa: E402
from npc_manifest_data import M  # noqa: E402
from npc_screen_decisions import EXCLUDE as STAGE2_EXCLUDE  # noqa: E402
from npc_modality import modalities  # noqa: E402

RUN_DATE = "2026-08-24"

HEADER = ["pmid", "doi", "year", "journal", "first_author", "last_author",
          "corresponding_affiliation", "modality", "tissue_primary", "tissue_LN_met",
          "tissue_distant_met", "tissue_normal_NP", "tissue_blood", "n_patients",
          "n_samples", "n_cells", "treatment_status", "longitudinal", "M_stage",
          "NK_depth", "NK_claim", "EBV_DNA_link", "accession",
          "accession_includes_blood", "source_field", "source_snippet", "note"]


def clean(s, n=200):
    if s is None:
        return "NOT_FOUND"
    s = re.sub(r"<[^>]+>", "", str(s))
    s = re.sub(r"[\t\r\n]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s[:n] if n else s


def authors(rec):
    al = (rec.get("authorList") or {}).get("author") or []
    if not al:
        s = rec.get("authorString") or ""
        parts = [p.strip() for p in s.rstrip(".").split(",") if p.strip()]
        if not parts:
            return "NOT_FOUND", "NOT_FOUND"
        return parts[0], parts[-1]
    def nm(a):
        return a.get("fullName") or ((a.get("lastName") or "") + " " + (a.get("initials") or "")).strip() or "NOT_FOUND"
    return nm(al[0]), nm(al[-1])


def journal_of(rec):
    ji = (rec.get("journalInfo") or {}).get("journal") or {}
    return (rec.get("journalTitle") or ji.get("medlineAbbreviation")
            or ji.get("title") or "NOT_FOUND")


def main():
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    all_dedup = json.load(open(os.path.join(RAW, "deduped.json")))
    dropped = json.load(open(os.path.join(RAW, "stage1_dropped.json")))
    cand = json.load(open(os.path.join(RAW, "candidate_ids.json")))
    counts_dd = json.load(open(os.path.join(RAW, "dedup_counts.json")))

    included = [i for i in cand if i not in VERDICT]
    secondary = [i for i, v in VERDICT.items() if v[0] == "secondary"]
    excl_stage3 = [i for i, v in VERDICT.items() if v[0] == "excluded"]

    # ---------------- manifest -------------------------------------------------
    rows = []
    for rid in sorted(included, key=lambda x: -(int(recs[x].get("pubYear") or 0))):
        r, m = recs[rid], M[rid]
        fa, la = authors(r)
        rows.append({
            "pmid": r.get("pmid") or "NOT_FOUND",
            "doi": r.get("doi") or "NOT_FOUND",
            "year": r.get("pubYear") or "NOT_FOUND",
            "journal": clean(journal_of(r), 80),
            "first_author": clean(fa, 60), "last_author": clean(la, 60),
            "corresponding_affiliation": clean(r.get("affiliation") or "NOT_FOUND", 300),
            "modality": m["modality"],
            "tissue_primary": m["primary"], "tissue_LN_met": m["ln"],
            "tissue_distant_met": m["dist"], "tissue_normal_NP": m["normal"],
            "tissue_blood": m["blood"],
            "n_patients": m["n_pat"], "n_samples": m["n_samp"], "n_cells": m["n_cells"],
            "treatment_status": m["treat"], "longitudinal": m["longit"],
            "M_stage": m["mstage"], "NK_depth": m["nk_depth"],
            "NK_claim": clean(m["nk_claim"], 200), "EBV_DNA_link": m["ebv"],
            "accession": m["acc"], "accession_includes_blood": m["acc_blood"],
            "source_field": m["src"], "source_snippet": clean(m["snip"], 200),
            "note": clean(m.get("note"), 400),
        })
    with open(os.path.join(OUT, "npc_sc_manifest.tsv"), "w") as f:
        f.write("# NPC single-cell occupancy audit -- included records\n")
        f.write("# run date: %s | source: Europe PMC REST + NCBI E-utilities\n" % RUN_DATE)
        f.write("# NOT_FOUND = not stated in the retrieved text; never an assertion of absence\n")
        f.write("\t".join(HEADER) + "\n")
        for row in rows:
            f.write("\t".join(str(row[h]) for h in HEADER) + "\n")

    # ---------------- secondary ------------------------------------------------
    with open(os.path.join(OUT, "npc_sc_secondary.tsv"), "w") as f:
        f.write("# Public-data re-analyses of NPC single-cell data (no new single-cell samples)\n")
        f.write("# run date: %s\n" % RUN_DATE)
        f.write("pmid\tdoi\tyear\tjournal\ttitle\treason\tsource_snippet\n")
        for rid in sorted(secondary, key=lambda x: -(int(recs[x].get("pubYear") or 0))):
            r = recs[rid]
            v = VERDICT[rid]
            f.write("\t".join([r.get("pmid") or "NOT_FOUND", r.get("doi") or "NOT_FOUND",
                               r.get("pubYear") or "NOT_FOUND", clean(journal_of(r), 60),
                               clean(r.get("title"), 160), v[1], clean(v[2], 200)]) + "\n")

    # ---------------- excluded -------------------------------------------------
    with open(os.path.join(OUT, "npc_sc_excluded.tsv"), "w") as f:
        f.write("# Excluded records with reason (spec v1.0 section 2)\n")
        f.write("# run date: %s\n" % RUN_DATE)
        f.write("stage\tpmid\tdoi\tyear\ttitle\treason\tbasis\n")
        # stage 1: no single-cell modality term anywhere in retrieved title/abstract
        for r in dropped:
            f.write("\t".join(["1_no_modality_term", r.get("pmid") or "NOT_FOUND",
                               r.get("doi") or "NOT_FOUND", r.get("pubYear") or "NOT_FOUND",
                               clean(r.get("title"), 160), r["_drop_reason"],
                               "no single-cell/CyTOF/CITE-seq/spatial term in retrieved title or abstract"]) + "\n")
        # stage 2: read and failed an inclusion criterion
        pool = json.load(open(os.path.join(RAW, "reviewpool_ids.json")))
        idx_by_pos = {}
        surv = [r for r in json.load(open(os.path.join(RAW, "survivors.json")))
                if r["_flags"]["has_nasopharyngeal_term"]]
        surv = [r for r in surv if modalities((r.get("title") or "") + " " + (r.get("abstractText") or ""))]
        surv.sort(key=lambda r: (-(int(r.get("pubYear") or 0)), r.get("pmid") or ""))
        for i, r in enumerate(surv, 1):
            idx_by_pos[i] = r
        for i, (code, basis) in sorted(STAGE2_EXCLUDE.items()):
            r = idx_by_pos[i]
            f.write("\t".join(["2_screening", r.get("pmid") or "NOT_FOUND",
                               r.get("doi") or "NOT_FOUND", r.get("pubYear") or "NOT_FOUND",
                               clean(r.get("title"), 160), code, clean(basis, 250)]) + "\n")
        # stage 3: failed on full-text adjudication
        for rid in excl_stage3:
            r = recs[rid]
            v = VERDICT[rid]
            f.write("\t".join(["3_adjudication", r.get("pmid") or "NOT_FOUND",
                               r.get("doi") or "NOT_FOUND", r.get("pubYear") or "NOT_FOUND",
                               clean(r.get("title"), 160), v[1], clean(v[2], 250)]) + "\n")
        # records with a nasopharyngeal term but no single-cell modality term
        nomod = [r for r in json.load(open(os.path.join(RAW, "survivors.json")))
                 if r["_flags"]["has_nasopharyngeal_term"]
                 and not modalities((r.get("title") or "") + " " + (r.get("abstractText") or ""))]
        for r in nomod:
            f.write("\t".join(["2_no_true_modality", r.get("pmid") or "NOT_FOUND",
                               r.get("doi") or "NOT_FOUND", r.get("pubYear") or "NOT_FOUND",
                               clean(r.get("title"), 160), "no_single_cell_genomics_modality",
                               "'single cell' in retrieved text refers to a comet assay, cell suspension, clonal subline, Raman/impedance single-cell measurement or FNA smear pattern"]) + "\n")
        # pan-cancer pool with no retrievable evidence of NPC samples
        ft = json.load(open(os.path.join(RAW, "fulltext_index.json")))
        for r in json.load(open(os.path.join(RAW, "survivors.json"))):
            if r["_flags"]["has_nasopharyngeal_term"] or r["id"] in cand:
                continue
            pm = r.get("pmcid")
            has = ft.get(pm, {}).get("npc_full") if pm else None
            reason = ("pan_cancer_record_full_text_has_no_NPC_mention" if has is False
                      else "pan_cancer_record_no_retrievable_evidence_of_NPC_samples")
            f.write("\t".join(["2_pan_cancer_no_NPC", r.get("pmid") or "NOT_FOUND",
                               r.get("doi") or "NOT_FOUND", r.get("pubYear") or "NOT_FOUND",
                               clean(r.get("title"), 160), reason,
                               "retrieved from Q7/Q8 (anchor-repair queries); retrieved text does not evidence human NPC samples at single-cell resolution"]) + "\n")

    print("manifest rows      :", len(rows))
    print("secondary rows     :", len(secondary))
    print("stage1 dropped     :", len(dropped))
    print("stage3 excluded    :", len(excl_stage3))
    json.dump({"included": included, "secondary": secondary, "excluded_stage3": excl_stage3,
               "run_date": RUN_DATE, "dedup": counts_dd, "deduped_total": len(all_dedup)},
              open(os.path.join(RAW, "final_sets.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
