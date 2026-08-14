#!/usr/bin/env python3
"""Phase C.2 — the preregistered power table for GSE302113.

One row per donor x compartment, which is the unit the preregistered rule
consumes ("NK cells per donor per compartment < 30 -> that cell is unusable"),
so libraries belonging to the same compartment of the same donor are summed
(G6: the summary is computed on the same basis as the rule that consumes it).

Informative variants are counted PER DONOR by union across that donor's
libraries, not per library. That is both the preregistered mgatk definition
(strand correlation >= 0.65, VMR >= 0.01, detected in >= 5 cells) and the
procedure the dataset's own authors use for cross-tissue clone calling: a
variant passing filters in any one sample of a donor is treated as informative
in all samples of that donor. Per-library counts would be the wrong denominator
for a cross-compartment question.
"""
import collections
import glob
import gzip
import json
import os
import statistics

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PC = os.path.join(HERE, "phaseC")

COMPARTMENT = {
    "Lung tumor": "tumour",
    "Lung normal": "adjacent_normal",
    "Ovarian tumor": "tumour",
    "Omentum met": "metastasis",
    "PBMC": "blood",
}
MIN_NK = 30          # preregistered
MIN_DONORS = 3       # preregistered
MIN_COV = 20.0       # preregistered


def load_samples():
    import re
    path = os.path.join(PC, "meta", "GSE302113_series_matrix.txt.gz")
    fields = collections.defaultdict(list)
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_"):
                key, _, rest = line.partition("\t")
                fields[key.strip()].append(re.findall(r'"([^"]*)"', rest))
    gsms = fields["!Sample_geo_accession"][0]
    titles = fields["!Sample_title"][0]
    ch = fields["!Sample_characteristics_ch1"]
    pick = lambda p: next(c for c in ch if c[0].startswith(p))  # noqa: E731
    t, s, d, dg = pick("tissue:"), pick("sorting:"), pick("donor:"), pick("diagnosis:")
    return [{"gsm": g, "library": titles[i], "tissue": t[i].split(": ", 1)[1],
             "sorting": s[i].split(": ", 1)[1], "donor": d[i].split(": ", 1)[1],
             "diagnosis": dg[i].split(": ", 1)[1]} for i, g in enumerate(gsms)]


def informative_by_library(gsm, lib):
    """Set of variants passing the preregistered mgatk filters in this library."""
    p = os.path.join(PC, "dl", f"{gsm}_{lib}.variant_stats.tsv.gz")
    keep, cov = set(), {}
    with gzip.open(p, "rt") as f:
        idx = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            q = line.rstrip("\n").split("\t")
            pos = int(q[idx["position"]])
            cov[pos] = float(q[idx["mean_coverage"]])
            sc = q[idx["strand_correlation"]]
            try:
                sc = float(sc) if sc else -1.0
            except ValueError:
                sc = -1.0
            if sc >= 0.65 and float(q[idx["vmr"]] or 0) >= 0.01 and int(q[idx["n_cells_conf_detected"]] or 0) >= 5:
                # the authors exclude chrM:307-314, a homopolymeric region that
                # creates spurious clonal links; excluded here for the same reason
                if not (307 <= pos <= 314):
                    keep.add(q[idx["variant"]])
    return keep, statistics.median(cov.values())


def main():
    samples = load_samples()
    calls_dir = os.path.join(PC, "calls")

    nk_by_lib, typed_by_lib = {}, {}
    for p in glob.glob(os.path.join(calls_dir, "*.calls.tsv")):
        g = os.path.basename(p).split(".")[0]
        c = collections.Counter()
        with open(p) as f:
            for line in f:
                c[line.split("\t")[1]] += 1
        nk_by_lib[g] = c.get("NK", 0)
        typed_by_lib[g] = sum(c.values())

    per_donor_vars = collections.defaultdict(set)
    lib_cov = {}
    for s in samples:
        try:
            keep, med = informative_by_library(s["gsm"], s["library"])
        except FileNotFoundError:
            continue
        per_donor_vars[s["donor"]] |= keep
        lib_cov[s["gsm"]] = med

    groups = collections.defaultdict(list)
    for s in samples:
        groups[(s["donor"], COMPARTMENT.get(s["tissue"], s["tissue"]))].append(s)

    rows = []
    for (donor, comp), ss in sorted(groups.items()):
        gsms = [s["gsm"] for s in ss]
        typed = sum(typed_by_lib.get(g, 0) for g in gsms)
        nk = sum(nk_by_lib.get(g, 0) for g in gsms)
        covs = [lib_cov[g] for g in gsms if g in lib_cov]
        untyped = [g for g in gsms if g not in typed_by_lib]
        rows.append({
            "dataset_id": "GSE302113",
            "donor_id": donor,
            "diagnosis": ss[0]["diagnosis"],
            "compartment": comp,
            "n_libraries": len(gsms),
            "sorting": ";".join(sorted({s["sorting"] for s in ss})),
            "n_cells_typed": typed,
            "n_NK_cells": nk if not untyped else "PENDING",
            "median_chrM_coverage": round(statistics.median(covs), 1) if covs else "NA",
            "min_chrM_coverage": round(min(covs), 1) if covs else "NA",
            "n_informative_variants": len(per_donor_vars[donor]),
            "coverage_pass": "TRUE" if covs and min(covs) >= MIN_COV else "FALSE",
            "nk_pass": ("PENDING" if untyped else ("TRUE" if nk >= MIN_NK else "FALSE")),
            "libraries": ";".join(gsms),
        })

    cols = ["dataset_id", "donor_id", "diagnosis", "compartment", "n_libraries", "sorting",
            "n_cells_typed", "n_NK_cells", "median_chrM_coverage", "min_chrM_coverage",
            "n_informative_variants", "coverage_pass", "nk_pass", "libraries"]
    out = os.path.join(HERE, "out", "T4_power.tsv")
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {out} ({len(rows)} donor x compartment rows)")

    # preregistered dataset-level verdicts, per claim-relevant compartment pairing
    def verdict(name, pair, donors):
        ok = []
        for d in donors:
            a = next((r for r in rows if r["donor_id"] == d and r["compartment"] == pair[0]), None)
            b = next((r for r in rows if r["donor_id"] == d and r["compartment"] == pair[1]), None)
            if not a or not b:
                continue
            if "PENDING" in (a["nk_pass"], b["nk_pass"]):
                return f"{name}: PENDING (cell typing incomplete)"
            if a["nk_pass"] == "TRUE" and b["nk_pass"] == "TRUE":
                ok.append(d)
        return (f"{name}: {len(ok)} donors with both compartments >= {MIN_NK} NK "
                f"({'PASS' if len(ok) >= MIN_DONORS else 'FAIL'}, need >= {MIN_DONORS}) {ok}")

    lung = sorted({r["donor_id"] for r in rows if r["donor_id"].startswith("SU-L")})
    ovar = sorted({r["donor_id"] for r in rows if r["donor_id"].startswith("SU-O")})
    print()
    print(verdict("X1 adjacent_normal vs tumour (lung)", ("adjacent_normal", "tumour"), lung))
    print(verdict("X2 blood vs tumour (lung)", ("blood", "tumour"), lung))
    print(verdict("X2 blood vs tumour (ovarian)", ("blood", "tumour"), ovar))


if __name__ == "__main__":
    main()
