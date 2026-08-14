#!/usr/bin/env python3
"""Phase C.1/C.2 — pull the small mtDNA tables for every GSE302113 library.

Deliberately does NOT pull the 36 GB of fragments: the per-library
`variant_stats` (mgatk output over all 16,569 chrM positions) and
`cell_heteroplasmic_df` (per-cell heteroplasmy at the retained variants)
carry everything the preregistered power check needs except the NK cell
counts, which require cell typing.
"""
import gzip
import os
import statistics
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DL = os.path.join(HERE, "phaseC", "dl")
UA = {"User-Agent": "nk-lineage-screen/1.0 (research use)"}

# GSM -> (library, tissue, sorting, donor) parsed from the series matrix
SAMPLES = []


def load_matrix():
    import re
    path = os.path.join(HERE, "phaseC", "meta", "GSE302113_series_matrix.txt.gz")
    fields = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.partition("\t")
            vals = re.findall(r'"([^"]*)"', rest)
            fields.setdefault(key.strip(), []).append(vals)
    gsms = fields["!Sample_geo_accession"][0]
    titles = fields["!Sample_title"][0]
    chars = fields["!Sample_characteristics_ch1"]
    tissue = next(c for c in chars if c[0].startswith("tissue:"))
    sorting = next(c for c in chars if c[0].startswith("sorting:"))
    donor = next(c for c in chars if c[0].startswith("donor:"))
    diag = next(c for c in chars if c[0].startswith("diagnosis:"))
    out = []
    for i, g in enumerate(gsms):
        out.append({
            "gsm": g, "library": titles[i],
            "tissue": tissue[i].split(": ", 1)[1],
            "sorting": sorting[i].split(": ", 1)[1],
            "donor": donor[i].split(": ", 1)[1],
            "diagnosis": diag[i].split(": ", 1)[1],
        })
    return out


def fetch(gsm, fname, tries=4):
    dest = os.path.join(DL, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return dest
    stub = gsm[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/suppl/{fname}"
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as fh:
                fh.write(r.read())
            return dest
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                sys.stderr.write(f"FAIL {fname}: {e}\n")
                return None
            time.sleep(2 ** i)
    return None


def stats_for(vs_path, het_path):
    """Return median chrM coverage, informative-variant count, retained variants, n cells."""
    cov = {}
    inf = 0
    with gzip.open(vs_path, "rt") as f:
        idx = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            p = line.rstrip("\n").split("\t")
            cov[int(p[idx["position"]])] = float(p[idx["mean_coverage"]])
            sc = p[idx["strand_correlation"]]
            try:
                sc = float(sc) if sc else -1.0
            except ValueError:
                sc = -1.0
            vmr = float(p[idx["vmr"]] or 0)
            nd = int(p[idx["n_cells_conf_detected"]] or 0)
            # preregistered informative-variant definition (mgatk defaults)
            if sc >= 0.65 and vmr >= 0.01 and nd >= 5:
                inf += 1
    ncells = nvars = 0
    if het_path and os.path.exists(het_path):
        with gzip.open(het_path, "rt") as f:
            nvars = len(f.readline().rstrip("\n").split("\t")) - 1
            ncells = sum(1 for _ in f)
    c = sorted(cov.values())
    return {
        "median_chrM_coverage": round(statistics.median(c), 1),
        "p10_chrM_coverage": round(c[len(c) // 10], 1),
        "chrM_positions": len(cov),
        "n_informative_variants": inf,
        "n_retained_variants_deposited": nvars,
        "n_cells": ncells,
    }


def main():
    os.makedirs(DL, exist_ok=True)
    samples = load_matrix()
    out = os.path.join(HERE, "out", "T4_power_GSE302113_libraries.tsv")
    cols = ["dataset_id", "gsm", "library", "donor", "diagnosis", "tissue", "sorting",
            "n_cells", "median_chrM_coverage", "p10_chrM_coverage", "chrM_positions",
            "n_informative_variants", "n_retained_variants_deposited", "chrM_data_present"]
    rows = []
    for i, s in enumerate(samples, 1):
        vs = fetch(s["gsm"], f"{s['gsm']}_{s['library']}.variant_stats.tsv.gz")
        het = fetch(s["gsm"], f"{s['gsm']}_{s['library']}.cell_heteroplasmic_df.tsv.gz")
        if not vs:
            rows.append({**s, "dataset_id": "GSE302113", "chrM_data_present": "FETCH_FAILED"})
            continue
        st = stats_for(vs, het)
        rows.append({**s, **st, "dataset_id": "GSE302113", "chrM_data_present": "TRUE"})
        sys.stderr.write(f"[{i}/{len(samples)}] {s['gsm']} {s['donor']:10s} {s['tissue']:14s} "
                         f"cov={st['median_chrM_coverage']:6.1f} inf={st['n_informative_variants']:4d} "
                         f"cells={st['n_cells']}\n")
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
