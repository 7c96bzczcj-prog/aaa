#!/usr/bin/env python3
"""Phase C.2 — per-library NK cell counts for GSE302113.

No cell-type annotations were deposited with GSE302113, so the NK counts the
preregistered power check needs have to be derived. This streams each library's
fragments file straight from GEO (nothing is written to disk), counts fragments
over marker gene bodies + 2 kb promoters for the authors' called cells, and
classifies each cell by which lineage's accessibility is most enriched.

The classifier is validated against the dataset's own sorted libraries:
GSM9096509 is CD56+ sorted PBMC and GSM9096510 is its CD56- counterpart, so the
NK fraction must be high in the former and low in the latter. That control is
reported by `--control` and its result is recorded in DECISIONS.

Usage:
  phaseC_genescores.py --gsm GSM9096509 [GSM...] [--jobs 3]
  phaseC_genescores.py --control
"""
import argparse
import collections
import glob
import gzip
import json
import math
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PC = os.path.join(HERE, "phaseC")
GENOME_BP = 3.1e9
MIN_FRAGS = 1000          # cells below this are not typed (too shallow)
MIN_MARGIN = 1.0          # argmax z must beat runner-up by this much


def libname(gsm):
    for f in glob.glob(os.path.join(PC, "dl", f"{gsm}_*.variant_stats.tsv.gz")):
        return os.path.basename(f)[len(gsm) + 1:-len(".variant_stats.tsv.gz")]
    raise SystemExit(f"no downloaded tables for {gsm}")


def barcodes(gsm, lib):
    out = os.path.join(PC, "bc", f"{gsm}.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        src = os.path.join(PC, "dl", f"{gsm}_{lib}.cell_heteroplasmic_df.tsv.gz")
        with gzip.open(src, "rt") as f, open(out, "w") as o:
            f.readline()
            for line in f:
                o.write(line.split("\t", 1)[0] + "\n")
    with open(out) as f:
        return out, sum(1 for _ in f)


def stream_counts(gsm, lib):
    """Stream fragments from GEO through the awk counter. Returns raw awk stdout path."""
    dest = os.path.join(PC, "gs", f"{gsm}.counts.tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return dest
    bcf, _ = barcodes(gsm, lib)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
           f"{gsm}_{lib}_fragments.tsv.gz")
    cmd = (f"curl -sS --retry 4 --retry-delay 3 --max-time 3000 '{url}' | gunzip -c | "
           f"awk -v BCFILE='{bcf}' -v BEDFILE='{PC}/meta/markers.bed' "
           f"-f '{HERE}/scripts/genescore.awk' > '{dest}.part' && mv '{dest}.part' '{dest}'")
    r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"[{gsm}] FAILED: {r.stderr[-400:]}\n")
        return None
    return dest


def classify(counts_path, markers):
    lin_genes = collections.defaultdict(list)
    lin_bp = collections.defaultdict(int)
    for g, v in markers.items():
        lin_genes[v["lineage"]].append(g)
        lin_bp[v["lineage"]] += v["end"] - v["start"]
    lineages = sorted(lin_genes)

    tot = {}
    mark = collections.defaultdict(dict)
    with open(counts_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "T":
                tot[p[1]] = int(p[2])
            else:
                mark[p[1]][p[2]] = int(p[3])

    cells = [c for c, t in tot.items() if t >= MIN_FRAGS]
    enr = {L: {} for L in lineages}
    for c in cells:
        m = mark.get(c, {})
        for L in lineages:
            raw = sum(m.get(g, 0) for g in lin_genes[L])
            exp = tot[c] * lin_bp[L] / GENOME_BP
            enr[L][c] = (raw / exp) if exp > 0 else 0.0

    # z-score each lineage across cells so lineage-level baseline accessibility
    # differences do not decide the argmax
    z = {L: {} for L in lineages}
    for L in lineages:
        vals = [enr[L][c] for c in cells]
        mu = statistics.mean(vals) if vals else 0.0
        sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        for c in cells:
            z[L][c] = ((enr[L][c] - mu) / sd) if sd > 0 else 0.0

    calls = {}
    for c in cells:
        ranked = sorted(((z[L][c], L) for L in lineages), reverse=True)
        top, second = ranked[0], ranked[1]
        calls[c] = top[1] if (top[0] - second[0]) >= MIN_MARGIN and top[0] > 0 else "ambiguous"
    return calls, tot, len(cells)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gsm", nargs="*", default=[])
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--control", action="store_true")
    a = ap.parse_args()

    markers = json.load(open(os.path.join(PC, "meta", "markers.json")))
    gsms = a.gsm or (["GSM9096509", "GSM9096510"] if a.control else [])
    if not gsms:
        raise SystemExit("nothing to do")

    os.makedirs(os.path.join(PC, "gs"), exist_ok=True)
    # streaming is network bound, so run a few at once
    from concurrent.futures import ThreadPoolExecutor
    libs = {g: libname(g) for g in gsms}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        list(ex.map(lambda g: stream_counts(g, libs[g]), gsms))

    rows = []
    for g in gsms:
        p = os.path.join(PC, "gs", f"{g}.counts.tsv")
        if not os.path.exists(p):
            sys.stderr.write(f"[{g}] no counts\n")
            continue
        calls, tot, n_typed = classify(p, markers)
        cnt = collections.Counter(calls.values())
        rows.append((g, len(tot), n_typed, cnt))
        frac = cnt["NK"] / n_typed if n_typed else 0
        print(f"{g}  cells={len(tot):6d} typed={n_typed:6d}  NK={cnt['NK']:5d} ({frac:5.1%})  "
              + " ".join(f"{k}={v}" for k, v in cnt.most_common() if k != "NK"))
        with open(os.path.join(PC, "gs", f"{g}.calls.tsv"), "w") as f:
            for c, L in sorted(calls.items()):
                f.write(f"{c}\t{L}\t{tot[c]}\n")

    if a.control and len(rows) == 2:
        pos = rows[0][3]["NK"] / max(1, rows[0][2])
        neg = rows[1][3]["NK"] / max(1, rows[1][2])
        print(f"\nCONTROL  CD56+ NK fraction = {pos:.1%}   CD56- NK fraction = {neg:.1%}   "
              f"ratio = {pos / neg if neg else math.inf:.1f}x")


if __name__ == "__main__":
    main()
