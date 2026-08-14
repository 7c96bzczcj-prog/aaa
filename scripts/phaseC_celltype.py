#!/usr/bin/env python3
"""Phase C.2 — NK cell counts per GSE302113 library.

Two stages, deliberately separated so the annotation rule can be fixed without
re-streaming 36 GB of fragments:

  stage 1 (`--cluster`): stream a library's fragments from GEO, build 5 kb tiles,
      TF-IDF -> LSI -> KMeans, and cache per-cluster marker counts. Nothing is
      written to disk except the small per-cluster summary.
  stage 2 (`--annotate`): pool every cluster from every library, z-score each
      lineage's enrichment ACROSS THE POOL, and label each cluster by argmax.

The pooling in stage 2 is the point. Scoring a cluster against the other clusters
of its own library assumes that library is a mixture of lineages; for a sorted
library it is not, and that assumption is what made the first two attempts fail
the control (DEC-07, DEC-08). The pooled cluster set does span all lineages.

Acceptance test (`--control`): GSM9096509 is CD56+ sorted PBMC and GSM9096510 is
the matched CD56- fraction, same donor. CD56+ must come out NK-high.
"""
import argparse
import collections
import glob
import gzip
import json
import os
import subprocess
import sys
from array import array

import numpy as np
from scipy.sparse import coo_matrix
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PC = os.path.join(HERE, "phaseC")
CL = os.path.join(PC, "clusters")
TILE = 5000
MIN_FRAGS = 1000
MIN_TILE_CELLS = 0.01
MAX_TILE_CELLS = 0.95
N_LSI = 30
CELLS_PER_CLUSTER = 120
MAX_CLUST = 20


def libname(gsm):
    for f in glob.glob(os.path.join(PC, "dl", f"{gsm}_*.variant_stats.tsv.gz")):
        return os.path.basename(f)[len(gsm) + 1:-len(".variant_stats.tsv.gz")]
    return None


def barcode_file(gsm, lib):
    out = os.path.join(PC, "bc", f"{gsm}.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        src = os.path.join(PC, "dl", f"{gsm}_{lib}.cell_heteroplasmic_df.tsv.gz")
        with gzip.open(src, "rt") as f, open(out, "w") as o:
            f.readline()
            for line in f:
                o.write(line.split("\t", 1)[0] + "\n")
    return out


AWK = r'''
BEGIN{FS="\t"; ci=0
  while((getline l < BCFILE)>0){ if(!(l in cid)){ cid[l]=ci; ci++ } }
  n=0
  while((getline l < BEDFILE)>0){ split(l,a,"\t"); n++; ms[n]=a[2]+0; me[n]=a[3]+0; mg[n]=a[4]
    b0=int(ms[n]/100000); b1=int(me[n]/100000)
    for(b=b0;b<=b1;b++){ k=a[1]":"b; bin[k]=bin[k]" "n } }
}
/^#/{next}
{ bc=$4; if(!(bc in cid)) next
  c=cid[bc]; fs=$2+0; fe=$3+0
  print "F\t" c "\t" $1 "\t" int(fs/TILEBP)
  k1=$1":"int(fs/100000); k2=$1":"int(fe/100000); hit=""
  if(k1 in bin) hit=bin[k1]
  if(k2!=k1 && (k2 in bin)) hit=hit" "bin[k2]
  if(hit=="") next
  m=split(hit,idx," "); delete seen
  for(j=1;j<=m;j++){ i=idx[j]+0; if(i==0||(i in seen)) continue; seen[i]=1
    if(fs<me[i] && fe>ms[i]) print "M\t" c "\t" mg[i] }
}'''


def stream(gsm, lib, bcf):
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
           f"{gsm}_{lib}_fragments.tsv.gz")
    cmd = (f"curl -sS --retry 4 --retry-delay 3 --max-time 4000 '{url}' | gunzip -c | "
           f"awk -v BCFILE='{bcf}' -v BEDFILE='{PC}/meta/markers.bed' -v TILEBP={TILE} '{AWK}'")
    proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 22)
    rows, cols = array("i"), array("i")
    tile_id = {}
    totals = collections.Counter()
    marker = collections.defaultdict(collections.Counter)
    for line in proc.stdout:
        p = line.split("\t")
        if p[0] == "F":
            c = int(p[1])
            key = p[2] + p[3]
            t = tile_id.get(key)
            if t is None:
                t = len(tile_id)
                tile_id[key] = t
            rows.append(c)
            cols.append(t)
            totals[c] += 1
        else:
            marker[int(p[1])][p[2].rstrip("\n")] += 1
    proc.stdout.close()
    err = proc.stderr.read()
    if proc.wait() != 0:
        sys.stderr.write(f"[{gsm}] stream failed: {err[-300:]}\n")
        return None
    return rows, cols, len(tile_id), totals, marker


def cluster_one(gsm, lib):
    dest = os.path.join(CL, f"{gsm}.json")
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return json.load(open(dest))
    bcf = barcode_file(gsm, lib)
    bcs = [l.strip() for l in open(bcf)]
    got = stream(gsm, lib, bcf)
    if got is None:
        return None
    rows, cols, n_tiles, totals, marker = got
    n_cells = len(bcs)
    M = coo_matrix((np.ones(len(rows), dtype=np.float32),
                    (np.frombuffer(rows, dtype=np.int32), np.frombuffer(cols, dtype=np.int32))),
                   shape=(n_cells, n_tiles)).tocsr()
    depth = np.asarray([totals.get(i, 0) for i in range(n_cells)])
    keep = np.where(depth >= MIN_FRAGS)[0]
    if len(keep) < 60:
        r = {"gsm": gsm, "status": "TOO_FEW_CELLS", "n_cells": n_cells, "n_typed": int(len(keep))}
        json.dump(r, open(dest, "w"))
        return r
    B = (M[keep] > 0).astype(np.float32)
    per_tile = np.asarray(B.sum(axis=0)).ravel()
    kt = np.where((per_tile >= MIN_TILE_CELLS * len(keep)) &
                  (per_tile <= MAX_TILE_CELLS * len(keep)))[0]
    B = B[:, kt]
    tf = normalize(B, norm="l1", axis=1)
    idf = np.log(1 + B.shape[0] / (np.asarray(B.sum(axis=0)).ravel() + 1))
    X = tf.multiply(idf).tocsr()
    X.data = np.log1p(X.data * 1e4)
    k = min(N_LSI, min(X.shape) - 1)
    L = TruncatedSVD(n_components=k, random_state=0).fit_transform(X)[:, 1:]
    L = normalize(L)
    nclust = int(min(MAX_CLUST, max(2, len(keep) // CELLS_PER_CLUSTER)))
    lab = KMeans(n_clusters=nclust, n_init=10, random_state=0).fit_predict(L)

    cl_tot = collections.Counter()
    cl_mark = collections.defaultdict(collections.Counter)
    cl_n = collections.Counter()
    for j, ci in enumerate(keep):
        c = int(lab[j])
        cl_n[c] += 1
        cl_tot[c] += totals.get(int(ci), 0)
        for g, v in marker.get(int(ci), {}).items():
            cl_mark[c][g] += v
    r = {"gsm": gsm, "status": "OK", "n_cells": n_cells, "n_typed": int(len(keep)),
         "n_clusters": nclust,
         "clusters": {str(c): {"n": cl_n[c], "frags": cl_tot[c], "markers": dict(cl_mark[c])}
                      for c in range(nclust)},
         "cells": [[bcs[int(ci)], int(lab[j])] for j, ci in enumerate(keep)]}
    os.makedirs(CL, exist_ok=True)
    json.dump(r, open(dest, "w"))
    return r


def annotate(markers, gsms):
    lin_genes = collections.defaultdict(list)
    lin_bp = collections.Counter()
    for g, v in markers.items():
        lin_genes[v["lineage"]].append(g)
        lin_bp[v["lineage"]] += v["end"] - v["start"]
    lineages = sorted(lin_genes)

    pool, keyed = [], []
    data = {}
    for g in gsms:
        p = os.path.join(CL, f"{g}.json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        data[g] = d
        if d.get("status") != "OK":
            continue
        for c, cd in d["clusters"].items():
            e = []
            for L in lineages:
                raw = sum(cd["markers"].get(x, 0) for x in lin_genes[L])
                exp = cd["frags"] * lin_bp[L] / 3.1e9
                e.append(raw / exp if exp > 0 else 0.0)
            pool.append(e)
            keyed.append((g, c))
    if not pool:
        return {}, {}
    E = np.asarray(pool)
    # log first: enrichment is multiplicative, and a few very pure clusters would
    # otherwise dominate the scale for their lineage
    Z = np.log1p(E)
    Z = (Z - Z.mean(axis=0)) / (Z.std(axis=0) + 1e-9)
    calls = {}
    for i, (g, c) in enumerate(keyed):
        li = int(np.argmax(Z[i]))
        calls[(g, c)] = lineages[li] if Z[i, li] > 0.25 else "ambiguous"
    return calls, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gsm", nargs="*", default=[])
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--annotate-only", action="store_true")
    a = ap.parse_args()

    markers = json.load(open(os.path.join(PC, "meta", "markers.json")))
    allg = sorted({os.path.basename(f).split("_")[0]
                   for f in glob.glob(os.path.join(PC, "dl", "*variant_stats.tsv.gz"))})
    if a.control:
        gsms = ["GSM9096509", "GSM9096510"]
    elif a.all:
        gsms = allg
    else:
        gsms = a.gsm or allg

    os.makedirs(CL, exist_ok=True)
    if not a.annotate_only:
        from concurrent.futures import ThreadPoolExecutor

        def run(g):
            lib = libname(g)
            if not lib:
                return
            try:
                r = cluster_one(g, lib)
                if r:
                    sys.stderr.write(f"[{g}] clustered typed={r.get('n_typed')} "
                                     f"k={r.get('n_clusters')}\n")
            except Exception as e:  # noqa: BLE001
                sys.stderr.write(f"[{g}] ERROR {type(e).__name__}: {e}\n")

        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            list(ex.map(run, gsms))

    # annotate against every cluster available, not just this run's libraries
    pool_gsms = sorted({os.path.basename(f)[:-5] for f in glob.glob(os.path.join(CL, "*.json"))})
    calls, data = annotate(markers, pool_gsms)
    os.makedirs(os.path.join(PC, "calls"), exist_ok=True)
    summary = {}
    for g in pool_gsms:
        d = data.get(g)
        if not d or d.get("status") != "OK":
            continue
        cc = {c: calls.get((g, c), "ambiguous") for c in d["clusters"]}
        cnt = collections.Counter()
        with open(os.path.join(PC, "calls", f"{g}.calls.tsv"), "w") as f:
            for bc, cl in d["cells"]:
                lab = cc[str(cl)]
                cnt[lab] += 1
                f.write(f"{bc}\t{lab}\t{cl}\n")
        summary[g] = {"n_typed": d["n_typed"], "nk": cnt.get("NK", 0), "counts": dict(cnt)}
        if g in gsms:
            print(f"{g} typed={d['n_typed']:6d} NK={cnt.get('NK',0):6d} "
                  f"({cnt.get('NK',0)/max(1,d['n_typed']):5.1%})  " +
                  " ".join(f"{k}={v}" for k, v in cnt.most_common()))
    json.dump(summary, open(os.path.join(PC, "calls", "summary.json"), "w"), indent=1)

    if a.control and all(g in summary for g in ("GSM9096509", "GSM9096510")):
        p = summary["GSM9096509"]["nk"] / max(1, summary["GSM9096509"]["n_typed"])
        n = summary["GSM9096510"]["nk"] / max(1, summary["GSM9096510"]["n_typed"])
        print(f"\nCONTROL CD56+ NK={p:.1%}  CD56- NK={n:.1%}  "
              f"ratio={p/n if n else float('inf'):.1f}x")
        print("PASS" if p > 0.5 and p > 3 * n else "FAIL")


if __name__ == "__main__":
    main()
