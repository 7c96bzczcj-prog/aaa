#!/usr/bin/env python3
"""GEO admission scanner for the BM-PB-NK preflight.

For every candidate GSE it fetches the series `filelist.txt` (small) and
classifies whether an UNFILTERED droplet matrix is deposited -- admission
criterion 1, which is a hard gate for task A.

Classification is deliberately conservative: a series only counts as
carrying raw droplets if a file NAME says so (raw_feature_bc_matrix /
raw_gene_bc_matrices / *raw*.h5) or a barcodes file is large enough that
it cannot be a called-cell list (>500 kB gzipped ~ >100k barcodes).
Everything else is recorded as filtered/other, with the evidence kept.
"""
import json, os, re, sys, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/scan"

# gzipped barcodes bigger than this cannot be a called-cell list
RAW_BARCODE_GZ_BYTES = 500_000
RAW_NAME_RE = re.compile(r"raw[_\-.]?(feature|gene)?[_\-]?bc|unfiltered|raw_matri|_raw\.h5|raw_feature", re.I)
BARCODE_RE = re.compile(r"barcode", re.I)


def get(url, tries=4):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if a == tries - 1:
                return ""
            time.sleep(2 ** a)
    return ""


def esearch(term, retmax=500):
    u = f"{EUTILS}/esearch.fcgi?db=gds&retmax={retmax}&term={urllib.parse.quote(term)}"
    return re.findall(r"<Id>(\d+)</Id>", get(u))


def esummary(ids):
    out = {}
    for i in range(0, len(ids), 100):
        chunk = ",".join(ids[i:i + 100])
        txt = get(f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&id={chunk}")
        try:
            res = json.loads(txt).get("result", {})
        except Exception:
            print(f"[warn] esummary chunk {i} unparseable ({len(txt)}B)", file=sys.stderr)
            time.sleep(1.0)
            continue
        for uid in res.get("uids", []):
            f = res.get(uid, {})
            acc = f.get("accession", "")
            if acc.startswith("GSE"):
                out[acc] = {"title": f.get("title", ""),
                            "n_samples": f.get("n_samples", ""),
                            "taxon": f.get("taxon", ""),
                            "gdsType": f.get("gdstype", f.get("gdsType", "")),
                            "pdat": f.get("pdat", ""),
                            "summary": f.get("summary", "")[:400]}
        time.sleep(0.5)
    return out


def filelist(gse):
    num = gse[3:]
    pre = "GSE" + (num[:-3] if len(num) > 3 else "") + "nnn"
    return get(f"{FTP}/{pre}/{gse}/suppl/filelist.txt")


def classify(gse):
    txt = filelist(gse)
    if not txt:
        return {"accession": gse, "filelist_ok": False, "raw_droplets": "UNKNOWN",
                "evidence": "filelist.txt unreachable", "files": 0}
    rows = []
    for line in txt.splitlines()[1:]:
        p = line.split("\t")
        if len(p) >= 5 and p[0] == "File":
            try:
                sz = int(p[3])
            except ValueError:
                sz = 0
            rows.append((p[1], sz, p[4]))
    named = [n for n, s, t in rows if RAW_NAME_RE.search(n)]
    bigbc = [(n, s) for n, s, t in rows if BARCODE_RE.search(n) and s > RAW_BARCODE_GZ_BYTES]
    smallbc = [(n, s) for n, s, t in rows if BARCODE_RE.search(n) and s <= RAW_BARCODE_GZ_BYTES]
    if named:
        verdict, ev = "YES_BY_NAME", "; ".join(named[:3])
    elif bigbc:
        verdict, ev = "YES_BY_SIZE", "; ".join(f"{n}={s}B" for n, s in bigbc[:3])
    elif smallbc:
        verdict, ev = "NO_FILTERED", "largest barcodes file %s" % max(
            (f"{n}={s}B" for n, s in smallbc), key=lambda x: int(x.split("=")[-1][:-1]))
    else:
        kinds = sorted({t for _, _, t in rows})
        verdict, ev = "NO_BARCODE_FILES", "file types: " + ",".join(kinds[:6])
    return {"accession": gse, "filelist_ok": True, "raw_droplets": verdict,
            "evidence": ev, "files": len(rows),
            "example_files": [n for n, _, _ in rows[:6]]}


def main():
    queries = json.load(open(sys.argv[1]))
    ids = []
    per_query = {}
    for name, term in queries.items():
        got = esearch(term)
        per_query[name] = len(got)
        ids += got
        time.sleep(0.4)
    ids = sorted(set(ids))
    print(f"[esearch] {per_query} -> {len(ids)} unique uids", file=sys.stderr)
    summ = esummary(ids)
    gses = sorted(summ)
    print(f"[esummary] {len(gses)} GSE", file=sys.stderr)
    with ThreadPoolExecutor(max_workers=4) as ex:
        cls = list(ex.map(classify, gses))
    for c in cls:
        c.update(summ.get(c["accession"], {}))
    os.makedirs(OUT, exist_ok=True)
    tag = sys.argv[2] if len(sys.argv) > 2 else "scan"
    with open(f"{OUT}/{tag}.json", "w") as fh:
        json.dump(cls, fh, indent=1)
    n_raw = sum(1 for c in cls if c["raw_droplets"].startswith("YES"))
    print(f"[done] {len(cls)} series, {n_raw} with raw droplet matrices -> {OUT}/{tag}.json",
          file=sys.stderr)


if __name__ == "__main__":
    main()
