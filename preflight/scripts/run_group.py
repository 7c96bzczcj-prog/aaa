#!/usr/bin/env python3
"""Annotate every donor of the named index files, skipping outputs that exist."""
import collections, json, os, subprocess, sys

PY = "/home/user/aaa/.venv/bin/python"
JOBS = {
    "hca": [("", "index/hca.json")],
    "rest": [("GSE233304_", "index/gse233304.json"),
             ("GSE181543_", "index/gse181543.json"),
             ("GSE120221_", "index/gse120221.json")],
}
for pre, idx in JOBS[sys.argv[1]]:
    ent = json.load(open(idx))
    by = collections.defaultdict(list)
    for e in ent:
        by[e["donor"]].append(e)
    for donor, es in sorted(by.items()):
        out = f"labels/{pre}{donor}.npz"
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"[skip] {pre}{donor}", flush=True)
            continue
        spec = [[e["library"], e["path"], e["kind"]] for e in es]
        r = subprocess.run([PY, "scripts/annotate.py", "--donor", pre + donor,
                            "--files", json.dumps(spec), "--out", out])
        if r.returncode:
            print("FAIL", pre, donor, flush=True)
print("ALL_DONE")
