#!/bin/bash
set -u
PY=/home/user/aaa/.venv/bin/python
cd /home/user/aaa/preflight
$PY - <<'PYX'
import json, subprocess, collections, os
PY="/home/user/aaa/.venv/bin/python"
jobs=[("", "index/hca.json", lambda d: f"labels/{d}.npz"),
      ("GSE233304_", "index/gse233304.json", lambda d: f"labels/GSE233304_{d}.npz"),
      ("GSE181543_", "index/gse181543.json", lambda d: f"labels/GSE181543_{d}.npz"),
      ("GSE120221_", "index/gse120221.json", lambda d: f"labels/GSE120221_{d}.npz")]
for pre, idx, outf in jobs:
    ent=json.load(open(idx))
    by=collections.defaultdict(list)
    for e in ent: by[e["donor"]].append(e)
    for donor, es in sorted(by.items()):
        spec=[[e["library"], e["path"], e["kind"]] for e in es]
        r=subprocess.run([PY,"scripts/annotate.py","--donor",pre+donor,
                          "--files",json.dumps(spec),"--out",outf(donor)])
        if r.returncode: print("FAIL",pre,donor,flush=True)
print("ALL_DONE")
PYX
