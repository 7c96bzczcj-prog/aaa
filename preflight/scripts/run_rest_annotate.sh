#!/bin/bash
set -u
PY=/home/user/aaa/.venv/bin/python
cd /home/user/aaa/preflight
while ! grep -q ALL_DONE scan/annotate_gse233304.log 2>/dev/null; do sleep 20; done
$PY - <<'PYX'
import json, subprocess, collections, os
PY="/home/user/aaa/.venv/bin/python"
for tag, idx in [("GSE181543","index/gse181543.json"), ("GSE120221","index/gse120221.json")]:
    ent=json.load(open(idx))
    by=collections.defaultdict(list)
    for e in ent: by[e["donor"]].append(e)
    for donor, es in sorted(by.items()):
        out=f"labels/{tag}_{donor}.npz"
        spec=[[e["library"], e["path"], e["kind"]] for e in es]
        r=subprocess.run([PY,"scripts/annotate.py","--donor",f"{tag}_{donor}",
                          "--files",json.dumps(spec),"--out",out])
        if r.returncode: print("FAIL",tag,donor,flush=True)
PYX
echo ALL_DONE
