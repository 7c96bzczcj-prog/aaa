#!/usr/bin/env python3
"""Download the HCA Census-of-Immune-Cells raw droplet matrices we admit:
MantonBM1-8 (bone marrow) and MantonBL1-8 (adult peripheral blood).
Cord blood (MantonCB*) and the pooled libraries (*P) are not downloaded.
"""
import json, os, sys, time, urllib.parse, urllib.request, re
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = "/home/user/aaa/data/hca_ica"
B = "https://service.azul.data.humancellatlas.org/index/files"
PROJ = "cc95ff89-2e68-4a08-a234-480eca21ce79"
KEEP = re.compile(r"^Manton(BM|BL)[1-8]_")


def listing():
    filt = json.dumps({"projectId": {"is": [PROJ]}, "fileFormat": {"is": ["h5"]}})
    url = f"{B}?size=75&filters={urllib.parse.quote(filt)}"
    rows = []
    while url:
        with urllib.request.urlopen(url, timeout=120) as r:
            d = json.load(r)
        for h in d["hits"]:
            f = h["files"][0]
            dn = (h.get("donorOrganisms") or [{}])[0]
            cs = (h.get("cellSuspensions") or [{}])[0]
            pr = h.get("protocols") or []
            chem = next((p["libraryConstructionApproach"][0] for p in pr
                         if p.get("libraryConstructionApproach")), None)
            rows.append({"name": f["name"], "size": f["size"], "url": f["azul_url"],
                         "sha256": f["sha256"],
                         "organ": (h.get("specimens") or [{}])[0].get("organ"),
                         "cell_type": cs.get("selectedCellType"),
                         "target_cells": cs.get("totalCells"),
                         "age": dn.get("organismAge"), "sex": dn.get("biologicalSex"),
                         "chemistry": chem})
        url = d["pagination"].get("next")
        time.sleep(0.2)
    return rows


def fetch(row):
    out = os.path.join(DEST, row["name"])
    if os.path.exists(out) and os.path.getsize(out) == row["size"]:
        return "cached"
    for a in range(5):
        try:
            urllib.request.urlretrieve(row["url"], out + ".part")
            if os.path.getsize(out + ".part") == row["size"]:
                os.rename(out + ".part", out)
                return "ok"
        except Exception:
            pass
        time.sleep(2 ** a)
    return "FAIL " + row["name"]


if __name__ == "__main__":
    os.makedirs(DEST, exist_ok=True)
    rows = listing()
    json.dump(rows, open(os.path.join(HERE, "scan", "hca_ica_files.json"), "w"), indent=1)
    keep = [r for r in rows if KEEP.match(r["name"])]
    print(f"{len(rows)} h5 total, {len(keep)} admitted (MantonBM1-8 + MantonBL1-8), "
          f"{sum(r['size'] for r in keep)/1e9:.2f} GB", flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = list(ex.map(fetch, keep))
    bad = [r for r in res if r.startswith("FAIL")]
    print(f"downloaded ok={res.count('ok')} cached={res.count('cached')} failed={len(bad)}")
    for b in bad:
        print(b)
