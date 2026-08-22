"""Task P-2 §2: verify the human/mouse orthology of the whole panel, gene by gene.

The spec forbids assuming a one-to-one correspondence, because topic 1 was
burned by assuming `XCL2` exists in mouse.  This queries Ensembl Compara for
each gene and records the homology *type*, which is what the assumption
actually turns on.  `XCL1`/`XCL2` are queried as a control on the procedure:
if the method is sound it must reproduce the known failure.

`CTGF` is queried under both its old and current HGNC symbols (`CTGF`,
`CCN2`) because a symbol that no longer resolves would otherwise look like
a missing ortholog rather than a renamed gene.

Run: python3 scripts/p2_orthologs.py     Writes: results/p2_orthologs.csv
"""

from __future__ import annotations

import json
import os
import subprocess
import time

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")

PANEL = {
    "channel": ["PIEZO1", "PIEZO2", "TRPV2", "TRPV4", "TRPM7", "TRPC1",
                "TRPC6", "TMEM63A", "TMEM63B", "KCNK2", "KCNK4"],
    "mechano_downstream": ["YAP1", "WWTR1", "ANKRD1", "CCN2", "CTGF"],
    "state_activation": ["IFNG", "GZMB", "MKI67"],
    "state_residency": ["ITGA1", "CD69", "ITGAE"],
    "state_egress": ["KLF2", "S1PR1", "SELL"],
    "procedure_control": ["XCL1", "XCL2"],
}


def get(url: str, tries: int = 4):
    for i in range(tries):
        r = subprocess.run(["curl", "-sS", "--max-time", "60", url],
                           capture_output=True, text=True).stdout
        try:
            return json.loads(r)
        except Exception:
            time.sleep(3 * (i + 1))
    return None


def main():
    rows = []
    for role, genes in PANEL.items():
        for g in genes:
            j = get(f"https://rest.ensembl.org/homology/symbol/human/{g}"
                    "?target_species=mus_musculus;type=orthologues;"
                    "format=condensed;content-type=application/json")
            if j is None:
                rows.append({"role": role, "human": g, "mouse_ensembl": "",
                             "mouse_symbol": "", "homology_type": "FETCH_FAILED",
                             "one_to_one": False})
                continue
            data = j.get("data") or []
            homs = data[0].get("homologies", []) if data else []
            if not homs:
                rows.append({"role": role, "human": g, "mouse_ensembl": "",
                             "mouse_symbol": "", "homology_type": "NO_MOUSE_ORTHOLOG",
                             "one_to_one": False})
            for h in homs:
                rows.append({"role": role, "human": g, "mouse_ensembl": h.get("id"),
                             "mouse_symbol": "", "homology_type": h.get("type"),
                             "one_to_one": h.get("type") == "ortholog_one2one"})
            time.sleep(1)

    ids = sorted({r["mouse_ensembl"] for r in rows if r["mouse_ensembl"]})
    if ids:
        r = subprocess.run(
            ["curl", "-sS", "--max-time", "90", "-X", "POST",
             "https://rest.ensembl.org/lookup/id",
             "-H", "Content-Type:application/json", "-H", "Accept:application/json",
             "-d", json.dumps({"ids": ids})], capture_output=True, text=True).stdout
        try:
            look = json.loads(r)
        except Exception:
            look = {}
        for row in rows:
            v = look.get(row["mouse_ensembl"]) or {}
            row["mouse_symbol"] = v.get("display_name") or ""

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "p2_orthologs.csv"), index=False)
    print(df.to_string(index=False))
    bad = df[(~df.one_to_one) & (df.role != "procedure_control")]
    print(f"\nnot one-to-one outside the control: {bad.human.tolist() or 'none'}")
    ctrl = df[df.role == "procedure_control"]
    print("procedure control (XCL1/XCL2) mouse symbols: "
          f"{sorted(set(ctrl.mouse_symbol))}, one_to_one={ctrl.one_to_one.tolist()}")


if __name__ == "__main__":
    main()
