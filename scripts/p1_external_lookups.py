"""Task P-1: the two external lookups the task requires, made reproducible.

Neither of these is a dataset download.  They are the two checks §2 asks
for by name, and both are small, addressable queries rather than files:

  1. Ortholog verification (§1 species note).  The task forbids assuming a
     one-to-one human/mouse correspondence, because topic 1 was burned by
     assuming `XCL2` exists in mouse.  `Ensembl REST /homology/symbol`
     answers this per gene and reports the homology *type*, which is what
     the assumption actually turns on.  `XCL1`/`XCL2` are queried as a
     control on the procedure itself: if the method is sound it must
     reproduce the known failure.

  2. An orthogonal transcript reference for Q1/Q2.  The Human Protein
     Atlas consensus blood-cell RNA-seq measures sorted NK cells in bulk,
     without droplets, without 3' capture, and without dropout.  It cannot
     confirm a single-cell detection rate -- different quantity -- but it
     can say whether the channel panel's ordering survives a change of
     assay, which is exactly what a floor-band ranking needs.

Run: python3 scripts/p1_external_lookups.py
Writes: results/p1_orthologs.csv, results/p1_hpa_blood_rna.csv
"""

from __future__ import annotations

import json
import os
import subprocess
import time
import urllib.parse

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")

PANEL = ["PIEZO1", "PIEZO2", "TRPV4", "TRPV2", "TRPM7", "TRPC1",
         "TMEM63A", "TMEM63B", "KCNK2", "KCNK4"]
Q6 = ["ITGA1", "ITGAE", "CD69", "KLF2", "S1PR1", "S1PR5", "SELL", "KLRD1", "NKG7"]
PROCEDURE_CONTROL = ["XCL1", "XCL2"]


def get(url: str, tries: int = 4):
    for i in range(tries):
        r = subprocess.run(["curl", "-sS", "--max-time", "60", url],
                           capture_output=True, text=True).stdout
        try:
            return json.loads(r)
        except Exception:
            time.sleep(3 * (i + 1))
    return None


def orthologs(genes) -> pd.DataFrame:
    rows = []
    for g in genes:
        j = get(f"https://rest.ensembl.org/homology/symbol/human/{g}"
                "?target_species=mus_musculus;type=orthologues;format=condensed;"
                "content-type=application/json")
        if j is None:
            rows.append({"human": g, "mouse_ensembl": "", "mouse_symbol": "",
                         "homology_type": "FETCH_FAILED", "one_to_one": False})
            continue
        data = j.get("data") or []
        homs = data[0].get("homologies", []) if data else []
        if not homs:
            rows.append({"human": g, "mouse_ensembl": "", "mouse_symbol": "",
                         "homology_type": "NO_MOUSE_ORTHOLOG", "one_to_one": False})
        for h in homs:
            rows.append({"human": g, "mouse_ensembl": h.get("id"),
                         "mouse_symbol": "", "homology_type": h.get("type"),
                         "one_to_one": h.get("type") == "ortholog_one2one"})
        time.sleep(1)

    ids = sorted({r["mouse_ensembl"] for r in rows if r["mouse_ensembl"]})
    if ids:
        r = subprocess.run(
            ["curl", "-sS", "--max-time", "90", "-X", "POST",
             "https://rest.ensembl.org/lookup/id", "-H", "Content-Type:application/json",
             "-H", "Accept:application/json", "-d", json.dumps({"ids": ids})],
            capture_output=True, text=True).stdout
        try:
            look = json.loads(r)
        except Exception:
            look = {}
        for row in rows:
            v = look.get(row["mouse_ensembl"]) or {}
            row["mouse_symbol"] = v.get("display_name") or ""
    return pd.DataFrame(rows)


def hpa_blood(genes, cell_types=("NK-cell", "T-reg", "neutrophil", "basophil")) -> pd.DataFrame:
    """HPA consensus blood-cell nTPM.

    Only single-token cell-type names resolve through this endpoint; the
    multi-word ones ("classical monocyte") return the gene with no value,
    so they are left out rather than reported as zero.
    """
    rows = []
    for g in genes:
        rec = {"gene": g}
        for ct in cell_types:
            cols = urllib.parse.quote(f"g,blood_RNA_{ct}", safe=",")
            j = get(f"https://www.proteinatlas.org/api/search_download.php"
                    f"?search={urllib.parse.quote(g)}&format=json&columns={cols}"
                    "&compress=no")
            val = None
            if j:
                hit = [x for x in j if x.get("Gene") == g]
                if hit:
                    val = hit[0].get(f"Blood RNA - {ct} [nTPM]")
            rec[f"nTPM_{ct}"] = float(val) if val not in (None, "") else None
        j = get("https://www.proteinatlas.org/api/search_download.php"
                f"?search={urllib.parse.quote(g)}&format=json&columns=g,rnabcs&compress=no")
        hit = [x for x in (j or []) if x.get("Gene") == g]
        rec["blood_cell_specificity"] = (hit[0].get("RNA blood cell specificity")
                                         if hit else "NOT_RETURNED")
        rows.append(rec)
        print(rec, flush=True)
    return pd.DataFrame(rows)


def main():
    orth = orthologs(PANEL + Q6 + PROCEDURE_CONTROL)
    orth.to_csv(os.path.join(OUT, "p1_orthologs.csv"), index=False)
    bad = orth[~orth.one_to_one]
    print("\n=== orthology ===")
    print(orth.to_string(index=False))
    print(f"not one-to-one: {bad.human.tolist() or 'none'}")

    hpa = hpa_blood(PANEL + Q6)
    hpa.to_csv(os.path.join(OUT, "p1_hpa_blood_rna.csv"), index=False)
    print("\n=== HPA consensus blood RNA (nTPM) ===")
    print(hpa.to_string(index=False))
    p = hpa[hpa.gene.isin(PANEL)].sort_values("nTPM_NK-cell", ascending=False)
    print("\npanel order in sorted blood NK cells: "
          + " > ".join(f"{r.gene}({r._asdict()['nTPM_NK-cell']})"
                       for r in p.itertuples()))


if __name__ == "__main__":
    main()
