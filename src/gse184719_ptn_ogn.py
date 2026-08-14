"""Check 2 and check 3, on the one dataset that can answer both.

GSE184719 is the bulk RNA-seq behind Du et al., Front Immunol 2022
("Human-Induced CD49a+ NK Cells Promote Fetal Growth") — the same group as
Fu et al., Immunity 2017.  20 libraries:

  pNK1-4        freshly sorted peripheral blood NK
  dNK1-4        freshly sorted decidual NK          <- check 2 asks about these
  CB/BM/pNK-iNK 12 in vitro induced NK cultures     <- check 3 point 4

Two questions, one table:

**Check 2 (contamination).** Sorting purity is measured on surface markers and
cannot exclude stromal *transcripts* riding along on carried-over cells or
free RNA.  The decisive quantity is whether stromal genes an NK cell cannot
transcribe — COL1A1, COL1A2, COL3A1, DCN, LUM — appear in the sorted dNK
libraries, and at what level relative to PTN/OGN.

**Check 3 point 4 (the in vitro arm).** The iNK cultures are feeder-free and
contain no stroma.  PTN/OGN mRNA there cannot be explained by contamination.
Absence there is equally informative: it means the paper's "high expression of
growth-promoting factors" does not rest on these two genes.

VIM is reported but is *not* a contamination control: leukocytes transcribe
vimentin.  PTPRC (CD45) is the reciprocal control — it should be high in every
library, and near zero in a stromal cell.

No scanpy (R10).  Gene ids are Ensembl with a Symbol column (R11).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

COUNTS = Path("/tmp/claude-0/-home-user-aaa/45313b32-c7b9-5c12-8800-46f354fcf387"
              "/scratchpad/gse184719_counts.txt.gz")
OUT = Path("out/GSE184719")

TARGETS = ["PTN", "OGN", "SPP1"]
STROMAL = ["COL1A1", "COL1A2", "COL3A1", "COL6A2", "DCN", "LUM", "PDGFRB", "ACTA2"]
DECIDUA_OTHER = ["PAEP", "KRT7", "HLA-G", "EPCAM"]
NOT_A_CONTROL = ["VIM"]
NK_IDENTITY = ["PTPRC", "NCAM1", "KLRD1", "NKG7", "GNLY", "ITGA1", "CD3E"]

GROUPS = {
    "pNK": ["pNK1", "pNK2", "pNK3", "pNK4"],
    "dNK": ["dNK1", "dNK2", "dNK3", "dNK4"],
    "CB_iNK": ["CB_iNK1", "CB_iNK2", "CB_iNK3", "CB_iNK4"],
    "BM_iNK": ["BM_iNK1", "BM_iNK2", "BM_iNK3", "BM_iNK4"],
    "pNK_iNK": ["pNK_iNK1", "pNK_iNK2", "pNK_iNK3", "pNK_iNK4"],
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(COUNTS, sep="\t")
    raw = raw.rename(columns={c: c.replace("_count", "") for c in raw.columns})
    samples = [s for g in GROUPS.values() for s in g]
    missing_s = [s for s in samples if s not in raw.columns]
    if missing_s:
        raise SystemExit(f"missing sample columns: {missing_s}")

    counts = raw[samples].astype(np.int64)
    lib = counts.sum(axis=0)
    cpm = counts / lib * 1e6
    cpm["Symbol"] = raw["Symbol"].astype(str)

    panel = TARGETS + STROMAL + DECIDUA_OTHER + NOT_A_CONTROL + NK_IDENTITY
    present = sorted(set(cpm["Symbol"]) & set(panel))
    rate = len(present) / len(panel)
    print(f"panel {len(present)}/{len(panel)} matched ({rate:.1%}); "
          f"missing {[g for g in panel if g not in present]}")
    if rate < 0.90:
        raise SystemExit(f"panel match rate {rate:.1%} < 90% (R11)")

    # a symbol can map to several Ensembl ids; sum them, which is what a
    # symbol-level statement means
    tab = cpm[cpm["Symbol"].isin(panel)].groupby("Symbol")[samples].sum()
    tab = tab.reindex([g for g in panel if g in tab.index])

    grp = pd.DataFrame(
        {name: tab[cols].mean(axis=1) for name, cols in GROUPS.items()}
    )
    grp["iNK_all_mean"] = tab[GROUPS["CB_iNK"] + GROUPS["BM_iNK"]
                              + GROUPS["pNK_iNK"]].mean(axis=1)
    grp["iNK_all_max"] = tab[GROUPS["CB_iNK"] + GROUPS["BM_iNK"]
                             + GROUPS["pNK_iNK"]].max(axis=1)

    print("\nlibrary size (M reads):")
    print((lib / 1e6).round(2).to_string())

    print("\n=== CPM by group (mean of 4 libraries) ===")
    print(grp.round(2).to_string())

    print("\n=== per-library CPM, the genes the two checks turn on ===")
    key = [g for g in ["PTN", "OGN", "SPP1", "COL1A1", "COL1A2", "COL3A1",
                       "DCN", "LUM", "PAEP", "VIM", "PTPRC"] if g in tab.index]
    print(tab.loc[key].round(2).to_string())

    print("\n=== raw counts, same genes (a CPM is not evidence if n reads = 0) ===")
    rawtab = raw[raw["Symbol"].isin(key)].groupby("Symbol")[samples].sum()
    print(rawtab.reindex(key).to_string())

    tab.round(4).to_csv(OUT / "panel_cpm_per_library.csv")
    grp.round(4).to_csv(OUT / "panel_cpm_by_group.csv")
    rawtab.reindex(key).to_csv(OUT / "panel_rawcounts_per_library.csv")

    # --- the two readings, stated as numbers rather than as prose ----------
    print("\n=== reading ===")
    for g in ["PTN", "OGN"]:
        if g not in tab.index:
            continue
        d = grp.loc[g, "dNK"]
        i_mean, i_max = grp.loc[g, "iNK_all_mean"], grp.loc[g, "iNK_all_max"]
        p = grp.loc[g, "pNK"]
        print(f"{g}: dNK {d:.2f} CPM | pNK {p:.2f} | iNK mean {i_mean:.2f}, "
              f"max {i_max:.2f} (12 stroma-free libraries)")
    col = [g for g in ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"] if g in tab.index]
    print(f"\nstromal sum in dNK: {grp.loc[col, 'dNK'].sum():.2f} CPM "
          f"(pNK {grp.loc[col, 'pNK'].sum():.2f}, "
          f"iNK {grp.loc[col, 'iNK_all_mean'].sum():.2f})")
    print("per-library stromal sum in dNK:")
    print(tab.loc[col, GROUPS['dNK']].sum().round(2).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
