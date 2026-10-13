"""Corrected reference integration (replaces the broken zero_check).

What was wrong before: the old script never opened the reference files and never used the candidate's own
CASP1 value — the verdict was a pure function of the lineage label. This script reads the unified reference
table and the candidates' own numbers, and reports a genuine side-by-side.

Per spec §9.5 this is PRESENTATION, not adjudication: no statistical test is run, and no candidate is
deleted. The evidence tier is an annotation that the reader can re-cut; the thresholds below are declared
engineering choices with no biological authority, exactly like the terciles in §8.

Ratio rule: IL18:CASP1 is computed ONLY within one source + one unit + one state, never across sources
(units are not comparable). Values at a source's zero floor are treated as "near-absent".
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import write_tsv

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ext" / "unified_reference.tsv"

# declared engineering thresholds (no biological authority)
RATIO_STRONG = 20.0   # IL18 >= 20x CASP1 within one source/unit -> CASP1 "near-absent" relative to IL18
RATIO_WEAK = 3.0      # IL18 >= 3x CASP1 -> CASP1 "low relative to IL18"
MIN_SOURCES_FOR_A = 2  # tier A needs >= 2 independent sources agreeing

# map each screen cell type to the reference table's cell_type_screen label (explicit, no regex matching)
SCREEN_TO_REF = {
    ("human", "liver", "basophil"): "liver basophil",
    ("human", "lung", "mast cell"): "lung mast cell",
    ("human", "lung", "granulocyte"): "lung mast cell",   # this lineage row is the mast cells re-aggregated
    ("human", "lung", "dendritic cell"): "lung pDC",
    ("human", "lung", "multiciliated columnar cell of tracheobronchial tree"): "lung multiciliated/ciliated epithelial cell",
    ("human", "lung", "ciliated cell"): "lung multiciliated/ciliated epithelial cell",
    ("mouse", "lung", "alveolar macrophage"): "lung alveolar macrophage",
    ("mouse", "lung", "resident_macrophage"): "lung alveolar macrophage",
    ("mouse", "brain", "astrocyte"): "brain astrocyte",
    ("mouse", "brain", "oligodendrocyte"): "brain oligodendrocyte",
    ("mouse", "brain", "glia"): "brain astrocyte",        # lineage row = astrocyte + oligodendrocyte
    ("mouse", "brain", "Cajal-Retzius cell"): "brain Cajal-Retzius neuron",
    ("mouse", "brain", "central nervous system macrophage"): "brain microglia",
    ("mouse", "brain", "resident_macrophage"): "brain microglia",
    ("mouse", "pancreas", "pancreatic D cell"): "islet delta cell",
    ("mouse", "pancreas", "pancreatic ductal cell"): "pancreas ductal cell",
}
# which candidate rows are re-aggregations of another row (so counts are not inflated)
DUPLICATE_OF = {
    ("human", "lung", "granulocyte"): "mast cell",
    ("mouse", "lung", "resident_macrophage"): "alveolar macrophage",
    ("mouse", "brain", "resident_macrophage"): "central nervous system macrophage",
    ("mouse", "brain", "glia"): "astrocyte + oligodendrocyte",
}


def norm_source(s):
    """Collapse variants of the same underlying resource so independence is not over-counted."""
    s = str(s)
    low = s.lower()
    if "protein atlas" in low or low.startswith("hpa"):
        return "HumanProteinAtlas"
    if "immgen" in low:
        return "ImmGen"
    if "cellxgene" in low or "cz_" in low:
        return "CELLxGENE"
    if "brainrnaseq" in low or "zhang" in low:
        return "brainrnaseq/Zhang"
    if "allen" in low:
        return "AllenBrainMap"
    if "gtex" in low:
        return "GTEx"
    return re.split(r"[ (\[]", s)[0][:40]


def numeric(v):
    try:
        f = float(v)
        return f if np.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def zero_floor(unit):
    """Value at or below which a source effectively reports zero."""
    u = str(unit).lower()
    if "immgen" in u:
        return 1.0          # ImmGen: 1.00 == zero reads
    if "fpkm" in u or "rpkm" in u:
        return 0.1          # Barres pipeline floors at 0.1
    return 0.0


def ratios_for(ref_rows):
    """Per (source, unit, state): median IL18 and CASP1, and how CASP1 sits relative to the source's zero
    floor. Groups where IL18 itself is at the floor are dropped as uninformative (they cannot support
    'IL18 present'). Exact matches are preferred; close matches are used only if no exact row exists."""
    if (ref_rows.match_quality == "exact").any():
        ref_rows = ref_rows[ref_rows.match_quality == "exact"]
    out = []
    for key, g in ref_rows.groupby(["source", "unit", "state"], dropna=False):
        iv = [x for x in (numeric(v) for v in g[g.gene.str.upper() == "IL18"].value) if x is not None]
        cv = [x for x in (numeric(v) for v in g[g.gene.str.upper() == "CASP1"].value) if x is not None]
        if not iv or not cv:
            continue
        i_med, c_med = float(np.median(iv)), float(np.median(cv))
        floor = zero_floor(key[1])
        if i_med <= floor:          # IL18 itself undetectable here -> this source cannot inform us
            continue
        casp1_absent = c_med <= floor
        c_eff = max(c_med - floor, 0.0)
        ratio = (i_med - floor) / c_eff if c_eff > 0 else np.inf
        out.append({"source": norm_source(key[0]), "unit": key[1], "state": key[2], "il18": i_med, "casp1": c_med,
                    "ratio": ratio, "casp1_absent": casp1_absent,
                    "match": ",".join(sorted(set(g.match_quality)))})
    return out


def tier(rs, has_ref):
    """A = CASP1 at/near the zero floor (or >=RATIO_STRONG below IL18) in >=2 independent sources, with no
    source showing CASP1 at the same order as IL18. B = the same in exactly one source, i.e. a genuine
    possibility that one reference supports. C = CASP1 clearly present alongside IL18. D = no matched
    reference. Thresholds are declared engineering choices, not biology."""
    if not has_ref:
        return "E_no_matched_reference", "参考层中没有对应群体，无法判读"
    steady = [r for r in rs if not any(k in str(r["state"]).lower()
                                       for k in ("challenge", "disease", "tumour", "tumor", "psoria"))]
    use = steady or rs
    if not use:
        return "E_no_matched_reference", "无可用的静息态参考值（IL18 在该源亦不可测）"
    good = [r for r in use if r["casp1_absent"] or r["ratio"] >= RATIO_STRONG]
    contra = [r for r in use if r["ratio"] < RATIO_WEAK]
    n_src = len({r["source"] for r in good})
    if n_src >= MIN_SOURCES_FOR_A and not contra:
        return "A_supported", f"{n_src} 个独立来源显示 CASP1 处于零底线或比 IL18 低 >={RATIO_STRONG:.0f} 倍，无反向来源"
    if n_src >= 1:
        extra = f"，但另有 {len(contra)} 个来源显示 CASP1 与 IL18 同量级" if contra else "（单一来源）"
        return "B_potential", f"{n_src} 个来源支持 CASP1 near-absent/低{extra}"
    mid = [r for r in use if RATIO_WEAK <= r["ratio"] < RATIO_STRONG]
    if mid:
        lo, hi = min(r["ratio"] for r in mid), max(r["ratio"] for r in mid)
        return "C_casp1_present", (f"CASP1 可测到，比 IL18 低 {lo:.3g}–{hi:.3g} 倍：IL18 仍是主导转录本，"
                                   f"但不满足“加工酶缺失”前提")
    inv = [r for r in use if r["ratio"] < 1.0]
    if inv:
        return "D_refuted", f"CASP1 >= IL18（最低比值 {min(r['ratio'] for r in inv):.2g}），与筛选前提方向相反"
    return "D_refuted", (f"CASP1 与 IL18 同量级（比值 {min(r['ratio'] for r in use):.2g}–"
                         f"{max(r['ratio'] for r in use):.2g}）")


def main():
    cand = pd.read_csv(ROOT / "results" / "candidates.tsv", sep="\t", comment="#")
    ref = pd.read_csv(REF, sep="\t")
    ref = ref[ref.match_quality.isin(["exact", "close"])]
    rows = []
    for _, c in cand.iterrows():
        key = (c.species, c.tissue, str(c.celltype))
        ref_label = SCREEN_TO_REF.get(key)
        sub = ref[(ref.species == c.species) & (ref.cell_type_screen == ref_label)] if ref_label else ref.iloc[0:0]
        rs = ratios_for(sub) if len(sub) else []
        t, why = tier(rs, bool(len(sub)))
        best = max(rs, key=lambda r: r["ratio"]) if rs else None
        rows.append({
            "dataset_id": c.dataset_id, "level": c.level, "species": c.species, "tissue": c.tissue,
            "suspension_type": c.suspension_type, "celltype": c.celltype,
            "duplicate_of": DUPLICATE_OF.get(key, ""),
            "sc_pct_IL18": c.pct_IL18, "sc_pct_CASP1": c.pct_CASP1,
            "sc_logcpm_IL18": c.logcpm_IL18, "sc_logcpm_CASP1": c.logcpm_CASP1,
            "sc_CASP1_detected": bool(c.pct_CASP1 > 0),
            "ref_label": ref_label or "", "n_ref_sources": len(rs),
            "ref_sources": "; ".join(f"{r['source'][:34]}[{r['unit'][:14]}] IL18={r['il18']:.3g}/CASP1={r['casp1']:.3g} ratio={r['ratio']:.3g}" for r in rs),
            "best_ratio_IL18_over_CASP1": (f"{best['ratio']:.3g}" if best else ""),
            "evidence_tier": t, "tier_reason": why,
        })
    out = pd.DataFrame(rows).sort_values(["evidence_tier", "species", "tissue"])
    write_tsv(out, ROOT / "results" / "il18_casp1_evidence.tsv")
    print(out.evidence_tier.value_counts().to_string())
    print("\n--- distinct populations (duplicates excluded) ---")
    d = out[out.duplicate_of == ""]
    for _, r in d.iterrows():
        print(f"  {r.evidence_tier:24} {r.species:6} {r.tissue:9} {str(r.celltype)[:38]:38} "
              f"IL18检出={r.sc_pct_IL18:.2f} CASP1检出={r.sc_pct_CASP1:.3f} ratio={r.best_ratio_IL18_over_CASP1}")


if __name__ == "__main__":
    main()
