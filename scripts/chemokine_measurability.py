"""Can this dataset measure NK chemokine axes at all, and across how many compartments?

The 四格挖掘协议 record ends on a falsification: three NK-specific
technical effects make the four-quadrant framework unreadable on
GSE154826.  A chemokine cross-compartment protocol inherits that
dataset, that soup, and those confounders, so the first question is not
"which axis is interesting" but "which axis is measurable, and between
which compartments".  This script answers both from data already on
disk, before any hypothesis is committed to.

Three parts:

  1. COMPARTMENT INVENTORY.  How many individuals carry each combination
     of tissues in GSE154826.  A cross-compartment contrast is a paired
     design, so its n is the number of individuals holding *both*
     compartments, not the number of libraries.

  2. MEASURABILITY.  For the full chemokine system (receptors, ligands,
     the retention/egress module) plus the same positive and negative
     controls used by the complosome premise check: detection against
     the protocol's own floor, mean CPM, and the ambient soup fraction,
     per lineage.  A receptor that fails the detection floor cannot be
     tested by any downstream phase -- that is D10's lesson, where the
     protocol's positive control turned out to be floor-saturated.

  3. EFFECT SIZES, DESCRIPTIVE ONLY.  For genes that clear both gates,
     the already-computed paired tumour-vs-normal statistics are joined
     in.  They are printed as a feasibility readout, not as findings:
     the quadrant reading of this dataset is retracted (README, "The
     result: the framework's core assumption is falsified").

The ambient scales are calibrated inside NK from genes of known truth,
exactly as in scripts/complosome_premise.py, because a raw soup fraction
cannot be read at face value in this dataset.
"""

from __future__ import annotations

import os
import sys
import collections
import csv

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
META = os.path.join(ROOT, "phase1_registry", "raw_meta",
                    "GSE154826_sample_annots.csv")

# ---- the panel, pre-registered here rather than chosen after looking --
# Grouped by the role each gene would play in a cross-compartment
# argument, so that a group failing wholesale is visible as such.
PANEL = {
    "receptor_inflammatory": ["CCR1", "CCR2", "CCR3", "CCR5", "CXCR1", "CXCR2",
                              "CXCR3", "CXCR6", "CX3CR1"],
    "receptor_homing": ["CCR4", "CCR6", "CCR7", "CCR8", "CCR9", "CCR10",
                        "CXCR4", "CXCR5", "XCR1"],
    "receptor_atypical": ["ACKR1", "ACKR2", "ACKR3", "ACKR4", "CCRL2"],
    "retention_egress": ["S1PR1", "S1PR5", "SELL", "SELPLG", "CD69",
                         "ITGAE", "ITGA1", "KLRG1"],
    "ligand_lymphoid": ["CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5",
                        "XCL1", "XCL2"],
    "ligand_myeloid_stromal": ["CCL2", "CCL7", "CCL8", "CCL13", "CCL14",
                               "CCL15", "CCL17", "CCL18", "CCL19", "CCL20",
                               "CCL21", "CCL22", "CCL23", "CCL24", "CCL28",
                               "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6",
                               "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL12",
                               "CXCL13", "CXCL14", "CXCL16", "CXCL17",
                               "CX3CL1"],
}
POS_CTRL = ["KLRD1", "NKG7", "GNLY"]      # genuinely NK
NEG_CTRL = ["C1QA", "LYZ", "IGKC"]        # myeloid / B, known to leak

MIN_COUNT, MIN_FRAC = 10, 0.5             # the protocol's detection floor


def compartment_inventory(log):
    """Part 1: how many individuals hold each combination of tissues."""
    rows = [r for r in csv.DictReader(open(META)) if r["Species"] == "human"]
    by_patient = collections.defaultdict(set)
    prep = collections.defaultdict(set)
    for r in rows:
        by_patient[r["patient_ID"]].add(r["tissue"])
        prep[r["tissue"]].add(r["prep"])

    combos = collections.Counter(
        tuple(sorted(v)) for v in by_patient.values())
    log(f"\n=== 1. compartment inventory ({len(by_patient)} human individuals, "
        f"{len(rows)} libraries) ===")
    inv = []
    for combo, n in sorted(combos.items(), key=lambda kv: -kv[1]):
        log(f"  {'+'.join(combo):22s} {n:3d} individuals")
        inv.append({"tissues": "+".join(combo), "n_individuals": n})

    def n_with(*t):
        return sum(1 for v in by_patient.values() if set(t) <= v)

    pairs = [("Tumor", "Normal"), ("Tumor", "PBMC"), ("Normal", "PBMC")]
    log("\n  paired n by contrast (individuals holding BOTH compartments):")
    for a, b in pairs:
        log(f"    {a} vs {b:8s} n = {n_with(a, b)}")
    log(f"    all three          n = {n_with('Tumor', 'Normal', 'PBMC')}")
    log("\n  enrichment protocol per tissue (a cross-tissue contrast "
        "inherits any difference here):")
    for t in sorted(prep):
        log(f"    {t:8s} {sorted(x or '(unrecorded)' for x in prep[t])}")

    pd.DataFrame(inv).to_csv(
        os.path.join(OUT, "chemokine_compartment_inventory.csv"), index=False)
    return {f"{a}|{b}": n_with(a, b) for a, b in pairs} | {
        "Tumor|Normal|PBMC": n_with("Tumor", "Normal", "PBMC")}


def main():
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    paired_n = compartment_inventory(log)

    # ---- part 2: measurability ---------------------------------------
    d = np.load(os.path.join(OUT, "pseudobulk_uncorrected_matched.npz"),
                allow_pickle=True)
    ps = PseudobulkSet(d["counts"].astype(float), d["genes"].astype(str),
                       d["patient"], d["condition"], d["lineage"],
                       np.ones(len(d["patient"]))).complete_pairs()
    sp = np.load(os.path.join(OUT, "soup_profile_mean.npz"), allow_pickle=True)
    soup_all = dict(zip(sp["genes"].astype(str), sp["soup"]))
    rho_lin = (pd.read_csv(os.path.join(OUT, "rho_variants.csv"))
               .groupby("lineage").rho_per_sample.median().to_dict())

    # Two separate reasons a panel gene can be untestable, kept apart
    # because they mean different things.  Gating genes are barred by
    # protocol 2.1 whatever their expression; the relaxed detection
    # filter (NK plus at least one witness lineage, docs/DEVIATIONS.md
    # D7) drops genes that are measurable but NK-restricted.
    excluded = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    in_universe = set(ps.genes[detection_filter(ps, require_nk=True,
                                                min_lineages=2)])

    genes = list(ps.genes)
    panel_flat = [(g, axis) for axis, gl in PANEL.items() for g in gl]
    panel_flat += [(g, "control_positive") for g in POS_CTRL]
    panel_flat += [(g, "control_negative") for g in NEG_CTRL]
    axis_of = dict(panel_flat)
    idx = {g: genes.index(g) for g, _ in panel_flat if g in genes}
    missing = [g for g, _ in panel_flat if g not in idx]
    if missing:
        log(f"\nnot present in the matrix at all: {missing}")

    rows = []
    for lin in LINEAGES:
        s = ps.subset_lineage(lin)
        lib = s.counts.sum(axis=0)
        cpm = s.counts / np.maximum(lib, 1) * 1e6
        total = float(lib.mean())
        rho = float(rho_lin.get(lin, 0.1))
        for g, i in idx.items():
            obs_mean = float(s.counts[i].mean())
            f = rho * total * soup_all.get(g, 0.0) / max(obs_mean, 1e-9)
            rows.append({
                "gene": g, "axis": axis_of[g], "lineage": lin,
                "detected": bool((s.counts[i] >= MIN_COUNT).mean() >= MIN_FRAC),
                "frac_ge_floor": float((s.counts[i] >= MIN_COUNT).mean()),
                "detection_rate": float((s.counts[i] > 0).mean()),
                "mean_CPM": float(cpm[i].mean()),
                "soup_frac": float(min(f, 1.0)),
            })
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(OUT, "chemokine_measurability.csv"), index=False)

    n_pat = len(set(ps.patient))
    nk_lib = float(ps.subset_lineage("NK").counts.sum(axis=0).mean())
    nkr = r[r.lineage == "NK"].set_index("gene")
    myr = r[r.lineage == "Myeloid"].set_index("gene")

    # ---- calibration, identical in construction to the complosome check
    soup_ceiling = float(nkr.loc[[g for g in NEG_CTRL if g in idx],
                                 "soup_frac"].min())
    real_floor = float(nkr.loc[[g for g in POS_CTRL if g in idx],
                               "soup_frac"].max())
    band = [float(nkr.loc[g, "mean_CPM"] / max(myr.loc[g, "mean_CPM"], 1e-9))
            for g in ("C1QA", "LYZ") if g in idx]
    log("\n=== 2. ambient calibration inside NK ===")
    log(f"  genuinely-NK genes read soup_frac up to     {real_floor:.3f}")
    log(f"  known-ambient genes read soup_frac down to  {soup_ceiling:.3f} "
        f"(their true value is 1.0)")
    log(f"  NK/Myeloid CPM for myeloid-only genes:      "
        f"{'  '.join(f'{b:.3f}' for b in band)}  (pure pickup)")

    verdicts = []
    for g, _ in panel_flat:
        if g not in idx:
            verdicts.append({"gene": g, "axis": axis_of[g],
                             "verdict": "ABSENT FROM MATRIX"})
            continue
        nk, my = nkr.loc[g], myr.loc[g]
        ratio = float(nk.mean_CPM / max(my.mean_CPM, 1e-9))
        # The NK/Myeloid ratio only discriminates where myeloid is the
        # dominant source; for broadly-expressed genes it flags ambient by
        # construction (LEADS_CLOSED.md, CTSL/C5AR1).  Record that as
        # untestable-this-way rather than as a negative.
        ratio_applies = bool(my.mean_CPM > nk.mean_CPM)
        if not nk.detected:
            v = "BELOW DETECTION FLOOR"
        elif nk.soup_frac >= soup_ceiling:
            v = "AMBIENT-INDISTINGUISHABLE"
        elif ratio_applies and ratio <= max(band) * 1.5:
            v = "AMBIENT-INDISTINGUISHABLE (ratio)"
        else:
            v = "MEASURABLE"
        verdicts.append({
            "gene": g, "axis": axis_of[g],
            "NK_detected": bool(nk.detected),
            "NK_frac_ge_floor": round(float(nk.frac_ge_floor), 3),
            "NK_detection_rate": round(float(nk.detection_rate), 3),
            "NK_CPM": round(float(nk.mean_CPM), 2),
            "Myeloid_CPM": round(float(my.mean_CPM), 2),
            "NK_soup_frac": round(float(nk.soup_frac), 3),
            "NK_over_Myeloid": round(ratio, 3),
            "ratio_test_applies": ratio_applies,
            "gating_excluded": g in excluded,
            "passes_witness_filter": g in in_universe,
            "verdict": v,
        })
    v = pd.DataFrame(verdicts)

    log(f"\n=== 3. measurability in NK (floor: >={MIN_COUNT} counts in "
        f">={MIN_FRAC:.0%} of NK pseudobulk samples) ===")
    for axis in list(PANEL) + ["control_positive", "control_negative"]:
        sub = v[v.axis == axis]
        ok = sub[sub.verdict == "MEASURABLE"]
        log(f"\n  {axis}  ({len(ok)}/{len(sub)} measurable)")
        for _, row in sub.iterrows():
            if row.verdict == "ABSENT FROM MATRIX":
                log(f"    {row.gene:10s} {'-':>7} {'-':>8} {'-':>7}   "
                    f"ABSENT FROM MATRIX")
                continue
            log(f"    {row.gene:10s} det={row.NK_frac_ge_floor:5.2f} "
                f"CPM={row.NK_CPM:8.2f} soup={row.NK_soup_frac:5.3f} "
                f"NK/Mye={row.NK_over_Myeloid:7.3f}"
                f"{'' if row.ratio_test_applies else '*'}   {row.verdict}")
    log("\n  det = fraction of NK pseudobulk samples reaching the floor "
        f"(>={MIN_COUNT} counts); at the mean NK library size of "
        f"{nk_lib:,.0f} counts the floor is ~{10 / nk_lib * 1e6:.0f} CPM.")
    log("  * NK/Myeloid ratio does not discriminate for this gene "
        "(NK >= Myeloid); the soup fraction is the only scale that applies.")

    # ---- part 3: effect sizes for what survives, descriptive only ----
    qpath = os.path.join(OUT, "quadrants_relaxed_d0.5_p95.csv")
    q = pd.read_csv(qpath)
    keep = [c for c in q.columns if c in ("gene", "quadrant") or
            c.endswith("_log2FC") or c.endswith("_SE")]
    v = v.merge(q[keep], on="gene", how="left")
    v.to_csv(os.path.join(OUT, "chemokine_panel_verdict.csv"), index=False)

    meas = (v[(v.verdict == "MEASURABLE") & v.axis.isin(PANEL)]
            .sort_values("NK_log2FC"))
    log(f"\n=== 4. paired tumour-vs-normal in NK, for the "
        f"{len(meas)} measurable panel genes ===")
    log("    DESCRIPTIVE ONLY -- the quadrant reading of this dataset is "
        "retracted; three")
    log("    NK-specific technical effects make any NK-vs-others contrast "
        "unreadable here.")
    log(f"    {'gene':10s} {'NK_log2FC':>10s} {'NK_SE':>7s} "
        f"{'CD8T_log2FC':>12s} {'quadrant':>22s}")
    for _, row in meas.iterrows():
        if pd.isna(row.get("NK_log2FC")):
            why = ("barred as a gating gene (2.1)" if row.gating_excluded
                   else "NK-restricted: no witness lineage (D7)"
                   if not row.passes_witness_filter
                   else "measurable here, absent from the 27-patient "
                   "Phase 5 table")
            log(f"    {row.gene:10s}   -- untested: {why}")
            continue
        log(f"    {row.gene:10s} {row.NK_log2FC:10.3f} {row.NK_SE:7.3f} "
            f"{row.CD8T_log2FC:12.3f} {str(row.quadrant):>22s}")

    # ---- the headline counts -----------------------------------------
    n_panel = int((v.axis.isin(PANEL)).sum())
    n_meas = int(((v.verdict == "MEASURABLE") & v.axis.isin(PANEL)).sum())
    rec_axes = ["receptor_inflammatory", "receptor_homing",
                "receptor_atypical"]
    n_rec = int(v.axis.isin(rec_axes).sum())
    n_rec_meas = int(((v.verdict == "MEASURABLE") &
                      v.axis.isin(rec_axes)).sum())
    n_lig = int(v.axis.isin(["ligand_lymphoid",
                             "ligand_myeloid_stromal"]).sum())
    n_lig_meas = int(((v.verdict == "MEASURABLE") &
                      v.axis.isin(["ligand_lymphoid",
                                   "ligand_myeloid_stromal"])).sum())
    log("\n=== SUMMARY ===")
    log(f"  panel genes                     {n_panel}")
    log(f"  measurable in NK                {n_meas}")
    log(f"    of which receptors            {n_rec_meas}/{n_rec}")
    log(f"    of which ligands              {n_lig_meas}/{n_lig}")
    log(f"  tumour vs normal, paired n      {paired_n['Tumor|Normal']} "
        f"individuals in the metadata; {n_pat} carry a complete "
        f"five-lineage pair here")
    log(f"  blood vs tissue, paired n       {paired_n['Tumor|PBMC']} "
        f"individuals")
    log("DONE")

    with open(os.path.join(OUT, "chemokine_measurability.log"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
