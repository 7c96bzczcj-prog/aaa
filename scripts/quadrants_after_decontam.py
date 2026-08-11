"""Quadrant assignment before and after ambient correction, matched A/B.

Both arms use the identical aggregation (all called cells per
patient x condition x lineage, summed across that patient's libraries),
so the only difference between them is the soup subtraction. The
previously reported analysis used balanced subsampling on top of that;
this pair deliberately does not, because mixing the two changes at once
would make the comparison uninterpretable.

Consequences of dropping balancing are stated rather than hidden: NK's
standard errors will be larger here than in the balanced analysis, so
absolute quadrant counts are not comparable to the earlier numbers.
Only the corrected-vs-uncorrected contrast is.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nkmine.de import paired_de  # noqa: E402
from nkmine.null_calibration import calibrate  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter, se_balance_report  # noqa: E402
from nkmine.quadrant import LINEAGES, classify_table  # noqa: E402
from s4_diagnosis import nk_saturation  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
DELTA = 0.5
IG = ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"]
CYTO_IN_B = ["GZMB", "GZMA", "PRF1", "CTSW", "GZMH"]


def load(name: str) -> PseudobulkSet:
    d = np.load(os.path.join(OUT, name), allow_pickle=True)
    ps = PseudobulkSet(d["counts"].astype(float), d["genes"], d["patient"],
                       d["condition"], d["lineage"], np.ones(len(d["patient"])))
    return ps.complete_pairs()


def run(ps: PseudobulkSet, label: str, shared_genes) -> tuple:
    keep = np.isin(ps.genes, shared_genes)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                       ps.condition, ps.lineage, ps.n_cells)
    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    rep = se_balance_report(de)
    null = calibrate(ps, case_label="Tumor", n_perm=100, rng=np.random.default_rng(1))
    sat = nk_saturation(ps.genes)
    res = classify_table(de, delta=DELTA, null_stats=null, equiv_key="p95",
                         saturated={l: (sat if l == "NK" else np.zeros(len(sat), bool))
                                    for l in LINEAGES})
    res.to_csv(os.path.join(OUT, f"quadrants_{label}.csv"), index=False)

    print(f"\n--- {label} ---", flush=True)
    print(f"  genes {len(res)}  patients {len(set(ps.patient))}", flush=True)
    print(f"  median SE: " + "  ".join(
        f"{l}={np.median(de[l].se):.3f}" for l in LINEAGES), flush=True)
    print(f"  S3: {rep['verdict'][:60]}", flush=True)
    print(f"  quadrants: {res.quadrant.value_counts().to_dict()}", flush=True)

    # the two known-zero controls
    def med(genes, lin):
        m = res[res.gene.isin(genes)]
        return float(np.median(np.abs(m[f"{lin}_log2FC"]))) if len(m) else np.nan
    print(f"  |log2FC| immunoglobulin in NK      : {med(IG,'NK'):.3f}", flush=True)
    print(f"  |log2FC| granzyme/perforin in B    : {med(CYTO_IN_B,'B'):.3f}", flush=True)

    # Q2 funnel
    s1 = res.NK_changed
    s2 = s1 & res.CD8T_changed
    s3 = s2 & (np.sign(res.NK_log2FC) == np.sign(res.CD8T_log2FC))
    s4 = s3 & res.B_equivalent
    s5 = s4 & res.Myeloid_equivalent
    print(f"  Q2 funnel: NK={int(s1.sum())} +CD8T={int(s2.sum())} "
          f"+concordant={int(s3.sum())} +B eq={int(s4.sum())} +Mye eq={int(s5.sum())}",
          flush=True)
    return res, de


def main():
    cor = load("pseudobulk_decontaminated.npz")
    raw = load("pseudobulk_uncorrected_matched.npz")
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())

    # one shared gene set, chosen on the UNCORRECTED data so that correction
    # cannot change which genes are tested (that would confound the A/B)
    tmp = PseudobulkSet(raw.counts, raw.genes, raw.patient, raw.condition,
                        raw.lineage, raw.n_cells)
    ok = ~np.isin(tmp.genes, list(excl))
    tmp = PseudobulkSet(tmp.counts[ok], tmp.genes[ok], tmp.patient,
                        tmp.condition, tmp.lineage, tmp.n_cells)
    keep = detection_filter(tmp, require_nk=True, min_lineages=2)
    shared = tmp.genes[keep]
    print(f"shared gene panel (fixed on uncorrected data): {len(shared)}", flush=True)

    res_raw, _ = run(raw, "uncorrected_matched", shared)
    res_cor, _ = run(cor, "decontaminated", shared)

    # Third arm, and the scientifically correct pipeline: filter AFTER
    # correcting.  Fixing the panel on uncorrected data (arms 1-2) isolates
    # the correction's effect, but it also lets genes that are pure soup
    # survive into the corrected analysis, where their near-zero residuals
    # give wild ratios.  Re-deriving the detection floor on corrected counts
    # is what removes them.
    ok2 = ~np.isin(cor.genes, list(excl))
    c2 = PseudobulkSet(cor.counts[ok2], cor.genes[ok2], cor.patient,
                       cor.condition, cor.lineage, cor.n_cells)
    keep2 = detection_filter(c2, require_nk=True, min_lineages=2)
    shared2 = c2.genes[keep2]
    lost = sorted(set(shared) - set(shared2))
    print(f"\npanel re-derived on corrected counts: {len(shared2)} genes "
          f"({len(set(shared)-set(shared2))} dropped, {len(set(shared2)-set(shared))} added)",
          flush=True)
    print(f"  immunoglobulin genes dropped: "
          f"{[g for g in IG if g in set(shared) and g not in set(shared2)]}", flush=True)
    print(f"  cytotoxic genes dropped: "
          f"{[g for g in CYTO_IN_B if g in set(shared) and g not in set(shared2)]}", flush=True)
    res_cor2, _ = run(cor, "decontaminated_refiltered", shared2)

    m = res_raw[["gene", "quadrant"]].merge(
        res_cor[["gene", "quadrant"]], on="gene", suffixes=("_raw", "_cor"))
    print("\n=== quadrant migration, uncorrected -> decontaminated ===", flush=True)
    print(pd.crosstab(m.quadrant_raw, m.quadrant_cor).to_string(), flush=True)
    m.to_csv(os.path.join(OUT, "quadrant_migration_decontam.csv"), index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
