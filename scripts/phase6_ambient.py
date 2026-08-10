"""Phase 6.1 - the shared-soup contrast.

GSE154826 contains both designs at once, which is what makes this
possible without any decontamination tool:

  * 8 patients have their tumour and adjacent-normal cells pooled into
    ONE droplet emulsion and separated afterwards by hashing.  Both
    conditions therefore sit in the same ambient soup, and whatever the
    soup contributes is added to both sides of the paired contrast and
    largely cancels.

  * 19 patients have their two conditions built as separate libraries.
    Each condition has its own soup, so the ambient contribution does
    NOT cancel -- and because ambient is added to every lineage alike,
    it pushes all five lineages in the same direction at once.

The protocol's prediction is directional and falsifiable: cross-lineage
concordance should be **higher** in the separately-built set, and the
gap estimates ambient's contribution.

The comparison is power-matched.  The separate group has 19 patients
against the shared group's 8, and correlation between lineages rises
as estimates get less noisy, so an unmatched comparison would find the
predicted gap whether or not ambient exists.  The separate group is
therefore repeatedly subsampled to 8 patients.
"""

from __future__ import annotations

import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.de import paired_de  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
# Patient 581 contributes BOTH a shared-emulsion library and two separate
# libraries, so its pseudobulk mixes the two designs and it cannot serve as
# a clean member of either group.  Excluded, leaving n=7 vs n=19.
SHARED = {"584", "593", "596", "626", "630", "695", "706"}
MIXED_DESIGN = {"581"}
N_DRAWS = 25
DELTA = 0.5


def subset_patients(ps: PseudobulkSet, keep: set) -> PseudobulkSet:
    m = np.isin(ps.patient, list(keep))
    return PseudobulkSet(ps.counts[:, m], ps.genes, ps.patient[m],
                         ps.condition[m], ps.lineage[m], ps.n_cells[m])


def fit(ps: PseudobulkSet) -> dict:
    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    return de


def concordance(de: dict) -> dict:
    """Cross-lineage sharing metrics.

    `mean_pairwise_r` needs no null calibration and no classification,
    so it is the primary metric: it is a direct readout of how much the
    five lineages move together.
    """
    lfc = {l: de[l].log2fc for l in LINEAGES}
    rs = []
    for a, b in itertools.combinations(LINEAGES, 2):
        rs.append(np.corrcoef(lfc[a], lfc[b])[0, 1])

    # how often do >=4 lineages move together, by CI alone
    sig = {l: ((de[l].ci_lo_95 > 0) | (de[l].ci_hi_95 < 0))
              & (np.abs(de[l].log2fc) >= DELTA) for l in LINEAGES}
    sign = {l: np.sign(de[l].log2fc) for l in LINEAGES}
    n_sig = np.sum([sig[l] for l in LINEAGES], axis=0)
    n_up = np.sum([sig[l] & (sign[l] > 0) for l in LINEAGES], axis=0)
    n_dn = np.sum([sig[l] & (sign[l] < 0) for l in LINEAGES], axis=0)
    concordant4 = (np.maximum(n_up, n_dn) >= 4)

    # NK-vs-rest: does NK move with the pack?
    others = np.mean([lfc[l] for l in LINEAGES if l != "NK"], axis=0)
    r_nk_rest = float(np.corrcoef(lfc["NK"], others)[0, 1])

    return {
        "mean_pairwise_r": float(np.mean(rs)),
        "min_pairwise_r": float(np.min(rs)),
        "max_pairwise_r": float(np.max(rs)),
        "r_NK_vs_rest": r_nk_rest,
        "frac_ge4_concordant": float(concordant4.mean()),
        "n_genes": int(len(lfc["NK"])),
        "median_abs_lfc_NK": float(np.median(np.abs(lfc["NK"]))),
    }


def main():
    d = np.load(os.path.join(OUT, "pseudobulk_raw.npz"), allow_pickle=True)
    ps = PseudobulkSet(d["counts"], d["genes"], d["patient"], d["condition"],
                       d["lineage"], np.ones(len(d["patient"])))
    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    keep = detection_filter(ps, require_nk=True, min_lineages=2)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                       ps.condition, ps.lineage, ps.n_cells)

    pats = sorted(set(ps.patient.tolist()))
    shared = sorted(p for p in pats if p in SHARED)
    separate = sorted(p for p in pats if p not in SHARED and p not in MIXED_DESIGN)
    print(f"genes {ps.counts.shape[0]} | shared-soup n={len(shared)} "
          f"| separate-library n={len(separate)} "
          f"| excluded for mixed design: {sorted(MIXED_DESIGN)}", flush=True)

    # Confounder check: the two groups must not differ systematically in the
    # things that also drive cross-lineage correlation.
    diag = []
    for nm, grp in (("shared", shared), ("separate", separate)):
        m = np.isin(ps.patient, list(grp))
        tot = ps.counts[:, m].sum(axis=0)
        diag.append({"group": nm, "n_patients": len(grp),
                     "median_pseudobulk_total_counts": float(np.median(tot)),
                     "min": float(tot.min()), "max": float(tot.max())})
    dg = pd.DataFrame(diag)
    dg.to_csv(os.path.join(OUT, "phase6_group_depth_check.csv"), index=False)
    print(dg.to_string(index=False), flush=True)

    rows = []
    m_shared = concordance(fit(subset_patients(ps, set(shared))))
    m_shared.update(group="shared_soup", draw=-1, n_patients=len(shared))
    rows.append(m_shared)
    print(f"shared soup   : r={m_shared['mean_pairwise_r']:.4f} "
          f"ge4concordant={m_shared['frac_ge4_concordant']:.4f}", flush=True)

    m_full = concordance(fit(subset_patients(ps, set(separate))))
    m_full.update(group="separate_full", draw=-1, n_patients=len(separate))
    rows.append(m_full)
    print(f"separate full : r={m_full['mean_pairwise_r']:.4f} "
          f"ge4concordant={m_full['frac_ge4_concordant']:.4f}", flush=True)

    # power-matched: draw 8 of the 19 separate-library patients, repeatedly
    rng = np.random.default_rng(0)
    for k in range(N_DRAWS):
        pick = set(rng.choice(separate, size=len(shared), replace=False).tolist())
        m = concordance(fit(subset_patients(ps, pick)))
        m.update(group="separate_matched", draw=k, n_patients=len(shared))
        rows.append(m)
        if (k + 1) % 5 == 0:
            print(f"  matched draw {k+1}/{N_DRAWS}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase6_shared_soup_contrast.csv"), index=False)

    matched = df[df.group == "separate_matched"]
    r_s = m_shared["mean_pairwise_r"]
    r_m = matched.mean_pairwise_r
    # one-sided empirical p: how often does a power-matched separate draw
    # fail to exceed the shared-soup value?
    p_emp = float((r_m <= r_s).mean())
    print("\n=== POWER-MATCHED COMPARISON (matched on n) ===", flush=True)
    print(f"  shared soup            mean pairwise r = {r_s:.4f}", flush=True)
    print(f"  separate (matched n)    mean pairwise r = {r_m.mean():.4f} "
          f"[{r_m.min():.4f}, {r_m.max():.4f}]", flush=True)
    print(f"  difference (separate - shared)         = {r_m.mean()-r_s:+.4f}", flush=True)
    print(f"  draws not exceeding shared soup        = {p_emp:.3f}", flush=True)
    f_s = m_shared["frac_ge4_concordant"]
    f_m = matched.frac_ge4_concordant
    print(f"  frac >=4 concordant: shared {f_s:.4f} vs separate {f_m.mean():.4f} "
          f"({f_m.mean()-f_s:+.4f})", flush=True)
    print(f"  r(NK, mean of others): shared {m_shared['r_NK_vs_rest']:.4f} vs "
          f"separate {matched.r_NK_vs_rest.mean():.4f}", flush=True)
    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
