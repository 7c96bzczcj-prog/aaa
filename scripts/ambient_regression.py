"""Regress the ambient component out of log2FC instead of subtracting counts.

Why this is not another round of tuning. The two failed arms both
subtracted counts, and both failed through the same two mechanisms:
per-sample rho noise injected into every gene, and unstable ratios of
near-zero residuals. This approach touches neither. It never modifies a
count, so nothing is floored at zero and no per-sample quantity enters;
it operates on the gene-level summary statistics that are already
computed.

MODEL
-----
Ambient adds `rho_c * total_c * soup_g` to the observed counts of gene g
in condition c. Its contribution to the observed log2FC is therefore
governed by how much of that gene's signal is soup rather than by the
soup's abundance alone: a gene the lineage expresses strongly is barely
moved, a gene it does not express at all is moved by the full ambient
difference.

So the regressor is the gene's *soup fraction* in that lineage,

    f_g = rho * total * soup_g / observed_g   (clipped to [0, 1])

and the correction is the residual of

    observed log2FC ~ s(f_g)

fitted across genes within a lineage. `s` is a monotone step-wise fit
over bins of f, which avoids assuming a functional form.

ACCEPTANCE, in the same units as the claim
------------------------------------------
Immunoglobulin genes in NK, CD8T, CD4T and myeloid, and granzyme genes
in B, have a true log2FC of exactly zero. They should sit *on* the
fitted curve, and their residuals should be near zero afterwards. They
are a small minority of genes, so they do not drive the fit, but to be
strict the fit is also repeated with them excluded and the residuals
re-checked out of sample.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nkmine.de import paired_de  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, detection_filter  # noqa: E402
from nkmine.quadrant import LINEAGES  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")

ZERO_TRUTH = {
    "NK": ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"],
    "CD8T": ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"],
    "CD4T": ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"],
    "Myeloid": ["IGKC", "IGHG1", "IGHG3", "IGHA1", "IGLC2", "IGLC3", "MZB1", "JCHAIN"],
    "B": ["GZMB", "GZMA", "PRF1", "CTSW", "GZMH", "NKG7", "KLRD1"],
}
N_BINS = 25


def binned_fit(f: np.ndarray, y: np.ndarray, mask_fit: np.ndarray) -> np.ndarray:
    """Median of y within bins of f, interpolated back to every gene."""
    edges = np.quantile(f[mask_fit], np.linspace(0, 1, N_BINS + 1))
    edges = np.unique(edges)
    if len(edges) < 3:
        return np.zeros_like(y)
    centres, meds = [], []
    for i in range(len(edges) - 1):
        sel = mask_fit & (f >= edges[i]) & (f <= edges[i + 1])
        if sel.sum() >= 10:
            centres.append(float(np.median(f[sel])))
            meds.append(float(np.median(y[sel])))
    if len(centres) < 3:
        return np.zeros_like(y)
    return np.interp(f, np.array(centres), np.array(meds))


def main():
    d = np.load(os.path.join(OUT, "pseudobulk_uncorrected_matched.npz"),
                allow_pickle=True)
    ps = PseudobulkSet(d["counts"].astype(float), d["genes"], d["patient"],
                       d["condition"], d["lineage"],
                       np.ones(len(d["patient"]))).complete_pairs()
    sp = np.load(os.path.join(OUT, "soup_profile_mean.npz"), allow_pickle=True)
    soup_all = dict(zip(sp["genes"].astype(str), sp["soup"]))
    rho_tab = pd.read_csv(os.path.join(OUT, "rho_variants.csv"))
    rho_lin = rho_tab.groupby("lineage").rho_per_sample.median().to_dict()

    excl = set(open(os.path.join(ROOT, "excluded_genes.txt")).read().split())
    ok = ~np.isin(ps.genes, list(excl))
    ps = PseudobulkSet(ps.counts[ok], ps.genes[ok], ps.patient, ps.condition,
                       ps.lineage, ps.n_cells)
    keep = detection_filter(ps, require_nk=True, min_lineages=2)
    ps = PseudobulkSet(ps.counts[keep], ps.genes[keep], ps.patient,
                       ps.condition, ps.lineage, ps.n_cells)
    soup = np.array([soup_all.get(g, 0.0) for g in ps.genes])
    print(f"genes {ps.counts.shape[0]}  patients {len(set(ps.patient))}", flush=True)

    rows, out = [], {"gene": ps.genes}
    for lin in LINEAGES:
        s = ps.subset_lineage(lin)
        de = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
        obs_mean = s.counts.mean(axis=1)
        total = float(s.counts.sum(axis=0).mean())
        rho = float(rho_lin.get(lin, 0.1))
        f = np.clip(rho * total * soup / np.maximum(obs_mean, 1e-9), 0.0, 1.0)

        y = de.log2fc
        truth = np.isin(ps.genes, ZERO_TRUTH[lin])
        # strict: fit excluding the known-zero genes, so their residual is
        # genuinely out of sample
        fit_out = binned_fit(f, y, mask_fit=~truth)
        resid = y - fit_out

        rows.append({
            "lineage": lin, "rho_used": rho, "n_zero_truth": int(truth.sum()),
            "soup_frac_median": float(np.median(f)),
            "soup_frac_zero_truth": float(np.median(f[truth])) if truth.any() else np.nan,
            "before_median_abs": float(np.median(np.abs(y[truth]))) if truth.any() else np.nan,
            "after_median_abs": float(np.median(np.abs(resid[truth]))) if truth.any() else np.nan,
            "background_before": float(np.median(np.abs(y[~truth]))),
            "background_after": float(np.median(np.abs(resid[~truth]))),
            "sd_before": float(np.std(y)), "sd_after": float(np.std(resid)),
        })
        out[f"{lin}_log2FC_raw"] = y
        out[f"{lin}_log2FC_corrected"] = resid
        out[f"{lin}_soup_frac"] = f
        out[f"{lin}_SE"] = de.se

    r = pd.DataFrame(rows)
    pd.DataFrame(out).to_csv(os.path.join(OUT, "ambient_regression_genes.csv"),
                             index=False)
    r.to_csv(os.path.join(OUT, "ambient_regression_summary.csv"), index=False)

    print("\n=== known-zero controls, |log2FC| before vs after "
          "(fit excluded them, so this is out of sample) ===", flush=True)
    print(r[["lineage", "soup_frac_zero_truth", "before_median_abs",
             "after_median_abs", "background_before", "background_after"]]
          .round(3).to_string(index=False), flush=True)

    imp = 1 - r.after_median_abs / r.before_median_abs
    print("\nreduction in the known-zero controls: " +
          "  ".join(f"{l}={v:.0%}" for l, v in zip(r.lineage, imp)), flush=True)
    bg_change = r.background_after / r.background_before - 1
    print("change in background |log2FC|: " +
          "  ".join(f"{l}={v:+.1%}" for l, v in zip(r.lineage, bg_change)), flush=True)
    print("\nSE is untouched by construction (counts are never modified): " +
          "  ".join(f"{l}={v:.3f}" for l, v in
                    zip(r.lineage, [float(np.median(out[f'{l}_SE'])) for l in r.lineage])),
          flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
