"""Task P-1 (topic 2, Q1-Q6): the PIEZO1 expression baseline in NK cells.

WHAT THIS SCRIPT CAN AND CANNOT DO
----------------------------------
The task specifies three datasets.  Only one of them leaves a usable
trace in this repository:

  GSE154826 (human NSCLC, paired tumour / adjacent normal)
      committed as *derived pseudobulks* -- results/pseudobulk_raw.npz,
      33,668 genes x 270 samples (27 paired patients x 2 conditions x
      5 lineages), balanced-subsampled raw counts.  Usable.

  Dominguez Conde cross-tissue atlas          absent from this repo
  E-MTAB-10176 (mouse Kaede photoconversion)  absent from this repo

`data/` is git-ignored and this container was cloned fresh, so no
cell-level matrix survives for any dataset.  That has one consequence
that shapes everything below: **per-cell detection rate cannot be
recomputed directly**, because it needs the cells.

THE SUBSTITUTE, AND WHY IT IS TRUSTWORTHY HERE
----------------------------------------------
Balanced pseudobulk stores, for each (patient x condition x lineage),
the summed counts over exactly `n_min` cells, where `n_min` is the
size of the rarest lineage in that group (scripts/phase3_5_run.py,
Phase 2.4).  `n_min` is recoverable from results/cell_counts_by_group.csv.
So the per-cell mean UMI of a gene, mu = count / n_min, IS recoverable.

Under any conditionally-Poisson model of per-cell counts -- which
includes negative binomial and every Poisson-mixture used in scRNA-seq --

    P(count > 0) = 1 - E[exp(-lambda)] <= 1 - exp(-E[lambda]) = 1 - exp(-mu)

by Jensen's inequality.  `1 - exp(-mu)` is therefore a strict UPPER
BOUND on detection rate, never an underestimate.  A gene whose bound
falls below the 5% floor is floor-band with certainty; the bound cannot
rescue it.

This is a self-supplied threshold in the sense of the topic-1 rule, so
it is calibrated against ground truth rather than asserted: the repo
records TOX's true cell-level NK detection as 4.3% (tumour) / 4.9%
(normal) (docs/DEVIATIONS.md D10), measured from the cells before they
were lost.  `calibration()` below recomputes TOX from pseudobulk alone
and prints both numbers side by side.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")
LINEAGES = ["NK", "CD8T", "CD4T", "B", "Myeloid"]

# Measurability bands, carried over unchanged from topic 1.
FLOOR, CEILING = 0.05, 0.85

# Q2 panel.  Ten mechanosensitive channels, PIEZO1 first.
PANEL = ["PIEZO1", "PIEZO2", "TRPV4", "TRPV2", "TRPM7", "TRPC1",
         "TMEM63A", "TMEM63B", "KCNK2", "KCNK4"]

# Q6 partners: residency side, egress side, lineage controls.
Q6_GENES = {
    "residency": ["ITGA1", "CD69", "ITGAE"],
    "egress": ["KLF2", "S1PR1", "S1PR5", "SELL"],
    "control": ["KLRD1", "NKG7"],
}

# docs/DEVIATIONS.md D10 / README: measured on the cells, before they were lost.
TOX_TRUTH = {"Tumor": 0.043, "Normal": 0.049}


def band(det: float) -> str:
    if not np.isfinite(det):
        return "INSUFFICIENT"
    if det < FLOOR:
        return "floor"
    if det > CEILING:
        return "ceiling"
    return "measurable"


class Bulk:
    """The committed GSE154826 pseudobulk plus the recovered cell counts."""

    def __init__(self):
        t0 = time.time()
        d = np.load(os.path.join(OUT, "pseudobulk_raw.npz"), allow_pickle=True)
        self.genes = np.array([str(x) for x in d["genes"]])
        self.counts = d["counts"].astype(float)
        self.patient = np.array([str(x) for x in d["patient"]])
        self.condition = np.array([str(x) for x in d["condition"]])
        self.lineage = np.array([str(x) for x in d["lineage"]])
        self.total = self.counts.sum(axis=0)
        self.idx = {g: i for i, g in enumerate(self.genes)}

        cc = pd.read_csv(os.path.join(OUT, "cell_counts_by_group.csv"))
        cc["patient"] = cc.patient.astype(str)
        cc = cc[cc.lineage.isin(LINEAGES)]
        wide = (cc.groupby(["patient", "condition", "lineage"])["n_cells"]
                  .sum().unstack("lineage")[LINEAGES])
        self.n_min = wide.min(axis=1)
        self.n_avail = wide
        self.ncell = np.array([self.n_min.get((p, c), np.nan)
                               for p, c in zip(self.patient, self.condition)])
        self.read_seconds = time.time() - t0

        if not np.isfinite(self.ncell).all():
            raise RuntimeError("cell counts missing for some pseudobulk sample")
        # Sanity: recovered n_min must give a believable per-cell UMI depth.
        upc = self.total / self.ncell
        if not (500 < np.median(upc) < 20000):
            raise RuntimeError(f"implausible UMI/cell ({np.median(upc):.0f}); "
                               "n_min mapping is wrong")

    def has(self, gene: str) -> bool:
        return gene in self.idx

    def mask(self, lineage=None, condition=None):
        m = np.ones(len(self.patient), bool)
        if lineage:
            m &= self.lineage == lineage
        if condition:
            m &= self.condition == condition
        return m

    def per_donor(self, gene: str, lineage: str, condition: str) -> pd.DataFrame:
        """One row per donor: counts, cells, CPM, mean UMI/cell, detection bound."""
        i = self.idx[gene]
        m = self.mask(lineage, condition)
        cnt = self.counts[i, m]
        n = self.ncell[m]
        mu = cnt / n
        return pd.DataFrame({
            "patient": self.patient[m], "condition": self.condition[m],
            "lineage": lineage, "gene": gene, "count": cnt, "n_cells": n,
            "cpm": 1e6 * cnt / self.total[m], "mu_umi_per_cell": mu,
            "det_ub": 1 - np.exp(-mu),
            "umi_per_cell": self.total[m] / n,
        })


class MarkerBulk:
    """Fallback for the 26 genes topic 1 barred from its own testing.

    `excluded_genes.txt` holds the RNA panel used for lineage gating, and
    `pseudobulk_raw.npz` drops those rows -- KLRD1 and NKG7, the Q6 control
    genes, among them.  They survive only in the rho-variant pseudobulks,
    which differ in two ways that are stated wherever a number from here is
    used: every gated cell contributes (no Phase 2.4 balancing), and the
    counts are ambient-subtracted (`decontaminate`, rho ~= 0.151 in NK).
    Ambient subtraction lowers a count, so a detection figure from here is
    if anything conservative for a ceiling call.
    """

    def __init__(self, n_avail: pd.DataFrame):
        d = np.load(os.path.join(OUT, "pseudobulk_per_sample.npz"), allow_pickle=True)
        self.genes = np.array([str(x) for x in d["genes"]])
        self.counts = np.clip(d["counts"].astype(float), 0, None)
        self.patient = np.array([str(x) for x in d["patient"]])
        self.condition = np.array([str(x) for x in d["condition"]])
        self.lineage = np.array([str(x) for x in d["lineage"]])
        self.total = self.counts.sum(axis=0)
        self.idx = {g: i for i, g in enumerate(self.genes)}
        self.ncell = np.array([n_avail["NK"].get((p, c), np.nan) if l == "NK"
                               else n_avail[l].get((p, c), np.nan)
                               for p, c, l in zip(self.patient, self.condition,
                                                  self.lineage)])

    def per_donor(self, gene: str, lineage: str, condition: str) -> pd.DataFrame:
        i = self.idx[gene]
        m = (self.lineage == lineage) & (self.condition == condition)
        m &= np.isfinite(self.ncell)
        cnt, n = self.counts[i, m], self.ncell[m]
        mu = cnt / n
        return pd.DataFrame({
            "patient": self.patient[m], "condition": condition,
            "lineage": lineage, "gene": gene, "count": cnt, "n_cells": n,
            "cpm": 1e6 * cnt / self.total[m], "mu_umi_per_cell": mu,
            "det_ub": 1 - np.exp(-mu), "umi_per_cell": self.total[m] / n})


def pooled_detection(df: pd.DataFrame, n_boot: int = 2000, seed: int = 0):
    """Cell-weighted detection bound, with a donor-resampling bootstrap.

    The statistical unit stays the donor: the bootstrap resamples donors,
    never cells, so the interval reflects between-donor variation.
    """
    c, n = df["count"].to_numpy(), df["n_cells"].to_numpy()
    est = 1 - np.exp(-c.sum() / n.sum())
    rng = np.random.default_rng(seed)
    k = len(c)
    if k < 2:
        return est, np.nan, np.nan
    picks = rng.integers(0, k, size=(n_boot, k))
    boot = 1 - np.exp(-c[picks].sum(axis=1) / n[picks].sum(axis=1))
    return est, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def summarise(bulk, gene: str, lineage: str, condition: str,
              source: str = "pseudobulk_raw (balanced, raw counts)") -> dict:
    d = bulk.per_donor(gene, lineage, condition)
    est, lo, hi = pooled_detection(d)
    q = d["det_ub"].quantile([0.25, 0.5, 0.75])
    qc = d["cpm"].quantile([0.25, 0.5, 0.75])
    return {
        "dataset": "GSE154826", "gene": gene, "lineage": lineage,
        "condition": condition, "n_donors": int(len(d)),
        "n_cells_total": int(d["n_cells"].sum()),
        "det_pooled": est, "det_pooled_lo95": lo, "det_pooled_hi95": hi,
        "det_donor_median": float(q[0.5]), "det_donor_q1": float(q[0.25]),
        "det_donor_q3": float(q[0.75]),
        "cpm_donor_median": float(qc[0.5]), "cpm_donor_q1": float(qc[0.25]),
        "cpm_donor_q3": float(qc[0.75]),
        "donors_zero_count": int((d["count"] == 0).sum()),
        "umi_per_cell_median": float(d["umi_per_cell"].median())
        if "umi_per_cell" in d else np.nan,
        "band_pooled": band(est), "band_donor_median": band(float(q[0.5])),
        "donors_in_floor_band": int((d["det_ub"] < FLOOR).sum()),
        "source": source,
    }


def calibration(bulk: Bulk, log) -> pd.DataFrame:
    """Recover TOX from pseudobulk and check it against the recorded truth."""
    rows = []
    for cond, truth in TOX_TRUTH.items():
        d = bulk.per_donor("TOX", "NK", cond)
        est, lo, hi = pooled_detection(d)
        rows.append({"gene": "TOX", "condition": cond,
                     "recorded_cell_level_detection": truth,
                     "poisson_bound_from_pseudobulk": est,
                     "absolute_error_pp": 100 * (est - truth),
                     "source_of_truth": "docs/DEVIATIONS.md D10 / README"})
    df = pd.DataFrame(rows)
    log("\n=== calibration of the detection estimator (§Q1 self-supplied rule) ===")
    log(df.round(4).to_string(index=False))
    worst = float(df.absolute_error_pp.abs().max())
    log(f"worst absolute error: {worst:.2f} percentage points "
        f"({'TIGHT' if worst < 1 else 'LOOSE'} at the floor)")
    return df


def head_to_head(bulk: Bulk, a: str, b: str, condition: str,
                 lineage: str = "NK") -> dict:
    """Is gene `a` above gene `b` in the SAME cells? Unit = donor.

    Both genes are counted in one pseudobulk column, so sequencing depth,
    cell number, dissociation and ambient burden are identical for the two
    and cancel in the within-donor difference.  This is what licenses a
    ranking among genes that are all in the floor band: the band says the
    absolute rate is unreliable, the paired contrast does not depend on it.
    """
    da = bulk.per_donor(a, lineage, condition).set_index("patient")
    db = bulk.per_donor(b, lineage, condition).set_index("patient")
    sh = sorted(set(da.index) & set(db.index))
    x = np.log2(da.loc[sh, "cpm"].to_numpy() + 1) - np.log2(db.loc[sh, "cpm"].to_numpy() + 1)
    n_hi = int((da.loc[sh, "count"].to_numpy() > db.loc[sh, "count"].to_numpy()).sum())
    n_lo = int((da.loc[sh, "count"].to_numpy() < db.loc[sh, "count"].to_numpy()).sum())
    if np.allclose(x, 0):
        return {"condition": condition, "gene_a": a, "gene_b": b, "n_donors": len(sh),
                "mean_log2_diff": 0.0, "paired_t_p": np.nan, "wilcoxon_p": np.nan,
                "donors_a_gt_b": n_hi, "donors_b_gt_a": n_lo,
                "verdict": "NO_VARIATION (both genes are zero in every donor)"}
    tt = stats.ttest_rel(np.log2(da.loc[sh, "cpm"].to_numpy() + 1),
                         np.log2(db.loc[sh, "cpm"].to_numpy() + 1))
    try:
        wp = float(stats.wilcoxon(x).pvalue)
    except ValueError:
        wp = np.nan
    d = float(x.mean())
    v = ("a_above_b" if d > 0 and tt.pvalue < 0.05 else
         "b_above_a" if d < 0 and tt.pvalue < 0.05 else "not_resolved")
    return {"condition": condition, "gene_a": a, "gene_b": b, "n_donors": len(sh),
            "mean_log2_diff": d, "paired_t_p": float(tt.pvalue), "wilcoxon_p": wp,
            "donors_a_gt_b": n_hi, "donors_b_gt_a": n_lo, "verdict": v}


def paired_test(bulk: Bulk, gene: str, lineage: str = "NK") -> dict:
    """Within-patient tumour vs adjacent normal, on log2 CPM. Unit = patient."""
    t = bulk.per_donor(gene, lineage, "Tumor").set_index("patient")
    n = bulk.per_donor(gene, lineage, "Normal").set_index("patient")
    shared = sorted(set(t.index) & set(n.index))
    lt = np.log2(t.loc[shared, "cpm"].to_numpy() + 1)
    ln = np.log2(n.loc[shared, "cpm"].to_numpy() + 1)
    diff = lt - ln
    if np.allclose(diff, 0):
        det0 = 1 - np.exp(-t.loc[shared, "count"].sum() / t.loc[shared, "n_cells"].sum())
        det1 = 1 - np.exp(-n.loc[shared, "count"].sum() / n.loc[shared, "n_cells"].sum())
        return {"gene": gene, "lineage": lineage, "n_pairs": len(shared),
                "mean_log2FC_tumour_vs_normal": 0.0, "ci95_lo": 0.0, "ci95_hi": 0.0,
                "paired_t_p": np.nan, "wilcoxon_p": np.nan,
                "det_tumour": det0, "det_normal": det1,
                "band_tumour": band(det0), "band_normal": band(det1),
                "measurable": False,
                "note": "NO_VARIATION - zero counts in every donor, no test defined"}
    tt = stats.ttest_rel(lt, ln)
    with np.errstate(invalid="ignore"):
        try:
            wp = float(stats.wilcoxon(diff).pvalue)
        except ValueError:
            wp = np.nan
    se = diff.std(ddof=1) / np.sqrt(len(diff))
    crit = stats.t.ppf(0.975, len(diff) - 1)
    det_t = 1 - np.exp(-t.loc[shared, "count"].sum() / t.loc[shared, "n_cells"].sum())
    det_n = 1 - np.exp(-n.loc[shared, "count"].sum() / n.loc[shared, "n_cells"].sum())
    return {
        "gene": gene, "lineage": lineage, "n_pairs": len(shared),
        "mean_log2FC_tumour_vs_normal": float(diff.mean()),
        "ci95_lo": float(diff.mean() - crit * se),
        "ci95_hi": float(diff.mean() + crit * se),
        "paired_t_p": float(tt.pvalue), "wilcoxon_p": wp,
        "det_tumour": det_t, "det_normal": det_n,
        "band_tumour": band(det_t), "band_normal": band(det_n),
        "measurable": band(det_t) == "measurable" or band(det_n) == "measurable",
        "note": "",
    }


def ambient_sensitivity(bulk: Bulk, genes: list[str]) -> pd.DataFrame:
    """How much of the NK signal could be ambient rather than NK-intrinsic?

    Not a correction -- topic 1 established that subtractive decontamination
    fails its own acceptance test on this dataset
    (docs/AMBIENT_CORRECTION_ATTEMPT.md).  This is the weaker, one-directional
    statement that survives that failure: with a fraction rho of an NK cell's
    UMIs drawn from the soup, a gene making up share `s` of the soup receives
    rho * depth * s ambient UMIs per cell whether or not NK transcribes it.
    Subtracting that can only LOWER a detection figure, so it can only
    strengthen a floor call and can only weaken a measurable one.
    """
    sp = np.load(os.path.join(OUT, "soup_profile_mean.npz"), allow_pickle=True)
    share = dict(zip([str(x) for x in sp["genes"]], sp["soup"]))
    rv = pd.read_csv(os.path.join(OUT, "rho_variants.csv"))
    rho = (rv[rv.lineage == "NK"].groupby("condition")["rho_per_sample"].median())
    rows = []
    for g in genes:
        for cond in ("Tumor", "Normal"):
            d = bulk.per_donor(g, "NK", cond)
            depth = d["cpm"].to_numpy() * 0 + (bulk.total[bulk.mask("NK", cond)]
                                               / bulk.ncell[bulk.mask("NK", cond)])
            mu_obs = d["count"].sum() / d["n_cells"].sum()
            mu_amb = float(rho[cond]) * float(np.mean(depth)) * share.get(g, 0.0)
            mu_int = max(mu_obs - mu_amb, 0.0)
            rows.append({"gene": g, "condition": cond,
                         "rho_NK_median": float(rho[cond]),
                         "soup_share": share.get(g, np.nan),
                         "mu_observed": mu_obs, "mu_ambient_expected": mu_amb,
                         "ambient_fraction_of_signal": mu_amb / mu_obs if mu_obs else np.nan,
                         "det_observed": 1 - np.exp(-mu_obs),
                         "det_ambient_adjusted": 1 - np.exp(-mu_int),
                         "band_observed": band(1 - np.exp(-mu_obs)),
                         "band_ambient_adjusted": band(1 - np.exp(-mu_int))})
    return pd.DataFrame(rows)


def independent_artifact_check(bulk: Bulk, genes: list[str]) -> pd.DataFrame:
    """Re-derive the same figures from the other committed pseudobulk.

    pseudobulk_per_sample.npz is built by a different script, over ALL gated
    NK cells rather than the balanced subsample, and with ambient already
    subtracted.  Agreement between the two is not independent biology, but it
    does rule out the balancing step or the recovered `n_min` being what
    produces the floor call.
    """
    alt = MarkerBulk(bulk.n_avail)
    rows = []
    for g in genes:
        for cond in ("Tumor", "Normal"):
            a = summarise(bulk, g, "NK", cond)
            b = summarise(alt, g, "NK", cond, source="pseudobulk_per_sample")
            rows.append({"gene": g, "condition": cond,
                         "det_balanced_raw": a["det_pooled"],
                         "cells_balanced": a["n_cells_total"],
                         "band_balanced_raw": a["band_pooled"],
                         "det_allcells_decontam": b["det_pooled"],
                         "cells_allcells": b["n_cells_total"],
                         "band_allcells_decontam": b["band_pooled"],
                         "same_band": a["band_pooled"] == b["band_pooled"]})
    return pd.DataFrame(rows)


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx])]
    m = len(order)
    adj = p[order] * m / (np.arange(m) + 1)
    out[order] = np.minimum.accumulate(adj[::-1])[::-1].clip(0, 1)
    return out


def main():
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(str(s))

    bulk = Bulk()
    log("=== Task P-1: PIEZO1 expression baseline (GSE154826 arm) ===")
    log(f"read pseudobulk_raw.npz: {bulk.counts.shape[0]} gene rows x "
        f"{bulk.counts.shape[1]} sample columns in {bulk.read_seconds:.2f}s")
    log(f"donors: {len(set(bulk.patient))} paired patients; "
        f"lineages: {sorted(set(bulk.lineage))}")
    log(f"targeted genes read: {len(set(PANEL + sum(Q6_GENES.values(), []) + ['TOX']))}"
        "  (no full-matrix rescan)")
    log(f"median UMI/cell implied by recovered n_min: "
        f"{np.median(bulk.total / bulk.ncell):.0f}")
    pats = sorted(set(bulk.patient))
    avail = bulk.n_avail.loc[[(p_, c) for p_ in pats for c in ("Tumor", "Normal")
                             if (p_, c) in bulk.n_avail.index], "NK"]
    used = {c: int(bulk.ncell[bulk.mask("NK", c)].sum()) for c in ("Tumor", "Normal")}
    log(f"NK cells gated in these 27 paired patients: {int(avail.sum())}; "
        f"NK cells entering the balanced pseudobulk: "
        f"{used['Tumor']} tumour + {used['Normal']} normal = "
        f"{used['Tumor'] + used['Normal']}")

    cal = calibration(bulk, log)
    cal.to_csv(os.path.join(OUT, "p1_calibration.csv"), index=False)

    # ---------------- Q1 ------------------------------------------------
    log("\n=== Q1: PIEZO1 in NK, and in every other gated lineage ===")
    rows, per_donor = [], []
    for lin in LINEAGES:
        for cond in ("Tumor", "Normal"):
            rows.append(summarise(bulk, "PIEZO1", lin, cond))
            per_donor.append(bulk.per_donor("PIEZO1", lin, cond))
    q1 = pd.DataFrame(rows)
    q1.to_csv(os.path.join(OUT, "p1_q1_baseline.csv"), index=False)
    pd.concat(per_donor).to_csv(os.path.join(OUT, "p1_q1_per_donor.csv"), index=False)
    show = ["lineage", "condition", "n_donors", "n_cells_total",
            "umi_per_cell_median", "cpm_donor_median",
            "det_pooled", "det_pooled_lo95", "det_pooled_hi95", "band_pooled",
            "donors_in_floor_band"]
    log(q1[show].round(4).to_string(index=False))

    nk = q1[(q1.lineage == "NK")]
    log("\nNote on the cross-lineage rows: CPM is depth-normalised and "
        "detection is not.\n  PIEZO1's CPM is flat across the five lineages, "
        "so the detection spread below tracks\n  library depth (UMI/cell), "
        "not a difference in relative abundance.")
    log(f"\nVERDICT Q1 (GSE154826): PIEZO1 in NK is "
        f"{nk.band_pooled.unique().tolist()} in both conditions "
        f"({100*nk.det_pooled.min():.2f}-{100*nk.det_pooled.max():.2f}% detection).")

    log("\n--- Q1 sensitivity 1: how much of it could be ambient ---")
    amb = ambient_sensitivity(bulk, ["PIEZO1"] + [g for g in PANEL if g != "PIEZO1"])
    amb.to_csv(os.path.join(OUT, "p1_ambient_sensitivity.csv"), index=False)
    log(amb[amb.gene.isin(["PIEZO1", "TRPV2", "TRPM7", "TMEM63A"])]
        [["gene", "condition", "rho_NK_median", "ambient_fraction_of_signal",
          "det_observed", "det_ambient_adjusted", "band_ambient_adjusted"]]
        .round(4).to_string(index=False))

    log("\n--- Q1 sensitivity 2: the same call from the other pseudobulk artifact ---")
    ind = independent_artifact_check(bulk, ["PIEZO1", "TRPV2", "TOX"])
    ind.to_csv(os.path.join(OUT, "p1_artifact_crosscheck.csv"), index=False)
    log(ind.round(4).to_string(index=False))

    # ---------------- Q2 ------------------------------------------------
    log("\n=== Q2: ten mechanosensitive channels, side by side in NK ===")
    rows = []
    for g in PANEL:
        for cond in ("Tumor", "Normal"):
            rows.append(summarise(bulk, g, "NK", cond))
    q2 = pd.DataFrame(rows)
    q2["rank_in_panel"] = (q2.groupby("condition")["det_pooled"]
                             .rank(ascending=False, method="min").astype(int))
    q2.to_csv(os.path.join(OUT, "p1_q2_panel.csv"), index=False)
    log(q2[["gene", "condition", "cpm_donor_median", "det_pooled",
            "det_pooled_lo95", "det_pooled_hi95", "band_pooled",
            "rank_in_panel"]].round(4).to_string(index=False))
    for cond in ("Tumor", "Normal"):
        s = q2[q2.condition == cond].sort_values("det_pooled", ascending=False)
        log(f"  {cond}: highest = {s.gene.iloc[0]} "
            f"({100*s.det_pooled.iloc[0]:.2f}%), PIEZO1 rank "
            f"{int(q2[(q2.condition==cond)&(q2.gene=='PIEZO1')].rank_in_panel.iloc[0])}/10")

    rows = []
    for g in PANEL:
        for lin in LINEAGES:
            for cond in ("Tumor", "Normal"):
                rows.append(summarise(bulk, g, lin, cond))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "p1_q2_panel_all_lineages.csv"),
                              index=False)

    log("\n--- Q2 head-to-head: PIEZO1 against each other channel, "
        "within donor, same cells ---")
    h2h = pd.DataFrame([head_to_head(bulk, "PIEZO1", g, c)
                        for c in ("Tumor", "Normal")
                        for g in PANEL if g != "PIEZO1"])
    h2h["q_BH"] = bh(h2h.paired_t_p)
    h2h.to_csv(os.path.join(OUT, "p1_q2_head_to_head.csv"), index=False)
    log(h2h[["condition", "gene_b", "n_donors", "mean_log2_diff", "paired_t_p",
             "q_BH", "donors_a_gt_b", "donors_b_gt_a", "verdict"]]
        .round(4).to_string(index=False))
    beaten = h2h[(h2h.verdict == "b_above_a") & (h2h.q_BH < 0.05)]
    log(f"  channels significantly ABOVE PIEZO1 (BH q<0.05): "
        f"{sorted(set(beaten.gene_b)) or 'none'}")

    # ---------------- Q4 (human arm only) --------------------------------
    log("\n=== Q4 (human arm): tumour vs paired adjacent normal, NK ===")
    tests = [paired_test(bulk, g) for g in PANEL]
    q4 = pd.DataFrame(tests)
    q4["paired_t_q_BH"] = bh(q4.paired_t_p)
    q4.to_csv(os.path.join(OUT, "p1_q4_paired_human.csv"), index=False)
    log(q4[["gene", "n_pairs", "mean_log2FC_tumour_vs_normal", "ci95_lo", "ci95_hi",
            "paired_t_p", "paired_t_q_BH", "band_tumour", "band_normal",
            "measurable"]].round(4).to_string(index=False))

    # ---------------- Q6 (measurability arm only) ------------------------
    log("\n=== Q6: measurability of the correlation partners (the "
        "correlation itself needs cells) ===")
    marker = MarkerBulk(bulk.n_avail)
    rows = []
    for side, gs in Q6_GENES.items():
        for g in gs:
            for cond in ("Tumor", "Normal"):
                if bulk.has(g):
                    r = summarise(bulk, g, "NK", cond)
                else:
                    r = summarise(marker, g, "NK", cond,
                                  source="pseudobulk_per_sample "
                                         "(all cells, ambient-subtracted; "
                                         "gene barred from pseudobulk_raw)")
                r["side"] = side
                rows.append(r)
    q6 = pd.DataFrame(rows)
    q6.to_csv(os.path.join(OUT, "p1_q6_measurability.csv"), index=False)
    log(q6[["side", "gene", "condition", "cpm_donor_median", "det_pooled",
            "det_pooled_lo95", "det_pooled_hi95", "band_pooled",
            "source"]].round(4).to_string(index=False))
    unt = sorted(set(q6[(q6.side == "residency") & (q6.band_pooled == "floor")].gene))
    log(f"  residency-side genes in the floor band: {unt or 'none'}")

    # Appendix: the between-donor correlation.  NOT the Q6 statistic --
    # the task forbids exactly this pooling, so it is filed as an appendix
    # and never as an answer.
    log("\n--- APPENDIX (NOT Q6): between-donor Spearman, "
        "confounded by donor as the task states ---")
    rows = []
    for cond in ("Tumor", "Normal"):
        p = bulk.per_donor("PIEZO1", "NK", cond).set_index("patient")["cpm"]
        for side, gs in Q6_GENES.items():
            for g in gs:
                src = bulk if bulk.has(g) else marker
                o = src.per_donor(g, "NK", cond).set_index("patient")["cpm"]
                sh = sorted(set(p.index) & set(o.index))
                rho, pv = stats.spearmanr(p.loc[sh], o.loc[sh])
                rows.append({"condition": cond, "side": side, "gene": g,
                             "n_donors": len(sh), "spearman_rho": rho,
                             "p": pv, "statistic_is_the_Q6_one": False})
    ap = pd.DataFrame(rows)
    ap.to_csv(os.path.join(OUT, "p1_q6_donorlevel_corr_APPENDIX.csv"), index=False)
    log(ap.round(3).to_string(index=False))

    with open(os.path.join(OUT, "p1_piezo1_baseline.log"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    log("\nDONE")


if __name__ == "__main__":
    main()
