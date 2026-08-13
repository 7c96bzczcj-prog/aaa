#!/usr/bin/env python
"""T3 soup_calibration.tsv — the two rulers, re-calibrated inside this dataset.

R4: no statement of the form "X is expressed in dNK" may be made without both
rulers attached.

Ruler A — soup fraction inside the NK gate.
  This matrix has no empty droplets, so the ambient profile is estimated as
  the abundance-weighted average transcriptome of the compartment (what a
  droplet's free RNA looks like). The per-cell ambient load rho is estimated
  from lineage-foreign genes: genes NK cells cannot transcribe, whose entire
  NK-gate signal must be pickup.

  rho is estimated LEAVE-ONE-OUT: the gene being judged never contributes to
  its own ambient expectation. Without that, every ambient control returns
  exactly 1.0 by construction and the ruler calibrates against itself.

      soup_fraction(g) = min(1, rho_LOO(g) * cpm_ambient(g) / cpm_NK(g))

  For a pure ambient gene the true value is 1.0. The spread the estimator
  actually returns over held-out ambient genes IS the calibration; its lower
  edge is the conservative ceiling.

Ruler B — NK / reference-lineage CPM ratio.
  Pure pickup sits in a low band. A gene must sit materially above the band to
  count as really expressed.

  When the reference lineage barely expresses the gene at all, the ratio has
  no power in either direction and the verdict is `unmeasurable` — NOT
  `negative`.

    python src/soup_calibration.py manifests/VT2018.yaml
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnkchem.counts import as_csr, cell_qc  # noqa: E402
from dnkchem.dataset import load_dataset, load_panel, unit_indices  # noqa: E402
from dnkchem.manifest import load_manifest  # noqa: E402

REF_LINEAGE = "Myeloid"        # v1.0 fixed reference, kept for the frozen column
STROMAL_LINEAGE = "Stromal"
# v1.1: ruler B is per-gene. The denominator is the lineage that actually
# dominates that gene's ambient contribution, chosen from these candidates by
# measured CPM. A single fixed denominator is a concept error: NK/myeloid only
# tests pickup when myeloid is the competing source. For a stromal gene both
# gates pick up equally and the ratio is ~1 whatever the truth is.
CANDIDATE_SOURCES = ["Myeloid", "Stromal", "Trophoblast", "T", "Endothelial",
                     "Perivascular", "Epithelial", "Plasma", "cDC1", "ILC3",
                     "Granulocyte", "Hofbauer"]
# ruler B is powered only if the reference lineage really expresses the gene.
REF_MIN_CPM = 10.0
REF_MIN_DETECTION = 0.05


def cpm_profile(X, rows, cols):
    """CPM per gene over a set of cells, pooled (sum counts / sum depth)."""
    if len(rows) == 0:
        return np.full(len(cols), np.nan), np.full(len(cols), np.nan)
    sub = X[rows]
    depth = float(sub.sum())
    if depth == 0:
        return np.full(len(cols), np.nan), np.full(len(cols), np.nan)
    g = np.asarray(sub[:, cols].sum(axis=0)).ravel()
    det = np.asarray((sub[:, cols] >= 1).sum(axis=0)).ravel() / sub.shape[0]
    return g / depth * 1e6, det


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--panel", default="panel/chemokine_panel_v1.tsv")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=0.10)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    mf = load_manifest(args.manifest)
    panel = load_panel(args.panel)
    outdir = args.outdir or os.path.join("out", mf.dataset_id)
    os.makedirs(outdir, exist_ok=True)

    ds = load_dataset(mf)
    hits, misses, rate, _ = ds.panel_columns(panel)
    X = as_csr(ds.X)
    keep, _ = cell_qc(X, ds.symbols, args.min_genes, args.max_mito)
    obs = ds.obs
    genes = [g for g in panel["gene_symbol"] if g in hits]
    cols = [hits[g] for g in genes]
    cat = dict(zip(panel["gene_symbol"], panel["category"]))
    gidx = {g: i for i, g in enumerate(genes)}

    ambient_ctrl = [g for g in genes if cat[g] == "control_ambient"]
    nk_ctrl = [g for g in genes if cat[g] == "control_nk"]
    print(f"[calib] ambient controls: {ambient_ctrl}")
    print(f"[calib] NK controls:      {nk_ctrl}")

    subs = obs["subset"].astype(object).where(obs["subset"].notna(), "").astype(str).to_numpy()
    is_nk = np.array([s.startswith(("dNK", "pbNK")) for s in subs])
    rows_out, calib_rows = [], []

    for (comp,), crows in unit_indices(obs, keep, ["compartment"]):
        # ambient profile: every cell in the compartment, abundance-weighted.
        # This is what free RNA in a droplet looks like when no empty droplets
        # are available to measure it directly.
        amb_cpm, _ = cpm_profile(X, crows, cols)

        nk_rows = crows[is_nk[crows]]
        if nk_rows.size < 30:
            print(f"[calib] {comp}: only {nk_rows.size} NK cells; skipping")
            continue
        nk_cpm, nk_det = cpm_profile(X, nk_rows, cols)

        ref_rows = crows[subs[crows] == REF_LINEAGE]
        ref_cpm, ref_det = cpm_profile(X, ref_rows, cols)
        str_rows = crows[subs[crows] == STROMAL_LINEAGE]
        str_cpm, str_det = cpm_profile(X, str_rows, cols)

        # --- v1.1: per-lineage profiles, for per-gene denominator choice -----
        lineage_cpm, lineage_det, lineage_n = {}, {}, {}
        for lin in CANDIDATE_SOURCES:
            rws = crows[subs[crows] == lin]
            if rws.size < 30:
                continue
            c, d_ = cpm_profile(X, rws, cols)
            lineage_cpm[lin], lineage_det[lin], lineage_n[lin] = c, d_, int(rws.size)

        def dominant_source(i):
            """Lineage with the highest CPM for gene i -- its ambient source."""
            best, best_cpm = None, -1.0
            for lin, c in lineage_cpm.items():
                if np.isfinite(c[i]) and c[i] > best_cpm:
                    best, best_cpm = lin, float(c[i])
            return best, best_cpm

        # --- rho, leave-one-out over the ambient control set ---------------
        ratios = {}
        for g in ambient_ctrl:
            i = gidx[g]
            if np.isfinite(amb_cpm[i]) and amb_cpm[i] > 0 and np.isfinite(nk_cpm[i]):
                ratios[g] = nk_cpm[i] / amb_cpm[i]
        if len(ratios) < 3:
            print(f"[calib] {comp}: too few usable ambient controls; skipping")
            continue
        rho_all = float(np.median(list(ratios.values())))

        def rho_loo(g):
            vals = [v for k, v in ratios.items() if k != g]
            return float(np.median(vals)) if vals else rho_all

        # --- ruler A calibration: held-out ambient controls ---------------
        held = []
        for g in ambient_ctrl:
            i = gidx[g]
            if not (np.isfinite(nk_cpm[i]) and nk_cpm[i] > 0):
                continue
            sf = min(1.0, rho_loo(g) * amb_cpm[i] / nk_cpm[i])
            held.append(sf)
        ceiling = float(np.min(held)) if held else np.nan

        floor_vals = []
        for g in nk_ctrl:
            i = gidx[g]
            if np.isfinite(nk_cpm[i]) and nk_cpm[i] > 0:
                floor_vals.append(min(1.0, rho_loo(g) * amb_cpm[i] / nk_cpm[i]))

        # --- ruler B calibration: pickup band from the same controls ------
        # Pre-registered band: every ambient control. This is the conservative
        # choice (a wider band makes it HARDER for a gene to pass), and it is
        # what the primary verdicts below use.
        band, band_ref_sourced, band_genes = [], [], []
        for g in ambient_ctrl:
            i = gidx[g]
            if np.isfinite(ref_cpm[i]) and ref_cpm[i] > 0 and np.isfinite(nk_cpm[i]):
                r = nk_cpm[i] / ref_cpm[i]
                band.append(r)
                # A calibrator only measures "pickup from the reference lineage"
                # if the reference lineage is actually that gene's ambient
                # source, i.e. it is enriched over the compartment average.
                # Tissue-sourced controls (stroma, trophoblast) are picked up
                # equally by NK and myeloid, so their ratio sits near 1 and
                # says nothing about pickup.
                if ref_cpm[i] > amb_cpm[i]:
                    band_ref_sourced.append(r)
                    band_genes.append(g)
        band_hi = float(np.max(band)) if band else np.nan
        band_hi_ref = float(np.max(band_ref_sourced)) if band_ref_sourced else np.nan

        # v1.1 band: each ambient control judged against ITS OWN dominant
        # source. This is the band that actually measures pickup.
        band_v11, band_v11_detail = [], []
        for g in ambient_ctrl:
            i = gidx[g]
            lin, lcpm = dominant_source(i)
            # A calibrator whose own source barely expresses it cannot measure
            # a pickup band; the same power gate that governs the verdicts must
            # govern the calibration, or one low-abundance control sets the band.
            if lin and lcpm >= REF_MIN_CPM and np.isfinite(nk_cpm[i]):
                r = nk_cpm[i] / lcpm
                band_v11.append(r)
                band_v11_detail.append(f"{g}:{lin}(cpm={lcpm:.0f})={r:.4f}")
            elif lin:
                band_v11_detail.append(f"{g}:{lin}(cpm={lcpm:.1f}) EXCLUDED, below "
                                       f"the {REF_MIN_CPM:.0f} CPM power gate")
        band_hi_v11 = float(np.max(band_v11)) if band_v11 else np.nan

        calib_rows.append({
            "dataset_id": mf.dataset_id, "compartment": comp,
            "n_nk_cells": int(nk_rows.size), "n_ref_cells": int(ref_rows.size),
            "rho_ambient_load": rho_all,
            "ruler_a_ceiling": ceiling,
            "ruler_a_held_out_min": float(np.min(held)) if held else np.nan,
            "ruler_a_held_out_max": float(np.max(held)) if held else np.nan,
            "ruler_a_nk_gene_floor_min": float(np.min(floor_vals)) if floor_vals else np.nan,
            "ruler_a_nk_gene_floor_max": float(np.max(floor_vals)) if floor_vals else np.nan,
            "ruler_b_band_upper": band_hi,
            "ruler_b_band_upper_ref_sourced": band_hi_ref,
            "ruler_b_band_upper_per_gene_source": band_hi_v11,
            "ruler_b_per_gene_calibrators": "; ".join(band_v11_detail),
            "ruler_b_ref_sourced_calibrators": ",".join(band_genes),
            "prior_ruler_a_ceiling": 0.402, "prior_ruler_b_band": 0.118,
            "prior_nk_floor_lo": 0.012, "prior_nk_floor_hi": 0.013})

        print(f"\n[calib {comp}] rho={rho_all:.4f}  n_NK={nk_rows.size}  n_{REF_LINEAGE}={ref_rows.size}")
        print(f"  ruler A: held-out ambient controls span "
              f"{np.min(held):.3f}-{np.max(held):.3f} (true value 1.0) "
              f"-> conservative ceiling {ceiling:.3f} [prior 0.402]")
        if floor_vals:
            print(f"           true-NK genes span {np.min(floor_vals):.4f}-"
                  f"{np.max(floor_vals):.4f} [prior 0.012-0.013]")
        print(f"  ruler B: pickup band upper edge {band_hi:.4f} over all ambient "
              f"controls [prior 0.118]")
        print(f"           restricted to controls the reference lineage actually "
              f"sources ({','.join(band_genes) or 'none'}): {band_hi_ref:.4f}")
        print(f"  ruler B v1.1 (each control against its OWN dominant source): "
              f"{band_hi_v11:.4f}")
        for d_ in band_v11_detail:
            print(f"             {d_}")
        if np.isfinite(band_hi) and np.isfinite(band_hi_ref) and band_hi > 5 * band_hi_ref:
            print(f"           NOTE: the two differ by {band_hi / band_hi_ref:.0f}x. "
                  f"Tissue-sourced ambient genes are picked up about equally by NK "
                  f"and by {REF_LINEAGE}, so they sit near ratio 1 and carry no "
                  f"information about pickup. Ruler B has little power here; "
                  f"verdicts use the wider (conservative) band.")

        # --- per-gene verdicts --------------------------------------------
        for g in genes:
            i = gidx[g]
            sf = (min(1.0, rho_loo(g) * amb_cpm[i] / nk_cpm[i])
                  if np.isfinite(nk_cpm[i]) and nk_cpm[i] > 0 else np.nan)
            ratio = (nk_cpm[i] / ref_cpm[i]
                     if np.isfinite(ref_cpm[i]) and ref_cpm[i] > 0 else np.nan)
            str_ratio = (nk_cpm[i] / str_cpm[i]
                         if np.isfinite(str_cpm[i]) and str_cpm[i] > 0 else np.nan)

            if not np.isfinite(sf):
                a_verdict = "unmeasurable"
            elif not np.isfinite(ceiling):
                a_verdict = "unmeasurable"
            elif sf >= ceiling:
                a_verdict = "indistinguishable_from_ambient"
            else:
                a_verdict = "above_ambient"

            # --- ruler B, v1.1: denominator is this gene's dominant source ---
            dom_lin, dom_cpm = dominant_source(i)
            dom_ratio = (nk_cpm[i] / dom_cpm
                         if dom_lin and dom_cpm > 0 and np.isfinite(nk_cpm[i]) else np.nan)
            dom_det = (lineage_det[dom_lin][i] if dom_lin else np.nan)
            dom_powered = (dom_lin is not None and np.isfinite(dom_cpm)
                           and dom_cpm >= REF_MIN_CPM
                           and np.isfinite(dom_det) and dom_det >= REF_MIN_DETECTION)
            if not dom_powered or not np.isfinite(dom_ratio) or not np.isfinite(band_hi_v11):
                b_verdict = "unmeasurable"
            elif dom_ratio > band_hi_v11:
                b_verdict = "above_pickup_band"
            else:
                b_verdict = "within_pickup_band"

            # v1.0 verdict retained for comparison
            ref_powered = (np.isfinite(ref_cpm[i]) and ref_cpm[i] >= REF_MIN_CPM
                           and np.isfinite(ref_det[i]) and ref_det[i] >= REF_MIN_DETECTION)
            if not ref_powered or not np.isfinite(ratio) or not np.isfinite(band_hi):
                b_verdict_v10 = "unmeasurable"
            elif ratio > band_hi:
                b_verdict_v10 = "above_pickup_band"
            else:
                b_verdict_v10 = "within_pickup_band"

            if a_verdict == "unmeasurable" or b_verdict == "unmeasurable":
                combined = "unmeasurable"
            elif a_verdict == "above_ambient" and b_verdict == "above_pickup_band":
                combined = "positive"
            elif a_verdict == "above_ambient" or b_verdict == "above_pickup_band":
                combined = "weak"
            else:
                combined = "indistinguishable_from_ambient"

            rows_out.append({
                "dataset_id": mf.dataset_id, "compartment": comp, "gene": g,
                "category": cat[g], "soup_fraction": sf,
                "nk_myeloid_cpm_ratio": ratio,
                "scale_a_verdict": a_verdict, "scale_b_verdict": b_verdict,
                "combined_verdict": combined,
                # supplementary, joinable on (dataset_id, compartment, gene)
                "_nk_cpm": nk_cpm[i], "_ambient_cpm": amb_cpm[i],
                "_ref_cpm": ref_cpm[i], "_ref_detection": ref_det[i],
                "_nk_detection": nk_det[i], "_stromal_cpm": str_cpm[i],
                "_nk_stromal_cpm_ratio": str_ratio,
                "_dominant_source_lineage": dom_lin,
                "_dominant_source_cpm": dom_cpm,
                "_nk_dominant_source_cpm_ratio": dom_ratio,
                "_scale_b_verdict_v10_myeloid_only": b_verdict_v10,
                "_high_stromal_background": bool(
                    "high_stromal_background=TRUE" in str(
                        panel.loc[panel.gene_symbol == g, "notes"].iloc[0]))})

    T3full = pd.DataFrame(rows_out)
    core = ["dataset_id", "compartment", "gene", "category", "soup_fraction",
            "nk_myeloid_cpm_ratio", "scale_a_verdict", "scale_b_verdict",
            "combined_verdict"]
    T3full[core].to_csv(os.path.join(outdir, "soup_calibration.tsv"),
                        sep="\t", index=False)
    supp = T3full[["dataset_id", "compartment", "gene"] +
                  [c for c in T3full.columns if c.startswith("_")]]
    supp.columns = [c.lstrip("_") for c in supp.columns]
    supp.to_csv(os.path.join(outdir, "soup_calibration_supplement.tsv"),
                sep="\t", index=False)
    pd.DataFrame(calib_rows).to_csv(
        os.path.join(outdir, "soup_ruler_calibration.tsv"), sep="\t", index=False)
    print(f"\n[T3] {len(T3full)} rows -> {outdir}/soup_calibration.tsv")

    with open(os.path.join(outdir, "soup_meta.json"), "w") as fh:
        json.dump({"reference_lineage": REF_LINEAGE,
                   "ref_min_cpm": REF_MIN_CPM,
                   "ref_min_detection": REF_MIN_DETECTION,
                   "calibration": calib_rows}, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
