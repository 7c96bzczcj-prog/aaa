#!/usr/bin/env python3
"""
Cross-organ check §5 aggregation, §7 admission, §8 table, §6 plots.

Statistical unit is the DONOR (§5): each donor gives one detection rate,
then the donor-level median and IQR are reported. Cells are never pooled
across donors.

§7 marks a stratum INSUFFICIENT when donors < 3, NK < 200, one donor
supplies > half the cells, or the gene sits in the < 5% floor band.

§4 forbids merging or ranking across chemistries, so every stratum is
(tissue x assay), never tissue alone.

No cross-organ test, fold change or ranking is computed. §6's four
pre-written interpretation rows are evaluated from the plots' diagnostics.
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

OUT = "results"
PLOTS = "COVARIATE_PLOTS"
TARGETS = ["CCL3", "CCL4", "CCL5", "XCL1", "XCL2"]
MIN_DONORS, MIN_NK, FLOOR = 3, 200, 0.05


def main():
    D = pd.read_csv(f"{OUT}/organ_per_donor.csv")
    os.makedirs(PLOTS, exist_ok=True)
    rows = []

    for (defn, tis, asy), g in D.groupby(["definition", "tissue", "assay"],
                                         observed=True):
        n_don = g.donor.nunique()
        n_nk = int(g.n_nk.sum())
        top_share = g.n_nk.max() / n_nk if n_nk else 1.0
        reasons = []
        if n_don < MIN_DONORS:
            reasons.append(f"donors {n_don}<{MIN_DONORS}")
        if n_nk < MIN_NK:
            reasons.append(f"NK {n_nk}<{MIN_NK}")
        if top_share > 0.5:
            reasons.append(f"one donor {top_share:.0%} of cells")

        rec = dict(definition=defn, tissue=tis, assay=asy, n_donors=n_don,
                   n_nk=n_nk, top_donor_share=top_share,
                   median_umi=float(g.median_umi.median()),
                   bright_frac=float(g.bright_frac.median()),
                   stress=float(g.stress.median()),
                   insufficient="; ".join(reasons) if reasons else "")
        for gene in TARGETS:
            v = g[gene].dropna()
            vb = g[gene + "_bright"].dropna()
            med = float(v.median()) if len(v) else np.nan
            rec[gene] = med
            rec[gene + "_iqr"] = (float(v.quantile(.75) - v.quantile(.25))
                                  if len(v) else np.nan)
            rec[gene + "_bright"] = float(vb.median()) if len(vb) else np.nan
            rec[gene + "_floor"] = bool(med < FLOOR) if not np.isnan(med) else False
        rows.append(rec)

    T = pd.DataFrame(rows).sort_values(["definition", "assay", "tissue"])
    T.to_csv(f"{OUT}/organ_table.csv", index=False)
    ok = T[T.insufficient == ""]
    print(f"strata total {len(T)}, passing §7: {len(ok)}")
    print(ok[["definition", "tissue", "assay", "n_donors", "n_nk",
              "median_umi", "bright_frac", "stress"] + TARGETS].to_string(index=False))

    # ---- §6 plots: detection vs median UMI, and vs CD56bright share ----
    for gene in TARGETS:
        for xvar, xlab, fname in (("median_umi", "median UMI per NK cell", "umi"),
                                  ("bright_frac", "CD56bright share of NK", "bright")):
            fig, ax = plt.subplots(figsize=(7.2, 5.0))
            for defn, mk in (("author", "o"), ("uniform", "^")):
                s = ok[ok.definition == defn]
                for asy, col in zip(sorted(ok.assay.unique()),
                                    ["#2563eb", "#dc2626", "#16a34a"]):
                    ss = s[s.assay == asy]
                    if ss.empty:
                        continue
                    ax.scatter(ss[xvar], ss[gene], marker=mk, s=70, color=col,
                               alpha=.85, edgecolor="k", linewidth=.5,
                               label=f"{defn} · {asy}")
                    for _, r in ss.iterrows():
                        ax.annotate(r.tissue[:14], (r[xvar], r[gene]),
                                    fontsize=6.5, alpha=.75,
                                    xytext=(3, 3), textcoords="offset points")
            ax.axhline(FLOOR, ls=":", c="grey", lw=1)
            ax.text(ax.get_xlim()[0], FLOOR, " 5% floor", fontsize=7,
                    va="bottom", c="grey")
            ax.set_xlabel(xlab)
            ax.set_ylabel(f"{gene} detection rate (donor median)")
            ax.set_title(f"{gene}: detection vs {xlab}\n"
                         "strata passing §7 only; NOT a cross-organ comparison",
                         fontsize=10)
            ax.legend(fontsize=6.5, frameon=False, ncol=2)
            fig.tight_layout()
            fig.savefig(f"{PLOTS}/{gene}_vs_{fname}.png", dpi=150)
            plt.close(fig)
    print(f"\nwrote plots to {PLOTS}/ ({len(TARGETS)*2} files)")

    # ---- §6 diagnostics feeding the four pre-written rows -------------
    print("\n" + "=" * 78)
    print("§6 DIAGNOSTICS (per chemistry; §4 bars cross-chemistry ranking)")
    print("=" * 78)
    diag = []
    for defn in ["author", "uniform"]:
        for asy in sorted(ok.assay.unique()):
            s = ok[(ok.definition == defn) & (ok.assay == asy)]
            if len(s) < 3:
                print(f"\n{defn} / {asy}: only {len(s)} organ strata -> "
                      "too few to evaluate §6")
                continue
            print(f"\n--- {defn} / {asy}  ({len(s)} organs: "
                  f"{', '.join(sorted(s.tissue))}) ---")
            for gene in TARGETS:
                v = s[gene].dropna()
                vb = s[gene + "_bright"].dropna()
                if len(v) < 3 or v.min() <= 0:
                    print(f"   {gene:5s}  spread not computable "
                          f"(n={len(v)}, min={v.min() if len(v) else np.nan})")
                    continue
                spread = v.max() / v.min()
                spread_b = (vb.max() / vb.min()) if (len(vb) >= 3 and vb.min() > 0) else np.nan
                r_umi = s[[gene, "median_umi"]].corr(method="spearman").iloc[0, 1]
                r_br = s[[gene, "bright_frac"]].corr(method="spearman").iloc[0, 1]
                print(f"   {gene:5s}  spread(max/min) all-NK {spread:5.2f}x  "
                      f"bright-only {spread_b if np.isnan(spread_b) else round(spread_b,2)}  "
                      f"rho(depth) {r_umi:+.2f}  rho(bright%) {r_br:+.2f}"
                      + ("   [FLOOR in >=1 organ]" if (v < FLOOR).any() else ""))
                diag.append(dict(definition=defn, assay=asy, gene=gene,
                                 n_organs=len(s), spread_all=spread,
                                 spread_bright=spread_b, rho_depth=r_umi,
                                 rho_bright=r_br,
                                 any_floor=bool((v < FLOOR).any())))
    pd.DataFrame(diag).to_csv(f"{OUT}/organ_diagnostics.csv", index=False)


if __name__ == "__main__":
    main()
