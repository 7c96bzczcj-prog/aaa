#!/usr/bin/env python3
"""
A-3 items 2 and 3 — re-report E-MTAB-10176 under A3.1, and issue the
A3.3 / A3.4 / A3.5 decomposition verdicts.

A3.1  ceiling genes (>85%) NEVER enter a module mean. They are reported
      individually as named observations, with start and end detection
      rates and a lower-bound note. Module means use 5-85% genes only.
      Consequence: R1 mean = Ccl3 + Ccl4. Ccl5 is listed separately.
A3.3  a subgroup with zero cells on one side makes standardisation across
      that composition UNDEFINED, not merely noisy. Distinct from
      UNIDENTIFIABLE, and unrescuable by any effective-n threshold.
A3.4  identifiability is anchored to the dataset's OWN measured replicate
      noise floor. Where no replicates exist the floor is UNMEASURABLE and
      the criterion cannot be evaluated.
A3.5  subset labels must be built independently of the contrast. Joint
      PCA + graph + leiden over both sides is a construction violation.

Record only; the dataset fails §4 admission.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)

OUT = "results"
R1_ALL = ["Ccl3", "Ccl4", "Ccl5"]
CEIL = 0.85
FLOOR = 0.05


def main():
    D = pd.read_csv(f"{OUT}/a2_gene_bands.csv")
    print("=" * 78)
    print("A3.1 RE-REPORT — E-MTAB-10176 (record only; §4-inadmissible)")
    print("=" * 78)

    def get(lin, hr, g, col):
        r = D[(D.lineage == lin) & (D.hours == hr) & (D.gene == g)]
        return float(r[col].iloc[0])

    print("\n--- R1 module mean, A3.1 compliant (5-85% genes only) ---")
    rows = []
    for lin in ["NK", "CD8T"]:
        for hr in [24, 72]:
            keep, ceil_genes = [], []
            for g in R1_ALL:
                pg, pr = get(lin, hr, g, "green"), get(lin, hr, g, "red")
                if pg > CEIL or pr > CEIL:
                    ceil_genes.append((g, pg, pr))
                elif pg < FLOOR or pr < FLOOR:
                    pass
                else:
                    keep.append(g)
            mg = np.mean([get(lin, hr, g, "green") for g in keep]) if keep else np.nan
            mr = np.mean([get(lin, hr, g, "red") for g in keep]) if keep else np.nan
            status = "OK" if len(keep) >= 2 else f"UNCOMPUTABLE ({len(keep)} gene)"
            print(f"\n{lin} {hr}h")
            print(f"   R1 mean over {keep}  [{status}]")
            print(f"      Green={mg:.4f}  Red={mr:.4f}   {100*(mr-mg):+.1f} pp")
            for g, pg, pr in ceil_genes:
                print(f"   [ceiling, EXCLUDED from mean, reported alone] {g}: "
                      f"{pg:.4f} -> {pr:.4f}  {100*(pr-pg):+.1f} pp  "
                      f"(start {'at/above' if pg > CEIL else 'below'} ceiling; "
                      f"decline is a LOWER BOUND)")
            rows.append(dict(lineage=lin, hours=hr, r1_genes=",".join(keep),
                             r1_status=status, r1_green=mg, r1_red=mr,
                             r1_delta_pp=100 * (mr - mg),
                             ceiling_excluded=";".join(g for g, _, _ in ceil_genes)))
    pd.DataFrame(rows).to_csv(f"{OUT}/a3_r1_reported.csv", index=False)

    # ---- A3.5 construction check ------------------------------------
    print("\n" + "=" * 78)
    print("A3.5 LABEL-CONSTRUCTION COMPLIANCE")
    print("=" * 78)
    a = ad.read_h5ad("data/raw/DW_T_NK.h5ad")
    single_pca = "X_pca" in a.obsm
    single_graph = "connectivities" in a.obsp
    print(f"   one X_pca over all {a.n_obs} cells (both sides): {single_pca}")
    print(f"   one neighbour graph over all cells:              {single_graph}")
    print(f"   leiden run over all cells:                       "
          f"{'leiden' in a.obs.columns}")
    print("\n   allowed constructions per A3.5:")
    print("     1. cluster ONE side, project the other  -> NOT used")
    print("     2. canonical marker panel, no unsupervised clustering -> NOT used")
    print("   >> VERDICT: CONSTRUCTION VIOLATION. Labels were built by joint")
    print("      PCA + graph + leiden across Green and Red, so they are partly")
    print("      DEFINED by the contrast. Per A3.5 retroactive clause, any")
    print("      decomposition on these labels is barred from RESULTS.md.")

    # ---- A3.3 / A3.4 verdicts ---------------------------------------
    print("\n" + "=" * 78)
    print("A3.3 / A3.4 DECOMPOSITION VERDICTS")
    print("=" * 78)
    I = pd.read_csv(f"{OUT}/a3_label_independence.csv")
    obs = a.obs
    verdicts = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            s = obs[(obs.main_celltype == lin) & (obs.hours == hr)]
            tab = pd.crosstab(s["celltype_cluster"], s["colour"])
            tab = tab.loc[tab.sum(axis=1) > 0]
            zero_side = [(str(cl), int(tab.loc[cl, "Green"]), int(tab.loc[cl, "Red"]))
                         for cl in tab.index
                         if tab.loc[cl, "Green"] == 0 or tab.loc[cl, "Red"] == 0]
            print(f"\n--- {lin} {hr}h ---")
            if zero_side:
                for cl, ng, nr in zero_side:
                    which = "GREEN" if ng == 0 else "RED"
                    print(f"   A3.3: subgroup '{cl}' has ZERO cells on the {which} "
                          f"side (green={ng}, red={nr})")
                print(f"   >> **UNDEFINED** — standardisation across this composition")
                print(f"      has an infinite weight. Not rescuable by any effective-n")
                print(f"      threshold. Dropping the subgroup is barred (A3.3).")
                v = "UNDEFINED"
            else:
                print("   A3.3: no zero-cell subgroup -> not UNDEFINED")
                v = "n/a"
            print(f"   A3.4: replicate-noise floor for this dataset is "
                  f"UNMEASURABLE (no same-condition replicate libraries),")
            print(f"         so the standardised estimator's SE cannot be compared "
                  f"to it -> **UNIDENTIFIABLE** (criterion not evaluable).")
            verdicts.append(dict(lineage=lin, hours=hr,
                                 a3_3=v,
                                 zero_side_subgroups=";".join(c for c, _, _ in zero_side),
                                 a3_4="UNIDENTIFIABLE (noise floor unmeasurable)",
                                 a3_5="CONSTRUCTION VIOLATION (joint clustering)"))
    pd.DataFrame(verdicts).to_csv(f"{OUT}/a3_decomposition_verdicts.csv", index=False)

    print("\n" + "=" * 78)
    print("NET: the Kitagawa decomposition is NOT run on this dataset, for three")
    print("independent reasons — UNDEFINED support (NK), UNIDENTIFIABLE (no")
    print("measurable floor), and a label-CONSTRUCTION violation (all four cells).")
    print("=" * 78)


if __name__ == "__main__":
    main()
