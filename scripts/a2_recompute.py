#!/usr/bin/env python3
"""
Amendment A-2 rev 2 — recompute on E-MTAB-10176 under the new definitions.

RECORD ONLY. This dataset FAILS §4 admission (n = 1 per side, colour
confounded with library, no same-condition replicates). Per the handoff these
numbers go to EXCLUDED.md and must NOT enter RESULTS.md.

Implements:
  A2.1  R1 = Ccl3/Ccl4/Ccl5 ; R2 = Xcl1 (single gene) ; K = Prf1/Gzma/Gzmb/Gzmk/Nkg7
        main readout R2/R1 ; (R1+R2)/K withdrawn
        rule 1  module means only over §5-passing genes
        rule 2  multi-gene module with <2 passing genes -> UNCOMPUTABLE
        rule 2' single-gene module (R2): (i) 5-85% band, (ii) must exceed the
                dataset's OWN measured replicate-noise floor, (iii) raw R2
                detection reported alongside every R2/R1
        rule 4  leave-one-out mandatory; sign flip -> UNSTABLE. R2 is one gene,
                so leave-one-out runs on R1 only.
  A2.2  §5 ceiling band is DIRECTIONAL:
          <5%   floor    -> unusable either direction
          5-85%          -> measurable
          >85%  ceiling  -> upward NOT measurable; downward MEASURABLE, reported
                            as a LOWER BOUND with the starting detection rate
  A2.4  NK-side subcluster composition reported alongside CD8T's

Units: detection rate from raw.X (log-normalised, zero-preserving).
adata.X is scaled/centred and zero-destroying; never used.
"""
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

warnings.filterwarnings("ignore")
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 60)

H5AD = "data/raw/DW_T_NK.h5ad"
OUT = "results"

R1 = ["Ccl3", "Ccl4", "Ccl5"]
R2 = ["Xcl1"]
K = ["Prf1", "Gzma", "Gzmb", "Gzmk", "Nkg7"]
FLOOR, CEIL = 0.05, 0.85


def classify(pg, pr):
    """A2.2 directional band verdict for a Green->Red comparison."""
    if np.isnan(pg) or np.isnan(pr):
        return "ABSENT"
    if pg < FLOOR or pr < FLOOR:
        return "FLOOR_UNMEASURED"
    if pg <= CEIL and pr <= CEIL:
        return "MEASURABLE"
    # some compartment sits above the ceiling
    if pr < pg:
        return "CEILING_DOWN_LOWERBOUND"   # A2.2: downward IS measurable
    return "CEILING_UP_UNMEASURED"         # A2.2: upward is NOT


USABLE = {"MEASURABLE", "CEILING_DOWN_LOWERBOUND"}


def main():
    a = ad.read_h5ad(H5AD)
    R = a.raw.X
    R = R.tocsr() if sp.issparse(R) else sp.csr_matrix(np.asarray(R))
    genes = pd.Index(a.raw.var_names)
    obs = a.obs

    print("=" * 78)
    print("A-2 rev 2 RECOMPUTE — E-MTAB-10176")
    print("RECORD ONLY: dataset fails §4 admission (n=1/side). -> EXCLUDED.md")
    print("=" * 78)

    gi = {g: genes.get_loc(g) for g in R1 + R2 + K if g in set(genes)}
    missing = [g for g in R1 + R2 + K if g not in gi]
    print(f"genes absent from matrix: {missing if missing else 'none'}")

    def det(lin, hr, col, g):
        if g not in gi:
            return np.nan
        m = ((obs.main_celltype == lin) & (obs.hours == hr) & (obs.colour == col)).values
        if m.sum() == 0:
            return np.nan
        v = np.asarray(R[:, gi[g]][m].todense()).ravel()
        return float((v > 0).mean())

    # ---- A2.4 composition, NK ALONGSIDE CD8T -----------------------------
    print("\n" + "-" * 78)
    print("A2.4  SUBCLUSTER COMPOSITION, Green -> Red  (NK reported alongside CD8T)")
    print("-" * 78)
    comp_rows = []
    for lin in ["NK", "CD8T"]:
        sub = obs[obs.main_celltype == lin]
        for hr in ["24", "72"]:
            s = sub[sub.hours == hr]
            tab = pd.crosstab(s["celltype_cluster"], s["colour"], normalize="columns")
            tab = tab.loc[tab.sum(axis=1) > 0]
            if "Green" not in tab or "Red" not in tab:
                continue
            tvd = 0.5 * (tab["Red"] - tab["Green"]).abs().sum()
            print(f"\n--- {lin} {hr}h ---")
            for cl in tab.index:
                g, r = tab.loc[cl, "Green"], tab.loc[cl, "Red"]
                print(f"   {str(cl):22s} {g:6.3f} -> {r:6.3f}   {100*(r-g):+6.1f} pp")
                comp_rows.append(dict(lineage=lin, hours=hr, subcluster=str(cl),
                                      frac_green=g, frac_red=r, delta_pp=100 * (r - g)))
            print(f"   >> total variation distance (0=identical, 1=disjoint): {tvd:.3f}")
            comp_rows.append(dict(lineage=lin, hours=hr, subcluster="__TVD__",
                                  frac_green=np.nan, frac_red=np.nan, delta_pp=tvd))
    pd.DataFrame(comp_rows).to_csv(f"{OUT}/a2_composition.csv", index=False)

    # ---- per-gene detection + A2.2 directional band ----------------------
    print("\n" + "-" * 78)
    print("A2.2  PER-GENE DETECTION AND DIRECTIONAL BAND VERDICT")
    print("-" * 78)
    rows = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            print(f"\n--- {lin} {hr}h ---")
            print(f"   {'gene':7s} {'mod':4s} {'Green':>7s} {'Red':>7s} {'delta_pp':>9s}  verdict")
            for mod, gl in (("R1", R1), ("R2", R2), ("K", K)):
                for g in gl:
                    pg, pr = det(lin, hr, "Green", g), det(lin, hr, "Red", g)
                    v = classify(pg, pr)
                    rows.append(dict(lineage=lin, hours=hr, gene=g, module=mod,
                                     green=pg, red=pr,
                                     delta_pp=100 * (pr - pg) if not np.isnan(pg) else np.nan,
                                     verdict=v, usable=v in USABLE))
                    print(f"   {g:7s} {mod:4s} {pg:7.4f} {pr:7.4f} {100*(pr-pg):+9.1f}  {v}")
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/a2_gene_bands.csv", index=False)

    # ---- module means and R2/R1 ------------------------------------------
    print("\n" + "-" * 78)
    print("A2.1  MODULE VALUES AND MAIN READOUT R2/R1")
    print("A2.1 rule 2'(iii): raw R2 detection reported alongside every ratio")
    print("-" * 78)

    def modmean(lin, hr, col, gl, drop=None):
        sub = D[(D.lineage == lin) & (D.hours == hr) & (D.module.isin(["R1", "R2", "K"]))]
        keep = [g for g in gl if g != drop and
                sub[(sub.gene == g)]["usable"].iloc[0]]
        if not keep:
            return np.nan, []
        vals = [det(lin, hr, col, g) for g in keep]
        return float(np.mean(vals)), keep

    res = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            r1g, r1keep = modmean(lin, hr, "Green", R1)
            r1r, _ = modmean(lin, hr, "Red", R1)
            r2g, r2keep = modmean(lin, hr, "Green", R2)
            r2r, _ = modmean(lin, hr, "Red", R2)
            kg, kkeep = modmean(lin, hr, "Green", K)
            kr, _ = modmean(lin, hr, "Red", K)

            # rule 2 : multi-gene module needs >=2 passing genes
            r1_status = "OK" if len(r1keep) >= 2 else f"UNCOMPUTABLE (only {len(r1keep)} gene passes)"
            k_status = "OK" if len(kkeep) >= 2 else f"UNCOMPUTABLE (only {len(kkeep)} gene passes)"
            # rule 2' : R2 single gene
            xg, xr = det(lin, hr, "Green", "Xcl1"), det(lin, hr, "Red", "Xcl1")
            band_ok = (FLOOR <= xg <= CEIL) and (FLOOR <= xr <= CEIL)
            r2_status = "band OK" if band_ok else "band FAIL"

            ratio_g = r2g / r1g if (len(r1keep) >= 2 and r1g and not np.isnan(r2g)) else np.nan
            ratio_r = r2r / r1r if (len(r1keep) >= 2 and r1r and not np.isnan(r2r)) else np.nan
            delta = ratio_r - ratio_g

            print(f"\n{lin} {hr}h")
            print(f"   R1 genes used {r1keep}  -> {r1_status}")
            print(f"      R1  Green={r1g:.4f}  Red={r1r:.4f}   {100*(r1r-r1g):+.1f} pp")
            print(f"   R2 (Xcl1, single gene) -> rule 2'(i) {r2_status}")
            print(f"      R2  Green={xg:.4f}  Red={xr:.4f}   {100*(xr-xg):+.1f} pp   [rule 2'(iii) raw]")
            print(f"   K  genes used {kkeep}  -> {k_status}")
            if kkeep:
                print(f"      K   Green={kg:.4f}  Red={kr:.4f}   {100*(kr-kg):+.1f} pp")
            print(f"   >> R2/R1  Green={ratio_g:.4f}  Red={ratio_r:.4f}  delta={delta:+.4f}")

            res.append(dict(lineage=lin, hours=hr, r1_genes=",".join(r1keep),
                            r1_status=r1_status, r1_green=r1g, r1_red=r1r,
                            r2_green=xg, r2_red=xr, r2_status=r2_status,
                            k_genes=",".join(kkeep), k_status=k_status,
                            k_green=kg, k_red=kr,
                            ratio_green=ratio_g, ratio_red=ratio_r, delta=delta))
    Rr = pd.DataFrame(res)
    Rr.to_csv(f"{OUT}/a2_module_values.csv", index=False)

    # ---- rule 4 leave-one-out, R1 only -----------------------------------
    print("\n" + "-" * 78)
    print("A2.1 rule 4  LEAVE-ONE-OUT on R1 only (R2 is a single gene)")
    print("-" * 78)
    lo = []
    for lin in ["NK", "CD8T"]:
        for hr in ["24", "72"]:
            base = Rr[(Rr.lineage == lin) & (Rr.hours == hr)]["delta"].iloc[0]
            if np.isnan(base):
                continue
            print(f"\n{lin} {hr}h baseline delta={base:+.4f}")
            _, keep = modmean(lin, hr, "Green", R1)
            for g in keep:
                r1g, k2 = modmean(lin, hr, "Green", R1, drop=g)
                r1r, _ = modmean(lin, hr, "Red", R1, drop=g)
                if len(k2) < 2:
                    print(f"   drop {g:6s} -> R1 UNCOMPUTABLE (<2 genes; rule 2)")
                    lo.append(dict(lineage=lin, hours=hr, dropped=g, delta=np.nan, flips=None))
                    continue
                d = det(lin, hr, "Red", "Xcl1") / r1r - det(lin, hr, "Green", "Xcl1") / r1g
                flip = np.sign(d) != np.sign(base)
                print(f"   drop {g:6s} -> delta={d:+.4f}  {'** FLIPS -> UNSTABLE **' if flip else 'stable'}")
                lo.append(dict(lineage=lin, hours=hr, dropped=g, delta=d, flips=bool(flip)))
    pd.DataFrame(lo).to_csv(f"{OUT}/a2_leave_one_out.csv", index=False)

    # ---- rule 2'(ii) / A2.3(g) noise floor -------------------------------
    print("\n" + "-" * 78)
    print("A2.1 rule 2'(ii) + A2.3(g)  REPLICATE-NOISE FLOOR")
    print("-" * 78)
    libs = obs.groupby(["colour", "hours"], observed=True)["sample"].nunique()
    print(libs)
    print(f"\nmax libraries in any one (colour,hours) cell: {libs.max()}")
    if libs.max() < 2:
        print("-> noise floor UNMEASURABLE in this dataset (A2.3(g) fires)")
        print("-> rule 2'(ii) CANNOT be satisfied: R2's change is not qualifiable")
        print("-> per A2.3(g) this dataset yields NO difference conclusion at all")

    # ---- hazard flag: selection on realised outcome ----------------------
    print("\n" + "-" * 78)
    print("HAZARD FLAG — A2.2 admission can select on the realised outcome")
    print("-" * 78)
    ceil_genes = D[D.verdict.str.startswith("CEILING")]
    if len(ceil_genes):
        print("Ceiling-band genes, and whether A2.2 admits them (admission depends on")
        print("the DIRECTION each gene happened to move):")
        for _, r in ceil_genes.iterrows():
            print(f"   {r.lineage:5s} {r.hours:>3s}h {r.gene:6s} "
                  f"{r.green:.4f}->{r.red:.4f} {r.delta_pp:+6.1f} pp  "
                  f"{'ADMITTED' if r.usable else 'EXCLUDED'}  ({r.verdict})")
        nadm = int(ceil_genes.usable.sum())
        nexc = int((~ceil_genes.usable).sum())
        print(f"\n   admitted {nadm}, excluded {nexc} — all excluded ones were "
              f"non-falling by construction.")
        print("   => a module mean over ceiling genes is biased DOWNWARD by construction,")
        print("      because only the fallers can enter. Affects K, not R1/R2 here.")

    print("\n" + "=" * 78)
    print("Wrote a2_composition.csv, a2_gene_bands.csv, a2_module_values.csv,")
    print("      a2_leave_one_out.csv")
    print("=" * 78)


if __name__ == "__main__":
    main()
