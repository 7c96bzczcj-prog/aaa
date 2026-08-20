#!/usr/bin/env python3
"""
v2.0 handoff items 2 and 3, on GSE154826.

Item 2 (§5.3) — per-gene x per-condition donor counts in each of the three
bands. NO donor is dropped for being out of band. The table is the
deliverable; it is what makes the band judgement auditable.

Item 3 (§8) — two noise floors, reported side by side:
  technical floor  : spread between same-condition REPLICATE LIBRARIES
                     within one donor. Only patient 522 has any, so this
                     is a single-donor floor and per §8 is a LOWER bound
                     that may not serve as the criterion on its own.
  inter-donor floor: dispersion ACROSS DONORS within one condition. Per
                     §8 this is the criterion, because n is counted by
                     donor (§4d).

Scope: this reports dispersion WITHIN each condition and band counts. It
does NOT compute any tumour-versus-normal difference, ratio, or module
value — those belong to the §6 main computation, which is gated behind
G-KILL-3 ① and ③ (§2, gate order is hard).
"""
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)

OUT = "results"
GENES = ["CCL3", "CCL4", "CCL5", "XCL1"]
FLOOR, CEIL = 0.05, 0.85
MIN_NK = 30


def band(p):
    return "floor" if p < FLOOR else ("ceiling" if p > CEIL else "measurable")


def main():
    L = pd.read_csv(f"{OUT}/gse154826_nk_detection_per_library.csv")

    print("=" * 78)
    print("v2.0 §5.3 + §8 on GSE154826")
    print("=" * 78)

    # donor-level values; a donor with >1 library per condition is averaged
    def donor_table(df):
        return df.groupby(["patient", "condition"])[GENES].mean().reset_index()

    for label, sub in (("NK >= 30 per library (as used for (f))", L[L.n_nk >= MIN_NK]),
                       ("ALL libraries, no cell-count filter", L)):
        print("\n" + "#" * 78)
        print(f"# {label}")
        print("#" * 78)
        P = donor_table(sub)

        # ---- §5.3 transparency table --------------------------------
        print("\n" + "-" * 78)
        print("§5.3  BAND MEMBERSHIP BY DONOR  (no donor is dropped)")
        print("-" * 78)
        rows = []
        for cond in ["normal", "tumor"]:
            s = P[P.condition == cond]
            print(f"\n--- {cond}  ({len(s)} donors) ---")
            print(f"   {'gene':6s} {'floor':>6s} {'meas':>6s} {'ceil':>6s}   "
                  f"{'median':>8s} {'min':>8s} {'max':>8s}   rule-5 (>half floor?)")
            for g in GENES:
                b = s[g].map(band)
                nf = int((b == "floor").sum())
                nm = int((b == "measurable").sum())
                nc = int((b == "ceiling").sum())
                unusable = nf > len(s) / 2
                print(f"   {g:6s} {nf:6d} {nm:6d} {nc:6d}   "
                      f"{s[g].median():8.4f} {s[g].min():8.4f} {s[g].max():8.4f}   "
                      f"{'UNUSABLE' if unusable else 'ok'}")
                rows.append(dict(filter=label, condition=cond, gene=g,
                                 n_donors=len(s), n_floor=nf, n_measurable=nm,
                                 n_ceiling=nc, median=s[g].median(),
                                 min=s[g].min(), max=s[g].max(),
                                 rule5_unusable=bool(unusable)))
        pd.DataFrame(rows).to_csv(
            f"{OUT}/gse154826_band_table_{'filtered' if 'NK >=' in label else 'all'}.csv",
            index=False)

        if "NK >=" not in label:
            continue

        # ---- §5.3 rule 4: bound direction ---------------------------
        print("\n" + "-" * 78)
        print("§5.3 rule 4  BOUND DIRECTION for cross-band genes")
        print("-" * 78)
        for g in GENES:
            n = P[P.condition == "normal"][g]
            t = P[P.condition == "tumor"][g]
            nc_n = int((n > CEIL).sum())
            nc_t = int((t > CEIL).sum())
            nf_n = int((n < FLOOR).sum())
            nf_t = int((t < FLOOR).sum())
            hi = "normal" if n.median() > t.median() else "tumor"
            lo = "tumor" if hi == "normal" else "normal"
            note = []
            if nc_n or nc_t:
                where = "normal" if nc_n >= nc_t else "tumor"
                note.append(f"ceiling donors ({nc_n} normal / {nc_t} tumor) sit mostly on "
                            f"{where}; higher-valued side is {hi} -> "
                            + ("difference is a LOWER BOUND" if where == hi
                               else "ceiling is on the LOWER side; rule 4 does not cover this"))
            if nf_n or nf_t:
                where = "normal" if nf_n >= nf_t else "tumor"
                note.append(f"floor donors ({nf_n} normal / {nf_t} tumor) sit mostly on "
                            f"{where}; lower-valued side is {lo}")
            print(f"\n   {g}: " + ("; ".join(note) if note else "no cross-band donors"))
        print("\n   NOTE: 'concentrated' is not given a threshold in §5.3 rule 4.")
        print("   The counts are reported as-is; no difference is claimed here.")

        # ---- §8 inter-donor floor -----------------------------------
        print("\n" + "-" * 78)
        print("§8  INTER-DONOR FLOOR  (dispersion ACROSS donors, within a condition)")
        print("      -- this is the CRITERION per §8 --")
        print("-" * 78)
        fl = []
        for cond in ["normal", "tumor"]:
            s = P[P.condition == cond]
            print(f"\n--- {cond} (n = {len(s)} donors) ---")
            print(f"   {'gene':6s} {'SD':>8s} {'IQR':>8s} {'SE(mean)':>9s} "
                  f"{'CV':>7s}")
            for g in GENES:
                v = s[g].dropna()
                sd = v.std(ddof=1)
                iqr = v.quantile(.75) - v.quantile(.25)
                se = sd / np.sqrt(len(v))
                cv = sd / v.mean() if v.mean() else np.nan
                print(f"   {g:6s} {sd:8.4f} {iqr:8.4f} {se:9.4f} {cv:7.3f}")
                fl.append(dict(level="inter_donor", condition=cond, gene=g,
                               n=len(v), sd=sd, iqr=iqr, se_mean=se, cv=cv))

        # ---- §8 technical floor (patient 522 only) ------------------
        print("\n" + "-" * 78)
        print("§8  TECHNICAL FLOOR  (replicate libraries within one donor)")
        print("      -- LOWER BOUND, single donor, NOT the criterion per §8 --")
        print("-" * 78)
        rep = sub.groupby(["patient", "condition"]).filter(lambda d: len(d) > 1)
        if rep.empty:
            print("   no same-condition replicate libraries")
        else:
            for (pt, cond), d in rep.groupby(["patient", "condition"]):
                print(f"\n--- patient {pt}, {cond}: {len(d)} libraries "
                      f"(batches {sorted(d.batch.tolist())}) ---")
                for g in GENES:
                    vals = d[g].values
                    spread = float(np.max(vals) - np.min(vals))
                    print(f"   {g:6s} " + "  ".join(f"{v:.4f}" for v in vals) +
                          f"   |spread| = {100*spread:5.1f} pp")
                    fl.append(dict(level="technical_single_donor", condition=cond,
                                   gene=g, n=len(d), sd=np.std(vals, ddof=1),
                                   iqr=np.nan, se_mean=np.nan, cv=np.nan,
                                   spread_pp=100 * spread, patient=pt))
        F = pd.DataFrame(fl)
        F.to_csv(f"{OUT}/gse154826_noise_floors.csv", index=False)

        # ---- the comparison §8 demands ------------------------------
        print("\n" + "-" * 78)
        print("§8  SIDE BY SIDE — is the single-donor technical floor a LOWER bound?")
        print("-" * 78)
        print(f"   {'gene':6s} {'cond':7s} {'technical (pp)':>15s} "
              f"{'inter-donor SD (pp)':>21s}  ratio")
        for cond in ["normal", "tumor"]:
            for g in GENES:
                tech = F[(F.level == "technical_single_donor") &
                         (F.condition == cond) & (F.gene == g)]
                inter = F[(F.level == "inter_donor") &
                          (F.condition == cond) & (F.gene == g)]
                if tech.empty or inter.empty:
                    continue
                tp = float(tech.spread_pp.iloc[0])
                ip = float(inter.sd.iloc[0]) * 100
                print(f"   {g:6s} {cond:7s} {tp:15.1f} {ip:21.1f}  "
                      f"{tp/ip:.2f}x" + ("  <== technical SMALLER, as §8 predicts"
                                          if tp < ip else "  <== technical LARGER"))
        print("\n   §8 predicts the single-donor technical floor understates true")
        print("   technical variation, so it should not be used as the criterion.")

    print("\n" + "=" * 78)
    print("NOTE: no tumour-versus-normal difference, ratio or module value was")
    print("computed. §6 main computation is gated behind G-KILL-3 ① and ③.")
    print("=" * 78)


if __name__ == "__main__":
    main()
