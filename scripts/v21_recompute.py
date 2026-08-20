#!/usr/bin/env python3
"""
v2.1 handoff items 2 and 3 on GSE154826.

Item 2 — recompute §6.8's arithmetic from RAW DONOR-LEVEL data.
The hand calculation combined two independent SEs. But the donors are
PAIRED (the same patient supplies tumour and adjacent normal), and a
paired contrast removes between-donor variance. If within-patient
correlation is positive, the paired SE is SMALLER, power is HIGHER, and
the MDPC verdict could flip. Both are computed and reported.

Item 3 — §5.3 rule 4 sensitivity controls: for every gene with any
out-of-band donor, report the main estimate over ALL donors alongside a
version excluding out-of-band donors. The exclusion version is a
SENSITIVITY control only and may never be the main estimate. Direction
mismatch -> UNINTERPRETABLE.

Scope note: this computes a tumour-versus-normal contrast because §6.8
scoring requires it, and the handoff authorises exactly that. It is used
to establish POWER (can a proportional change be detected), not to assert
a mechanism: per §7.5 only observational-layer wording is permitted, and
per §2 the §6 main computation remains gated behind G-KILL-3 ①.
"""
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 210)

OUT = "results"
GENES = ["CCL3", "CCL4", "CCL5", "XCL1"]
FLOOR, CEIL = 0.05, 0.85
MIN_NK = 30


def main():
    L = pd.read_csv(f"{OUT}/gse154826_nk_detection_per_library.csv")
    K = L[L.n_nk >= MIN_NK]
    P = K.groupby(["patient", "condition"])[GENES].mean().reset_index()

    wide = P.pivot(index="patient", columns="condition", values=GENES)
    paired_ids = [p for p in wide.index
                  if not wide.loc[p].isna().any()]
    print("=" * 78)
    print("v2.1 §6.8 RECOMPUTE FROM RAW DONOR-LEVEL DATA — GSE154826")
    print("=" * 78)
    print(f"donors with normal: {(P.condition=='normal').sum()}   "
          f"with tumour: {(P.condition=='tumor').sum()}   "
          f"PAIRED (both): {len(paired_ids)}")

    rows = []
    for g in GENES:
        n = P[P.condition == "normal"].set_index("patient")[g]
        t = P[P.condition == "tumor"].set_index("patient")[g]
        base = n.mean()                      # baseline = normal mean

        # ---- unpaired (the hand calculation) ------------------------
        sd_n, sd_t = n.std(ddof=1), t.std(ddof=1)
        se_n, se_t = sd_n / np.sqrt(len(n)), sd_t / np.sqrt(len(t))
        se_unpaired = np.sqrt(se_n ** 2 + se_t ** 2)
        diff_unpaired = t.mean() - n.mean()

        # ---- paired (correct for this design) -----------------------
        nn = n.loc[paired_ids]
        tt = t.loc[paired_ids]
        d = (tt - nn).values
        npair = len(d)
        sd_d = d.std(ddof=1)
        se_paired = sd_d / np.sqrt(npair)
        diff_paired = d.mean()
        r = np.corrcoef(nn.values, tt.values)[0, 1]
        tstat, pval = stats.ttest_rel(tt.values, nn.values)
        tcrit = stats.t.ppf(0.975, npair - 1)

        prop_paired = diff_paired / base
        # MDPC as the amendment defines it: SE / baseline (1-SE bar)
        mdpc_1se_paired = se_paired / base
        mdpc_1se_unpaired = se_unpaired / base
        # and at a real detection bar
        mdpc_t_paired = tcrit * se_paired / base

        rows.append(dict(gene=g, baseline_normal=base,
                         diff_pp_unpaired=100 * diff_unpaired,
                         se_pp_unpaired=100 * se_unpaired,
                         diff_pp_paired=100 * diff_paired,
                         se_pp_paired=100 * se_paired,
                         sd_of_diff_pp=100 * sd_d,
                         within_patient_r=r, n_pairs=npair,
                         t=tstat, p=pval,
                         prop_change_paired=prop_paired,
                         mdpc_1se_paired=mdpc_1se_paired,
                         mdpc_1se_unpaired=mdpc_1se_unpaired,
                         mdpc_t_paired=mdpc_t_paired))

    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/v21_recompute.csv", index=False)

    print("\n" + "-" * 78)
    print("PAIRED vs UNPAIRED  (pp)")
    print("-" * 78)
    print(f"{'gene':6s} {'base':>6s} {'d(unpair)':>10s} {'SE(unpair)':>11s} "
          f"{'d(paired)':>10s} {'SE(paired)':>11s} {'r':>6s} {'SE ratio':>9s}")
    for _, x in R.iterrows():
        print(f"{x.gene:6s} {x.baseline_normal:6.3f} {x.diff_pp_unpaired:10.2f} "
              f"{x.se_pp_unpaired:11.2f} {x.diff_pp_paired:10.2f} "
              f"{x.se_pp_paired:11.2f} {x.within_patient_r:6.2f} "
              f"{x.se_pp_paired/x.se_pp_unpaired:9.2f}")

    print("\n" + "-" * 78)
    print("§6.8 SCORING, PAIRED (the design's correct analysis)")
    print("-" * 78)
    r1prop = R[R.gene.isin(["CCL3", "CCL4", "CCL5"])]["prop_change_paired"]
    r1med = float(np.median(r1prop))
    print(f"R1 observed proportional change (median of CCL3/4/5): {100*r1med:+.1f}%")
    for _, x in R.iterrows():
        det = abs(x.diff_pp_paired) / x.se_pp_paired
        print(f"\n{x.gene}:")
        print(f"   paired diff = {x.diff_pp_paired:+.2f} pp "
              f"({100*x.prop_change_paired:+.1f}%),  SE = {x.se_pp_paired:.2f} pp"
              f"   -> |d|/SE = {det:.2f}")
        print(f"   paired t({x.n_pairs-1}) = {x.t:.3f},  p = {x.p:.4g}")
        print(f"   MDPC (1-SE bar)  = {100*x.mdpc_1se_paired:.1f}%   "
              f"[unpaired would give {100*x.mdpc_1se_unpaired:.1f}%]")
        print(f"   MDPC (t.975 bar) = {100*x.mdpc_t_paired:.1f}%")
        if x.gene == "XCL1":
            need_pp = abs(r1med) * x.baseline_normal * 100
            print(f"   >> if XCL1 fell by R1's {100*r1med:+.1f}%, that is "
                  f"{need_pp:.2f} pp = {need_pp/x.se_pp_paired:.2f} x SE")
            verdict = ("DETECTABLE" if need_pp / x.se_pp_paired >= 1.96
                       else "MARGINAL" if need_pp / x.se_pp_paired >= 1.0
                       else "NOT DETECTABLE")
            print(f"   >> §6.7 rule 3: MDPC(1SE)={100*x.mdpc_1se_paired:.1f}% vs "
                  f"R1 observed {abs(100*r1med):.1f}%  -> "
                  f"{'SATISFIED' if x.mdpc_1se_paired < abs(r1med) else 'NOT SATISFIED'}")
            print(f"   >> a same-proportion fall would be: {verdict}")

    # ---- item 3: §5.3 rule 4 sensitivity controls -------------------
    print("\n" + "=" * 78)
    print("§5.3 rule 4  SENSITIVITY CONTROLS (exclude out-of-band donors)")
    print("main estimate = ALL donors; exclusion version is sensitivity ONLY")
    print("=" * 78)
    srows = []
    for g in GENES:
        n = P[P.condition == "normal"].set_index("patient")[g]
        t = P[P.condition == "tumor"].set_index("patient")[g]
        oob = sorted(set(n[(n < FLOOR) | (n > CEIL)].index) |
                     set(t[(t < FLOOR) | (t > CEIL)].index))
        keep = [p for p in paired_ids if p not in oob]
        d_all = (t.loc[paired_ids] - n.loc[paired_ids]).mean()
        d_sen = (t.loc[keep] - n.loc[keep]).mean() if keep else np.nan
        agree = (np.sign(d_all) == np.sign(d_sen)) if keep else None
        print(f"\n{g}: out-of-band donors = {oob if oob else 'none'}")
        print(f"   main (all {len(paired_ids)} pairs)      : {100*d_all:+.2f} pp")
        if keep:
            print(f"   sensitivity ({len(keep)} pairs kept) : {100*d_sen:+.2f} pp")
            print(f"   direction agrees: {agree}"
                  + ("" if agree else "   <== **UNINTERPRETABLE**"))
        else:
            print("   sensitivity: no pairs left")
        srows.append(dict(gene=g, n_oob=len(oob), oob=";".join(map(str, oob)),
                          main_pp=100 * d_all,
                          sensitivity_pp=100 * d_sen if keep else np.nan,
                          n_pairs_main=len(paired_ids), n_pairs_sens=len(keep),
                          direction_agrees=agree,
                          verdict="UNINTERPRETABLE" if agree is False else "ok"))
    pd.DataFrame(srows).to_csv(f"{OUT}/v21_sensitivity.csv", index=False)

    print("\n" + "=" * 78)
    print("Wrote v21_recompute.csv, v21_sensitivity.csv")
    print("=" * 78)


if __name__ == "__main__":
    main()
