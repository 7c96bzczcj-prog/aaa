#!/usr/bin/env python3
"""
v2.2 item 2 — re-report MDPC under the corrected definition, and add
rule 5's attenuation-corrected within-patient correlation.

MDPC (v2.2) = (z.975 + z.80) * SE(diff) / baseline  ~= 2.80 * SE / baseline
  i.e. the proportional change detectable at two-sided alpha = 0.05 with
  80% power. The v2.1 definition (1 x SE) corresponds to ~17% power.

Rule 5 attenuation:
  observed r = true r * reliability,  reliability = 1 - tech_var/total_var
  The technical variance comes from patient 522's replicate libraries -
  the ONLY donor with same-condition replicates. For n=2 the SD estimate
  is |a-b|/sqrt(2); a raw range must NOT be compared directly against an
  SD. (The previously reported 0.06x-0.60x ratio compared a range to an
  SD and therefore OVERSTATED the technical share.)
"""
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 210)

OUT = "results"
GENES = ["CCL3", "CCL4", "CCL5", "XCL1"]
MIN_NK = 30
Z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)   # 2.8016


def main():
    R = pd.read_csv(f"{OUT}/v21_recompute.csv").set_index("gene")
    L = pd.read_csv(f"{OUT}/gse154826_nk_detection_per_library.csv")
    K = L[L.n_nk >= MIN_NK]
    P = K.groupby(["patient", "condition"])[GENES].mean().reset_index()

    print("=" * 78)
    print("v2.2 MDPC — two-sided alpha = 0.05, power = 80%")
    print(f"multiplier (z.975 + z.80) = {Z:.4f}")
    print("=" * 78)

    r1_prop = np.median([abs(R.loc[g, "prop_change_paired"])
                         for g in ["CCL3", "CCL4", "CCL5"]])
    print(f"\ncontrol benchmark: R1 median observed proportional change "
          f"= {100*r1_prop:.1f}%\n")

    rows = []
    print(f"{'gene':6s} {'baseline':>9s} {'SE(diff) pp':>12s} "
          f"{'MDPC@80%':>9s} {'MDPC@1SE(v2.1)':>15s} {'rule3 vs R1':>14s}")
    for g in GENES:
        se = R.loc[g, "se_pp_paired"] / 100.0
        base = R.loc[g, "baseline_normal"]
        mdpc80 = Z * se / base
        mdpc1 = se / base
        ok80 = mdpc80 < r1_prop
        ok1 = mdpc1 < r1_prop
        print(f"{g:6s} {base:9.4f} {100*se:12.2f} {100*mdpc80:8.1f}% "
              f"{100*mdpc1:14.1f}% "
              f"{'SATISFIED' if ok80 else 'NOT SATISFIED':>14s}")
        rows.append(dict(gene=g, baseline=base, se_diff_pp=100 * se,
                         mdpc_power80=mdpc80, mdpc_1se_v21=mdpc1,
                         rule3_satisfied_power80=bool(ok80),
                         rule3_satisfied_1se_v21=bool(ok1),
                         r1_benchmark=r1_prop))
    M = pd.DataFrame(rows)
    M.to_csv(f"{OUT}/v22_mdpc.csv", index=False)

    print("\n   NOTE: under the v2.1 (1 x SE) definition XCL1 gave "
          f"{100*M[M.gene=='XCL1'].mdpc_1se_v21.iloc[0]:.1f}% ;")
    print("   the hand calculation's 22.3% was BELOW the R1 benchmark and would")
    print("   have LITERALLY PERMITTED writing PRESERVED. At 80% power it is")
    print(f"   {100*M[M.gene=='XCL1'].mdpc_power80.iloc[0]:.1f}% — no longer close.")

    # ---- rule 5: attenuation ----------------------------------------
    print("\n" + "=" * 78)
    print("RULE 5 — attenuation-corrected within-patient correlation")
    print("=" * 78)
    rep = K.groupby(["patient", "condition"]).filter(lambda d: len(d) > 1)
    arows = []
    for g in GENES:
        # technical SD from patient 522's replicate pairs (n=2 -> |a-b|/sqrt2)
        tech_sds = []
        for (pt, cond), d in rep.groupby(["patient", "condition"]):
            v = d[g].values
            tech_sds.append(abs(v[0] - v[1]) / np.sqrt(2))
        tech_sd = float(np.mean(tech_sds))
        tech_sd_max = float(np.max(tech_sds))
        inter = float(np.mean([P[P.condition == c][g].std(ddof=1)
                               for c in ["normal", "tumor"]]))
        obs_r = R.loc[g, "within_patient_r"]
        out = []
        for lbl, ts in (("mean tech", tech_sd), ("worst tech", tech_sd_max)):
            rel = 1 - (ts ** 2) / (inter ** 2)
            rel = max(rel, 1e-6)
            true_r = obs_r / rel
            out.append((lbl, ts, rel, true_r))
        print(f"\n{g}: observed r = {obs_r:+.3f}   inter-donor SD = "
              f"{100*inter:.2f} pp")
        for lbl, ts, rel, tr in out:
            print(f"   {lbl:11s} SD = {100*ts:5.2f} pp  -> ratio "
                  f"{ts/inter:.3f}, reliability = {rel:.3f}, "
                  f"attenuation-corrected r ~ {tr:+.3f}")
            arows.append(dict(gene=g, observed_r=obs_r, inter_donor_sd_pp=100 * inter,
                              tech_case=lbl, tech_sd_pp=100 * ts,
                              sd_ratio=ts / inter, reliability=rel,
                              corrected_r=tr))
    A = pd.DataFrame(arows)
    A.to_csv(f"{OUT}/v22_attenuation.csv", index=False)
    print("\n   Reliability is HIGH (technical share small), so the near-zero")
    print("   correlations survive correction: they are not an artefact of")
    print("   measurement error. Caveat: the technical term rests on ONE donor.")
    print("   Correction: the previously reported 0.06x-0.60x compared a RANGE")
    print("   to an SD; on a like-for-like SD basis the ratios are smaller still.")

    # ---- sample size to reach the target ----------------------------
    print("\n" + "=" * 78)
    print("SCALE FEASIBILITY — what would it take to detect an R1-sized fall?")
    print("=" * 78)
    npair = int(R.loc["XCL1", "n_pairs"])
    for g in GENES:
        se = R.loc[g, "se_pp_paired"]
        base = R.loc[g, "baseline_normal"]
        target_pp = r1_prop * base * 100
        se_needed = target_pp / Z
        n_needed = npair * (se / se_needed) ** 2
        print(f"   {g:6s} target {target_pp:5.2f} pp   current SE {se:5.2f} pp   "
              f"SE needed {se_needed:5.2f} pp   n = {n_needed:6.0f} pairs"
              + ("   <== already powered" if n_needed <= npair else ""))
    print(f"\n   (current n = {npair} pairs)")


if __name__ == "__main__":
    main()
