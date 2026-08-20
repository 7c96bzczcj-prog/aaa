#!/usr/bin/env python3
"""
A-3 item 4 — the ONLY authorised computation this round.

Determine A2.3 criterion (f) on GSE154826: are CCL3 / CCL4 / CCL5 / XCL1
measurable in NK cells, per §5's bands, in BOTH compartments?

ADMISSION ONLY. This script deliberately does NOT compute R1, R2, any
ratio, or any tumour-versus-normal difference. It reports per-compartment
detection rates and the band verdict, nothing else.

n is counted by DONOR, not by library (§4(d)): the headline detection
rate per compartment is the mean over patients of each patient's
detection rate, so a patient with many cells does not dominate. The
cell-pooled value is reported alongside for transparency only.

Lineage gating is the project's existing RNA gate (gate_lineages_rna),
whose marker genes are all in excluded_genes.txt. None of the four target
genes is a gating marker, so there is no circularity.
"""
import os
import re
import sys
import time
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from nkmine.gating import read_10x_tar, gate_lineages_rna  # noqa: E402

RAW = "data/raw"
OUT = "results"
TARGETS = ["CCL3", "CCL4", "CCL5", "XCL1"]
FLOOR, CEIL = 0.05, 0.85


def parse_condition(inner: str):
    """(patient, condition) from the tarball's internal filename.

    Returns condition in {'tumor','normal','hashed_both','unknown'}.
    Hashed libraries (…TN…) carry BOTH conditions in one emulsion and
    need HTO demux, so they are not used for this admission check.
    """
    m = re.match(r"^(\d+)_patient_([0-9]+)-(.+?)_(?:features|barcodes|matrix)\.tsv|"
                 r"^(\d+)_patient_([0-9]+)-(.+?)_features\.tsv", inner)
    mm = re.match(r"^(\d+)_patient_([0-9]+)-(.+?)_features\.tsv", inner)
    if not mm:
        return None, "unknown", inner
    batch, patient, tag = mm.group(1), mm.group(2), mm.group(3)
    t = tag.lower()
    # hashed tumour+normal in one emulsion
    if re.search(r"\dtn\d?$", t) or t.endswith("tn"):
        return patient, "hashed_both", tag
    if "tumor" in t or "tumour" in t:
        return patient, "tumor", tag
    if "normal" in t:
        return patient, "normal", tag
    if re.search(r"t$", t) and not re.search(r"test$", t):
        return patient, "tumor", tag
    if re.search(r"n$", t):
        return patient, "normal", tag
    return patient, "unknown", tag


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0

    man = pd.read_csv(f"{OUT}/gse154826_manifest_raw.tsv", sep="\t",
                      header=None, names=["batch", "inner"])
    rows = []
    for _, r in man.iterrows():
        p, c, tag = parse_condition(str(r.inner))
        rows.append(dict(batch=int(r.batch), patient=p, condition=c, tag=tag))
    M = pd.DataFrame(rows)
    M.to_csv(f"{OUT}/gse154826_manifest.csv", index=False)

    print("=" * 78)
    print("GSE154826 — criterion (f) admission check  [ADMISSION ONLY]")
    print("=" * 78)
    print("\nlibraries by parsed condition:")
    print(M.condition.value_counts().to_string())

    usable = M[M.condition.isin(["tumor", "normal"])].copy()
    pats = usable.groupby("patient")["condition"].nunique()
    paired = sorted(pats[pats == 2].index)
    print(f"\nunambiguous single-condition libraries: {len(usable)}")
    print(f"patients with BOTH tumour and normal among them: {len(paired)}")
    print(f"  -> {paired}")
    print("hashed (…TN) libraries are excluded here: they need HTO demux, "
          "which is\n  beyond an admission check.")

    todo = usable[usable.patient.isin(paired)].sort_values(["patient", "condition"])
    if limit:
        todo = todo.head(limit)
    print(f"\nlibraries to read: {len(todo)}")

    recs = []
    t0 = time.time()
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        path = f"{RAW}/GSE154826_amp_batch_ID_{r.batch}.tar.gz"
        if not os.path.exists(path):
            print(f"  [{i}/{len(todo)}] batch {r.batch} MISSING, skipped")
            continue
        ts = time.time()
        try:
            lib = read_10x_tar(path)
            lineage, _ = gate_lineages_rna(lib)
        except Exception as e:
            print(f"  [{i}/{len(todo)}] batch {r.batch} FAILED: {type(e).__name__}: {e}")
            continue
        nk = lineage == "NK"
        n_nk = int(nk.sum())
        idx = {g: k for k, g in enumerate(lib.gene_names)}
        rec = dict(batch=int(r.batch), patient=r.patient, condition=r.condition,
                   n_cells=int(lineage.size), n_nk=n_nk)
        for g in TARGETS:
            if g in idx and n_nk > 0:
                v = np.asarray(lib.rna[idx[g], :][:, nk].todense()).ravel()
                rec[g] = float((v > 0).mean())
                rec[f"{g}_npos"] = int((v > 0).sum())
            else:
                rec[g] = np.nan
                rec[f"{g}_npos"] = 0
        recs.append(rec)
        el = time.time() - ts
        print(f"  [{i}/{len(todo)}] batch {r.batch:>3} pt {r.patient} "
              f"{r.condition:6s} cells={rec['n_cells']:6d} NK={n_nk:5d}  "
              f"({el:.1f}s)  " +
              "  ".join(f"{g}={rec[g]:.3f}" if not np.isnan(rec[g]) else f"{g}=NA"
                        for g in TARGETS))
        pd.DataFrame(recs).to_csv(f"{OUT}/gse154826_nk_detection_per_library.csv",
                                  index=False)
    print(f"\ntotal read time {time.time()-t0:.0f}s for {len(recs)} libraries")

    if not recs:
        print("no libraries read; aborting")
        return
    L = pd.DataFrame(recs)

    # ---- donor-level aggregation (§4(d): n by donor, not library) -----
    print("\n" + "-" * 78)
    print("§4(d)  DONOR-LEVEL DETECTION RATE  (per-patient mean, equal donor weight)")
    print("-" * 78)
    MIN_NK = 30
    keep = L[L.n_nk >= MIN_NK]
    print(f"libraries with NK >= {MIN_NK}: {len(keep)} of {len(L)}")
    per_pat = (keep.groupby(["patient", "condition"])[TARGETS]
               .mean().reset_index())
    summary = []
    for cond in ["normal", "tumor"]:
        s = per_pat[per_pat.condition == cond]
        pooled = keep[keep.condition == cond]
        print(f"\n--- {cond}  (n donors = {s.patient.nunique()}, "
              f"n libraries = {len(pooled)}, NK cells = {int(pooled.n_nk.sum())}) ---")
        for g in TARGETS:
            donor_mean = s[g].mean()
            # cell-pooled, for transparency only
            tot_pos = pooled[f"{g}_npos"].sum()
            tot_nk = pooled.n_nk.sum()
            cellpool = tot_pos / tot_nk if tot_nk else np.nan
            band = ("FLOOR" if donor_mean < FLOOR else
                    "CEILING" if donor_mean > CEIL else "MEASURABLE")
            print(f"   {g:5s} donor-mean={donor_mean:.4f}  "
                  f"[cell-pooled {cellpool:.4f}]   band={band}")
            summary.append(dict(condition=cond, gene=g, donor_mean=donor_mean,
                                cell_pooled=cellpool, band=band,
                                n_donors=int(s.patient.nunique())))
    S = pd.DataFrame(summary)
    S.to_csv(f"{OUT}/gse154826_f_verdict.csv", index=False)

    # ---- (f) verdict ---------------------------------------------------
    print("\n" + "=" * 78)
    print("CRITERION (f) VERDICT")
    print("=" * 78)
    r1 = ["CCL3", "CCL4", "CCL5"]
    ok_r1 = [g for g in r1
             if all(S[(S.gene == g) & (S.condition == c)]["band"].iloc[0] == "MEASURABLE"
                    for c in ["normal", "tumor"])]
    x = S[S.gene == "XCL1"]
    ok_r2 = all(x[x.condition == c]["band"].iloc[0] == "MEASURABLE"
                for c in ["normal", "tumor"])
    print(f"R1: genes MEASURABLE in BOTH compartments = {ok_r1}  "
          f"({len(ok_r1)} of 3; need >= 2)  -> {'PASS' if len(ok_r1) >= 2 else 'FAIL'}")
    print(f"R2 (XCL1, rule 2'(i) needs 5-85% in both) -> "
          f"{'PASS' if ok_r2 else 'FAIL'}")
    verdict = "PASS" if (len(ok_r1) >= 2 and ok_r2) else "FAIL"
    print(f"\n>> criterion (f) on GSE154826: **{verdict}**")
    print("\nNOTE: no ratio, no module value, and no tumour-versus-normal")
    print("difference has been computed. Both compartments are shown because")
    print("§5 requires the band to hold in each; that is admission, not comparison.")


if __name__ == "__main__":
    main()
