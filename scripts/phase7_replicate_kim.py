"""Phase 7 - replication of the GSE154826 quadrant calls in GSE131907.

This is the protocol's hard gate: a Q4 candidate that does not reproduce
in an independent cohort does not advance, however artefact-proof it
looked in the first dataset.

Two things about this particular replication have to be said up front,
because they bound what a negative result would mean.

1. It is underpowered by construction. GSE131907 reaches the A1 minimum
   of 8 individuals only at a relaxed 20-cell threshold (see D8), so a
   Q4 gene failing here is ambiguous between "the original call was
   false" and "there was never enough power to see it". That is why the
   Q1 and Q3 replication rates are reported alongside Q4 -- stop rule S6
   is a *relative* comparison for exactly this reason.

2. The two datasets are annotated by different means: lineages here come
   from the authors' own published cell-type labels, not from this
   repository's gating. That makes the replication genuinely independent
   of the gating code, at the cost of the two lineage definitions not
   being identical.

Library sizes are computed from the retained gene panel rather than the
full transcriptome, because the source matrix is a 390 MB gzipped text
file and only the panel rows are parsed. Filtering genes before
normalisation is standard practice, but it is an approximation and is
noted here rather than hidden.
"""

from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from nkmine.de import paired_de  # noqa: E402
from nkmine.null_calibration import calibrate  # noqa: E402
from nkmine.pseudobulk import PseudobulkSet, se_balance_report  # noqa: E402
from nkmine.quadrant import LINEAGES, classify_table  # noqa: E402
from nkmine.replication import (  # noqa: E402
    evaluate_stop_rules_s5_s6, replicate, replication_rate_by_quadrant,
)

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
KIM = os.path.join(ROOT, "data", "kim")
MIN_CELLS = 20          # relaxed threshold; see docs/DEVIATIONS.md D8
N_REPS = 5
DELTA = 0.5


def lineage_of(row) -> str | None:
    st, ct = str(row.Cell_subtype), str(row.Cell_type)
    if st == "NK":
        return "NK"
    if st.startswith("CD8"):
        return "CD8T"
    if st.startswith("CD4") or st == "Treg":
        return "CD4T"
    if ct == "B lymphocytes":
        return "B"
    if ct == "Myeloid cells":
        return "Myeloid"
    return None


def build_groups():
    ann = pd.read_csv(os.path.join(KIM, "cellann.txt.gz"), sep="\t")
    ann = ann[ann.Sample_Origin.isin(["tLung", "nLung"])].copy()
    ann["pid"] = ann.Sample.str.extract(r"LUNG_[NT](\d+)")
    ann["lin"] = ann.apply(lineage_of, axis=1)
    ann = ann[ann.lin.notna()]
    ann["cond"] = ann.Sample_Origin.map({"tLung": "Tumor", "nLung": "Normal"})

    n = ann.groupby(["pid", "cond", "lin"]).size().reset_index(name="n")
    mn = n.groupby(["pid", "cond"]).n.min().reset_index(name="m")
    ok = mn[mn.m >= MIN_CELLS]
    both = ok.groupby("pid").cond.nunique()
    keep_pids = sorted(both[both == 2].index)
    ann = ann[ann.pid.isin(keep_pids)]
    print(f"patients usable at >={MIN_CELLS} cells in all five lineages: "
          f"{len(keep_pids)} -> {keep_pids}", flush=True)

    nmin = (ann.groupby(["pid", "cond", "lin"]).size()
               .groupby(level=[0, 1]).min().to_dict())
    return ann, nmin, keep_pids


def main():
    panel = pd.read_csv(os.path.join(OUT, "quadrants_relaxed_d0.5_p95.csv"))
    want = set(panel.gene)
    ann, nmin, pids = build_groups()
    if len(pids) < 8:
        print(f"STOP: only {len(pids)} individuals, below A1 minimum of 8", flush=True)

    # balanced subsample: index -> (rep, group) membership
    rng = np.random.default_rng(0)
    cell_group = {}       # barcode index label -> list of (rep, pid, cond, lin)
    for (pid, cond), m in nmin.items():
        for lin in LINEAGES:
            ids = ann[(ann.pid == pid) & (ann.cond == cond) & (ann.lin == lin)].Index.values
            if len(ids) < m:
                continue
            for r in range(N_REPS):
                for cid in rng.choice(ids, size=m, replace=False):
                    cell_group.setdefault(cid, []).append((r, pid, cond, lin))
    print(f"cells participating in >=1 subsample: {len(cell_group)}", flush=True)

    path = os.path.join(KIM, "raw_umi.txt.gz")
    acc = {}
    lib = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = header[1:]
        col_idx = {c: i for i, c in enumerate(cols)}
        take_pos, take_keys = [], []
        for cid, groups in cell_group.items():
            if cid in col_idx:
                take_pos.append(col_idx[cid])
                take_keys.append(groups)
        take_pos = np.array(take_pos)
        print(f"matched {len(take_pos)} of {len(cell_group)} cells to matrix columns "
              f"({len(cols)} columns total)", flush=True)

        genes_kept = []
        rows_kept = []
        n_lines = 0
        for line in fh:
            n_lines += 1
            tab = line.find("\t")
            g = line[:tab].strip('"')
            if g not in want:
                continue
            vals = np.fromstring(line[tab + 1:], sep="\t")
            if vals.size != len(cols):
                continue
            genes_kept.append(g)
            rows_kept.append(vals[take_pos].astype(np.float32))
            if len(genes_kept) % 500 == 0:
                print(f"  parsed {len(genes_kept)} panel genes "
                      f"({n_lines} lines scanned)", flush=True)
    print(f"panel genes found in GSE131907: {len(genes_kept)} of {len(want)}", flush=True)

    mat = np.vstack(rows_kept)            # genes x selected cells
    genes = np.array(genes_kept)

    keys = sorted({k for gs in take_keys for k in gs})
    key_index = {k: i for i, k in enumerate(keys)}
    sums = np.zeros((mat.shape[0], len(keys)), dtype=np.float64)
    for j, groups in enumerate(take_keys):
        for k in groups:
            sums[:, key_index[k]] += mat[:, j]

    # median across the N_REPS balancing repeats
    final_keys = sorted({k[1:] for k in keys})
    counts = np.column_stack([
        np.median([sums[:, key_index[(r,) + k]] for r in range(N_REPS)
                   if (r,) + k in key_index], axis=0)
        for k in final_keys
    ])
    ps = PseudobulkSet(counts, genes,
                       np.array([k[0] for k in final_keys]),
                       np.array([k[1] for k in final_keys]),
                       np.array([k[2] for k in final_keys]),
                       np.ones(len(final_keys))).complete_pairs()
    print(f"GSE131907 pseudobulk: {ps.counts.shape[0]} genes x "
          f"{ps.counts.shape[1]} samples ({len(set(ps.patient))} patients)", flush=True)
    # persist so the synthetic-recall test can be repeated in THIS cohort's
    # design; recall measured in the discovery cohort says nothing about
    # whether the replication cohort could have detected anything.
    np.savez_compressed(os.path.join(OUT, "pseudobulk_kim.npz"),
                        counts=ps.counts, genes=ps.genes, patient=ps.patient,
                        condition=ps.condition, lineage=ps.lineage)

    de = {}
    for l in LINEAGES:
        s = ps.subset_lineage(l)
        de[l] = paired_de(s.counts, s.genes, s.patient, s.condition, "Tumor")
    rep = se_balance_report(de)
    print("S3:", rep["verdict"], flush=True)

    null = calibrate(ps, case_label="Tumor", n_perm=200, rng=np.random.default_rng(1))
    res = classify_table(de, delta=DELTA, null_stats=null, equiv_key="p95")
    res.to_csv(os.path.join(OUT, "quadrants_GSE131907_d0.5_p95.csv"), index=False)
    print("GSE131907 quadrants:", res.quadrant.value_counts().to_dict(), flush=True)

    rp = replicate(panel, res)
    rp.to_csv(os.path.join(OUT, "replication_GSE154826_vs_GSE131907.csv"), index=False)
    rates = replication_rate_by_quadrant(rp)
    print("\n=== replication by quadrant (GSE154826 -> GSE131907) ===", flush=True)
    print(rates.to_string(index=False), flush=True)
    ev = evaluate_stop_rules_s5_s6(rates)
    print("\nS5/S6:", ev["verdict"], flush=True)
    rates.to_csv(os.path.join(OUT, "replication_rates_by_quadrant.csv"), index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
