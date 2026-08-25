#!/usr/bin/env python3
"""Second-stage refinement of the NK / T boundary.

The first pass clusters a donor's whole library at leiden resolution 1.0.  That
is enough when a donor contributes 40,000+ cells, but not when it contributes
3,000: on GSE120221 donor A the entire NK population sat inside a 237-cell
cluster whose T panel outscored its NK panel (1.54 vs 1.08), so the donor
returned zero NK.  NK and cytotoxic CD8 T are transcriptionally adjacent, and at
low cell counts one resolution cannot separate them.

The refinement re-clusters only the lymphoid compartment -- every cell in a
first-pass cluster whose RAW panel argmax was NK or T -- and re-annotates the
sub-clusters with exactly the same rules.  Cells outside that compartment keep
their first-pass label.  Applied uniformly to every dataset and both
compartments, so it cannot introduce an asymmetry between the arms of the
contrast.
"""
from __future__ import annotations

import argparse, json, os, sys, warnings

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.join(HERE, "scripts"))
import pf_core as pf  # noqa: E402
from annotate import build, annotate, LEIDEN_RES, N_PCS, N_NEIGHBORS, HVG_N  # noqa: E402

warnings.filterwarnings("ignore")
LYMPHOID = {"NK", "T"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--files", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    z = np.load(args.labels, allow_pickle=True)
    clusters = json.load(open(args.labels.replace(".npz", "_clusters.json")))
    lymph_clusters = {c["cluster"] for c in clusters if c.get("raw_argmax") in LYMPHOID}
    leiden = z["leiden"].astype(str)
    in_lymph = np.isin(leiden, list(lymph_clusters)) & (z["lineage"] != "doublet")
    if in_lymph.sum() < 200:
        np.savez_compressed(args.out, **{k: z[k] for k in z.files})
        json.dump(clusters, open(args.out.replace(".npz", "_clusters.json"), "w"), indent=1)
        print(f"[refine] {os.path.basename(args.out)}: lymphoid compartment too small "
              f"({int(in_lymph.sum())}), first pass kept", flush=True)
        return

    spec = json.loads(args.files)
    libs = []
    for name, path, kind in spec:
        libs.append((name, pf.read_10x_h5(path) if kind == "h5"
                     else pf.read_10x_mtx(*path.split("|"), name)))
    a = build(libs)
    # build() calls cells in the same order as the first pass, so rows align
    assert list(a.obs["barcode"]) == list(z["barcode"]), "cell order changed"
    assert list(a.obs["library"]) == list(z["library"]), "library order changed"

    import scanpy as sc
    b = a[in_lymph].copy()
    b.layers["counts"] = b.X.copy()
    b.obs["is_doublet"] = False
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=min(HVG_N, b.shape[1] - 1))
    b.raw = b
    c = b[:, b.var.highly_variable].copy()
    sc.pp.scale(c, max_value=10)
    sc.tl.pca(c, n_comps=min(N_PCS, min(c.shape) - 1), svd_solver="arpack")
    sc.pp.neighbors(c, n_neighbors=N_NEIGHBORS, n_pcs=min(N_PCS, c.obsm["X_pca"].shape[1]))
    sc.tl.leiden(c, resolution=LEIDEN_RES, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    b.obs["leiden"] = ["L" + s for s in c.obs["leiden"].astype(str).values]
    lin_ratio, lin_cluster, lin_zero, detail = annotate(b)

    out = {k: z[k].copy() for k in z.files}
    for key, new in (("lineage", lin_ratio), ("lineage_cluster", lin_cluster),
                     ("lineage_cd3zero", lin_zero)):
        v = out[key].astype(object)
        v[in_lymph] = new
        out[key] = v.astype(str)
    lv = out["leiden"].astype(object)
    lv[in_lymph] = np.asarray(b.obs["leiden"])
    out["leiden"] = lv.astype(str)
    np.savez_compressed(args.out, **out)
    json.dump(clusters + detail, open(args.out.replace(".npz", "_clusters.json"), "w"),
              indent=1)
    import collections
    before = int((z["lineage"] == "NK").sum())
    after = int((out["lineage"] == "NK").sum())
    print(f"[refine] {os.path.basename(args.out)}: lymphoid={int(in_lymph.sum())} "
          f"NK {before} -> {after} | " +
          " ".join(f"{k}={v}" for k, v in collections.Counter(out["lineage"]).most_common(6)),
          flush=True)


if __name__ == "__main__":
    main()
