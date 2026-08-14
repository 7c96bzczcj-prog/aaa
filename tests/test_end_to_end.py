"""End-to-end smoke test on a synthetic dataset.

Builds a tiny dense-TSV dataset with a manifest, runs ingest -> validate ->
qc -> detection_rates -> donor_tests -> soup_calibration, and asserts the
tables come out with the right shape and the pre-registered flags set.

Run:  .venv/bin/python -m pytest tests/test_end_to_end.py -q
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from dnkchem.dataset import load_panel  # noqa: E402

PY = os.path.join(ROOT, ".venv", "bin", "python")
if not os.path.exists(PY):
    PY = sys.executable


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    """A 5-donor decidua + blood dataset with a planted dNK2 > dNK1 effect."""
    d = tmp_path_factory.mktemp("synth")
    panel = load_panel(os.path.join(ROOT, "panel", "chemokine_panel_v1.tsv"))
    genes = panel["gene_symbol"].tolist()
    ens = panel["ensembl_id"].tolist()
    # padding genes so the null distribution has somewhere to draw from
    extra = [f"FILLER{i}" for i in range(300)]
    genes = genes + extra
    ens = ens + [f"ENSGX{i:011d}" for i in range(300)]

    rng = np.random.default_rng(0)
    donors = ["S1", "S2", "S3", "S4", "S5"]
    rows = []
    meta = []
    plan = [("Decidua", "dNK1", 60), ("Decidua", "dNK2", 60), ("Decidua", "dNK3", 45),
            ("Decidua", "dM1", 60), ("Decidua", "dS1", 60), ("Decidua", "DC1", 40),
            ("Decidua", "Tcells", 50), ("Blood", "NK CD16+", 60)]
    gi = {g: i for i, g in enumerate(genes)}
    for di, donor in enumerate(donors):
        for comp, ct, n in plan:
            for k in range(n):
                v = rng.poisson(0.25, size=len(genes)).astype(int)
                v[gi["KLRD1"]] += rng.poisson(6)     # NK-ish everywhere, fine
                if ct.startswith("dNK") or ct.startswith("NK"):
                    v[gi["NKG7"]] += rng.poisson(8)
                    v[gi["XCL1"]] += rng.poisson(4 if ct == "dNK2" else 1)
                if ct == "dM1":
                    v[gi["LYZ"]] += rng.poisson(30)
                    v[gi["C1QA"]] += rng.poisson(20)
                if ct == "dS1":
                    v[gi["DCN"]] += rng.poisson(40)
                    v[gi["COL1A1"]] += rng.poisson(40)
                if ct == "DC1":
                    v[gi["XCR1"]] += rng.poisson(6)
                if ct == "Tcells":
                    for g in ("CD3E", "CD3D", "TRBC2"):
                        v[gi[g]] += rng.poisson(5)
                rows.append(v)
                meta.append((f"RUN{di}_{comp[:2]}_{ct.replace(' ', '')}_{k}", donor, comp, ct))
    M = np.vstack(rows)                      # cells x genes
    ids = [m[0] for m in meta]

    mtx = d / "counts.txt"
    with open(mtx, "w") as fh:
        fh.write("Gene\t" + "\t".join(ids) + "\n")
        for j, g in enumerate(genes):
            fh.write(f"{g}_{ens[j]}\t" + "\t".join(map(str, M[:, j])) + "\n")

    pd.DataFrame({"cell": ids, "donor": [m[1] for m in meta],
                  "loc": [m[2] for m in meta], "ann": [m[3] for m in meta]}
                 ).to_csv(d / "meta.txt", sep="\t", index=False)

    mdir = d / "manifests"
    mdir.mkdir()
    (mdir / "SYNTH.yaml").write_text(f"""
dataset_id: SYNTH
citation: "synthetic"
accession: "none"
local_path: "{d}/SYNTH.h5ad"
file_format: "h5ad"
species: "human"
source:
  path: "{mtx}"
  format: "dense_tsv_genes_by_cells"
  obs_table: "{d}/meta.txt"
gene_id:
  var_index_type: "symbol_ensembl_concat"
  id_separator: "_"
  symbol_field: 0
  ensembl_field: -1
  symbol_column: null
  genome_build: "GRCh38"
obs_columns:
  donor: "donor"
  compartment: "loc"
  celltype: "ann"
  library: null
fixed_values:
  library: "unspecified"
celltype_map:
  dNK1: "dNK1"
  dNK2: "dNK2"
  dNK3: "dNK3"
  "NK CD16+": "pbNK_CD56dim"
  dM1: "Myeloid"
  dS1: "Stromal"
  DC1: "cDC1"
  Tcells: "T"
counts_layer: "X"
raw_counts_verified: false
notes: |
  synthetic fixture
""")
    return d, str(mdir / "SYNTH.yaml")


def run(script, manifest, cwd, extra=()):
    env = dict(os.environ)
    r = subprocess.run([PY, os.path.join(ROOT, "src", script), manifest,
                        "--panel", os.path.join(ROOT, "panel", "chemokine_panel_v1.tsv"),
                        *extra],
                       cwd=cwd, capture_output=True, text=True, env=env)
    return r


def test_pipeline_end_to_end(synth):
    d, manifest = synth
    cwd = str(d)

    r = subprocess.run([PY, os.path.join(ROOT, "src", "ingest.py"), manifest],
                       cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "round-trip verified" in r.stdout

    r = run("validate_manifest.py", manifest, cwd)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "PASS" in r.stdout

    r = run("qc_report.py", manifest, cwd, extra=["--min-genes", "5"])
    assert r.returncode == 0, r.stdout + r.stderr

    r = run("detection_rates.py", manifest, cwd, extra=["--min-genes", "5", "--min-cells", "20"])
    assert r.returncode == 0, r.stdout + r.stderr

    T1 = pd.read_csv(os.path.join(cwd, "out", "SYNTH", "detection_rates.tsv"), sep="\t")
    assert {"dataset_id", "donor_id", "compartment", "subset", "gene", "n_cells",
            "n_detected", "detection_rate", "mean_cpm", "depth_floor",
            "cells_retained_frac", "qc_flag"} == set(T1.columns)
    assert T1["depth_floor"].nunique() >= 2, "sensitivity floors must be present"

    r = run("donor_tests.py", manifest, cwd, extra=["--min-cells", "20"])
    assert r.returncode == 0, r.stdout + r.stderr
    T2 = pd.read_csv(os.path.join(cwd, "out", "SYNTH", "donor_level_tests.tsv"), sep="\t")
    frozen_v10 = ["dataset_id", "comparison_id", "group_a", "group_b", "gene",
                  "n_donors", "mean_diff_pp", "ci_low", "ci_high", "p_raw", "q_bh",
                  "min_achievable_p", "on_test_floor", "uninterpretable_ceiling",
                  "uninterpretable_floor", "cross_compartment_soup_caveat"]
    # the v1.0 schema is frozen: same columns, same order, at the front
    assert list(T2.columns)[:len(frozen_v10)] == frozen_v10
    # v1.1 appends, never reorders or removes
    assert set(T2.columns) - set(frozen_v10) == {
        "ceiling_arm_a", "ceiling_arm_b", "magnitude_is_lower_bound"}

    # v1.1 A1: the ceiling rule is per arm
    both = T2.uninterpretable_ceiling
    assert (T2.loc[both, "ceiling_arm_a"] & T2.loc[both, "ceiling_arm_b"]).all(), \
        "uninterpretable_ceiling must require BOTH arms saturated"
    lb = T2[T2.magnitude_is_lower_bound & T2.n_donors.gt(0)]
    assert not lb["uninterpretable_ceiling"].any(), \
        "a one-sided ceiling is a lower bound, not uninterpretable"
    for _, r in lb.iterrows():
        if r.ceiling_arm_b and not r.ceiling_arm_a:
            assert r.mean_diff_pp > 0, "effect must point toward the saturated arm"
        elif r.ceiling_arm_a and not r.ceiling_arm_b:
            assert r.mean_diff_pp < 0

    # R2 travels with every test
    tested = T2[T2.n_donors > 0]
    assert tested["min_achievable_p"].notna().all()
    assert (tested["p_raw"] >= tested["min_achievable_p"] - 1e-12).all()

    # analysis B rows carry the caveat unconditionally, analysis A rows do not
    assert T2[T2.comparison_id.str.startswith("B_")]["cross_compartment_soup_caveat"].all()
    assert not T2[T2.comparison_id.str.startswith("A_")]["cross_compartment_soup_caveat"].any()

    # the planted dNK2 > dNK1 XCL1 effect is recovered with the right sign
    xcl1 = T2[(T2.comparison_id == "A_decidua_dNK1_vs_dNK2") & (T2.gene == "XCL1")]
    assert len(xcl1) == 1
    assert xcl1.iloc[0]["mean_diff_pp"] > 0, "dNK2 should exceed dNK1 on XCL1"

    # dNKp is never tested against the trio
    assert not T2.comparison_id.str.contains("dNKp").any()

    # v1.1 A2: the empirical null carries the inference and gets its own BH
    r = run("null_distribution.py", manifest, cwd,
            extra=["--min-genes", "5", "--min-cells", "20", "--n-null", "270",
                   "--min-null-per-target", "5"])
    assert r.returncode == 0, r.stdout + r.stderr
    T4 = pd.read_csv(os.path.join(cwd, "out", "SYNTH", "null_distribution.tsv"), sep="\t")

    # v1.3 A8: no FDR is reported, and no q column survives
    assert "q_bh_empirical" not in T4.columns
    assert "q_bh_empirical_per_comparison" not in T4.columns
    assert "empirical_p" not in T4.columns, "a rank p would be read as a verdict"
    assert (T4["fdr_estimable"] == False).all()  # noqa: E712
    assert T4["fdr_not_estimable_reason"].notna().all()
    assert "NO FDR REPORTED" in r.stdout

    # v1.3 A9: nulls are matched on baseline detection, and both the
    # standardised and the raw effect must be present so neither is read alone
    for c in ("baseline_detection_ref_arm", "empirical_z", "observed_diff_pp",
              "reference_arm", "rank_of_observed"):
        assert c in T4.columns, c
    # the reference arm of a contrast is its group_a
    assert (T4.loc[T4.comparison_id == "A_decidua_dNK1_vs_dNK2",
                   "reference_arm"] == "dNK1").all()

    # the planted XCL1 effect should stand out against its matched null
    x = T4[(T4.comparison_id == "A_decidua_dNK1_vs_dNK2") & (T4.gene == "XCL1")]
    assert len(x) == 1 and x.iloc[0]["empirical_z"] > 2

    r = run("soup_calibration.py", manifest, cwd, extra=["--min-genes", "5"])
    assert r.returncode == 0, r.stdout + r.stderr
    T3 = pd.read_csv(os.path.join(cwd, "out", "SYNTH", "soup_calibration.tsv"), sep="\t")
    assert list(T3.columns) == ["dataset_id", "compartment", "gene", "category",
                                "soup_fraction", "nk_myeloid_cpm_ratio",
                                "scale_a_verdict", "scale_b_verdict", "combined_verdict"]
    assert set(T3.combined_verdict) <= {"positive", "weak",
                                        "indistinguishable_from_ambient", "unmeasurable"}
    # v1.1 A3: ruler B denominator is per gene
    sup = pd.read_csv(os.path.join(cwd, "out", "SYNTH",
                                   "soup_calibration_supplement.tsv"), sep="\t")
    for c in ("dominant_source_lineage", "nk_dominant_source_cpm_ratio",
              "scale_b_verdict_v10_myeloid_only"):
        assert c in sup.columns, c
    dec = sup[sup.compartment == "Decidua"].set_index("gene")
    # DCN is planted in the stromal population only, so its source must be Stromal
    if "DCN" in dec.index:
        assert dec.loc["DCN", "dominant_source_lineage"] == "Stromal"
    if "LYZ" in dec.index:
        assert dec.loc["LYZ", "dominant_source_lineage"] == "Myeloid"

    # XCR1 is expressed only by cDC1 here, so NK pickup must not read as real
    x = T3[(T3.gene == "XCR1") & (T3.compartment == "Decidua")]
    if len(x):
        assert x.iloc[0]["combined_verdict"] != "positive"


def test_panel_mismatch_aborts(synth, tmp_path):
    """R11: a panel that matches almost nothing must abort, not print a table."""
    d, manifest = synth
    bad = tmp_path / "bad_panel.tsv"
    pd.DataFrame({"gene_symbol": [f"NOTAGENE{i}" for i in range(20)],
                  "ensembl_id": [""] * 20, "category": ["ligand"] * 20,
                  "role": ["primary_target"] * 20, "notes": [""] * 20}
                 ).to_csv(bad, sep="\t", index=False)
    r = subprocess.run([PY, os.path.join(ROOT, "src", "validate_manifest.py"),
                        manifest, "--panel", str(bad)],
                       cwd=str(d), capture_output=True, text=True)
    assert r.returncode != 0
    assert "panel match rate" in (r.stdout + r.stderr)
