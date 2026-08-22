"""Task P-2: the mechanosensor panel across six sorted-NK bulk RNA-seq datasets.

Six datasets, four comparison groups (spec §1).  Every comparison is made
INSIDE one dataset only -- §3.1 forbids comparing absolute values across
datasets, and nothing here pools or meta-analyses them.  Cross-dataset use is
limited to checking whether the sign of an effect agrees.

  group 1  CD56bright vs CD56dim
      GSE236394  6 donors, blood, sorted bright/dim x CD8a+/-, paired
      GSE133383  4 donors, 5 tissues, sorted CD56brightCD16- / CD56dimCD16+
  group 2  resting vs cytokine-activated
      GSE140035  human: 6 donors x {UNSTIM, IFNA, IL12IL18, IL2IL15}, paired
      GSE140035  mouse: 3 replicates x 8 cytokine combinations, unpaired
      GSE242941  human: 3 control vs 3 cytokine-treated, unpaired
  group 3  blood vs tissue
      GSE133383  blood vs BM/spleen/lung/lung-LN, donor+subset held fixed
      GSE200319  blood vs liver perfusate, subset held fixed, unpaired
  group 4  tumour-infiltrating vs blood
      GSE205492  4 sarcoma patients, matched blood NK vs tumour NK, paired

Statistics (§3.3).  Paired designs get a paired t-test on the log2 scale.
Designs that are within-donor but not two-column -- GSE133383's tissue axis,
where donors contribute different tissue sets -- get an ordinary least squares
fit with donor and subset as fixed effects, so the contrast is estimated
within donor rather than across donors.  Genuinely unpaired designs get
Welch's t-test.  Every row carries log2FC, a 95% interval, p, and a BH q
computed across the panel within that comparison.

Detection floor (§3.3).  A gene whose expression sits near a dataset's own
detection limit cannot support a difference claim.  Two flags are recorded
per gene per dataset: CPM/TPM below 1 (where the dataset's unit allows it)
and a rank in the bottom decile of that dataset's own expression
distribution.  Either one sets NEAR_FLOOR.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import subprocess
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data", "p2")

PANEL = {
    "channel": ["PIEZO1", "PIEZO2", "TRPV2", "TRPV4", "TRPM7", "TRPC1",
                "TRPC6", "TMEM63A", "TMEM63B", "KCNK2", "KCNK4"],
    "mechano_downstream": ["YAP1", "WWTR1", "ANKRD1", "CCN2"],
    "state_activation": ["IFNG", "GZMB", "MKI67"],
    "state_residency": ["ITGA1", "CD69", "ITGAE"],
    "state_egress": ["KLF2", "S1PR1", "SELL"],
}
GENES = [g for v in PANEL.values() for g in v]
ROLE = {g: r for r, v in PANEL.items() for g in v}
# CTGF was renamed CCN2; some older annotations still carry the old symbol.
ALIASES = {"CCN2": ["CCN2", "CTGF"], "WWTR1": ["WWTR1", "TAZ"]}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


LOG: list[str] = []


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if not len(idx):
        return out
    order = idx[np.argsort(p[idx])]
    m = len(order)
    adj = p[order] * m / (np.arange(m) + 1)
    out[order] = np.minimum.accumulate(adj[::-1])[::-1].clip(0, 1)
    return out


def cpm_log2(counts: pd.DataFrame) -> pd.DataFrame:
    """counts (genes x samples) -> log2(CPM + 1)."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.divide(lib, axis=1) * 1e6 + 1)


def collapse_symbols(df: pd.DataFrame) -> pd.DataFrame:
    """Sum duplicate gene-symbol rows so a symbol appears exactly once."""
    return df.groupby(level=0).sum() if df.index.duplicated().any() else df


def mouse_symbols() -> dict:
    """human symbol -> [mouse symbols], from the verified ortholog table.

    Kept as a list because the verification found one gene in the panel that
    is NOT one-to-one: human GZMB maps to both mouse Gzmb and Gzmc.  Reading
    the mouse arm through a one-to-one assumption would silently pick one of
    them, which is the error §2 exists to prevent.
    """
    f = os.path.join(OUT, "p2_orthologs.csv")
    d = pd.read_csv(f)
    out = {}
    for _, r in d.iterrows():
        if isinstance(r.mouse_symbol, str) and r.mouse_symbol:
            out.setdefault(r.human, []).append(r.mouse_symbol)
    return out


def pick(expr: pd.DataFrame, gene: str, smap: dict | None = None):
    """Row for `gene`, trying its aliases / orthologues; None if absent."""
    names = (smap or {}).get(gene) or ALIASES.get(gene, [gene])
    for name in names:
        if name in expr.index:
            return expr.loc[name], name
    return None, None


def floor_flags(expr_log2: pd.DataFrame, gene_row, unit: str) -> dict:
    """§3.3 detection-floor flags, computed inside this dataset only.

    Two percentiles are reported because they answer different questions.
    `pct_all` ranks the gene against every row the depositors left in the
    matrix, which is inflated by however many never-expressed genes that is.
    `pct_expressed` ranks it only against rows that are above the matrix's
    own floor value, and it is the one the NEAR_FLOOR flag uses -- it is the
    stricter of the two and does not depend on the depositors' filtering.
    """
    mean_log2 = float(np.mean(gene_row))
    all_means = expr_log2.mean(axis=1)
    pct_all = float((all_means < mean_log2).mean() * 100)
    floor_val = float(all_means.min())
    expressed = all_means[all_means > floor_val + 1e-9]
    pct_expr = float((expressed < mean_log2).mean() * 100) if len(expressed) else np.nan
    if unit.startswith("log2(CPM"):
        linear = 2 ** mean_log2 - 1
    elif unit.startswith("moderated log2 CPM"):
        linear = 2 ** mean_log2
    else:
        linear = np.nan
    below_one = bool(linear < 1) if np.isfinite(linear) else None
    near = bool((np.isfinite(pct_expr) and pct_expr < 10) or (below_one is True))
    return {"mean_log2": mean_log2, "linear_cpm": linear,
            "pct_all_genes": pct_all, "pct_expressed_genes": pct_expr,
            "below_1_cpm": below_one, "near_floor": near}


def paired_test(a: np.ndarray, b: np.ndarray) -> dict:
    """a - b, paired.  a and b are aligned log2 values."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n < 2 or np.allclose(d, 0):
        return {"n": n, "log2FC": float(d.mean()) if n else np.nan,
                "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan,
                "test": "paired_t", "note": "no variation or n<2"}
    se = d.std(ddof=1) / np.sqrt(n)
    crit = stats.t.ppf(0.975, n - 1)
    t = stats.ttest_rel(a, b)
    return {"n": n, "log2FC": float(d.mean()), "ci_lo": float(d.mean() - crit * se),
            "ci_hi": float(d.mean() + crit * se), "p": float(t.pvalue),
            "test": "paired_t", "note": ""}


def welch_test(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {"n": len(a) + len(b), "log2FC": float(np.mean(a) - np.mean(b)),
                "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan,
                "test": "welch_t", "note": "n<2 in a group"}
    diff = float(a.mean() - b.mean())
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return {"n": len(a) + len(b), "log2FC": diff, "ci_lo": np.nan,
                "ci_hi": np.nan, "p": np.nan, "test": "welch_t",
                "note": "zero variance in both arms (gene flat or absent)"}
    t = stats.ttest_ind(a, b, equal_var=False)
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / ((va ** 2) / (len(a) - 1) + (vb ** 2) / (len(b) - 1))
    crit = stats.t.ppf(0.975, df) if np.isfinite(df) and df > 0 else np.nan
    return {"n": len(a) + len(b), "log2FC": diff,
            "ci_lo": diff - crit * se if np.isfinite(crit) else np.nan,
            "ci_hi": diff + crit * se if np.isfinite(crit) else np.nan,
            "p": float(t.pvalue), "test": "welch_t", "note": ""}


def ols_contrast(y: np.ndarray, effect: np.ndarray, *covariates) -> dict:
    """log2 expression ~ effect + fixed-effect covariates.

    Used where the design is within-donor but unbalanced, so a two-column
    paired test would silently drop most of the data.  Donor (and subset)
    enter as dummies, which is what keeps the contrast within donor.
    """
    y = np.asarray(y, float)
    cols = [np.ones(len(y)), np.asarray(effect, float)]
    for cov in covariates:
        cov = np.asarray(cov)
        for lvl in sorted(set(cov))[1:]:
            cols.append((cov == lvl).astype(float))
    X = np.column_stack(cols)
    ok = np.isfinite(y)
    X, y = X[ok], y[ok]
    n, k = X.shape
    if n <= k:
        return {"n": n, "log2FC": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                "p": np.nan, "test": "ols", "note": "not enough residual df"}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - np.linalg.matrix_rank(X)
    s2 = resid @ resid / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(s2 * XtX_inv[1, 1])
    t = beta[1] / se if se > 0 else np.nan
    p = 2 * stats.t.sf(abs(t), dof) if np.isfinite(t) else np.nan
    crit = stats.t.ppf(0.975, dof)
    return {"n": n, "log2FC": float(beta[1]), "ci_lo": float(beta[1] - crit * se),
            "ci_hi": float(beta[1] + crit * se), "p": float(p),
            "test": f"ols(+{len(covariates)} fixed effects)", "note": ""}


# ---------------------------------------------------------------------------
# dataset loaders
# ---------------------------------------------------------------------------
def load_gse236394():
    f = os.path.join(DATA, "GSE236394_all.gene_moderated_log2cpm.tsv.gz")
    df = pd.read_csv(f, sep="\t")
    expr = df.set_index("external_gene_name").filter(regex=r"^sample\.")
    expr.columns = [c.replace("sample.", "") for c in expr.columns]
    expr = expr[~expr.index.isna()]
    expr = expr.groupby(level=0).max()          # duplicate symbols: keep the
    meta = pd.DataFrame({"sample": expr.columns})  # better-measured copy
    meta["donor"] = meta["sample"].str.split("_").str[0]
    meta["subset"] = np.where(meta["sample"].str.contains("cd56bright"),
                              "CD56bright", "CD56dim")
    meta["cd8"] = np.where(meta["sample"].str.contains("cd8pos"), "CD8+", "CD8-")
    return expr, meta, {
        "dataset": "GSE236394",
        "sorting": "FACS: CD56/CD16 gate for bright vs dim, then CD8a+/-; "
                   "dump gate excluding CD3/CD19/CD14",
        "material": "peripheral blood NK, rested overnight in 1 ng/mL IL-15 "
                    "before sorting",
        "platform": "GPL24676 Illumina NovaSeq 6000",
        "unit": "moderated log2 CPM (as deposited)",
        "n_samples": expr.shape[1], "n_donors": meta.donor.nunique(),
        "paired": "yes, within donor",
        "depth": "not deposited (counts not provided) [待核实]"}


def load_gse133383():
    f = os.path.join(DATA, "GSE133383_filtered_NKcounts_table.csv.gz")
    counts = pd.read_csv(f, index_col=0)
    counts = collapse_symbols(counts)
    expr = cpm_log2(counts)
    meta = pd.DataFrame({"sample": expr.columns})
    parts = meta["sample"].str.split("_")
    meta["donor"] = parts.str[0]
    meta["tissue"] = parts.str[1]
    meta["subset"] = np.where(parts.str[2] == "I", "CD56bright", "CD56dim")
    meta["depth"] = counts.sum(axis=0).values
    return expr, meta, {
        "dataset": "GSE133383",
        "sorting": "FACS: CD56brightCD16- and CD56dimCD16+ NK",
        "material": "blood, bone marrow, spleen, lung, lung lymph node from "
                    "4 organ donors",
        "platform": "GPL16791 Illumina HiSeq 2500",
        "unit": "log2(CPM+1) from deposited filtered counts",
        "n_samples": expr.shape[1], "n_donors": meta.donor.nunique(),
        "paired": "yes, within donor (unbalanced tissue coverage)",
        "depth": f"median {np.median(meta.depth):,.0f} counts "
                 f"({meta.depth.min():,.0f}-{meta.depth.max():,.0f}); "
                 f"{counts.shape[0]} genes after the depositors' filter"}


def _gse140035(species):
    f = os.path.join(DATA, f"GSE140035_RawCounts_RNA_CytoStims_{species}.txt.gz")
    df = pd.read_csv(f, sep="\t")
    sym = df["gene.name"].astype(str).str.split(",").str[0].str.strip()
    counts = df.filter(regex=r"^RNA_").set_index(pd.Index(sym, name="gene"))
    counts = collapse_symbols(counts)
    expr = cpm_log2(counts)
    meta = pd.DataFrame({"sample": expr.columns})
    core = meta["sample"].str.replace(f"RNA_{species}_", "", regex=False)
    meta["condition"] = core.str.replace(r"^(HD|WT)_", "", regex=True).str.rsplit(
        "_", n=1).str[0]
    meta["replicate"] = core.str.rsplit("_", n=1).str[1]
    meta["depth"] = counts.sum(axis=0).values
    return counts, expr, meta


def load_gse140035_human():
    counts, expr, meta = _gse140035("human")
    return expr, meta, {
        "dataset": "GSE140035 (human)",
        "sorting": "ex vivo human NK cells; sorting markers not stated in the "
                   "GEO record [待核实]",
        "material": "healthy-donor NK, stimulated 'with various cytokine "
                    "combinations'",
        "platform": "GPL16791 Illumina HiSeq 2500",
        "unit": "log2(CPM+1) from deposited raw counts",
        "n_samples": expr.shape[1], "n_donors": meta.replicate.nunique(),
        "paired": "yes, donor 4-9 present in every condition",
        "depth": f"median {np.median(meta.depth):,.0f} counts "
                 f"({meta.depth.min():,.0f}-{meta.depth.max():,.0f})"}


def load_gse140035_mouse():
    counts, expr, meta = _gse140035("mouse")
    return expr, meta, {
        "dataset": "GSE140035 (mouse)",
        "sorting": "ex vivo mouse NK cells; sorting markers not stated in the "
                   "GEO record [待核实]",
        "material": "WT mouse NK, same cytokine panel as the human arm",
        "platform": "GPL21103 Illumina HiSeq 4000",
        "unit": "log2(CPM+1) from deposited raw counts",
        "n_samples": expr.shape[1], "n_donors": "3 replicates per condition",
        "paired": "no (replicates not matched across conditions)",
        "depth": f"median {np.median(meta.depth):,.0f} counts "
                 f"({meta.depth.min():,.0f}-{meta.depth.max():,.0f})"}


def load_gse242941():
    files = {"control": ["GSM7775364_NK_control_rep_1.txt.gz",
                         "GSM7775365_NK_control_rep_2.txt.gz",
                         "GSM7775366_NK_control_rep_3.txt.gz"],
             "cytokine": ["GSM7775370_NK_cytokine_control_1.txt.gz",
                          "GSM7775371_NK_cytokine_control_2.txt.gz",
                          "GSM7775372_NK_cytokine_control_3.txt.gz"]}
    cols, rows = {}, None
    for cond, fs in files.items():
        for i, fn in enumerate(fs, 1):
            d = pd.read_csv(os.path.join(DATA, fn), sep="\t")
            s = d.set_index("gene_name").iloc[:, -1]
            s = s.groupby(level=0).sum()
            cols[f"{cond}_{i}"] = s
            rows = s.index if rows is None else rows
    counts = pd.DataFrame(cols)
    # Library-size QC.  This is a self-supplied threshold, reported as such:
    # one library (control_1) has 32,276 total counts against 1.1-2.2M for
    # the rest and correlates 0.14-0.40 with every other sample where the
    # others correlate 0.63-0.93.  It is a failed library, and leaving it in
    # put PIEZO1's interval at [-5.6, +8.4] purely from its own noise.  The
    # cut is set at 1e5 counts, an order of magnitude below the next-lowest
    # library, so it removes that one sample and nothing else.
    keep = counts.sum(axis=0) >= 1e5
    dropped = sorted(counts.columns[~keep])
    counts = counts.loc[:, keep]
    expr = cpm_log2(counts)
    meta = pd.DataFrame({"sample": expr.columns})
    meta["condition"] = meta["sample"].str.split("_").str[0]
    meta["depth"] = counts.sum(axis=0).values
    return expr, meta, {
        "dataset": "GSE242941 (bulk arm)",
        "sorting": "NK cells sorted from PBMC; gate not stated in the GEO "
                   "record [待核实]",
        "material": "blood NK, untreated vs cytokine-treated; whether the three "
                    "replicates are three donors is not stated [待核实]",
        "platform": "GPL24676 Illumina NovaSeq 6000",
        "unit": "log2(CPM+1) from deposited raw counts",
        "n_samples": expr.shape[1], "n_donors": "3 replicates per arm [待核实]",
        "paired": "no",
        "qc": f"dropped for library size < 1e5 counts: {dropped or 'none'} "
              "(self-supplied threshold, see the loader docstring)",
        "depth": f"median {np.median(meta.depth):,.0f} counts "
                 f"({meta.depth.min():,.0f}-{meta.depth.max():,.0f}); "
                 "an order of magnitude shallower than the other datasets here"}


def load_gse200319():
    f = os.path.join(DATA, "GSE200319_ProcessedGeneList.csv.gz")
    # Layout: 3 banner rows (filename / sample type / cell type), then a
    # header row whose sample columns are unnamed, then the data.  Columns
    # 0-5 are annotation, column 6 is a blank spacer, 7.. are the samples.
    with gzip.open(f, "rt") as fh:
        raw = fh.read().splitlines()
    fname = raw[0].split(",")
    stype = raw[1].split(",")
    ctype = raw[2].split(",")
    samples = [x.replace(".hs.bam", "") for x in fname[7:]]
    names = raw[3].split(",")[:6] + ["_spacer"] + samples
    body = pd.read_csv(io.StringIO("\n".join(raw[4:])), header=None, names=names)
    expr = body.set_index("Gene Symbol")[samples].apply(pd.to_numeric,
                                                        errors="coerce")
    expr = expr.groupby(level=0).max()
    meta = pd.DataFrame({"sample": samples,
                         "tissue": stype[7:], "subset": ctype[7:]})
    return expr, meta, {
        "dataset": "GSE200319",
        "sorting": "FACS: CD56hiCD16-, CD56hiCD16+, CD56loCD16+ NK subsets",
        "material": "liver perfusate (n=5 donors) vs healthy blood (n=5 donors)",
        "platform": "GPL18573 Illumina NextSeq 500; modified SMART-Seq2, "
                    "low input",
        "unit": "log2 normalised expression as deposited (not CPM)",
        "n_samples": expr.shape[1], "n_donors": "5 + 5, different individuals",
        "paired": "no (blood and liver come from different donors)",
        "depth": "not deposited [待核实]"}


def _ensembl_symbol_map(genes):
    """symbol -> human ENSG, cached, for the dataset that ships ENSG ids."""
    cache = os.path.join(OUT, "p2_symbol_to_ensg.csv")
    if os.path.exists(cache):
        d = pd.read_csv(cache)
        return dict(zip(d.symbol, d.ensg))
    out = {}
    for g in genes:
        r = subprocess.run(
            ["curl", "-sS", "--max-time", "60",
             f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}"
             "?content-type=application/json"], capture_output=True, text=True).stdout
        try:
            out[g] = json.loads(r).get("id", "")
        except Exception:
            out[g] = ""
        time.sleep(1)
    pd.DataFrame({"symbol": list(out), "ensg": list(out.values())}).to_csv(
        cache, index=False)
    return out


def load_gse205492():
    f = os.path.join(DATA, "GSE205492_rnaseq_workshop_normalized_counts.txt.gz")
    df = pd.read_csv(f, sep="\t", index_col=0)
    df.index = [i.split(".")[0] for i in df.index]
    smap = _ensembl_symbol_map(GENES + ["CTGF"])
    rev = {v: k for k, v in smap.items() if v}
    nk = [c for c in df.columns if "_NK_" in c]
    expr = np.log2(df[nk] + 1)
    expr.index = [rev.get(i, i) for i in expr.index]
    expr = expr.groupby(level=0).max()
    meta = pd.DataFrame({"sample": nk})
    meta["tissue"] = np.where(meta["sample"].str.contains("TUMOR"), "tumour", "blood")
    # blood libraries are named CCS24002/3/11 + SA1556; tumour ones X002/3/11 +
    # SA1556.  The trailing digits are the patient id in both schemes.
    meta["patient"] = (meta["sample"].str.replace("CCS24", "", regex=False)
                       .str.replace("^X", "", regex=True).str.split("_").str[0]
                       .str.lstrip("0"))
    return expr, meta, {
        "dataset": "GSE205492",
        "sorting": "FACS: CD3-CD56+ NK from matched blood and tumour",
        "material": "soft-tissue sarcoma patients (liposarcoma, "
                    "myxofibrosarcoma, UPS)",
        "platform": "GPL20301 Illumina HiSeq 4000",
        "unit": "log2(normalised counts + 1) as deposited (not CPM)",
        "n_samples": len(nk), "n_donors": meta.patient.nunique(),
        "paired": "yes, matched blood and tumour from the same patient",
        "depth": "not deposited (normalised counts only) [待核实]"}


# ---------------------------------------------------------------------------
# comparison runner
# ---------------------------------------------------------------------------
def run_comparison(expr, meta, info, name, group, fn, smap=None):
    """fn(gene_row, meta) -> dict from one of the *_test helpers."""
    rows = []
    for gene in GENES:
        row, used = pick(expr, gene, smap)
        if row is None:
            rows.append({"dataset": info["dataset"], "comparison": name,
                         "group": group, "gene": gene, "role": ROLE[gene],
                         "present": False, "note": "gene absent from this matrix"})
            continue
        res = fn(row, meta)
        fl = floor_flags(expr, row.values, info["unit"])
        rows.append({"dataset": info["dataset"], "comparison": name,
                     "group": group, "gene": gene, "role": ROLE[gene],
                     "symbol_used": used, "present": True, "unit": info["unit"],
                     "ortholog_note": ("human GZMB is one-to-many in mouse "
                                       "(Gzmb + Gzmc); this row reads only "
                                       f"{used}")
                     if (smap and gene == "GZMB") else "",
                     **res, **fl})
    df = pd.DataFrame(rows)
    if "p" in df:
        df["q_BH"] = bh(df["p"].to_numpy(dtype=float))
    if "ci_lo" in df:
        # The largest effect still compatible with the interval.  A null
        # result only means something if this bound is small; otherwise the
        # comparison did not have the power to see a moderate change and
        # "no difference" is a statement about the study, not the biology.
        df["max_effect_still_compatible"] = np.nanmax(
            np.abs(df[["ci_lo", "ci_hi"]].to_numpy(dtype=float)), axis=1)
        df["informative_null"] = (df["q_BH"] >= 0.05) & (
            df["max_effect_still_compatible"] < 0.585)   # < 1.5-fold
    return df


def main():
    os.makedirs(OUT, exist_ok=True)
    all_rows, metas = [], []
    log("=== Task P-2: mechanosensor panel in sorted-NK bulk RNA-seq ===")

    # ---------------- group 1 -------------------------------------------
    expr, meta, info = load_gse236394()
    metas.append(info)
    log(f"\n[{info['dataset']}] {expr.shape[0]} genes x {expr.shape[1]} samples, "
        f"{info['n_donors']} donors")

    def bright_vs_dim_236394(row, m):
        # one value per donor per subset: mean of the CD8+ and CD8- arms
        v = pd.DataFrame({"donor": m.donor.values, "subset": m.subset.values,
                          "y": row[m["sample"]].values})
        w = v.groupby(["donor", "subset"]).y.mean().unstack("subset")
        w = w.dropna()
        return paired_test(w["CD56bright"].values, w["CD56dim"].values)

    all_rows.append(run_comparison(expr, meta, info,
                                   "CD56bright vs CD56dim (blood)", "1",
                                   bright_vs_dim_236394))

    expr3, meta3, info3 = load_gse133383()
    metas.append(info3)
    log(f"[{info3['dataset']}] {expr3.shape[0]} genes x {expr3.shape[1]} samples, "
        f"{info3['n_donors']} donors, tissues {sorted(set(meta3.tissue))}")

    def bright_vs_dim_133383(row, m):
        y = row[m["sample"]].values
        eff = (m.subset == "CD56bright").astype(float).values
        return ols_contrast(y, eff, m.donor.values, m.tissue.values)

    all_rows.append(run_comparison(expr3, meta3, info3,
                                   "CD56bright vs CD56dim (5 tissues, donor+tissue fixed)",
                                   "1", bright_vs_dim_133383))

    # ---------------- group 2 -------------------------------------------
    exprH, metaH, infoH = load_gse140035_human()
    metas.append(infoH)
    log(f"[{infoH['dataset']}] {exprH.shape[0]} genes x {exprH.shape[1]} samples, "
        f"conditions {sorted(set(metaH.condition))}")

    for cond in ("IL2IL15", "IL12IL18", "IFNA"):
        def f(row, m, cond=cond):
            a = m[m.condition == cond].sort_values("replicate")
            b = m[m.condition == "UNSTIM"].sort_values("replicate")
            common = sorted(set(a.replicate) & set(b.replicate))
            av = [row[a[a.replicate == r]["sample"].iloc[0]] for r in common]
            bv = [row[b[b.replicate == r]["sample"].iloc[0]] for r in common]
            return paired_test(np.array(av), np.array(bv))
        all_rows.append(run_comparison(exprH, metaH, infoH,
                                       f"{cond} vs UNSTIM", "2", f))

    MOUSE = mouse_symbols()
    exprM, metaM, infoM = load_gse140035_mouse()
    metas.append(infoM)
    log(f"[{infoM['dataset']}] {exprM.shape[0]} genes x {exprM.shape[1]} samples, "
        f"conditions {sorted(set(metaM.condition))}")
    for cond in ("IL2IL15", "IL12IL18", "IFNA"):
        def fm(row, m, cond=cond):
            a = row[m[m.condition == cond]["sample"]].values
            b = row[m[m.condition == "UNSTIM"]["sample"]].values
            return welch_test(a, b)
        all_rows.append(run_comparison(exprM, metaM, infoM,
                                       f"{cond} vs UNSTIM (mouse)", "2", fm,
                                       smap=MOUSE))

    expr2, meta2, info2 = load_gse242941()
    metas.append(info2)
    log(f"[{info2['dataset']}] {expr2.shape[0]} genes x {expr2.shape[1]} samples")

    def cyto_242941(row, m):
        a = row[m[m.condition == "cytokine"]["sample"]].values
        b = row[m[m.condition == "control"]["sample"]].values
        return welch_test(a, b)

    all_rows.append(run_comparison(expr2, meta2, info2,
                                   "cytokine vs control", "2", cyto_242941))

    # ---------------- group 3 -------------------------------------------
    def tissue_vs_blood(row, m, tissue=None):
        sub = m if tissue is None else m[m.tissue.isin(["BL", tissue])]
        y = row[sub["sample"]].values
        eff = (sub.tissue != "BL").astype(float).values
        return ols_contrast(y, eff, sub.donor.values, sub.subset.values)

    all_rows.append(run_comparison(
        expr3, meta3, info3,
        "all tissues vs blood (donor+subset fixed)", "3", tissue_vs_blood))
    for t in ("BM", "SP", "LU", "LN"):
        all_rows.append(run_comparison(
            expr3, meta3, info3, f"{t} vs blood (donor+subset fixed)", "3",
            lambda row, m, t=t: tissue_vs_blood(row, m, t)))

    exprL, metaL, infoL = load_gse200319()
    metas.append(infoL)
    log(f"[{infoL['dataset']}] {exprL.shape[0]} genes x {exprL.shape[1]} samples, "
        f"tissues {sorted(set(metaL.tissue))}, subsets {sorted(set(metaL.subset))}")

    def liver_vs_blood(row, m):
        shared = {s for s in set(m.subset)
                  if len(set(m[m.subset == s].tissue)) == 2}
        sub = m[m.subset.isin(shared)]
        y = row[sub["sample"]].values
        eff = (sub.tissue == "Liver Perfusate").astype(float).values
        return ols_contrast(y, eff, sub.subset.values)

    all_rows.append(run_comparison(exprL, metaL, infoL,
                                   "liver perfusate vs blood (subset fixed)",
                                   "3", liver_vs_blood))

    # ---------------- group 4 -------------------------------------------
    exprT, metaT, infoT = load_gse205492()
    metas.append(infoT)
    log(f"[{infoT['dataset']}] {exprT.shape[0]} genes x {exprT.shape[1]} samples, "
        f"{infoT['n_donors']} patients, pairs "
        f"{sorted(set(metaT.patient))}")

    def tumour_vs_blood(row, m):
        w = pd.DataFrame({"patient": m.patient.values, "tissue": m.tissue.values,
                          "y": row[m["sample"]].values})
        w = w.pivot(index="patient", columns="tissue", values="y").dropna()
        return paired_test(w["tumour"].values, w["blood"].values)

    all_rows.append(run_comparison(exprT, metaT, infoT,
                                   "tumour NK vs blood NK (matched)", "4",
                                   tumour_vs_blood))

    res = pd.concat(all_rows, ignore_index=True)
    res.to_csv(os.path.join(OUT, "p2_panel_results.csv"), index=False)
    pd.DataFrame(metas).to_csv(os.path.join(OUT, "p2_dataset_metadata.csv"),
                               index=False)

    # ---------------- reporting -----------------------------------------
    log("\n=== PIEZO1, every comparison ===")
    p = res[res.gene == "PIEZO1"]
    log(p[["dataset", "comparison", "n", "log2FC", "ci_lo", "ci_hi", "p", "q_BH",
           "linear_cpm", "pct_expressed_genes", "near_floor",
           "max_effect_still_compatible", "informative_null"]]
        .round(3).to_string(index=False))
    log(f"\n  PIEZO1 comparisons with q<0.05: "
        f"{int((p.q_BH < 0.05).sum())} of {len(p)}")
    log(f"  PIEZO1 measurable (not near floor) in: "
        f"{int((~p.near_floor.astype(bool)).sum())} of {len(p)} comparisons; "
        f"CPM range {p.linear_cpm.min():.1f}-{p.linear_cpm.max():.1f} "
        "(datasets reporting CPM)")
    log(f"  informative nulls (CI excludes a 1.5-fold change): "
        f"{int(p.informative_null.sum())} of {len(p)}")

    log("\n=== channels reaching q<0.05 in any comparison ===")
    sig = res[(res.role == "channel") & (res.q_BH < 0.05)]
    if len(sig):
        log(sig[["dataset", "comparison", "gene", "log2FC", "q_BH",
                 "near_floor"]].round(4).to_string(index=False))
    else:
        log("  none")

    log("\n=== state controls: do the two arms actually differ? ===")
    st = res[res.role.str.startswith("state") & (res.q_BH < 0.05)]
    log(st[["dataset", "comparison", "gene", "log2FC", "q_BH"]]
        .round(3).to_string(index=False))

    with open(os.path.join(OUT, "p2_bulk_analysis.log"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    log("\nDONE")


if __name__ == "__main__":
    main()
