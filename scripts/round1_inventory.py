"""Round-1 inventory of two published NK atlases. Reads metadata only.

Deliverable 1: a field list for every NK-level object, plus the
required-field checklist.
Deliverable 2: sample / patient / NK counts by cancer type x tissue source,
with healthy donors in separate rows.

This script reads cell metadata (obs) and never touches expression values.
It computes no signature scores and no differential expression.

Local files (verified by md5 against the Zenodo API; see scripts/atlas_fetch.sh):
  data/atlasA/comb_CD56_CD16_NK.h5ad.gz        Zenodo 10163908 (plain HDF5 despite .gz)
  data/atlasA/comb_CD56_CD16_NK_blood.h5ad.gz  Zenodo 10163908 (plain HDF5 despite .gz)
  data/atlasB/adata_nk_tumor_query.h5ad        Zenodo 14178285 (concept 8434223)
  data/atlasB/adata_ref_nk.h5ad                Zenodo 14178285
Remote Atlas B objects: only obs is read, through HTTP range requests
(scripts/remote_h5.py): pb_12_donors, adata_all_nk_after_mapping,
adata_all_nk_milo, adata_ref_after_training.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from remote_h5 import RangeFile  # noqa: E402

OUT = ROOT / "results" / "round1"
OUT.mkdir(parents=True, exist_ok=True)
ZB = "https://zenodo.org/records/14178285/files/{}?download=1"

OBJECTS = {
    # label: (atlas, source, role)
    "A_main": ("A", ROOT / "data/atlasA/comb_CD56_CD16_NK.h5ad.gz",
               "Zenodo 10163908 comb_CD56_CD16_NK.h5ad.gz: all tissues, used in most of the paper"),
    "A_blood": ("A", ROOT / "data/atlasA/comb_CD56_CD16_NK_blood.h5ad.gz",
                "Zenodo 10163908 comb_CD56_CD16_NK_blood.h5ad.gz: circulating-NK analysis"),
    "B_milo": ("B", ZB.format("adata_all_nk_milo.h5ad"),
               "Zenodo 14178285 adata_all_nk_milo.h5ad: reference + TiNK, published nhood groups; matches paper 27,732 TrNK / 38,982 TiNK"),
    "B_after_mapping": ("B", ZB.format("adata_all_nk_after_mapping.h5ad"),
                        "Zenodo 14178285 adata_all_nk_after_mapping.h5ad: reference + TiNK, pre-Milo filter"),
    "B_tumor_query": ("B", ROOT / "data/atlasB/adata_nk_tumor_query.h5ad",
                      "Zenodo 14178285 adata_nk_tumor_query.h5ad: TiNK query"),
    "B_ref_nk": ("B", ROOT / "data/atlasB/adata_ref_nk.h5ad",
                 "Zenodo 14178285 adata_ref_nk.h5ad: reference NK used for mapping"),
    "B_ref_after_training": ("B", ZB.format("adata_ref_after_training.h5ad"),
                             "Zenodo 14178285 adata_ref_after_training.h5ad"),
    "B_pb_12_donors": ("B", ZB.format("pb_12_donors.h5ad"),
                       "Zenodo 14178285 pb_12_donors.h5ad: healthy PB-NK, 12 donors; matches paper 44,640"),
}

# The other h5ad objects in the Atlas B record are all-cell or PB-derived
# objects rather than the NK objects. Their obs is inventoried for
# Deliverable 1 only and is not used in the census.
B_EXTRA = ["lung_tumor", "lung_normal", "breast_tumor", "breast_normal", "glioblastoma", "melanoma",
           "sarcoma", "pancreas_tumor", "pancreas_normal", "prostate_tumor", "prostate_normal",
           "skin_normal", "pb_12_donors_bulk", "pb_12_donors_processed", "pb_12_donors_pseudotime"]
for _k in B_EXTRA:
    OBJECTS[f"Bx_{_k}"] = ("B", ZB.format(f"{_k}.h5ad"),
                           f"Zenodo 14178285 {_k}.h5ad: not an NK-only object; D1 field list only")

# Strings counted as missing, as well as real NaN/None.
MISSING_TOKENS = {"", "nan", "none", "na", "n/a", "null", "unknown", "notavailable", "not available"}


def open_h5(src):
    if isinstance(src, Path):
        return h5py.File(src, "r")
    return h5py.File(RangeFile(src), "r")


CACHE = ROOT / "data" / "obs_cache"  # under data/, which is git-ignored


def load(label):
    atlas, src, role = OBJECTS[label]
    hit = CACHE / f"{label}.pkl"
    if hit.exists():
        return pd.read_pickle(hit)
    out = _load(label, atlas, src, role)
    CACHE.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(out, hit)
    return out


def _load(label, atlas, src, role):
    f = open_h5(src)
    obs = read_elem(f["obs"])
    x = f["X"]
    shape = tuple(int(v) for v in (x.attrs["shape"] if "shape" in x.attrs else x.shape))
    layers = sorted(f["layers"].keys()) if "layers" in f else []
    info = dict(object=label, atlas=atlas, role=role, n_cells=shape[0], n_genes=shape[1],
                has_raw=("raw" in f), layers=";".join(layers),
                uns_keys=";".join(sorted(f["uns"].keys())) if "uns" in f else "")
    f.close()
    return obs, info


def is_missing(s: pd.Series) -> pd.Series:
    m = s.isna()
    if s.dtype == object or isinstance(s.dtype, pd.CategoricalDtype):
        m |= s.astype(str).str.strip().str.lower().isin(MISSING_TOKENS)
    return m


def field_table(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in obs.columns:
        s = obs[c]
        miss = is_missing(s)
        categorical = isinstance(s.dtype, pd.CategoricalDtype) or s.dtype == object or s.dtype == bool
        if categorical:
            vc = s[~miss].astype(str).value_counts()
            n = len(vc)
            shown = vc.head(20)
            vals = "; ".join(f"{k} ({v})" for k, v in shown.items())
            if n > 20:
                vals += f"; ... [first 20 of {n} values, by cell count]"
            dtype = "categorical" if isinstance(s.dtype, pd.CategoricalDtype) else str(s.dtype)
        else:
            n, dtype = s.nunique(), str(s.dtype)
            vals = "(numeric)"  # no summary statistics in round 1
        rows.append(dict(field=c, dtype=dtype, non_missing_frac=round(1 - miss.mean(), 4),
                         n_non_missing=int((~miss).sum()), n_values=int(n), values=vals))
    return pd.DataFrame(rows)


# ------------------------------------------------------------ keyword scan
# The scan only provides evidence for the checklist. Every checklist call
# below names an exact field or says 缺失 (missing).
KEYWORDS = {
    "metastasis": r"metast|\bm[_ ]?stage|\bstage\b|tnm|distant|\bmet\b|site",
    "survival": r"surviv|\bos\b|pfs|dfs|event|death|dead|alive|vital|follow|relapse|recur",
    "timepoint": r"time|pre|post|treat|therapy|baseline|on[_-]?treat|visit|day|\bd\d|cycle",
    "tissue": r"tissue|source|site|origin|compartment|organ|sample[_ ]?type",
    "age_sex": r"\bage\b|sex|gender|female|male",
    "batch_platform": r"batch|platform|tech|chem|protocol|10x|seq|library|kit|run|lane",
}


def keyword_scan(label, obs):
    rows = []
    for c in obs.columns:
        for k, pat in KEYWORDS.items():
            if re.search(pat, c, re.I):
                rows.append(dict(object=label, hit_in="field_name", field=c, keyword_group=k, match=c))
        s = obs[c]
        if isinstance(s.dtype, pd.CategoricalDtype) or s.dtype == object:
            vals = pd.Series(s.dropna().astype(str).unique())
            for k, pat in KEYWORDS.items():
                if k in ("batch_platform", "tissue"):
                    continue  # value-level scan kept to the clinically relevant groups
                hits = vals[vals.str.contains(pat, case=False, regex=True)]
                if len(hits):
                    rows.append(dict(object=label, hit_in="values", field=c, keyword_group=k,
                                     match=f"{len(hits)} values, e.g. " + "; ".join(hits.head(5))))
    return rows


# ------------------------------------------------------------ checklist
# Each call was made by reading the full field tables. A cell names the
# exact obs field only where that field directly encodes the item.
# Otherwise it says 缺失 (missing). No approximating field is substituted.
MISSING = "缺失"
DEDUP_LOG = []
HEME = ["Chronic Lymphocytic Leukemia(CLL)", "Multiple Myeloma(MM)",
        "Acute Lymphocytic Leukemia(ALL)", "Acute Myeloid Leukemia(AML)"]


ALL_B = ["B_milo", "B_after_mapping", "B_tumor_query", "B_ref_nk", "B_ref_after_training"]
B_REF = ["B_milo", "B_after_mapping", "B_ref_nk", "B_ref_after_training"]
B_TUM = ["B_milo", "B_after_mapping", "B_tumor_query"]
# item -> {object: (field, regex over the field's values selecting the level or None, note)}
ITEMS = [
    ("转移状态 — M 分期", {}),
    ("转移状态 — 是否远处转移", {}),
    ("转移状态 — 转移部位", {}),
    ("随访终点 — 总生存期 (OS) 及事件指示", {}),
    ("随访终点 — 无进展生存期 (PFS) 及事件指示", {}),
    ("随访终点 — 无病生存期 (DFS) 及事件指示", {}),
    ("外周血样本采样时点 (治疗前/中/后)", {}),
    ("组织来源标签 — 字段", {
        "A_main": ("meta_tissue_in_paper", None, ""),
        "A_blood": ("meta_tissue", None, "常数字段：全部为 Blood"),
        **{o: ("source", None, "器官与组织合一的标签；另有 reference 字段区分参照/肿瘤") for o in ALL_B}}),
    ("组织来源 — 肿瘤", {
        "A_main": ("meta_tissue_in_paper", r"^Tumor$", "注意：血液肿瘤（如 CLL 外周血 GSE111015）也标为 Tumor"),
        **{o: ("source", r"_tumor$|^glioblastoma$|^melanoma$|^sarcoma$", "") for o in B_TUM}}),
    ("组织来源 — 癌旁", {}),
    ("组织来源 — 外周血", {
        "A_main": ("meta_tissue_in_paper", r"^Blood$", ""),
        "A_blood": ("meta_tissue", r"^Blood$", ""),
        **{o: ("source", r"^PBMC$", "仅健康供者") for o in B_REF}}),
    ("组织来源 — 淋巴结", {}),
    ("组织来源 — 正常组织", {
        "A_main": ("meta_tissue_in_paper", r"^Normal$", "单一 Normal 取值，与癌旁不可区分"),
        **{o: ("source", r"_normal$", "参照组织；对象内未注明是否癌旁") for o in B_REF}}),
    ("组织来源 — 其他（未细分）", {
        "A_main": ("meta_tissue_in_paper", r"^Other tissue$",
                   "混合类：含淋巴结、转移灶、积液、腹水及部分正常组织，对象内不可拆分")}),
    ("健康供者年龄", {}),
    ("健康供者性别", {}),
    ("批次 — 作者整合批次键 (batch)", {
        "A_main": ("batch", None, "= 病人ID + 数据集 (BBKNN batch_key)，与病人一一对应，非技术批次"),
        "A_blood": ("batch", None, "= 病人ID + 数据集 (BBKNN batch_key)，与病人一一对应，非技术批次"),
        **{o: ("batch", None, "样本/文库级 ID，用作 scVI 批次键") for o in ALL_B + ["B_pb_12_donors"]}}),
    ("批次 — 技术批次 (测序批次/文库制备/化学版本)", {}),
    ("测序平台", {"A_main": ("meta_platform", None, "11 个拼写取值约对应 8 个平台；无 10x 化学版本")}),
]


# Extra B objects: level rows come from a regex over the values of their
# own tissue/source field. That field is `tissue` in the per-tissue objects
# and `source` in pb_12_donors_bulk.
EXTRA_LEVELS = {
    "组织来源标签 — 字段": None,
    # Only values that literally name a tumour count (tumor, cancer, melanoma,
    # sarcoma, ...). Abbreviations such as tLung or TEC_* are not decoded,
    # and bare organ names ("lung", "prostate") do not count.
    "组织来源 — 肿瘤": r"(?i)tumou?r|cancer|carcinoma|melanoma|sarcoma|glioblastoma",
    "组织来源 — 癌旁": r"(?i)adjacent|\badj",
    "组织来源 — 外周血": r"(?i)pbmc|blood",
    "组织来源 — 淋巴结": r"(?i)lymph|\bLN\b",
    "组织来源 — 正常组织": r"(?i)normal",
}


def _extra_spec(item, obs):
    if item.startswith("批次 — 作者整合批次键") and "batch" in obs.columns:
        return ("batch", None, "样本/文库级 ID")
    fld = "tissue" if "tissue" in obs.columns else "source" if "source" in obs.columns else None
    if fld and item in EXTRA_LEVELS:
        return (fld, EXTRA_LEVELS[item], "按取值正则匹配，取值原样列出")
    return None


def checklist(obs_by_obj):
    rows = []
    for item, present in ITEMS:
        for obj, obs in obs_by_obj.items():
            spec = present.get(obj) if not obj.startswith("Bx_") else _extra_spec(item, obs)
            if not spec:
                rows.append(dict(item=item, object=obj, status=MISSING, field="", values="", note=""))
                continue
            fld, level, note = spec
            assert fld in obs.columns, (obj, fld)
            vc = obs[fld].astype(str).value_counts()
            if level is not None:
                vc = vc[vc.index.str.contains(level, regex=True)]
                if vc.sum() == 0:
                    rows.append(dict(item=item, object=obj, status=MISSING, field="", values="", note=""))
                    continue
            shown = "; ".join(f"{k} ({v})" for k, v in vc.head(12).items())
            rows.append(dict(item=item, object=obj, status="存在", field=fld,
                             values=f"{len(vc)} 个取值: {shown}" + ("; ..." if len(vc) > 12 else ""),
                             note=note))
    return pd.DataFrame(rows)


# ------------------------------------------------------------ census
def per_group(df, sample_key, patient_key, cols):
    rows = []
    for key, g in df.groupby(cols, observed=True):
        per_sample = g.groupby(sample_key, observed=True).size()
        q1, med, q3 = np.percentile(per_sample.values, [25, 50, 75])
        rows.append(dict(zip(cols, key if isinstance(key, tuple) else (key,)),
                         n_samples=len(per_sample), n_patients=g.groupby(patient_key, observed=True).ngroups,
                         n_patient_ids_raw=g[patient_key[-1]].nunique(),
                         n_nk=len(g), nk_per_sample_median=med, nk_per_sample_q1=q1, nk_per_sample_q3=q3,
                         nk_per_sample_iqr=q3 - q1,
                         n_datasets=g[patient_key[0]].nunique(),
                         **({} if "cancer_type" in cols else {"n_cancer_types": g["cancer_type"].nunique()})))
    return pd.DataFrame(rows)


def census_A(obs, tissue_col, label):
    o = obs.copy()
    # Normalise dataset labels by dropping the ",part_N" suffix, so that one
    # submission split into parts counts once: "No_552,part_1" -> "No_552".
    # This changes the counting key only; no biological label is modified.
    o["dataset_key"] = o["datasets"].astype(str).str.replace(r",part_\d+$", "", regex=True)
    # Some cells were deposited twice, with the same patient and cellID
    # under two dataset labels. In A_main these are 98 cells of CID4471,
    # CID44971 and CID4513, in GSE176078 and 9_Cryopreservation. The copy in
    # the dataset holding fewer of that patient's cells is dropped.
    # Barcodes repeated inside one dataset are left alone: they can be real
    # collisions between samples.
    key = o["meta_patientID"].astype(str) + "|" + o["cellID"].astype(str)
    cross = key.map(o.groupby(key, observed=True)["dataset_key"].nunique()) > 1
    n_pt = o.groupby(["meta_patientID", "dataset_key"], observed=True).size()
    size = pd.Series(list(zip(o["meta_patientID"], o["dataset_key"])), index=o.index).map(n_pt)
    keep = ~cross | (size == size.groupby(key).transform("max"))
    DEDUP_LOG.append(dict(object=label, cells_in=len(o), cross_dataset_duplicate_cells_dropped=int((~keep).sum()),
                          patients=";".join(sorted(o.loc[~keep, "meta_patientID"].astype(str).unique())),
                          datasets_dropped_from=";".join(sorted(o.loc[~keep, "dataset_key"].unique()))))
    o = o[keep].copy()
    o["cancer_type"] = o["meta_histology"].astype(str)
    o["tissue"] = o[tissue_col].astype(str)
    sk, pk = ["dataset_key", "sampleID"], ["dataset_key", "meta_patientID"]
    o["group"] = np.where(o["cancer_type"] == "Healthy donor", "healthy_donor", "patient")
    tab = per_group(o, sk, pk, ["group", "cancer_type", "tissue"])
    tab.insert(0, "object", label)
    tot = per_group(o, sk, pk, ["group", "tissue"])
    tot.insert(0, "object", label)
    # Patients (dataset-scoped key) with NK cells in both Blood and Tumor in
    # this object. This is a count of paired designs, not a test.
    # The Tumor label also covers haematological malignancies, whose
    # samples are blood or bone marrow. A separate row keeps only solid
    # types. The list is explicit: the four meta_histology labels of
    # haematological cancers.
    solid = o[(o.group == "patient") & (o.tissue == "Tumor") & ~o.cancer_type.isin(HEME)].copy()
    solid["tissue"] = "Tumor (solid types only)"
    st = per_group(solid, sk, pk, ["group", "tissue"])
    st.insert(0, "object", label)
    tot = pd.concat([tot, st], ignore_index=True)
    pt = o[o.group == "patient"].groupby(pk, observed=True)["tissue"].agg(set)
    tot["patients_with_blood_and_tumor"] = int(pt.map(lambda t: {"Blood", "Tumor"} <= t).sum())
    return tab, tot


def census_B(milo, after_mapping, pb12):
    o = milo.join(after_mapping[["patient"]], how="left")
    assert o["patient"].notna().all(), "patient join from after_mapping incomplete"
    src = o["source"].astype(str)
    o["cancer_type"] = np.where(o["reference"].astype(str) == "tumor",
                                src.str.replace("_tumor", "", regex=False),
                                np.where(src == "PBMC", "Healthy donor", "—"))
    o["tissue"] = np.select(
        [o["reference"].astype(str) == "tumor", src == "PBMC"],
        ["Tumor (reference=tumor)", "PBMC (reference)"],
        default=src.radd("Reference tissue: "))
    o["group"] = np.where(o["tissue"] == "PBMC (reference)", "healthy_donor",
                          np.where(o["cancer_type"] == "—", "reference_tissue", "patient"))
    ob = o.rename(columns={"dataset": "dataset_key"})
    tab = per_group(ob, ["dataset_key", "sample"], ["dataset_key", "patient"], ["group", "cancer_type", "tissue"])
    tab.insert(0, "object", "B_milo")
    tot = per_group(ob, ["dataset_key", "sample"], ["dataset_key", "patient"], ["group", "tissue"])
    tot.insert(0, "object", "B_milo")
    # Sensitivity: in dataset "chan" the sample field is one ID per cell
    # (n_samples == n_cells), which drags the per-sample median to 1.
    # Detected within each (dataset, reference/tumor) part, because a dataset
    # can use per-cell IDs in one part only. Only that part is dropped.
    part = o.groupby(["dataset", "reference"], observed=True)
    flag = (part["sample"].nunique() == part.size()) & (part.size() > 5)
    per_cell = sorted(flag[flag].index)
    drop = pd.MultiIndex.from_frame(o[["dataset", "reference"]].astype(str)).isin(
        [tuple(map(str, k)) for k in per_cell])
    o2 = o[~drop].rename(columns={"dataset": "dataset_key"})
    sens = per_group(o2, ["dataset_key", "sample"], ["dataset_key", "patient"], ["group", "cancer_type", "tissue"])
    sens.insert(0, "object", f"B_milo excluding per-cell-sample datasets {per_cell}")
    sens_tot = per_group(o2, ["dataset_key", "sample"], ["dataset_key", "patient"], ["group", "tissue"])
    sens_tot.insert(0, "object", f"B_milo excluding {per_cell}")
    # Healthy PB donors, full 44,640-cell object. Batches *_sorted are sorted
    # subset fractions from donors malm1/malm2, so the donor key is the batch
    # with the suffix dropped. This is a counting key only.
    p = pb12.copy()
    p["donor"] = p["batch"].astype(str).str.replace(r"_(sorted|bulk)$", "", regex=True)
    # The object has no dataset or donor field. The cohort is read from the
    # sample-name prefix (amir/crinier/malm/yang) and the donor from the batch
    # name. Both are counting keys only.
    p["dataset_key"] = p["sample"].astype(str).str.replace(r"\d+$", "", regex=True)
    p["cancer_type"], p["tissue"] = "Healthy donor", "PBMC (pb_12_donors, bulk + sorted)"
    p["group"] = "healthy_donor"
    pbt = per_group(p, ["dataset_key", "sample"], ["dataset_key", "donor"], ["group", "cancer_type", "tissue"])
    pbt.insert(0, "object", "B_pb_12_donors")
    pbtot = per_group(p, ["dataset_key", "sample"], ["dataset_key", "donor"], ["group", "tissue"])
    pbtot.insert(0, "object", "B_pb_12_donors")
    return (pd.concat([tab, pbt], ignore_index=True), sens, per_cell,
            pd.concat([tot, sens_tot, pbtot], ignore_index=True))


def main():
    obs_by_obj, infos, fields, scans = {}, [], [], []
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(4) as ex:
        loaded = dict(zip(OBJECTS, ex.map(load, OBJECTS)))
    for label in OBJECTS:
        obs, info = loaded[label]
        obs_by_obj[label] = obs
        infos.append(info)
        ft = field_table(obs)
        ft.insert(0, "object", label)
        fields.append(ft)
        scans += keyword_scan(label, obs)
        print(f"{label}: {obs.shape}", flush=True)
    pd.DataFrame(infos).to_csv(OUT / "objects.csv", index=False)
    pd.concat(fields).to_csv(OUT / "D1_field_list.csv", index=False)
    pd.DataFrame(scans).to_csv(OUT / "D1_keyword_scan.csv", index=False)
    checklist(obs_by_obj).to_csv(OUT / "D1_required_fields.csv", index=False)

    a_main, a_main_tot = census_A(obs_by_obj["A_main"], "meta_tissue_in_paper", "A_main")
    a_blood, a_blood_tot = census_A(obs_by_obj["A_blood"], "meta_tissue", "A_blood")
    b_tab, b_sens, per_cell, b_tot = census_B(obs_by_obj["B_milo"], obs_by_obj["B_after_mapping"],
                                       obs_by_obj["B_pb_12_donors"])
    for name, t in [("D2_census_A_main", a_main), ("D2_census_A_blood", a_blood),
                    ("D2_census_B", b_tab), ("D2_census_B_sensitivity", b_sens)]:
        t.sort_values(["group", "tissue", "n_nk"], ascending=[True, True, False]).to_csv(
            OUT / f"{name}.csv", index=False, float_format="%.2f")

    # Compartment totals, counted directly from cells (never by summing rows).
    pd.concat([a_main_tot, a_blood_tot, b_tot], ignore_index=True).to_csv(
        OUT / "D2_compartment_totals.csv", index=False, float_format="%.2f")
    pd.DataFrame(DEDUP_LOG).to_csv(OUT / "D2_dedup_log.csv", index=False)
    print("per-cell sample datasets (B):", per_cell)
    write_markdown()


def _md(df):
    df = df.fillna("").astype(str).apply(lambda c: c.str.replace("|", "\\|", regex=False))
    head = "| " + " | ".join(df.columns) + " |\n|" + "---|" * len(df.columns) + "\n"
    return head + "\n".join("| " + " | ".join(r) + " |" for r in df.values) + "\n"


def write_markdown():
    """Render the CSVs as Markdown tables so they read directly on GitHub."""
    objs = pd.read_csv(OUT / "objects.csv")
    fl = pd.read_csv(OUT / "D1_field_list.csv")
    for atlas in ("A", "B"):
        parts = [f"# Deliverable 1 — Atlas {atlas} metadata field list\n\n"
                 "Generated by `scripts/round1_inventory.py` from obs only. Categorical fields list "
                 "values by descending cell count, with the count in parentheses. When a field has "
                 "more than 20 values, the first 20 are shown and the total is given.\n"]
        for _, o in objs[objs.atlas == atlas].iterrows():
            parts.append(f"\n## `{o.object}`: {o.n_cells:,} cells × {o.n_genes:,} genes\n\n"
                         f"{o.role}. raw: {o.has_raw}; layers: {o.layers if isinstance(o.layers, str) else '—'}\n\n")
            t = fl[fl.object == o.object].drop(columns=["object", "n_non_missing"])
            t.columns = ["字段名", "数据类型", "非缺失比例", "取值数", "取值 (前 20, 细胞数)"]
            parts.append(_md(t))
        (OUT / f"D1_field_list_{atlas}.md").write_text("".join(parts))
    ck = pd.read_csv(OUT / "D1_required_fields.csv")
    (OUT / "D1_required_fields.md").write_text(
        "# Deliverable 1 — required-field checklist\n\n" + _md(ck.rename(columns={
            "item": "项目", "object": "对象", "status": "状态", "field": "字段", "values": "取值",
            "note": "说明"})))
    cols = ["group", "cancer_type", "tissue", "n_samples", "n_patients", "n_patient_ids_raw", "n_nk",
            "nk_per_sample_median", "nk_per_sample_q1", "nk_per_sample_q3", "n_datasets"]
    zh = ["组别", "瘤种", "组织来源", "样本数", "病人数", "病人ID数(原始)", "NK 细胞数",
          "每样本 NK 中位数", "Q1", "Q3", "数据集数"]
    parts = ["# Deliverable 2 — sample census by cancer type × tissue source\n"]
    for name in ["D2_census_A_main", "D2_census_A_blood", "D2_census_B", "D2_census_B_sensitivity"]:
        t = pd.read_csv(OUT / f"{name}.csv")
        parts.append(f"\n## {name} ({', '.join(dict.fromkeys(t.object))})\n\n")
        t = t[cols].copy()
        for c in ["nk_per_sample_median", "nk_per_sample_q1", "nk_per_sample_q3"]:
            t[c] = t[c].map(lambda v: f"{v:.2f}".rstrip("0").rstrip("."))
        t.columns = zh
        parts.append(_md(t))
    (OUT / "D2_census.md").write_text("".join(parts))


if __name__ == "__main__":
    main()
