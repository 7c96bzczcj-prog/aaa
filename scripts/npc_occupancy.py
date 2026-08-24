#!/usr/bin/env python3
"""Answer Q-A .. Q-F from the manifest and write OCCUPANCY.md + COUNTS.txt."""
import csv, json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = os.path.join(ROOT, "raw"), os.path.join(ROOT, "out")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ftlib  # noqa: E402

RUN_DATE = "2026-08-24"
CAVEAT = ("本结论基于 %s 的检索，NOT_FOUND 表示未检索到证据，不等于不存在。" % RUN_DATE)


def load_manifest():
    rows = []
    with open(os.path.join(OUT, "npc_sc_manifest.tsv")) as f:
        lines = [l for l in f if not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        rows.append(r)
    return rows


def ident(r):
    return r["pmid"] if r["pmid"] != "NOT_FOUND" else ("doi:" + r["doi"])


def cascade_scan(finals, recs):
    """Q-E: which included records frame results by metastatic-cascade stages."""
    PAT = (r"metastatic\s+cascade|invasion[\s\-]metastasis\s+cascade|"
           r"step[\s\-]?wise\s+(?:process\s+of\s+)?metasta|"
           r"seed(?:ing)?\s+and\s+coloniz|"
           r"(?:lymphatic|haematogenous|hematogenous)\s+route|"
           r"evolutionary\s+route\s+of\s+.{0,20}metasta|"
           r"primary\s*(?:->|→|to|vs\.?)\s*(?:regional\s+)?lymph\s+node\s*(?:->|→|to)\s*distant")
    hits = {}
    for rid in finals["included"]:
        r = recs[rid]
        ti = re.sub(r"<[^>]+>", "", r.get("title") or "")
        ab = re.sub(r"<[^>]+>", " ", r.get("abstractText") or "")
        xml = ftlib.load_xml(r["pmcid"]) if r.get("pmcid") else None
        ft = ftlib.parse(xml)["plain"] if xml else ""
        blob = ti + " " + ab + " " + ft
        m = re.search(PAT, blob, re.I)
        if m:
            sents = ftlib.sentences_with(blob, PAT, maxn=2)
            hits[rid] = re.sub(r"\s+", " ", sents[0])[:200] if sents else m.group(0)
    return hits


def main():
    rows = load_manifest()
    finals = json.load(open(os.path.join(RAW, "final_sets.json")))
    recs = {r["id"]: r for r in json.load(open(os.path.join(RAW, "survivors.json")))}
    ds = json.load(open(os.path.join(RAW, "datasets_rows.json")))
    zen = json.load(open(os.path.join(RAW, "anchor15_zenodo_blood.json")))
    zall = json.load(open(os.path.join(RAW, "anchor15_zenodo_all.json")))
    anchors = json.load(open(os.path.join(RAW, "anchor_check.json")))
    dd = json.load(open(os.path.join(RAW, "dedup_counts.json")))

    DEEP = {"subclustered", "NK_is_theme"}

    qa = [r for r in rows if r["tissue_blood"] == "Y" and r["NK_depth"] in DEEP]
    qb = [r for r in rows if r["tissue_distant_met"].startswith("Y") or r["M_stage"] == "M1_included"]
    qb_blood = [r for r in qb if r["tissue_blood"] == "Y"]
    qc = [r for r in rows if r["longitudinal"].startswith("Y") and r["tissue_blood"] == "Y"]
    qc_all = [r for r in rows if r["longitudinal"].startswith("Y")]
    qd = [r for r in rows if r["EBV_DNA_link"] == "continuous"]
    qe = cascade_scan(finals, recs)
    qf = [d for d in ds if d["is_NPC"] == "Y" and d["is_single_cell"] == "Y" and d["includes_blood"] == "Y"]
    qf_ctrl = [d for d in ds if d["repository"] in ("GSA-Human", "NODE", "OMIX", "CNGB", "EGA")]

    L = []
    A = L.append
    A("# NPC 单细胞占位审计 — 占位问题答案")
    A("")
    A("**run date**: %s  " % RUN_DATE)
    A("**检索来源**: Europe PMC REST (Q1–Q5, Q7–Q8) + NCBI E-utilities db=gds (Q6)  ")
    A("**included 记录数**: %d  |  **secondary_analysis**: %d  |  **excluded**: %d"
      % (len(rows), len(finals["secondary"]), dd["deduped"] - len(rows) - len(finals["secondary"])))
    A("")
    A("> Q7/Q8 是 spec §5 强制的 query 修复：anchor #15 (PMID 37607536) 在 Europe PMC 中非开放获取")
    A("> (isOpenAccess=N, 无 PMCID)，全文未被索引，任何以 TITLE_ABS:\"nasopharyngeal\" 为轴的 query")
    A("> 都不可能召回它。Q7/Q8 改以 pan-cancer 单细胞图谱语言为轴，召回后再按 NPC 样本证据筛选。")
    A("")
    A("---")
    A("")

    # ---- Q-A
    A("## Q-A 有几篇做了 NPC 外周血的单细胞，且 NK_depth ≥ subclustered？")
    A("")
    A("**答：%d 篇。**" % len(qa))
    A("")
    A("| pmid / doi | year | modality | n_patients | NK_depth | NK_claim |")
    A("|---|---|---|---|---|---|")
    for r in qa:
        A("| %s | %s | %s | %s | %s | %s |" % (ident(r), r["year"], r["modality"],
                                               r["n_patients"], r["NK_depth"], r["NK_claim"][:90]))
    A("")
    A("PMID 列表：%s" % ", ".join(ident(r) for r in qa))
    A("")
    A("**关键限定**（决定这个答案是「没人做过」还是「图谱里已经有了」）：")
    A("")
    A("- 血 + NK 深度分析的记录里，只有 **PMID 41992060**(Br J Cancer 2026) 是以 NPC 外周血 NK 亚群")
    A("  为主题的原始研究；其摘要明确 `%s`。" % "CD16+CD57+ NK cells in blood correlated with better patient outcomes")
    A("- **PMID 37607536**(anchor #15, Cell 2023 泛癌 NK 图谱) 确实含 NPC 循环 NK：从 Zenodo record")
    A("  8275845 的 `comb_CD56_CD16_NK_blood.h5ad` 直接读出 —— `meta_histology` 取值")
    A("  `Nasopharyngeal Carcinoma(NPC)` 共 **%d 个细胞 / %d 位病人 / %d 个样本**，"
      % (zen["NPC_cells"], zen["NPC_patients"], zen["NPC_samples"]))
    A("  全部来自单一数据集 **%s**，即 anchor #1 (PMID 33531485) 的 10 对瘤-血配对样本。" % zen["NPC_source_datasets"][0])
    A("  这些细胞在图谱里被细分到 %d 个 NK 亚群（最大三个：%s）。"
      % (len(zen["NPC_cells_by_NK_subset"]),
         "、".join("%s %d" % (k, v) for k, v in list(zen["NPC_cells_by_NK_subset"].items())[:3])))
    A("- 因此：NPC 循环 NK 的**单细胞数据**在公开图谱里已经存在并已被亚群化，但其**原始样本**")
    A("  仍只有 GSE162025 这一份 10 人队列；没有第二份独立的 NPC 外周血 NK 单细胞队列。")
    A("")
    A(CAVEAT)
    A("")

    # ---- Q-B
    A("## Q-B 有几篇纳入了远处转移灶（M1）样本？其中几篇同时有血？")
    A("")
    A("**答：%d 篇纳入远处转移灶样本；其中 %d 篇同时有血。**" % (len(qb), len(qb_blood)))
    A("")
    A("| pmid / doi | year | tissue_distant_met | M_stage | tissue_blood | n_patients |")
    A("|---|---|---|---|---|---|")
    for r in qb:
        A("| %s | %s | %s | %s | %s | %s |" % (ident(r), r["year"], r["tissue_distant_met"][:60],
                                               r["M_stage"], r["tissue_blood"], r["n_patients"]))
    A("")
    A("PMID 列表（M1）：%s" % ", ".join(ident(r) for r in qb))
    A("PMID 列表（M1 且有血）：%s" % (", ".join(ident(r) for r in qb_blood) or "（无）"))
    A("")
    A("注：PMID 36739462 的 163 份配对血+组织样本是**基因组**层面（WES/WGS）；其单细胞部分只有")
    A("2 位病人的 11 个样本，血未进入单细胞。这是唯一一篇 M1 + 血同时出现的记录，但两者不在同一层。")
    A("")
    A(CAVEAT)
    A("")

    # ---- Q-C
    A("## Q-C 有几篇对同一病人做了纵向重复取血的单细胞？时间点分别是什么？")
    A("")
    A("**答：%d 篇。**（另有 %d 篇是纵向的但取的是组织而非血）" % (len(qc), len(qc_all) - len(qc)))
    A("")
    A("| pmid / doi | year | modality | 时间点（原文） |")
    A("|---|---|---|---|")
    for r in qc:
        A("| %s | %s | %s | %s |" % (ident(r), r["year"], r["modality"], r["longitudinal"][:110]))
    A("")
    A("PMID 列表：%s" % ", ".join(ident(r) for r in qc))
    A("")
    A("纵向但非血的记录：%s" % ", ".join(ident(r) for r in qc_all if r not in qc))
    A("")
    A(CAVEAT)
    A("")

    # ---- Q-D
    A("## Q-D 有几篇把免疫单细胞状态对连续血浆 EBV DNA 载量做了回归（而非阳性/阴性二分）？")
    A("")
    from npc_manifest_data import M as MD
    # key the lookup exactly the way ident() keys manifest rows, so records with no
    # PMID (preprints) resolve too
    by_pmid = {}
    for k, v in MD.items():
        pm = recs[k].get("pmid")
        by_pmid[pm if pm else ("doi:" + (recs[k].get("doi") or k))] = (k, v)
    A("**答：%d 篇。**" % len(qd))
    A("")
    if qd:
        A("| pmid / doi | year | EBV_DNA_link | EBV 变量类型 | 原文依据 |")
        A("|---|---|---|---|---|")
        for r in qd:
            rid, md = by_pmid[ident(r)]
            A("| %s | %s | %s | %s | %s |" % (ident(r), r["year"], r["EBV_DNA_link"],
                                              md.get("ebv_kind", "-"), md.get("ebv_snip", "-")[:130]))
        A("")
        A("PMID 列表：%s" % ", ".join(ident(r) for r in qd))
    A("")
    cat = [r for r in rows if r["EBV_DNA_link"] == "categorical"]
    A("其余 included 记录的 EBV_DNA_link 分布：categorical %d 篇、none %d 篇。"
      % (len(cat), sum(1 for r in rows if r["EBV_DNA_link"] == "none")))
    A("")
    A("**`categorical` 的两种子类型必须分开看** —— spec 的字段名是 `EBV_DNA_link`，但被判 categorical")
    A("的 6 篇里只有一部分真的用了**血浆 EBV DNA**，其余用的是**肿瘤组织 EBV 状态**(EBER ISH /")
    A("EBV-high vs EBV-low)，两者不是同一个变量：")
    A("")
    A("| pmid | EBV 变量类型 | 原文依据 |")
    A("|---|---|---|")
    for r in cat:
        rid, md = by_pmid[ident(r)]
        A("| %s | %s | %s |" % (ident(r), md.get("ebv_kind", "-"), md.get("ebv_snip", "-")[:130]))
    A("")
    plasma = [r for r in rows if by_pmid[ident(r)][1].get("ebv_kind") == "plasma_EBV_DNA"]
    A("即：**真正把单细胞结果连到血浆 EBV DNA 的只有 %d 篇**（%s），其中只有 %d 篇是连续量。"
      % (len(plasma), ", ".join(ident(r) for r in plasma), len(qd)))
    A("")
    A("注：被判为 `continuous` 的唯一一篇是 PMID 32686767（Cell Res 2020），其全文写明")
    A("`We further explored their correlations with plasma EBV DNA concentrations.`——把免疫细胞")
    A("signature 对**连续**血浆 EBV DNA 浓度作相关，而非 Sero+/Sero− 二分。作为对照，anchor #12")
    A("(PMID 36860875) 全篇以 EBV DNA **血清阳性/阴性**分组，属 `categorical`；但该记录在本轮")
    A("被判为 secondary_analysis（其 scRNA-seq 数据取自 GSE150430），故不在 included 表内。")
    A("")
    A("**没有任何一篇**把免疫单细胞状态对连续血浆 EBV DNA 拷贝数做**回归建模**；32686767 做的是")
    A("相关性探索。若把标准严格限定为「回归」，答案是 0 篇。")
    A("")
    A(CAVEAT)
    A("")

    # ---- Q-E
    A("## Q-E 有几篇明确用 metastatic cascade 分阶段的框架来组织免疫学结论？")
    A("")
    A("**答：严格口径 %d 篇。**（下面同时给出较宽的结构性口径）" % len(qe))
    A("")
    A("**严格口径** —— 检索文本里出现明确的 cascade / 分阶段路径措辞：")
    A("")
    A("| pmid / doi | year | 原文依据 |")
    A("|---|---|---|")
    for rid, sn in qe.items():
        r = recs[rid]
        A("| %s | %s | %s |" % (r.get("pmid") or ("doi:" + (r.get("doi") or rid)),
                                r.get("pubYear"), sn[:150]))
    A("")
    A("PMID 列表：%s" % ", ".join((recs[k].get("pmid") or ("doi:" + (recs[k].get("doi") or k)))
                                  for k in qe))
    A("")
    staged = [r for r in rows
              if sum(1 for k in ("tissue_primary", "tissue_LN_met", "tissue_distant_met")
                     if r[k].startswith("Y")) >= 2]
    A("**结构性口径** —— 不要求用 cascade 这个词，只要求同一研究在单细胞层面同时纳入")
    A("解剖学上 ≥2 个阶段（原发灶 / 区域淋巴结转移 / 远处转移），共 %d 篇：" % len(staged))
    A("")
    A("| pmid / doi | year | primary | LN_met | distant_met | blood |")
    A("|---|---|---|---|---|---|")
    for r in staged:
        A("| %s | %s | %s | %s | %s | %s |" % (ident(r), r["year"], r["tissue_primary"],
                                               r["tissue_LN_met"], r["tissue_distant_met"][:28],
                                               r["tissue_blood"]))
    A("")
    A("PMID 列表：%s" % ", ".join(ident(r) for r in staged))
    A("")
    A("两个口径的差距本身就是答案：**采样上跨越 ≥2 个转移阶段的有 %d 篇，但只有 %d 篇真的把结论"
      % (len(staged), len(qe)))
    A("按 cascade 阶段组织**；其余都是两点对比（原发 vs 转移），没有分阶段框架。")
    A("")
    A(CAVEAT)
    A("")

    # ---- Q-F
    A("## Q-F 公开可下载、且含外周血的 NPC 单细胞数据集有哪些？")
    A("")
    A("**答：%d 个。**" % len(qf))
    A("")
    A("| accession | repository | n_samples | 组织类型（据 GEO sample title） | 关联 PMID |")
    A("|---|---|---|---|---|")
    ORIGIN = {"GSE162025": "33531485 (anchor #1, Nat Commun 2021)",
              "GSE202604": "无关联论文（Europe PMC 检索 'GSE202604' 命中 0 条，属仅存 GEO 的数据集）"}
    for d in qf:
        st = d["sample_titles"]
        tt = "primary tumour + PBMC" if "PBMC" in st else d["title"][:60]
        cited = ",".join(sorted(set(x for x in d["cited_by_included"].split(",") if x and x != "-")))
        A("| %s | %s | %s | %s | originating: %s; re-used by: %s |"
          % (d["accession"], d["repository"], d["n_samples"], tt,
             ORIGIN.get(d["accession"], "-"), cited or "-"))
    A("")
    for d in qf:
        A("- **%s** — %s  " % (d["accession"], d["title"][:150]))
        A("  sample titles: `%s`" % d["sample_titles"][:260])
    A("")
    A("**注意 GSE202604 的血不是外周血**：其 series title 写明取的是 *right atrial blood*")
    A("（右心房血，术中/导管取样），不是肘静脉外周血。若把 Q-F 严格限定为**外周血**，")
    A("公开可下载的 NPC 单细胞数据集只有 **GSE162025 一个**。该数据集同时也是 anchor #15")
    A("泛癌 NK 图谱中全部 NPC 循环 NK 细胞的唯一来源 —— 也就是说，")
    A("**目前全世界公开可下载的 NPC 外周血单细胞数据，追溯到底只有这一份 10 人队列。**")
    A("")
    A("**同样含血、但不是自由下载的**（受控访问，需申请）：")
    A("")
    A("| accession | repository | 关联 PMID | 说明 |")
    A("|---|---|---|---|")
    A("| HRA006885 | GSA-Human | 39231979 | anchor #11 的 77 份活检+血 scRNA-seq 原始数据；GSA-Human 为受控访问 |")
    A("| HRA008590 | GSA-Human | 40315843 | TIL-ACT 试验的 scRNA/scTCR；含 PBMC 分离流程 |")
    A("")
    A("**重要更正**：GSE206245（anchor #11 的 GEO 条目）的 series summary 复述了论文摘要中的")
    A("「77 biopsy and blood samples」，但其**实际存放的 sample 全部是空间转录组**（NPC_ST1–ST19）")
    A("加一条 BCR 记录，**不含血**。本审计的 `includes_blood` 一律以 sample title / series title")
    A("判定，不采信 series summary，正是为避免这一类误判。")
    A("")
    A(CAVEAT)
    A("")
    A("---")
    A("")
    A("## anchor #15 的额外核实（spec §5 第 15 项）")
    A("")
    A("Zenodo record **8275845**（`%s`）:" % zen["zenodo_doi"])
    A("")
    A("- 文件 `comb_CD56_CD16_NK_blood.h5ad`：`meta_tissue` 唯一取值 `Blood`，共 **%d** 个循环 NK 细胞。"
      % zen["total_cells_in_object"])
    A("- 其中 **Nasopharyngeal Carcinoma(NPC) = %d 个细胞 / %d 位病人 / %d 个样本**，"
      % (zen["NPC_cells"], zen["NPC_patients"], zen["NPC_samples"]))
    A("  patient ID 为 `%s`。" % "`, `".join(zen["NPC_patient_ids"]))
    A("- NPC 细胞的来源数据集：**%s**（该文 24 个瘤种确含 nasopharyngeal cancer，此为直接证据）。"
      % ", ".join(zen["NPC_source_datasets"]))
    A("- NPC 循环 NK 在图谱中的亚群分布：")
    A("")
    A("| NK subset | cells |")
    A("|---|---|")
    for k, v in zen["NPC_cells_by_NK_subset"].items():
        A("| %s | %d |" % (k, v))
    A("")
    A("全部瘤种在该 blood 对象中的占比：")
    A("")
    A("| histology | cells | patients |")
    A("|---|---|---|")
    for k, v in sorted(zen["all_histologies"].items(), key=lambda x: -x[1]["cells"]):
        A("| %s | %d | %d |" % (k, v["cells"], v["patients"]))
    A("")
    A("同一 record 的另一个文件 `comb_CD56_CD16_NK.h5ad`（全部组织，%d 个细胞）中，NPC 共 **%d 个细胞 / %d 位病人**，按组织拆开："
      % (zall["total_cells"], zall["NPC_cells"], zall["NPC_patients"]))
    A("")
    A("| meta_tissue | cells | patients |")
    A("|---|---|---|")
    for k, v in sorted(zall["NPC_by_tissue"].items(), key=lambda x: -x[1]["cells"]):
        A("| %s | %d | %d |" % (k, v["cells"], v["patients"]))
    A("")
    A("即：这份泛癌 NK 图谱里 NPC 的 NK 细胞**绝大多数来自血**（%d/%d），"
      % (zall["NPC_by_tissue"].get("Blood", {}).get("cells", 0), zall["NPC_cells"]))
    A("瘤内 NK 只有 %d 个、分散在 %d 位病人 —— 这与 scRNA-seq 对瘤内 NK 捕获率低的已知现象一致，"
      % (zall["NPC_by_tissue"].get("Tumor", {}).get("cells", 0),
         zall["NPC_by_tissue"].get("Tumor", {}).get("patients", 0)))
    A("也说明「图谱里已经有了」这句话对**血**成立、对**瘤内 NK** 只在极稀疏的意义上成立。")
    A("")
    A("（两个文件的 NPC 血细胞数相差 1 个细胞：blood 专用对象 %d，全组织对象 %d。"
      % (zen["NPC_cells"], zall["NPC_by_tissue"].get("Blood", {}).get("cells", 0)))
    A("此处如实记录，未做调和。）")
    A("")
    A(CAVEAT)

    with open(os.path.join(OUT, "OCCUPANCY.md"), "w") as f:
        f.write("\n".join(L) + "\n")

    json.dump({"qa": [ident(r) for r in qa], "qb": [ident(r) for r in qb],
               "qb_blood": [ident(r) for r in qb_blood], "qc": [ident(r) for r in qc],
               "qc_all": [ident(r) for r in qc_all], "qd": [ident(r) for r in qd],
               "qe": list(qe.keys()), "qf": [d["accession"] for d in qf]},
              open(os.path.join(RAW, "occupancy_answers.json"), "w"), indent=1)
    print("Q-A %d | Q-B %d (blood %d) | Q-C %d | Q-D %d | Q-E %d | Q-F %d"
          % (len(qa), len(qb), len(qb_blood), len(qc), len(qd), len(qe), len(qf)))


if __name__ == "__main__":
    main()
