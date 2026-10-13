// Word report for the five-organ IL18-high / CASP1-low screen.
// Reads results/answer_table.json (scripts/answer_pack.py). Every detection rate shown here is the
// POOLED rate (positive cells / analysed cells) with its raw counts; the donor median is reported in a
// separate column and explicitly labelled as a censoring statistic.
const fs = require("fs"), d = require("docx");
const { Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, BorderStyle, ShadingType, ImageRun } = d;
const ROOT = "/home/user/aaa/il18_casp5organ";
const R = JSON.parse(fs.readFileSync(ROOT + "/results/answer_table.json", "utf8"));
const F = "Microsoft YaHei", INK = "1A1A1A", MU = "666666", HD = "EEF1F5";
const A = "1B7A43", B = "1F5FA8", B2 = "2F6FB4", C = "B07A00", D = "9A2B2B", E = "6B6B6B";
const COL = { A_supported: A, B_potential: B, B_potential_screen_missed: B2, C_casp1_present: C, D_refuted: D, E_no_matched_reference: E };
const one = (sp, name) => R.find(r => r.sp === sp && r.name.startsWith(name));
const pct = v => (v >= 1 ? v.toFixed(1) : v.toFixed(2)) + "%";

const P = (t, o = {}) => new Paragraph({
  spacing: { after: o.after ?? 120, before: o.before ?? 0, line: 288 },
  children: [].concat(t).map(x => typeof x === "string"
    ? new TextRun({ text: x, font: F, size: o.size ?? 21, color: o.color ?? INK, bold: o.bold }) : x)
});
const H = (t, l = HeadingLevel.HEADING_1) => new Paragraph({
  heading: l, spacing: { before: 260, after: 120 },
  children: [new TextRun({ text: t, font: F, bold: true, size: l === HeadingLevel.HEADING_1 ? 28 : 23, color: INK })]
});
const run = (t, o = {}) => new TextRun({ text: t, font: F, size: o.size ?? 18, color: o.color ?? INK, bold: o.bold });

function tbl(hdr, rows, w, colorCol) {
  const tot = w.reduce((a, b) => a + b, 0);
  const hr = new TableRow({
    tableHeader: true, children: hdr.map((h, i) => new TableCell({
      width: { size: w[i], type: WidthType.DXA }, shading: { type: ShadingType.CLEAR, fill: HD },
      margins: { top: 40, bottom: 40, left: 60, right: 60 },
      children: [P(run(h, { bold: true, size: 15 }), { after: 0 })]
    }))
  });
  const br = rows.map(r => new TableRow({
    children: r.cells.map((c, i) => new TableCell({
      width: { size: w[i], type: WidthType.DXA }, margins: { top: 30, bottom: 30, left: 60, right: 60 },
      children: [P(run(String(c), { size: 14, color: (colorCol !== undefined && i === colorCol) ? (COL[r.tier] || INK) : INK, bold: colorCol === i }), { after: 0 })]
    }))
  }));
  return new Table({
    width: { size: tot, type: WidthType.DXA }, columnWidths: w, rows: [hr, ...br],
    borders: ["top", "bottom", "left", "right", "insideHorizontal", "insideVertical"]
      .reduce((o, k) => (o[k] = { style: BorderStyle.SINGLE, size: 2, color: "D8DBE0" }, o), {})
  });
}

const K = [];
K.push(new Paragraph({ spacing: { after: 60 }, children: [run("跨器官 IL18 高 / CASP1 低 细胞类型筛选", { bold: true, size: 34 })] }));
K.push(P("人与小鼠 · 肺 肝 胰 脑 脊髓 · 静息基线 · 修订版（参考层扩展 + 判读更正 + 检出率口径更正）", { color: MU, size: 20, after: 60 }));
K.push(P("CELLxGENE Census LTS 2025-11-08 · 以供体为统计单位 · 2026-10-13", { color: MU, size: 16, after: 200 }));

// ---- 0. what the percentage is a percentage of ---------------------------------------------------
K.push(H("〇、检出率的定义：百分比是什么的百分比"));
K.push(P([run("检出率 = 分子 ÷ 分母。", { bold: true }),
run("分母 = 该细胞类型在某一分层（数据集 × 器官 × 单细胞/单核）内、每位供体至少 20 个细胞后进入分析的细胞总数；分子 = 其中该基因原始计数 ≥ 1 的细胞个数。即「有多少比例的细胞里测得到这个基因」——不是表达量占比，不是相对最高值的百分比，不是读段占比。")]));
K.push(P([run("本版的更正：改用合并检出率，不再用供体中位数。", { bold: true, color: D }),
run("上一版图表标注的是「每位供体各算一次检出率、再取中位数」。每位供体只有数十到上百个该类型细胞，当真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，中位数被截断为 0——这是统计量把低值压成零，不是真的没有检出。本版全部改用合并检出率（全部供体的阳性细胞数 ÷ 全部分析细胞数），并同时标注阳性细胞数与有阳性的供体数；供体中位数单列保留，供对照。")]));
const ex = [["人·肺 肥大细胞", one("人", "肥大细胞")], ["小鼠·脑 小胶质/CNS巨噬", one("小鼠", "小胶质")],
["人·肺 纤毛上皮（单核）", one("人", "纤毛上皮")], ["小鼠·胰 胰导管细胞", one("小鼠", "胰导管")],
["小鼠·脑 星形胶质", one("小鼠", "星形胶质")], ["小鼠·脑 少突胶质", one("小鼠", "少突胶质")]];
K.push(tbl(["细胞类型", "CASP1 阳性细胞 / 分析细胞", "合并检出率", "供体中位数", "有阳性的供体"],
  ex.map(([n, r]) => ({ tier: "x", cells: [n, `${r.casp1_pos} / ${r.n_cells}`, pct(r.casp1), r.med_casp1.toFixed(2) + "%", `${r.casp1_donors} / ${r.n_donors}`] })),
  [3100, 2600, 1400, 1300, 1500]));
K.push(P("这六行在上一版里全部显示为 0.0%，而真实检出率相差一百倍以上：肥大细胞 2.60%（138 个阳性细胞、21/44 供体有阳性）是实实在在的信号；少突胶质 0.02%（8987 个细胞里 2 个）才是在充分样本量下的真·近零。因此「0%」在上一版既可能是真低，也可能只是统计量的截断，必须分开看。", { size: 17, color: MU }));
K.push(P([run("IL18:CASP1 参考比值", { bold: true }), run(" = 参考层中同一来源、同一单位、同一状态下 IL18 与 CASP1 的中位值之比；不同来源单位不可比，故只在来源内部计算。inf 表示该来源中 CASP1 处于其零底线。")]));

// ---- 1. what was corrected -----------------------------------------------------------------------
K.push(H("一、本次修订更正了什么"));
K.push(P([run("（1）判读实现错误。", { bold: true }), run("前一版的「零值判读」脚本从未打开过参考文件，也从未使用候选自身的 CASP1 值——判读结果是细胞谱系标签的纯函数。因此星形胶质被判为可信，而肺泡巨噬被判为「dropout 伪影」，而后者恰恰是全部候选中 CASP1 检出最高的一个。该脚本已废弃，由 evidence_integrate.py 取代：真正读取参考值与候选自身数值，按来源内部比值分级。")]));
K.push(P([run("（2）不再以「你的零不可信」否定检出。", { bold: true, color: D }), run("按规格 §9.5，CASP1 零值只做并列呈现、不做统计检验，也不删除候选。本版被排除的候选全部给出正面理由——CASP1 在该细胞中确实被检出（C 级），或方向相反（D 级）。")]));
K.push(P([run("（3）参考层扩展。", { bold: true }), run("从 3 个来源扩展到 Human Protein Atlas（含抗体 IHC 蛋白层面）、GTEx（含脊髓）、Monaco 2019 / Haemopedia、ImmGen 扩展面板、Allen Brain Map、brainrnaseq/Zhang、CELLxGENE WMG，以及本次新检索到的 FACS 分选数据集（分选胰岛 α/β/δ、腺泡/导管、纯化肝细胞、急性纯化人小胶质）。统一参考表共 1624 行，带匹配质量标注（exact / close / tissue_level / extrapolated），只有 exact 与 close 参与判读。扩展后，前一版因「无分选对照」而无法判读的五个实质细胞群全部变为可判读。")]));
K.push(P([run("（4）检出率口径更正。", { bold: true }), run("见上一节。")]));

// ---- 2. the answer -------------------------------------------------------------------------------
K.push(H("二、明确答案：按证据强度分级"));
const tierName = {
  A_supported: "A 多源一致支持", B_potential: "B 潜在可能", B_potential_screen_missed: "B 潜在（相对规则漏掉）",
  C_casp1_present: "C CASP1 确实存在", D_refuted: "D 方向相反，排除", E_no_matched_reference: "E 仍无匹配参考"
};
K.push(tbl(["证据级别", "物种·器官", "类型", "细胞类型", "IL18 检出（阳性/分析细胞）", "CASP1 检出（阳性/分析细胞）", "有 CASP1 阳性的供体", "参考比值", "参考源数"],
  R.map(r => ({
    tier: r.tier, cells: [tierName[r.tier], r.sp + "·" + r.organ, r.susp, r.name,
    `${pct(r.il18)}（${r.il18_pos}/${r.n_cells}）`, `${pct(r.casp1)}（${r.casp1_pos}/${r.n_cells}）`,
    `${r.casp1_donors}/${r.n_donors}`, (r.ratio || "—"), String(r.src)]
  })),
  [1650, 780, 560, 1950, 1500, 1550, 820, 700, 540], 0));
K.push(P("检出率为单细胞/单核数据的合并值（分母见第〇节）；参考比值为参考层同源同单位的 IL18:CASP1。最后一行 Cajal-Retzius 神经元无匹配参考，仅作记录。", { size: 15, color: MU }));

K.push(new Paragraph({ spacing: { before: 200, after: 100 }, children: [run("图：候选总览（左 IL18、右 CASP1，同一分母；CASP1 为对数轴）", { bold: true, size: 17 })] }));
K.push(new Paragraph({
  spacing: { after: 160 },
  children: [new ImageRun({ data: fs.readFileSync(ROOT + "/figures/answer_v2.png"), transformation: { width: 620, height: 322 } })]
}));

// ---- 3. tier by tier -----------------------------------------------------------------------------
const ast = one("小鼠", "星形胶质"), oli = one("小鼠", "少突胶质"), mast = one("人", "肥大细胞"),
  amac = one("小鼠", "肺泡巨噬"), mg = one("小鼠", "小胶质"), at1 = R.filter(r => r.name.startsWith("AT1")),
  dIslet = R.find(r => r.name.indexOf("分选胰岛") > -1), duct = one("小鼠", "胰导管");
K.push(H("三、逐级说明"));
K.push(P([run("A 级 — 小鼠脑 星形胶质、少突胶质。", { bold: true, color: A }),
run(`Allen Brain Map、CELLxGENE、brainrnaseq/Zhang 三个独立来源一致：Il18 明确存在而 Casp1 处于零底线（Zhang：星形胶质 Il18 115.1 vs Casp1 1.13 FPKM，约 102 倍；少突胶质 24.7 vs 0.12）。本数据中星形胶质 IL18 ${pct(ast.il18)}（${ast.il18_pos}/${ast.n_cells} 细胞、${ast.n_donors}/${ast.n_donors} 供体），CASP1 仅 ${ast.casp1_pos}/${ast.n_cells}；少突胶质 IL18 ${pct(oli.il18)}（${oli.il18_pos}/${oli.n_cells}），CASP1 ${oli.casp1_pos}/${oli.n_cells}。样本量足够大，近零是真的近零，不是统计截断。`)]));
K.push(P([run("B 级 — 人肺肥大细胞。", { bold: true, color: B }),
run(`HPA 单细胞（n=1237）IL18 331.5 nCPM 为人肺所有细胞类型中最高，CASP1 仅 11.1 nCPM（约 30 倍）。本数据 IL18 ${pct(mast.il18)}（${mast.il18_pos}/${mast.n_cells}），CASP1 ${pct(mast.casp1)}（${mast.casp1_pos}/${mast.n_cells}，${mast.casp1_donors}/${mast.n_donors} 供体有阳性）——CASP1 低但确有检出，另有来源显示 CASP1 可测到，故列为潜在而非确证。`)]));
K.push(P([run("B 级 — 小鼠胰岛 D(δ) 细胞、胰导管细胞。", { bold: true, color: B }),
run(`前一版因「无分选对照」判为无法判读；本次检索到 FACS 分选数据（GSE80673）：Il18 13.16 RPKM vs Casp1 0.016 RPKM。两套独立单细胞数据集同向：分选胰岛 δ 细胞 IL18 ${pct(dIslet.il18)}（${dIslet.il18_pos}/${dIslet.n_cells}）而 CASP1 仅 ${dIslet.casp1_pos}/${dIslet.n_cells}；胰导管细胞 IL18 ${pct(duct.il18)} 而 CASP1 ${duct.casp1_pos}/${duct.n_cells}。`)]));
K.push(P([run("B 级 — 人肺 AT1（被筛选的相对规则漏掉）。", { bold: true, color: B2 }),
run(`AT1 在人肺单细胞分层 IL18 排名 3/39、单核分层排名 1/23，HPA 给出 IL18 214.1 vs CASP1 7.9 nCPM（27 倍）。它未进候选，仅因该分层另有三分之一的细胞类型 CASP1 比它更低——这是「分层内最低三分位」这一相对排序规则的结构性盲点，不是证据不足。需注意其自身 CASP1 检出并不低（${at1.map(r => r.susp + " " + pct(r.casp1) + "，" + r.casp1_pos + "/" + r.n_cells).join("；")}），故列为需复核的潜在项而非确证。`)]));
K.push(P([run("C 级 — CASP1 确实存在，不满足「加工酶缺失」前提。", { bold: true, color: C }),
run(`小鼠肺泡巨噬 IL18 检出 ${pct(amac.il18)}（${amac.il18_pos}/${amac.n_cells}）为全部候选最高，但其 CASP1 同样是全部候选中最高的 ${pct(amac.casp1)}（${amac.casp1_pos}/${amac.n_cells}，${amac.casp1_donors}/${amac.n_donors} 供体全部有阳性），ImmGen 分选 Casp1 = 511（Il18 2384，比值仅 4.7）。小鼠脑小胶质 CASP1 ${pct(mg.casp1)}（${mg.casp1_pos}/${mg.n_cells}），参考比值 3.4。人肺纤毛上皮参考比值 11.9–13.3。它们被排除，是因为 CASP1 被检出，不是因为测不准。`)]));
K.push(P([run("D 级 — 方向相反。", { bold: true, color: D }),
run("人肝嗜碱粒细胞参考层 CASP1 60–234 nTPM 而 IL18 仅 4.7–14.8（比值 0.085）；人肺树突状细胞参考比值 2.0。两者在参考层中 CASP1 不低于或接近 IL18，与筛选前提方向相反。")]));

// ---- 4. species ----------------------------------------------------------------------------------
K.push(H("四、必须注明的物种差异"));
K.push(P([run("A 级结论为小鼠特异，不可外推到人。", { bold: true, color: D }),
run("HPA 人脑单核数据显示 IL18 压倒性集中于小胶质（338.5 nCPM），星形胶质仅 3.7、少突胶质 2.6 nCPM。即小鼠 CNS 胶质的「IL18 高 / CASP1 低」在人脑并不成立：人脑该表型属于小胶质，而人小胶质 CASP1 为 41.2 nCPM（比值 8.2），不满足前提。人脑目前没有对应的 A 级候选。")]));

// ---- 5. limits -----------------------------------------------------------------------------------
K.push(H("五、限制"));
["环境 RNA 未校正：本机无 GPU（CellBender 不可行），Census 源矩阵已过滤、无空液滴（SoupX 不可行），全部候选标 ambient_uncorrected=TRUE，为待验证结果。环境污染会抬高低表达细胞的 CASP1，使保留的近零值偏保守，但对低表达的神经/胶质候选影响最大。",
  "单核与单细胞不可混比；跨器官只比候选名单，不比表达量；跨物种不比表达量。",
  "三分位、CASP3 最低 20% 下限、以及本报告的 20 倍 / 3 倍参考比值分界，均为工程取值，无生物学依据；完整数值已输出，可自行重切。",
  "小鼠脊髓数据（GSE103892）仅有基因符号、无 Ensembl ID，按规格 §1.2 不进入候选判定；人脑成年正常单细胞、人脊髓多供体为已知数据缺口，未用邻近组织顶替。",
  "IL18 转录本 ≠ pro-IL-18 蛋白储量；CASP3 mRNA ≠ caspase-3 酶活性；静息基线不能回答刺激后 CASP1 能否诱导。"
].forEach(t => K.push(P("• " + t)));

K.push(H("六、建议的后续验证优先级"));
K.push(P("1. 小鼠脑星形胶质与少突胶质（A 级，多源一致）——蛋白层面检测 pro-IL-18 与 pro-caspase-1。"));
K.push(P("2. 小鼠胰岛 D(δ) 细胞与胰导管细胞、人肺肥大细胞、人肺 AT1（B 级，单一来源支持或被相对规则漏掉）——先用带空液滴的原始矩阵做环境 RNA 校正复核，再进蛋白层面。"));
K.push(P("3. 人脑：无 A 级候选；若要在人体系统中找同类表型，需另行筛选而非沿用小鼠结论。"));

K.push(H("附：交付物"));
K.push(P("answer_table.tsv / answer_table.json（本报告的答案表，含合并检出率与供体中位数两套口径）、il18_casp1_evidence.tsv（参考层整合判读）、screen_missed_but_reference_supported.tsv（被相对规则漏掉者）、reference/ext/unified_reference.tsv（统一参考表，含匹配质量标注）、candidates.tsv、il18_casp_summary.tsv、il18_casp_celltype_table.tsv、census_scan.tsv、threshold_sensitivity.tsv、celltype_map.tsv、datasets.yaml、run_log.md、README.md、figures/。", { size: 16, color: MU }));

const doc = new Document({
  styles: { default: { document: { run: { font: F } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1080, right: 1080 } } }, children: K }]
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(ROOT + "/IL18_CASP_五器官筛选报告.docx", b); console.log("wrote docx"); });
