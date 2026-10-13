// Word report for the five-organ IL18-high / CASP1-low screen.
// Reads results/answer_table.json (scripts/answer_pack.py). Every detection rate here is the POOLED rate
// (positive cells / analysed cells) with its raw counts; the donor median is a separate column and is
// labelled as a censoring statistic.
const fs = require("fs"), d = require("docx");
const { Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, BorderStyle, ShadingType, ImageRun } = d;
const ROOT = "/home/user/aaa/il18_casp5organ";
const R = JSON.parse(fs.readFileSync(ROOT + "/results/answer_table.json", "utf8"));
const F = "Microsoft YaHei", INK = "1A1A1A", MU = "666666", HD = "EEF1F5";
const A = "1B7A43", B = "1F5FA8", D = "9A2B2B";
const COL = { "确定": A, "可能": B, "相反": D };   // the only three classes shown to the reader
const one = (sp, name) => R.find(r => r.sp === sp && r.name.startsWith(name));
const pick = (sp, name, ds) => R.find(r => r.sp === sp && r.name.startsWith(name) && r.dataset_cn === ds);
const pct = v => (v >= 1 ? v.toFixed(1) : v.toFixed(2)) + "%";
const nV = v => R.filter(r => r.verdict === v).length;

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
K.push(P("人与小鼠 · 肺 肝 胰 脑 脊髓 · 静息基线 · 29 个分层（扩充批次 v0.2）", { color: MU, size: 20, after: 60 }));
K.push(P("CELLxGENE Census LTS 2025-11-08 · 以供体为统计单位 · 2026-10-13", { color: MU, size: 16, after: 200 }));

// ---- 0. definition --------------------------------------------------------------------------------
K.push(H("〇、检出率的定义：百分比是什么的百分比"));
K.push(P([run("检出率 = 分子 ÷ 分母。", { bold: true }),
run("分母 = 该细胞类型在某一分层（数据集 × 器官 × 单细胞/单核）内、每位供体至少 20 个细胞后进入分析的细胞总数；分子 = 其中该基因原始计数 ≥ 1 的细胞个数。即「有多少比例的细胞里测得到这个基因」——不是表达量占比，不是相对最高值的百分比，不是读段占比。")]));
K.push(P([run("本报告一律用合并检出率，不用供体中位数。", { bold: true, color: D }),
run("每位供体只有数十到上百个该类型细胞，当真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，中位数被截断为 0。例如人肺肥大细胞 CASP1 合并 138/5299 = 2.60%（21/44 供体有阳性）与小鼠脑少突胶质 CASP1 合并 2/8987 = 0.02%，在旧口径下同为 0.00%，真实值差一百倍以上。图中每行的灰色短竖线是 1/n，即该行可检出的最低非零值；一个落在该线上的 CASP1 点，或一个 1/n 很高的行里的「0」，都是样本量不足的信号。")]));

// ---- 1. what the expansion did --------------------------------------------------------------------
K.push(H("一、这一版扩充了什么，以及为什么必须扩充"));
K.push(P([run("扩充前，人侧除肺以外的器官分层太窄，三分位规则跑不动：", { bold: true }),
run("人脑只有 5 个细胞类型达到「≥3 供体」（整个分层只有 3 位供体），三分位 = 1.7 个，最高与最低不可能是同一个；人胰单细胞 0 个细胞类型合格（供体仅 3 位，且该数据集以外分泌为主，胰岛 δ 细胞一个都没有）；人脊髓只有 1 位供体。实测确认：把 IL18 检出率门槛从 10% 降到 5% 或 2%，人侧一个新候选都不会多出来——卡住的是分层宽度与 CASP1 相对排名，不是检出率门槛。")]));
K.push(P("分层数由 15 增至 29，全部经本次会话内的 Census obs 查询逐条核实（物种、器官、suspension_type、normal、成年、供体数、细胞类型数、关键细胞类型的细胞数），记录在 datasets.yaml 的 verified_how 字段。新增："));
K.push(tbl(["器官 / 物种", "新增数据", "供体", "细胞", "细胞类型", "补的是什么缺口"],
  [["人 胰腺", "3 套人胰岛 scRNA（multiomics / 48-donor / cross-species）", "19 / 17 / 5", "47k / 89k / 26k", "10 / 12 / 4", "此前完全没有人胰岛 δ(D) 细胞可查"],
  ["人 脑", "Allen Human Variation Study（10x 3' v3，suspension=cell，真单细胞）", "78", "379k", "18", "此前被判定为「不存在」的成年正常人脑 scRNA"],
  ["人 脑", "SEA-AD MTG / DLPFC、PsychAD MSSM（单核，仅取 disease=normal）", "46 / 44 / 261", "764k / 797k / 1050k", "18 / 18 / 25", "人脑单核此前只有 3 位供体"],
  ["人 脊髓", "颈段脊髓白质单核图谱", "19", "20k", "15", "此前只有 1 位供体，整个器官无结论"],
  ["人 肝", "肝 full / 肝 global", "14 / 7", "219k / 13k", "63 / 34", "此前只有 11 个可分析细胞类型"],
  ["人 肺", "CellRef / Distal", "21 / 12", "122k / 72k", "46 / 77", "补细胞类型更全的两套"],
  ["小鼠 肺 / 肝", "Tabula Muris Senis Smart-seq2（FACS）", "14 / 9", "5k / 3k", "32 / 11", "补一条与 10x 互相独立的技术路线"]].map(c => ({ tier: "x", cells: c })),
  [1250, 3000, 1050, 1500, 1000, 2350]));
K.push(P([run("一处必须说明的例外：", { bold: true }), run("CELLxGENE 把 Allen HVS 这一个研究按细胞类型拆成 24 个 dataset_id 分发，每个只含 1 类细胞。若按 dataset_id 分层，每层只有一个细胞类型，分层内排序无从谈起。故把同一 collection 的 24 个 id 合并为一个分层——这是重组同一研究（同 collection、同供体、同 assay、同组织），不是跨数据集整合。该例外写在 datasets.yaml 与代码注释里。")], { size: 18 }));

// ---- 2. the answer --------------------------------------------------------------------------------
K.push(H("二、明确答案：确定 / 可能 / 相反"));
K.push(P([run("只分三类：", { bold: true }),
run(`确定 = ≥2 个独立参考来源一致显示 CASP1 处于零底线或比 IL18 低 ≥20 倍（${nV("确定")} 行）；可能 = 只有单一参考来源支持，或本数据一致但尚无匹配参考（${nV("可能")} 行）；相反 = 参考层显示 CASP1 确实存在，或方向与「CASP1 低」这一前提相反（${nV("相反")} 行）。末列「细分依据」保留六档内部判读，仅作留痕，不改变上述三分。`)]));
K.push(P([run("「×n」一列是重复性：", { bold: true }),
run("同一细胞类型在几个相互独立的分层（不同数据集、不同供体）里都成为候选。这是与参考层相互独立的另一条证据轴——在四套数据里重复出现，比一套数据里的注释偏差或批次效应更难解释掉。但重复不等于参考层支持，二者分开记。")], { size: 18 }));
K.push(tbl(["判定", "物种·器官", "细胞类型", "数据集", "类型", "IL18 检出（阳性/分析细胞）", "CASP1 检出（阳性/分析细胞）", "CASP1 阳性供体", "×n", "参考比值", "细分依据"],
  R.map(r => ({
    tier: r.verdict, cells: [r.verdict, r.sp + "·" + r.organ, r.name, r.dataset_cn, r.susp,
    `${pct(r.il18)}（${r.il18_pos}/${r.n_cells}）`, `${pct(r.casp1)}（${r.casp1_pos}/${r.n_cells}）`,
    `${r.casp1_donors}/${r.n_donors}`, r.n_strata > 1 ? "×" + r.n_strata : "", (r.ratio || "—"), r.tier_cn]
  })),
  [620, 700, 1500, 1200, 500, 1380, 1400, 760, 380, 620, 1380], 0));
K.push(P("检出率为单细胞/单核数据的合并值（分母见第〇节）；参考比值为参考层同源同单位的 IL18:CASP1，inf 表示该来源中 CASP1 处于零底线，— 表示无匹配参考。", { size: 15, color: MU }));

K.push(new Paragraph({ spacing: { before: 200, after: 100 }, children: [run("图：候选总览——每行一个细胞类型，实心点 = IL18、空心点 = CASP1，同一条对数检出率轴；灰色短竖线 = 该行可检出的最低非零值 1/n；人与小鼠分开绘制", { bold: true, size: 17 })] }));
K.push(new Paragraph({
  spacing: { after: 160 },
  children: [new ImageRun({ data: fs.readFileSync(ROOT + "/figures/answer_v2.png"), transformation: { width: 620, height: 484 } })]
}));

// ---- 3. verdict by verdict ------------------------------------------------------------------------
const ast = one("小鼠", "星形胶质"), oli = one("小鼠", "少突胶质"), amac = one("小鼠", "肺泡巨噬"),
  pv = pick("人", "PV 中间神经元", "PsychAD"), pvh = pick("人", "PV 中间神经元", "HVS 皮层"),
  l5 = pick("人", "L5 ET", "HVS 皮层"), mast = pick("人", "肥大细胞", "HLCA"),
  nk = one("人", "CD56bright");
K.push(H("三、三类逐条说明"));

K.push(H("确定", HeadingLevel.HEADING_2));
K.push(P([run("小鼠脑 星形胶质、少突胶质。", { bold: true, color: A }),
run(`Allen Brain Map、CELLxGENE、brainrnaseq/Zhang 三个独立来源一致：Il18 明确存在而 Casp1 处于零底线（Zhang：星形胶质 Il18 115.1 vs Casp1 1.13 FPKM，约 102 倍；少突胶质 24.7 vs 0.12）。本数据星形胶质 IL18 ${pct(ast.il18)}（${ast.il18_pos}/${ast.n_cells} 细胞、${ast.n_donors} 位供体全部有阳性），CASP1 仅 ${ast.casp1_pos}/${ast.n_cells}；少突胶质 IL18 ${pct(oli.il18)}（${oli.il18_pos}/${oli.n_cells}），CASP1 ${oli.casp1_pos}/${oli.n_cells}。样本量足够大，近零是真的近零，不是统计截断。`)]));

K.push(H("可能", HeadingLevel.HEADING_2));
K.push(P([run("人脑皮层 PV 中间神经元、L5 ET 投射神经元（本次扩充的最大发现）。", { bold: true, color: B }),
run(`PV 中间神经元在 Allen HVS 真单细胞数据（${pvh.n_donors} 供体）IL18 ${pct(pvh.il18)}、CASP1 ${pvh.casp1_pos}/${pvh.n_cells} = ${pct(pvh.casp1)}，在 PsychAD 单核数据（${pv.n_donors} 供体）IL18 ${pct(pv.il18)}、CASP1 ${pv.casp1_pos}/${pv.n_cells} = ${pct(pv.casp1)}；L5 ET 投射神经元在 HVS、SEA-AD MTG、SEA-AD DLPFC 三套数据集里重复，IL18 31.8%–${pct(l5.il18)}、CASP1 0.00%–0.22%。单细胞与单核两类技术、四套独立数据集、合计数百位供体同向。列为「可能」而非「确定」，是因为参考层对人 CNS 只解析到「神经元」这一级（没有按中间神经元亚类分选的 bulk），而该级参考中 HPA 给出 CASP1 = 0（比值 inf），Darmanis 2015 单细胞与 Zhang FPKM 两个来源却显示 CASP1 与 IL18 同量级。参考层本身不一致，故不升为确定。`)]));
K.push(P([run("人肺 AT1 肺泡上皮。", { bold: true, color: B }),
run("在 HLCA、CellRef、肺单核三套分层里都是 IL18 高而被同一条相对规则漏掉（分层内 IL18 排名 3/39、3/30、1/23），HPA 给出 IL18 214.1 vs CASP1 7.9 nCPM（27 倍）。但其自身 CASP1 检出 1.85%–5.49% 并不算低，故为需复核的可能项。")]));
K.push(P([run("人肺 肥大细胞。", { bold: true, color: B }),
run(`HPA 单细胞 IL18 331.5 nCPM 为人肺最高、CASP1 仅 11.1 nCPM（约 30 倍）。本数据 HLCA 分层 IL18 ${pct(mast.il18)}、CASP1 ${mast.casp1_pos}/${mast.n_cells} = ${pct(mast.casp1)}，CellRef 分层同向。单一参考来源支持。`)]));
K.push(P([run("人胰腺 腺泡细胞。", { bold: true, color: B }),
run("在两套人胰岛数据集里 IL18 检出 10.8% 与 30.3%、CASP1 0.25% 与 3.45%，参考比值 43.2；被相对规则漏掉。注意这两套数据以胰岛为主，腺泡细胞是伴随采到的，细胞数与注释质量不如专门的胰腺图谱。")]));
K.push(P([run("小鼠 胰岛 D(δ) 细胞、胰导管细胞。", { bold: true, color: B }),
run("FACS 分选数据（GSE80673）Il18 13.16 RPKM vs Casp1 0.016 RPKM；两套独立小鼠数据集同向。但见第四节——人侧不支持这一发现。")]));

K.push(H("相反", HeadingLevel.HEADING_2));
K.push(P([run("CASP1 确实存在，不满足「加工酶缺失」前提。", { bold: true, color: D }),
run(`小鼠肺泡巨噬 IL18 ${pct(amac.il18)}（${amac.il18_pos}/${amac.n_cells}）为全部最高，但其 CASP1 同样是全部中最高的 ${pct(amac.casp1)}（${amac.casp1_pos}/${amac.n_cells}，${amac.casp1_donors}/${amac.n_donors} 供体全部有阳性），ImmGen 分选 Casp1 = 511（Il18 2384，比值仅 4.7）。小鼠脑小胶质参考比值 3.4；人肺气道/纤毛上皮 11.9。它们被排除是因为 CASP1 被检出，不是因为测不准。`)]));
K.push(P([run("方向相反。", { bold: true, color: D }),
run(`人肝嗜碱粒细胞参考层 CASP1 60–234 nTPM 而 IL18 仅 4.7–14.8（比值 0.085）；人肺树突状细胞参考比值 2.0。人肝 CD56bright NK 细胞没有匹配参考，但其自身数据就是反的——CASP1 检出 ${pct(nk.casp1)}（${nk.casp1_pos}/${nk.n_cells}）高于 IL18 的 ${pct(nk.il18)}，故按自身数据判为相反而非可能。`)]));

// ---- 4. what the expansion answered ---------------------------------------------------------------
K.push(H("四、扩充数据回答了四个此前无法回答的问题"));
K.push(P([run("1. 人脑里 IL18 高的是神经元，不是胶质。", { bold: true }),
run("扩充前人脑只有 Siletti MTG 一套（3 位供体、5 个可分析细胞类型），什么也判不出。现在四套数据集一致指向 PV 中间神经元与 L5 ET 投射神经元。")]));
K.push(P([run("2. 小鼠的星形胶质/少突胶质结论是真实物种差异，不是单核测序灵敏度问题。", { bold: true, color: D }),
run("这是上一版明确标为「悬而未决」的问题。现在有了同技术的对照：同为真单细胞 scRNA，人皮层星形胶质 IL18 检出 4.1%（分层内排 16/18，78 位供体）、少突胶质 1.4%（17/18）；小鼠皮层星形胶质 55.4%（1/24，40 位供体）、少突胶质 27.4%（2/24）。差 13–20 倍，技术相同、供体数都大。小鼠的「确定」结论不可外推到人。")]));
K.push(P([run("3. 人胰岛 δ 细胞不支持小鼠的发现。", { bold: true }),
run("新增三套人胰岛 scRNA（19/17/5 供体，δ 细胞 1217/4677/2889 个），δ 细胞 IL18 检出率均为 0.0–0.5%。这是明确的阴性结论；扩充前人侧根本没有 δ 细胞可查，这个问题无法回答。")]));
K.push(P([run("4. 人脊髓从「无结论」变成「明确阴性」。", { bold: true }),
run("新增颈段脊髓白质单核图谱（19 供体、15 个细胞类型，此前只有 1 位供体）。IL18 几乎只在小胶质（检出 52.7%），而小胶质 CASP1 检出 10.5%、分层内排第 2；其余细胞类型 IL18 检出率 ≤3.3%。该器官无候选——这是结论，不是数据缺口。注意该数据集只有白质，不含灰质与运动神经元。")]));

// ---- 5. limits ------------------------------------------------------------------------------------
K.push(H("五、限制"));
["环境 RNA 未校正：本机无 GPU（CellBender 不可行），Census 源矩阵已过滤、无空液滴（SoupX 不可行），全部候选标 ambient_uncorrected=TRUE，为待验证结果。环境污染会抬高低表达细胞的 CASP1，使保留的近零值偏保守，但对低表达的神经元/胶质候选影响最大。",
  "单核与单细胞不可混比；跨器官只比候选名单，不比表达量；跨物种不比表达量。人脑的「可能」同时来自单细胞（HVS）与单核（SEA-AD、PsychAD）两类数据，二者各自成立，但不能当作互相验证。",
  "人脑神经元的参考层只解析到「CNS 神经元」这一级，没有按中间神经元亚类分选的 bulk；该级参考内部不一致（HPA 给 CASP1 = 0，Darmanis 2015 与 Zhang 给同量级），这正是它停在「可能」的原因。",
  "三分位、CASP3 最低 20% 下限、以及本报告的 20 倍 / 3 倍参考比值分界，均为工程取值，无生物学依据；完整数值已输出，可自行重切。",
  "SEA-AD 与 PsychAD 的「normal」是各自队列里的对照供体，取自以神经退行性与精神疾病为研究对象的队列；年龄偏大，不等同于健康年轻人群。",
  "小鼠脊髓数据（GSE103892）仅有基因符号、无 Ensembl ID，按 §1.2 不进入候选判定；Census 中仍无小鼠脊髓数据。",
  "IL18 转录本 ≠ pro-IL-18 蛋白储量；CASP3 mRNA ≠ caspase-3 酶活性；静息基线不能回答刺激后 CASP1 能否诱导。"
].forEach(t => K.push(P("• " + t)));

K.push(H("六、建议的后续验证优先级"));
K.push(P("1. 小鼠脑星形胶质与少突胶质（确定，多源一致）——蛋白层面检测 pro-IL-18 与 pro-caspase-1。"));
K.push(P("2. 人脑皮层 PV 中间神经元与 L5 ET 投射神经元（可能，四套独立数据集重复）——先补按神经元亚类分选的参考，这是把它们升为「确定」唯一缺的一环，再进蛋白层面。"));
K.push(P("3. 人肺 AT1 与肥大细胞、人胰腺腺泡细胞（可能）——先用带空液滴的原始矩阵做环境 RNA 校正复核。"));
K.push(P("4. 小鼠胰岛 D(δ) 细胞（可能，但人侧不支持）——若目标是人体系统，优先级应低于上面三项。"));

K.push(H("附：交付物"));
K.push(P("answer_table.tsv / answer_table.json（答案表，含合并检出率、供体中位数、重复数据集数）、il18_casp1_evidence.tsv（参考层整合判读）、screen_missed_but_reference_supported.tsv（被相对规则漏掉者）、reference/ext/unified_reference.tsv（统一参考表）、candidates.tsv、il18_casp_summary.tsv、il18_casp_celltype_table.tsv、census_scan.tsv、threshold_sensitivity.tsv、celltype_map.tsv（702 行逐条映射）、datasets.yaml（29 个分层，含 verified_how）、run_log.md、README.md、figures/。", { size: 16, color: MU }));

const doc = new Document({
  styles: { default: { document: { run: { font: F } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1080, right: 1080 } } }, children: K }]
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(ROOT + "/IL18_CASP_五器官筛选报告.docx", b); console.log("wrote docx"); });
