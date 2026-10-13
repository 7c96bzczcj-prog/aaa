// Deck for the five-organ IL18-high / CASP1-low screen. Minimal: no decoration, answer first.
// Reads results/answer_table.json (scripts/answer_pack.py) so every number on a slide is the pooled
// detection rate with its raw counts, never the donor median (which censors low rates to 0).
const pptxgen = require("pptxgenjs"), fs = require("fs");
const ROOT = "/home/user/aaa/il18_casp5organ", FIG = ROOT + "/figures";
const R = JSON.parse(fs.readFileSync(ROOT + "/results/answer_table.json", "utf8"));
const by = v => R.filter(r => r.verdict === v);          // 确定 / 可能 / 相反 — the only split shown
const uniq = a => a.filter((x, i) => a.indexOf(x) === i);
const names = (v, sp) => uniq(by(v).filter(r => !sp || r.sp === sp).map(r => r.sp + r.organ + " " + r.name)).join(" · ");
const one = (sp, name) => R.find(r => r.sp === sp && r.name.startsWith(name));

const p = new pptxgen(); p.layout = "LAYOUT_WIDE";
const W = 13.33, H = 7.5, F = "Microsoft YaHei";
const INK = "1A1A1A", MU = "6B6B6B", FA = "9AA0A6";
const A = "1B7A43", B = "1F5FA8", D = "9A2B2B";
const fit = (iw, ih, bw, bh) => { const r = Math.min(bw / iw, bh / ih); return { w: iw * r, h: ih * r }; };
const head = (s, t) => s.addText(t, { isTextBox: true, x: .7, y: .45, w: 12, h: .7, fontFace: F, fontSize: 24, color: INK, bold: true });
const cap = (s, t) => s.addText(t, { isTextBox: true, x: .7, y: H - .55, w: 12, h: .38, fontFace: F, fontSize: 10.5, color: MU });
const fig = (s, path, iw, ih, top, bot) => { const d = fit(iw, ih, W - 1.4, bot - top); s.addImage({ path, x: (W - d.w) / 2, y: top + (bot - top - d.h) / 2, w: d.w, h: d.h }); };
const pct = v => (v >= 1 ? v.toFixed(1) : v.toFixed(2)) + "%";

// ---- 1. answer -----------------------------------------------------------------------------------
let s = p.addSlide();
s.addText("哪些细胞 IL18 高且 CASP1 低", { isTextBox: true, x: .7, y: 1.75, w: 12, h: 1.0, fontFace: F, fontSize: 33, color: INK, bold: true, valign: "middle" });
s.addText("人与小鼠 · 肺 肝 胰 脑 脊髓 · 静息基线 · 29 个分层 · 以供体为统计单位", { isTextBox: true, x: .7, y: 2.75, w: 12, h: .5, fontFace: F, fontSize: 16, color: MU });
s.addText([
  { text: "答案分三类。", options: { bold: true } },
  { text: "确定：", options: { bold: true, color: A } },
  { text: "小鼠脑星形胶质、少突胶质。　", options: {} },
  { text: "可能：", options: { bold: true, color: B } },
  { text: "人脑皮层 PV 中间神经元与 L5 ET 投射神经元（3–4 套独立数据集重复）、人肺 AT1 与肥大细胞、人胰腺腺泡细胞、小鼠胰岛 D(δ) 细胞与胰导管细胞。　", options: {} },
  { text: "相反：", options: { bold: true, color: D } },
  { text: "小鼠肺泡巨噬与小胶质、人肺气道/纤毛上皮、人肝嗜碱粒细胞与 CD56bright NK、人肺树突状细胞——CASP1 在它们之中确实被检出。", options: {} },
], { isTextBox: true, x: .7, y: 3.85, w: 12, h: 2.0, fontFace: F, fontSize: 16.5, color: INK, lineSpacingMultiple: 1.3 });
s.addText("检出率 = 原始计数 ≥1 的细胞数 ÷ 分析细胞数（合并全部供体），每个数字都标注阳性细胞数。",
  { isTextBox: true, x: .7, y: 6.2, w: 12, h: .4, fontFace: F, fontSize: 12.5, color: FA });

// ---- 2. the three verdicts -----------------------------------------------------------------------
s = p.addSlide(); head(s, "三类答案");
const amac = one("小鼠", "肺泡巨噬"), ast = one("小鼠", "星形胶质"), oli = one("小鼠", "少突胶质");
const rows = [
  ["确定", A, "≥2 个独立参考来源一致", "小鼠脑 星形胶质 · 小鼠脑 少突胶质",
    "Allen Brain Map + CELLxGENE + brainrnaseq/Zhang 三源一致：Il18 明确存在而 Casp1 处于零底线（星形胶质 115.1 vs 1.13 FPKM）。本数据 Casp1 阳性细胞仅 "
    + ast.casp1_pos + "/" + ast.n_cells + " 与 " + oli.casp1_pos + "/" + oli.n_cells + "，样本量足够大，近零是真的近零。"],
  ["可能", B, "单一来源支持，或本数据一致但尚无匹配参考", names("可能", "人") + " · " + names("可能", "小鼠"),
    "人脑皮层 PV 中间神经元（HVS 单细胞 78 供体 + PsychAD 单核 234 供体）与 L5 ET 投射神经元（HVS + SEA-AD MTG + SEA-AD DLPFC）CASP1 合并检出 0.07–0.22%，在 3–4 套独立数据集里重复；"
    + "人肺 AT1 在 3 套数据集被同一条相对规则漏掉；小鼠胰岛 D(δ) 有 FACS 分选证据（13.16 vs 0.016 RPKM）。"],
  ["相反", D, "参考层显示 CASP1 确实存在，或方向与前提相反", names("相反", "人") + " · " + names("相反", "小鼠"),
    "小鼠肺泡巨噬 IL18 " + pct(amac.il18) + " 为全部最高，但 CASP1 也是最高的 " + pct(amac.casp1) + "（" + amac.casp1_pos + "/" + amac.n_cells + "，14/14 供体），ImmGen 分选 Casp1=511；"
    + "人肝 CD56bright NK 的 CASP1 检出 25.2% 反而高于 IL18 的 19.1%。排除理由是 CASP1 被检出，不是测不准。"],
];
let y = 1.5;
rows.forEach(([tag, col, def, who, why]) => {
  s.addText([{ text: tag + "\n", options: { fontSize: 19, bold: true } }, { text: def, options: { fontSize: 10.5, color: MU } }],
    { isTextBox: true, x: .7, y, w: 2.9, h: 1.6, fontFace: F, color: col, valign: "top", lineSpacingMultiple: 1.15 });
  s.addText([{ text: who + "\n", options: { bold: true } }, { text: why, options: { color: MU } }],
    { isTextBox: true, x: 3.85, y, w: 8.8, h: 1.85, fontFace: F, fontSize: 11.5, color: INK, valign: "top", lineSpacingMultiple: 1.15 });
  y += 1.87;
});

// ---- 3-4. per-species figures --------------------------------------------------------------------
s = p.addSlide(); head(s, "人：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_human.png`, 2890, 1550, 1.15, 7.0);
cap(s, "每行一个细胞类型：实心点 = IL18，空心点 = CASP1，两点距离即差距；同一细胞类型在多个数据集重复时取供体最多的一套，右侧标 ×n。");

s = p.addSlide(); head(s, "小鼠：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_mouse.png`, 2890, 908, 1.4, 6.3);
cap(s, "人与小鼠分开绘制，不作跨物种定量比较。");

// ---- 5. what the expansion changed ---------------------------------------------------------------
s = p.addSlide(); head(s, "扩充数据后，四个此前无法回答的问题有了答案");
[["人脑：IL18 高的是神经元，不是胶质",
  "新增 Allen HVS（真单细胞 scRNA，78 供体）与 SEA-AD MTG / DLPFC、PsychAD（单核，46 / 44 / 261 供体）。PV 中间神经元与 L5 ET 投射神经元 IL18 检出 21–43%、CASP1 0.07–0.22%，在 3–4 套独立数据集重复。"],
 ["小鼠的胶质结论是真实物种差异，不是单核测序灵敏度问题",
  "同为真单细胞 scRNA：人皮层星形胶质 IL18 检出 4.1%（分层内排 16/18，78 供体）、少突胶质 1.4%（17/18）；小鼠 55.4%（1/24）与 27.4%（2/24）。技术相同、供体数都大，差 13–20 倍。"],
 ["人胰岛 δ 细胞不支持小鼠的发现",
  "新增 3 套人胰岛 scRNA（19 / 17 / 5 供体，δ 细胞 1217 / 4677 / 2889 个）。三套数据里 δ 细胞 IL18 检出率均为 0.0–0.5%。这是明确的阴性，而扩充前人侧根本没有 δ 细胞可查。"],
 ["人脊髓：从无结论变成明确阴性",
  "新增颈段脊髓白质单核图谱（19 供体 15 类，此前只有 1 位供体）。IL18 几乎只在小胶质（检出 52.7%），而小胶质 CASP1 检出 10.5%、分层内排第 2；其余细胞类型 IL18 ≤3.3%。该器官无候选。"],
].forEach(([k, v], i) => {
  s.addText(k, { isTextBox: true, x: .8, y: 1.4 + i * 1.33, w: 3.5, h: 1.2, fontFace: F, fontSize: 13, color: INK, bold: true, valign: "top", lineSpacingMultiple: 1.15 });
  s.addText(v, { isTextBox: true, x: 4.5, y: 1.4 + i * 1.33, w: 8.1, h: 1.25, fontFace: F, fontSize: 12, color: MU, valign: "top", lineSpacingMultiple: 1.2 });
});

// ---- 6. the detection-rate definition ------------------------------------------------------------
s = p.addSlide(); head(s, "检出率的口径：百分比是什么的百分比");
s.addText([{ text: "检出率 = 该基因原始计数 ≥1 的细胞数 ÷ 分析细胞数（全部供体合并）。", options: { bold: true } },
{ text: "不是表达量占比，不是相对最高值的百分比，不是读段占比。本版不用「供体中位数」：每位供体只有数十到上百个该类型细胞，真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，中位数被截断为 0。", options: {} }],
  { isTextBox: true, x: .8, y: 1.45, w: 11.7, h: 1.0, fontFace: F, fontSize: 14.5, color: INK, lineSpacingMultiple: 1.25 });
const ex = [["人·脑 PV 中间神经元（PsychAD）", one("人", "PV 中间神经元")],
["人·肺 肥大细胞（HLCA）", R.find(r => r.name === "肥大细胞" && r.dataset_cn === "HLCA")],
["小鼠·脑 星形胶质", ast], ["小鼠·脑 少突胶质", oli]];
s.addTable(
  [[{ text: "细胞类型" }, { text: "CASP1 阳性细胞 / 分析细胞" }, { text: "合并检出率" }, { text: "供体中位数" }, { text: "有阳性的供体" }]
    .map(c => ({ ...c, options: { bold: true, fill: "EEF1F5" } }))].concat(
      ex.filter(e => e[1]).map(([n, r]) => [n, `${r.casp1_pos} / ${r.n_cells}`, pct(r.casp1), r.med_casp1.toFixed(2) + "%", `${r.casp1_donors} / ${r.n_donors}`])),
  { x: .8, y: 2.75, w: 11.7, colW: [3.6, 2.9, 1.8, 1.7, 1.7], fontFace: F, fontSize: 12.5, color: INK,
    border: { pt: .5, color: "D8DBE0" }, rowH: .42, valign: "middle" });
s.addText("同样显示「0.0%」的两行真实检出率可以差一百倍：肥大细胞 2.60%（138 个阳性细胞、21/44 供体）是真有信号；少突胶质 0.02%（8987 个细胞里 2 个）才是在充分样本下的真·近零。",
  { isTextBox: true, x: .8, y: 5.1, w: 11.7, h: .9, fontFace: F, fontSize: 13.5, color: INK, lineSpacingMultiple: 1.2 });
cap(s, "分母 = 该细胞类型在该分层内、每位供体 ≥20 个细胞后进入分析的细胞总数；图中灰色短竖线为该行可检出的最低非零值 1/n。");

// ---- 7. mouse brain dot plot ---------------------------------------------------------------------
s = p.addSlide(); head(s, "小鼠脑（单细胞）：胶质 IL18 可检、CASP1 近乎测不到");
fig(s, `${FIG}/mouse/brain_cell_dots.png`, 2312, 1122, 1.22, 6.92);
cap(s, "星形胶质 IL18 56.7% 而 CASP1 仅 2/3766；小胶质 IL18 20.5% 但 CASP1 15/818 = 1.83%，是全脑最高的一档。");

// ---- 8. limits -----------------------------------------------------------------------------------
s = p.addSlide(); head(s, "限制");
["环境 RNA 未校正：无 GPU、Census 源矩阵无空液滴；全部候选为待验证，对低表达的神经元/胶质候选影响最大。",
  "单核与单细胞不可混比；跨器官只比候选名单，跨物种不比表达量。人脑的「可能」来自单细胞与单核两类数据，二者各自成立，不互为验证。",
  "人脑神经元的参考层只到「CNS 神经元」这一级，没有按中间神经元亚类分选的 bulk；该级参考里另有两个来源显示 CASP1 与 IL18 同量级，故为「可能」而非「确定」。",
  "三分位、CASP3 最低 20%、20 倍 / 3 倍比值分界均为工程取值，无生物学依据；完整数值已输出，可自行重切。",
  "IL18 转录本 ≠ pro-IL-18 蛋白储量；CASP3 mRNA ≠ caspase-3 酶活性；静息基线不回答刺激后 CASP1 能否诱导。"
].forEach((t, i) => s.addText(t, { isTextBox: true, x: .8, y: 1.45 + i * .98, w: 11.8, h: .9, fontFace: F, fontSize: 13.5, color: INK, lineSpacingMultiple: 1.15, valign: "top" }));

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
