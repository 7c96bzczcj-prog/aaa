// Deck for the five-organ IL18-high / CASP1-low screen. Minimal: no decoration, answer first.
// Reads results/answer_table.json (produced by scripts/answer_pack.py) so every number on a slide is
// the pooled detection rate with its raw counts, never the donor median (which censors low rates to 0).
const pptxgen = require("pptxgenjs"), fs = require("fs");
const ROOT = "/home/user/aaa/il18_casp5organ", FIG = ROOT + "/figures";
const R = JSON.parse(fs.readFileSync(ROOT + "/results/answer_table.json", "utf8"));
const by = t => R.filter(r => r.tier === t);
const one = (sp, name) => R.find(r => r.sp === sp && r.name.startsWith(name));

const p = new pptxgen(); p.layout = "LAYOUT_WIDE";
const W = 13.33, H = 7.5, F = "Microsoft YaHei";
const INK = "1A1A1A", MU = "6B6B6B", FA = "9AA0A6";
const A = "1B7A43", B = "1F5FA8", B2 = "4C87C4", C = "B07A00", D = "9A2B2B", E = "8A8A8A";
const fit = (iw, ih, bw, bh) => { const r = Math.min(bw / iw, bh / ih); return { w: iw * r, h: ih * r }; };
const head = (s, t) => s.addText(t, { isTextBox: true, x: .7, y: .5, w: 12, h: .7, fontFace: F, fontSize: 25, color: INK, bold: true });
const cap = (s, t) => s.addText(t, { isTextBox: true, x: .7, y: H - .55, w: 12, h: .38, fontFace: F, fontSize: 10.5, color: MU });
const fig = (s, path, iw, ih, top, bot) => { const d = fit(iw, ih, W - 1.4, bot - top); s.addImage({ path, x: (W - d.w) / 2, y: top + (bot - top - d.h) / 2, w: d.w, h: d.h }); };
const pct = v => (v >= 1 ? v.toFixed(1) : v.toFixed(2)) + "%";

// ---- 1. answer -----------------------------------------------------------------------------------
let s = p.addSlide();
s.addText("哪些细胞 IL18 高且 CASP1 低", { isTextBox: true, x: .7, y: 1.9, w: 12, h: 1.0, fontFace: F, fontSize: 33, color: INK, bold: true, valign: "middle" });
s.addText("人与小鼠 · 肺 肝 胰 脑 脊髓 · 静息基线 · 以供体为统计单位", { isTextBox: true, x: .7, y: 2.9, w: 12, h: .5, fontFace: F, fontSize: 16, color: MU });
s.addText([
  { text: "答案：", options: { bold: true } },
  { text: "多源参考一致支持的是 ", options: {} },
  { text: "小鼠脑星形胶质与少突胶质", options: { bold: true, color: A } },
  { text: "；有单一来源支持、值得优先验证的潜在可能是 ", options: {} },
  { text: "人肺肥大细胞、小鼠胰岛 D(δ) 细胞、小鼠胰导管细胞，以及被相对规则漏掉的人肺 AT1", options: { bold: true, color: B } },
  { text: "。肺泡巨噬、小胶质等不入选，是因为 CASP1 在它们之中确实被检出（肺泡巨噬 178/1341 = 13.3%），不是因为测不准。", options: {} },
], { isTextBox: true, x: .7, y: 4.15, w: 12, h: 1.6, fontFace: F, fontSize: 17, color: INK, lineSpacingMultiple: 1.3 });
s.addText("检出率 = 原始计数 ≥1 的细胞数 ÷ 分析细胞数（合并全部供体）；本版全部标注阳性细胞数与阳性供体数。",
  { isTextBox: true, x: .7, y: 6.15, w: 12, h: .4, fontFace: F, fontSize: 12.5, color: FA });

// ---- 2. tiers ------------------------------------------------------------------------------------
s = p.addSlide(); head(s, "分级答案");
const amac = one("小鼠", "肺泡巨噬"), mast = one("人", "肥大细胞"), at1 = R.filter(r => r.name.startsWith("AT1"));
const rows = [
  ["A 多源一致支持", A, by("A_supported").map(r => r.sp + r.organ + " " + r.name).join(" · "),
    "Allen Brain Map + CELLxGENE + brainrnaseq/Zhang 三源一致：Il18 明确存在而 Casp1 处于零底线（星形胶质 115.1 vs 1.13 FPKM）。本数据中 Casp1 阳性细胞 2/3766 与 2/8987。"],
  ["B 潜在可能", B, by("B_potential").map(r => r.sp + r.organ + " " + r.name.replace(/（.*/, "")).filter((v, i, a) => a.indexOf(v) === i).join(" · "),
    "肥大细胞 HPA IL18 331.5 nCPM 为人肺最高、CASP1 11.1（约 30 倍），本数据 CASP1 " + mast.casp1_pos + "/" + mast.n_cells + " = " + pct(mast.casp1) + "；胰岛 δ 细胞有 FACS 分选证据（13.16 vs 0.016 RPKM），两套独立数据同向。"],
  ["B 潜在（相对规则漏掉）", B2, "人肺 AT1 肺泡上皮（单细胞 + 单核）",
    "IL18 分层内排名 3/39（单细胞）与 1/23（单核），HPA IL18 214.1 vs CASP1 7.9 nCPM（27 倍）。未进候选只因该分层另有三分之一细胞类型 CASP1 更低；其自身 CASP1 检出 " + at1.map(r => r.susp + " " + pct(r.casp1)).join(" / ") + "，并不算低，故列为需复核的潜在项。"],
  ["C CASP1 确实存在", C, by("C_casp1_present").map(r => r.sp + r.organ + " " + r.name).join(" · "),
    "肺泡巨噬 IL18 检出 " + pct(amac.il18) + " 为全部候选最高，但 ImmGen 分选 Casp1=511（比值 4.7），本数据 CASP1 " + amac.casp1_pos + "/" + amac.n_cells + " = " + pct(amac.casp1) + "，14/14 供体均有阳性。排除理由是 CASP1 被检出。"],
  ["D 方向相反，排除", D, by("D_refuted").map(r => r.sp + r.organ + " " + r.name).join(" · "),
    "人肝嗜碱粒细胞参考层 CASP1 60–234 nTPM vs IL18 4.7–14.8（比值 0.085）；人肺树突状细胞 2.0 倍。"],
];
let y = 1.6;
rows.forEach(([tag, col, who, why]) => {
  s.addText(tag, { isTextBox: true, x: .7, y, w: 2.75, h: 1.0, fontFace: F, fontSize: 13.5, color: col, bold: true, valign: "top" });
  s.addText([{ text: who + "\n", options: { bold: true } }, { text: why, options: { color: MU } }],
    { isTextBox: true, x: 3.6, y, w: 9.05, h: 1.1, fontFace: F, fontSize: 12, color: INK, valign: "top", lineSpacingMultiple: 1.12 });
  y += 1.1;
});

// ---- 3. figure -----------------------------------------------------------------------------------
s = p.addSlide(); head(s, "候选总览：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_v2.png`, 2890, 1394, 1.22, 6.92);
cap(s, "每行一个细胞类型：实心点 = IL18，空心点 = CASP1，两点之间的距离就是差距；人与小鼠分开，不作跨物种比较。");

// ---- 4. the 0% question --------------------------------------------------------------------------
s = p.addSlide(); head(s, "检出率的口径：为什么上一版这一页全是 0%");
s.addText([{ text: "上一版标注的是「供体中位数检出率」。", options: { bold: true } },
{ text: "每位供体只有数十到上百个该类型细胞，当真实检出率低于约 1/细胞数 时，半数以上供体恰好 0 个阳性细胞，中位数即为 0——这是统计量把低值截断，不是真的没有检出。本版改用合并检出率并标注阳性细胞数。", options: {} }],
  { isTextBox: true, x: .8, y: 1.5, w: 11.7, h: 1.0, fontFace: F, fontSize: 15, color: INK, lineSpacingMultiple: 1.25 });
const ex = [["人·肺 肥大细胞", mast], ["小鼠·脑 小胶质/CNS巨噬", one("小鼠", "小胶质")],
["小鼠·胰 胰导管细胞", one("小鼠", "胰导管")], ["小鼠·脑 少突胶质", one("小鼠", "少突胶质")]];
s.addTable(
  [[{ text: "细胞类型" }, { text: "CASP1 阳性细胞 / 分析细胞" }, { text: "合并检出率" }, { text: "供体中位数" }, { text: "有阳性的供体" }]
    .map(c => ({ ...c, options: { bold: true, fill: "EEF1F5" } }))].concat(
      ex.map(([n, r]) => [n, `${r.casp1_pos} / ${r.n_cells}`, pct(r.casp1), r.med_casp1.toFixed(2) + "%", `${r.casp1_donors} / ${r.n_donors}`])),
  { x: .8, y: 2.85, w: 11.7, colW: [3.3, 3.0, 1.9, 1.7, 1.8], fontFace: F, fontSize: 13, color: INK,
    border: { pt: .5, color: "D8DBE0" }, rowH: .42, valign: "middle" });
s.addText("同样显示「0.0%」的两行，真实检出率相差一百倍以上：肥大细胞 2.60%（138 个阳性细胞、21/44 供体）是真有信号；少突胶质 0.02%（8987 个细胞里 2 个）才是在充分样本下的真·近零。",
  { isTextBox: true, x: .8, y: 5.3, w: 11.7, h: .9, fontFace: F, fontSize: 14, color: INK, lineSpacingMultiple: 1.2 });
cap(s, "分母 = 该细胞类型在该分层内、每位供体 ≥20 个细胞后进入分析的细胞总数；分子 = 其中原始计数 ≥1 的细胞数。不是表达量占比、不是相对最高值的百分比、不是读段占比。");

// ---- 5. scatter ----------------------------------------------------------------------------------
s = p.addSlide(); head(s, "小鼠脑（单细胞）：胶质 IL18 可检、CASP1 处于零底线");
fig(s, `${FIG}/mouse/brain_cell.png`, 1856, 362, 2.3, 5.3);
cap(s, "点大小 = 检出率，颜色 = log2(CPM+1)，三基因各自色标。astrocyte / oligodendrocyte 归入 glia 谱系。");

// ---- 6. species ----------------------------------------------------------------------------------
s = p.addSlide(); head(s, "必须注明的物种差异");
s.addText("A 级结论为小鼠特异，不可外推到人。", { isTextBox: true, x: .8, y: 1.85, w: 11.7, h: .5, fontFace: F, fontSize: 17, color: D, bold: true });
[["小鼠脑（Zhang 分选）", "星形胶质 Il18 115.1 / Casp1 1.13 FPKM（约 102 倍）；少突胶质 24.7 / 0.12。"],
["人脑（HPA 单核）", "小胶质 IL18 338.5 nCPM 压倒性最高，星形胶质仅 3.7、少突胶质 2.6。"],
["含义", "人脑的「IL18 高」属于小胶质，而人小胶质 CASP1 为 41.2 nCPM（比值 8.2），不满足「加工酶缺失」前提。人脑中目前没有对应的 A 级候选。"]
].forEach(([k, v], i) => {
  s.addText(k, { isTextBox: true, x: .8, y: 2.7 + i * 1.05, w: 2.9, h: .9, fontFace: F, fontSize: 13.5, color: INK, bold: true, valign: "top" });
  s.addText(v, { isTextBox: true, x: 3.9, y: 2.7 + i * 1.05, w: 8.7, h: .9, fontFace: F, fontSize: 13.5, color: INK, valign: "top", lineSpacingMultiple: 1.15 });
});

// ---- 7. limits -----------------------------------------------------------------------------------
s = p.addSlide(); head(s, "限制");
["环境 RNA 未校正：无 GPU、Census 源矩阵无空液滴；全部候选为待验证，对低表达的神经/胶质候选影响最大。",
  "单核与单细胞不可混比；跨器官只比候选名单，跨物种不比表达量。",
  "三分位、CASP3 最低 20%、以及本次的 20 倍 / 3 倍参考比值分界均为工程取值，无生物学依据；完整数值已输出，可自行重切。",
  "IL18 转录本 ≠ pro-IL-18 蛋白储量；CASP3 mRNA ≠ caspase-3 酶活性；静息基线不回答刺激后 CASP1 能否诱导。"
].forEach((t, i) => s.addText(t, { isTextBox: true, x: .8, y: 1.85 + i * .95, w: 11.8, h: .8, fontFace: F, fontSize: 14.5, color: INK, lineSpacingMultiple: 1.15, valign: "top" }));

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
