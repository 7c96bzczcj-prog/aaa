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
const ast = one("小鼠", "星形胶质"), oli = one("小鼠", "少突胶质"), amac = one("小鼠", "肺泡巨噬");

// Slide order follows the user's own edit of the deck (6 slides, conclusions last), with their two
// comments applied: the short-IL-18 page is cut down to data sources only, and the "four questions"
// page is removed. A human three-gene panel is added immediately before the mouse one.

// ---- 1. short IL-18: data sources only (per comment on this page) --------------------------------
let s = p.addSlide(); head(s, "关于「短 IL-18」：转录本层面是明确的否定");
s.addText("人任何正常组织、任何正常细胞类型都不显著使用短 IL18 转录本。",
  { isTextBox: true, x: .8, y: 1.35, w: 11.7, h: .45, fontFace: F, fontSize: 17, color: D, bold: true });
s.addTable([[{ text: "数据源" }, { text: "内容" }, { text: "短异构体占比" }].map(c => ({ ...c, options: { bold: true, fill: "EEF1F5" } }))].concat([
  ["Ensembl REST", "IL18 转录本结构与异构体蛋白长度（人 193 aa vs 148 / 144 aa）", "—"],
  ["Snaptron（GTEx v8）", "逐样本剪接连接计数，54 个组织、6792 个样本", "0.09%"],
  ["Snaptron（recount / SRA v3h）", "逐样本剪接连接计数，50790 个样本", "0.20%"]]),
  { x: .8, y: 2.05, w: 11.7, colW: [3.0, 6.4, 2.3], fontFace: F, fontSize: 13, color: INK,
    border: { pt: .5, color: "D8DBE0" }, rowH: .44, valign: "middle" });
s.addText("本分析的全部参考数据库", { isTextBox: true, x: .8, y: 4.0, w: 11.7, h: .32, fontFace: F, fontSize: 13.5, color: INK, bold: true });
[["单细胞/单核", "CZ CELLxGENE Census LTS 2025-11-08（29 个分层）· Allen Human Variation Study · SEA-AD · PsychAD · Tabula Muris Senis · Tabula Sapiens · HLCA · CellRef"],
 ["分选 / 纯化 bulk", "ImmGen（RNA-seq + 微阵列）· brainrnaseq/Zhang · Pressl 2024 sFANS 分选核 · Allen Human Cortex SMART-seq · RiboTag / TRAP 系列 · LCM 分选 · 多套 GEO FACS 分选"],
 ["组织 bulk / 蛋白", "GTEx v8 · Human Protein Atlas v25（RNA + IHC 蛋白）· Monaco 2019 · Haemopedia"],
 ["剪接连接", "Snaptron（recount3 / SRA / GTEx）"]
].forEach(([k, v], i) => {
  s.addText(k, { isTextBox: true, x: .8, y: 4.4 + i * .62, w: 2.3, h: .55, fontFace: F, fontSize: 12, color: MU, bold: true, valign: "top" });
  s.addText(v, { isTextBox: true, x: 3.2, y: 4.4 + i * .62, w: 9.3, h: .58, fontFace: F, fontSize: 11.5, color: INK, valign: "top", lineSpacingMultiple: 1.1 });
});
s.addText("统一参考表 1714 行，每条带来源、单位、状态与匹配质量标注。", { isTextBox: true, x: .8, y: 6.95, w: 11.7, h: .3, fontFace: F, fontSize: 11, color: FA });

// ---- 2-3. per-species overview figures ----------------------------------------------------------
s = p.addSlide(); head(s, "人：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_human.png`, 2890, 1550, 1.15, 7.0);
cap(s, "每行一个细胞类型：实心点 = IL18，空心点 = CASP1，两点距离即差距；同一细胞类型在多个数据集重复时取供体最多的一套，右侧标 ×n。");

s = p.addSlide(); head(s, "小鼠：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_mouse.png`, 2890, 908, 1.4, 6.3);
cap(s, "人与小鼠分开绘制，不作跨物种定量比较。");

// ---- 4-5. three-gene panels, human first then mouse ---------------------------------------------
s = p.addSlide(); head(s, "人脑（单细胞）：IL18 / CASP1 / CASP3 三基因");
fig(s, `${FIG}/human/brain_cell_dots.png`, 3162, 1564, 1.15, 7.0);
cap(s, "Allen HVS 真单细胞、78 供体；该分层全部 18 个细胞类型一律列出。PV 与 L5 ET 的 CASP1 近零而 CASP3 为 27.9% 与 34.3%，即不满足「CASP3 最低档」这一排除项。");

s = p.addSlide(); head(s, "小鼠脑（单细胞）：IL18 / CASP1 / CASP3 三基因");
fig(s, `${FIG}/mouse/brain_cell_dots.png`, 2312, 1122, 1.3, 6.6);
cap(s, "星形胶质 IL18 56.7% 而 CASP1 仅 2/3766；小胶质 IL18 20.5% 但 CASP1 15/818 = 1.83%，是全脑最高的一档。");

// ---- 6. conclusions: two axes, factual labels -----------------------------------------------------
s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .4, w: 6, h: .6, fontFace: F, fontSize: 26, color: INK, bold: true });
s.addText("「方向」= 是否符合「IL18 高 + CASP1 低」这一前提；「证据」= 支持该方向的独立参考来源有多少、是否一致。两者分开看。",
  { isTextBox: true, x: .7, y: 1.0, w: 12, h: .35, fontFace: F, fontSize: 12, color: MU });
const G = (dir, evi) => uniq(R.filter(r => r.direction === dir && r.evidence === evi)
  .map(r => r.sp + r.organ + " " + r.name.replace(/（[^）]*）$/, "")));
const BLK = (y, dir, evi, col, note) => {
  const items = G(dir, evi);
  if (!items.length) return y;
  const one = R.find(r => r.direction === dir && r.evidence === evi);
  s.addText(`${dir} · ${evi}`, { isTextBox: true, x: .8, y, w: 3.1, h: .3, fontFace: F, fontSize: 13.5, color: col, bold: true });
  s.addText(note, { isTextBox: true, x: .8, y: y + .28, w: 3.1, h: .5, fontFace: F, fontSize: 9.5, color: MU, valign: "top", lineSpacingMultiple: 1.1 });
  s.addText(items.join(" · "), { isTextBox: true, x: 4.1, y: y - .04, w: 8.5, h: .72, fontFace: F, fontSize: 12.5, color: INK, valign: "top", lineSpacingMultiple: 1.15 });
  return y + .82;
};
let yy = 1.42;
yy = BLK(yy, "符合", "多源一致", A, "≥2 个独立来源，无一例外");
yy = BLK(yy, "符合", "多源有分歧", "2E8B57", "多数来源支持，但有来源测到 CASP1");
yy = BLK(yy, "符合", "单一来源", B, "仅 1 个来源支持");
yy = BLK(yy, "符合", "无参考可查", "4C87C4", "无分选参考，仅本数据方向相符");
yy = BLK(yy, "不符合", "多源一致", D, "参考层 CASP1 确实可测或 ≥ IL18");
yy = BLK(yy, "不符合", "无参考可查", "B5564E", "无参考，本数据自身即方向相反");
s.addText("两条必须一起说的边界：小鼠星形胶质的结论仅限皮层与海马（小脑 Bergmann 胶质与脊髓星形胶质 Casp1 明确表达）；人的 PV / L5 ET 神经元结论在小鼠没有对应物（小鼠这两类 Casp1 同样近零，但 Il18 不高）。",
  { isTextBox: true, x: .8, y: 6.25, w: 11.8, h: .6, fontFace: F, fontSize: 11.5, color: INK, lineSpacingMultiple: 1.2 });
s.addText("检出率 = 原始计数 ≥1 的细胞数 ÷ 分析细胞数（合并全部供体）。29 个分层，人与小鼠分开，静息基线；未做环境 RNA 校正。",
  { isTextBox: true, x: .8, y: 6.95, w: 11.8, h: .35, fontFace: F, fontSize: 10.5, color: FA });

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
