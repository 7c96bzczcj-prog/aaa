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
// captions run to three lines on the organ pages, so the box is sized for three and the figure stops
// above it; a one-line caption just leaves white space rather than colliding with the image
const cap = (s, t) => s.addText(t, { isTextBox: true, x: .7, y: H - 1.0, w: 11.95, h: .86, fontFace: F,
                                     fontSize: 10.5, color: MU, valign: "top", lineSpacingMultiple: 1.2 });
const fig = (s, path, iw, ih, top, bot) => { const d = fit(iw, ih, W - 1.4, bot - top); s.addImage({ path, x: (W - d.w) / 2, y: top + (bot - top - d.h) / 2, w: d.w, h: d.h }); };
const pct = v => (v >= 1 ? v.toFixed(1) : v.toFixed(2)) + "%";
const ast = one("小鼠", "星形胶质"), oli = one("小鼠", "少突胶质"), amac = one("小鼠", "肺泡巨噬");

// Slide order follows the user's own edit of the deck (6 slides, conclusions last), with their two
// comments applied: the short-IL-18 page is cut down to data sources only, and the "four questions"
// page is removed. A human three-gene panel is added immediately before the mouse one.

// ---- 1-2. the libraries consulted ----------------------------------------------------------------
// Every citation below was resolved against NCBI E-utilities in-session (PMID printed by
// scratchpad/cite*.py); none is written from memory.
const REFS = [
  ["单细胞 / 单核", [
    ["CZ CELLxGENE Discover（Census LTS 2025-11-08）", "Nucleic Acids Res. 2025 Jan 6;53(D1):D886-D900."],
    ["Allen Human Variation Study（人皮层单细胞）", "Science. 2023 Oct 13;382(6667):eadf2359."],
    ["Allen SEA-AD（MTG / DLPFC 单核）", "Nat Neurosci. 2024 Dec;27(12):2366-2383."],
    ["PsychAD（人前额叶单核）", "Nature. 2026 Sep;657(8133):1003-1015."],
    ["Siletti 人脑细胞图谱（MTG、脊髓）", "Science. 2023 Oct 13;382(6667):eadd7046."],
    ["Yao 小鼠全脑图谱（皮层+海马）", "Nature. 2023 Dec;624(7991):317-332."],
    ["Human Lung Cell Atlas", "Nat Med. 2023 Jun;29(6):1563-1577."],
    ["LungMAP CellRef（人肺参考）", "Nat Commun. 2023 Jul 29;14(1):4566."],
    ["Guilliams 人肝图谱", "Cell. 2022 Jan 20;185(2):379-396.e38."],
    ["Tabula Sapiens（人全器官）", "Science. 2022 May 13;376(6594):eabl4896."],
    ["Tabula Muris Senis（小鼠全器官）", "Nature. 2020 Jul;583(7817):590-595."],
    ["人胰岛多组学", "Nat Metab. 2022 Feb;4(2):284-299."],
    ["小鼠胰岛跨物种图谱", "Mol Metab. 2022 Dec;66:101595."],
    ["人颈段脊髓白质单核", "Acta Neuropathol Commun. 2023 May 22;11(1):84."],
    ["小鼠腰段脊髓单核（GSE103892）", "Cell Rep. 2018 Feb 20;22(8):2216-2225."],
  ]],
];
const REFS2 = [
  ["分选 / 纯化 bulk", [
    ["ImmGen（小鼠免疫细胞，GSE109125）", "Cell. 2019 Feb 7;176(4):897-912.e20."],
    ["brainrnaseq 小鼠皮层分选", "J Neurosci. 2014 Sep 3;34(36):11929-47."],
    ["brainrnaseq 人星形胶质分选", "Neuron. 2016 Jan 6;89(1):37-53."],
    ["Pressl sFANS 分选核（人皮层 L5a）", "Neuron. 2024 Mar 20;112(6):924-941.e10."],
    ["Allen 人皮层 SMART-seq", "Nature. 2019 Sep;573(7772):61-68."],
  ]],
  ["组织 bulk / 蛋白", [
    ["GTEx v8", "Science. 2020 Sep 11;369(6509):1318-1330."],
    ["Human Protein Atlas（RNA + IHC 蛋白）", "Science. 2015 Jan 23;347(6220):1260419."],
    ["Monaco 人免疫细胞绝对定量", "Cell Rep. 2019 Feb 5;26(6):1627-1640.e7."],
    ["Haemopedia（造血谱系）", "Nucleic Acids Res. 2019 Jan 8;47(D1):D780-D785."],
  ]],
  ["剪接连接", [
    ["Snaptron", "Bioinformatics. 2018 Jan 1;34(1):114-116."],
    ["recount3", "Genome Biol. 2021 Nov 29;22(1):323."],
  ]],
  ["基因模型", [
    ["Ensembl", "Nucleic Acids Res. 2025 Jan 6;53(D1):D948-D957."],
  ]],
];
const refPage = (groups, title, foot) => {
  const sl = p.addSlide(); head(sl, title);
  let y = 1.35;
  for (const [cat, items] of groups) {
    sl.addText(cat, { isTextBox: true, x: .75, y, w: 2.3, h: .3, fontFace: F, fontSize: 12.5, color: MU, bold: true });
    for (const [name, cit] of items) {
      sl.addText(name, { isTextBox: true, x: 3.05, y, w: 5.35, h: .28, fontFace: F, fontSize: 11.5, color: INK, valign: "top" });
      sl.addText(cit, { isTextBox: true, x: 8.45, y, w: 4.2, h: .28, fontFace: F, fontSize: 11.5, color: "3C3C3C", valign: "top" });
      y += .292;
    }
    y += .17;
  }
  if (foot) sl.addText(foot, { isTextBox: true, x: .75, y: H - .52, w: 11.9, h: .3, fontFace: F, fontSize: 10, color: FA });
};
let s;
refPage(REFS, "查阅的数据库（一）：单细胞 / 单核", "本筛选从 CELLxGENE Census 取 29 个分层。全部引文经 NCBI E-utilities 逐条核对。");
refPage(REFS2, "查阅的数据库（二）：参考层", "统一参考表 1714 行，每条带来源、单位、状态与匹配质量标注。");

// ---- 2-3. per-species overview figures ----------------------------------------------------------
s = p.addSlide(); head(s, "人：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_human.png`, 2890, 1550, 1.15, H - 1.05);
cap(s, "实心点 = IL18，空心点 = CASP1。同一细胞类型在多个数据集重复时取供体最多的一套。");

s = p.addSlide(); head(s, "小鼠：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_mouse.png`, 2890, 908, 1.4, 6.3);
cap(s, "人与小鼠分开绘制，不作跨物种定量比较。");

// ---- 4-13. three-gene panels: five organs × two species -------------------------------------------
// Earlier drafts carried brain only. That was an artefact of which page got built first, not a
// decision — the screen covers five organs, so every organ gets the same page, including the ones
// where nothing qualifies. Each page lists every cell type of one stratum, whatever its verdict.
const PANEL = [
  ["人肺（HLCA 单细胞）", "human/lung_dots.png",
   "肥大细胞是人肺唯一两项都成立的一档；IL18 最高的三类 CASP1 同样高。"],
  ["人肝（Guilliams 单细胞）", "human/liver_dots.png",
   "11 类无一符合。肝的 IL18 集中在髓系，而髓系的 CASP1 最高。"],
  ["人胰（48 供体胰岛单细胞）", "human/pancreas_dots.png",
   "胰腺腺泡细胞 IL18 30.3%、CASP1 3.45%，人胰唯一方向相符者。"],
  ["人脑（Allen HVS 单细胞）", "human/brain_cell_dots.png",
   "PV 与 L5 ET 的 CASP1 近零，CASP3 为 27.9% 与 34.3%，不在最低档。"],
  ["人脊髓（颈段白质单核）", "human/spinal_dots.png",
   "小胶质 IL18 52.3% 但 CASP1 11.8%，其余 10 类 IL18 均 ≤3.2%。"],
  ["小鼠肺（Tabula Muris Senis 10x）", "mouse/lung_dots.png",
   "肺泡巨噬 IL18 71.6% 全器官最高，CASP1 13.3%，参考层确认可测。"],
  ["小鼠肝（Tabula Muris Senis 10x）", "mouse/liver_dots.png",
   "肝细胞 IL18 11.2%、CASP1 0.34% 方向相符，但 IL18 未进最高三分位。"],
  ["小鼠胰（小鼠胰岛图谱单细胞）", "mouse/pancreas_dots.png",
   "胰岛 δ 细胞与胰导管细胞方向相符；人胰岛 δ 的 IL18 仅 0.17%。"],
  ["小鼠脑（Yao 皮层+海马单细胞）", "mouse/brain_cell_dots.png",
   "星形胶质 IL18 56.7%、CASP1 2/3766；小胶质 CASP1 1.83%，全脑最高。"],
  ["小鼠脊髓（GSE103892 单核）", "mouse/spinal_dots.png",
   "星形胶质 CASP1 0/404，但参考层 RiboTag 分选脊髓星形胶质测到 Casp1。"],
];
// placed 1:1 — answer_pack.py draws these at exactly 11.93 × 5.30 in so 8 pt in the figure is 8 pt here
PANEL.forEach(([title, f, c]) => {
  s = p.addSlide(); head(s, title + "：IL18 / CASP1 / CASP3 三基因");
  s.addImage({ path: `${FIG}/${f}`, x: (W - 11.93) / 2, y: 1.15, w: 11.93, h: 5.30 });
  cap(s, c);
});

// ---- 14. conclusion: the cell types the screen returns, ordered by reference support -------------
const HIT = R.filter(r => r.direction === "符合")
  .sort((a, b) => b.n_donors - a.n_donors)
  .filter((r, i, arr) => arr.findIndex(x => x.sp === r.sp && x.organ === r.organ && x.name === r.name) === i)
  .sort((a, b) => b.n_sup - a.n_sup || b.il18 - a.il18);
const pc = v => (v >= 10 ? v.toFixed(1) : v.toFixed(2)) + "%";

s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .36, w: 6, h: .5, fontFace: F, fontSize: 24, color: INK, bold: true });
s.addText("符合「IL18 高 + CASP1 低」的细胞类型，按参考层支持由强到弱。支持 / 分歧 = 参考层中显示 CASP1 近零 / 可测的独立来源数",
  { isTextBox: true, x: .7, y: .88, w: 11.9, h: .3, fontFace: F, fontSize: 12, color: MU });
const HEAD = ["物种 · 器官", "细胞类型", "IL18 检出", "CASP1 检出", "CASP3 检出", "供体", "支持", "分歧"];
s.addTable([HEAD.map(c => ({ text: c, options: { bold: true, fill: "EEF1F5", color: "3C3C3C" } }))].concat(
  HIT.map(r => [`${r.sp} · ${r.organ}`, r.name + (r.n_strata > 1 ? ` ×${r.n_strata}` : ""),
                `${pc(r.il18)}　${r.il18_pos}/${r.n_cells}`,
                `${pc(r.casp1)}　${r.casp1_pos}/${r.n_cells}`,
                `${pc(r.casp3)}　${r.casp3_pos}/${r.n_cells}${r.casp3_floor ? "  ▽" : ""}`, `${r.n_donors}`,
                r.n_sup ? `${r.n_sup}` : "—", r.n_dis ? `${r.n_dis}` : "—"])),
  { x: .7, y: 1.26, w: 11.93, colW: [1.38, 2.72, 2.14, 2.14, 2.14, .55, .46, .40], fontFace: F, fontSize: 10,
    color: INK, border: { pt: .5, color: "E2E5E9" }, rowH: .225, valign: "middle", autoPage: false });
s.addText([
  { text: "明确阴性：", options: { bold: true } },
  { text: "髓系 62 个「细胞类型 × 分层」0 个符合；人肝 11 类、人脊髓白质 11 类 0 个符合；短 IL-18 转录本占 GTEx v8 6792 个样本的 0.09%。" }],
  { isTextBox: true, x: .7, y: H - 1.24, w: 11.93, h: .4, fontFace: F, fontSize: 12, color: INK, valign: "top", lineSpacingMultiple: 1.2 });
s.addText([
  { text: "物种不互通：", options: { bold: true } },
  { text: "人皮层星形胶质 IL18 4.8%（782/16319）；小鼠皮层 PV 中间神经元 Casp1 同样近零，但 Il18 仅 2.2%（684/30434）。" }],
  { isTextBox: true, x: .7, y: H - .86, w: 11.93, h: .4, fontFace: F, fontSize: 12, color: INK, valign: "top", lineSpacingMultiple: 1.2 });
s.addText("检出率为合并口径（原始计数 ≥1 的细胞数 ÷ 分析细胞数）。同一细胞类型见于多个分层时取供体最多者，×n 为分层数。CASP3 不参与候选定义，仅用于排除分层内最低档，▽ 为处于该分层 CASP3 最低档者。静息基线，未做环境 RNA 校正。",
  { isTextBox: true, x: .7, y: H - .44, w: 11.93, h: .3, fontFace: F, fontSize: 10, color: FA });

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
