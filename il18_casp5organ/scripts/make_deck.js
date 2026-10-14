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

// ---- 14-15. conclusions and limits ----------------------------------------------------------------
const CJK = s_ => [...s_].reduce((n, ch) => n + (ch.codePointAt(0) > 0x2E80 ? 1 : 0.55), 0);
// 0.86: PowerPoint adds space at every CJK/Latin boundary, so a line holds fewer glyphs than width/size
const wraps = (t_, w, fs) => Math.max(1, Math.ceil(CJK(t_) / (w * 72 / fs * 0.86)));
const LIST = (sl, y, items, fs, w, x) => {
  for (const it of items) {
    const h = wraps(it, w, fs) * (fs / 50) + .10;
    sl.addText(it, { isTextBox: true, x, y, w, h, fontFace: F, fontSize: fs, color: INK,
                     valign: "top", lineSpacingMultiple: 1.2, bullet: { characterCode: "2022", indent: 15 } });
    y += h + .13;
  }
  return y;
};

s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .42, w: 6, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
LIST(s, 1.40, [
  "髓系细胞不满足筛选前提：两物种五器官 62 个「细胞类型 × 分层」中无一符合，IL18 检出率最高者同时也是 CASP1 检出率最高者。",
  "小鼠命中胶质：皮层与海马星形胶质 IL18 56.7%（2135/3766，40 供体）、CASP1 2/3766；少突胶质 27.3%、CASP1 2/8987。单核分层中星形胶质 CASP1 为 4.65%，不复现。",
  "人命中皮层神经元亚类：PV 中间神经元 IL18 26.8%（5286/19711，78 供体）、CASP1 24/19711，另一 234 供体队列复现；L5 ET 三个分层 IL18 31.8%–43.3%。",
  "两组命中不跨物种：人皮层星形胶质 IL18 4.8%，小鼠皮层 PV 中间神经元 Il18 2.2%。",
  "脑外 7 个细胞类型方向符合，证据均为单一来源或存在分歧；最高者人肺肥大细胞 IL18 34.9%（1851/5299，44 供体）、CASP1 2.60%。",
  "人肝与人脊髓白质各 11 类，无符合者。",
  "短 IL-18 转录本在人正常组织不被显著使用：GTEx v8 6792 个样本中占 0.09%。",
], 14.5, 11.5, .95);
s.addText("检出率为合并口径（原始计数 ≥1 的细胞数 ÷ 分析细胞数）。29 个分层，静息基线。完整判定见 results/answer_table.tsv。",
  { isTextBox: true, x: .95, y: H - .52, w: 11.7, h: .34, fontFace: F, fontSize: 10.5, color: FA });

s = p.addSlide();
s.addText("限制与后续", { isTextBox: true, x: .7, y: .42, w: 6, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
LIST(s, 1.40, [
  "静息基线，未做环境 RNA 校正。CASP1 为可诱导基因，静息态近零不预测刺激后；CASP3 的 mRNA 不反映 caspase-3 酶活性，仅用于排除最低档。",
  "物种之间、数据集之间、单核与单细胞之间仅比较候选名单，不比较数值。",
  "参考层覆盖不均：人肺 AT1 与人胰腺腺泡各 1 个参考来源，4 个细胞类型无来源，小鼠星形胶质 4 个。证据分级在相当程度上反映覆盖差异。",
  "零值标尺：静息小鼠小胶质 Casp1 约 50–61 CPM、17% 细胞。阳性对照 L5a IT 锥体神经元 8.34 CPM（25 个文库中 24 个非零）。",
  "小鼠星形胶质结论限皮层与海马：小脑 Bergmann 胶质 22.4 CPM、脊髓星形胶质 12.5–16.7 CPM 均测到 Casp1。",
  "三分位与「CASP3 最低 20%」为工程取值：改用四分位后候选数由 17 降至 12。",
  "待补测量：人肺肥大细胞、人肺 AT1、人胰腺腺泡各需一次原代分选后的 CASP1 定量；小鼠胰岛 δ 细胞需复核比值 5.64 的分歧来源。",
], 14.5, 11.5, .95);

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
