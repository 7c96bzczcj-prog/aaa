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
cap(s, "每行一个细胞类型：实心点 = IL18，空心点 = CASP1，两点距离即差距；同一细胞类型在多个数据集重复时取供体最多的一套，右侧标 ×n。");

s = p.addSlide(); head(s, "小鼠：同一条检出率轴上看 IL18 与 CASP1");
fig(s, `${FIG}/answer_mouse.png`, 2890, 908, 1.4, 6.3);
cap(s, "人与小鼠分开绘制，不作跨物种定量比较。");

// ---- 4-13. three-gene panels: five organs × two species -------------------------------------------
// Earlier drafts carried brain only. That was an artefact of which page got built first, not a
// decision — the screen covers five organs, so every organ gets the same page, including the ones
// where nothing qualifies. Each page lists every cell type of one stratum, whatever its verdict.
const PANEL = [
  ["人肺（HLCA 单细胞）", "human/lung_dots.png",
   "肥大细胞 IL18 34.9% 而 CASP1 2.60%（138/5299），是人肺唯一两项都成立的一档。IL18 最高的三类反而不符合：CD1c+ 树突状 52.8% 对 CASP1 32.3%、间皮 51.8% 对 18.3%、肺泡巨噬 50.7% 对 57.9%。该分层 39 个细胞类型，本页按 IL18 取前 20。"],
  ["人肝（Guilliams 单细胞）", "human/liver_dots.png",
   "11 个细胞类型无一符合。IL18 最高的巨噬 53.1% 对 CASP1 38.5%，浆细胞样树突状 37.7% 对 42.8%，单核 20.9% 对 52.3%——肝的 IL18 信号集中在髓系，而髓系正是 CASP1 最高的一类。"],
  ["人胰（48 供体胰岛单细胞）", "human/pancreas_dots.png",
   "胰腺腺泡细胞 IL18 30.3%（1301/4288）而 CASP1 3.45%，是人胰唯一方向相符者，且被分层内相对规则漏掉、由参考层补回。同页肥大细胞 IL18 34.9% 但此分层 CASP1 10.4%，与人肺的肥大细胞不是同一回事。"],
  ["人脑（Allen HVS 单细胞）", "human/brain_cell_dots.png",
   "78 供体真单细胞，该分层全部 18 类一律列出。PV 与 L5 ET 的 CASP1 近零而 CASP3 为 27.9% 与 34.3%，即不满足「CASP3 最低档」这一排除项。"],
  ["人脊髓（颈段白质单核）", "human/spinal_dots.png",
   "11 类全部不符合：小胶质 IL18 52.3% 但 CASP1 11.8%（215/1829），其余 10 类 IL18 均 ≤3.2%。人脊髓白质在转录本层面是明确的否定，不是「数据不足」。"],
  ["小鼠肺（Tabula Muris Senis 10x）", "mouse/lung_dots.png",
   "肺泡巨噬 IL18 71.6% 为全器官最高，但 CASP1 13.3%（178/1341），参考层亦确认可测，判为不符合；其余 16 类 IL18 均 ≤18.5%。"],
  ["小鼠肝（Tabula Muris Senis 10x）", "mouse/liver_dots.png",
   "Kupffer 细胞 IL18 74.2% 对 CASP1 51.9%。肝细胞 IL18 11.2% 而 CASP1 0.34% 方向相符，但 IL18 未进该分层最高三分位，按 §8 不入候选——这一条列在此处是为了让落选理由可见。"],
  ["小鼠胰（小鼠胰岛图谱单细胞）", "mouse/pancreas_dots.png",
   "胰岛 D(δ) 细胞与胰导管细胞是本器官的两个方向相符项。胰岛内分泌细胞 IL18 19.5% 而 CASP1 0.09%，CASP3 33.3%。同一位置的人胰岛 δ 细胞 IL18 仅 0.17%（8/4677），该发现在人不复现。"],
  ["小鼠脑（Yao 皮层+海马单细胞）", "mouse/brain_cell_dots.png",
   "星形胶质 IL18 56.7% 而 CASP1 仅 2/3766；小胶质 IL18 20.5% 但 CASP1 15/818 = 1.83%，是全脑最高的一档。"],
  ["小鼠脊髓（GSE103892 单核）", "mouse/spinal_dots.png",
   "本页不产候选。星形胶质 IL18 21.0% 而 CASP1 0/404、少突胶质 1/707，方向看似相符，但该分层矩阵只带基因符号（gene_id_unavailable，§1.2），细胞数又少；而参考层的 RiboTag 分选脊髓星形胶质明确测到 Casp1。故小鼠星形胶质的结论只覆盖皮层与海马，不含脊髓。"],
];
// placed 1:1 — answer_pack.py draws these at exactly 11.93 × 5.30 in so 8 pt in the figure is 8 pt here
PANEL.forEach(([title, f, c]) => {
  s = p.addSlide(); head(s, title + "：IL18 / CASP1 / CASP3 三基因");
  s.addImage({ path: `${FIG}/${f}`, x: (W - 11.93) / 2, y: 1.15, w: 11.93, h: 5.30 });
  cap(s, c);
});

// ---- 14-16. conclusions --------------------------------------------------------------------------
// The conclusion page states findings only: one claim per line, the measurement in parentheses, no
// rhetorical asides. Everything that qualifies a finding — reference-layer coverage, source
// independence, the within-species non-replication, the zero calibration — moves to p.15, and the
// §0.2 boundary conditions to p.16, so the conclusion itself stays readable in one pass.
const CJK = s_ => [...s_].reduce((n, ch) => n + (ch.codePointAt(0) > 0x2E80 ? 1 : 0.55), 0);
// 0.86: PowerPoint adds space at every CJK/Latin boundary, so a line holds fewer glyphs than width/size
const wraps = (t_, w, fs) => Math.max(1, Math.ceil(CJK(t_) / (w * 72 / fs * 0.86)));

// --- 14. findings ---------------------------------------------------------------------------------
s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .40, w: 6, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
s.addText("检出率均为合并口径（原始计数 ≥1 的细胞数 ÷ 分析细胞数，合并全部供体）。29 个分层，静息基线。",
  { isTextBox: true, x: .7, y: .97, w: 12, h: .3, fontFace: F, fontSize: 11, color: MU });
const FIND = [
  ["髓系细胞不满足筛选前提。",
   "巨噬、单核、树突状、小胶质与粒细胞在两物种五器官共 62 个「细胞类型 × 分层」中无一方向符合；IL18 检出率最高的细胞类型同时也是 CASP1 检出率最高的。"],
  ["小鼠的命中集中于胶质。",
   "皮层与海马星形胶质 IL18 56.7%（2135/3766，40 供体）、CASP1 2/3766；少突胶质 27.3%（2452/8987，55 供体）、CASP1 2/8987。"],
  ["人的命中集中于皮层神经元亚类。",
   "PV 中间神经元 IL18 26.8%（5286/19711，78 供体）、CASP1 24/19711，另一队列（234 供体）复现；L5 ET 投射神经元三个分层 IL18 31.8%–43.3%、CASP1 均 ≤5 个阳性细胞。"],
  ["两组命中不跨物种成立。",
   "人皮层星形胶质 IL18 仅 4.8%（782/16319，78 供体）；小鼠皮层 PV 中间神经元 Casp1 同样近零，但 Il18 仅 2.2%（684/30434，52 供体）。"],
  ["脑外共 7 个细胞类型方向符合，证据均为单一来源或存在分歧。",
   "其中检出率最高者为人肺肥大细胞 IL18 34.9%（1851/5299，44 供体）、CASP1 2.60%。"],
  ["人肝与人脊髓无方向符合者。",
   "Guilliams 肝图谱 11 类、颈段脊髓白质 11 类，逐类核对后均不符合。"],
  ["短 IL-18 转录本在人正常组织中不被显著使用。",
   "GTEx v8 6792 个样本中占 0.09%（707/791471 跨接读段）。"],
];
let yy = 1.52;
FIND.forEach(([lead, rest], i) => {
  const h = wraps(lead + rest, 11.15, 14) * .27 + .10;
  s.addText(`${i + 1}`, { isTextBox: true, x: .72, y: yy, w: .42, h: .3, fontFace: F, fontSize: 12, color: FA, align: "right" });
  s.addText([{ text: lead, options: { bold: true, color: INK } }, { text: rest, options: { color: "3C3C3C" } }],
    { isTextBox: true, x: 1.30, y: yy - .045, w: 11.15, h, fontFace: F, fontSize: 14, valign: "top", lineSpacingMultiple: 1.22 });
  yy += h + .12;
});
s.addText("32 行判定按同一细胞类型跨分层折叠为 22 个类型：15 个方向符合（人 10、小鼠 5，两物种零重合），7 个不符合。完整判定见 results/answer_table.tsv。",
  { isTextBox: true, x: .75, y: H - .50, w: 11.9, h: .34, fontFace: F, fontSize: 10.5, color: FA });

// --- 15-16. what qualifies the findings, and the boundary conditions -------------------------------
const LW = 12, GUT = 3.35, BODY_W = 12.58 - GUT;
const SEC = (sl, y, head_, note, lines, col) => {
  sl.addText(head_, { isTextBox: true, x: .75, y, w: GUT - .25, h: .3, fontFace: F, fontSize: 12.5, color: col, bold: true });
  if (note) sl.addText(note, { isTextBox: true, x: .75, y: y + .26, w: GUT - .25, h: .6, fontFace: F,
                               fontSize: 9.5, color: MU, valign: "top", lineSpacingMultiple: 1.12 });
  let y2 = y - .035;
  for (const t_ of lines) {
    const h = wraps(t_, BODY_W, LW) * .232 + .07;
    sl.addText(t_, { isTextBox: true, x: GUT, y: y2, w: BODY_W, h, fontFace: F, fontSize: LW,
                     color: INK, valign: "top", lineSpacingMultiple: 1.16 });
    y2 += h;
  }
  return Math.max(y2 + .18, y + (note ? .92 : .5));
};

s = p.addSlide();
s.addText("判定依据与证据强度", { isTextBox: true, x: .7, y: .40, w: 9, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
let zz = 1.22;
zz = SEC(s, zz, "参考层覆盖不均", "影响证据分级，不影响检出率", [
  "8 行「符合 × 多源一致」全部位于脑。人肺 AT1 与人胰腺腺泡各仅 1 个参考来源，另有 4 个细胞类型无参考来源（n_ref_sources = 0），而小鼠星形胶质为 4 个。证据分级在相当程度上反映参考层覆盖，而非检出率差异。",
  "人肺 AT1 与人胰腺腺泡的 screen_status 为「未进候选（CASP1 未进最低三分位）」，二者系由参考层补回，未通过分层内相对规则。",
], B);
zz = SEC(s, zz, "来源独立性", "", [
  "小鼠星形胶质的 4 个参考来源中，3 个为独立纯化测量（GSE100329 RiboTag、Tabula Muris Senis FACS-SmartSeq2、brainrnaseq/Zhang），第 4 个为 CELLxGENE 跨语料汇总，与本筛选同属一个数据家族。",
  "人肺肥大细胞的 4 个来源中，支持 CASP1 近零者为 HPA 组织图谱（比值 29.9），两个 ex vivo 分选来源均测到 CASP1（比值 3.82 与 4.48）。",
], B);
zz = SEC(s, zz, "分层内不复现", "按 §3 应予报告", [
  "小鼠脑另一分层（单核，10 供体）中星形胶质 IL18 仍居首位（14579/21878 = 66.6%），但 CASP1 为 1018/21878 = 4.65%，不入候选。单核与单细胞不作定量比较，候选集成员这一层级的分歧仍须报告。",
  "人脑 CGE 源中间神经元为全表证据最弱的一条「符合」：IL18 12.1%（618/5096，44 供体），5 个参考来源中 1 个支持、3 个反对，其中 2 个显示 CASP1 与 IL18 处于同一量级。",
], D);
zz = SEC(s, zz, "零值标尺", "判断「近零」的参照", [
  "静息小鼠小胶质 Casp1 约 50–61 CPM、17% 细胞；小鼠 CNS 内「髓系高」与零底线相差约 50 倍。阳性对照：同一分选、同一供体、同一建库下 L5a IT 锥体神经元 CASP1 = 8.34 CPM（25 个文库中 24 个非零）。",
], A);
zz = SEC(s, zz, "后续测量", "将单一来源推进为多源一致所需", [
  "人肺肥大细胞、人肺 AT1、人胰腺腺泡各需一次新鲜原代分选后的 CASP1 定量；小鼠胰岛 δ 细胞需复核比值为 5.64 的分歧来源。",
], MU);

s = p.addSlide();
s.addText("适用边界与阈值敏感性", { isTextBox: true, x: .7, y: .40, w: 11, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
zz = 1.22;
zz = SEC(s, zz, "四条适用边界", "规格 §0.2", [
  "IL18 转录本为 pro-IL-18 的转录信号；pro-IL-18 属储备型蛋白，转录水平与细胞内前体储量不成线性关系。",
  "CASP1 为可诱导基因。本分析为静息基线筛选，静息态近零不预测刺激后能否诱导，该限制只能由带刺激条件的数据或实验解决。",
  "CASP3 的 mRNA 不反映 caspase-3 酶活性（酶原活化与 XIAP 抑制均发生于转录后）。本分析仅用其排除最低档，不构成三重阳性门。",
  "未做环境 RNA 校正；物种之间、数据集之间、单核与单细胞之间仅比较候选名单，不比较数值。小鼠星形胶质的结论仅覆盖皮层与海马：小脑 Bergmann 胶质 22.4 CPM、脊髓星形胶质 12.5–16.7 CPM 均测到 Casp1。",
], D);
zz = SEC(s, zz, "阈值敏感性", "§9.6", [
  "三分位与「CASP3 最低 20%」为工程取值，无生物学依据。改用四分位后候选数由 17 降至 12；CASP3 下限改为 10% 与 30% 时分别为 17 与 16。",
], MU);
zz = SEC(s, zz, "用途", "", [
  "本名单用于为蛋白与功能层面的验证排定器官与细胞类型的优先级。列入名单不表示 pro-IL-18 确被加工，未列入亦不表示不会。",
], MU);

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
