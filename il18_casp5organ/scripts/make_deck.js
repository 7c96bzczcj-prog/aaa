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

// ---- 14-15. conclusions ---------------------------------------------------------------------------
// The previous version listed every cell type under six near-identical bucket names and carried no
// numbers: a classification table, not a conclusion — the reader learned which bucket each row fell
// into but not what the screen found. Replaced with declarative sentences, each one checkable on its
// own. Two slides, because the boundary conditions (§0.2) appeared nowhere else in the deck.
//
// Lead with the myeloid rule rather than the hit list: it is a co-expression pattern, not a near-zero
// claim, so it survives the two limitations that bite every other line here — uncorrected ambient RNA
// and the detection floor.
const CJK = s_ => [...s_].reduce((n, ch) => n + (ch.codePointAt(0) > 0x2E80 ? 1 : 0.55), 0);
const LW = 12.5, GUT = 3.35, BODY_W = 12.58 - GUT;   // gutter label | body column
// a line's height from how many times it wraps, so blocks never overlap however long a sentence runs
const linesOf = t => Math.max(1, Math.ceil(CJK(t) / (BODY_W * 72 / LW - 0.5)));
const SEC = (sl, y, head, note, lines, col) => {
  sl.addText(head, { isTextBox: true, x: .75, y, w: GUT - .25, h: .3, fontFace: F, fontSize: 13, color: col, bold: true });
  if (note) sl.addText(note, { isTextBox: true, x: .75, y: y + .27, w: GUT - .25, h: .62, fontFace: F,
                               fontSize: 9.5, color: MU, valign: "top", lineSpacingMultiple: 1.12 });
  let yy2 = y - .035;
  for (const t of lines) {
    const h = linesOf(t) * .238 + .07;
    sl.addText(t, { isTextBox: true, x: GUT, y: yy2, w: BODY_W, h, fontFace: F, fontSize: LW,
                    color: INK, valign: "top", lineSpacingMultiple: 1.16 });
    yy2 += h;
  }
  return Math.max(yy2 + .17, y + (note ? .95 : .5));
};

s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .38, w: 6, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
s.addText("每行都是可自行核对的观察，不是判定标签。括号内为「阳性细胞 / 分析细胞」，均为合并全部供体的检出率。",
  { isTextBox: true, x: .7, y: .95, w: 12, h: .3, fontFace: F, fontSize: 11, color: MU });
let yy = 1.38;
yy = SEC(s, yy, "一条通则", "不依赖环境 RNA 校正与检出下限", [
  "髓系（巨噬 / 单核 / 树突状 / 小胶质 / 粒细胞）在两物种五器官共 62 个「细胞类型 × 分层」，方向符合的是 0 个：IL18 检出最高的细胞，同时也是 CASP1 检出最高的。",
  "这是共表达模式而不是「近零」断言，因此是本页最不吃边界条件的一条。",
], D);
yy = SEC(s, yy, "小鼠：胶质在顶", "Yao 皮层+海马，真单细胞", [
  "星形胶质 IL18 2135/3766 = 56.7%（40 供体，分层内 1/24），CASP1 仅 2/3766，4 个参考来源无分歧；少突胶质 27.3%（2452/8987，55 供体），CASP1 2/8987。",
  "同为真单细胞的人皮层（78 供体）里这两类只有 4.8%（782/16319）与 1.4%（838/57855）——小鼠这一组在人不成立。",
], A);
yy = SEC(s, yy, "人：神经元在顶", "Allen HVS 单细胞 + 三套单核复现", [
  "PV 中间神经元 IL18 5286/19711 = 26.8%（78 供体），CASP1 24/19711；另一套 234 供体数据复现 24.9%。L5 ET 在三套数据 31.8%–43.3%，CASP1 为 0/201、3/1358、5/2443。",
  "小鼠同名神经元 Casp1 同样近零，但 Il18 只有 2.2%（684/30434）与 5.5%（1333/24181）——人这一组在小鼠不成立。两条物种差异方向相反。",
], B);
yy = SEC(s, yy, "脑外 7 个", "证据都只到单一来源或有分歧", [
  "人肺肥大细胞 34.9%（1851/5299，44 供体）、CASP1 2.60%：4 个参考来源里只有 HPA 组织图谱支持近零（29.9 倍），两个 ex vivo 分选来源都测到 CASP1（3.82 与 4.48 倍）。",
  "人肺 AT1 与人胰腺腺泡各只有 1 个参考来源，且都写明「未进候选（CASP1 未进最低三分位）」——靠参考层补回，不是通过了筛选。小鼠胰岛 δ 细胞 6 个来源中 4 支持、1 分歧；人胰岛 δ 细胞 IL18 仅 8/4677，不复现。",
  "另 4 类（人脑血管软脑膜、人肝巨核-红系祖、人肺泡外膜成纤维、小鼠 Cajal-Retzius）参考层查不到（n_ref_sources = 0）。",
], "4C87C4");
s.addText("32 行判定按「同一细胞类型跨分层」折叠为 22 个类型：15 个方向符合（人 10 / 小鼠 5，两边零重合），7 个不符合。完整 32 行见 results/answer_table.tsv。",
  { isTextBox: true, x: .75, y: H - .52, w: 11.9, h: .34, fontFace: F, fontSize: 10.5, color: FA });

s = p.addSlide();
s.addText("这份名单怎么用，以及它不是什么", { isTextBox: true, x: .7, y: .38, w: 9, h: .55, fontFace: F, fontSize: 26, color: INK, bold: true });
let zz = 1.18;
zz = SEC(s, zz, "弱在哪里", "是没人分选测过，不是检出率低", [
  "8 行「符合 × 多源一致」全部在脑，不等于脑的信号更强：人肺 AT1 与人胰腺腺泡各只有 1 个参考来源，另有 4 类一个都没有，而小鼠星形胶质有 4 个。",
  "要把一格从单一来源推到多源一致，缺的是具体一次测量：人肺肥大细胞 / AT1 / 人胰腺腺泡各缺一次新鲜原代分选后的 CASP1 定量；小鼠胰岛 δ 缺的是复核那一个分歧来源（比值 5.64）。",
], B);
zz = SEC(s, zz, "「近零」相对于什么", "同一测法下的标尺与阳性对照", [
  "阳性对照：同一分选、同一供体、同一建库下，L5a IT 锥体神经元 CASP1 = 8.34 CPM（25 个文库 24 个非零）——不是测不出来，是亚类特异的。",
  "标尺：静息小鼠小胶质 Casp1 约 50–61 CPM / 17% 细胞。小鼠 CNS 里「髓系高」与「零底线」相差约 50 倍，不是 1000 倍。",
], A);
zz = SEC(s, zz, "四条不能省的边界", "规格 §0.2", [
  "IL18 转录本即 pro-IL-18 的转录信号，但 pro-IL-18 是储备型蛋白，转录水平与细胞内前体储量不成线性关系。",
  "CASP1 是可诱导基因。本分析是静息基线筛选，静息近零不预测刺激后能否诱导——这一条只能靠带刺激的数据或实验解决。",
  "CASP3 的 mRNA 不反映 caspase-3 酶活性（酶原活化与 XIAP 抑制都在转录后）。本分析仅用它排除最低档，不是三重阳性门。",
  "未做环境 RNA 校正；人鼠之间、数据集之间、单核与单细胞之间只比名单不比数值。小鼠星形胶质的结论只覆盖皮层与海马：小脑 Bergmann 胶质 22.4 CPM、脊髓星形胶质 12.5–16.7 CPM 都测到 Casp1。",
], D);
zz = SEC(s, zz, "名单有多硬", "阈值与用途", [
  "三分位与「CASP3 最低 20%」是工程取值、无生物学依据：换成四分位，候选数从 17 掉到 12。",
  "短 IL-18 转录本在人正常组织不被显著使用：GTEx v8 6792 个样本 707/791471 跨接读段 = 0.09%（SRA v3h 的 0.20% 来自全人类 RNA-seq 汇编，含癌系与剪接因子敲低样本，不能称作正常组织）。",
  "本名单用于给蛋白与功能层面的验证排优先级。进名单不代表 pro-IL-18 真被加工，不进也不代表不会。",
], MU);

p.writeFile({ fileName: ROOT + "/IL18_CASP_五器官筛选.pptx" }).then(f => console.log("wrote", f));
