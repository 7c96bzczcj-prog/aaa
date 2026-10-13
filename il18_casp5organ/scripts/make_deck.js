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

// ---- 1. one-line conclusions ---------------------------------------------------------------------
let s = p.addSlide();
s.addText("结论", { isTextBox: true, x: .7, y: .45, w: 12, h: .7, fontFace: F, fontSize: 26, color: INK, bold: true });
const L = (y, color, txt, sub) => {
  s.addText(txt, { isTextBox: true, x: .8, y, w: 11.7, h: .36, fontFace: F, fontSize: 15.5, color, bold: true, valign: "top" });
  if (sub) s.addText(sub, { isTextBox: true, x: .8, y: y + .33, w: 11.7, h: .3, fontFace: F, fontSize: 11, color: MU, valign: "top" });
};
s.addText("确定 —— 两个以上独立参考来源一致", { isTextBox: true, x: .8, y: 1.25, w: 11.7, h: .3, fontFace: F, fontSize: 12.5, color: A, bold: true });
L(1.62, INK, "人脑皮层 PV 中间神经元  IL18 高、CASP1 近零", "HVS 单细胞 78 供体 + PsychAD 单核 234 供体；参考层 sFANS 分选核 566 倍、Allen SMART-seq 358 倍");
L(2.28, INK, "人脑皮层 L5 ET 投射神经元  IL18 高、CASP1 近零", "HVS + SEA-AD MTG + SEA-AD DLPFC 三套；参考层 Betz 细胞 275 倍、von Economo 神经元 384 倍");
L(2.94, INK, "小鼠皮层与海马 星形胶质  Il18 高、Casp1 近零", "4 个来源；但仅限皮层与海马——小脑 Bergmann 胶质与脊髓星形胶质 Casp1 明确表达");
L(3.60, INK, "小鼠脑 少突胶质  Il18 高、Casp1 近零", "3 个来源；RiboTag 中 Casp1 比配对全组织更低，是去除而非带入");
L(4.26, INK, "小鼠胰岛 D(δ) 细胞  Il18 可检、Casp1 近零", "6 个来源；GSE269616 九个文库 Casp1 零条读段。但人胰岛 δ 细胞不成立");
s.addText("可能 —— 单一来源支持，或尚无匹配参考", { isTextBox: true, x: .8, y: 5.0, w: 11.7, h: .3, fontFace: F, fontSize: 12.5, color: B, bold: true });
s.addText("人肺 AT1 肺泡上皮 · 人肺 肥大细胞 · 人胰腺 腺泡细胞 · 小鼠胰 导管细胞 · 人脑 血管软脑膜细胞 · 小鼠脑 Cajal-Retzius 神经元",
  { isTextBox: true, x: .8, y: 5.32, w: 11.7, h: .4, fontFace: F, fontSize: 13, color: INK });
s.addText("相反 —— CASP1 在其中确实被检出，不是测不准", { isTextBox: true, x: .8, y: 5.85, w: 11.7, h: .3, fontFace: F, fontSize: 12.5, color: D, bold: true });
s.addText("小鼠肺 肺泡巨噬 · 小鼠脑 小胶质 · 人肺 气道/纤毛上皮 · 人肝 嗜碱粒细胞 · 人肝 CD56bright NK · 人肺 树突状细胞",
  { isTextBox: true, x: .8, y: 6.17, w: 11.7, h: .4, fontFace: F, fontSize: 13, color: INK });
s.addText("检出率 = 原始计数 ≥1 的细胞数 ÷ 分析细胞数（合并全部供体）。29 个分层，人与小鼠分开，静息基线。",
  { isTextBox: true, x: .8, y: 6.75, w: 11.7, h: .35, fontFace: F, fontSize: 11, color: FA });

// ---- 2. short IL-18 ------------------------------------------------------------------------------
s = p.addSlide(); head(s, "关于「短 IL-18」：转录本层面是明确的否定");
s.addText([{ text: "人任何正常组织、任何正常细胞类型都不显著使用短 IL18 转录本。", options: { bold: true, color: D } }],
  { isTextBox: true, x: .8, y: 1.35, w: 11.7, h: .4, fontFace: F, fontSize: 16, color: INK });
s.addText("长短异构体差一个内含框盒式外显子 E4（chr11:112150072-112150206，135 bp = 45 个密码子）：跳过它 193 aa → 148 aa。所有相关剪接连接共用同一受体位点，因此可直接算同一位点上的竞争比例。",
  { isTextBox: true, x: .8, y: 1.9, w: 11.7, h: .7, fontFace: F, fontSize: 13, color: MU, lineSpacingMultiple: 1.2 });
s.addTable([[{ text: "数据源" }, { text: "短异构体占比" }, { text: "说明" }].map(c => ({ ...c, options: { bold: true, fill: "EEF1F5" } }))].concat([
  ["GTEx 全组织（6792 样本）", "0.09%", "最高的组织也只有 0.65%，无任何组织富集"],
  ["SRA（50790 样本）", "0.20%", "占比最高的几个研究全是剪接机器被破坏的实验或癌细胞系"],
  ["U2AF1 敲低（HFF 原代）", "11.1%", "U2AF1 是 3' 剪接位点识别的核心因子"],
  ["Tra2 敲低（MDA-MB-231）", "13.5%", "Tra2 为剪接因子"],
  ["CLK2 抑制剂 CC-671（Cal-51）", "4.2%", "125 个样本、84 个检出；CLK 磷酸化 SR 蛋白"]]),
  { x: .8, y: 2.75, w: 11.7, colW: [4.0, 2.2, 5.5], fontFace: F, fontSize: 12.5, color: INK,
    border: { pt: .5, color: "D8DBE0" }, rowH: .4, valign: "middle" });
s.addText([{ text: "所以短 IL18 转录本是剪接噪声，不是受调控的细胞类型特异异构体。", options: { bold: true } },
{ text: "若要追「短 IL-18」，正确的落点是蛋白层面的非经典切割（caspase-3、糜蛋白酶、蛋白酶-3、granzyme B 产物与 caspase-1 不同），转录组看不见，只能做蛋白实验。", options: {} }],
  { isTextBox: true, x: .8, y: 5.55, w: 11.7, h: .9, fontFace: F, fontSize: 13.5, color: INK, lineSpacingMultiple: 1.2 });
cap(s, "数据源 Snaptron（服务 recount/SRA 与 GTEx 的逐样本剪接连接计数）；坐标由 Ensembl 转录本结构推出并与 Snaptron 标注连接核对。");

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
