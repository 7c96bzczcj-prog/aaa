// 方向汇报 PPT — built to spec v1.0
// HARD RULE (spec §0): no simulated, example, or invented data anywhere.
// Any figure whose data is not in this environment becomes a spec'd placeholder box.
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');

// ---- palette (spec §3). Hex without '#', no alpha. ----
const C = {
  bg: 'FFFFFF',
  title: '1F3864',
  sub: '595959',
  lab: '2E5C9A',      // 本实验室数据
  cite: 'A6A6A6',     // 引用数据
  accent: '0F4C81',   // only P4 and P7, once each
  grid: 'D9D9D9',
};
const FONT = 'Microsoft YaHei';
const M = 0.6;              // margin
const W = 13.333, H = 7.5;
const CW = W - 2 * M;       // content width

const pres = new PptxGenJS();
pres.layout = 'LAYOUT_WIDE';   // MUST precede any addSlide
pres.author = 'Richard Tang';
pres.title = '在肿瘤到达之前装好 NK';

const placeholders = [];

function titleBar(s, text) {
  s.addText(text, { x: M, y: 0.62, w: CW, h: 0.72, fontSize: 32, bold: true,
                    color: C.title, fontFace: FONT, align: 'left', valign: 'middle' });
}

function sourceTag(s, text, isLab) {
  s.addText(text, { x: W - M - 4.6, y: H - 0.95, w: 4.6, h: 0.3, fontSize: 10,
                    color: isLab ? C.lab : C.sub, fontFace: FONT, align: 'right',
                    italic: !isLab });
}

// spec'd placeholder: light-grey dashed frame, NO data marks
function placeholder(s, o) {
  placeholders.push({ slide: o.slideNo, title: o.title, need: o.files });
  s.addShape(pres.ShapeType.rect, {
    x: o.x, y: o.y, w: o.w, h: o.h, fill: { color: 'FFFFFF' },
    line: { color: C.grid, width: 1.25, dashType: 'dash' },
  });
  const lines = [
    { text: o.title + '\n', options: { fontSize: 14, bold: true, color: C.sub, fontFace: FONT } },
    { text: '横轴：' + o.xaxis + '\n', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
    { text: '纵轴：' + o.yaxis + '\n', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
    { text: '所需数据：' + o.files + '\n', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
    { text: '所需统计量：' + o.stats + '\n\n', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
    { text: '[待补：数据不在本环境]', options: { fontSize: 12, bold: true, color: '8C8C8C', fontFace: FONT } },
  ];
  s.addText(lines, { x: o.x + 0.3, y: o.y + 0.25, w: o.w - 0.6, h: o.h - 0.5,
                     align: 'left', valign: 'top', lineSpacingMultiple: 1.25 });
}

// ============================== P1 封面 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  s.addText('在肿瘤到达之前装好 NK', { x: M, y: 2.45, w: CW, h: 1.1, fontSize: 44, bold: true,
    color: C.title, fontFace: FONT, align: 'left' });
  s.addText('研究方向汇报', { x: M, y: 3.65, w: CW, h: 0.45, fontSize: 18,
    color: C.sub, fontFace: FONT, align: 'left' });
  s.addText('唐泽宇 (Richard Tang) · 沈俊辰组 · 2026-08', { x: M, y: H - 1.15, w: CW, h: 0.35,
    fontSize: 12, color: C.sub, fontFace: FONT, align: 'left' });
  s.addNotes('这次讲一条线：NK 在实体瘤里之所以没用，可能不是分子问题，而是到达顺序问题。' +
    '十分钟讲现状和转折，五分钟讲提案。');
}

// ============================== P2 NK 极少 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, 'NK 在实体瘤中极少');
  s.addText([
    { text: '4.0%', options: { fontSize: 60, bold: true, color: C.lab, fontFace: FONT } },
  ], { x: M, y: 1.68, w: 3.6, h: 1.42, align: 'left' });
  s.addText('肿瘤内 NK 占免疫细胞比例\n癌旁 7.4%（CD45⁺ 富集后）', {
    x: M, y: 3.02, w: 3.6, h: 0.9, fontSize: 15, color: C.title, fontFace: FONT,
    align: 'left', lineSpacingMultiple: 1.3 });
  s.addText('未分选组织中约 1%（GSE178341）', {
    x: M, y: 3.98, w: 3.6, h: 0.4, fontSize: 11, color: C.sub, fontFace: FONT, align: 'left' });

  // real values computed from results/cell_counts_by_group.csv (single-mode libraries, 28 patients)
  const cats = ['髓系', 'B', '未分配', 'CD4T', 'CD8T', 'NK'];
  const vals = [28.7, 20.0, 18.8, 17.8, 10.7, 4.0];
  s.addChart(pres.ChartType.bar, [{ name: '肿瘤内占比', labels: cats, values: vals }], {
    x: 4.6, y: 1.55, w: CW - 4.0, h: 4.5,
    barDir: 'col', showLegend: false, showValue: true,
    dataLabelColor: C.title, dataLabelFontSize: 11, dataLabelFontFace: FONT,
    dataLabelFormatCode: '0.0"%"',
    chartColors: [C.cite, C.cite, C.cite, C.cite, C.cite, C.lab],
    catAxisLabelColor: C.sub, catAxisLabelFontSize: 12, catAxisLabelFontFace: FONT,
    valAxisLabelColor: C.sub, valAxisLabelFontSize: 11, valAxisLabelFontFace: FONT,
    valAxisMaxVal: 35, valGridLine: { color: C.grid, style: 'solid', size: 1 },
    catAxisLineShow: true, valAxisLineShow: false,
    catAxisLineColor: C.grid, valAxisLineColor: C.grid, serAxisLineColor: C.grid, barGapWidthPct: 60,
  });
  sourceTag(s, '本实验室分析 · GSE154826 (Leader NSCLC)，n = 28 例配对病人', true);
  s.addNotes('这是我们自己再分析的 GSE154826，CD45 阳性富集数据。' +
    'NK 只占肿瘤内免疫细胞的 4%，癌旁 7.4%，进瘤以后掉一半。' +
    '注意这已经是富集过的；未分选组织里 NK 大概只有 1%。' +
    '“未分配”那一条是没能明确归类的细胞，如实列出。');
}

// ============================== P3 进瘤后快速失能 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '进入肿瘤后 24–72 小时内同时丢失杀伤与趋化因子');
  s.addText('不是慢性耗竭，是快速、同步的功能塌陷。', {
    x: M, y: 1.48, w: CW, h: 0.4, fontSize: 16, color: C.title, fontFace: FONT });
  placeholder(s, {
    slideNo: 3,
    title: 'Dean 2024 时间轴：5 h / 24 h / 72 h',
    xaxis: '进入肿瘤后时间（5 h / 24 h / 72 h）',
    yaxis: 'CD11b⁺CD49a⁻ → CD11b⁻CD49a⁺ 比例；Ccl3 / Ccl4 / Ccl5 与细胞毒同步下降',
    files: 'Dean I, et al. Nat Commun 2024;15:683 图内数值（本环境无）',
    stats: '各时间点两群比例 ± SD，n（动物）',
    x: M, y: 1.95, w: CW, h: 3.7,
  });
  sourceTag(s, '引用 · Dean I, et al. Nat Commun 2024;15:683', false);
  s.addNotes('这一页是引用 Dean 2024，不是我们的数据，所以是灰色标签。' +
    '他们做的是活体追踪：NK 进入肿瘤后 24 到 72 小时内，表型从 CD11b 阳性转成 CD49a 阳性，' +
    '同时趋化因子和杀伤能力一起掉。关键是“同时”——这排除了单纯的杀伤耗竭解释。' +
    '数值我没有从原文取，所以这里是占位框，补图时按框里列的要求填。');
}

// ============================== P4 转折页 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '从已成型肿瘤中清除 NK，肿瘤生长不变');
  s.addText('NK 到得太晚，晚到有没有它都一样。', {
    x: M, y: 1.45, w: CW, h: 0.55, fontSize: 22, bold: true, color: C.accent, fontFace: FONT });
  placeholder(s, {
    slideNo: 4,
    title: '植瘤前敲除 vs 植瘤后敲除，各配对照',
    xaxis: '时间（天）或终点分组',
    yaxis: '肿瘤体积（mm³）／终点重量（g）',
    files: '本实验室未发表数据（本环境无）',
    stats: '分组、时间点、肿瘤体积/重量、n（动物）；配对对照',
    x: M, y: 2.05, w: CW, h: 3.55,
  });
  s.addText('Dean 2024 在多个模型中得到一致结果', {
    x: M, y: 5.75, w: 7.0, h: 0.3, fontSize: 10, color: C.sub, fontFace: FONT, italic: true });
  sourceTag(s, '本实验室未发表数据', true);
  s.addNotes('这是全场第一个支点。如果在肿瘤已经长起来之后再敲除 NK，生长曲线不动；' +
    '而在植瘤之前敲除则有差别。这说明晚到的 NK 不承担控制作用。' +
    'Dean 2024 在多个模型里也是一致的，所以不是我们模型的特例。' +
    '数据文件不在分析环境里，这里是占位框，需要分组/时间点/体积/n。');
}

// ============================== P5 蜕膜 NK ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '蜕膜 NK：不杀伤，但趋化因子满载');
  s.addText('低杀伤，但保有招募功能。', {
    x: M, y: 1.48, w: CW, h: 0.4, fontSize: 16, color: C.title, fontFace: FONT });

  const genes = ['CCL4', 'CCL3', 'XCL2', 'XCL1', 'CCL5'];
  const rates = [94.4, 81.7, 80.1, 75.6, 58.5];
  s.addChart(pres.ChartType.bar, [{ name: '检出率', labels: genes, values: rates }], {
    x: M, y: 2.0, w: 7.4, h: 3.5,
    barDir: 'bar', showLegend: false, showValue: true,
    dataLabelColor: C.title, dataLabelFontSize: 11, dataLabelFontFace: FONT,
    dataLabelFormatCode: '0.0"%"',
    chartColors: [C.lab],
    catAxisLabelColor: C.sub, catAxisLabelFontSize: 12, catAxisLabelFontFace: FONT,
    valAxisLabelColor: C.sub, valAxisLabelFontSize: 11, valAxisLabelFontFace: FONT,
    valAxisMaxVal: 100, valGridLine: { color: C.grid, style: 'solid', size: 1 },
    catAxisLineShow: true, valAxisLineShow: false,
    catAxisLineColor: C.grid, valAxisLineColor: C.grid, serAxisLineColor: C.grid, barGapWidthPct: 55,
  });
  s.addText([
    { text: '此处为跨亚群均值\n\n', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
    { text: '未作拟时序／轨迹图\n', options: { fontSize: 12, bold: true, color: C.title, fontFace: FONT } },
    { text: '（原论文无排序声明）', options: { fontSize: 11, color: C.sub, fontFace: FONT } },
  ], { x: 8.3, y: 2.1, w: CW - 7.7, h: 3.2, align: 'left', valign: 'top', lineSpacingMultiple: 1.3 });
  sourceTag(s, '本实验室分析 · VT2018 / E-MTAB-6701', true);
  s.addNotes('蜕膜 NK 是这个状态天然稳定存在的证明：它几乎不杀伤，但趋化因子基因检出率很高，' +
    'CCL4 到 94%。所以“低杀伤”和“无功能”不是一回事。' +
    '按亚群分列的值我这里没有，画的是跨亚群均值，已在页面上注明。' +
    '刻意没有画拟时序图——定义这三个亚群的原文根本没做过排序声明，后面三篇的排序互相矛盾。');
}

// ============================== P6 经得起检验 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '这条结论经得起簇内检验');
  s.addText('效应与反号对照来自同一批细胞。', {
    x: M, y: 1.48, w: CW, h: 0.4, fontSize: 16, color: C.title, fontFace: FONT });

  placeholder(s, {
    slideNo: 6,
    title: '左：按供者配对的 XCL1 检出率',
    xaxis: 'dNK1 → dNK2（配对）',
    yaxis: 'XCL1 检出率（%）',
    files: '6 名供者的按亚群检出率（本环境无）',
    stats: '每供者一条连线，n = 6（供者）',
    x: M, y: 1.95, w: 5.9, h: 3.6,
  });

  const items = ['XCL1', 'XCL2', 'CXCR4'];
  const dpp = [6.7, 5.7, -3.5];
  s.addChart(pres.ChartType.bar, [{ name: '簇内差值', labels: items, values: dpp }], {
    x: 7.0, y: 1.95, w: CW - 6.4, h: 3.6,
    barDir: 'col', showLegend: false, showValue: true,
    dataLabelColor: C.title, dataLabelFontSize: 12, dataLabelFontFace: FONT,
    dataLabelFormatCode: '+0.0;-0.0',
    chartColors: [C.lab, C.lab, C.cite],
    catAxisLabelColor: C.sub, catAxisLabelFontSize: 12, catAxisLabelFontFace: FONT,
    valAxisLabelColor: C.sub, valAxisLabelFontSize: 11, valAxisLabelFontFace: FONT,
    valAxisMinVal: -6, valAxisMaxVal: 9,
    valGridLine: { color: C.grid, style: 'solid', size: 1 },
    catAxisLineShow: true, valAxisLineShow: false,
    catAxisLineColor: C.grid, valAxisLineColor: C.grid, serAxisLineColor: C.grid, barGapWidthPct: 70,
  });
  s.addText('dNK1 簇内 CD39⁺ → CD39⁻，百分点差，3/3 供者一致；CXCR4 反号内对照', {
    x: 7.0, y: 5.62, w: CW - 6.4, h: 0.35, fontSize: 10, color: C.sub, fontFace: FONT });
  sourceTag(s, '本实验室分析', true);
  s.addNotes('右边是关键：在 dNK1 这一个簇内部再分 CD39 阳性和阴性，XCL1 涨 6.7 个百分点、' +
    'XCL2 涨 5.7，三个供者方向都一致；而 CXCR4 在同一批细胞上是反号的，掉 3.5 个百分点。' +
    '反号对照来自同一批细胞，这排除了整体检出深度差异造成的假象。' +
    '左边按供者配对的图数据不在环境里，是占位框。');
}

// ============================== P7 关键页 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '差别在顺序');
  s.addText('不是分子不同，是顺序不同。', {
    x: M, y: 1.45, w: CW, h: 0.7, fontSize: 26, bold: true, color: C.accent, fontFace: FONT });

  const axY = [2.7, 4.6];
  const labels = ['蜕膜', '肿瘤'];
  const shades = ['2E5C9A', '7FA3CC'];   // same blue family, two depths
  const steps = [
    ['uNK 群体先建立', '滋养层侵袭到达', 'NK 保持功能'],
    ['肿瘤已成型', 'NK 逐个进入', '24–72 h 内失能'],
  ];
  for (let i = 0; i < 2; i++) {
    s.addText(labels[i], { x: M, y: axY[i] - 0.05, w: 1.0, h: 0.4, fontSize: 15, bold: true,
      color: C.title, fontFace: FONT, align: 'left' });
    s.addShape(pres.ShapeType.line, { x: M + 1.05, y: axY[i] + 0.16, w: CW - 1.35, h: 0,
      line: { color: shades[i], width: 2.5, endArrowType: 'triangle' } });
    for (let k = 0; k < 3; k++) {
      const cx = M + 1.55 + k * 3.5;
      s.addShape(pres.ShapeType.ellipse, { x: cx - 0.09, y: axY[i] + 0.07, w: 0.18, h: 0.18,
        fill: { color: shades[i] }, line: { color: shades[i], width: 1 } });
      s.addText(steps[i][k], { x: cx - 1.5, y: axY[i] + 0.34, w: 3.0, h: 0.5, fontSize: 12,
        color: C.sub, fontFace: FONT, align: 'center' });
    }
  }
  s.addText('时间 →', { x: W - M - 1.2, y: 5.65, w: 1.2, h: 0.3, fontSize: 11,
    color: C.sub, fontFace: FONT, align: 'right' });
  s.addNotes('这是第二个支点，也是整条线的转轴。' +
    '蜕膜里，NK 群体是在滋养层侵袭之前就建好的；肿瘤里，NK 是在肿瘤已经成型之后一个一个进去的。' +
    '同样的分子环境，先到和后到结果完全不同。' +
    '这一页是示意图，不是数据图，我不会拿它当证据。');
}

// ============================== P8 提案 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '提案：在转移靶器官预先装好这个状态');
  const steps = ['手术切除原发灶', '在转移靶器官（肝／肺）\n预装驻留 NK 状态', '肿瘤细胞到达时\n遭遇已就位的招募型 NK'];
  for (let i = 0; i < 3; i++) {
    const x = M + i * 4.15;
    s.addShape(pres.ShapeType.roundRect, { x, y: 2.2, w: 3.6, h: 1.9,
      fill: { color: 'FFFFFF' }, line: { color: C.lab, width: 1.5 }, rectRadius: 0.08 });
    s.addText(steps[i], { x: x + 0.2, y: 2.4, w: 3.2, h: 1.5, fontSize: 14, color: C.title,
      fontFace: FONT, align: 'center', valign: 'middle', lineSpacingMultiple: 1.3 });
    if (i < 2) s.addShape(pres.ShapeType.line, { x: x + 3.72, y: 3.15, w: 0.3, h: 0,
      line: { color: C.lab, width: 2, endArrowType: 'triangle' } });
  }
  s.addText([
    { text: '目标不是杀伤，是分泌 CCL5 / XCL1 招募 cDC1。\n', options: { fontSize: 17, bold: true, color: C.title, fontFace: FONT } },
    { text: '装的是：趋化因子满载、抗 TGF-β 沉默的驻留 NK 状态。', options: { fontSize: 14, color: C.sub, fontFace: FONT } },
  ], { x: M, y: 4.6, w: CW, h: 1.0, align: 'left', lineSpacingMultiple: 1.35 });
  s.addNotes('提案是把顺序倒过来：不在肿瘤里补 NK，而是在肿瘤还没到的靶器官里先把状态装好。' +
    '要特别说清楚，我们要的不是杀伤能力，是招募能力——分泌 CCL5 和 XCL1 把 cDC1 拉进来。' +
    '措辞统一为“趋化因子满载、抗 TGF-β 沉默的组织驻留 NK 状态”，不用任何借代说法。');
}

// ============================== P9 两条 Aim ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '两条 Aim 是同一个问题的正反两面');
  const boxes = [
    { t: 'Aim 1', b: '一个被诱导的状态，\n靠什么在信号撤除后持久？' },
    { t: 'Aim 2', b: '要做成这件事，\n需要一个持久的被诱导状态。' },
  ];
  for (let i = 0; i < 2; i++) {
    const x = M + i * 6.33;
    s.addShape(pres.ShapeType.roundRect, { x, y: 2.3, w: 5.4, h: 2.3,
      fill: { color: 'FFFFFF' }, line: { color: C.lab, width: 1.5 }, rectRadius: 0.08 });
    s.addText(boxes[i].t, { x: x + 0.3, y: 2.55, w: 4.8, h: 0.45, fontSize: 17, bold: true,
      color: C.lab, fontFace: FONT });
    s.addText(boxes[i].b, { x: x + 0.3, y: 3.1, w: 4.8, h: 1.2, fontSize: 15, color: C.title,
      fontFace: FONT, lineSpacingMultiple: 1.35 });
  }
  s.addShape(pres.ShapeType.line, { x: M + 5.5, y: 3.45, w: 0.9, h: 0,
    line: { color: C.sub, width: 2, beginArrowType: 'triangle', endArrowType: 'triangle' } });
  s.addText('同一个问题，符号相反。', { x: M, y: 5.0, w: CW, h: 0.5, fontSize: 18, bold: true,
    color: C.title, fontFace: FONT, align: 'center' });
  s.addNotes('两条 Aim 不是并列的两个课题，是一个问题的两面：' +
    'Aim 1 问诱导态为什么能持久，Aim 2 需要一个能持久的诱导态。机制和应用互为条件。');
}

// ============================== P10 第一个实验 ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, '第一个实验');
  const flow = ['预处理靶器官 NK', '尾静脉注射肿瘤细胞', '读出'];
  for (let i = 0; i < 3; i++) {
    const x = M + i * 4.15;
    s.addShape(pres.ShapeType.roundRect, { x, y: 1.85, w: 3.6, h: 1.0,
      fill: { color: 'FFFFFF' }, line: { color: C.lab, width: 1.5 }, rectRadius: 0.08 });
    s.addText(flow[i], { x: x + 0.15, y: 1.95, w: 3.3, h: 0.8, fontSize: 15, color: C.title,
      fontFace: FONT, align: 'center', valign: 'middle' });
    if (i < 2) s.addShape(pres.ShapeType.line, { x: x + 3.72, y: 2.35, w: 0.3, h: 0,
      line: { color: C.lab, width: 2, endArrowType: 'triangle' } });
  }
  const reads = [
    ['结局', '转移结节计数'],
    ['机制 · 状态', '驻留 NK 的 CCL5 / XCL1 保留'],
    ['机制 · 下游', '瘤内 cDC1 数量'],
  ];
  for (let i = 0; i < 3; i++) {
    const y = 3.35 + i * 0.85;
    s.addText(reads[i][0], { x: M, y, w: 2.2, h: 0.6, fontSize: 12, color: C.sub,
      fontFace: FONT, valign: 'middle' });
    s.addText(reads[i][1], { x: M + 2.5, y, w: 8.0, h: 0.6, fontSize: 15, color: C.title,
      fontFace: FONT, valign: 'middle' });
  }
  s.addText('一个实验同时测机制与结局。', { x: M, y: 6.05, w: CW, h: 0.4, fontSize: 15,
    bold: true, color: C.title, fontFace: FONT });
  s.addNotes('第一个实验刻意设计成一次同时给机制和结局：' +
    '结节数是结局，CCL5/XCL1 保留度说明状态有没有装住，cDC1 数量说明招募链路有没有通。' +
    '如果结节少了但 cDC1 没变，那机制解释就不成立，这个设计能让我们知道。');
}

// ============================== Backup ==============================
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, 'Backup B1 · CXCR4 的亚群梯度在簇内反号');
  s.addText('反号对照排除检出深度假象。', {
    x: M, y: 1.45, w: CW, h: 0.4, fontSize: 15, color: C.title, fontFace: FONT });
  placeholder(s, {
    slideNo: 'B1',
    title: 'CXCR4：亚群间梯度 vs 簇内差值',
    xaxis: '亚群间（dNK1/2/3）／dNK1 簇内（CD39⁺ → CD39⁻）',
    yaxis: 'CXCR4 检出率（%）与百分点差',
    files: '亚群间 CXCR4 梯度（本环境无）；簇内值已知 = −3.5 pp',
    stats: '亚群均值 ± SD；簇内配对差，n = 3（供者）',
    x: M, y: 1.95, w: CW, h: 3.7,
  });
  sourceTag(s, '本实验室分析', true);
  s.addNotes('被问到“会不会是检出深度差异”时翻这一页。CXCR4 在亚群之间是一个方向，' +
    '在 dNK1 簇内部是反方向，−3.5 个百分点。同一批细胞、同一套深度，方向相反，' +
    '所以不是深度造成的。');
}
{
  const s = pres.addSlide();
  s.background = { color: C.bg };
  titleBar(s, 'Backup B2 · 可测性分带');
  s.addText('27 个候选中 14 个在该深度不可测。', {
    x: M, y: 1.45, w: CW, h: 0.45, fontSize: 16, color: C.title, fontFace: FONT });
  s.addChart(pres.ChartType.bar, [{ name: '基因数', labels: ['可测', '该深度下不可测'], values: [13, 14] }], {
    x: M, y: 2.1, w: 7.0, h: 3.3,
    barDir: 'col', showLegend: false, showValue: true,
    dataLabelColor: C.title, dataLabelFontSize: 13, dataLabelFontFace: FONT,
    chartColors: [C.lab, C.cite],
    catAxisLabelColor: C.sub, catAxisLabelFontSize: 12, catAxisLabelFontFace: FONT,
    valAxisLabelColor: C.sub, valAxisLabelFontSize: 11, valAxisLabelFontFace: FONT,
    valAxisMaxVal: 20, valGridLine: { color: C.grid, style: 'solid', size: 1 },
    catAxisLineShow: true, valAxisLineShow: false,
    catAxisLineColor: C.grid, valAxisLineColor: C.grid, serAxisLineColor: C.grid, barGapWidthPct: 80,
  });
  s.addText('不可测 = 对照组检出率贴近 0%。\n报为“不可测”，不报为“无差异”。', {
    x: 7.9, y: 2.4, w: CW - 7.3, h: 2.0, fontSize: 14, color: C.sub, fontFace: FONT,
    lineSpacingMultiple: 1.35 });
  sourceTag(s, '本实验室分析', true);
  s.addNotes('被问“为什么只报这几个基因”时翻这一页。27 个候选里有 14 个在这个深度下贴地板，' +
    '不是没差异，是测不出来。我们把它们标成不可测，不标成阴性——这两件事不一样。');
}

pres.writeFile({ fileName: '方向汇报_NK预装.pptx' }).then(() => {
  console.log('WROTE 方向汇报_NK预装.pptx');
  console.log('placeholders:', JSON.stringify(placeholders, null, 1));
  fs.writeFileSync('/tmp/placeholders.json', JSON.stringify(placeholders, null, 1));
});
