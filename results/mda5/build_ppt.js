const pptxgen = require('pptxgenjs');
const fs = require('fs');
const FIG = '/home/user/aaa/results/mda5/figs/';
const BLACK = '000000', GREY = '6E6E6E', RULE = 'C8CCD2';
const CN = 'Microsoft YaHei', EN = 'Calibri';

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';                  // 13.333 x 7.5 英寸
pres.title = 'MDA5（IFIH1）在人脑与多发性硬化中的表达';

/* 读 PNG 的 IHDR 取宽高比，图片按页面可用区等比铺满，不拉伸 */
const DIMS = {};
const sizeOf = p => {
  const b = fs.readFileSync(p);
  let i = 8;
  while (i < b.length) {
    const len = b.readUInt32BE(i);
    if (b.toString('ascii', i + 4, i + 8) === 'IHDR')
      return { w: b.readUInt32BE(i + 8), h: b.readUInt32BE(i + 12) };
    i += len + 12;
  }
  return null;
};
['fig_q1_celltype.png', 'fig_q1_region.png', 'figB_null.png', 'fig_q2_wm.png',
 'fig_q2_gm.png', 'fig_q2_interaction.png', 'fig_confound.png',
 'fig_umap_wm.png'].forEach(n => {
  try { const s = sizeOf(FIG + n); if (s) DIMS[n] = s.w / s.h; } catch (e) {}
});

let pageNo = 0;
function plain() {
  const s = pres.addSlide();
  pageNo++;
  s.addText(String(pageNo), { x: 12.5, y: 6.95, w: 0.5, h: 0.3, fontSize: 10,
    color: GREY, align: 'right', fontFace: EN, isTextBox: true, margin: 0 });
  return s;
}
function title(s, t, sub) {
  s.addText(t, { x: 0.6, y: 0.30, w: 12.1, h: 0.52, fontSize: 25, bold: true,
    color: BLACK, fontFace: CN, isTextBox: true, margin: 0 });
  if (sub) s.addText(sub, { x: 0.6, y: 0.86, w: 12.1, h: 0.36, fontSize: 13.5,
    color: GREY, fontFace: CN, isTextBox: true, margin: 0 });
}
function figure(s, name, top) {
  const a = DIMS[name];
  if (!a) { console.warn('缺图', name); return; }
  const bx = 0.45, by = top, bw = 12.45, bh = 7.5 - top - 0.42;
  let w = bw, h = w / a;
  if (h > bh) { h = bh; w = h * a; }
  s.addImage({ path: FIG + name, x: bx + (bw - w) / 2, y: by + (bh - h) / 2, w, h });
}
function table(s, rows, colW, y, rowH = 0.42) {
  s.addTable(rows, { x: 0.6, y, w: 12.1, colW, fontSize: 12.5, fontFace: CN, color: BLACK,
    border: [{ type: 'none' }, { type: 'none' }, { type: 'solid', color: RULE, pt: 1 }, { type: 'none' }],
    rowH: rowH, valign: 'middle', autoPage: false });
}
const B = t => ({ text: t, options: { bold: true } });

/* ---------------------------------------------------------------- 1 标题 */
let s = pres.addSlide();
s.addText('MDA5（IFIH1）在人脑与多发性硬化中的表达',
  { x: 0.9, y: 2.72, w: 11.5, h: 0.8, fontSize: 34, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });
s.addText('单核转录组水平的细胞类型特异性分析',
  { x: 0.9, y: 3.64, w: 11.5, h: 0.5, fontSize: 17, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0 });

/* ------------------------------------------------------ 2 两个问题与数据 */
s = plain();
title(s, '两个问题');
s.addText('一　正常人脑中 IFIH1 由哪些细胞类型表达？两个解剖区域是否不同？',
  { x: 0.75, y: 1.06, w: 12, h: 0.40, fontSize: 17, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });
s.addText('二　MS 病变中 IFIH1 的细胞分布如何改变？两个区域的改变是否相同？',
  { x: 0.75, y: 1.56, w: 12, h: 0.40, fontSize: 17, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });
table(s, [
  [B('队列'), B('区室'), B('供体'), B('核数'), B('承担')],
  ['CELLxGENE Census', '正常脑', '783', '18 551 076', '细胞类型表达谱'],
  ['GSE180759  Absinta 2021', '皮层下白质', '3 对照 + 5 MS', '66 432', '白质，含病灶分期'],
  ['GSE279180  Lerma-Martin 2024', '皮层下白质', '6 对照 + 7 MS', '103 780', '白质，主队列'],
  ['Schirmer 2019', '皮层标本', '9 对照 + 12 MS', '48 919', '皮层'],
  ['GSE118257  Jäkel 2019', '白质', '5 对照 + 4 MS', '17 799', '白质，验证'],
], [3.3, 3.1, 2.1, 1.6, 2.0], 2.20, 0.64);

/* ------------------------------------------- 3 问题一：哪些细胞表达 */
s = plain();
title(s, 'IFIH1 在人脑各类细胞中的表达');
figure(s, 'fig_q1_celltype.png', 1.05);

/* ------------------------------------------- 4 问题一：灰质与白质 */
s = plain();
title(s, '正常脑：两个解剖区域无差异');
figure(s, 'fig_q1_region.png', 1.05);

/* ------------------------------------------- 5 第 4 页的批次对照 */
s = plain();
title(s, '跨研究偏移的稳健性检验');
figure(s, 'figB_null.png', 1.05);

/* ------------------------------------------- 6 问题二：白质中的改变 */
s = plain();
title(s, 'MS 白质：少突胶质细胞 IFIH1 阳性率升高');
figure(s, 'fig_q2_wm.png', 1.05);

/* ------------------------------------------- 7 问题二：灰质中的改变 */
s = plain();
title(s, 'MS 皮层：未见同等幅度的变化');
figure(s, 'fig_q2_gm.png', 1.05);

/* ------------------------------------------- 8 问题二的落点 */
s = plain();
title(s, '两个区域变化幅度之差');
figure(s, 'fig_q2_interaction.png', 1.05);

/* ------------------------------------------- 9 皮层块的混入 */
s = plain();
title(s, '皮层标本含皮层下白质，MS 取材更深');
figure(s, 'fig_confound.png', 1.05);

/* ------------------------------------------- 10 白质图谱 */
s = plain();
title(s, '白质单核图谱与脱髓鞘病灶分期');
figure(s, 'fig_umap_wm.png', 1.05);

/* ------------------------------------------- 10 与文献的关系 */
s = plain();
title(s, '与已发表结果的关系');
table(s, [
  [B('本研究结果'), B('文献状态')],
  ['少突胶质在白质病灶上调', 'Lerma-Martin 2024 补充表 6 已报道：CA log2FC +0.905（padj 0.012），CI +1.095（padj 0.005）'],
  ['同一队列髓系未达显著', '与原作者一致，其补充表髓系 padj = 0.116'],
  ['bulk 白质病灶同向上调', 'GSE138614 活动病灶 FDR 0.021、慢性活动 FDR 0.007；GSE283092 泡沫样活动病灶 1.79 倍'],
  ['皮层组织块的 IFIH1', 'Schirmer 2019 补充表中 IFIH1 低于其检出阈值，本次为原创测量'],
  ['GSE180759 中的 IFIH1', '原文 13 个补充表中未出现'],
  ['两个区室改变的直接对比', '未见任何已发表研究做过'],
  ['MDA5 蛋白', '人类 MS 脑组织的免疫组化与蛋白质组测量为零篇'],
], [3.4, 8.7], 1.05, 0.72);

/* ------------------------------------------- 11 限制 */
s = plain();
title(s, '局限');
const lim = [
  ['皮层标本含白质成分', 'Schirmer 取材跨越皮层全层并带有其下的皮层下白质，MS 标本取材更深，疾病对比与取材部位混杂。'],
  ['区域交互仅两类可估计', '小胶质细胞、神经元、内皮细胞因一侧样本量不足，少突胶质前体细胞因对照组未检出。'],
  ['白质结论主要由一个队列支持', 'GSE279180 十三例供体，3.5 倍，区间 [1.5, 4.9]；GSE180759 仅三例对照，1.5 倍，区间包含 1。'],
  ['测序深度残余影响', '已限定在每核 1500 至 5000 UMI 的共同深度区间内比较，区间内仍存在残余差异。'],
  ['GSE180759 存在平台混杂', '病灶周围白质文库全为 NovaSeq，病灶边缘文库全为 HiSeq；图谱已按文库做 Harmony 整合。'],
  ['以阳性核比例而非表达量', 'Schirmer 仅提供对数归一化矩阵，标度因子不可恢复，跨队列表达量不可比。'],
  ['仅限转录本水平', '人类 MS 脑组织尚无 MDA5 蛋白水平的测量。'],
];
let yy = 1.10;
lim.forEach(L => {
  s.addText(L[0], { x: 0.6, y: yy, w: 3.9, h: 0.40, fontSize: 13.5, bold: true,
    color: BLACK, fontFace: CN, isTextBox: true, margin: 0 });
  s.addText(L[1], { x: 4.65, y: yy, w: 8.1, h: 0.62, fontSize: 12, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0, lineSpacing: 16 });
  yy += 0.80;
});

pres.writeFile({ fileName: '/home/user/aaa/results/mda5/MDA5_IFIH1_人脑与MS.pptx' })
  .then(f => console.log('written', f, '|', pageNo + 1, '页'));
