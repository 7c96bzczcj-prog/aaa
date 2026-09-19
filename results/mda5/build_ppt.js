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
 'fig_q2_gm.png', 'fig_q2_interaction.png', 'fig_umap_wm.png'].forEach(n => {
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
s.addText('全篇测量同一个量：IFIH1 的表达。细胞类型、灰质与白质、正常与 MS，只是分格的三条轴。',
  { x: 0.9, y: 3.64, w: 11.5, h: 0.5, fontSize: 16, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0 });

/* ------------------------------------------------------ 2 两个问题与数据 */
s = plain();
title(s, '两个问题');
s.addText('一　正常人脑中 IFIH1 由哪些细胞表达？灰质与白质是否不同？',
  { x: 0.75, y: 1.06, w: 12, h: 0.40, fontSize: 17, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });
s.addText('二　MS 中 IFIH1 的分布如何改变？灰质与白质的改变是否相同？',
  { x: 0.75, y: 1.56, w: 12, h: 0.40, fontSize: 17, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });
s.addText('第一问的区域对比跨研究，需要批次对照；第二问的「改变」在各研究内部相减，批次自行抵消。',
  { x: 0.75, y: 2.06, w: 12, h: 0.36, fontSize: 13, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0 });
table(s, [
  [B('队列'), B('区室'), B('供体'), B('核数'), B('承担')],
  ['CELLxGENE Census', '正常脑', '783', '18 551 076', '问题一：细胞类型次序'],
  ['GSE180759  Absinta 2021', '皮层下白质', '3 对照 + 5 MS', '66 432', '白质基线与改变'],
  ['GSE279180  Lerma-Martin 2024', '皮层下白质', '6 对照 + 7 MS', '103 780', '白质基线与改变（主力）'],
  ['Schirmer 2019', '皮层灰质', '9 对照 + 12 MS', '48 919', '灰质基线与改变'],
  ['GSE118257  Jäkel 2019', '白质', '5 对照 + 4 MS', '17 799', '白质验证'],
], [3.5, 2.0, 2.3, 1.7, 2.6], 2.62, 0.62);

/* ------------------------------------------- 3 问题一：哪些细胞表达 */
s = plain();
title(s, '问题一：IFIH1 由哪些细胞表达',
  '正常脑内 IFIH1 最高的是内皮，其次小胶质；神经元最低。左为 783 位供体的伪整体排序，右为白质单核层面的检出。');
figure(s, 'fig_q1_celltype.png', 1.30);

/* ------------------------------------------- 4 问题一：灰质与白质 */
s = plain();
title(s, '问题一：正常人的灰质与白质',
  '只取对照供体。星形胶质在白质高于灰质，其余细胞类型无差别；小胶质与内皮因灰质对照供体不足而无法判定。');
figure(s, 'fig_q1_region.png', 1.30);

/* ------------------------------------------- 5 第 4 页的批次对照 */
s = plain();
title(s, '第 4 页的批次对照：把跨研究偏移扣掉',
  '偏移同等作用于所有基因，故全基因组 log2（白质/灰质）的中位即偏移量。扣掉它之后，只有星形胶质在两个白质队列中都高于中位。');
figure(s, 'figB_null.png', 1.44);

/* ------------------------------------------- 6 问题二：白质中的改变 */
s = plain();
title(s, '问题二：白质中，MS 相对对照的改变',
  '两个白质队列分开计算，每点一位供体。少突胶质在两队列同向升高，但区间不跨零的只有 GSE279180。');
figure(s, 'fig_q2_wm.png', 1.30);

/* ------------------------------------------- 7 问题二：灰质中的改变 */
s = plain();
title(s, '问题二：灰质中，MS 相对对照的改变',
  '皮层，Schirmer 队列。少突胶质与星形胶质的区间跨零。OPC 区间不跨零，但对照中位贴近零，倍数变化在此不稳。');
figure(s, 'fig_q2_gm.png', 1.30);

/* ------------------------------------------- 8 问题二的落点 */
s = plain();
title(s, '问题二：两个区室的改变是否相同',
  '画的是自助 95% 区间，不是显著性星号。区间跨零意味着判不出，而不是两侧相同。');
figure(s, 'fig_q2_interaction.png', 1.30);

/* ------------------------------------------- 9 白质图谱 */
s = plain();
title(s, '白质内部：细胞图谱与病灶分期',
  'GSE180759，66 432 个核，已按文库做 Harmony 整合。对照白质几乎不表达，慢性活动边缘在小胶质与少突胶质上点亮。');
figure(s, 'fig_umap_wm.png', 1.30);

/* ------------------------------------------- 10 与文献的关系 */
s = plain();
title(s, '与已发表结果的关系');
table(s, [
  [B('本次结果'), B('文献状态')],
  ['少突胶质在白质病灶上调', 'Lerma-Martin 2024 补充表 6 已报道：CA log2FC +0.905（padj 0.012），CI +1.095（padj 0.005）'],
  ['同一队列髓系未达显著', '与原作者一致，其补充表髓系 padj = 0.116'],
  ['bulk 白质病灶同向上调', 'GSE138614 活动病灶 FDR 0.021、慢性活动 FDR 0.007；GSE283092 泡沫样活动病灶 1.79 倍'],
  ['皮层灰质的 IFIH1', 'Schirmer 2019 补充表中 IFIH1 低于其检出阈值，本次为原创测量'],
  ['GSE180759 中的 IFIH1', '原文 13 个补充表中未出现'],
  ['灰质与白质改变的直接对比', '未见任何已发表研究做过'],
  ['MDA5 蛋白', '人类 MS 脑组织的免疫组化与蛋白质组测量为零篇'],
], [3.4, 8.7], 1.22, 0.70);

/* ------------------------------------------- 11 限制 */
s = plain();
title(s, '限制');
const lim = [
  ['交互项全部跨零', '现有供体数判不出灰白质的改变是否不同。这是功效不足，不是「两侧相同」。'],
  ['白质结论主要由一个队列承担', 'GSE279180 十三位供体，留一法 +1.16 到 +1.41；GSE180759 仅三位对照，去掉一位即从 +1.09 降到 +0.25。'],
  ['灰质对照几乎没采到小胶质', 'Schirmer 九位对照中七位的小胶质核数不超过 11 个，故小胶质的交互项无法计算。'],
  ['OPC 的灰质结果不稳', '对照中位仅 0.002 CP10k，倍数变化对伪计数敏感；置换检验与自助区间在此结论相反。'],
  ['灰质与白质来自不同研究', '影响第一问的区域对比，已用全基因组分布作对照；不影响第二问的研究内相减。'],
  ['GSE180759 存在文库批次', '文库 18 单独形成少突胶质孤岛，其检出率为其余少突胶质的三倍；图谱已做 Harmony 整合。'],
  ['全部为转录本层面', '人类 MS 脑组织无任何 MDA5 蛋白测量。'],
];
let yy = 1.22;
lim.forEach(L => {
  s.addText(L[0], { x: 0.6, y: yy, w: 3.9, h: 0.40, fontSize: 13.5, bold: true,
    color: BLACK, fontFace: CN, isTextBox: true, margin: 0 });
  s.addText(L[1], { x: 4.65, y: yy, w: 8.1, h: 0.62, fontSize: 12, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0, lineSpacing: 16 });
  yy += 0.80;
});

pres.writeFile({ fileName: '/home/user/aaa/results/mda5/MDA5_IFIH1_人脑与MS.pptx' })
  .then(f => console.log('written', f, '|', pageNo + 1, '页'));
