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
s.addText('全篇测量同一个量：单个核是否检出 IFIH1。细胞类型、区室、正常与 MS，只是分格的三条轴。',
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
s.addText('口径：每核是否检出 IFIH1，限定在三队列共有的测序深度区间（每核 1500 至 5000 UMI）内比较。该量在任何归一化下都精确，故三队列可比。',
  { x: 0.75, y: 2.06, w: 12, h: 0.36, fontSize: 13, color: GREY,
    fontFace: CN, isTextBox: true, margin: 0 });
table(s, [
  [B('队列'), B('区室'), B('供体'), B('核数'), B('承担')],
  ['CELLxGENE Census', '正常脑', '783', '18 551 076', '问题一：细胞类型次序'],
  ['GSE180759  Absinta 2021', '皮层下白质', '3 对照 + 5 MS', '66 432', '白质基线与改变'],
  ['GSE279180  Lerma-Martin 2024', '皮层下白质', '6 对照 + 7 MS', '103 780', '白质基线与改变（主力）'],
  ['Schirmer 2019', '皮层组织块（含皮层下白质）', '9 对照 + 12 MS', '48 919', '皮层侧基线与改变'],
  ['GSE118257  Jäkel 2019', '白质', '5 对照 + 4 MS', '17 799', '白质验证'],
], [3.3, 3.1, 2.1, 1.6, 2.0], 2.62, 0.62);

/* ------------------------------------------- 3 问题一：哪些细胞表达 */
s = plain();
title(s, '问题一：IFIH1 由哪些细胞表达',
  '正常脑内 IFIH1 最高的是内皮，其次小胶质；神经元最低。左为 783 位供体的伪整体排序，右为白质单核层面的检出。');
figure(s, 'fig_q1_celltype.png', 1.30);

/* ------------------------------------------- 4 问题一：灰质与白质 */
s = plain();
title(s, '问题一：正常人的两个区室',
  '只取对照供体。六类细胞在皮层下白质与皮层组织块之间均无显著差别。');
figure(s, 'fig_q1_region.png', 1.30);

/* ------------------------------------------- 5 第 4 页的批次对照 */
s = plain();
title(s, '第 4 页的批次对照：把跨研究偏移扣掉',
  '偏移同等作用于所有基因，故全基因组 log2（白质/灰质）的中位即偏移量。扣掉它之后，只有星形胶质在两个白质队列中都高于中位。');
figure(s, 'figB_null.png', 1.44);

/* ------------------------------------------- 6 问题二：白质中的改变 */
s = plain();
title(s, '问题二：皮层下白质中，MS 相对对照的改变',
  '两个队列分开计算，每点一位供体，括号内为倍数变化的自助 95% 区间。少突胶质两队列同向升高，GSE279180 达 3.5 倍且区间不跨 1。');
figure(s, 'fig_q2_wm.png', 1.30);

/* ------------------------------------------- 7 问题二：灰质中的改变 */
s = plain();
title(s, '问题二：皮层组织块中，MS 相对对照的改变',
  'Schirmer 队列。六类细胞的倍数变化区间全部跨 1，无一类达到白质中少突胶质的幅度。');
figure(s, 'fig_q2_gm.png', 1.30);

/* ------------------------------------------- 8 问题二的落点 */
s = plain();
title(s, '问题二：两个区室的改变是否相同',
  '画的是自助 95% 区间，不是显著性星号。六类中只有两类两侧都可估计，其余因供体不足或对照中位为零而不可判定。');
figure(s, 'fig_q2_interaction.png', 1.30);

/* ------------------------------------------- 9 皮层块的混入 */
s = plain();
title(s, '为什么皮层那一侧只能说到这里',
  'Schirmer 的组织块按原文方法跨越整个皮层并带有其下的皮层下白质。MS 块取得更深，白质成分显著多于对照块，故该侧的疾病对比与解剖深度混杂。');
figure(s, 'fig_confound.png', 1.44);

/* ------------------------------------------- 10 白质图谱 */
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
  ['皮层组织块的 IFIH1', 'Schirmer 2019 补充表中 IFIH1 低于其检出阈值，本次为原创测量'],
  ['GSE180759 中的 IFIH1', '原文 13 个补充表中未出现'],
  ['两个区室改变的直接对比', '未见任何已发表研究做过'],
  ['MDA5 蛋白', '人类 MS 脑组织的免疫组化与蛋白质组测量为零篇'],
], [3.4, 8.7], 1.22, 0.70);

/* ------------------------------------------- 11 限制 */
s = plain();
title(s, '限制');
const lim = [
  ['皮层侧不是纯灰质', 'Schirmer 组织块含皮层下白质，且 MS 块的白质成分显著多于对照块，该侧的疾病对比与解剖深度混杂。'],
  ['交互项只有两类可估计', '小胶质、神经元、内皮因一侧供体不足，OPC 因对照中位为零。两类可估计者区间均触及零。'],
  ['白质结论主要由一个队列承担', 'GSE279180 十三位供体，3.5 倍，区间 [1.5, 4.9]；GSE180759 仅三位对照，1.5 倍，区间跨 1。'],
  ['检出率仍随深度变化', '已限定在每核 1500 至 5000 UMI 的共同区间内比较，但区间内仍有残余深度差异。'],
  ['GSE180759 存在平台混杂', '斑块周围文库全为 NovaSeq，病灶边缘文库全为 HiSeq；图谱已按文库做 Harmony 整合。'],
  ['不用表达量而用检出率', 'Schirmer 仅有对数归一化矩阵，尺度因子无法反解，跨队列的表达量不可比；检出与否不受归一化影响。'],
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
