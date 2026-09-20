/* Nestin（NES）单页。数据来源放在右下角，是可编辑文本框，不烧进图片。 */
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const FIG = '/home/user/aaa/results/mda5/figs/';
const BLACK = '000000', GREY = '6E6E6E';
const CN = 'Microsoft YaHei';

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';                  // 13.333 x 7.5 英寸
pres.title = 'Nestin（NES）在神经系统各类细胞中的表达';

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

const s = pres.addSlide();
s.addText('Nestin（NES）在神经系统各类细胞中的表达',
  { x: 0.6, y: 0.30, w: 12.1, h: 0.52, fontSize: 25, bold: true, color: BLACK,
    fontFace: CN, isTextBox: true, margin: 0 });

/* 图按可用区等比铺满，底部给来源留出空间 */
const name = 'fig_nes_celltype.png';
const d = sizeOf(FIG + name);
const a = d.w / d.h;
const bx = 0.45, by = 1.05, bw = 12.45, bh = 7.5 - by - 0.95;
let w = bw, h = w / a;
if (h > bh) { h = bh; w = h * a; }
s.addImage({ path: FIG + name, x: bx + (bw - w) / 2, y: by + (bh - h) / 2, w, h });

/* 数据来源：右下角，可编辑 */
s.addText(
  [{ text: 'CELLxGENE Census · CZI Cell Science Program, et al. Nucleic Acids Res. 2025;53(D1):D886-D900',
     options: { breakLine: true } },
   { text: 'GSE180759 · Absinta M, et al. Nature. 2021;597(7878):709-714' }],
  { x: 4.6, y: 6.74, w: 8.15, h: 0.52, fontSize: 9, color: GREY, fontFace: CN,
    align: 'right', isTextBox: true, margin: 0, lineSpacing: 13 });

pres.writeFile({ fileName: '/home/user/aaa/results/mda5/Nestin_NES_细胞类型表达.pptx' })
  .then(f => console.log('written', f));
