/* English edition, mirroring the nine-slide Chinese deck. Figures come from figs_en/. */
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const FIG = '/home/user/aaa/results/mda5/figs_en/';
const BLACK = '000000', GREY = '6E6E6E', RULE = 'C8CCD2';
const EN = 'Calibri';

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';                  // 13.333 x 7.5 in
pres.title = 'MDA5 (IFIH1) and Nestin (NES) in human brain and multiple sclerosis';

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
['fig_q1_celltype.png', 'fig_nes_celltype.png', 'fig_q1_region.png', 'fig_q2_wm.png',
 'fig_q2_gm.png', 'fig_umap_wm.png', 'fig_lesion.png'].forEach(n => {
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
function title(s, t) {
  s.addText(t, { x: 0.6, y: 0.30, w: 12.1, h: 0.56, fontSize: 24, bold: true,
    color: BLACK, fontFace: EN, isTextBox: true, margin: 0 });
}
function figure(s, name, top, bottomPad = 0.42) {
  const a = DIMS[name];
  if (!a) { console.warn('missing figure', name); return; }
  const bx = 0.45, by = top, bw = 12.45, bh = 7.5 - top - bottomPad;
  let w = bw, h = w / a;
  if (h > bh) { h = bh; w = h * a; }
  s.addImage({ path: FIG + name, x: bx + (bw - w) / 2, y: by + (bh - h) / 2, w, h });
}
const B = t => ({ text: t, options: { bold: true } });

/* --------------------------------------------------------------- 1 title */
let s = pres.addSlide();
s.addText('MDA5 (IFIH1) and Nestin (NES) expression in the human brain and in multiple sclerosis',
  { x: 0.9, y: 2.60, w: 11.5, h: 1.0, fontSize: 30, bold: true, color: BLACK,
    fontFace: EN, isTextBox: true, margin: 0 });
s.addText('Cell-type-resolved single-nucleus transcriptomic analysis',
  { x: 0.9, y: 3.72, w: 11.5, h: 0.5, fontSize: 17, color: GREY,
    fontFace: EN, isTextBox: true, margin: 0 });

/* ------------------------------------------------- 2 questions and data */
s = plain();
title(s, 'Two questions');
s.addText('1   Which cell types express IFIH1 in the normal human brain, and does this differ between grey and white matter?',
  { x: 0.75, y: 0.98, w: 12, h: 0.70, fontSize: 15, bold: true, color: BLACK,
    fontFace: EN, isTextBox: true, margin: 0, lineSpacing: 21 });
s.addText('2   How does the cellular distribution of IFIH1 change in MS lesions, and is the change the same in the two regions?',
  { x: 0.75, y: 1.70, w: 12, h: 0.70, fontSize: 15, bold: true, color: BLACK,
    fontFace: EN, isTextBox: true, margin: 0, lineSpacing: 21 });
s.addTable([
  [B('Cohort'), B('Compartment'), B('Donors'), B('Nuclei'), B('Role')],
  ['CELLxGENE Census', 'Normal brain', '783', '18,551,076', 'Cell-type expression profile'],
  ['GSE180759  Absinta 2021', 'Subcortical white matter', '3 ctrl + 5 MS', '66,432', 'White matter, lesion stages'],
  ['GSE279180  Lerma-Martin 2024', 'Subcortical white matter', '6 ctrl + 7 MS', '103,780', 'White matter, primary cohort'],
  ['Schirmer 2019', 'Cortical block', '9 ctrl + 12 MS', '48,919', 'Cortex'],
  ['Macnair 2025', 'Cortex + subcortical WM', '26 ctrl + 54 MS', '632,375', 'Lesion stages, within one study'],
], { x: 0.6, y: 2.50, w: 12.1, colW: [3.0, 3.1, 2.1, 1.7, 2.2], fontSize: 12, fontFace: EN,
     color: BLACK, rowH: 0.56, valign: 'middle', autoPage: false,
     border: [{ type: 'none' }, { type: 'none' }, { type: 'solid', color: RULE, pt: 1 }, { type: 'none' }] });
s.addText(
  [{ text: 'Nucleic Acids Res. 2025 Jan 6;53(D1):D886-D900;', options: { breakLine: true } },
   { text: 'Nature. 2021 Sep;597(7878):709-714;', options: { breakLine: true } },
   { text: 'Nat Neurosci. 2024 Dec;27(12):2354-2365;', options: { breakLine: true } },
   { text: 'Nature. 2019 Sep;573(7772):75-82;', options: { breakLine: true } },
   { text: 'Neuron. 2025 Feb 5;113(3):396-410.e9.' }],
  { x: 0.6, y: 5.95, w: 12.1, h: 1.10, fontSize: 10, color: GREY, fontFace: EN,
    isTextBox: true, margin: 0, lineSpacing: 14 });

/* --------------------------------------------------- 3 IFIH1 cell types */
s = plain();
title(s, 'IFIH1 expression across human brain cell types');
figure(s, 'fig_q1_celltype.png', 1.00);

/* ----------------------------------------------------- 4 Nestin */
s = plain();
title(s, 'Nestin (NES) expression across CNS cell types');
figure(s, 'fig_nes_celltype.png', 1.00, 0.95);
s.addText(
  [{ text: 'CELLxGENE Census · CZI Cell Science Program, et al. Nucleic Acids Res. 2025;53(D1):D886-D900',
     options: { breakLine: true } },
   { text: 'GSE180759 · Absinta M, et al. Nature. 2021;597(7878):709-714' }],
  { x: 4.4, y: 6.68, w: 7.95, h: 0.52, fontSize: 9, color: GREY, fontFace: EN,
    align: 'right', isTextBox: true, margin: 0, lineSpacing: 13 });

/* ------------------------------------------------- 5 normal, two regions */
s = plain();
title(s, 'Normal brain: no difference between the two anatomical regions');
figure(s, 'fig_q1_region.png', 1.00);

/* ------------------------------------------------------ 6 MS white matter */
s = plain();
title(s, 'MS white matter: more IFIH1-positive oligodendrocytes');
figure(s, 'fig_q2_wm.png', 1.00);

/* ------------------------------------------------------------ 7 MS cortex */
s = plain();
title(s, 'MS cortex: no change of comparable magnitude');
figure(s, 'fig_q2_gm.png', 1.00);

/* ----------------------------------------------------------- 8 atlas */
s = plain();
title(s, 'IFIH1 across white-matter nuclei, by lesion stage');
figure(s, 'fig_umap_wm.png', 1.00);

/* ------------------------------------------------- 9 lesion stage testing */
s = plain();
title(s, 'Each demyelinating lesion stage against its own control');
figure(s, 'fig_lesion.png', 1.00);

pres.writeFile({ fileName: '/home/user/aaa/results/mda5/MDA5_IFIH1_NES_human_brain_MS_EN.pptx' })
  .then(f => console.log('written', f, '|', pageNo + 1, 'slides'));
