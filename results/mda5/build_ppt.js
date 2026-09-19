const pptxgen = require('pptxgenjs');
const fs=require('fs');
const FIG='/home/user/aaa/results/mda5/figs/';
const BLACK='000000', GREY='6E6E6E', RULE='C8CCD2';
const CN='Microsoft YaHei', EN='Calibri';
const pres=new pptxgen();
pres.layout='LAYOUT_WIDE';
pres.title='MDA5（IFIH1）在人脑与多发性硬化中的表达';

const DIMS={};
const sizeOf=(p)=>{const b=fs.readFileSync(p);let i=8;while(i<b.length){const len=b.readUInt32BE(i);const t=b.toString('ascii',i+4,i+8);if(t==='IHDR')return{w:b.readUInt32BE(i+8),h:b.readUInt32BE(i+12)};i+=len+12;}return null;};
['figA_2x2.png','fig0_umap_atlas.png','fig1_baseline.png','fig3_gm_vs_wm.png',
 'fig4_lesion_stage.png','fig5_controls.png','fig6_umap_180759.png'].forEach(n=>{
  try{const s=sizeOf(FIG+n); if(s) DIMS[n]=s.w/s.h;}catch(e){}
});

let pageNo=0;
function plain(){const s=pres.addSlide();pageNo++;const n=pageNo;
  s.addText(String(n),{x:12.5,y:6.95,w:0.5,h:0.3,fontSize:10,color:GREY,align:'right',fontFace:EN,isTextBox:true,margin:0});
  return s;}
function title(s,t){s.addText(t,{x:0.6,y:0.34,w:12.1,h:0.62,fontSize:26,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});}
function figure(s,name,top){
  const a=DIMS[name]; if(!a) return;
  const bx=0.45,by=top,bw=12.45,bh=7.5-top-0.42;
  let w=bw,h=w/a; if(h>bh){h=bh;w=h*a;}
  s.addImage({path:FIG+name,x:bx+(bw-w)/2,y:by+(bh-h)/2,w,h});
}

/* 1 */
let s=pres.addSlide();
s.addText('MDA5（IFIH1）在人脑与多发性硬化中的表达',{x:0.9,y:2.7,w:11.5,h:0.8,fontSize:34,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
s.addText('正常与 MS ｜ 灰质与白质',{x:0.9,y:3.62,w:11.5,h:0.5,fontSize:20,color:GREY,fontFace:CN,isTextBox:true,margin:0});

/* 2 问题与数据 */
s=plain(); title(s,'两个问题');
s.addText('一　正常人脑中 MDA5 由哪些细胞表达？灰质与白质是否不同？',
  {x:0.75,y:1.22,w:12,h:0.44,fontSize:18,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
s.addText('二　MS 中这一分布如何改变？灰质与白质的改变是否相同？',
  {x:0.75,y:1.80,w:12,h:0.44,fontSize:18,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
const rows=[
 [{text:'队列',options:{bold:true}},{text:'区室',options:{bold:true}},{text:'供体',options:{bold:true}},{text:'核数',options:{bold:true}},{text:'承担',options:{bold:true}}],
 ['CELLxGENE Census','正常脑（灰质区域）','783','18 551 076','细胞类型次序'],
 ['GSE180759  Absinta 2021','皮层下白质','5 MS + 3 对照','66 432','白质 · 病灶分期'],
 ['GSE279180  Lerma-Martin 2024','皮层下白质','7 MS + 6 对照','103 794','白质 · 主力'],
 ['GSE118257  Jäkel 2019','白质','4 MS + 5 对照','17 799','白质 · 验证'],
 ['Schirmer 2019','皮层灰质','12 MS + 9 对照','48 919','灰质'],
];
s.addTable(rows,{x:0.6,y:2.72,w:12.1,colW:[3.6,3.0,2.1,1.9,1.5],fontSize:12.5,fontFace:CN,color:BLACK,
  border:[{type:'none'},{type:'none'},{type:'solid',color:RULE,pt:1},{type:'none'}],
  rowH:0.46,valign:'middle',autoPage:false});

/* 3 问题一 细胞类型 */
s=plain(); title(s,'问题一：哪些细胞表达 MDA5'); figure(s,'fig1_baseline.png',1.22);

/* 4 核心：区域对比 */
s=plain(); title(s,'问题一与二：灰质与白质，先看正常人，再看病人'); figure(s,'figA_2x2.png',1.22);

/* 5 疾病效应按区室 */
s=plain(); title(s,'同一批数据的另一种读法：MS 相对对照，分区室'); figure(s,'fig3_gm_vs_wm.png',1.22);

/* 6 UMAP */
s=plain(); title(s,'白质中 IFIH1 的细胞分布'); figure(s,'fig0_umap_atlas.png',1.22);

/* 7 白质内部 */
s=plain(); title(s,'白质内部：病灶分期与小胶质状态'); figure(s,'fig4_lesion_stage.png',1.22);

/* 8 对照 */
s=plain(); title(s,'对照'); figure(s,'fig5_controls.png',1.22);

/* 9 文献 */
s=plain(); title(s,'与已发表结果的关系');
const r2=[
 [{text:'结果',options:{bold:true}},{text:'文献状态',options:{bold:true}}],
 ['少突胶质在白质病灶上调','Lerma-Martin 2024 补充表 6 已报道：CA log2FC +0.905（padj 0.012），CI +1.095（padj 0.005）'],
 ['同一队列髓系未达显著','与原作者一致，其补充表髓系 padj = 0.116'],
 ['bulk 白质病灶同向上调','GSE138614 活动病灶 FDR 0.021、慢性活动 FDR 0.007；GSE283092 泡沫样活动病灶 1.79 倍'],
 ['皮层灰质为阴性','Schirmer 2019 补充表中 IFIH1 低于其检出阈值，本次为原创测量'],
 ['GSE180759 中的 IFIH1','原文 13 个补充表中未出现'],
 ['灰质与白质的直接对比','未见任何已发表研究做过'],
 ['MDA5 蛋白','人类 MS 脑组织的免疫组化与蛋白质组测量为零篇'],
];
s.addTable(r2,{x:0.6,y:1.30,w:12.1,colW:[3.6,8.5],fontSize:12.5,fontFace:CN,color:BLACK,
  border:[{type:'none'},{type:'none'},{type:'solid',color:RULE,pt:1},{type:'none'}],
  rowH:0.48,valign:'middle',autoPage:false});

/* 10 限制 */
s=plain(); title(s,'限制');
const lim=[
 ['灰质与白质来自不同研究','无任何公开队列同时含两者；Macnair 2025 有，但在 EGA 受控访问'],
 ['正常灰质的小胶质仅 2 个供体','因此“正常时两区室无差别”这一条在小胶质上无法判定'],
 ['两个白质队列存在系统偏移','GSE180759 的绝对值约为 GSE279180 的一半，方向一致'],
 ['测序深度与疾病状态混杂','GSE279180 的 MS 小胶质深度为对照的 3.84 倍，已全程深度归一化'],
 ['GSE180759 平台完全混杂','4 个斑块周围文库全为 NovaSeq，11 个病灶边缘文库全为 HiSeq'],
 ['全部为转录本层面','人类 MS 脑组织无任何 MDA5 蛋白测量'],
];
let yy=1.42;
lim.forEach(L=>{
  s.addText(L[0],{x:0.6,y:yy,w:4.3,h:0.42,fontSize:14,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
  s.addText(L[1],{x:5.0,y:yy,w:7.7,h:0.6,fontSize:12.5,color:GREY,fontFace:CN,isTextBox:true,margin:0,lineSpacing:17});
  yy+=0.88;
});

pres.writeFile({fileName:'/home/user/aaa/results/mda5/MDA5_IFIH1_人脑与MS.pptx'}).then(f=>console.log('written',f));
