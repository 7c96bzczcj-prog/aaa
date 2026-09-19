const pptxgen = require('pptxgenjs');
const FIG='/home/user/aaa/results/mda5/figs/';
const BLACK='000000', GREY='6E6E6E', RULE='C8CCD2';
const CN='Microsoft YaHei', EN='Calibri';
const pres=new pptxgen();
pres.layout='LAYOUT_WIDE';
pres.title='MDA5（IFIH1）在人脑与多发性硬化中的表达';

const DIMS={};               // filled from disk
const fs=require('fs');
const sizeOf=(p)=>{const b=fs.readFileSync(p);let i=8;while(i<b.length){const len=b.readUInt32BE(i);const type=b.toString('ascii',i+4,i+8);if(type==='IHDR')return{w:b.readUInt32BE(i+8),h:b.readUInt32BE(i+12)};i+=len+12;}return null;};
['fig0_umap_atlas.png','fig1_baseline.png','fig2_ms_vs_ctrl.png','fig3_gm_vs_wm.png',
 'fig4_lesion_stage.png','fig5_controls.png','fig6_umap_180759.png'].forEach(n=>{
  try{const s=sizeOf(FIG+n); if(s) DIMS[n]=s.w/s.h;}catch(e){}
});

let pageNo=0;
function plain(){ const s=pres.addSlide(); pageNo++; const n=pageNo;
  s.addText(String(n),{x:12.5,y:6.95,w:0.5,h:0.3,fontSize:10,color:GREY,align:'right',fontFace:EN,isTextBox:true,margin:0});
  return s; }
function title(s,t){ s.addText(t,{x:0.6,y:0.36,w:12.1,h:0.62,fontSize:27,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0}); }
function figure(s,name,top){
  const a=DIMS[name]; if(!a) return;
  const boxX=0.45,boxY=top,boxW=12.45,boxH=7.5-top-0.42;
  let w=boxW,h=w/a; if(h>boxH){h=boxH;w=h*a;}
  s.addImage({path:FIG+name,x:boxX+(boxW-w)/2,y:boxY+(boxH-h)/2,w,h});
}

/* 1 title */
let s=pres.addSlide();
s.addText('MDA5（IFIH1）在人脑与多发性硬化中的表达',{x:0.9,y:2.75,w:11.5,h:0.8,fontSize:36,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
s.addText('CELLxGENE Census · GSE180759 · GSE279180 · GSE118257 · Schirmer 2019',
  {x:0.9,y:3.75,w:11.5,h:0.4,fontSize:14,color:GREY,fontFace:CN,isTextBox:true,margin:0});

/* 2 datasets */
s=plain(); title(s,'数据');
const rows=[
 [{text:'队列',options:{bold:true}},{text:'组织',options:{bold:true}},{text:'供体',options:{bold:true}},{text:'核数',options:{bold:true}},{text:'用途',options:{bold:true}}],
 ['CELLxGENE Census','正常人脑','783','18 551 076','基线'],
 ['GSE180759  Absinta 2021','皮层下白质，MRI 分期','5 MS + 3 对照','66 432','病灶分期'],
 ['GSE279180  Lerma-Martin 2024','皮层下白质','7 MS + 6 对照','103 794','主力队列'],
 ['GSE118257  Jäkel 2019','白质','4 MS + 5 对照','17 799','方向验证'],
 ['Schirmer 2019','皮层灰质','12 MS + 9 对照','48 919','灰质与神经元'],
];
s.addTable(rows,{x:0.6,y:1.35,w:12.1,colW:[3.5,3.0,2.1,1.9,1.6],fontSize:13,fontFace:CN,color:BLACK,
  border:[{type:'none'},{type:'none'},{type:'solid',color:RULE,pt:1},{type:'none'}],
  rowH:0.5,valign:'middle',autoPage:false});
s.addText('灰质与白质分属不同研究：Macnair 2025 为唯一同时含两者的队列，存于 EGA 受控访问。',
  {x:0.6,y:5.1,w:12.1,h:0.4,fontSize:13,color:GREY,fontFace:CN,isTextBox:true,margin:0});

/* 3 baseline */
s=plain(); title(s,'基线：MDA5 在正常人脑各细胞类型中的分布'); figure(s,'fig1_baseline.png',1.25);

/* 4 UMAP 180759 (added later if present) */
if(DIMS['fig6_umap_180759.png']){ s=plain(); title(s,'GSE180759：细胞类型与 IFIH1 按病理分区的分布'); figure(s,'fig6_umap_180759.png',1.25); }

/* 5 UMAP atlas */
s=plain(); title(s,'GSE279180：细胞类型与 IFIH1 按病灶类型的分布'); figure(s,'fig0_umap_atlas.png',1.25);

/* 6 MS vs ctrl */
s=plain(); title(s,'MS 与正常对照的差异'); figure(s,'fig2_ms_vs_ctrl.png',1.25);

/* 7 GM vs WM */
s=plain(); title(s,'灰质与白质'); figure(s,'fig3_gm_vs_wm.png',1.25);

/* 8 lesion stage */
s=plain(); title(s,'病灶分期与小胶质状态'); figure(s,'fig4_lesion_stage.png',1.25);

/* 9 controls */
s=plain(); title(s,'对照'); figure(s,'fig5_controls.png',1.25);

/* 10 prior art */
s=plain(); title(s,'与已发表结果的关系');
const r2=[
 [{text:'结果',options:{bold:true}},{text:'文献状态',options:{bold:true}}],
 ['少突胶质在白质病灶上调','Lerma-Martin 2024 补充表 6 已报道：CA log2FC +0.905（padj 0.012），CI +1.095（padj 0.005）'],
 ['同一队列髓系未达显著','与原作者一致，其补充表髓系 padj = 0.116'],
 ['bulk 白质病灶同向上调','GSE138614 活动病灶 FDR 0.021、慢性活动 FDR 0.007；GSE283092 泡沫样活动病灶 1.79 倍'],
 ['RIG-I 不随 MDA5 上升','GSE138614 中 DDX58 与 DHX58 在任何病灶类型均未显著'],
 ['GSE180759 中的 IFIH1','原文 13 个补充表中未出现'],
 ['皮层灰质为阴性','Schirmer 2019 补充表中 IFIH1 低于其检出阈值'],
 ['MDA5 蛋白','人类 MS 脑组织的免疫组化与蛋白质组测量为零篇'],
];
s.addTable(r2,{x:0.6,y:1.35,w:12.1,colW:[3.6,8.5],fontSize:12.5,fontFace:CN,color:BLACK,
  border:[{type:'none'},{type:'none'},{type:'solid',color:RULE,pt:1},{type:'none'}],
  rowH:0.48,valign:'middle',autoPage:false});

/* 11 limits */
s=plain(); title(s,'限制');
const lim=[
 ['测序深度与疾病状态混杂','GSE279180 的 MS 小胶质深度为对照的 3.84 倍；配平后效应量由 4.13 降至 1.64'],
 ['GSE180759 平台完全混杂','4 个斑块周围文库全为 NovaSeq，11 个病灶边缘文库全为 HiSeq'],
 ['对照脑数量过少','GSE180759 仅 3 个对照脑，可用于髓系检验者 2 个，p 值下限 0.095'],
 ['灰白质比较跨研究','灰质来自 Schirmer，白质来自另外三项研究'],
 ['全部为转录本层面','人类 MS 脑组织无任何 MDA5 蛋白测量'],
 ['环境 RNA 未排除','GSE180759 中小胶质与少突胶质的 IFIH1 检出率跨文库相关，Spearman ρ = 0.63'],
];
let yy=1.45;
lim.forEach(L=>{
  s.addText(L[0],{x:0.6,y:yy,w:3.8,h:0.4,fontSize:14,bold:true,color:BLACK,fontFace:CN,isTextBox:true,margin:0});
  s.addText(L[1],{x:4.5,y:yy,w:8.2,h:0.55,fontSize:12.5,color:GREY,fontFace:CN,isTextBox:true,margin:0,lineSpacing:17});
  yy+=0.87;
});

pres.writeFile({fileName:'/home/user/aaa/results/mda5/MDA5_IFIH1_人脑与MS.pptx'}).then(f=>console.log('written',f));
