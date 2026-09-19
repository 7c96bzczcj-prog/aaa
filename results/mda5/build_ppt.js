const pptxgen = require('pptxgenjs');
const path = require('path');
const FIG = '/home/user/aaa/results/mda5/figs/';
const NAVY='16233F', BLUE='2A78D6', ORANGE='EB6834', INK='1A1A1A', MUT='5A5A5A', LINE='E4E6EA', WHITE='FFFFFF';
const CN='Microsoft YaHei', EN='Calibri';
const pres = new pptxgen();
pres.layout='LAYOUT_WIDE';            // 13.33 x 7.5
pres.author='MDA5 analysis'; pres.title='MDA5（IFIH1）在人脑与多发性硬化中的表达';

const DIMS={'fig1_baseline.png':2.654,'fig2_ms_vs_ctrl.png':2.065,'fig3_gm_vs_wm.png':1.993,
            'fig4_lesion_stage.png':2.49,'fig5_controls.png':2.401};
function fit(name, boxX, boxY, boxW, boxH){
  const a=DIMS[name]; let w=boxW, h=w/a;
  if(h>boxH){ h=boxH; w=h*a; }
  return {path:FIG+name, x:boxX+(boxW-w)/2, y:boxY+(boxH-h)/2, w:w, h:h};
}
function head(s, num, title, reading){
  s.addText(String(num).padStart(2,'0'), {x:0.55,y:0.30,w:0.7,h:0.42,fontSize:15,bold:true,color:BLUE,fontFace:EN,isTextBox:true,margin:0});
  s.addText(title, {x:1.15,y:0.26,w:11.7,h:0.5,fontSize:25,bold:true,color:INK,fontFace:CN,isTextBox:true,margin:0});
  if(reading) s.addText(reading, {x:1.15,y:0.82,w:11.7,h:0.42,fontSize:13.5,color:MUT,fontFace:CN,isTextBox:true,margin:0});
}
function src(s, t){
  s.addText(t, {x:0.55,y:6.98,w:12.3,h:0.32,fontSize:10,color:'8A8F98',fontFace:CN,isTextBox:true,margin:0});
}

/* ---------- 1 封面 ---------- */
let s=pres.addSlide(); s.background={color:NAVY};
s.addText('MDA5（IFIH1）', {x:0.9,y:2.05,w:11.5,h:0.85,fontSize:44,bold:true,color:WHITE,fontFace:CN,isTextBox:true,margin:0});
s.addText('在人脑与多发性硬化中的表达', {x:0.9,y:2.92,w:11.5,h:0.8,fontSize:36,bold:true,color:'CADCFC',fontFace:CN,isTextBox:true,margin:0});
s.addText('四个人类单核转录组队列 · 21 个脑 · 全部深度归一化后重新分析', {x:0.9,y:4.0,w:11.5,h:0.42,fontSize:16,color:'AFC0DE',fontFace:CN,isTextBox:true,margin:0});
s.addText('GSE180759 · GSE279180 · GSE118257 · Schirmer 2019 · CELLxGENE Census（1 851 万个细胞）',
  {x:0.9,y:6.25,w:11.5,h:0.36,fontSize:12,color:'8FA3C8',fontFace:CN,isTextBox:true,margin:0});
s.addNotes('本页为封面。全部数据均由原始数据重新计算，未直接引用原文结论。');

/* ---------- 2 数据集 ---------- */
s=pres.addSlide();
head(s,1,'用了哪些数据，各能回答什么','每个队列只在它能回答的问题上被使用；不能回答的，明确标注为空白。');
const rows=[
 [{text:'队列',options:{bold:true}},{text:'组织与分区',options:{bold:true}},{text:'供体',options:{bold:true}},{text:'在本分析中承担的问题',options:{bold:true}}],
 ['CELLxGENE Census','正常人脑，144 个数据集','783','基线：各细胞类型的表达次序'],
 ['GSE180759  Absinta 2021','皮层下白质，MRI 分期病灶','5 MS + 3 对照','白质的病灶分期梯度'],
 ['GSE279180  Lerma-Martin 2024','皮层下白质，慢性活动/非活动','7 MS + 6 对照','白质主力队列；小胶质亚型'],
 ['GSE118257  Jäkel 2019','白质','4 MS + 5 对照','方向一致性（功效不足）'],
 ['Schirmer 2019','皮层灰质（四个脑区）','12 MS + 9 对照','灰质；神经元'],
];
s.addTable(rows,{x:0.55,y:1.55,w:12.25,colW:[2.9,3.5,2.0,3.85],
  fontSize:12.5,fontFace:CN,color:INK,border:{type:'solid',color:LINE,pt:1},
  fill:{color:WHITE},rowH:0.46,valign:'middle',autoPage:false});
s.addShape(pres.ShapeType.roundRect,{x:0.55,y:4.78,w:12.25,h:1.25,fill:{color:'FFF4EE'},line:{color:'F7CDB8',pt:1},rectRadius:0.06});
s.addText([{text:'一个结构性限制：',options:{bold:true,color:'B4542C'}},
  {text:'没有任何一个可公开获取处理后数据的队列同时采了灰质和白质。为此设计的 Macnair 2025（63.2 万个核，54 MS + 26 对照，含灰质与白质病灶）存放于 EGA 受控访问，无法取用。因此灰质与白质的比较是跨研究的，这是本分析最主要的弱点。',options:{color:'6B4A38'}}],
  {x:0.85,y:4.95,w:11.65,h:0.95,fontSize:12.5,fontFace:CN,isTextBox:true,margin:0,lineSpacing:19});
src(s,'供体数为该队列在本分析中实际使用的数量。');
s.addNotes('强调：区室比较跨研究，是最主要的弱点。');

/* ---------- 3 基线 ---------- */
s=pres.addSlide();
head(s,2,'基线：MDA5 在正常人脑各细胞类型中的分布','三种彼此独立的技术给出同一次序，从最高到最低约 60 倍。');
s.addImage(fit('fig1_baseline.png',0.55,1.40,12.25,5.35));
src(s,'A：CELLxGENE Census 2025-11-08，单核子集。B：人类蛋白图谱单核脑图谱。C：Zhang 2016 免疫纯化细胞 bulk 测序（GSE73721），非单核技术。');
s.addNotes('C 是 bulk 纯化细胞，技术上与单核完全独立，用来排除单核测序的低检出偏倚。');

/* ---------- 4 MS vs 对照 ---------- */
s=pres.addSlide();
head(s,3,'MS 与正常对照：哪些细胞类型发生变化','小胶质/巨噬与少突胶质各约上升 2 倍；星形胶质、OPC、内皮不变。');
s.addImage(fit('fig2_ms_vs_ctrl.png',0.55,1.40,12.25,5.35));
src(s,'A：GSE180759 单核分布，纵轴为深度归一化后的 log1p(CP10k)。B：GSE180759 + GSE279180 共 21 个脑，负二项模型，供体为分析单位。');
s.addNotes('B 用供体作单位，避免同一供体内核之间的伪重复。');

/* ---------- 5 GM vs WM ---------- */
s=pres.addSlide();
head(s,4,'灰质与白质：变化只出现在白质','少突胶质在白质上升、在灰质不动；神经元只有灰质测得到，结果为阴性。');
s.addImage(fit('fig3_gm_vs_wm.png',0.55,1.40,12.25,5.30));
src(s,'A：同一深度配平统计量。灰质少突胶质每组 3 070 个核、神经元 10 134 个核，阴性不是功效不足；灰质小胶质仅 159 个核，无法判定。');
s.addNotes('灰质小胶质 159 个核，9 个对照供体中只有 2 个达到 20 个核，所以灰质不能验证小胶质结果。');

/* ---------- 6 分期与状态 ---------- */
s=pres.addSlide();
head(s,5,'白质内部：随病灶分期与小胶质状态递增','对照 0.118 → 病灶核心 0.386；慢性活动态与病灶边缘态小胶质最高。');
s.addImage(fit('fig4_lesion_stage.png',0.55,1.40,12.25,5.35));
src(s,'A：GSE180759 小胶质/巨噬，按 MRI 分区。B、C：GSE279180 小胶质亚型标注为原作者所给。');
s.addNotes('B 的慢性活动态与病灶边缘态正是 Absinta 2021 描述的驱动病灶扩张的那一群。');

/* ---------- 7 对照 ---------- */
s=pres.addSlide();
head(s,6,'三项必须先做的对照','深度混杂、深度校正的后果、以及 MDA5 相对 RIG-I 的选择性。');
s.addImage(fit('fig5_controls.png',0.55,1.40,12.25,5.35));
src(s,'A、B：GSE279180。C：三个彼此独立的设定。GSE138614 为 73 个 MS 白质区对 25 个对照白质区的 bulk 测序。');
s.addNotes('A 是本分析最重要的方法学发现：不做深度配平，小胶质效应会被高估 2.5 倍。');

/* ---------- 8 与已发表结果的关系 ---------- */
s=pres.addSlide();
head(s,7,'与已发表结果的关系','哪些是重现别人已报道的，哪些是这批数据里原本没有人写过的。');
const r2=[
 [{text:'本分析中的结果',options:{bold:true}},{text:'文献中的状态',options:{bold:true}}],
 ['少突胶质在白质病灶上调','已报道。Lerma-Martin 2024 补充表 6：慢性活动 log2FC +0.905（padj 0.012），慢性非活动 +1.095（padj 0.005）'],
 ['同一队列中髓系未达显著','与原作者一致。其补充表中髓系 padj = 0.116'],
 ['bulk 白质病灶同向上调','已报道。GSE138614 活动病灶 FDR 0.021、慢性活动 FDR 0.007；GSE283092 泡沫样活动病灶 1.79 倍'],
 ['RIG-I 不随 MDA5 上升','已报道。GSE138614 中 DDX58 与 DHX58 在任何病灶类型均未显著'],
 ['GSE180759 中的 IFIH1','原文 13 个补充表中从未出现，本次为该数据集的原创测量'],
 ['皮层灰质为阴性','Schirmer 2019 补充表中 IFIH1 低于其检出阈值，本次为原创测量'],
 ['MDA5 蛋白','人类 MS 脑组织的免疫组化、蛋白质组测量为零篇'],
];
s.addTable(r2,{x:0.55,y:1.58,w:12.25,colW:[3.7,8.55],fontSize:12,fontFace:CN,color:INK,
  border:{type:'solid',color:LINE,pt:1},fill:{color:WHITE},rowH:0.42,valign:'middle',autoPage:false});
src(s,'文献核查覆盖 8 个方向共 145 条论断，全部经对抗式复核：84 条确认、57 条判定为夸大、3 条被推翻。');
s.addNotes('少突胶质这一条不是新发现，是独立重现了作者自己的补充表。');

/* ---------- 9 限制 ---------- */
s=pres.addSlide(); s.background={color:NAVY};
s.addText('读这批数据时必须同时记住的六条限制', {x:0.9,y:0.72,w:11.6,h:0.62,fontSize:28,bold:true,color:WHITE,fontFace:CN,isTextBox:true,margin:0});
const lim=[
 ['测序深度与疾病状态混杂','GSE279180 的 MS 小胶质深度是对照的 3.84 倍。不配平，效应量从 1.64 被高估到 4.13。'],
 ['GSE180759 存在平台完全混杂','4 个斑块周围文库全部是 NovaSeq，11 个病灶边缘文库全部是 HiSeq，二者不可直接比较。'],
 ['对照脑数量过少','GSE180759 只有 3 个对照脑，其中可用于髓系检验的只有 2 个，其 p 值不可能低于 0.095。'],
 ['灰白质比较跨研究','灰质来自 Schirmer，白质来自另外三项研究，区室差异与研究批次无法分离。'],
 ['全部为转录本层面','人类 MS 脑组织没有任何 MDA5 蛋白测量；正常脑的抗体数据反而指向神经元而非胶质。'],
 ['环境 RNA 未排除','GSE180759 中小胶质与少突胶质的 IFIH1 检出率跨文库相关（Spearman ρ=0.63）。'],
];
let yy=1.62;
lim.forEach((L,i)=>{
  s.addShape(pres.ShapeType.ellipse,{x:0.95,y:yy+0.03,w:0.30,h:0.30,fill:{color:i<3?'E05A4B':'3D6FB8'}});
  s.addText(String(i+1),{x:0.95,y:yy+0.03,w:0.30,h:0.30,fontSize:12,bold:true,color:WHITE,align:'center',valign:'middle',fontFace:EN,isTextBox:true,margin:0});
  s.addText(L[0],{x:1.42,y:yy,w:3.5,h:0.34,fontSize:14,bold:true,color:'CADCFC',fontFace:CN,isTextBox:true,margin:0});
  s.addText(L[1],{x:5.0,y:yy-0.02,w:7.5,h:0.62,fontSize:12,color:'B6C4DC',fontFace:CN,isTextBox:true,margin:0,lineSpacing:17});
  yy+=0.82;
});
s.addNotes('前三条是技术性的，后三条是设计性的，都不能靠重新分析消除。');

pres.writeFile({fileName:'/home/user/aaa/results/mda5/MDA5_IFIH1_人脑与MS.pptx'}).then(f=>console.log('written',f));
