import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
from scipy import stats
OUT='/home/user/aaa/results/mda5/figs/'

# ================= 图1 基线 =================
cons=pd.read_csv('/home/user/aaa/results/mda5/out_baseline_consolidated.csv')
cons=cons.rename(columns={cons.columns[0]:'grp'})
NAME={'endothelial':'内皮细胞','microglia/CNS macrophage':'小胶质/巨噬','fibroblast':'成纤维',
      'mural/perivascular':'周细胞/血管壁','lymphocyte/leukocyte':'淋巴/白细胞',
      'ependymal/choroid':'室管膜/脉络丛','astrocyte':'星形胶质','OPC':'OPC',
      'oligodendrocyte':'少突胶质','neuron':'神经元'}
cons['cn']=cons.grp.map(NAME); cons=cons.dropna(subset=['cn']).sort_values('CP10k',ascending=True)
FIVE=['内皮细胞','小胶质/巨噬','星形胶质','少突胶质','神经元']
bulk={'内皮细胞':6.5,'小胶质/巨噬':5.13,'星形胶质':1.70,'少突胶质':1.18,'神经元':0.8}

fig=plt.figure(figsize=(13.2,4.5))
gs=fig.add_gridspec(1,3,width_ratios=[1.55,1,1],wspace=0.52)
ax=fig.add_subplot(gs[0])
cols=[BLUE if g not in ('神经元',) else GREY for g in cons.cn]
ax.barh(cons.cn,cons.CP10k,color=cols,height=0.68)
for y,(v,n) in enumerate(zip(cons.CP10k,cons.n)):
    ax.text(v+0.022,y,f'{v:.3f}',va='center',fontsize=9.5,color=MUT)
ax.set_xlabel('IFIH1 表达量  CP10k（每万计数）'); ax.set_xlim(0,1.32)
ax.set_title('A  单核转录组｜CELLxGENE Census\n1 851 万个细胞 · 783 名供体 · 144 个数据集',fontsize=11.5,loc='left',pad=10)

ax2=fig.add_subplot(gs[1])
h=cons.dropna(subset=['nCPM']).sort_values('nCPM')
ax2.barh(h.cn,h.nCPM,color=[BLUE if g!='神经元' else GREY for g in h.cn],height=0.68)
for y,v in enumerate(h.nCPM): ax2.text(v+2.3,y,f'{v:.1f}',va='center',fontsize=9.5,color=MUT)
ax2.set_xlabel('IFIH1  nCPM'); ax2.set_xlim(0,143)
ax2.set_title('B  单核转录组｜人类蛋白图谱\n（Siletti 人脑图谱）',fontsize=11.5,loc='left',pad=10)

ax3=fig.add_subplot(gs[2])
bs=pd.Series(bulk).sort_values()
ax3.barh(bs.index,bs.values,color=[BLUE if g!='神经元' else GREY for g in bs.index],height=0.6)
for y,v in enumerate(bs.values): ax3.text(v+0.13,y,f'{v:.2f}',va='center',fontsize=9.5,color=MUT)
ax3.set_xlabel('IFIH1  FPKM'); ax3.set_xlim(0,8.2)
ax3.set_title('C  免疫纯化细胞 bulk 测序\n（Zhang 2016，GSE73721，非单核技术）',fontsize=11.5,loc='left',pad=10)
fig.text(0.008,-0.045,'三种彼此独立的技术给出同一次序：内皮最高，神经元最低。灰色为神经元。',fontsize=10.5,color=MUT)
plt.savefig(OUT+'fig1_baseline.png'); plt.close()
print('fig1 ok')

# ================= 图2 MS vs 对照 =================
A=pd.read_parquet('df_180759.parquet')
A['cp10k']=1e4*A.IFIH1/A.total
STAGES=['control_white_matter','MS_periplaque_white_matter','chronic_inactive_MS_lesion_edge',
        'chronic_active_MS_lesion_edge','MS_lesion_core']
CELLS=['immune','oligodendrocytes','astrocytes','neurons']
fig=plt.figure(figsize=(13.2,5.6))
gs=fig.add_gridspec(1,2,width_ratios=[2.25,1],wspace=0.3)
axA=fig.add_subplot(gs[0])
pos=0; xt=[]; xl=[]; sep=[]
for ci,c in enumerate(CELLS):
    for s in STAGES:
        v=np.log1p(A[(A.ct==c)&(A.path==s)].cp10k.values)
        rng=np.random.default_rng(int(pos*10))
        if len(v)>3:
            p=axA.violinplot([v],positions=[pos],widths=0.78,showextrema=False)
            for b in p['bodies']: b.set_facecolor(STAGE[s]); b.set_edgecolor('none'); b.set_alpha(0.38)
        nz=v[v>0]
        if len(nz):
            ss=nz if len(nz)<=220 else rng.choice(nz,220,replace=False)
            axA.scatter(pos+rng.uniform(-0.19,0.19,len(ss)),ss,s=3.6,color=STAGE[s],alpha=0.8,linewidths=0)
        if len(v): axA.hlines(np.mean(v),pos-0.30,pos+0.30,color='#333333',lw=1.6,zorder=6)
        axA.text(pos,-0.30,f'{len(v)}',ha='center',fontsize=7.4,color=MUT,rotation=90)
        pos+=1
    xt.append(pos-3); xl.append(CT_CN[c]); pos+=0.9
    if ci<len(CELLS)-1: sep.append(pos-1.45)
for x in sep: axA.axvline(x,color='#e0e0e0',lw=1)
axA.set_xticks(xt); axA.set_xticklabels(xl,fontsize=12)
axA.set_ylabel('IFIH1 表达量  log1p(CP10k)'); axA.set_ylim(-0.42,3.6)
hs=[plt.Line2D([],[],marker='s',ls='',ms=8,color=STAGE[s],label=STAGE_CN[s]) for s in STAGES]
axA.legend(handles=hs,loc='upper left',ncol=3,fontsize=9.5,handletextpad=0.4,columnspacing=1.1)
axA.set_title('A  GSE180759（Absinta 2021）｜按病理分区的单核表达分布',fontsize=11.5,loc='left',pad=10)

axB=fig.add_subplot(gs[1])
m=pd.read_csv('/home/user/aaa/results/mda5/out_meta_nb_glm.csv')
CN={'microglia':'小胶质/巨噬','oligodendrocyte':'少突胶质','astrocyte':'星形胶质','OPC':'OPC','endothelial':'内皮/血管'}
m['cn']=m.cell_type.map(CN); m=m.sort_values('fold_change')
y=np.arange(len(m))
sig=m.p<0.05
axB.hlines(y,m.ci_low,m.ci_high,color=[BLUE if s else GREY for s in sig],lw=2.2)
axB.scatter(m.fold_change,y,s=64,color=[BLUE if s else GREY for s in sig],zorder=5)
axB.axvline(1,color='#b0b0b0',ls='--',lw=1.1)
axB.set_yticks(y); axB.set_yticklabels(m.cn)
for i,r in enumerate(m.itertuples()):
    axB.text(5.1,i,f'{r.fold_change:.2f}×',va='center',fontsize=10,ha='left',
             color=INK if r.p<0.05 else MUT,fontweight='bold' if r.p<0.05 else 'normal')
    axB.text(9.8,i,f'p={r.p:.3f}',va='center',fontsize=9.3,ha='left',
             color=INK if r.p<0.05 else MUT)
axB.set_xscale('log'); axB.set_xlim(0.3,20)
axB.set_xlabel('MS / 对照  倍数变化（95% 置信区间）')
import matplotlib.ticker as mt
axB.xaxis.set_major_locator(mt.FixedLocator([0.5,1,2,4]))
axB.xaxis.set_major_formatter(mt.FixedFormatter(['0.5','1','2','4']))
axB.xaxis.set_minor_locator(mt.NullLocator())
axB.spines['bottom'].set_bounds(0.3,4.4)
for sp in ['right','top']: axB.spines[sp].set_visible(False)
axB.set_title('B  21 个脑的供体水平合并分析\n负二项模型，深度归一化',fontsize=11.5,loc='left',pad=10)
fig.text(0.008,-0.035,'A：每根小提琴下方数字为该组核数；深色短横为该组均值。B：虚线为无变化（1.0×），蓝色 p<0.05，灰色不显著。',fontsize=10.5,color=MUT)
plt.savefig(OUT+'fig2_ms_vs_ctrl.png'); plt.close()
print('fig2 ok')
