import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
OUT='/home/user/aaa/results/mda5/figs/'
CN={'microglia':'小胶质/巨噬','oligodendrocyte':'少突胶质','astrocyte':'星形胶质',
    'OPC':'OPC','endothelial':'内皮/血管','neuron':'神经元'}

# ================= 图3 灰质 vs 白质 =================
T=pd.read_csv('/home/user/aaa/results/mda5/out_gm_wm_depthmatched.csv')
T['cn']=T.cell_type.map(CN)
order=['少突胶质','小胶质/巨噬','星形胶质','神经元','OPC','内皮/血管']
fig=plt.figure(figsize=(13.4,5.4)); gs=fig.add_gridspec(1,2,width_ratios=[1.62,1],wspace=0.3)
ax=fig.add_subplot(gs[0])
SLOTS=[('white matter','GSE180759 Absinta','#7fb0e8'),
       ('white matter','GSE279180 Lerma-Martin',BLUE),
       ('cortex / grey matter','Schirmer 2019',ORANGE)]
h=0.24; LABX=4.45
for i,cn in enumerate(order):
    for k,(comp,coh,col) in enumerate(SLOTS):
        r=T[(T.cn==cn)&(T.compartment==comp)&(T.cohort==coh)]
        y=i+(k-1)*h
        if r.empty:
            ax.text(0.06,y,'该区室无足够细胞',va='center',fontsize=8.2,color='#c0c0c0'); continue
        r=r.iloc[0]
        c=col if r.p<0.05 else ('#c9dcf2' if comp=='white matter' else '#f7cdb8')
        ax.barh(y,min(r.OR,4.3),height=h*0.86,color=c)
        if r.OR>4.3: ax.text(4.33,y,'▸',va='center',fontsize=9,color=c)
        ptxt='p<0.005' if r.p<0.005 else f'p={r.p:.2f}'
        ax.text(LABX,y,f'{r.OR:.2f}',va='center',ha='right',fontsize=9.4,
                color=INK if r.p<0.05 else MUT,fontweight='bold' if r.p<0.05 else 'normal')
        ax.text(LABX+0.25,y,ptxt,va='center',ha='left',fontsize=8.7,color=INK if r.p<0.05 else MUT)
        ax.text(LABX+1.95,y,f'n={int(r.n_per_group)}',va='center',ha='right',fontsize=8.7,color=MUT)
ax.axvline(1,color='#b0b0b0',ls='--',lw=1.1)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=11.5)
ax.set_xlim(0,6.5); ax.set_xlabel('MS / 对照  检出率优势比（深度十分位配平后）')
ax.invert_yaxis()
import matplotlib.ticker as mt
ax.xaxis.set_major_locator(mt.FixedLocator([0,1,2,3,4]))
ax.spines['bottom'].set_bounds(0,4.3)
hs=[plt.Line2D([],[],marker='s',ls='',ms=9,color='#7fb0e8',label='白质 GSE180759'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color=BLUE,label='白质 GSE279180'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color=ORANGE,label='皮层灰质 Schirmer'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color='#d8d8d8',label='浅色 = 不显著')]
ax.legend(handles=hs,loc='upper center',bbox_to_anchor=(0.42,-0.13),ncol=2,fontsize=9.3)
ax.set_title('A  同一统计量在两个区室的比较',fontsize=11.5,loc='left',pad=10)

ax2=fig.add_subplot(gs[1])
B=pd.read_csv('/home/user/aaa/results/mda5/out_gm_vs_wm_baseline.csv')
B['cn']=B.cell_type.map(CN); B=B.set_index('cn').loc[['内皮/血管','小胶质/巨噬','星形胶质','少突胶质','OPC','神经元']]
x=np.arange(len(B)); w=0.27
ax2.bar(x-w,B.cortex_GM_rel_astro,w,color=ORANGE,label='皮层灰质')
ax2.bar(x,B.WM_Lerma_rel_astro,w,color=BLUE,label='白质（Lerma）')
ax2.bar(x+w,B.WM_Absinta_rel_astro,w,color='#7fb0e8',label='白质（Absinta）')
ax2.axhline(1,color='#b0b0b0',ls='--',lw=1.1)
ax2.set_xticks(x); ax2.set_xticklabels(B.index,rotation=32,ha='right',fontsize=10.5)
ax2.set_ylabel('相对同一数据集星形胶质的倍数'); ax2.legend(fontsize=9.5)
ax2.set_title('B  基线次序在两个区室一致',fontsize=11.5,loc='left',pad=10)
fig.text(0.008,-0.155,'A：白质每个细胞类型有两条（GSE180759 与 GSE279180），灰质一条（Schirmer 2019）。少突胶质在白质显著上升、在灰质不动，\n    且灰质每组有 3 070 个核，不是功效不足。神经元只有灰质才测得到，结果为阴性（n=10 134）。',fontsize=10.2,color=MUT)
plt.savefig(OUT+'fig3_gm_vs_wm.png'); plt.close(); print('fig3 ok')

# ================= 图4 病灶分期与小胶质状态 =================
A=pd.read_parquet('df_180759.parquet'); A['cp10k']=1e4*A.IFIH1/A.total
STAGES=['control_white_matter','MS_periplaque_white_matter','chronic_inactive_MS_lesion_edge',
        'chronic_active_MS_lesion_edge','MS_lesion_core']
fig=plt.figure(figsize=(13.2,4.6)); gs=fig.add_gridspec(1,3,width_ratios=[1,1.05,1.1],wspace=0.42)
ax=fig.add_subplot(gs[0])
vals=[];ns=[]
for s in STAGES:
    d=A[(A.ct=='immune')&(A.path==s)]
    vals.append(1e4*d.IFIH1.sum()/d.total.sum()); ns.append(len(d))
ax.bar(range(5),vals,color=[STAGE[s] for s in STAGES],width=0.68)
for i,(v,n) in enumerate(zip(vals,ns)):
    ax.text(i,v+0.020,f'{v:.3f}',ha='center',fontsize=9.6,color=INK)
    ax.text(i,v+0.006,f'n={n}',ha='center',fontsize=8.1,color=MUT)
ax.set_xticks(range(5)); ax.set_xticklabels([STAGE_CN[s] for s in STAGES],rotation=32,ha='right',fontsize=10)
ax.set_ylabel('IFIH1  CP10k'); ax.set_ylim(0,0.47)
ax.set_title('A  小胶质/巨噬：随病灶分期递增\nGSE180759',fontsize=11.5,loc='left',pad=10)

ax2=fig.add_subplot(gs[1])
S=pd.read_csv('/home/user/aaa/results/mda5/out_lerma_microglia_subtype.csv')
SUB={'MG_CA':'慢性活动态','MG_Rim':'病灶边缘态','MG_Dis':'疾病相关态','MG_Homeo3':'稳态3',
     'MG_NA':'未分类','MG_Phago1':'吞噬态','MG_Homeo1':'稳态1','MG_Homeo2':'稳态2'}
S['cn']=S.subtype.map(SUB); S=S.sort_values('cp10k')
cols=[RED if s in ('MG_CA','MG_Rim') else BLUE if s=='MG_Dis' else GREY for s in S.subtype]
ax2.barh(S.cn,S.cp10k,color=cols,height=0.66)
for y,(v,n) in enumerate(zip(S.cp10k,S.n)): ax2.text(v+0.012,y,f'{v:.3f}  n={int(n)}',va='center',fontsize=8.8,color=MUT)
ax2.set_xlabel('IFIH1  CP10k'); ax2.set_xlim(0,0.99)
ax2.set_title('B  小胶质状态梯度\nGSE279180（深度归一化）',fontsize=11.5,loc='left',pad=10)

ax3=fig.add_subplot(gs[2])
import anndata as ad, scipy.sparse as sp
a=ad.read_h5ad('lerma/GSE279180_ctype_MG.h5ad',backed='r')
um=np.asarray(a.obsm['X_umap'])
j=int(np.where(a.var_names=='IFIH1')[0][0])
D=pd.read_parquet('lerma_ifih1.parquet'); D=D[D.ctype=='microglia']
v=np.log1p(1e4*D.IFIH1.values/D.total.values)
zero=v<=0; pos=~zero
ax3.scatter(um[zero,0],um[zero,1],c='#e2e2e2',s=2.2,linewidths=0)
o=np.argsort(v[pos])
vmax=np.percentile(v[pos],97) if pos.sum() else 1
ax3.scatter(um[pos][o,0],um[pos][o,1],c=v[pos][o],s=7.0,cmap='YlOrRd',vmin=0,vmax=vmax,linewidths=0)
ax3.set_xticks([]); ax3.set_yticks([])
for sp_ in ['left','bottom']: ax3.spines[sp_].set_visible(False)
sm=plt.cm.ScalarMappable(cmap='YlOrRd',norm=plt.Normalize(0,vmax))
cb=plt.colorbar(sm,ax=ax3,fraction=0.045,pad=0.02); cb.set_label('log1p(CP10k)',fontsize=9.5); cb.outline.set_visible(False)
ax3.set_title(f'C  小胶质 UMAP：IFIH1 阳性核的分布\nGSE279180，9 239 个核中 {int(pos.sum())} 个阳性',fontsize=11.5,loc='left',pad=10)
fig.text(0.008,-0.06,'A：柱上方为该组核数。B：红色为慢性活动态与病灶边缘态，即驱动病灶扩张的那一群，灰色为稳态小胶质。C：灰点为 IFIH1 未检出的核。',fontsize=10.2,color=MUT)
plt.savefig(OUT+'fig4_lesion_stage.png'); plt.close(); print('fig4 ok')

# ================= 图5 关键对照 =================
fig=plt.figure(figsize=(13.2,4.5)); gs=fig.add_gridspec(1,3,width_ratios=[1.05,1,1.05],wspace=0.42)
depth={'小胶质/巨噬':(927,3564),'少突胶质':(2358,3816),'星形胶质':(4033,4334),'OPC':(4158,5416),'内皮/血管':(3605,5790),'神经元':(25058,6959)}
ax=fig.add_subplot(gs[0]); k=list(depth); x=np.arange(len(k)); w=0.36
ax.bar(x-w/2,[depth[i][0] for i in k],w,color=GREY,label='对照')
ax.bar(x+w/2,[depth[i][1] for i in k],w,color=BLUE,label='MS')
for i,kk in enumerate(k):
    r=depth[kk][1]/depth[kk][0]
    ax.text(i,max(depth[kk])+700,f'{r:.2f}×',ha='center',fontsize=9,color=RED if r>1.5 or r<0.7 else MUT)
ax.set_xticks(x); ax.set_xticklabels(k,rotation=32,ha='right',fontsize=10)
ax.set_ylabel('每核中位测序计数'); ax.legend(fontsize=9.5)
ax.set_title('A  测序深度与疾病状态混杂\nGSE279180',fontsize=11.5,loc='left',pad=10)

ax2=fig.add_subplot(gs[1])
raw={'小胶质/巨噬':(4.13,1.64),'少突胶质':(2.84,2.40),'星形胶质':(1.22,1.14),'内皮/血管':(1.51,0.83),'神经元':(0.53,0.88),'OPC':(1.10,0.72)}
k2=list(raw); x=np.arange(len(k2))
for i,kk in enumerate(k2):
    ax2.plot([0,1],[raw[kk][0],raw[kk][1]],'-o',color=RED if kk=='小胶质/巨噬' else (BLUE if kk=='少突胶质' else '#c8c8c8'),
             ms=6,lw=2 if kk in ('小胶质/巨噬','少突胶质') else 1.3,zorder=5 if kk in ('小胶质/巨噬','少突胶质') else 1)
    dy={'神经元':0.075,'内皮/血管':-0.075,'OPC':-0.02}.get(kk,0)
    ax2.text(1.04,raw[kk][1]+dy,kk,fontsize=9.2,va='center',
             color=RED if kk=='小胶质/巨噬' else (BLUE if kk=='少突胶质' else MUT))
ax2.axhline(1,color='#b0b0b0',ls='--',lw=1.1)
ax2.set_xticks([0,1]); ax2.set_xticklabels(['未校正','深度配平后'],fontsize=10.5); ax2.set_xlim(-0.12,1.75)
ax2.set_ylabel('MS / 对照  优势比')
ax2.set_title('B  深度校正前后的效应量\nGSE279180',fontsize=11.5,loc='left',pad=10)

ax3=fig.add_subplot(gs[2])
sets=['GSE180759\n小胶质','GSE279180\n少突胶质','GSE138614\nbulk 活动病灶']
ifih=[2.45,2.45,1.47]; ddx=[0.65,1.42,1.00]
x=np.arange(3); w=0.34
ax3.bar(x-w/2,ifih,w,color=BLUE,label='IFIH1（MDA5）')
ax3.bar(x+w/2,ddx,w,color=GREY,label='DDX58（RIG-I）')
for i in range(3):
    ax3.text(i-w/2,ifih[i]+0.05,f'{ifih[i]:.2f}',ha='center',fontsize=9.2,color=INK)
    ax3.text(i+w/2,ddx[i]+0.05,f'{ddx[i]:.2f}'+('*' if i==2 else ''),ha='center',fontsize=9.2,color=MUT)
ax3.axhline(1,color='#b0b0b0',ls='--',lw=1.1)
ax3.set_xticks(x); ax3.set_xticklabels(sets,fontsize=9.6)
ax3.set_ylabel('MS / 对照  倍数变化'); ax3.set_ylim(0,3.0); ax3.legend(fontsize=9.5)
ax3.set_title('C  MDA5 相对 RIG-I 的选择性',fontsize=11.5,loc='left',pad=10)
fig.text(0.008,-0.075,'A：MS 的小胶质测序深度是对照的 3.84 倍，检出率类指标必须先配平。B：配平把小胶质的效应从 4.13 降到 1.64，少突胶质基本不变。\nC：三个独立设定中 IFIH1 上升而其旁系 RIG-I 不升；* 表示该基因在原研究中任何病灶类型均未达显著。',fontsize=10.2,color=MUT)
plt.savefig(OUT+'fig5_controls.png'); plt.close(); print('fig5 ok')
