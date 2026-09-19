import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
from scipy import stats
D=pd.read_parquet('twoway.parquet')
CN={'oligodendrocyte':'少突胶质','microglia':'小胶质/巨噬','astrocyte':'星形胶质',
    'OPC':'OPC','neuron':'神经元','endothelial':'内皮/血管'}
WM='#2A78D6'; GM='#EB6834'
CELLS=['oligodendrocyte','microglia','astrocyte','OPC','neuron','endothelial']
fig,axes=plt.subplots(2,3,figsize=(13.6,6.6))
POS=[('白质','Control',0,WM),('灰质','Control',1,GM),('白质','MS',2.55,WM),('灰质','MS',3.55,GM)]
for ax,g in zip(axes.ravel(),CELLS):
    rng=np.random.default_rng(1); ns={}
    vmax=0
    for comp,cond,x,c in POS:
        s=D[(D.细胞类型==g)&(D.区室==comp)&(D.状态==cond)]
        if not len(s): continue
        med=s.值.median(); vmax=max(vmax,s.值.max())
        ax.bar(x,med,width=0.74,color=c,alpha=0.30,zorder=1)
        ax.hlines(med,x-0.37,x+0.37,color=c,lw=2.4,zorder=4)
        ax.scatter(x+rng.uniform(-0.19,0.19,len(s)),s.值,s=26,color=c,alpha=0.9,
                   linewidths=0.6,edgecolors='white',zorder=5)
        ns[x]=len(s)
    # p 值：同一状态下 白质 vs 灰质
    for cond,xc,lab in [('Control',0.5,'正常'),('MS',3.05,'MS')]:
        w=D[(D.细胞类型==g)&(D.区室=='白质')&(D.状态==cond)].值.values
        m=D[(D.细胞类型==g)&(D.区室=='灰质')&(D.状态==cond)].值.values
        if len(w)>=3 and len(m)>=3:
            p=stats.mannwhitneyu(w,m,alternative='two-sided').pvalue
            t=f'p={p:.3f}' if p>=0.001 else 'p<0.001'
            bold=p<0.05
        else:
            t='供体不足'; bold=False
        y=vmax*1.06
        ax.plot([xc-0.5,xc+0.5],[y,y],color='#555' if bold else '#bbb',lw=1.0)
        ax.text(xc,y*1.02,t,ha='center',fontsize=9,color=INK if bold else MUT,
                fontweight='bold' if bold else 'normal')
    ax.set_xticks([0,1,2.55,3.55])
    ax.set_xticklabels([f'白质\n{ns.get(0,0)} 人',f'灰质\n{ns.get(1,0)} 人',
                        f'白质\n{ns.get(2.55,0)} 人',f'灰质\n{ns.get(3.55,0)} 人'],fontsize=9.5)
    ax.set_xlim(-0.65,4.2); ax.set_ylim(0,vmax*1.22)
    ax.set_title(CN[g],fontsize=12.5,loc='left',pad=8)
    ax.text(0.5,-0.30,'正常',ha='center',fontsize=11,transform=ax.get_xaxis_transform())
    ax.text(3.05,-0.30,'MS',ha='center',fontsize=11,transform=ax.get_xaxis_transform())
    ax.axvline(1.78,color='#e2e2e2',lw=1)
for ax in axes[:,0]: ax.set_ylabel('IFIH1  伪整体 CPM')
plt.subplots_adjust(hspace=0.72,wspace=0.26)
plt.savefig('/home/user/aaa/results/mda5/figs/figA_2x2.png',dpi=225,bbox_inches='tight')
print('ok')
