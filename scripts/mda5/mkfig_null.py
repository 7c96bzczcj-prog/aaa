import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
N=pd.read_csv('/home/user/aaa/results/mda5/out_gm_wm_genomewide_null.csv')
CELLS=['少突胶质','小胶质/巨噬','星形胶质','OPC','神经元','内皮/血管']
COH={'GSE180759':'#7fb0e8','GSE279180':'#2A78D6'}
fig,ax=plt.subplots(figsize=(10.4,5.6))
for i,cn in enumerate(CELLS):
    for coh,c in COH.items():
        s=N[(N.细胞类型==cn)&(N.白质队列==coh)]
        ctl=s[s.状态=='Control']; ms=s[s.状态=='MS']
        if not len(ctl) or not len(ms): continue
        y0=float(ctl.百分位.iloc[0]); y1=float(ms.百分位.iloc[0])
        x=i+(-0.13 if coh=='GSE180759' else 0.13)
        up=y1>y0
        ax.annotate('',xy=(x,y1),xytext=(x,y0),
                    arrowprops=dict(arrowstyle='-|>',lw=2.4,color=c if up else '#c4c4c4',
                                    shrinkA=0,shrinkB=0,mutation_scale=15))
        ax.scatter([x],[y0],s=38,color='white',edgecolors=c if up else '#c4c4c4',linewidths=1.8,zorder=5)
        ax.scatter([x],[y1],s=46,color=c if up else '#c4c4c4',zorder=5)
ax.axhline(50,color='#b0b0b0',ls='--',lw=1.1)
ax.text(5.62,51.5,'全基因组中位',fontsize=9.5,color=MUT,ha='right')
ax.set_xticks(range(len(CELLS))); ax.set_xticklabels(CELLS,fontsize=12)
ax.set_ylim(0,102); ax.set_ylabel('IFIH1 的 白质/灰质 比值\n在全基因组同类比值中的百分位')
ax.set_xlim(-0.55,5.7)
hs=[plt.Line2D([],[],marker='o',ls='',ms=8,mfc='white',mec='#555',label='正常（空心）'),
    plt.Line2D([],[],marker='o',ls='',ms=8,color='#555',label='MS（实心）'),
    plt.Line2D([],[],lw=2.4,color=COH['GSE180759'],label='GSE180759'),
    plt.Line2D([],[],lw=2.4,color=COH['GSE279180'],label='GSE279180'),
    plt.Line2D([],[],lw=2.4,color='#c4c4c4',label='未上升')]
ax.legend(handles=hs,loc='upper center',bbox_to_anchor=(0.5,-0.10),ncol=5,fontsize=9.6)
plt.savefig('/home/user/aaa/results/mda5/figs/figB_null.png',dpi=225,bbox_inches='tight')
print('ok')
