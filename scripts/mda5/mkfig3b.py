import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
import matplotlib.ticker as mt
T=pd.read_csv('/home/user/aaa/results/mda5/out_gm_wm_depthmatched.csv')
CN={'microglia':'小胶质/巨噬','oligodendrocyte':'少突胶质','astrocyte':'星形胶质',
    'OPC':'OPC','endothelial':'内皮/血管','neuron':'神经元'}
T['cn']=T.cell_type.map(CN)
order=['少突胶质','小胶质/巨噬','星形胶质','神经元','OPC','内皮/血管']
SLOTS=[('white matter','GSE180759 Absinta','#7fb0e8'),
       ('white matter','GSE279180 Lerma-Martin','#2A78D6'),
       ('cortex / grey matter','Schirmer 2019','#EB6834')]
fig,ax=plt.subplots(figsize=(9.6,5.4))
h=0.24; LABX=4.45
for i,cn in enumerate(order):
    for k,(comp,coh,col) in enumerate(SLOTS):
        r=T[(T.cn==cn)&(T.compartment==comp)&(T.cohort==coh)]
        y=i+(k-1)*h
        if r.empty: continue
        r=r.iloc[0]
        c=col if r.p<0.05 else ('#c9dcf2' if comp=='white matter' else '#f7cdb8')
        ax.barh(y,min(r.OR,4.3),height=h*0.86,color=c)
        if r.OR>4.3: ax.text(4.33,y,'▸',va='center',fontsize=9,color=c)
        ptxt='p<0.005' if r.p<0.005 else f'p={r.p:.2f}'
        ax.text(LABX,y,f'{r.OR:.2f}',va='center',ha='right',fontsize=9.6,
                color=INK if r.p<0.05 else MUT,fontweight='bold' if r.p<0.05 else 'normal')
        ax.text(LABX+0.28,y,ptxt,va='center',ha='left',fontsize=8.9,color=INK if r.p<0.05 else MUT)
        ax.text(LABX+2.05,y,f'{int(r.n_per_group)} 核',va='center',ha='right',fontsize=8.9,color=MUT)
ax.axvline(1,color='#b0b0b0',ls='--',lw=1.1)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=12)
ax.set_xlim(0,6.7); ax.set_xlabel('MS / 对照  检出率优势比（深度十分位配平后）')
ax.invert_yaxis()
ax.xaxis.set_major_locator(mt.FixedLocator([0,1,2,3,4]))
ax.spines['bottom'].set_bounds(0,4.3)
hs=[plt.Line2D([],[],marker='s',ls='',ms=9,color='#7fb0e8',label='白质 GSE180759'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color='#2A78D6',label='白质 GSE279180'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color='#EB6834',label='灰质 Schirmer 2019'),
    plt.Line2D([],[],marker='s',ls='',ms=9,color='#d8d8d8',label='浅色 = 不显著')]
ax.legend(handles=hs,loc='upper center',bbox_to_anchor=(0.40,-0.12),ncol=2,fontsize=9.6)
plt.savefig('/home/user/aaa/results/mda5/figs/fig3_gm_vs_wm.png',dpi=225,bbox_inches='tight')
print('ok')
