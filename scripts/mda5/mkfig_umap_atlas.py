import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np
from matplotlib.colors import LinearSegmentedColormap
D=pd.read_parquet('atlas_ifih1.parquet')
D['cp10k']=1e4*D.IFIH1/D.total
D['lg']=np.log1p(D.cp10k)
CT={'OL':'少突胶质','AS':'星形胶质','MG':'小胶质','NEU':'神经元','OPC':'OPC',
    'EC':'内皮','TC':'T 细胞','SC':'基质','BC':'B 细胞'}
COL={'OL':'#E8A33D','AS':'#E4572E','MG':'#2E933C','NEU':'#17A398','OPC':'#2A78D6',
     'EC':'#7B4FA8','TC':'#C64191','SC':'#8C6239','BC':'#5B7C99'}
LES=[('Ctrl','对照白质'),('CA','慢性活动病灶'),('CI','慢性非活动病灶')]
PUR=LinearSegmentedColormap.from_list('pur',['#D9D2E9','#8E76C1','#4B2E83','#2A1352'])

fig=plt.figure(figsize=(14.2,5.2))
gs=fig.add_gridspec(1,2,width_ratios=[1.06,1.55],wspace=0.06)

# ---- left: compartment UMAP ----
ax=fig.add_subplot(gs[0])
rng=np.random.default_rng(0); idx=rng.permutation(len(D))
ax.scatter(D.umap1.values[idx],D.umap2.values[idx],c=[COL[c] for c in D.celltype.values[idx]],
           s=1.0,linewidths=0,rasterized=True)
import matplotlib.patheffects as pe
OFF={'BC':(-2.6,0.9),'SC':(2.4,-1.1),'TC':(0,-0.9)}
for c,g in D.groupby('celltype',observed=True):
    dx,dy=OFF.get(c,(0,0))
    ax.text(np.median(g.umap1)+dx,np.median(g.umap2)+dy,CT[c],fontsize=10.5,ha='center',va='center',
            color='#111111',fontweight='bold',
            path_effects=[pe.withStroke(linewidth=3.4,foreground='white')])
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
x0,y0=D.umap1.min()-0.6,D.umap2.min()-0.6
ax.annotate('',xy=(x0+4.2,y0),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>',color='#333',lw=1.1))
ax.annotate('',xy=(x0,y0+4.2),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>',color='#333',lw=1.1))
ax.text(x0+2.1,y0-1.3,'UMAP1',fontsize=9,ha='center',color='#333')
ax.text(x0-1.3,y0+2.1,'UMAP2',fontsize=9,va='center',rotation=90,color='#333')
ax.set_title('103 794 个核，九类细胞',fontsize=12,loc='left',pad=14)

# ---- right: IFIH1 feature plots by lesion type ----
sub=gs[1].subgridspec(1,3,wspace=0.04)
vmax=np.percentile(D.lg[D.lg>0],97)
for k,(code,name) in enumerate(LES):
    axk=fig.add_subplot(sub[k])
    axk.scatter(D.umap1,D.umap2,c='#E6E6E6',s=0.7,linewidths=0,rasterized=True)
    d=D[(D.lesion_type==code)&(D.lg>0)].sort_values('lg')
    axk.scatter(d.umap1,d.umap2,c=d.lg,cmap=PUR,vmin=0,vmax=vmax,s=3.4,linewidths=0,rasterized=True)
    axk.set_xticks([]); axk.set_yticks([])
    for s in axk.spines.values(): s.set_visible(False)
    n=int((D.lesion_type==code).sum()); npos=len(d)
    axk.set_title(f'{name}\n{n:,} 个核 · {100*npos/n:.1f}% 检出'.replace(',',' '),fontsize=11,pad=6)
sm=plt.cm.ScalarMappable(cmap=PUR,norm=plt.Normalize(0,vmax))
cax=fig.add_axes([0.915,0.22,0.009,0.5])
cb=plt.colorbar(sm,cax=cax); cb.set_label('IFIH1  log1p(CP10k)',fontsize=9.5); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8.5)
plt.savefig('/home/user/aaa/results/mda5/figs/fig0_umap_atlas.png',dpi=230,bbox_inches='tight')
print('ok')
