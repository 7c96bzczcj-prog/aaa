import sys; sys.path.insert(0,'.')
from figstyle import *
import pandas as pd, numpy as np, pickle, csv
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
um=np.load('umap_coords.npy')
hdr=open('header.txt').read().strip().split(',')
ann={r['nucleus_barcode']:r for r in csv.DictReader(open('GSE180759_annotation.txt'),delimiter='\t')}
genes=pickle.load(open('genes.pkl','rb')); lib=np.load('libsize.npy')
D=pd.DataFrame({'ct':[ann[b]['cell_type'] for b in hdr],'path':[ann[b]['pathology'] for b in hdr],
                'umap1':um[:,0],'umap2':um[:,1],'IFIH1':genes['IFIH1'],'total':lib})
D['lg']=np.log1p(1e4*D.IFIH1/D.total)
CT={'oligodendrocytes':'少突胶质','astrocytes':'星形胶质','immune':'小胶质/巨噬','neurons':'神经元',
    'opc':'OPC','vascular_cells':'血管/内皮','lymphocytes':'淋巴细胞'}
COL={'oligodendrocytes':'#E8A33D','astrocytes':'#E4572E','immune':'#2E933C','neurons':'#17A398',
     'opc':'#2A78D6','vascular_cells':'#7B4FA8','lymphocytes':'#C64191'}
ORD=['control_white_matter','MS_periplaque_white_matter','chronic_inactive_MS_lesion_edge',
     'chronic_active_MS_lesion_edge','MS_lesion_core']
PUR=LinearSegmentedColormap.from_list('pur',['#D9D2E9','#8E76C1','#4B2E83','#2A1352'])

fig=plt.figure(figsize=(14.2,5.4))
gs=fig.add_gridspec(1,2,width_ratios=[1.0,1.62],wspace=0.05)
ax=fig.add_subplot(gs[0])
rng=np.random.default_rng(0); idx=rng.permutation(len(D))
ax.scatter(D.umap1.values[idx],D.umap2.values[idx],c=[COL[c] for c in D.ct.values[idx]],
           s=1.1,linewidths=0,rasterized=True)
OFF={'lymphocytes':(0,1.2),'vascular_cells':(0,-1.0)}
for c,g in D.groupby('ct',observed=True):
    dx,dy=OFF.get(c,(0,0))
    ax.text(np.median(g.umap1)+dx,np.median(g.umap2)+dy,CT[c],fontsize=10.5,ha='center',va='center',
            color='#111111',fontweight='bold',path_effects=[pe.withStroke(linewidth=3.4,foreground='white')])
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
x0,y0=D.umap1.min()-0.6,D.umap2.min()-0.6
ax.annotate('',xy=(x0+4.0,y0),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>',color='#333',lw=1.1))
ax.annotate('',xy=(x0,y0+4.0),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>',color='#333',lw=1.1))
ax.text(x0+2.0,y0-1.4,'UMAP1',fontsize=9,ha='center',color='#333')
ax.text(x0-1.4,y0+2.0,'UMAP2',fontsize=9,va='center',rotation=90,color='#333')
ax.set_title('66 432 个核，七类细胞',fontsize=12,loc='left',pad=14)

sub=gs[1].subgridspec(2,3,wspace=0.04,hspace=0.30)
vmax=np.percentile(D.lg[D.lg>0],97)
for k,p in enumerate(ORD):
    r,c=divmod(k,3); axk=fig.add_subplot(sub[r,c])
    axk.scatter(D.umap1,D.umap2,c='#E6E6E6',s=0.6,linewidths=0,rasterized=True)
    d=D[(D.path==p)&(D.lg>0)].sort_values('lg')
    axk.scatter(d.umap1,d.umap2,c=d.lg,cmap=PUR,vmin=0,vmax=vmax,s=4.0,linewidths=0,rasterized=True)
    axk.set_xticks([]); axk.set_yticks([])
    for s in axk.spines.values(): s.set_visible(False)
    n=int((D.path==p).sum())
    axk.set_title(f'{STAGE_CN[p]}\n{n:,} 个核 · {100*len(d)/n:.1f}% 检出'.replace(',',' '),fontsize=10.5,pad=4)
axlast=fig.add_subplot(sub[1,2]); axlast.axis('off')
sm=plt.cm.ScalarMappable(cmap=PUR,norm=plt.Normalize(0,vmax))
cax=fig.add_axes([0.80,0.13,0.010,0.26])
cb=plt.colorbar(sm,cax=cax); cb.set_label('IFIH1  log1p(CP10k)',fontsize=9.5); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8.5)
plt.savefig('/home/user/aaa/results/mda5/figs/fig6_umap_180759.png',dpi=230,bbox_inches='tight')
print('ok')
