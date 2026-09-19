import pandas as pd, numpy as np, gc, time
t0=time.time()
keep_names=[]; keep_rows=[]
NCELL=66432
LO,HI=1000,40000            # detected in 1.5% .. 60% of nuclei
n=0
for ch in pd.read_csv('GSE180759_expression_matrix.csv.gz',index_col=0,chunksize=500):
    V=ch.values.astype(np.float32)
    nnz=(V>0).sum(axis=1)
    sel=np.where((nnz>=LO)&(nnz<=HI))[0]
    for i in sel:
        keep_names.append(ch.index[i]); keep_rows.append(V[i].copy())
    n+=ch.shape[0]
    if n%5000==0: print(n,'genes scanned |',len(keep_names),'kept |',round(time.time()-t0),'s',flush=True)
    del V; gc.collect()
X=np.vstack(keep_rows).T          # cells x genes
del keep_rows; gc.collect()
print('matrix',X.shape,'time',round(time.time()-t0),'s',flush=True)
np.save('umap_X.npy',X); np.save('umap_genes.npy',np.array(keep_names,dtype=object))
print('SAVED_MATRIX',flush=True)

import scanpy as sc, anndata as ad
a=ad.AnnData(X=X)
a.var_names=[str(g) for g in keep_names]
hdr=open('header.txt').read().strip().split(',')
a.obs_names=hdr
import csv
ann={r['nucleus_barcode']:r for r in csv.DictReader(open('GSE180759_annotation.txt'),delimiter='\t')}
a.obs['cell_type']=[ann[b]['cell_type'] for b in hdr]
a.obs['pathology']=[ann[b]['pathology'] for b in hdr]
a.obs['donor']=[ann[b]['NBB_case'] for b in hdr]
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000); a=a[:,a.var.highly_variable].copy()
sc.pp.scale(a,max_value=10)
sc.tl.pca(a,n_comps=50,svd_solver='arpack')
import scanpy.external as sce
sc.pp.neighbors(a,n_neighbors=15,n_pcs=40)
sc.tl.umap(a,min_dist=0.3)
np.save('umap_coords.npy',a.obsm['X_umap'])
a.obs[['cell_type','pathology','donor']].to_csv('umap_obs.csv')
print('UMAP_DONE',round(time.time()-t0),'s',flush=True)
