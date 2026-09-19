import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
a=ad.read_h5ad('lerma/GSE279180_sn_atlas.h5ad',backed='r')
j=int(np.where(a.var_names=='IFIH1')[0][0])
n=a.shape[0]; tot=np.zeros(n); g=np.zeros(n); i=0
for ch in a.chunked_X(8000):
    C=sp.csr_matrix(ch[0] if isinstance(ch,tuple) else ch); m=C.shape[0]
    tot[i:i+m]=np.asarray(C.sum(axis=1)).ravel()
    g[i:i+m]=np.asarray(C[:,j].todense()).ravel(); i+=m
o=a.obs[['celltype','subtype','lesion_type','condition','patient_id','sample_id']].copy().reset_index(drop=True)
o['total']=tot; o['IFIH1']=g
um=np.asarray(a.obsm['X_umap']); o['umap1']=um[:,0]; o['umap2']=um[:,1]
o.to_parquet('atlas_ifih1.parquet')
print('done',len(o),'IFIH1+',int((g>0).sum()))
print(o.groupby(['celltype','lesion_type'],observed=True).size().unstack(fill_value=0).to_string())
