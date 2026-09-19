import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, time, csv
t0=time.time()
# ---------- GSE279180 atlas ----------
a=ad.read_h5ad('lerma/GSE279180_sn_atlas.h5ad',backed='r')
o=a.obs[['celltype','condition']].copy().reset_index(drop=True)
grp=(o.celltype.astype(str)+'|'+o.condition.astype(str)).values
keys=sorted(set(grp)); kidx={k:i for i,k in enumerate(keys)}
gi=np.array([kidx[g] for g in grp])
G=a.shape[1]; M=np.zeros((G,len(keys)))
i=0
for ch in a.chunked_X(8000):
    C=sp.csr_matrix(ch[0] if isinstance(ch,tuple) else ch); m=C.shape[0]
    sub=gi[i:i+m]
    for k in range(len(keys)):
        sel=np.where(sub==k)[0]
        if len(sel): M[:,k]+=np.asarray(C[sel].sum(axis=0)).ravel()
    i+=m
    print('atlas',i,round(time.time()-t0),'s',flush=True)
pd.DataFrame(M,index=[str(v) for v in a.var_names],columns=keys).to_parquet('lerma_pseudobulk.parquet')
print('ATLAS_DONE',M.shape,flush=True)

# ---------- GSE180759 ----------
X=np.load('umap_X.npy',mmap_mode='r')
gn=np.load('umap_genes.npy',allow_pickle=True)
hdr=open('header.txt').read().strip().split(',')
ann={r['nucleus_barcode']:r for r in csv.DictReader(open('GSE180759_annotation.txt'),delimiter='\t')}
ct=np.array([ann[b]['cell_type'] for b in hdr])
cond=np.array(['Control' if ann[b]['pathology']=='control_white_matter' else 'MS' for b in hdr])
g2=np.char.add(np.char.add(ct.astype(str),'|'),cond.astype(str))
k2=sorted(set(g2)); ki={k:i for i,k in enumerate(k2)}
gi2=np.array([ki[g] for g in g2])
M2=np.zeros((X.shape[1],len(k2)))
for k in range(len(k2)):
    sel=np.where(gi2==k)[0]
    for s in range(0,len(sel),4000):
        M2[:,k]+=X[sel[s:s+4000]].sum(axis=0)
    print('180759',k2[k],round(time.time()-t0),'s',flush=True)
pd.DataFrame(M2,index=[str(g) for g in gn],columns=k2).to_parquet('abs_pseudobulk.parquet')
print('ALL_DONE',M2.shape,round(time.time()-t0),'s',flush=True)
