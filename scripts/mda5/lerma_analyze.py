import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, glob, os
GENES=['IFIH1','DDX58','DHX58','MAVS','STAT1','MX1','ISG15','IFIT3','OAS1','BST2','IRF7','XAF1','EIF2AK2']
NAME={'MG':'microglia','AS':'astrocytes','OL':'oligodendrocytes','OPC':'OPC','NEU':'neurons',
      'EC':'endothelial/vascular','BC':'B cells','TC':'T cells','SC':'stromal/other'}
out=[]
for f in sorted(glob.glob('lerma/GSE279180_ctype_*.h5ad')):
    key=os.path.basename(f).split('_ctype_')[1].replace('.h5ad','')
    a=ad.read_h5ad(f,backed='r')
    idx={g:int(np.where(a.var_names==g)[0][0]) for g in GENES if g in a.var_names}
    n=a.shape[0]; tot=np.zeros(n); vals={g:np.zeros(n) for g in idx}
    i=0
    for ch in a.chunked_X(4000):
        C=sp.csr_matrix(ch[0] if isinstance(ch,tuple) else ch)
        m=C.shape[0]
        tot[i:i+m]=np.asarray(C.sum(axis=1)).ravel()
        for g,j in idx.items(): vals[g][i:i+m]=np.asarray(C[:,j].todense()).ravel()
        i+=m
    o=a.obs.copy().reset_index(drop=True); o['ctype']=NAME.get(key,key); o['total']=tot
    for g in idx: o[g]=vals[g]
    out.append(o); print(f"{key:<4} n={n:<7} cells done",flush=True)
D=pd.concat(out,ignore_index=True)
D.to_parquet('lerma_ifih1.parquet')
print("\nTOTAL cells",len(D))
print(D.groupby(['ctype','condition'],observed=True).size().unstack(fill_value=0).to_string())
print("\npatients:",D.patient_id.nunique(),"| controls:",sorted(D[D.condition=='Control'].patient_id.unique()))
print("MS:",sorted(D[D.condition=='MS'].patient_id.unique()))
