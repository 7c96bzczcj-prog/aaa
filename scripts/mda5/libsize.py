import pandas as pd, numpy as np
tot=None; ngenes=0; names=[]
for ch in pd.read_csv('GSE180759_expression_matrix.csv.gz',index_col=0,chunksize=500):
    v=ch.values.sum(axis=0,dtype=np.int64)
    tot = v if tot is None else tot+v
    ngenes+=ch.shape[0]; names.append(ch.shape[0])
    if len(names)%20==0: print("genes so far",ngenes,flush=True)
np.save('libsize.npy',tot)
print("DONE genes",ngenes,"cells",len(tot),"median",np.median(tot),"min",tot.min(),"total UMI",tot.sum())
