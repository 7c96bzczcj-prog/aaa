import pandas as pd, numpy as np, gzip, time
t0=time.time()
m=pd.read_parquet('schirmer_ifih1.parquet')
grp=(m.grp.astype(str)+'|'+m.cond.astype(str)).values
keys=sorted(set(grp)); kidx={k:i for i,k in enumerate(keys)}
gi=np.array([kidx[g] for g in grp])
names=[]; mat=[]
with gzip.open('schirmer/exprMatrix.tsv.gz','rt') as f:
    f.readline()
    for n,line in enumerate(f):
        p=line.rstrip('\n').split('\t')
        v=np.array(p[1:],dtype=np.float32)
        if (v>0).sum()<50: continue
        s=np.bincount(gi,weights=np.expm1(v).astype(np.float64),minlength=len(keys))
        names.append(p[0].split('|')[1]); mat.append(s)
        if n%4000==0: print(n,'rows |',len(names),'kept |',round(time.time()-t0),'s',flush=True)
M=np.vstack(mat)
pd.DataFrame(M,index=names,columns=keys).to_parquet('schirmer_pseudobulk.parquet')
print('DONE',M.shape,round(time.time()-t0),'s',flush=True)
