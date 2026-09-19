import numpy as np, pandas as pd, pickle, csv
from scipy import stats
hdr=open('header.txt').read().strip().split(',')
genes=pickle.load(open('genes.pkl','rb'))
ann={r['nucleus_barcode']:r for r in csv.DictReader(open('GSE180759_annotation.txt'),delimiter='\t')}
lib=np.load('libsize.npy')
df=pd.DataFrame({'ct':[ann[b]['cell_type'] for b in hdr],'donor':[ann[b]['NBB_case'] for b in hdr],
                 'path':[ann[b]['pathology'] for b in hdr],'clus':[int(ann[b]['seurat_cluster']) for b in hdr],
                 'libr':[int(b.split('-')[-1]) for b in hdr],'total':lib})
for g in genes: df[g]=genes[g]
df['cond']=np.where(df.path=='control_white_matter','Control','MS')
ORD=['oligodendrocytes','astrocytes','immune','neurons','opc','vascular_cells','lymphocytes']
PORD=['control_white_matter','MS_periplaque_white_matter','chronic_inactive_MS_lesion_edge',
      'chronic_active_MS_lesion_edge','MS_lesion_core']
print(f"cells {len(df)} | median total counts {np.median(lib):.0f} | total UMI {lib.sum():.3g}")

print("\n########## A. DEPTH by cell type and condition (the confounder to rule out) ##########")
print(f"{'cell_type':<20}{'ctrl_median':>13}{'MS_median':>11}{'MS/ctrl':>9}")
for c in ORD:
    s=df[df.ct==c]; a=s[s.cond=='MS'].total.median(); b=s[s.cond=='Control'].total.median()
    print(f"{c:<20}{b:>13.0f}{a:>11.0f}{a/b:>9.2f}")

print("\n########## B. BASELINE IFIH1 by cell type (depth-normalised CP10k) ##########")
print(f"{'cell_type':<20}{'n_all':>8}{'CP10k_all':>11}{'pct_pos':>9}   |{'n_ctrl':>8}{'CP10k_ctrl':>12}")
rows=[]
for c in ORD:
    s=df[df.ct==c]; sc=s[s.cond=='Control']
    cp=1e4*s[c=='x'].sum() if False else 1e4*s.IFIH1.sum()/s.total.sum()
    cpc=1e4*sc.IFIH1.sum()/sc.total.sum() if len(sc) else np.nan
    print(f"{c:<20}{len(s):>8}{cp:>11.3f}{100*(s.IFIH1>0).mean():>8.2f}%   |{len(sc):>8}{cpc:>12.3f}")
    rows.append((c,len(s),cp,100*(s.IFIH1>0).mean(),len(sc),cpc))
pd.DataFrame(rows,columns=['cell_type','n','cp10k_all','pct_pos','n_ctrl','cp10k_ctrl']).to_csv('out_180759_baseline.csv',index=False)

print("\n########## C. MS vs control, DONOR-LEVEL, depth-normalised CP10k ##########")
print(f"{'cell_type':<20}{'n_c_don':>8}{'n_MS_don':>9}{'ctrl_med':>10}{'MS_med':>9}{'ratio':>7}{'MW_p':>8}{'sep':>10}")
res=[]
for c in ORD:
    s=df[df.ct==c]
    per=s.groupby(['donor','cond'],observed=True).apply(
        lambda d: pd.Series({'n':len(d),'cp10k':1e4*d.IFIH1.sum()/d.total.sum()}),include_groups=False).reset_index()
    per=per[per.n>=20]
    ms=per[per.cond=='MS'].cp10k.values; cs=per[per.cond=='Control'].cp10k.values
    if len(ms)<3 or len(cs)<2:
        print(f"{c:<20}{len(cs):>8}{len(ms):>9}   -- too few donors --"); continue
    p=stats.mannwhitneyu(ms,cs,alternative='two-sided').pvalue
    sep='COMPLETE' if ms.min()>cs.max() else 'overlap'
    print(f"{c:<20}{len(cs):>8}{len(ms):>9}{np.median(cs):>10.3f}{np.median(ms):>9.3f}{np.median(ms)/max(np.median(cs),1e-9):>7.2f}{p:>8.4f}{sep:>10}")
    res.append((c,len(cs),len(ms),np.median(cs),np.median(ms),p,sep))
pd.DataFrame(res,columns=['cell_type','n_ctrl_don','n_MS_don','cp10k_ctrl','cp10k_MS','MW_p','separation']).to_csv('out_180759_cp10k_donor.csv',index=False)

print("\n########## D. DEPTH-MATCHED re-test (immune cells) ##########")
M=df[df.ct=='immune']; c=M[M.cond=='Control']; m=M[M.cond=='MS']
bins=np.percentile(np.concatenate([c.total,m.total]),np.arange(0,101,10))
tc=tm=pc=pm=0
for i in range(len(bins)-1):
    cc=c[(c.total>=bins[i])&(c.total<bins[i+1])]; mm=m[(m.total>=bins[i])&(m.total<bins[i+1])]
    k=min(len(cc),len(mm))
    if k<5: continue
    tc+=k; tm+=k; pc+=(cc.sample(k,random_state=1).IFIH1>0).sum(); pm+=(mm.sample(k,random_state=1).IFIH1>0).sum()
orr,p=stats.fisher_exact([[int(pm),tm-int(pm)],[int(pc),tc-int(pc)]])
print(f"  depth-matched n={tc}/group: control {100*pc/tc:.2f}%  MS {100*pm/tm:.2f}%  OR={orr:.2f}  p={p:.2e}")

print("\n########## E. IFIH1 CP10k by cell type x lesion region ##########")
print(f"{'cell_type':<20}"+"".join(f"{p[:14]:>16}" for p in PORD))
for c in ORD:
    row=[]
    for p in PORD:
        s=df[(df.ct==c)&(df.path==p)]
        row.append(1e4*s.IFIH1.sum()/s.total.sum() if s.total.sum()>0 else np.nan)
    print(f"{c:<20}"+"".join(f"{v:>16.3f}" for v in row))

print("\n########## F. SPECIFICITY: RLR/ISG panel in immune cells, CP10k, MS vs control ##########")
print(f"{'gene':<9}{'ctrl_CP10k':>12}{'MS_CP10k':>11}{'ratio':>8}")
for g in ['IFIH1','DDX58','DHX58','MAVS','STAT1','MX1','ISG15','IFIT3','OAS1','BST2','IRF7','XAF1','EIF2AK2']:
    if g not in df: continue
    a=M[M.cond=='MS']; b=M[M.cond=='Control']
    ca=1e4*a[g].sum()/a.total.sum(); cb=1e4*b[g].sum()/b.total.sum()
    print(f"{g:<9}{cb:>12.3f}{ca:>11.3f}{(ca/cb if cb>0 else np.inf):>8.2f}")

print("\n########## G. microglia subclusters (c5/c10/c17) CP10k by condition ##########")
for cl in [5,10,17,16,11,13]:
    s=df[df.clus==cl]; a=s[s.cond=='MS']; b=s[s.cond=='Control']
    lbl=s.ct.iloc[0]
    ca=1e4*a.IFIH1.sum()/a.total.sum() if a.total.sum()>0 else np.nan
    cb=1e4*b.IFIH1.sum()/b.total.sum() if b.total.sum()>0 else np.nan
    print(f"  c{cl:<3} {lbl:<16} n={len(s):>5}  ctrl CP10k {cb:>7.3f} (n={len(b)})   MS CP10k {ca:>7.3f} (n={len(a)})")
df[['ct','donor','path','clus','libr','total','IFIH1','DDX58','MX1','BST2','STAT1']].to_parquet('df_180759.parquet')

print("\n########## H. COMPOSITION CONTROL: is the rise within true microglia, or a myeloid shift? ##########")
print("Alternative explanation: MS recruits CD163+/MRC1+ macrophages, which may carry more IFIH1.")
print("Test restricted to P2RY12-positive nuclei (homeostatic-microglia-defining, not macrophage).")
for label,mask in [('P2RY12+ only', (df.ct=='immune')&(df.P2RY12>0)),
                   ('CD163+ or MRC1+', (df.ct=='immune')&((df.CD163>0)|(df.MRC1>0))),
                   ('all immune', df.ct=='immune')]:
    s=df[mask]; a=s[s.cond=='MS']; b=s[s.cond=='Control']
    if len(b)<10: print(f"  {label:<18} control n={len(b)} too few"); continue
    ca=1e4*a.IFIH1.sum()/a.total.sum(); cb=1e4*b.IFIH1.sum()/b.total.sum()
    orr,p=stats.fisher_exact([[int((a.IFIH1>0).sum()),int((a.IFIH1==0).sum())],
                              [int((b.IFIH1>0).sum()),int((b.IFIH1==0).sum())]])
    print(f"  {label:<18} ctrl n={len(b):>4} CP10k {cb:>6.3f} ({100*(b.IFIH1>0).mean():>5.2f}%) | MS n={len(a):>5} CP10k {ca:>6.3f} ({100*(a.IFIH1>0).mean():>5.2f}%) | OR={orr:>5.2f} p={p:.2e}")

print("\n########## I. baseline IFIH1 vs macrophage markers (is IFIH1 a macrophage gene here?) ##########")
im=df[df.ct=='immune']
for grp,m in [('P2RY12+ microglia', im.P2RY12>0), ('CD163+/MRC1+ macrophage', (im.CD163>0)|(im.MRC1>0))]:
    s=im[m]; print(f"  {grp:<26} n={len(s):>5}  IFIH1 CP10k {1e4*s.IFIH1.sum()/s.total.sum():>6.3f}  pct+ {100*(s.IFIH1>0).mean():>5.2f}%")

import os
if os.path.exists('gene_rows2.csv'):
    g2={}
    for line in open('gene_rows2.csv'):
        p=line.rstrip('\n').split(','); g2[p[0]]=np.array(p[1:],dtype=np.int32)
    for k,v in g2.items(): df[k]=v
    print("\n########## J. ADAR and the wider dsRNA-sensing panel in immune cells (CP10k) ##########")
    print("ADAR1 edits endogenous dsRNA and is the brake on MDA5; loss of that brake activates MDA5.")
    print(f"{'gene':<10}{'ctrl_CP10k':>12}{'MS_CP10k':>11}{'ratio':>8}{'fisher_p':>11}")
    a=df[(df.ct=='immune')&(df.cond=='MS')]; b=df[(df.ct=='immune')&(df.cond=='Control')]
    for g in ['IFIH1','ADAR','ADARB1','ZBP1','TLR3','DDX60','IFI16','OAS2','OAS3','RNASEL','NLRP3',
              'TREM2','APOE','SPP1','GPNMB','CD68','HLA-DRA','CD74','CXCL10','IL1B','P2RY13','SALL1','MERTK','ITGAX']:
        if g not in df: continue
        ca=1e4*a[g].sum()/a.total.sum(); cb=1e4*b[g].sum()/b.total.sum()
        orr,p=stats.fisher_exact([[int((a[g]>0).sum()),int((a[g]==0).sum())],[int((b[g]>0).sum()),int((b[g]==0).sum())]])
        print(f"{g:<10}{cb:>12.3f}{ca:>11.3f}{(ca/cb if cb>0 else np.inf):>8.2f}{p:>11.2e}")
