import sys, csv, os, statistics
import numpy as np
d=sys.argv[1]; out=sys.argv[2]
labels={'B_CELL_NAIVE':'B cell, naive','MONOCYTES':'Monocyte, classical','M2':'Monocyte, non-classical','NK':'NK cell, CD56dim CD16+','TREG_MEM':'T cell, CD4, memory TREG','CD4_NAIVE':'T cell, CD4, naive','CD4_STIM':'T cell, CD4, naive [activated]','TREG_NAIVE':'T cell, CD4, naive TREG','TFH':'T cell, CD4, TFH','TH1':'T cell, CD4, TH1','THSTAR':'T cell, CD4, TH1/17','TH17':'T cell, CD4, TH17','TH2':'T cell, CD4, TH2','CD8_NAIVE':'T cell, CD8, naive','CD8_STIM':'T cell, CD8, naive [activated]'}
cols="source accession_or_url species population tissue n_reps values summary_value unit gene_id_used q1 q3 iqr frac_donors_tpm_gt1 n_donors_tpm_eq0".split()
rows=[]
for k,lab in labels.items():
    fn=os.path.join(d,k+'_TPM.csv')
    with open(fn) as f:
        r=csv.reader(f); hdr=next(r)
        hit=[row for row in r if row[0].split('.')[0]=='ENSG00000150782']
    assert len(hit)==1, (k,len(hit))
    row=hit[0]; v=np.array([float(x) for x in row[3:] if x!=''])
    q1,med,q3=np.percentile(v,[25,50,75])
    rows.append(['DICE','https://dice-database.org/download/%s_TPM.csv'%k,'Homo sapiens',lab,'peripheral blood (sorted)',len(v),';'.join('%.4g'%x for x in v),'%.3f'%med,'TPM (median across donors)',row[0]+' ('+row[2]+')','%.3f'%q1,'%.3f'%q3,'%.3f'%(q3-q1),'%.3f'%np.mean(v>1),int((v==0).sum())])
    print(lab,len(v),round(q1,2),round(med,2),round(q3,2),round(np.mean(v>1),3),int((v==0).sum()),round(v.min(),3),round(v.max(),2),file=sys.stderr)
with open(out,'w') as f:
    f.write('\t'.join(cols)+'\n')
    for r in rows: f.write('\t'.join(map(str,r))+'\n')
