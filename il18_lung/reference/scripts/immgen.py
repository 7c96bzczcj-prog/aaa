import sys, re, csv, gzip
import numpy as np

dl = sys.argv[1]; out = sys.argv[2]

def load(fn, sep):
    op = gzip.open if fn.endswith('.gz') else open
    with op(fn, 'rt') as f:
        r = csv.reader(f, delimiter=sep); hdr = next(r)
        rows = {}; want = {'Il18', 'Actb', 'Ptprc'}
        tot = np.zeros(len(hdr) - 1)
        for row in r:
            v = np.array([float(x) if x not in ('', 'NA') else np.nan for x in row[1:]])
            tot += np.nan_to_num(v)
            if row[0] in want: rows[row[0]] = v
    return hdr[1:], rows, tot

sets = {
    'GSE109125': (dl + '/GSE109125_Normalized_Gene_count_table.csv', dl + '/GSE109125_Genes_count_table.tsv.gz', '\t', lambda s: s.rsplit('#', 1)[0]),
    'GSE122108': (dl + '/GSE122108_Normalized_Gene_count_table.csv', dl + '/GSE122108_Gene_count_table.csv.gz', '\t', lambda s: re.sub(r'\.\d+$', '', s)),
}
pick = {
    'GSE109125': ['MF.Alv.Lu', 'MF.pIC.Alv.Lu', 'MAIT.Lu', 'NKT.19-8-TCRb+CD1daGalCerTet+.Lu', 'Baso.Nb.B7hi.Lung', 'MC.Nb.B7hi.Lung',
                  'Mo.6C+II-.Bl', 'Mo.6C-II-.Bl', 'GN.BM', 'GN.Sp', 'Eo.Sp', 'NK.27+11b-.Sp', 'NK.27+11b+.Sp', 'NK.27-11b+.Sp',
                  'T.4.Nve.Sp', 'T.8.Nve.Sp', 'Treg.4.25hi.Sp', 'Tgd.Sp', 'B.Fo.Sp', 'B.MZ.Sp', 'DC.4+.Sp', 'DC.8+.Sp', 'DC.pDC.Sp',
                  'MF.RP.Sp', 'MF.PC', 'BEC.SLN', 'LEC.SLN', 'FRC.CD140a+.Madcam-.CD35-.SLN', 'IAP.SLN', 'Ep.MECHi.Th'],
    'GSE122108': None,  # all lung / BAL / lung-LN populations + blood monocytes
}
tissue_map = [('BAL', 'lung (bronchoalveolar lavage)'), ('LuLN', 'lung-draining lymph node'), ('Lung', 'lung'), ('Lu', 'lung'),
              ('Neo', 'lung (neonatal)'), ('Bl', 'blood'), ('Sp', 'spleen'), ('BM', 'bone marrow'), ('SLN', 'skin-draining lymph node'),
              ('Th', 'thymus'), ('PC', 'peritoneal cavity')]

def tissue(p):
    last = p.split('.')[-1]
    for k, v in tissue_map:
        if last == k: return v
    return last

def note(p):
    n = []
    if 'LPS' in p: n.append('LPS-challenged')
    if 'pIC' in p: n.append('poly(I:C)-challenged')
    if re.search(r'E1\d\.5', p): n.append('embryonic')
    if 'Neo' in p: n.append('neonatal')
    if re.search(r'^(BEC|LEC|FRC|IAP|Ep)\.', p): n.append('non-immune stromal/epithelial, NOT lung')
    return '; '.join(n) if n else 'baseline'

cols = ("source accession_or_url species population tissue n_reps values summary_value unit gene_id_used "
        "replicate_ids raw_counts_values raw_cpm_values note").split()
out_rows = []; sig = {}
for g, (nf, rf, sep, grp) in sets.items():
    nh, nr, _ = load(nf, ',')
    rh, rr, rtot = load(rf, sep)
    assert nh == rh, (g, 'header mismatch')
    nv = nr['Il18']; rv = rr['Il18']
    # verify normalized = raw * per-sample factor + 1
    k = (nr['Actb'] - 1) / rr['Actb']; k2 = (nr['Ptprc'] - 1) / rr['Ptprc']
    print(g, 'max rel diff of per-sample factor Actb vs Ptprc:', np.nanmax(np.abs(k - k2) / k), file=sys.stderr)
    pops = {}
    for i, s in enumerate(nh): pops.setdefault(grp(s), []).append(i)
    sel = pick[g] if pick[g] is not None else [p for p in pops if re.search(r'\.(Lu|BAL|LuLN)(\.Neo)?$', p) or p in ('Mo.6Cp.Bl', 'Mo.6Cn.Bl')]
    for p in sel:
        idx = pops[p]
        for i in idx:
            sig.setdefault((rr['Ptprc'][i], rv[i]), set()).add(g + ':' + p)
        vals = nv[idx]; raw = rv[idx]; cpm = raw / rtot[idx] * 1e6
        out_rows.append(['ImmGen ULI RNA-seq (Smart-seq2, ~1000 sorted cells)',
                         'https://sharehost.hms.harvard.edu/immgen/%s/%s_Normalized_Gene_count_table.csv (GEO %s)' % (g, g, g),
                         'Mus musculus', p, tissue(p), len(idx), ';'.join('%.2f' % x for x in vals), '%.2f' % np.mean(vals),
                         'ImmGen normalized count = raw reads x per-sample scale factor + 1 (verified from raw table); 1.00 = 0 raw reads; summary = mean',
                         'symbol Il18 (tables carry symbols only, no Ensembl IDs)', ';'.join(nh[i] for i in idx),
                         ';'.join('%d' % x for x in raw), ';'.join('%.2f' % x for x in cpm), note(p)])

for r in out_rows:
    me = r[1].split('GEO ')[1].rstrip(')') + ':' + r[3]
    dup = sorted(set(x for v in sig.values() if me in v for x in v if x != me))
    if dup:
        r[13] += '; replicate(s) identical (identical Ptprc and Il18 raw counts, Actb within 1-2 reads; GEO annotations differ, vM12 vs vM16) to ' + ','.join(dup) + ' = same samples deposited in both series'

with open(out, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in out_rows: f.write('\t'.join(map(str, r)) + '\n')
