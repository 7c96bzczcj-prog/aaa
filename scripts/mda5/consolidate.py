"""把 Census 的细胞本体标签合并成十个大类，并算伪整体 CP10k。

伪整体 = Σ该类全部核的基因计数 / Σ该类全部核的总计数 × 1e4，
即按测序深度加权，不是每核比值的平均。
映射先用 IFIH1 复现已有的 out_baseline_consolidated.csv 作为校验。
"""
import numpy as np, pandas as pd, re

GROUPS = {
    'neuron': [r'neuron', r'interneuron', r'neuroblast', r'Purkinje', r'granule cell',
               r'medium spiny', r'Cajal-Retzius', r'unipolar brush'],
    'oligodendrocyte': [r'^oligodendrocyte$', r'^mature oligodendrocyte'],
    'OPC': [r'oligodendrocyte precursor', r'committed oligodendrocyte precursor'],
    'astrocyte': [r'astrocyte', r'Bergmann glial'],
    'microglia/CNS macrophage': [r'microglial', r'macrophage', r'central nervous system macrophage'],
    'endothelial': [r'endothelial'],
    'mural/perivascular': [r'pericyte', r'mural cell', r'smooth muscle',
                           r'perivascular cell', r'vascular leptomeningeal',
                           r'brain vascular cell'],
    'fibroblast': [r'fibroblast'],
    'lymphocyte/leukocyte': [r'T cell', r'B cell', r'natural killer', r'leukocyte',
                             r'lymphocyte'],
    'ependymal/choroid': [r'ependymal', r'choroid plexus'],
}
# 归入 neuron 前先排除这些（它们是胶质或前体，名字里带 glial/progenitor）
EXCLUDE_NEURON = [r'glial', r'progenitor', r'radial glia', r'glioblast', r'neural cell']


def assign(ct):
    s = str(ct)
    for g, pats in GROUPS.items():
        if g == 'neuron':
            continue
        for p in pats:
            if re.search(p, s, re.I):
                return g
    if any(re.search(p, s, re.I) for p in EXCLUDE_NEURON):
        return None
    for p in GROUPS['neuron']:
        if re.search(p, s, re.I):
            return 'neuron'
    return None


def consolidate(parquet, col, suspension='nucleus', disease='normal'):
    d = pd.read_parquet(parquet, columns=['cell_type', 'disease', 'suspension_type',
                                          'raw_sum', col])
    d = d[(d.suspension_type == suspension) & (d.disease == disease)]
    d['grp'] = d.cell_type.map(assign)
    d = d.dropna(subset=['grp'])
    g = (d.groupby('grp', observed=True)
           .agg(n=('raw_sum', 'size'), 计数=(col, 'sum'), 总计数=('raw_sum', 'sum'),
                阳性核=(col, lambda v: int((v > 0).sum())))
           .reset_index())
    g['CP10k'] = g.计数 / g.总计数 * 1e4
    g['阳性率'] = g.阳性核 / g.n
    return g.sort_values('CP10k', ascending=False).reset_index(drop=True)


if __name__ == '__main__':
    import sys
    gene = sys.argv[1] if len(sys.argv) > 1 else 'ifih1'
    pq = f'census_{gene}.parquet'
    out = consolidate(pq, gene)
    pd.set_option('display.width', 200)
    print(out.to_string(index=False, float_format=lambda v: f'{v:12.4f}'))
    if gene == 'ifih1':
        ref = pd.read_csv('/home/user/aaa/results/mda5/out_baseline_consolidated.csv')
        j = out.merge(ref[['grp', 'CP10k', 'n']], on='grp', suffixes=('_新', '_原'))
        j['比值'] = j.CP10k_新 / j.CP10k_原
        print('\n=== 与原表比对（比值应为 1）===')
        print(j[['grp', 'n_新', 'n_原', 'CP10k_新', 'CP10k_原', '比值']]
              .to_string(index=False, float_format=lambda v: f'{v:10.4f}'))
    else:
        out.to_csv(f'/home/user/aaa/results/mda5/out_{gene}_baseline.csv', index=False)
        print(f'\n写入 out_{gene}_baseline.csv')
