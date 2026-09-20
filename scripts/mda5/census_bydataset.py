"""把 Census 拆成每数据集每细胞大类的伪整体表达量。

为什么不能只报合并值：合并值按总计数加权，一个测序深、细胞多或背景高的
数据集可以独占某一行。NES 的淋巴细胞就是如此，96% 的计数来自单个数据集，
该数据集阳性率 64%，其余 17 个数据集的中位只有 0.013。
按数据集画点、取中位，这一行会自动落到应有的位置。
"""
import pandas as pd, numpy as np, sys
from consolidate import assign

MIN_CELLS = 100     # 每数据集每细胞类型的最少核数


def build(gene):
    d = pd.read_parquet(f'census_{gene}.parquet',
                        columns=['cell_type', 'disease', 'suspension_type',
                                 'raw_sum', gene, 'dataset_id'])
    d = d[(d.suspension_type == 'nucleus') & (d.disease == 'normal')]
    d['grp'] = d.cell_type.map(assign)
    d = d.dropna(subset=['grp'])
    t = (d.groupby(['grp', 'dataset_id'], observed=True)
           .agg(核数=(gene, 'size'), 计数=(gene, 'sum'), 总计数=('raw_sum', 'sum'),
                阳性核=(gene, lambda v: int((v > 0).sum())))
           .reset_index())
    t = t[t.核数 >= MIN_CELLS]
    t['CP10k'] = t.计数 / t.总计数 * 1e4
    t['阳性率'] = t.阳性核 / t.核数
    out = f'/home/user/aaa/results/mda5/out_{gene}_bydataset.csv'
    t.to_csv(out, index=False)
    s = (t.groupby('grp')
           .agg(数据集数=('CP10k', 'size'), 中位CP10k=('CP10k', 'median'),
                四分一=('CP10k', lambda v: v.quantile(.25)),
                四分三=('CP10k', lambda v: v.quantile(.75)),
                核数合计=('核数', 'sum'))
           .reset_index().sort_values('中位CP10k', ascending=False))
    s.to_csv(f'/home/user/aaa/results/mda5/out_{gene}_bydataset_summary.csv', index=False)
    print(f'--- {gene} ---')
    print(s.to_string(index=False, float_format=lambda v: f'{v:10.4f}'))
    return t


if __name__ == '__main__':
    for g in (sys.argv[1:] or ['ifih1', 'nes']):
        build(g)
        print()
