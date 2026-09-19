"""按供体 × 病灶分期 × 细胞类型汇总 Macnair 2025 的 IFIH1。

GSE180759 只有 3 位对照，按分期逐个比时双侧最小可能 p 为 0.057，
单侧为 0.029，任何星号都会正好卡在这个下限上，不构成证据。
Macnair 白质有 13 位对照、34 位 MS 供体，且病灶分期字段完整
（WM 对照白质，NAWM 正常表现白质，AL 活动性，CAL 慢性活动性，
CIL 慢性非活动性，RL 再髓鞘化），因此按分期比较才有统计力。
"""
import gzip, numpy as np, pandas as pd

COLS = ['individual_id_anon', 'type_broad', 'lesion_type', 'matter',
        'diagnosis', 'log10_counts', 'exclude_pseudobulk']

meta = pd.read_csv('geo/coldata.txt.gz', usecols=COLS)
print('元数据', meta.shape)

iv = pd.read_csv('geo/ifih1_macnair.tsv', sep='\t', header=None,
                 names=['idx', 'count'])
print('IFIH1 非零核', len(iv), '索引范围', iv.idx.min(), iv.idx.max())

# 判定索引基数：R 导出通常为 1 基
base = 1 if iv.idx.min() >= 1 and iv.idx.max() <= len(meta) else 0
pos = iv.idx.values - base
assert pos.min() >= 0 and pos.max() < len(meta), (pos.min(), pos.max())
print('索引基数', base)

cnt = np.zeros(len(meta), np.int32)
cnt[pos] = iv['count'].values
meta['ifih1'] = cnt
meta['umi'] = np.power(10.0, meta.log10_counts.values)

keep = meta.exclude_pseudobulk.astype(str).str.upper() != 'TRUE'
meta = meta[keep]
print('排除混合簇后', len(meta))

g = (meta.groupby(['individual_id_anon', 'matter', 'lesion_type',
                   'type_broad', 'diagnosis'], observed=True)
         .agg(核数=('ifih1', 'size'), 计数=('ifih1', 'sum'),
              UMI=('umi', 'sum'), 阳性核=('ifih1', lambda v: int((v > 0).sum())))
         .reset_index())
g['CPM'] = g.计数 / g.UMI * 1e6
g['阳性率'] = g.阳性核 / g.核数
g['状态'] = np.where(g.diagnosis.values == 'CTR', 'Control', 'MS')
g.to_parquet('macnair_lesion.parquet')
g.to_csv('/home/user/aaa/results/mda5/out_macnair_lesion.csv', index=False)

print()
print('白质各分期的供体数（核数 >= 30）：')
w = g[(g.matter == 'WM') & (g.核数 >= 30)]
print(w.pivot_table(index='type_broad', columns='lesion_type',
                    values='individual_id_anon', aggfunc='nunique').to_string())
print()
print('灰质各分期的供体数（核数 >= 30）：')
m = g[(g.matter == 'GM') & (g.核数 >= 30)]
print(m.pivot_table(index='type_broad', columns='lesion_type',
                    values='individual_id_anon', aggfunc='nunique').to_string())
