"""把三个队列的细胞级 IFIH1 统一到 CP10k。

GSE180759 与 GSE279180 存的是原始计数，直接除以该核总计数。
Schirmer 存的是 Seurat 对数归一化值，expm1 后还差一个每细胞尺度因子 S。
S 在同一组内是常数，可由该组的供体级 CP10k 反解：
    CP10k_组 = (Σ_cells expm1) / (n_cells · S) · 1e4
    ⇒ S = (Σ_cells expm1) · 1e4 / (n_cells · CP10k_组)
因此不需要重读表达矩阵，也不含任何拟合。
"""
import pandas as pd, numpy as np

CT = {'oligodendrocytes': 'oligodendrocyte', 'astrocytes': 'astrocyte', 'immune': 'microglia',
      'neurons': 'neuron', 'opc': 'OPC', 'vascular_cells': 'endothelial',
      'lymphocytes': 'lymphocyte', 'microglia': 'microglia', 'OPC': 'OPC',
      'endothelial/vascular': 'endothelial', 'T cells': 'lymphocyte',
      'B cells': 'lymphocyte', 'stromal/other': 'stromal'}


def load_cells():
    out = []

    d = pd.read_parquet('df_180759.parquet')
    out.append(pd.DataFrame({
        '区室': '白质', '队列': 'GSE180759',
        '细胞类型': d.ct.map(lambda x: CT.get(x, x)).values,
        '供体': d.donor.astype(str).values,
        '状态': np.where(d.path.values == 'control_white_matter', 'Control', 'MS'),
        'CP10k': d.IFIH1.values / np.maximum(d.total.values, 1) * 1e4}))

    d = pd.read_parquet('lerma_ifih1.parquet')
    out.append(pd.DataFrame({
        '区室': '白质', '队列': 'GSE279180',
        '细胞类型': d.ctype.map(lambda x: CT.get(x, x)).values,
        '供体': d.patient_id.astype(str).values,
        '状态': np.where(d.condition.values == 'Control', 'Control', 'MS'),
        'CP10k': d.IFIH1.values / np.maximum(d.total.values, 1) * 1e4}))

    s = pd.read_parquet('schirmer_ifih1.parquet')
    raw = np.expm1(s.IFIH1.values.astype(np.float64))
    g = pd.DataFrame({'细胞类型': s.grp.map(lambda x: CT.get(x, x)).values,
                      '供体': s['sample'].astype(str).values,
                      '状态': np.where(s.cond.values == 'Control', 'Control', 'MS'),
                      'raw': raw})
    tw = pd.read_parquet('twoway.parquet')
    tw = tw[tw.区室 == '灰质'][['细胞类型', '供体', '值']].rename(columns={'值': 'CPM组'})
    agg = g.groupby(['细胞类型', '供体']).agg(S_raw=('raw', 'sum'), n=('raw', 'size')).reset_index()
    agg = agg.merge(tw, on=['细胞类型', '供体'], how='left')
    # CP10k = CPM / 100
    agg['S'] = agg.S_raw * 1e4 / (agg.n * (agg['CPM组'] / 100.0))
    g = g.merge(agg[['细胞类型', '供体', 'S']], on=['细胞类型', '供体'], how='left')
    g['CP10k'] = g.raw / g.S * 1e4
    g = g[np.isfinite(g.CP10k)]
    out.append(pd.DataFrame({'区室': '灰质', '队列': 'Schirmer 2019',
                             '细胞类型': g.细胞类型.values, '供体': g.供体.values,
                             '状态': g.状态.values, 'CP10k': g.CP10k.values}))

    return pd.concat(out, ignore_index=True)


if __name__ == '__main__':
    c = load_cells()
    print(c.groupby(['区室', '队列']).size().to_string())
    print()
    # 校验：细胞级均值应当复现供体级 CP10k
    tw = pd.read_parquet('twoway.parquet')
    m = c.groupby(['区室', '队列', '细胞类型', '供体']).CP10k.mean().reset_index()
    j = m.merge(tw, on=['区室', '队列', '细胞类型', '供体'], suffixes=('_细胞均值', ''))
    j['供体CP10k'] = j['值'] / 100.0
    j['比值'] = j.CP10k / j['供体CP10k']
    print('细胞级均值 / 供体级 CP10k 的比值（应接近 1）：')
    print(j.groupby(['区室', '队列'])['比值'].describe()[['count', '50%', 'min', 'max']].to_string())
