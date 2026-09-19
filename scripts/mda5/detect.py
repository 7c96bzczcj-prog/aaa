"""统一口径重算：每个核是否检出 IFIH1。

为什么换这个量：
  白质两队列存的是原始计数，灰质（Schirmer, UCSC）只有 Seurat 对数归一化值，
  尺度因子无法反解，两侧的伪整体 CPM 因此不在同一标度上。
  "该核是否检出 IFIH1" 在任何单调归一化下都不变，三个队列精确可比。

检出率随测序深度升高，因此只在三个队列共有的深度区间内比较，不用模型。
"""
import numpy as np, pandas as pd

BAND = (1500, 5000)          # 每核总 UMI 的共同区间
MIN_NUC = 30                 # 区间内每供体每细胞类型的最少核数

CT = {'oligodendrocytes': 'oligodendrocyte', 'astrocytes': 'astrocyte', 'immune': 'microglia',
      'neurons': 'neuron', 'opc': 'OPC', 'vascular_cells': 'endothelial',
      'lymphocytes': 'lymphocyte', 'microglia': 'microglia', 'OPC': 'OPC',
      'endothelial/vascular': 'endothelial', 'T cells': 'lymphocyte',
      'B cells': 'lymphocyte', 'stromal/other': 'stromal'}


def nuclei():
    out = []
    d = pd.read_parquet('df_180759.parquet')
    out.append(pd.DataFrame({
        '区室': '皮层下白质', '队列': 'GSE180759',
        '细胞类型': d.ct.map(lambda x: CT.get(x, x)).values,
        '供体': d.donor.astype(str).values,
        '状态': np.where(d.path.values == 'control_white_matter', 'Control', 'MS'),
        'UMI': d.total.values.astype(float),
        '检出': (d.IFIH1.values > 0).astype(int)}))

    d = pd.read_parquet('lerma_ifih1.parquet')
    out.append(pd.DataFrame({
        '区室': '皮层下白质', '队列': 'GSE279180',
        '细胞类型': d.ctype.map(lambda x: CT.get(x, x)).values,
        '供体': d.patient_id.astype(str).values,
        '状态': np.where(d.condition.values == 'Control', 'Control', 'MS'),
        'UMI': d.total.values.astype(float),
        '检出': (d.IFIH1.values > 0).astype(int)}))

    s = pd.read_parquet('schirmer_ifih1.parquet')
    out.append(pd.DataFrame({
        '区室': '皮层组织块', '队列': 'Schirmer 2019',
        '细胞类型': s.grp.map(lambda x: CT.get(x, x)).values,
        '供体': s['sample'].astype(str).values,
        '状态': np.where(s.cond.values == 'Control', 'Control', 'MS'),
        'UMI': s.UMIs.values.astype(float),
        '检出': (s.IFIH1.values > 0).astype(int)}))
    return pd.concat(out, ignore_index=True)


def donor_rates(n, band=BAND, min_nuc=MIN_NUC):
    """深度区间内，每供体每细胞类型的检出率。"""
    b = n[(n.UMI >= band[0]) & (n.UMI <= band[1])]
    g = (b.groupby(['区室', '队列', '细胞类型', '供体', '状态'])
          .agg(核数=('检出', 'size'), 阳性=('检出', 'sum'), UMI中位=('UMI', 'median'))
          .reset_index())
    g['检出率'] = g.阳性 / g.核数
    return g[g.核数 >= min_nuc].reset_index(drop=True)


def wm_content(n):
    """Schirmer 每个供体组织块中少突胶质占比，作为白质混入量的代理。"""
    s = n[n.队列 == 'Schirmer 2019']
    tot = s.groupby('供体').size().rename('全部')
    oli = s[s.细胞类型 == 'oligodendrocyte'].groupby('供体').size().rename('少突')
    w = pd.concat([tot, oli], axis=1).fillna(0)
    w['白质代理'] = w.少突 / w.全部
    st = s.groupby('供体').状态.first()
    return w.join(st).reset_index()


if __name__ == '__main__':
    n = nuclei()
    n.to_parquet('nuclei_detect.parquet')
    print('核数合计', len(n))
    print(n.groupby(['区室', '队列']).agg(核=('检出', 'size'), UMI中位=('UMI', 'median')).to_string())
    b = n[(n.UMI >= BAND[0]) & (n.UMI <= BAND[1])]
    print(f'\n深度区间 {BAND[0]}–{BAND[1]} 内保留 {len(b)} 核（{100*len(b)/len(n):.0f}%）')
    print(b.groupby(['队列']).agg(核=('检出', 'size'), UMI中位=('UMI', 'median')).to_string())

    g = donor_rates(n)
    g.to_csv('/home/user/aaa/results/mda5/out_detect_donor.csv', index=False)
    print('\n每供体记录数', len(g))
    piv = (g.groupby(['细胞类型', '队列', '状态'])
            .agg(供体=('供体', 'nunique'), 检出率中位=('检出率', 'median'))
            .reset_index())
    print(piv.to_string(index=False))

    w = wm_content(n)
    w.to_csv('/home/user/aaa/results/mda5/out_schirmer_wm_content.csv', index=False)
    print('\nSchirmer 组织块中少突胶质占比（白质混入代理）：')
    print(w.groupby('状态')['白质代理'].describe()[['count', '25%', '50%', '75%']].to_string())
    from scipy.stats import mannwhitneyu
    a = w[w.状态 == 'MS'].白质代理.values; c = w[w.状态 == 'Control'].白质代理.values
    print('MS vs 对照  Mann-Whitney p = %.4f' % mannwhitneyu(a, c).pvalue)
