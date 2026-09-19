"""第二问的区间估计：白质改变、灰质改变、以及两者之差。

供体层自助法，在每个 区室 x 队列 x 状态 的格子里有放回重抽供体，
重算统计量，取 2.5 / 97.5 百分位。
不显著不等于相同，所以图上画区间，不画星号。
"""
import numpy as np, pandas as pd

C = 0.25
NBOOT = 20000
MIN_N = 3
rng = np.random.default_rng(11)

d = pd.read_parquet('twoway.parquet')
d['stratum'] = d['区室'] + '|' + d['队列']
WM = sorted([s for s in d.stratum.unique() if s.startswith('白质')])
GM = sorted([s for s in d.stratum.unique() if s.startswith('灰质')])
CT_CN = {'oligodendrocyte': '少突胶质', 'microglia': '小胶质/巨噬', 'astrocyte': '星形胶质',
         'OPC': 'OPC', 'neuron': '神经元', 'endothelial': '内皮/血管'}
ORDER = ['oligodendrocyte', 'microglia', 'astrocyte', 'OPC', 'neuron', 'endothelial']


def lfc(ms, ct):
    return np.log2(np.median(ms) + C) - np.log2(np.median(ct) + C)


def cells(sub, strata):
    """取出每层的 (MS 值数组, 对照值数组)，供体数不足的层剔除。"""
    out = []
    for s in strata:
        g = sub[sub.stratum == s]
        ms = g[g.状态 == 'MS'].值.values
        ct = g[g.状态 == 'Control'].值.values
        if len(ms) >= MIN_N and len(ct) >= MIN_N:
            out.append((s, ms, ct))
    return out


def combine(parts, draw=False):
    """按供体数加权平均各层 log2FC；draw=True 时先重抽供体。"""
    num = w = 0.0
    for _, ms, ct in parts:
        if draw:
            ms = rng.choice(ms, len(ms), replace=True)
            ct = rng.choice(ct, len(ct), replace=True)
        n = len(ms) + len(ct)
        num += n * lfc(ms, ct); w += n
    return num / w if w else np.nan


rows = []
for key in ORDER:
    sub = d[d.细胞类型 == key]
    pw, pg = cells(sub, WM), cells(sub, GM)
    r = dict(细胞类型=CT_CN[key], 键=key,
             白质可算=len(pw) > 0, 灰质可算=len(pg) > 0,
             白质层='; '.join(f'{s.split("|")[1]}(MS {len(m)}/对照 {len(c)})' for s, m, c in pw),
             灰质层='; '.join(f'{s.split("|")[1]}(MS {len(m)}/对照 {len(c)})' for s, m, c in pg))
    dW = combine(pw) if pw else np.nan
    dG = combine(pg) if pg else np.nan
    r['ΔWM'] = dW; r['ΔGM'] = dG
    r['交互'] = dW - dG if (pw and pg) else np.nan

    bw = np.array([combine(pw, True) for _ in range(NBOOT)]) if pw else None
    bg = np.array([combine(pg, True) for _ in range(NBOOT)]) if pg else None
    if pw:
        r['ΔWM_lo'], r['ΔWM_hi'] = np.percentile(bw, [2.5, 97.5])
    if pg:
        r['ΔGM_lo'], r['ΔGM_hi'] = np.percentile(bg, [2.5, 97.5])
    if pw and pg:
        bi = bw - bg
        r['交互_lo'], r['交互_hi'] = np.percentile(bi, [2.5, 97.5])
    rows.append(r)

    # 各队列单独
    for s, ms, ct in pw + pg:
        rows.append(dict(细胞类型=CT_CN[key], 键=key, 层=s,
                         ΔWM=np.nan, ΔGM=np.nan, 交互=np.nan,
                         单层log2FC=lfc(ms, ct), nMS=len(ms), n对照=len(ct),
                         中位MS=np.median(ms), 中位对照=np.median(ct),
                         **dict(zip(['单层_lo', '单层_hi'],
                                    np.percentile([lfc(rng.choice(ms, len(ms), True),
                                                       rng.choice(ct, len(ct), True))
                                                   for _ in range(NBOOT)], [2.5, 97.5])))))

r = pd.DataFrame(rows)
r.to_csv('/home/user/aaa/results/mda5/out_q2_ci.csv', index=False)

main = r[r.层.isna()] if '层' in r.columns else r
pd.set_option('display.width', 250)
print('=== 主估计（自助 95% 区间）===')
for _, x in main.iterrows():
    def f(a, lo, hi):
        if not np.isfinite(x.get(a, np.nan)): return '      不可算      '
        return f'{x[a]:+.2f} [{x[lo]:+.2f},{x[hi]:+.2f}]'
    print(f'{x.细胞类型:<8} 白质 {f("ΔWM","ΔWM_lo","ΔWM_hi"):<22} '
          f'灰质 {f("ΔGM","ΔGM_lo","ΔGM_hi"):<22} 之差 {f("交互","交互_lo","交互_hi")}')
print()
print('=== 分层 ===')
per = r[r.层.notna()] if '层' in r.columns else pd.DataFrame()
for _, x in per.iterrows():
    print(f'{x.细胞类型:<8} {x.层:<22} log2FC {x.单层log2FC:+.2f} '
          f'[{x.单层_lo:+.2f},{x.单层_hi:+.2f}]  MS {int(x.nMS)} / 对照 {int(x.n对照)}')
