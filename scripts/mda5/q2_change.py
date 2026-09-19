"""问题二：MS 相对对照的改变，白质与灰质是否相同。

关键点：MS 与对照在白质里同属一个研究，在灰质里同属另一个研究，
因此 log2FC 是研究内的量，跨研究的系统偏移在相减时抵消。
交互项 = ΔWM - ΔGM，同样只由研究内的量构成。

统计：分层置换检验（在 区室×队列 内打乱疾病标签），不作分布假设。
"""
import numpy as np, pandas as pd

C = 0.25          # CPM 伪计数
NPERM = 50000
MIN_N = 3         # 每组最少供体数

d = pd.read_parquet('twoway.parquet')
d['stratum'] = d['区室'] + '|' + d['队列']
WM_STRATA = [s for s in d.stratum.unique() if s.startswith('白质')]
GM_STRATA = [s for s in d.stratum.unique() if s.startswith('灰质')]

def lfc(v_ms, v_ct):
    return np.log2(np.median(v_ms) + C) - np.log2(np.median(v_ct) + C)

def delta(sub, strata, is_ms):
    """分层 log2FC：各层按供体数加权平均。层内任一组不足 MIN_N 则该层不计。"""
    num = 0.0; wsum = 0.0; used = []
    for s in strata:
        m = sub.stratum.values == s
        if not m.any():
            continue
        y = sub.值.values[m]; ms = is_ms[m]
        n1, n0 = int(ms.sum()), int((~ms).sum())
        if n1 < MIN_N or n0 < MIN_N:
            continue
        w = n1 + n0
        num += w * lfc(y[ms], y[~ms]); wsum += w
        used.append((s, n1, n0))
    if wsum == 0:
        return np.nan, used
    return num / wsum, used

rows = []
for ct, sub in d.groupby('细胞类型'):
    sub = sub.reset_index(drop=True)
    is_ms = (sub.状态.values == 'MS')

    dWM, uWM = delta(sub, WM_STRATA, is_ms)
    dGM, uGM = delta(sub, GM_STRATA, is_ms)
    if not np.isfinite(dWM) or not np.isfinite(dGM):
        rows.append(dict(细胞类型=ct, ΔWM=dWM, ΔGM=dGM, 交互=np.nan,
                         p_WM=np.nan, p_GM=np.nan, p_交互=np.nan,
                         n白质=str(uWM), n灰质=str(uGM)))
        continue
    inter = dWM - dGM

    # 分层置换：在每个 区室×队列 层内独立打乱疾病标签
    rng = np.random.default_rng(20240919)
    strat_idx = [np.where(sub.stratum.values == s)[0] for s in sub.stratum.unique()]
    cWM = cGM = cIN = 0
    for _ in range(NPERM):
        perm = is_ms.copy()
        for idx in strat_idx:
            perm[idx] = rng.permutation(perm[idx])
        pWM, _ = delta(sub, WM_STRATA, perm)
        pGM, _ = delta(sub, GM_STRATA, perm)
        if np.isfinite(pWM) and abs(pWM) >= abs(dWM) - 1e-12: cWM += 1
        if np.isfinite(pGM) and abs(pGM) >= abs(dGM) - 1e-12: cGM += 1
        if np.isfinite(pWM) and np.isfinite(pGM) and abs(pWM - pGM) >= abs(inter) - 1e-12: cIN += 1

    rows.append(dict(细胞类型=ct, ΔWM=dWM, ΔGM=dGM, 交互=inter,
                     p_WM=(cWM + 1) / (NPERM + 1),
                     p_GM=(cGM + 1) / (NPERM + 1),
                     p_交互=(cIN + 1) / (NPERM + 1),
                     n白质=str(uWM), n灰质=str(uGM)))

r = pd.DataFrame(rows)

def bh(p):
    p = np.asarray(p, float); ok = np.isfinite(p); q = np.full(p.shape, np.nan)
    pv = p[ok]; n = pv.size
    if n == 0: return q
    o = np.argsort(pv); s = pv[o] * n / (np.arange(n) + 1)
    s = np.minimum.accumulate(s[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(s, 1)
    q[ok] = out; return q

for c in ['WM', 'GM', '交互']:
    r['q_' + c] = bh(r['p_' + c])

order = ['少突胶质', '小胶质/巨噬', '星形胶质', 'OPC', '神经元', '内皮/血管']
r['_o'] = r.细胞类型.map({k: i for i, k in enumerate(order)})
r = r.sort_values('_o').drop(columns='_o')
r.to_csv('/home/user/aaa/results/mda5/out_q2_change_interaction.csv', index=False)
pd.set_option('display.width', 220)
print(r[['细胞类型', 'ΔWM', 'p_WM', 'q_WM', 'ΔGM', 'p_GM', 'q_GM', '交互', 'p_交互', 'q_交互']].to_string(
    index=False, float_format=lambda x: f'{x:7.4f}'))
print()
print('白质层：', [s for s in WM_STRATA], ' 灰质层：', [s for s in GM_STRATA])

# 每个队列单独的 log2FC，供图与稳健性检查
per = []
for ct, sub in d.groupby('细胞类型'):
    for s, g in sub.groupby('stratum'):
        ms = g[g.状态 == 'MS'].值.values; ct0 = g[g.状态 == 'Control'].值.values
        if len(ms) < MIN_N or len(ct0) < MIN_N:
            per.append(dict(细胞类型=ct, 层=s, log2FC=np.nan, nMS=len(ms), n对照=len(ct0),
                            中位MS=np.median(ms) if len(ms) else np.nan,
                            中位对照=np.median(ct0) if len(ct0) else np.nan))
            continue
        per.append(dict(细胞类型=ct, 层=s, log2FC=lfc(ms, ct0), nMS=len(ms), n对照=len(ct0),
                        中位MS=np.median(ms), 中位对照=np.median(ct0)))
p = pd.DataFrame(per)
p['_o'] = p.细胞类型.map({k: i for i, k in enumerate(order)})
p = p.sort_values(['_o', '层']).drop(columns='_o')
p.to_csv('/home/user/aaa/results/mda5/out_q2_change_bycohort.csv', index=False)
print()
print(p.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
