"""统一口径的第一问与第二问统计。

量：深度匹配区间内每供体的 IFIH1 检出率。
效应：log2（MS 检出率 / 对照检出率），供体自助 95% 区间。
不用伪计数，不用跨队列标度，因此不受 Schirmer 归一化尺度未知的影响。
"""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from detect import nuclei, donor_rates

NBOOT = 20000
EPS = 1e-4                    # 检出率的下限，避免零率取对数；远小于任何可观测率
MIN_D = 3
rng = np.random.default_rng(5)

ORDER = ['oligodendrocyte', 'microglia', 'astrocyte', 'OPC', 'neuron', 'endothelial']
import figstyle as _F
CN = {k: _F.CT_CN[k] for k in ORDER}
WM_COH = ['GSE180759', 'GSE279180']
GM_COH = ['Schirmer 2019']

g = donor_rates(nuclei())


def rates(ct, coh, st):
    x = g[(g.细胞类型 == ct) & (g.队列 == coh) & (g.状态 == st)]
    return x.检出率.values, x.供体.values


def lr(ms, ct):
    return np.log2(np.median(ms) + EPS) - np.log2(np.median(ct) + EPS)


def effect(ct, cohorts, draw=False):
    """按供体数加权合并各队列的 log2 比值。"""
    num = w = 0.0
    for coh in cohorts:
        m, _ = rates(ct, coh, 'MS'); c, _ = rates(ct, coh, 'Control')
        if len(m) < MIN_D or len(c) < MIN_D:
            continue
        if draw:
            m = rng.choice(m, len(m), True); c = rng.choice(c, len(c), True)
        n = len(m) + len(c)
        num += n * lr(m, c); w += n
    return num / w if w else np.nan


rows = []
for ct in ORDER:
    r = {'细胞类型': CN[ct], '键': ct}
    for tag, cohs in (('WM', WM_COH), ('GM', GM_COH)):
        e = effect(ct, cohs)
        r['Δ' + tag] = e
        if np.isfinite(e):
            b = np.array([effect(ct, cohs, True) for _ in range(NBOOT)])
            b = b[np.isfinite(b)]
            r['Δ%s_lo' % tag], r['Δ%s_hi' % tag] = np.percentile(b, [2.5, 97.5])
            # 供体级 Mann-Whitney（多队列时按队列合并秩，用分层近似：取各队列 p 的 Fisher 合并）
            ps = []
            for coh in cohs:
                m, _ = rates(ct, coh, 'MS'); c, _ = rates(ct, coh, 'Control')
                if len(m) >= MIN_D and len(c) >= MIN_D:
                    ps.append(mannwhitneyu(m, c, alternative='two-sided').pvalue)
            if ps:
                chi = -2 * np.sum(np.log(ps))
                from scipy.stats import chi2
                r['p_' + tag] = chi2.sf(chi, 2 * len(ps))
    if np.isfinite(r.get('ΔWM', np.nan)) and np.isfinite(r.get('ΔGM', np.nan)):
        r['交互'] = r['ΔWM'] - r['ΔGM']
        bi = np.array([effect(ct, WM_COH, True) - effect(ct, GM_COH, True)
                       for _ in range(NBOOT)])
        bi = bi[np.isfinite(bi)]
        r['交互_lo'], r['交互_hi'] = np.percentile(bi, [2.5, 97.5])
    rows.append(r)

R = pd.DataFrame(rows)
R.to_csv('/home/user/aaa/results/mda5/out_detect_effects.csv', index=False)

# 分队列明细
per = []
for ct in ORDER:
    for coh in WM_COH + GM_COH:
        m, _ = rates(ct, coh, 'MS'); c, _ = rates(ct, coh, 'Control')
        d = dict(细胞类型=CN[ct], 键=ct, 队列=coh, nMS=len(m), n对照=len(c),
                 率MS=np.median(m) if len(m) else np.nan,
                 率对照=np.median(c) if len(c) else np.nan)
        if len(m) >= MIN_D and len(c) >= MIN_D:
            d['log2比值'] = lr(m, c)
            b = np.array([lr(rng.choice(m, len(m), True), rng.choice(c, len(c), True))
                          for _ in range(NBOOT)])
            d['lo'], d['hi'] = np.percentile(b, [2.5, 97.5])
            d['p'] = mannwhitneyu(m, c, alternative='two-sided').pvalue
        per.append(d)
P = pd.DataFrame(per)
P.to_csv('/home/user/aaa/results/mda5/out_detect_bycohort.csv', index=False)

pd.set_option('display.width', 220)
print('=== 主估计：log2（MS 检出率 / 对照检出率），自助 95% 区间 ===')
for _, x in R.iterrows():
    def f(a):
        if not np.isfinite(x.get(a, np.nan)): return '        不可算        '
        return f"{x[a]:+.2f} [{x[a+'_lo']:+.2f},{x[a+'_hi']:+.2f}]"
    print(f"{x.细胞类型:<8} 白质 {f('ΔWM'):<24} 皮层块 {f('ΔGM'):<24} 之差 {f('交互')}")
print()
print('=== 分队列 ===')
print(P.to_string(index=False, float_format=lambda v: f'{v:8.4f}'))
