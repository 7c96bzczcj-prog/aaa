"""以 Macnair 2025 为主队列重算两个问题。

为什么换主队列：此前皮层与白质分属不同研究，区域与研究完全混同，
任何在 log2 比值内部的推理都无法把两者分开。Macnair 在同一研究、
同一建库与分析流程下同时采集皮层与皮层下白质，80 位供体，63 万个核，
因此区域对比第一次成为研究内的量。

统计：供体层 CPM，log2 比值用中位，95% 区间用供体自助，
p 用供体标签置换（区域对比在状态内置换，疾病对比在区域内置换）。
"""
import numpy as np, pandas as pd

NBOOT = 20000
NPERM = 20000
MIN_NUC = 30
MIN_D = 5
EPS = 0.05          # CPM 下限，远小于任何可观测值
rng = np.random.default_rng(17)

CT = {'Oligodendrocytes': '少突胶质细胞', 'Microglia': '小胶质细胞',
      'Astrocytes': '星形胶质细胞', 'OPCs + COPs': '少突胶质前体细胞',
      'Excitatory neurons': '兴奋性神经元', 'Inhibitory neurons': '抑制性神经元',
      'Endo + Peri': '内皮/周细胞'}
ORDER = list(CT)

d = pd.read_parquet('macnair_twoway.parquet')
d = d[(d.ct.isin(ORDER)) & (d.nuclei >= MIN_NUC)].copy()


def vals(ct, matter, status):
    x = d[(d.ct == ct) & (d.matter == matter) & (d.status == status)]
    return x.cpm.values


def lr(a, b):
    return np.log2(np.median(a) + EPS) - np.log2(np.median(b) + EPS)


def boot(a, b, n=NBOOT):
    out = np.empty(n)
    for i in range(n):
        out[i] = lr(rng.choice(a, len(a), True), rng.choice(b, len(b), True))
    return out


def perm_p(a, b, obs, n=NPERM):
    """两组合并后随机划分，看 |统计量| 是否常常达到观测值。"""
    pool = np.concatenate([a, b]); k = len(a); c = 0
    for _ in range(n):
        p = rng.permutation(pool)
        if abs(lr(p[:k], p[k:])) >= abs(obs) - 1e-12:
            c += 1
    return (c + 1) / (n + 1)


rows = []
for key in ORDER:
    nm = CT[key]
    r = {'细胞类型': nm, '键': key}

    # 问题一：正常供体，白质 vs 灰质
    w0, g0 = vals(key, 'WM', 'Control'), vals(key, 'GM', 'Control')
    if len(w0) >= MIN_D and len(g0) >= MIN_D:
        e = lr(w0, g0); b = boot(w0, g0)
        r.update({'正常_白比灰': e, '正常_lo': np.percentile(b, 2.5),
                  '正常_hi': np.percentile(b, 97.5), '正常_p': perm_p(w0, g0, e),
                  'n正常白': len(w0), 'n正常灰': len(g0)})

    # 问题二：各区室内部 MS vs 对照
    for mt, tag in (('WM', '白质'), ('GM', '灰质')):
        c0, c1 = vals(key, mt, 'Control'), vals(key, mt, 'MS')
        if len(c0) >= MIN_D and len(c1) >= MIN_D:
            e = lr(c1, c0); b = boot(c1, c0)
            r.update({f'Δ{tag}': e, f'Δ{tag}_lo': np.percentile(b, 2.5),
                      f'Δ{tag}_hi': np.percentile(b, 97.5),
                      f'p{tag}': perm_p(c1, c0, e),
                      f'n{tag}对照': len(c0), f'n{tag}MS': len(c1)})

    # 交互：白质的改变 减 灰质的改变
    if np.isfinite(r.get('Δ白质', np.nan)) and np.isfinite(r.get('Δ灰质', np.nan)):
        r['交互'] = r['Δ白质'] - r['Δ灰质']
        w1, wc = vals(key, 'WM', 'MS'), vals(key, 'WM', 'Control')
        g1, gc = vals(key, 'GM', 'MS'), vals(key, 'GM', 'Control')
        bi = np.empty(NBOOT)
        for i in range(NBOOT):
            bi[i] = (lr(rng.choice(w1, len(w1), True), rng.choice(wc, len(wc), True))
                     - lr(rng.choice(g1, len(g1), True), rng.choice(gc, len(gc), True)))
        r['交互_lo'], r['交互_hi'] = np.percentile(bi, [2.5, 97.5])
        # 置换：在各区室内部独立打乱疾病标签
        obs = r['交互']; c = 0
        pw = np.concatenate([w1, wc]); kw = len(w1)
        pg = np.concatenate([g1, gc]); kg = len(g1)
        for _ in range(NPERM):
            a = rng.permutation(pw); b2 = rng.permutation(pg)
            if abs(lr(a[:kw], a[kw:]) - lr(b2[:kg], b2[kg:])) >= abs(obs) - 1e-12:
                c += 1
        r['p交互'] = (c + 1) / (NPERM + 1)
    rows.append(r)

R = pd.DataFrame(rows)


def bh(p):
    p = np.asarray(p, float); ok = np.isfinite(p); q = np.full(p.shape, np.nan)
    pv = p[ok]; n = pv.size
    if n == 0: return q
    o = np.argsort(pv); s = pv[o] * n / (np.arange(n) + 1)
    s = np.minimum.accumulate(s[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(s, 1); q[ok] = out
    return q


for c in ['正常_p', 'p白质', 'p灰质', 'p交互']:
    if c in R:
        R['q' + c.lstrip('p').lstrip('_')] = bh(R[c])
R.to_csv('/home/user/aaa/results/mda5/out_macnair_effects.csv', index=False)

pd.set_option('display.width', 250)


def f(x, a):
    if not np.isfinite(x.get(a, np.nan)):
        return '        —        '
    lo, hi = (a + '_lo', a + '_hi') if a.startswith('Δ') or a == '交互' else ('正常_lo', '正常_hi')
    return f'{x[a]:+.2f} [{x[lo]:+.2f},{x[hi]:+.2f}]'


print('=== 问题一：正常供体，白质 / 灰质 ===')
for _, x in R.iterrows():
    if np.isfinite(x.get('正常_白比灰', np.nan)):
        print(f"{x.细胞类型:<10} {f(x,'正常_白比灰'):<24} p={x['正常_p']:.3f}  "
              f"n={int(x.n正常白)}/{int(x.n正常灰)}")
print()
print('=== 问题二：MS 相对对照，及两区室之差 ===')
for _, x in R.iterrows():
    print(f"{x.细胞类型:<10} 白质 {f(x,'Δ白质'):<24} 灰质 {f(x,'Δ灰质'):<24} "
          f"之差 {f(x,'交互'):<24} p交互={x.get('p交互', float('nan')):.3f}")
