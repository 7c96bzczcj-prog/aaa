"""按病灶分期，每一期各自与本区室的对照比较（Macnair 2025）。

白质：对照白质 → 正常表现白质 → 活动性 → 慢性活动性 → 慢性非活动性 → 再髓鞘化
灰质：对照灰质 → 正常表现灰质 → 皮层病灶
每一期对对照做双侧 Mann-Whitney，区室内按细胞类型 BH 校正。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import figstyle as F

MIN_NUC = 30
g = pd.read_parquet('macnair_lesion.parquet')
g = g[g.核数 >= MIN_NUC]

WM = [('WM', '对照白质'), ('NAWM', '正常表现白质'), ('AL', '活动性病灶'),
      ('CAL', '慢性活动性病灶'), ('CIL', '慢性非活动性病灶'), ('RL', '再髓鞘化病灶')]
GM = [('GM', '对照灰质'), ('NAGM', '正常表现灰质'), ('GML', '皮层病灶')]
CTS = [('Oligodendrocytes', '少突胶质细胞'), ('Microglia', '小胶质细胞'),
       ('Astrocytes', '星形胶质细胞'), ('OPCs + COPs', '少突胶质前体细胞')]


def vals(ct, matter, lt):
    x = g[(g.type_broad == ct) & (g.matter == matter) & (g.lesion_type == lt)]
    return x.CPM.values


def bh(p):
    p = np.asarray(p, float); ok = np.isfinite(p); q = np.full(p.shape, np.nan)
    pv = p[ok]; n = pv.size
    if n == 0: return q
    o = np.argsort(pv); s = pv[o] * n / (np.arange(n) + 1)
    s = np.minimum.accumulate(s[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(s, 1); q[ok] = out
    return q


rows = []
for matter, stages in (('WM', WM), ('GM', GM)):
    for key, nm in CTS:
        ref = vals(key, matter, stages[0][0])
        for lt, lab in stages[1:]:
            v = vals(key, matter, lt)
            p = (mannwhitneyu(v, ref, alternative='two-sided').pvalue
                 if len(v) >= 3 and len(ref) >= 3 else np.nan)
            rows.append(dict(区室=matter, 细胞类型=nm, 分期=lt, 标签=lab.replace('\n', ''),
                             n=len(v), n对照=len(ref),
                             中位=np.median(v) if len(v) else np.nan,
                             中位对照=np.median(ref) if len(ref) else np.nan, p=p))
S = pd.DataFrame(rows)
S['q'] = np.nan
for (mt, ct), idx in S.groupby(['区室', '细胞类型']).groups.items():
    S.loc[idx, 'q'] = bh(S.loc[idx, 'p'].values)
S.to_csv('/home/user/aaa/results/mda5/out_macnair_lesion_tests.csv', index=False)

fig, axes = plt.subplots(2, len(CTS), figsize=(15.6, 7.0),
                         gridspec_kw=dict(hspace=0.72, wspace=0.30))
for r, (matter, stages) in enumerate((('WM', WM), ('GM', GM))):
    ramp = ['#9a9a9a', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
    for c, (key, nm) in enumerate(CTS):
        ax = axes[r, c]
        tops = []
        for j, (lt, lab) in enumerate(stages):
            v = vals(key, matter, lt)
            if not len(v):
                continue
            col = ramp[j] if matter == 'WM' else ['#9a9a9a', '#f7c3ab', '#eb6834'][j]
            rng = np.random.default_rng(r * 40 + c * 8 + j)
            ax.scatter(j + rng.uniform(-0.19, 0.19, len(v)), v, s=30, facecolor='white',
                       edgecolor=col, linewidths=1.4, zorder=5)
            ax.hlines(np.median(v), j - 0.32, j + 0.32, color=col, lw=2.8, zorder=6)
            tops.extend(v)
        if not tops:
            continue
        hi = np.percentile(tops, 97)
        ax.set_ylim(-hi * 0.06, hi * 1.46)
        ax.set_xticks(range(len(stages)))
        ax.set_xticklabels([l for _, l in stages], fontsize=8.4,
                           rotation=34, ha='right')
        ax.set_xlim(-0.7, len(stages) - 0.3)
        ax.set_title(nm, fontsize=11.5, pad=6, color=F.INK)
        ax.tick_params(labelsize=8.5)
        if c == 0:
            ax.set_ylabel(('白质' if matter == 'WM' else '灰质') +
                          '\nIFIH1  每百万转录本计数', fontsize=9.5)
        for j, (lt, lab) in enumerate(stages[1:], start=1):
            row = S[(S.区室 == matter) & (S.细胞类型 == nm) & (S.分期 == lt)]
            if row.empty or not np.isfinite(row.iloc[0].q):
                continue
            sg = F.stars(row.iloc[0].q)
            ax.text(j, hi * 1.10, sg, ha='center', va='bottom', fontsize=10.5,
                    color=F.INK if sg not in ('ns', '') else F.MUT)
        ax.text(0, hi * 1.10, '对照', ha='center', va='bottom', fontsize=8.2, color=F.MUT)
F.star_key(fig, y=-0.004)
fig.savefig('/home/user/aaa/results/mda5/figs/fig_lesion.png', dpi=240,
            bbox_inches='tight', facecolor='white')
print('ok')
pd.set_option('display.width', 200)
print(S[['区室', '细胞类型', '标签', 'n', 'n对照', '中位', '中位对照', 'p', 'q']]
      .to_string(index=False, float_format=lambda v: f'{v:8.4f}'))
