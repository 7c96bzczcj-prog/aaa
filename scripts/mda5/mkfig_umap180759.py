"""GSE180759 白质图谱：左侧细胞类型 UMAP，右侧 IFIH1 按病灶分期的特征图。

与参考图同一语法：类群名直接标在图上，特征图零表达为浅灰底，
表达细胞按值叠在上层，五个分期共用一个色标。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import figstyle as F

obs = pd.read_csv('umap_obs_harmony.csv')
x = obs.UMAP1.values; y = obs.UMAP2.values
v = obs.IFIH1_lognorm.values

PAL = {'oligodendrocytes': '#e8a33d', 'astrocytes': '#e0585b', 'immune': '#2fa36b',
       'neurons': '#1f9c9c', 'opc': '#17a2b8', 'vascular_cells': '#4a63c8',
       'lymphocytes': '#c46ec4'}
STAGES = ['control_white_matter', 'MS_periplaque_white_matter',
          'chronic_inactive_MS_lesion_edge', 'chronic_active_MS_lesion_edge',
          'MS_lesion_core']
STAGES = [s for s in STAGES if (obs.pathology.values == s).any()]

fig = plt.figure(figsize=(16.0, 7.4))
gs = GridSpec(2, 1 + len(STAGES), width_ratios=[1.75] + [1] * len(STAGES),
              height_ratios=[1.32, 1.0], wspace=0.06, hspace=0.42, figure=fig)

ax0 = fig.add_subplot(gs[0, 0])
F.umap_annotated(ax0, x, y, obs.cell_type.values, PAL, s=1.5, label_fs=10.5)
ax0.set_title('细胞类型', fontsize=13, color=F.INK, pad=8)

vmax = np.quantile(v[v > 0], 0.99) if (v > 0).any() else 1.0
for j, st in enumerate(STAGES):
    ax = fig.add_subplot(gs[0, j + 1])
    m = obs.pathology.values == st
    # 灰底为该分期之外的全部细胞，保持轮廓可见
    ax.scatter(x[~m], y[~m], s=1.2, c='#f2f2f2', linewidths=0, rasterized=True)
    F.umap_feature(ax, x[m], y[m], v[m], s=1.5, vmax=vmax)
    ax.set_xlim(ax0.get_xlim()); ax.set_ylim(ax0.get_ylim())
    ax.set_title(F.STAGE_CN[st], fontsize=11.5, color=F.INK, pad=8)
    ax.text(0.5, -0.045, f'n = {int(m.sum()):,}', transform=ax.transAxes,
            ha='center', va='top', fontsize=8.5, color=F.MUT)

# ---- 下排：按病灶分期的供体级数值，每点一位供体的一个样本 ----
d = pd.read_parquet('df_180759.parquet').reset_index(drop=True)
d['stage'] = obs.pathology.values
d['ct'] = obs.cell_type.values
CTS = [('oligodendrocytes', '少突胶质细胞'), ('immune', '小胶质细胞'),
       ('astrocytes', '星形胶质细胞'), ('opc', '少突胶质前体细胞'),
       ('vascular_cells', '内皮细胞')]
SHORT = {'control_white_matter': '对照', 'MS_periplaque_white_matter': '斑块周围',
         'chronic_inactive_MS_lesion_edge': '慢性非活动',
         'chronic_active_MS_lesion_edge': '慢性活动', 'MS_lesion_core': '病灶核心'}
from matplotlib.gridspec import GridSpecFromSubplotSpec
gb = GridSpecFromSubplotSpec(1, len(CTS), subplot_spec=gs[1, :], wspace=0.40)
rng = np.random.default_rng(7)
for k, (key, nm) in enumerate(CTS):
    ax = fig.add_subplot(gb[0, k])
    sub = d[d.ct == key]
    tops = []; byst = {}
    for j, st in enumerate(STAGES):
        g = sub[sub.stage == st]
        vals = (g.groupby('donor')
                 .apply(lambda t: t.IFIH1.sum() / max(t.total.sum(), 1) * 1e4,
                        include_groups=False)
                 .values)
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            continue
        tops.extend(vals); byst[st] = vals
        jit = rng.uniform(-0.17, 0.17, len(vals)) if len(vals) > 1 else np.zeros(1)
        ax.scatter(j + jit, vals, s=34, facecolor='white', edgecolor=F.STAGE[st],
                   linewidths=1.6, zorder=5, clip_on=False)
        ax.hlines(np.median(vals), j - 0.3, j + 0.3, color=F.STAGE[st], lw=2.6, zorder=6)
    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels([SHORT[s] for s in STAGES], fontsize=8.6, rotation=32, ha='right')
    ax.set_xlim(-0.65, len(STAGES) - 0.35)
    if tops:
        ax.set_ylim(-max(tops) * 0.07, max(tops) * 1.30)
    ax.set_title(nm, fontsize=11.5, pad=6, color=F.INK)
    ax.tick_params(labelsize=8.5)
    if k == 0:
        ax.set_ylabel('IFIH1 表达量\n每万转录本中的计数', fontsize=9.5)

    # 对照 vs 全部 MS 分期合并。逐个分期比不可能显著：
    # 3 位对照对 4 位供体，Mann-Whitney 双侧最小 p 为 0.057。
    from scipy.stats import mannwhitneyu
    ctrl = np.array(byst.get(STAGES[0], []), float)
    msv = np.concatenate([np.asarray(byst.get(st, []), float) for st in STAGES[1:]]) \
        if len(STAGES) > 1 else np.array([])
    if len(ctrl) >= 2 and len(msv) >= 3 and tops:
        pv = mannwhitneyu(msv, ctrl, alternative='two-sided').pvalue
        sg = F.stars(pv)
        F.bracket(ax, 0, len(STAGES) - 1, max(tops) * 1.14,
                  f'对照 vs MS 各期   {sg}', color=F.INK if sg not in ('ns', '') else F.MUT,
                  fs=8.6)

cax = fig.add_axes([0.925, 0.60, 0.007, 0.25])
cb = fig.colorbar(ScalarMappable(norm=Normalize(0, vmax), cmap=F.FEAT), cax=cax)
cb.set_label('IFIH1 表达水平（对数归一化）', fontsize=9.5, color=F.INK)
cb.ax.tick_params(labelsize=8.5, length=2)
cb.outline.set_visible(False)

F.star_key(fig, y=-0.005)
fig.savefig('/home/user/aaa/results/mda5/figs/fig_umap_wm.png',
            dpi=240, bbox_inches='tight', facecolor='white')
print('ok  vmax', round(float(vmax), 3))
print(obs.groupby('pathology').size().to_string())
