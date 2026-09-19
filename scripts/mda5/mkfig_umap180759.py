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

obs = pd.read_csv('umap_obs.csv')
x = obs.UMAP1.values; y = obs.UMAP2.values
v = obs.IFIH1_lognorm.values

PAL = {'oligodendrocytes': '#e8a33d', 'astrocytes': '#e0585b', 'immune': '#2fa36b',
       'neurons': '#1f9c9c', 'opc': '#17a2b8', 'vascular_cells': '#4a63c8',
       'lymphocytes': '#c46ec4'}
STAGES = ['control_white_matter', 'MS_periplaque_white_matter',
          'chronic_inactive_MS_lesion_edge', 'chronic_active_MS_lesion_edge',
          'MS_lesion_core']
STAGES = [s for s in STAGES if (obs.pathology.values == s).any()]

fig = plt.figure(figsize=(16.4, 4.5))
gs = GridSpec(1, 1 + len(STAGES), width_ratios=[1.75] + [1] * len(STAGES),
              wspace=0.06, figure=fig)

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

cax = fig.add_axes([0.915, 0.30, 0.008, 0.36])
cb = fig.colorbar(ScalarMappable(norm=Normalize(0, vmax), cmap=F.FEAT), cax=cax)
cb.set_label('IFIH1  log 归一化表达', fontsize=9.5, color=F.INK)
cb.ax.tick_params(labelsize=8.5, length=2)
cb.outline.set_visible(False)

fig.savefig('/home/user/aaa/results/mda5/figs/fig_umap_wm.png',
            dpi=240, bbox_inches='tight', facecolor='white')
print('ok  vmax', round(float(vmax), 3))
print(obs.groupby('pathology').size().to_string())
