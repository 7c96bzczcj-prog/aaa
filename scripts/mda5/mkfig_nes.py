"""Nestin（NES）在神经系统各类细胞中的表达，版式同 IFIH1 那一页。

左：CELLxGENE Census 正常人脑，按细胞大类的伪整体表达量。
右：正常白质中检出 NES 的核的表达水平，横轴下标注阳性率。
两张图合起来才是完整分布：右图的小提琴是"表达者表达多少"，
阳性率是"多少细胞在表达"，左图是两者的乘积。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import figstyle as F

FIG = '/home/user/aaa/results/mda5/figs/'
RES = '/home/user/aaa/results/mda5/'
GENE = 'NES'

NAME = {'endothelial': '内皮细胞', 'fibroblast': '成纤维细胞',
        'mural/perivascular': '壁细胞/血管周细胞',
        'lymphocyte/leukocyte': '淋巴细胞', 'microglia/CNS macrophage': '小胶质细胞',
        'ependymal/choroid': '室管膜/脉络丛上皮', 'astrocyte': '星形胶质细胞',
        'OPC': '少突胶质前体细胞', 'oligodendrocyte': '少突胶质细胞', 'neuron': '神经元'}
CT180 = {'oligodendrocytes': '少突胶质细胞', 'astrocytes': '星形胶质细胞',
         'immune': '小胶质细胞', 'neurons': '神经元', 'opc': '少突胶质前体细胞',
         'vascular_cells': '内皮细胞', 'lymphocytes': '淋巴细胞'}
# 淋巴细胞必须在列：它是左图的第一名，读者会到右图核对。
ORDER180 = ['lymphocytes', 'vascular_cells', 'opc', 'astrocytes',
            'oligodendrocytes', 'neurons', 'immune']

base = pd.read_csv(RES + 'out_nes_baseline.csv')
base = base[base.grp.isin(NAME)].copy()
base['名'] = base.grp.map(NAME)
base = base.sort_values('CP10k')

d = pd.read_parquet('nes_180759.parquet')
d['CP10k'] = d.NES / np.maximum(d.total, 1) * 1e4
# 用全部核，不限对照：对照白质仅 5291 核，血管 65 个、神经元 34 个，撑不起分布。
# NES 在 MS 与对照之间无差异（六类中五类供体级 p > 0.36），故合并不引入疾病偏倚。
ctrl = d

fig = plt.figure(figsize=(15.8, 5.5))
gs = GridSpec(1, 2, width_ratios=[0.92, 1.30], wspace=0.24, figure=fig)

# ---- 左：Census 伪整体排序 ----
ax = fig.add_subplot(gs[0, 0])
bd = pd.read_csv(RES + 'out_nes_bydataset.csv')
F.census_panel(ax, bd, 'NES', highlight=('endothelial',))
ax.set_title('正常人脑，每点一个数据集，横线为中位', fontsize=12.5,
             pad=10, color=F.INK)

# ---- 右：白质单核层面 ----
ax2 = fig.add_subplot(gs[0, 1])
rng = np.random.default_rng(3)
tops = []; xlab = []
order = [k for k in ORDER180 if (ctrl.ct == k).sum() >= 200]
for i, k in enumerate(order):
    v = ctrl[ctrl.ct == k].CP10k.values
    nz = v[v > 0]
    col = F.BLUE if k in ('vascular_cells', 'opc', 'lymphocytes') else '#9fb3c8'
    if len(nz) >= 3:
        p = ax2.violinplot([np.log10(nz)], positions=[i], widths=0.74,
                           showextrema=False, showmedians=False)
        for b in p['bodies']:
            b.set_facecolor(col); b.set_edgecolor('none'); b.set_alpha(0.30)
        sm = nz if len(nz) <= 500 else rng.choice(nz, 500, replace=False)
        ax2.scatter(i + rng.uniform(-0.15, 0.15, len(sm)), np.log10(sm),
                    s=2.6, color=col, alpha=0.5, linewidths=0, rasterized=True)
        ax2.hlines(np.log10(np.median(nz)), i - 0.3, i + 0.3, color=col, lw=2.8, zorder=6)
        tops.append(np.log10(nz).max())
        if len(nz) < 30:
            ax2.text(i, np.log10(nz).max() + 0.13, f'{len(nz)} 核',
                     ha='center', va='bottom', fontsize=8, color=F.MUT)
    xlab.append(f'{CT180[k]}\n阳性率 {100 * len(nz) / len(v):.1f}%'
                .replace('少突胶质前体细胞', '少突胶质\n前体细胞'))
ax2.set_xticks(range(len(order))); ax2.set_xticklabels(xlab, fontsize=8.4)
ax2.set_xlim(-0.7, len(order) - 0.3)
ax2.set_ylabel(f'{GENE}（每万转录本中的计数）', fontsize=10.5)
ax2.set_yticks([0, 1, 2]); ax2.set_yticklabels(['1', '10', '100'])
ax2.set_ylim(-0.55, max(tops) + 0.35)
ax2.set_title(f'阳性核的 {GENE} 表达水平（白质，66 432 核）', fontsize=12.5, pad=10, color=F.INK)
ax2.tick_params(labelsize=9.5)

# 来源不烧进图片：由 PPT 右下角的可编辑文本框承担
fig.savefig(FIG + 'fig_nes_celltype.png', dpi=240, bbox_inches='tight', facecolor='white')
print('fig_nes_celltype ok')
print(base[['名', 'n', 'CP10k', '阳性率']].sort_values('CP10k', ascending=False)
      .to_string(index=False, float_format=lambda v: f'{v:10.4f}'))
