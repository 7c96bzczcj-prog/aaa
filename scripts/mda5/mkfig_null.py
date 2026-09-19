"""第 4 页的批次对照（只针对正常供体）。

第 4 页把白质（Absinta、Lerma）与灰质（Schirmer）相比，跨了研究。
跨研究的系统偏移会同等作用于所有基因，因此该偏移量可由全基因组
log2(白质/灰质) 的中位估计出来。把它扣掉后，看 IFIH1 落在同一分布的
第几百分位；这个百分位对偏移免疫。

左：偏移有多大（灰点为全基因组中位，色点为 IFIH1）。
右：扣掉偏移后 IFIH1 的百分位。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as F

d = pd.read_csv('/home/user/aaa/results/mda5/out_gm_wm_genomewide_null.csv')
d = d[d.状态 == 'Control']
ORDER = ['少突胶质', '小胶质/巨噬', '星形胶质', 'OPC', '神经元', '内皮/血管']
COH = {'GSE180759': '#a9cbf4', 'GSE279180': '#2a78d6'}

rows = []
for ct in ORDER:
    for coh in COH:
        r = d[(d.细胞类型 == ct) & (d.白质队列 == coh)]
        rows.append((ct, coh, None if r.empty else r.iloc[0]))

fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0),
                         gridspec_kw=dict(width_ratios=[1.15, 1.0], wspace=0.24))
ypos = {ct: len(ORDER) - 1 - i for i, ct in enumerate(ORDER)}
OFF = {'GSE180759': +0.17, 'GSE279180': -0.17}

ax = axes[0]
ax.axvline(0, color='#c8c8c8', lw=1.0)
for ct, coh, r in rows:
    if r is None:
        continue
    y = ypos[ct] + OFF[coh]
    ax.plot([r.全基因组中位, r.IFIH1_log2_WM_GM], [y, y], lw=1.6, color='#c8c8c8', zorder=2)
    ax.scatter([r.全基因组中位], [y], s=46, marker='s', color='#9a9a9a',
               zorder=3, linewidths=0)
    ax.scatter([r.IFIH1_log2_WM_GM], [y], s=64, color=COH[coh],
               edgecolor='white', linewidths=1.2, zorder=4)
ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=11)
ax.set_xlabel('log2（白质 / 灰质）', fontsize=10.5)
ax.set_title('跨研究偏移有多大', fontsize=12.5, pad=10, color=F.INK)
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
ax.set_ylim(-0.55, len(ORDER) - 0.45)
ax.legend(handles=[Line2D([], [], marker='s', ls='', color='#9a9a9a', ms=7,
                          label='全基因组中位（即偏移量）'),
                   Line2D([], [], marker='o', ls='', color='#2a78d6', ms=8, label='IFIH1')],
          loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=10)

ax = axes[1]
ax.axvline(50, color='#9a9a9a', lw=1.0, ls='--')
for ct, coh, r in rows:
    if r is None:
        y = ypos[ct] + OFF[coh]
        ax.text(22, y, '该队列缺对照供体', ha='center', va='center',
                fontsize=8.5, color=F.MUT, style='italic')
        continue
    y = ypos[ct] + OFF[coh]
    ax.plot([50, r.百分位], [y, y], lw=2.2, color=COH[coh], solid_capstyle='round', zorder=3)
    ax.scatter([r.百分位], [y], s=66, color=COH[coh], edgecolor='white',
               linewidths=1.2, zorder=4)
    ax.text(r.百分位 + (2.6 if r.百分位 >= 50 else -2.6), y, f'{r.百分位:.0f}',
            ha='left' if r.百分位 >= 50 else 'right', va='center',
            fontsize=9, color=F.INK)
ax.set_xlim(0, 100)
ax.set_yticks(list(ypos.values())); ax.set_yticklabels(['' for _ in ypos])
ax.set_xlabel('扣掉偏移后，IFIH1 在全基因组分布中的百分位', fontsize=10.5)
ax.set_title('扣掉偏移之后', fontsize=12.5, pad=10, color=F.INK)
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
ax.set_ylim(-0.55, len(ORDER) - 0.45)
ax.text(0.5, -0.17, '只有星形胶质在两个白质队列中都高于全基因组中位',
        transform=ax.transAxes, ha='center', va='top', fontsize=10, color=F.MUT)
ax.legend(handles=[Line2D([], [], color=COH[c], lw=2.4, marker='o', ms=7, label=c)
                   for c in COH],
          loc='upper center', bbox_to_anchor=(0.5, -0.26), ncol=2, fontsize=10)

fig.savefig('/home/user/aaa/results/mda5/figs/figB_null.png',
            dpi=240, bbox_inches='tight', facecolor='white')
print('ok')
print(d[['细胞类型', '白质队列', 'IFIH1_log2_WM_GM', '全基因组中位', '百分位']].to_string(index=False))
