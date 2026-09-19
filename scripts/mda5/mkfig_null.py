"""正常供体中，IFIH1 的白质/皮层比值在全基因组分布中的百分位。

白质与皮层标本来自不同研究，系统性批次偏移会同等作用于所有基因，
因此全基因组 log2(白质/皮层) 的中位即偏移量。扣掉后取 IFIH1 的百分位，
该百分位不受偏移影响。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as F

NEW = {'少突胶质': '少突胶质细胞', '小胶质/巨噬': '小胶质细胞', '星形胶质': '星形胶质细胞',
       'OPC': '少突胶质前体细胞', '神经元': '神经元', '内皮/血管': '内皮细胞'}
ORDER = list(NEW.values())
COH = {'GSE180759': '#a9cbf4', 'GSE279180': '#2a78d6'}

d = pd.read_csv('/home/user/aaa/results/mda5/out_gm_wm_genomewide_null.csv')
d = d[d.状态 == 'Control'].copy()
d['名'] = d.细胞类型.map(NEW)

fig, ax = plt.subplots(figsize=(9.6, 4.8))
ypos = {ct: len(ORDER) - 1 - i for i, ct in enumerate(ORDER)}
OFF = {'GSE180759': +0.17, 'GSE279180': -0.17}

ax.axvline(50, color='#9a9a9a', lw=1.0, ls='--')
for ct in ORDER:
    for coh in COH:
        r = d[(d.名 == ct) & (d.白质队列 == coh)]
        y = ypos[ct] + OFF[coh]
        if r.empty:
            ax.text(24, y, '该队列无对照样本', ha='center', va='center',
                    fontsize=8.8, color=F.MUT, style='italic')
            continue
        v = float(r.iloc[0].百分位)
        ax.plot([50, v], [y, y], lw=2.3, color=COH[coh], solid_capstyle='round', zorder=3)
        ax.scatter([v], [y], s=68, color=COH[coh], edgecolor='white', linewidths=1.2, zorder=4)
        ax.text(v + (2.8 if v >= 50 else -2.8), y, f'{v:.0f}',
                ha='left' if v >= 50 else 'right', va='center', fontsize=9.2, color=F.INK)

ax.set_xlim(0, 100)
ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=11)
ax.set_xlabel('IFIH1 的白质/皮层比值在全基因组分布中的百分位', fontsize=10.5)
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
ax.set_ylim(-0.55, len(ORDER) - 0.45)
ax.legend(handles=[Line2D([], [], color=COH[c], lw=2.4, marker='o', ms=7, label=c) for c in COH],
          loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=10.5)
fig.savefig('/home/user/aaa/results/mda5/figs/figB_null.png',
            dpi=240, bbox_inches='tight', facecolor='white')
print('ok')
