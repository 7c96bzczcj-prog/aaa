"""第一问与第二问的全部图。

图型的选择按信息量决定，不按好看决定：
  细胞级小提琴只用在"哪类细胞表达"那页，那里有十几万个核，分布形状有意义；
  区室与疾病的比较一律画供体点，因为统计量就是供体级的，
  且 IFIH1 只在百分之一到五的核里检出，细胞级小提琴在那些页只会是零点上的一根钉。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import figstyle as F
from celllevel import load_cells

FIG = '/home/user/aaa/results/mda5/figs/'
RES = '/home/user/aaa/results/mda5/'
ORDER_EN = ['oligodendrocyte', 'microglia', 'astrocyte', 'OPC', 'neuron', 'endothelial']
CN = {e: F.CT_CN[e] for e in ORDER_EN}

tw = pd.read_parquet('twoway.parquet')
cells = load_cells()
ci = pd.read_csv(RES + 'out_q2_ci.csv')
ci_main = ci[ci.层.isna()] if '层' in ci.columns else ci
ci_str = ci[ci.层.notna()] if '层' in ci.columns else pd.DataFrame()


def donors(ct, comp, coh=None, status=None):
    g = tw[(tw.细胞类型 == ct) & (tw.区室 == comp)]
    if coh: g = g[g.队列 == coh]
    if status: g = g[g.状态 == status]
    return g.值.values / 100.0          # CPM -> CP10k


def dotgroup(ax, pos, vals, color, width=0.34, seed=0, ms=52):
    """一组供体：抖动点 + 中位横线。点是数据，横线是统计量。"""
    v = np.asarray(vals, float)
    if not len(v):
        return
    rng = np.random.default_rng(seed)
    jit = rng.uniform(-width * 0.45, width * 0.45, len(v)) if len(v) > 1 else np.zeros(1)
    ax.scatter(pos + jit, v, s=ms, facecolor='white', edgecolor=color,
               linewidths=1.7, zorder=5, clip_on=False)
    ax.hlines(np.median(v), pos - width, pos + width, color=color, lw=3.0, zorder=6)


def tidy(ax, labels, positions, ylab=None, title=None):
    ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    if ylab: ax.set_ylabel(ylab, fontsize=10.5)
    if title: ax.set_title(title, fontsize=12.5, pad=9, color=F.INK)
    ax.tick_params(labelsize=9.5)


# ==========================================================================
# 图一（第 3 页）  哪些细胞表达 IFIH1
# ==========================================================================
def fig_celltype():
    base = pd.read_csv(RES + 'out_baseline_consolidated.csv')
    NAME = {'endothelial': '内皮', 'fibroblast': '成纤维', 'mural/perivascular': '周细胞/血管周',
            'lymphocyte/leukocyte': '淋巴/白细胞', 'microglia/CNS macrophage': '小胶质/巨噬',
            'ependymal/choroid': '室管膜/脉络丛', 'astrocyte': '星形胶质', 'OPC': 'OPC',
            'oligodendrocyte': '少突胶质', 'neuron': '神经元'}
    base = base[base.grp.isin(NAME)].copy()
    base['名'] = base.grp.map(NAME)
    base = base.sort_values('CP10k')

    fig = plt.figure(figsize=(14.6, 5.5))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.15], wspace=0.28, figure=fig)

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(base))
    for yi, (_, r) in zip(y, base.iterrows()):
        c = F.VIOLET if r.grp == 'endothelial' else (
            F.BLUE if r.grp == 'microglia/CNS macrophage' else '#9fb3c8')
        ax.hlines(yi, 0, r.CP10k, color=c, lw=2.2, alpha=0.55)
        ax.scatter([r.CP10k], [yi], s=74, color=c, zorder=4, linewidths=0)
        ax.text(r.CP10k * 1.09, yi, f'{r.CP10k:.2f}', va='center', ha='left',
                fontsize=9.5, color=F.INK)
    ax.set_yticks(y); ax.set_yticklabels(base.名.values, fontsize=11)
    ax.set_xlim(0, base.CP10k.max() * 1.28)
    ax.set_xlabel('IFIH1  CP10k（伪整体）', fontsize=10.5)
    ax.set_title('正常人脑，783 位供体，1850 万个核', fontsize=12.5, pad=10, color=F.INK)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color=F.HAIR, lw=0.7, alpha=0.7); ax.set_axisbelow(True)

    # 右：正常白质里，检出 IFIH1 的核的表达水平；检出率单列标注。
    # 只画检出核，因为九成以上的核为零，含零的小提琴只会是零点上的一根钉；
    # 零的信息由上方的检出率给出，两者合起来才是完整的分布。
    ax2 = fig.add_subplot(gs[0, 1])
    wm = cells[(cells.区室 == '白质') & (cells.状态 == 'Control')]
    order = [e for e in ORDER_EN if (wm.细胞类型 == e).sum() > 200]
    rng = np.random.default_rng(3)
    tops = []; xlab = []
    for i, e in enumerate(order):
        v = wm[wm.细胞类型 == e].CP10k.values
        nz = v[v > 0]
        col = F.BLUE if e in ('endothelial', 'microglia') else '#9fb3c8'
        if len(nz) > 5:
            p = ax2.violinplot([np.log10(nz)], positions=[i], widths=0.74,
                               showextrema=False, showmedians=False)
            for b in p['bodies']:
                b.set_facecolor(col); b.set_edgecolor('none'); b.set_alpha(0.30)
            s = nz if len(nz) <= 500 else rng.choice(nz, 500, replace=False)
            ax2.scatter(i + rng.uniform(-0.15, 0.15, len(s)), np.log10(s),
                        s=2.6, color=col, alpha=0.5, linewidths=0, rasterized=True)
            ax2.hlines(np.log10(np.median(nz)), i - 0.3, i + 0.3, color=col, lw=2.8, zorder=6)
            tops.append(np.log10(nz).max())
        xlab.append(f'{CN[e]}\n检出 {100 * len(nz) / len(v):.1f}%')
    tidy(ax2, xlab, list(range(len(order))), None,
         '正常白质，检出 IFIH1 的核的表达水平')
    ax2.set_ylabel('IFIH1  CP10k', fontsize=10.5)
    ticks = np.array([0, 1, 2])
    ax2.set_yticks(ticks); ax2.set_yticklabels(['1', '10', '100'])
    ax2.set_ylim(-0.55, max(tops) + 0.35)
    ax2.tick_params(axis='x', rotation=0, labelsize=9.5)
    fig.savefig(FIG + 'fig_q1_celltype.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_q1_celltype ok')


# ==========================================================================
# 图二（第 4 页）  正常人：白质对灰质
# ==========================================================================
def fig_region_normal():
    fig, axes = plt.subplots(2, 3, figsize=(13.6, 7.2))
    g2x2 = pd.read_csv(RES + 'out_gm_wm_2x2.csv')
    for ax, e in zip(axes.ravel(), ORDER_EN):
        w = donors(e, '白质', status='Control')
        g = donors(e, '灰质', status='Control')
        dotgroup(ax, 0, w, F.WM_MS, seed=1); dotgroup(ax, 1, g, F.GM_MS, seed=2)
        lab = [f'白质\n{len(w)} 人', f'灰质\n{len(g)} 人']
        tidy(ax, lab, [0, 1], 'IFIH1  CP10k', CN[e])
        top = max([*w, *g, 0.01]) if len(w) + len(g) else 1
        ax.set_ylim(-top * 0.06, top * 1.30)
        row = g2x2[(g2x2.细胞类型 == CN[e]) & (g2x2.状态 == '正常')]
        if len(w) < 3 or len(g) < 3:
            txt = '供体不足'
        elif len(row):
            txt = f'p = {float(row.iloc[0].p):.3f}'
        else:
            txt = ''
        if txt:
            F.bracket(ax, 0, 1, top * 1.09, txt,
                      color=F.MUT if ('不足' in txt or float(row.iloc[0].p) >= .05) else F.INK,
                      fs=10.5)
    fig.suptitle('', y=1)
    fig.tight_layout(h_pad=2.4, w_pad=2.0)
    fig.savefig(FIG + 'fig_q1_region.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_q1_region ok')


# ==========================================================================
# 图三、四（第 6、7 页）  各区室内部，MS 相对对照
# ==========================================================================
def _change_panel(ax, e, comp, cohorts, cmap_ct, cmap_ms):
    pos = 0; ticks = []; labs = []; seps = []; ann = []
    for k, coh in enumerate(cohorts):
        c0 = donors(e, comp, coh, 'Control'); c1 = donors(e, comp, coh, 'MS')
        dotgroup(ax, pos, c0, cmap_ct, seed=10 + k); ticks.append(pos)
        labs.append(f'对照\n{len(c0)} 人')
        dotgroup(ax, pos + 1, c1, cmap_ms, seed=20 + k); ticks.append(pos + 1)
        labs.append(f'MS\n{len(c1)} 人')
        ann.append((pos, pos + 1, coh, len(c0), len(c1)))
        pos += 2
        if k < len(cohorts) - 1:
            seps.append(pos - 0.5); pos += 0.6
    return ticks, labs, seps, ann


def fig_change(comp, cohorts, fname, title_note):
    fig, axes = plt.subplots(2, 3, figsize=(14.4 if len(cohorts) > 1 else 12.6, 7.4))
    ct_c = F.WM_CT if comp == '白质' else F.GM_CT
    ms_c = F.WM_MS if comp == '白质' else F.GM_MS
    for ax, e in zip(axes.ravel(), ORDER_EN):
        ticks, labs, seps, ann = _change_panel(ax, e, comp, cohorts, ct_c, ms_c)
        tidy(ax, labs, ticks, 'IFIH1  CP10k', CN[e])
        allv = [v for coh in cohorts for st in ('Control', 'MS')
                for v in donors(e, comp, coh, st)]
        top = max(allv) if allv else 1.0
        ax.set_ylim(-top * 0.06, top * 1.42)
        for s in seps:
            ax.axvline(s, color=F.HAIR, lw=1.0)
        for x0, x1, coh, n0, n1 in ann:
            row = ci_str[(ci_str.细胞类型 == CN[e]) & (ci_str.层 == f'{comp}|{coh}')]
            if len(row):
                r = row.iloc[0]
                txt = f'log2FC {r.单层log2FC:+.2f}\n[{r.单层_lo:+.2f}, {r.单层_hi:+.2f}]'
                solid = (r.单层_lo > 0) or (r.单层_hi < 0)
            else:
                txt = '供体不足'; solid = False
            F.bracket(ax, x0, x1, top * 1.06, txt, color=F.INK if solid else F.MUT, fs=9.2)
        if len(cohorts) > 1:
            for x0, x1, coh, _, _ in ann:
                ax.text((x0 + x1) / 2, -0.205, coh, transform=ax.get_xaxis_transform(),
                        ha='center', va='top', fontsize=9.2, color=F.MUT)
    fig.tight_layout(h_pad=3.4 if len(cohorts) > 1 else 2.6, w_pad=2.0)
    fig.savefig(FIG + fname, dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print(fname, 'ok', title_note)


# ==========================================================================
# 图五（第 8 页）  两个区室的改变是否相同
# ==========================================================================
def _interval(ax, x0, x1, est, y, col, xlo, xhi, lw=2.6, ms=62):
    """画区间；超出坐标范围的一端改成箭头，避免个别宽区间把横轴拉散。"""
    a, b = max(x0, xlo), min(x1, xhi)
    if a < b:
        ax.plot([a, b], [y, y], lw=lw, color=col, solid_capstyle='round', zorder=3)
    if x1 > xhi:
        ax.annotate('', xy=(xhi + (xhi - xlo) * 0.028, y), xytext=(xhi, y),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', lw=lw, color=col, mutation_scale=11))
    if x0 < xlo:
        ax.annotate('', xy=(xlo - (xhi - xlo) * 0.028, y), xytext=(xlo, y),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', lw=lw, color=col, mutation_scale=11))
    if xlo <= est <= xhi:
        ax.scatter([est], [y], s=ms, facecolor=col, edgecolor='white',
                   linewidths=1.3, zorder=4)


def fig_interaction():
    rows = [CN[e] for e in ORDER_EN]
    m = ci_main.set_index('细胞类型')
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 4.9),
                             gridspec_kw=dict(width_ratios=[1.05, 1.0], wspace=0.26))
    ypos = np.arange(len(rows))[::-1]

    ax = axes[0]
    XL, XH = -1.8, 4.2
    ax.axvline(0, color='#9a9a9a', lw=1.0, ls='--')
    for i, nm in enumerate(rows):
        r = m.loc[nm]
        for off, a, lo, hi, col in ((+0.18, 'ΔWM', 'ΔWM_lo', 'ΔWM_hi', F.WM_MS),
                                    (-0.18, 'ΔGM', 'ΔGM_lo', 'ΔGM_hi', F.GM_MS)):
            if not np.isfinite(r.get(a, np.nan)):
                ax.text(XL + (XH - XL) * 0.012, ypos[i] + off, '该区室供体不足',
                        ha='left', va='center', fontsize=9, color=F.MUT, style='italic')
                continue
            _interval(ax, r[lo], r[hi], r[a], ypos[i] + off, col, XL, XH)
    ax.set_xlim(XL, XH)
    ax.set_yticks(ypos); ax.set_yticklabels(rows, fontsize=11)
    ax.set_xlabel('MS 相对对照的改变   log2 倍数变化', fontsize=10.5)
    ax.set_title('各区室内部的改变', fontsize=12.5, pad=10, color=F.INK)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.legend(handles=[Line2D([], [], color=F.WM_MS, lw=2.6, marker='o', ms=7, label='白质'),
                       Line2D([], [], color=F.GM_MS, lw=2.6, marker='o', ms=7, label='灰质')],
              loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=10.5)

    ax = axes[1]
    XL2, XH2 = -5.2, 4.2
    ax.axvline(0, color='#9a9a9a', lw=1.0, ls='--')
    for i, nm in enumerate(rows):
        r = m.loc[nm]
        if not np.isfinite(r.get('交互', np.nan)):
            ax.text(0, ypos[i], '灰质侧供体不足，不可判定', ha='center', va='center',
                    fontsize=9.5, color=F.MUT, style='italic')
            continue
        _interval(ax, r['交互_lo'], r['交互_hi'], r['交互'], ypos[i], '#3a3a3a',
                  XL2, XH2, lw=2.8, ms=66)
    ax.set_xlim(XL2, XH2)
    ax.set_yticks(ypos); ax.set_yticklabels(['' for _ in rows])
    ax.set_xlabel('白质的改变  减去  灰质的改变', fontsize=10.5)
    ax.set_title('两个区室的改变之差', fontsize=12.5, pad=10, color=F.INK)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.text(0.5, -0.20, '四类可判定的细胞，区间全部跨过零：现有供体数判不出两侧改变是否不同',
            transform=ax.transAxes, ha='center', va='top', fontsize=10, color=F.MUT)
    fig.savefig(FIG + 'fig_q2_interaction.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_q2_interaction ok')


if __name__ == '__main__':
    fig_celltype()
    fig_region_normal()
    fig_change('白质', ['GSE180759', 'GSE279180'], 'fig_q2_wm.png', '白质')
    fig_change('灰质', ['Schirmer 2019'], 'fig_q2_gm.png', '灰质')
    fig_interaction()
