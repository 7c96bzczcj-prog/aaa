"""统一口径下的第一问与第二问各图。

纵轴一律是 IFIH1 检出率（深度匹配区间内，阳性核占比）。
这个量在任何归一化下都精确，三个队列可比。
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as F
from detect import nuclei, donor_rates, wm_content, BAND

FIG = '/home/user/aaa/results/mda5/figs/'
RES = '/home/user/aaa/results/mda5/'
ORDER = ['oligodendrocyte', 'microglia', 'astrocyte', 'OPC', 'neuron', 'endothelial']
CN = {k: F.CT_CN[k] for k in ORDER}
WM_COH = ['GSE180759', 'GSE279180']

N = nuclei()
g = donor_rates(N)
eff = pd.read_csv(RES + 'out_detect_effects.csv')
per = pd.read_csv(RES + 'out_detect_bycohort.csv')
PCT = lambda v: v * 100          # 画成百分数


def dots(ax, pos, vals, color, width=0.34, seed=0, ms=54):
    v = np.asarray(vals, float)
    if not len(v):
        return
    rng = np.random.default_rng(seed)
    jit = rng.uniform(-width * 0.45, width * 0.45, len(v)) if len(v) > 1 else np.zeros(1)
    ax.scatter(pos + jit, PCT(v), s=ms, facecolor='white', edgecolor=color,
               linewidths=1.7, zorder=5, clip_on=False)
    ax.hlines(PCT(np.median(v)), pos - width, pos + width, color=color, lw=3.0, zorder=6)


def rates(ct, coh, st):
    x = g[(g.细胞类型 == ct) & (g.队列 == coh) & (g.状态 == st)]
    return x.检出率.values


def tidy(ax, labels, positions, title, ylab='IFIH1 阳性核比例 %'):
    ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=9.6)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.set_ylabel(ylab, fontsize=10.5)
    ax.set_title(title, fontsize=12.5, pad=9, color=F.INK)
    ax.tick_params(labelsize=9.5)


# ======================================================== 疾病效应，按区室
def fig_change(comp_label, cohorts, fname, ct_col, ms_col):
    fig, axes = plt.subplots(2, 3, figsize=(14.4 if len(cohorts) > 1 else 12.4, 7.4))
    for ax, ct in zip(axes.ravel(), ORDER):
        pos = 0; ticks = []; labs = []; seps = []; ann = []
        for k, coh in enumerate(cohorts):
            c0 = rates(ct, coh, 'Control'); c1 = rates(ct, coh, 'MS')
            dots(ax, pos, c0, ct_col, seed=10 + k); ticks.append(pos)
            labs.append(f'对照\n{len(c0)} 例')
            dots(ax, pos + 1, c1, ms_col, seed=20 + k); ticks.append(pos + 1)
            labs.append(f'MS\n{len(c1)} 例')
            ann.append((pos, pos + 1, coh))
            pos += 2
            if k < len(cohorts) - 1:
                seps.append(pos - 0.5); pos += 0.6
        tidy(ax, labs, ticks, CN[ct])
        allv = [v for coh in cohorts for st in ('Control', 'MS') for v in rates(ct, coh, st)]
        top = PCT(max(allv)) if allv else 1.0
        ax.set_ylim(-top * 0.06, top * 1.44)
        for s in seps:
            ax.axvline(s, color=F.HAIR, lw=1.0)
        for x0, x1, coh in ann:
            row = per[(per.键 == ct) & (per.队列 == coh)]
            if len(row) and np.isfinite(row.iloc[0].get('log2比值', np.nan)):
                r = row.iloc[0]
                if not (r.率对照 > 0):
                    txt = '对照组未检出'; solid = False
                else:
                    sg = F.stars(r.get('p', float('nan')))
                    txt = f"×{2**r['log2比值']:.1f}  {sg}\n[{2**r.lo:.1f}, {2**r.hi:.1f}]"
                    solid = sg not in ('ns', '')
            else:
                txt = '样本量不足'; solid = False
            F.bracket(ax, x0, x1, top * 1.08, txt, color=F.INK if solid else F.MUT, fs=9.4)
        if len(cohorts) > 1:
            for x0, x1, coh in ann:
                ax.text((x0 + x1) / 2, -0.205, coh, transform=ax.get_xaxis_transform(),
                        ha='center', va='top', fontsize=9.2, color=F.MUT)
    fig.tight_layout(h_pad=3.4 if len(cohorts) > 1 else 2.6, w_pad=2.2)
    F.star_key(fig, y=-0.012)
    fig.savefig(FIG + fname, dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print(fname, 'ok')


# ======================================================== 正常人的区室对比
def fig_region():
    fig, axes = plt.subplots(2, 3, figsize=(13.6, 7.2))
    for ax, ct in zip(axes.ravel(), ORDER):
        w = np.concatenate([rates(ct, c, 'Control') for c in WM_COH]) if True else []
        m = rates(ct, 'Schirmer 2019', 'Control')
        dots(ax, 0, w, F.WM_MS, seed=1); dots(ax, 1, m, F.GM_MS, seed=2)
        tidy(ax, [f'皮层下白质\n{len(w)} 例', f'皮层标本\n{len(m)} 例'], [0, 1], CN[ct])
        top = PCT(max([*w, *m, 1e-4]))
        ax.set_ylim(-top * 0.06, top * 1.34)
        if len(w) >= 3 and len(m) >= 3:
            from scipy.stats import mannwhitneyu
            p = mannwhitneyu(w, m, alternative='two-sided').pvalue
            F.bracket(ax, 0, 1, top * 1.10, f'{F.stars(p)}   p = {p:.3f}',
                      color=F.INK if p < .05 else F.MUT, fs=10.5)
        else:
            F.bracket(ax, 0, 1, top * 1.10, '样本量不足', color=F.MUT, fs=10.5)
    fig.tight_layout(h_pad=2.6, w_pad=2.2)
    F.star_key(fig, y=-0.012)
    fig.savefig(FIG + 'fig_q1_region.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_q1_region ok')


# ======================================================== 交互项
def _iv(ax, lo, hi, est, y, col, xlo, xhi, lw=2.6, ms=62):
    a, b = max(lo, xlo), min(hi, xhi)
    if a < b:
        ax.plot([a, b], [y, y], lw=lw, color=col, solid_capstyle='round', zorder=3)
    for edge, beyond in ((xhi, hi > xhi), (xlo, lo < xlo)):
        if beyond:
            d = (xhi - xlo) * 0.028 * (1 if edge == xhi else -1)
            ax.annotate('', xy=(edge + d, y), xytext=(edge, y), annotation_clip=False,
                        arrowprops=dict(arrowstyle='-|>', lw=lw, color=col, mutation_scale=11))
    if xlo <= est <= xhi:
        ax.scatter([est], [y], s=ms, facecolor=col, edgecolor='white', linewidths=1.3, zorder=4)


def _zero_ctrl(ct, cohorts):
    """该区室里是否有任一队列的对照中位为零；为零则比值无定义。"""
    for coh in cohorts:
        row = per[(per.键 == ct) & (per.队列 == coh)]
        if len(row) and np.isfinite(row.iloc[0].get('log2比值', np.nan)) \
                and not (row.iloc[0].率对照 > 0):
            return True
    return False


def fig_interaction():
    rows = [CN[c] for c in ORDER]
    m = eff.set_index('细胞类型')
    BAD = {CN[c]: (_zero_ctrl(c, WM_COH), _zero_ctrl(c, ['Schirmer 2019'])) for c in ORDER}
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 4.9),
                             gridspec_kw=dict(width_ratios=[1.05, 1.0], wspace=0.26))
    ypos = np.arange(len(rows))[::-1]

    ax = axes[0]; XL, XH = -2.4, 3.2
    ax.axvline(0, color='#9a9a9a', lw=1.0, ls='--')
    for i, nm in enumerate(rows):
        r = m.loc[nm]
        for kk, (off, a, col) in enumerate(((+0.18, 'ΔWM', F.WM_MS), (-0.18, 'ΔGM', F.GM_MS))):
            if not np.isfinite(r.get(a, np.nan)):
                ax.text(XL + (XH - XL) * 0.012, ypos[i] + off, '样本量不足',
                        ha='left', va='center', fontsize=9, color=F.MUT, style='italic')
                continue
            if BAD[nm][kk]:
                ax.text(XL + (XH - XL) * 0.012, ypos[i] + off, '对照组未检出',
                        ha='left', va='center', fontsize=9, color=F.MUT, style='italic')
                continue
            _iv(ax, r[a + '_lo'], r[a + '_hi'], r[a], ypos[i] + off, col, XL, XH)
    ax.set_xlim(XL, XH)
    ax.set_yticks(ypos); ax.set_yticklabels(rows, fontsize=11)
    ax.set_xlabel('IFIH1 阳性核比例的 log2 倍数变化', fontsize=10.5)
    ax.set_title('各区域内 MS 相对对照', fontsize=12.5, pad=10, color=F.INK)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.legend(handles=[Line2D([], [], color=F.WM_MS, lw=2.6, marker='o', ms=7, label='皮层下白质'),
                       Line2D([], [], color=F.GM_MS, lw=2.6, marker='o', ms=7, label='皮层标本')],
              loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=10.5)

    ax = axes[1]; XL2, XH2 = -3.0, 3.6
    ax.axvline(0, color='#9a9a9a', lw=1.0, ls='--')
    for i, nm in enumerate(rows):
        r = m.loc[nm]
        if not np.isfinite(r.get('交互', np.nan)) or any(BAD[nm]):
            msg = '对照组未检出，无法判定' if any(BAD[nm]) else '样本量不足，无法判定'
            ax.text(0.3, ypos[i], msg, ha='center', va='center',
                    fontsize=9.5, color=F.MUT, style='italic')
            continue
        _iv(ax, r['交互_lo'], r['交互_hi'], r['交互'], ypos[i], '#3a3a3a', XL2, XH2, lw=2.8, ms=66)
    ax.set_xlim(XL2, XH2)
    ax.set_yticks(ypos); ax.set_yticklabels(['' for _ in rows])
    ax.set_xlabel('白质倍数变化  减去  皮层倍数变化', fontsize=10.5)
    ax.set_title('两区域变化幅度之差', fontsize=12.5, pad=10, color=F.INK)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.text(0.5, -0.20, '仅两类细胞可判定，置信区间均包含零',
            transform=ax.transAxes, ha='center', va='top', fontsize=10, color=F.MUT)
    fig.savefig(FIG + 'fig_q2_interaction.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_q2_interaction ok')


# ======================================================== 皮层块的白质混入
def fig_confound():
    """皮层组织块的细胞组成：MS 块取得更深，白质成分显著更多。"""
    from scipy.stats import mannwhitneyu
    sch = N[N.队列 == 'Schirmer 2019']
    tot = sch.groupby('供体').size()
    st = sch.groupby('供体').状态.first()
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    CTS = [('neuron', '神经元\n皮层标志'), ('oligodendrocyte', '少突胶质细胞\n白质标志'),
           ('microglia', '小胶质细胞')]
    for i, (key, nm) in enumerate(CTS):
        frac = (sch[sch.细胞类型 == key].groupby('供体').size()
                .reindex(tot.index).fillna(0) / tot)
        for j, cond in enumerate(['Control', 'MS']):
            v = frac[st == cond].values * 100
            col = F.GM_CT if cond == 'Control' else F.GM_MS
            x = i * 2.6 + j
            rng = np.random.default_rng(i * 7 + j)
            ax.scatter(x + rng.uniform(-0.17, 0.17, len(v)), v, s=52, facecolor='white',
                       edgecolor=col, linewidths=1.7, zorder=5)
            ax.hlines(np.median(v), x - 0.33, x + 0.33, color=col, lw=3.0, zorder=6)
        pv = mannwhitneyu(frac[st == 'MS'], frac[st == 'Control']).pvalue
        ax.text(i * 2.6 + 0.5, 112, f'{nm}\np = {pv:.3f}', ha='center', va='top',
                fontsize=10.5, color=F.INK if pv < .05 else F.MUT)
    ax.set_xticks([i * 2.6 + j for i in range(3) for j in (0, 1)])
    ax.set_xticklabels(['对照\n9 例', 'MS\n12 例'] * 3, fontsize=9.6)
    ax.set_xlim(-0.85, 2 * 2.6 + 1.85); ax.set_ylim(-3, 113)
    ax.set_ylabel('占该标本核数的百分比', fontsize=10.5)
    fig.savefig(FIG + 'fig_confound.png', dpi=240, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('fig_confound ok')


if __name__ == '__main__':
    fig_region()
    fig_change('白质', WM_COH, 'fig_q2_wm.png', F.WM_CT, F.WM_MS)
    fig_change('皮层块', ['Schirmer 2019'], 'fig_q2_gm.png', F.GM_CT, F.GM_MS)
    fig_interaction()
    fig_confound()
