"""统一作图样式。

配色约定，全篇不变：
  区室用色相      白质 = 蓝，灰质 = 橙
  疾病状态用深浅  对照 = 浅，MS = 深
所以一张图里蓝橙之分永远是灰白质，深浅之分永远是对照与 MS。
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patheffects import withStroke
import numpy as np

font_manager.fontManager.addfont('/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc')
plt.rcParams.update({
    'font.family': 'WenQuanYi Zen Hei', 'font.size': 11, 'axes.unicode_minus': False,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#4a4a4a', 'axes.linewidth': 0.9,
    'xtick.color': '#4a4a4a', 'ytick.color': '#4a4a4a', 'text.color': '#1a1a1a',
    'axes.labelcolor': '#1a1a1a', 'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.dpi': 240, 'savefig.bbox': 'tight', 'legend.frameon': False,
})

INK = '#1a1a1a'; MUT = '#6b6b6b'; HAIR = '#d8d8d8'
BLUE = '#2a78d6'; ORANGE = '#eb6834'; GREY = '#9a9a9a'
VIOLET = '#4a3aa7'; RED = '#d03b3b'

# 区室 x 状态
WM_MS = '#2a78d6'; WM_CT = '#a9cbf4'
GM_MS = '#eb6834'; GM_CT = '#f7c3ab'
def cc(comp, ms):
    return (WM_MS if ms else WM_CT) if comp == '白质' else (GM_MS if ms else GM_CT)

# 特征图色标：近零为极浅灰，高表达为紫
FEAT = LinearSegmentedColormap.from_list('feat', ['#ececec', '#cfc9e6', '#8b7ec8', '#4a3aa7', '#2b1f6b'])

STAGE = {'control_white_matter': '#9a9a9a', 'MS_periplaque_white_matter': '#a9cbf4',
         'chronic_inactive_MS_lesion_edge': '#5598e7', 'chronic_active_MS_lesion_edge': '#2a78d6',
         'MS_lesion_core': '#104281'}
STAGE_CN = {'control_white_matter': '对照白质',
            'MS_periplaque_white_matter': '病灶周围白质',
            'chronic_inactive_MS_lesion_edge': '慢性非活动性病灶边缘',
            'chronic_active_MS_lesion_edge': '慢性活动性病灶边缘',
            'MS_lesion_core': '病灶核心'}
CT_CN = {'oligodendrocytes': '少突胶质细胞', 'astrocytes': '星形胶质细胞',
         'immune': '小胶质细胞', 'neurons': '神经元', 'opc': '少突胶质前体细胞',
         'vascular_cells': '内皮细胞', 'lymphocytes': '淋巴细胞',
         'oligodendrocyte': '少突胶质细胞', 'astrocyte': '星形胶质细胞',
         'microglia': '小胶质细胞', 'neuron': '神经元', 'OPC': '少突胶质前体细胞',
         'endothelial': '内皮细胞', 'lymphocyte': '淋巴细胞',
         'stromal': '基质细胞', 'mural': '壁细胞'}
CT_ORDER = ['少突胶质细胞', '小胶质细胞', '星形胶质细胞', '少突胶质前体细胞', '神经元', '内皮细胞']


def bare(ax):
    """去掉坐标轴，供 UMAP 使用。"""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def umap_axes_cue(ax, x, y, lab1='UMAP1', lab2='UMAP2', frac=0.17):
    """左下角的短箭头坐标提示，替代完整坐标轴。"""
    x0, x1 = np.min(x), np.max(x); y0, y1 = np.min(y), np.max(y)
    px = x0 - (x1 - x0) * 0.04; py = y0 - (y1 - y0) * 0.04
    lx = (x1 - x0) * frac; ly = (y1 - y0) * frac
    ka = dict(arrowprops=dict(arrowstyle='-|>', lw=1.1, color='#4a4a4a', mutation_scale=9))
    ax.annotate('', xy=(px + lx, py), xytext=(px, py), **ka)
    ax.annotate('', xy=(px, py + ly), xytext=(px, py), **ka)
    ax.text(px + lx / 2, py - (y1 - y0) * 0.035, lab1, ha='center', va='top', fontsize=7.5, color=MUT)
    ax.text(px - (x1 - x0) * 0.028, py + ly / 2, lab2, ha='right', va='center',
            fontsize=7.5, color=MUT, rotation=90)


def umap_annotated(ax, x, y, labels, palette, s=1.6, label_fs=10, cue=True):
    """按细胞类型着色的 UMAP，类群名直接标在质心上。"""
    labels = np.asarray(labels)
    for k in palette:
        m = labels == k
        if not m.any():
            continue
        ax.scatter(x[m], y[m], s=s, c=palette[k], linewidths=0, alpha=0.85, rasterized=True)
    for k in palette:
        m = labels == k
        if m.sum() < 20:
            continue
        cx, cy = np.median(x[m]), np.median(y[m])
        ax.text(cx, cy, CT_CN.get(k, k), ha='center', va='center', fontsize=label_fs,
                color='#141414', fontweight='bold',
                path_effects=[withStroke(linewidth=3.2, foreground='white')])
    bare(ax)
    if cue:
        umap_axes_cue(ax, x, y)


def umap_feature(ax, x, y, v, s=1.6, vmax=None, cmap=FEAT, bg='#e9e9e9', cue=False):
    """IFIH1 特征图：零表达为浅灰底，表达细胞按值叠在上层，高值在最上。"""
    v = np.asarray(v, float)
    nz = v > 0
    ax.scatter(x[~nz], y[~nz], s=s, c=bg, linewidths=0, rasterized=True)
    if nz.any():
        o = np.argsort(v[nz])
        if vmax is None:
            vmax = np.quantile(v[nz], 0.99) or v[nz].max()
        ax.scatter(x[nz][o], y[nz][o], s=s * 1.5, c=v[nz][o], cmap=cmap,
                   vmin=0, vmax=vmax, linewidths=0, rasterized=True)
    bare(ax)
    if cue:
        umap_axes_cue(ax, x, y)
    return vmax


def violin_cells_donors(ax, pos, cell_vals, donor_vals, color, width=0.62,
                        seed=0, maxpts=900, donor_ms=34):
    """一个格子：底层是细胞级分布的小提琴，上层是供体均值的大点。

    统计只用供体点，细胞级仅用来显示分布形状。
    """
    rng = np.random.default_rng(seed)
    v = np.asarray(cell_vals, float)
    if len(v) > 5 and np.ptp(v) > 0:
        p = ax.violinplot([v], positions=[pos], widths=width, showextrema=False, showmedians=False)
        for b in p['bodies']:
            b.set_facecolor(color); b.set_edgecolor('none'); b.set_alpha(0.28)
    nz = v[v > 0]
    if len(nz):
        sm = nz if len(nz) <= maxpts else rng.choice(nz, maxpts, replace=False)
        ax.scatter(pos + rng.uniform(-width * 0.22, width * 0.22, len(sm)), sm,
                   s=2.0, color=color, alpha=0.45, linewidths=0, rasterized=True)
    d = np.asarray(donor_vals, float)
    if len(d):
        ax.scatter(pos + rng.uniform(-width * 0.16, width * 0.16, len(d)), d,
                   s=donor_ms, facecolor='white', edgecolor=color, linewidths=1.5, zorder=6)
        ax.hlines(np.median(d), pos - width * 0.36, pos + width * 0.36,
                  color=color, lw=2.6, zorder=7)


def bracket(ax, x0, x1, y, text, color=INK, fs=10, tick=None):
    """两组之间的横线与标注。"""
    if tick is None:
        tick = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.018
    ax.plot([x0, x0, x1, x1], [y - tick, y, y, y - tick], lw=1.0, color=color, clip_on=False)
    ax.text((x0 + x1) / 2, y + tick * 0.5, text, ha='center', va='bottom',
            fontsize=fs, color=color, clip_on=False)


def forest(ax, labels, est, lo, hi, colors, xlab, ref=0.0, untestable=None, fs=11):
    """区间图。untestable 里的行只画一条说明，不画点。"""
    n = len(labels)
    ypos = np.arange(n)[::-1]
    ax.axvline(ref, color='#9a9a9a', lw=1.0, ls='--', zorder=1)
    for i, y in enumerate(ypos):
        if untestable and labels[i] in untestable:
            ax.text(ref, y, untestable[labels[i]], ha='center', va='center',
                    fontsize=fs - 1.5, color=MUT, style='italic')
            continue
        ax.plot([lo[i], hi[i]], [y, y], lw=2.4, color=colors[i], solid_capstyle='round', zorder=3)
        ax.scatter([est[i]], [y], s=64, facecolor=colors[i], edgecolor='white',
                   linewidths=1.4, zorder=4)
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=fs)
    ax.set_xlabel(xlab)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.7, n - 0.3)


def stars(p):
    """显著性标记。ns 表示不显著，不是"无差异"。"""
    if p is None or p != p:
        return ''
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def star_key(fig, y=0.005, fs=9):
    fig.text(0.5, y, '* p<0.05    ** p<0.01    *** p<0.001    ns 不显著',
             ha='center', va='bottom', fontsize=fs, color=MUT)


# 数据来源，NLM 格式；图左下角标注用
CITE = {
    'census':  'CELLxGENE Census · CZI Cell Science Program, et al. Nucleic Acids Res. 2025;53(D1):D886-D900',
    'absinta': 'GSE180759 · Absinta M, et al. Nature. 2021;597(7878):709-714',
    'lerma':   'GSE279180 · Lerma-Martin C, et al. Nat Neurosci. 2024;27(12):2354-2365',
    'schirmer':'Schirmer L, et al. Nature. 2019;573(7772):75-82',
    'macnair': 'Macnair W, et al. Neuron. 2025;113(3):396-410.e9',
}


def source(fig, keys, y=-0.002, fs=7.8, x=0.0):
    """图左下角的数据来源。keys 为 CITE 的键，按给定次序逐行排列。"""
    if isinstance(keys, str):
        keys = [keys]
    txt = '\n'.join(CITE[k] for k in keys)
    fig.text(x, y, txt, ha='left', va='top', fontsize=fs, color=MUT, linespacing=1.5)
