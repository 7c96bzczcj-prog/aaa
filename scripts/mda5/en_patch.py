"""英文版：运行时拦截 matplotlib 的文字接口做翻译，输出到 figs_en/。

不改动原脚本，因此 DataFrame 的中文列名（区室、队列、细胞类型、状态…）
完全不受影响。只有真正被画到图上的字符串会经过 set_title / set_xlabel /
set_ylabel / set_xticklabels / text / Line2D(label=) 这几个入口，逐一翻译。

用法：  python -c "import en_patch; import mkfig_nes"
"""
import re, matplotlib
matplotlib.use('Agg')
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.colorbar import Colorbar

TR = [
    ('少突胶质\n前体细胞', 'OPCs'), ('少突胶质前体细胞', 'OPCs'), ('少突胶质细胞', 'Oligodendrocytes'),
    ('小胶质/巨噬', 'Microglia'), ('小胶质细胞', 'Microglia'),
    ('星形胶质细胞', 'Astrocytes'), ('壁细胞/血管周细胞', 'Mural / perivascular'),
    ('室管膜/脉络丛上皮', 'Ependymal / choroid'),
    ('内皮/周细胞', 'Endothelial + pericytes'), ('内皮/血管', 'Endothelial'),
    ('内皮细胞', 'Endothelial cells'), ('兴奋性神经元', 'Excitatory neurons'),
    ('抑制性神经元', 'Inhibitory neurons'), ('淋巴细胞', 'Lymphocytes'),
    ('成纤维细胞', 'Fibroblasts'), ('基质细胞', 'Stromal'), ('神经元', 'Neurons'),

    ('慢性非活动性病灶边缘', 'Chronic inactive edge'),
    ('慢性活动性病灶边缘', 'Chronic active edge'),
    ('慢性非活动性病灶', 'Chronic inactive lesion'),
    ('慢性活动性病灶', 'Chronic active lesion'),
    ('再髓鞘化病灶', 'Remyelinating lesion'), ('活动性病灶', 'Active lesion'),
    ('病灶周围白质', 'Periplaque WM'), ('正常表现白质', 'NAWM'),
    ('正常表现灰质', 'NAGM'), ('皮层病灶', 'Cortical lesion'),
    ('病灶核心', 'Lesion core'), ('对照白质', 'Control WM'), ('对照灰质', 'Control GM'),
    ('慢性非活动', 'Chronic inactive'), ('慢性活动', 'Chronic active'),
    ('斑块周围', 'Periplaque'),

    ('正常人脑，每点一个数据集，横线为中位',
     'Normal human brain · one point per dataset · bar = median'),
    ('正常人脑，783 例供体，1850 万个核',
     'Normal human brain · 783 donors · 18.55 M nuclei'),
    ('阳性核的 IFIH1 表达水平（正常白质）',
     'IFIH1 in positive nuclei (normal WM)'),
    ('阳性核的 NES 表达水平（白质，66 432 核）',
     'NES in positive nuclei (WM, 66,432 nuclei)'),
    ('各区域内 MS 相对对照', 'MS vs control within each region'),
    ('两区域变化幅度之差', 'Difference between regions'),
    ('IFIH1 阳性核比例的 log2 倍数变化', 'log2 fold change, % IFIH1-positive nuclei'),
    ('白质倍数变化  减去  皮层倍数变化', 'WM fold change  minus  cortical fold change'),
    ('IFIH1 阳性核比例 %', '% nuclei IFIH1-positive'),
    ('IFIH1 表达量（每万转录本中的计数）', 'IFIH1 (counts per 10,000 transcripts)'),
    ('NES 表达量（每万转录本中的计数）', 'NES (counts per 10,000 transcripts)'),
    ('IFIH1（每万转录本中的计数）', 'IFIH1 (counts per 10,000)'),
    ('NES（每万转录本中的计数）', 'NES (counts per 10,000)'),
    ('IFIH1  每百万转录本计数', 'IFIH1 (counts per million)'),
    ('IFIH1 表达水平（对数归一化）', 'IFIH1 (log-normalised)'),
    ('细胞类型', 'Cell type'),

    ('仅两类细胞可判定，置信区间均包含零',
     'Only two cell types estimable; both intervals include zero'),
    ('对照组未检出，无法判定', 'not estimable'),
    ('样本量不足，无法判定', 'too few donors'),
    ('对照组未检出\n比值无定义', 'not detected\nin controls'),
    ('对照组未检出', 'not detected in controls'),
    ('样本量不足', 'too few donors'),
    ('* p<0.05    ** p<0.01    *** p<0.001    ns 不显著',
     '* p<0.05    ** p<0.01    *** p<0.001    ns not significant'),
    ('不显著', 'not significant'),
    ('每点一位供体', 'one point per donor'),
    ('阳性率', 'positive'), ('对照', 'Control'),

    ('皮层下白质', 'Subcortical WM'), ('皮层组织块', 'Cortical block'),
    ('皮层标本', 'Cortical block'),
    ('白质', 'WM'), ('灰质', 'GM'), ('皮层', 'Cortex'),
    ('核数', 'nuclei'),
]

_NUM_UNIT = [(re.compile(r'(\d[\d ,]*)\s*例'), r'\1 donors'),
             (re.compile(r'(\d[\d ,]*)\s*核'), r'\1 nuclei')]


def tr(s):
    if not isinstance(s, str) or not s:
        return s
    for a, b in TR:
        s = s.replace(a, b)
    for rx, rep in _NUM_UNIT:
        s = rx.sub(rep, s)
    return s


def _wrap1(cls, name, argidx=0):
    orig = getattr(cls, name)
    def f(self, *a, **k):
        a = list(a)
        if len(a) > argidx:
            a[argidx] = tr(a[argidx])
        if 'label' in k: k['label'] = tr(k['label'])
        if 's' in k: k['s'] = tr(k['s'])
        return orig(self, *a, **k)
    setattr(cls, name, f)


def _wrap_list(cls, name):
    orig = getattr(cls, name)
    def f(self, labels=None, *a, **k):
        if labels is not None:
            labels = [tr(x) if isinstance(x, str) else
                      (tr(x.get_text()) if hasattr(x, 'get_text') else x) for x in labels]
        return orig(self, labels, *a, **k)
    setattr(cls, name, f)


for nm in ('set_title', 'set_xlabel', 'set_ylabel'):
    _wrap1(Axes, nm, 0)
_wrap1(Axes, 'text', 2)          # ax.text(x, y, s)
_wrap1(Figure, 'text', 2)        # fig.text(x, y, s)
_wrap1(Figure, 'suptitle', 0)
_wrap1(Colorbar, 'set_label', 0)
_wrap_list(Axes, 'set_xticklabels')
_wrap_list(Axes, 'set_yticklabels')

_line_init = Line2D.__init__
def _li(self, *a, **k):
    if 'label' in k: k['label'] = tr(k['label'])
    return _line_init(self, *a, **k)
Line2D.__init__ = _li

_ann = Axes.annotate
def _an(self, text, *a, **k):
    return _ann(self, tr(text), *a, **k)
Axes.annotate = _an

# 输出改到 figs_en/
_sf = Figure.savefig
def _save(self, fname, *a, **k):
    if isinstance(fname, str):
        fname = fname.replace('/results/mda5/figs/', '/results/mda5/figs_en/')
    return _sf(self, fname, *a, **k)
Figure.savefig = _save

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'DejaVu Sans'
print('[en_patch] 已启用英文翻译与 figs_en 输出')
