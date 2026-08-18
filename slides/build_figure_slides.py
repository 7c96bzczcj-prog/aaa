#!/usr/bin/env python3
"""在原 PPT 中插入论文数据原图（截图）页面。"""
import copy
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR
from pptx.oxml.ns import qn
from PIL import Image

SRC = 'deck.pptx'
OUT = 'deck_out.pptx'
CROPS = 'figures/'

SW = 13.3333            # 幻灯片宽（英寸）
NAVY = RGBColor(0x1F, 0x38, 0x64)
BLUE = RGBColor(0x2E, 0x5C, 0x9A)
BODY = RGBColor(0x59, 0x59, 0x59)
MUTE = RGBColor(0xA6, 0xA6, 0xA6)
LINE = RGBColor(0xD9, 0xD9, 0xD9)
FONT = 'Microsoft YaHei'


def set_font(run, size, color, bold=False):
    f = run.font
    f.size = Pt(size)
    f.bold = bold
    f.color.rgb = color
    f.name = FONT
    rPr = f._rPr
    for tag in ('a:ea', 'a:cs'):
        el = rPr.makeelement(qn(tag), {'typeface': FONT})
        rPr.append(el)


def textbox(slide, x, y, w, h, parts, anchor=MSO_ANCHOR.MIDDLE, line_pts=None):
    """parts: [(text, size, color, bold), ...] —— 同一段内的多个 run。"""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    p = tf.paragraphs[0]
    if line_pts:
        pPr = p._p.get_or_add_pPr()
        lnSpc = pPr.makeelement(qn('a:lnSpc'), {})
        spcPts = pPr.makeelement(qn('a:spcPts'), {'val': str(int(line_pts * 100))})
        lnSpc.append(spcPts)
        pPr.insert(0, lnSpc)
    for text, size, color, bold in parts:
        r = p.add_run()
        r.text = text
        set_font(r, size, color, bold)
    return tb


def title(slide, text):
    textbox(slide, 0.70, 0.45, 11.90, 0.80, [(text, 26, NAVY, True)])


def cite(slide, text):
    textbox(slide, 0.70, 6.80, 11.90, 0.30, [(text, 10, MUTE, False)])


def lead(slide, text, y=1.28, h=0.34, size=12.5):
    textbox(slide, 0.70, y, 11.90, h, [(text, size, BODY, False)], line_pts=17)


def accent(slide, x, y, h=0.30):
    """细竖条，沿用原 deck 的强调条。"""
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(0.06), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = BLUE
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def aspect(path):
    with Image.open(CROPS + path) as im:
        return im.width / im.height


def layout_row(paths, avail_w, max_h, gap, x_center=SW / 2):
    """等高排布一行图片，返回 [(path, x, y_offset_w, w, h)] 的宽高与起始 x。"""
    ars = [aspect(p) for p in paths]
    n = len(paths)
    h = min(max_h, (avail_w - gap * (n - 1)) / sum(ars))
    widths = [h * a for a in ars]
    total = sum(widths) + gap * (n - 1)
    x = x_center - total / 2
    out = []
    for p, w in zip(paths, widths):
        out.append((p, x, w, h))
        x += w + gap
    return out


def picture(slide, path, x, y, w, h):
    pic = slide.shapes.add_picture(CROPS + path, Inches(x), Inches(y),
                                   Inches(w), Inches(h))
    pic.line.color.rgb = LINE
    pic.line.width = Pt(0.75)
    return pic


def panel_label(slide, x, y, w, text, sub=None):
    parts = [(text, 12, NAVY, True)]
    if sub:
        parts.append(('　' + sub, 11.5, BODY, False))
    textbox(slide, x, y, w, 0.28, parts, anchor=MSO_ANCHOR.BOTTOM)


def new_slide(prs):
    s = prs.slides.add_slide(prs.slide_layouts[0])
    # 显式白底，和原有页面一致
    cSld = s._element.find(qn('p:cSld'))
    bg = cSld.makeelement(qn('p:bg'), {})
    bgPr = bg.makeelement(qn('p:bgPr'), {})
    fill = bgPr.makeelement(qn('a:solidFill'), {})
    clr = fill.makeelement(qn('a:srgbClr'), {'val': 'FFFFFF'})
    fill.append(clr)
    bgPr.append(fill)
    bgPr.append(bgPr.makeelement(qn('a:effectLst'), {}))
    bg.append(bgPr)
    cSld.insert(0, bg)
    return s


# ---------------------------------------------------------------- 各页内容

def slide_fig3(prs):
    s = new_slide(prs)
    title(s, '原图（一）：进入肿瘤 24 h 内，效应分子同步下降')
    lead(s, '每个点是一只小鼠。绿 = 标记之后才进入肿瘤（Green 24hrs，"新兵"）；'
            '红 = 标记时已在瘤内（Red 24hrs，"老兵"）。同一张图里比的是同一群细胞的先后两个状态。')
    rows = [
        [('f3_A.png', 'A · CCL5', '招募 DC 的趋化因子'),
         ('f3_G.png', 'G · 穿孔素 Perforin', '杀伤机器的核心组件')],
        [('f3_D.png', 'D · 颗粒酶 A', '与穿孔素同向下降'),
         ('f3_H.png', 'H · CD107a', '脱颗粒 —— 实际动手的读数')],
    ]
    y = 1.72
    for row in rows:
        placed = layout_row([r[0] for r in row], 11.2, 1.94, 0.70)
        for (path, x, w, h), (_, lab, sub) in zip(placed, row):
            panel_label(s, x, y, w, lab, sub)
            picture(s, path, x, y + 0.30, w, h)
        y += 0.30 + placed[0][3] + 0.36
    cite(s, 'Dean et al., Nat Commun 2024;15:683 — Fig. 3A / 3G / 3D / 3H（CC BY 4.0，原图截取）')


def slide_supp13(prs):
    s = new_slide(prs)
    title(s, '原图（二）：清除已成型肿瘤中的 NK，生长曲线不动')
    lead(s, '第 6 天起注射 anti-NK1.1（曲线上方三角为给药时点）。清除效率经流式确认，'
            '但两条生长曲线在 MC38 与 B16-F10 中都完全重合 —— 在缺乏 T / B 细胞的 RAG 小鼠中同样如此。')
    row = [('s13_A.png', 'A · MC38', '两条曲线重合'),
           ('s13_D.png', 'D · 清除效率', '瘤内 NK 已清干净'),
           ('s13_F.png', 'F · B16-F10', 'WT 与 RAG 都无差别')]
    placed = layout_row([r[0] for r in row], 12.0, 3.35, 0.50)
    y = 1.86
    for (path, x, w, h), (_, lab, sub) in zip(placed, row):
        panel_label(s, x, y, w, lab, sub)
        picture(s, path, x, y + 0.30, w, h)
    ny = y + 0.30 + placed[0][3] + 0.40
    accent(s, 0.70, ny, 0.62)
    textbox(s, 1.00, ny, 11.30, 0.62,
            [('三角 = 给药时点；A 为 mean ± SEM，D / F 为逐鼠散点 + mean ± SD。'
              'NK 被清得很干净（D），曲线却一点没动 —— 说明限制变量不在"瘤内 NK 够不够多"。',
              12.5, BODY, False)], anchor=MSO_ANCHOR.MIDDLE, line_pts=17)
    cite(s, 'Dean et al., Nat Commun 2024;15:683 — Supplementary Fig. 13A / 13D / 13F（CC BY 4.0，原图截取）')


def slide_method(prs):
    s = new_slide(prs)
    title(s, '方法原图：光转换到底测到了什么')
    lead(s, '① 405 nm 照射后，瘤内细胞由绿变红且不可逆；② 照射当时瘤内 98.4% 转红，'
            '引流淋巴结只有 0.35%（本底）；③ 因此，此后在淋巴结里出现的红色细胞只可能来自肿瘤。')
    r1 = layout_row(['star_f1.png'], 11.0, 2.00, 0)
    y = 1.72
    for path, x, w, h in r1:
        panel_label(s, x, y, w, '① 光转换原理')
        picture(s, path, x, y + 0.30, w, h)
    y2 = y + 0.30 + r1[0][3] + 0.34
    row = [('star_f4C.png', '② 照射当时（0 h）的实测'),
           ('star_f4E.png', '③ 红 / 绿的判读规则')]
    placed = layout_row([r[0] for r in row], 11.6, 1.76, 0.85)
    for (path, x, w, h), (_, lab) in zip(placed, row):
        panel_label(s, x, y2, w, lab)
        picture(s, path, x, y2 + 0.30, w, h)
    cite(s, 'Dean et al., STAR Protoc 2024;5:102956 — Fig. 1 与 Fig. 4C / 4E（CC BY 4.0，原图截取）')


def slide_supp6(prs):
    s = new_slide(prs)
    title(s, '原图：Supplementary Fig. 6 —— "NK 出不来"的全部证据')
    placed = layout_row(['s6_full.png'], 6.30, 5.05, 0, x_center=0.70 + 6.30 / 2)
    path, x, w, h = placed[0]
    picture(s, path, x, 1.42, w, h)
    tx = x + w + 0.55
    tw = SW - 0.70 - tx
    blocks = [
        ('A / D：代表性流式图',
         '全篇唯一直接反映"出境量"的读数，只有一只代表性动物。红门内比例：脾 0.1% / 0.4%，'
         '引流淋巴结 1.5% / 1.7% —— 没有逐鼠散点，也没有统计检验。'),
        ('B / C / E / F：做了统计的四个面板',
         '量的是回收到的红色 NK 中 CD49a / CD11b 的构成比例，不是出境量本身；'
         '分母还是淋巴结里的 NK 总数，而非瘤内被光转换的红色起始池。'),
        ('全图没有同管的 T 细胞对照',
         '抗体面板里含 CD3 / NKp46，T 细胞被染上了，但在设门时排除掉了。'),
    ]
    ty = 1.55
    for head, sub in blocks:
        accent(s, tx, ty, 1.42)
        textbox(s, tx + 0.24, ty, tw - 0.24, 0.34, [(head, 13, NAVY, True)],
                anchor=MSO_ANCHOR.TOP)
        textbox(s, tx + 0.24, ty + 0.40, tw - 0.24, 1.02, [(sub, 11.5, BODY, False)],
                anchor=MSO_ANCHOR.TOP, line_pts=16)
        ty += 1.75
    cite(s, 'Dean et al., Nat Commun 2024;15:683 — Supplementary Fig. 6（CC BY 4.0，原图截取；'
            '图注与坐标轴为原文）')


def slide_pnas(prs):
    s = new_slide(prs)
    title(s, '原图：耳部肿瘤 —— 同一分母下，T 与 NK 的出境构成')
    lead(s, '同一批动物、同一个分母。左边是瘤内构成，右边是走到引流淋巴结的红色细胞构成 —— '
            '正文只报告了 T 细胞的 44%，NK 那一格从未被提及。')
    row = [('pnas_A.png', 'A · 瘤内构成', '分母＝瘤内 Kaede 阳性细胞'),
           ('pnas_F.png', 'F · 出境到 dLN 的构成', 'T≈44%，NK≈6.5%')]
    placed = layout_row([r[0] for r in row], 11.2, 3.35, 0.80)
    y = 1.86
    for (path, x, w, h), (_, lab, sub) in zip(placed, row):
        panel_label(s, x, y, w, lab, sub)
        picture(s, path, x, y + 0.30, w, h)
    ny = y + 0.30 + placed[0][3] + 0.34
    accent(s, 0.70, ny, 0.66)
    textbox(s, 1.00, ny, 11.30, 0.66,
            [('NK 那一列在 2017 年就印在图上：柱高约 6.5%，8 只小鼠全部大于 0。'
              '本页百分数为图上目测（±2–3 pp）。', 12.5, BODY, False)],
            anchor=MSO_ANCHOR.MIDDLE, line_pts=17)
    cite(s, 'Torcellan et al., PNAS 2017;114:5677 — Fig. 1A / 1F（原图截取）')


def slide_jem(prs):
    s = new_slide(prs)
    title(s, '分子候选的原始证据：S1PR5 决定 NK / ILC1 的组织去留')
    lead(s, 'T 细胞靠 CD69 降解 S1PR1 留在组织里；NK 出境走的是 S1PR5。'
            '"S1PR5 不与 CD69 结合"一句在该文中是引用 Jenne et al. 2009，本文自己的数据是下面这两张。')
    placed = layout_row(['jem_AB.png'], 11.4, 4.05, 0)
    path, x, w, h = placed[0]
    picture(s, path, x, 1.86, w, h)
    ly = 1.86 + h + 0.20
    half = w / 2 - 0.15
    textbox(s, x, ly, half, 0.62,
            [('A　', 12, NAVY, True),
             ('循环型 cNK 高表达 S1pr1 / S1pr5；组织驻留的 ILC1 则低表达。',
              11.5, BODY, False)], anchor=MSO_ANCHOR.TOP, line_pts=16)
    textbox(s, x + w / 2 + 0.15, ly, half, 0.62,
            [('B　', 12, NAVY, True),
             ('WT 与 S1pr5 缺失细胞的骨髓嵌合体 —— 缺 S1PR5 的细胞在组织里过度蓄积'
              '（SI-IEL、唾液腺 ILC1），在血液中减少。',
              11.5, BODY, False)], anchor=MSO_ANCHOR.TOP, line_pts=16)
    cite(s, 'Evrard et al., J Exp Med 2022;219:e20210116 — Fig. 6A / 6B（该文未使用肿瘤模型）')


# ---------------------------------------------------------------- 组装

def main():
    prs = Presentation(SRC)
    n_before = len(prs.slides._sldIdLst)

    builders = [
        (2, slide_fig3),      # 插到原第 2 页之后
        (2, slide_supp13),
        (3, slide_method),    # 原第 3 页之后
        (4, slide_supp6),
        (5, slide_pnas),
        (8, slide_jem),
    ]
    # 先全部追加到末尾，记录新页在 sldIdLst 中的元素
    added = []
    for after, fn in builders:
        fn(prs)
        added.append((after, list(prs.slides._sldIdLst)[-1]))

    # 按"原页码"重新插入：从后往前处理可避免下标漂移
    lst = prs.slides._sldIdLst
    orig = list(lst)[:n_before]
    for after, el in reversed(added):
        lst.remove(el)
        anchor = orig[after - 1]
        anchor.addnext(el)

    prs.save(OUT)
    print('saved', OUT, 'slides =', len(list(lst)))


if __name__ == '__main__':
    main()
