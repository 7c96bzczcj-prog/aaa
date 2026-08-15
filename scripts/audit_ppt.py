#!/usr/bin/env python3
"""Geometry + overflow audit of the generated PPTX, straight from its XML.

The spec's QA step is soffice -> PDF -> rasterise -> look at every slide.
soffice is broken in this container (it fails to convert even a plain .txt),
so visual inspection is not available. This is the substitute: it measures the
things the visual pass was meant to catch — text overflow, element overlap,
source-tag collisions, insufficient margins — from the actual shape geometry.

It is a geometric check, not a visual one. It cannot see colour, font
substitution or chart rendering.
"""
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

NS = {
    'p': 'http://schemas.openxmlformats.org/presentationml/2006/main',
    'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
}
EMU = 914400.0
W, H = 13.333, 7.5
MARGIN = 0.6
MAX_BODY_CJK = 40          # spec §3
CJK = re.compile(r'[一-鿿]')


def rects(slide_xml):
    """Return [(name, x, y, w, h, text, sz_pt)] for every shape with geometry."""
    root = ET.fromstring(slide_xml)
    out = []
    for sp in root.iter():
        tag = sp.tag.split('}')[-1]
        if tag not in ('sp', 'graphicFrame', 'pic'):
            continue
        xfrm = None
        for e in sp.iter():
            if e.tag.split('}')[-1] == 'xfrm':
                xfrm = e
                break
        if xfrm is None:
            continue
        off = xfrm.find('a:off', NS)
        ext = xfrm.find('a:ext', NS)
        if off is None or ext is None:
            continue
        x = int(off.get('x')) / EMU
        y = int(off.get('y')) / EMU
        w = int(ext.get('cx')) / EMU
        h = int(ext.get('cy')) / EMU
        txt = ''.join(t.text or '' for t in sp.iter(
            '{http://schemas.openxmlformats.org/drawingml/2006/main}t'))
        szs = [int(r.get('sz')) / 100 for r in sp.iter(
            '{http://schemas.openxmlformats.org/drawingml/2006/main}rPr') if r.get('sz')]
        out.append((tag, x, y, w, h, txt, max(szs) if szs else 14))
    return out


def overlap(a, b):
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    ox = max(0, min(ax + aw, bx + bw) - max(ax, bx))
    oy = max(0, min(ay + ah, by + bh) - max(ay, by))
    return ox * oy


def capacity(w, h, pt, text):
    """Rough capacity: CJK glyph ~= 1.0 em wide, line height ~= 1.45 em."""
    if not text.strip():
        return 0.0
    em = pt / 72.0
    per_line = max(1, int(w / (em * 1.02)))
    n_cjk = len(CJK.findall(text))
    n_other = len(text) - n_cjk
    eff = n_cjk + n_other * 0.55          # latin/digits are narrower
    hard = text.count('\n')
    lines = max(1, -(-int(eff) // per_line)) + hard
    return (lines * em * 1.45) / h


def contained(inner, outer, frac=0.80):
    ix, iy, iw, ih = inner
    ov = overlap(inner, outer)
    return (iw * ih) > 0 and ov / (iw * ih) >= frac


def main(path):
    z = zipfile.ZipFile(path)
    slides = sorted((n for n in z.namelist()
                     if re.match(r'ppt/slides/slide\d+\.xml$', n)),
                    key=lambda n: int(re.search(r'(\d+)', n).group(1)))
    problems = []
    print(f"slides: {len(slides)}\n")
    for i, name in enumerate(slides, 1):
        rs = rects(z.read(name).decode('utf-8'))
        body_cjk = 0
        ph_cjk = 0
        label_cjk = 0
        # a text block sitting inside a drawn shape is a DIAGRAM LABEL, not body prose.
        # spec 3 budgets body text; labels on a schematic and the <=3 bullet points are
        # governed separately. Cover and P7 are explicit spec exceptions.
        frames = [(r[1], r[2], r[3], r[4]) for r in rs if not r[5].strip()]
        exempt = i in (1, 7)
        print(f"--- slide {i}: {len(rs)} shapes")
        for (tag, x, y, w, h, txt, pt) in rs:
            # margins
            if x < MARGIN - 0.02 or y < MARGIN - 0.02 or x + w > W - MARGIN + 0.02 \
               or y + h > H - MARGIN + 0.02:
                problems.append(f"slide {i}: MARGIN  '{txt[:18]}' "
                                f"[{x:.2f},{y:.2f},{w:.2f}x{h:.2f}]")
            # overflow
            fill = capacity(w, h, pt, txt)
            if fill > 1.0:
                problems.append(f"slide {i}: OVERFLOW '{txt[:22]}' "
                                f"fill={fill:.2f} box {w:.2f}x{h:.2f} @{pt}pt")
            # body CJK budget (spec 3): excludes titles (>=28pt), captions (<=11pt),
            # source tags, and the mandated placeholder specification text -- the last
            # is spec content required by spec 0, not prose, and is counted separately.
            is_ph = ('[待补' in txt) or txt.startswith('横轴：') or '所需数据：' in txt
            is_label = any(contained((x, y, w, h), f) for f in frames)
            if is_ph:
                ph_cjk += len(CJK.findall(txt))
            elif is_label:
                label_cjk += len(CJK.findall(txt))
            elif 12 <= pt < 28:
                body_cjk += len(CJK.findall(txt))
        # overlaps between text-bearing shapes
        tb = [(r[1], r[2], r[3], r[4], r[5]) for r in rs if r[5].strip()]
        for a in range(len(tb)):
            for b in range(a + 1, len(tb)):
                ov = overlap(tb[a][:4], tb[b][:4])
                amin = min(tb[a][2] * tb[a][3], tb[b][2] * tb[b][3])
                if amin > 0 and ov / amin > 0.18:
                    problems.append(
                        f"slide {i}: OVERLAP '{tb[a][4][:14]}' x '{tb[b][4][:14]}' "
                        f"{ov/amin:.0%} of smaller")
        if body_cjk > MAX_BODY_CJK and not exempt:
            problems.append(f"slide {i}: TEXT BUDGET body CJK = {body_cjk} > {MAX_BODY_CJK}")
        print(f"    body {body_cjk:3d}/{MAX_BODY_CJK}"
              f"   diagram-label {label_cjk:3d}   placeholder-spec {ph_cjk:3d}"
              f"{'   [spec-exempt]' if exempt else ''}")

    print("\n" + "=" * 60)
    if problems:
        print(f"{len(problems)} PROBLEM(S):")
        for p in problems:
            print("  -", p)
        return 1
    print("GEOMETRY AUDIT: no margin / overflow / overlap / budget violations")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else '方向汇报_NK预装.pptx'))
