# 汇报幻灯片：肿瘤中的 NK 的流动性

`NK_流动性_方向汇报_v3_带原图.pptx` —— 在原 v2 稿（10 页）基础上，插入 6 页论文数据原图，共 16 页。
原有 10 页一字未改，新页沿用同一套版式（Microsoft YaHei / `1F3864` 深蓝标题 / `2E5C9A` 强调条 / `595959` 正文 / `A6A6A6` 出处）。

## 新增页面与图源

| 位置 | 页面 | 图源 |
|---|---|---|
| 原第 2 页后 | 原图（一）：进入肿瘤 24 h 内，效应分子同步下降 | Dean et al., *Nat Commun* 2024;15:683 — Fig. 3A / 3G / 3D / 3H |
| 同上 | 原图（二）：清除已成型肿瘤中的 NK，生长曲线不动 | 同上 — Supplementary Fig. 13A / 13D / 13F |
| 原第 3 页后 | 方法原图：光转换到底测到了什么 | Dean et al., *STAR Protoc* 2024;5:102956 — Fig. 1 与 Fig. 4C / 4E |
| 原第 4 页后 | 原图：Supplementary Fig. 6 —— “NK 出不来”的全部证据 | Dean et al., *Nat Commun* 2024;15:683 — Supplementary Fig. 6（整图） |
| 原第 5 页后 | 原图：耳部肿瘤 —— 同一分母下，T 与 NK 的出境构成 | Torcellan et al., *PNAS* 2017;114:5677 — Fig. 1A / 1F |
| 原第 8 页后 | 分子候选的原始证据：S1PR5 决定 NK / ILC1 的组织去留 | Evrard et al., *J Exp Med* 2022;219:e20210116 — Fig. 6A / 6B |

## 授权

- Dean 2024（Nat Commun）与 Dean 2024（STAR Protoc）为 **CC BY 4.0**，可在注明出处的前提下复用；每页脚注已标明。
- Torcellan 2017（PNAS）与 Evrard 2022（JEM）非 CC 授权，此处按学术引用方式截图使用并注明出处；若要对外发表需另行取得许可。

## 两处与原稿的事实校正

- **S1PR5 与 CD69**：“S1PR5 不受 CD69 调控”在 Evrard 2022 中是引用 Jenne et al. 2009，并非该文自己的数据；该文自己的证据是 Fig. 6（S1PR5 缺失改变 NK / ILC1 的组织分布）。新页已如实标注。
- **Supplementary Fig. 6 的分母**：原图 B/C/E/F 四个面板的分母是脾 / 引流淋巴结中的 NK 总数，量的是回收到的 Kaede Red⁺ NK 的 CD49a / CD11b 构成 —— 与原稿第 4 页的判断一致，新页把原图放出来供当场核对。

## 复现

```bash
pip install python-pptx pillow pymupdf
python build_figure_slides.py          # 需要 deck.pptx（v2 原稿）与 papers/crops/ 下的截图
```

`figures/` 中为从各论文原图裁出的面板截图，文件名对应上表；`build_figure_slides.py` 负责排版与插页。
