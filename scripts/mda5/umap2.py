"""GSE180759 UMAP，分块处理，避免上一次的内存崩溃。

上一版把 3.1 GB 稠密矩阵交给 AnnData，normalize / log1p / HVG 各复制一份，
峰值超过物理内存。这一版用 mmap 读入，两遍流式处理：
  第一遍：逐块算 log 归一化后的每基因均值与方差
  第二遍：只把方差最大的 2000 个基因抽成一个小矩阵
之后的 PCA / 近邻 / UMAP 都在 66432 x 2000 上做，约 0.5 GB。
"""
import numpy as np, pandas as pd, time, gc, csv

t0 = time.time()
def log(*a): print(*a, '|', round(time.time() - t0), 's', flush=True)

X = np.load('umap_X.npy', mmap_mode='r')          # cells x genes, 原始计数
genes = np.array([str(g) for g in np.load('umap_genes.npy', allow_pickle=True)])
NC, NG = X.shape
log('matrix', X.shape)

CH = 4000
# ---- 第一遍：流式均值 / 方差（log1p(CP10k) 尺度）----
s1 = np.zeros(NG, np.float64); s2 = np.zeros(NG, np.float64)
for i in range(0, NC, CH):
    B = np.array(X[i:i + CH], dtype=np.float32, copy=True)
    tot = B.sum(axis=1, keepdims=True); tot[tot == 0] = 1
    np.multiply(B, 1e4 / tot, out=B); np.log1p(B, out=B)
    s1 += B.sum(axis=0); s2 += (B.astype(np.float64) ** 2).sum(axis=0)
    del B; gc.collect()
    if (i // CH) % 5 == 0: log('pass1', i)
mu = s1 / NC
var = np.maximum(s2 / NC - mu ** 2, 0)

NTOP = 2000
hv = np.argsort(var)[::-1][:NTOP]
# IFIH1 单独保留，便于特征图直接用同一矩阵
ifih1_col = int(np.where(genes == 'IFIH1')[0][0])
log('HVG selected; IFIH1 in HVG:', ifih1_col in set(hv.tolist()))

# ---- 第二遍：抽取 HVG 子矩阵 ----
hv_sorted = np.sort(hv)
Y = np.empty((NC, NTOP), np.float32)
ifih1 = np.empty(NC, np.float32)
for i in range(0, NC, CH):
    B = np.array(X[i:i + CH], dtype=np.float32, copy=True)
    tot = B.sum(axis=1, keepdims=True); tot[tot == 0] = 1
    np.multiply(B, 1e4 / tot, out=B); np.log1p(B, out=B)
    Y[i:i + B.shape[0]] = B[:, hv_sorted]
    ifih1[i:i + B.shape[0]] = B[:, ifih1_col]
    del B; gc.collect()
    if (i // CH) % 5 == 0: log('pass2', i)
del X; gc.collect()
log('subset', Y.shape)

# ---- 标准化并降维 ----
Y -= Y.mean(axis=0)
sd = Y.std(axis=0); sd[sd == 0] = 1
Y /= sd
np.clip(Y, -10, 10, out=Y)

from sklearn.decomposition import PCA
pca = PCA(n_components=50, svd_solver='randomized', random_state=0)
P = pca.fit_transform(Y).astype(np.float32)
del Y; gc.collect()
np.save('umap_pca.npy', P)
log('PCA done', P.shape)

import scanpy as sc, anndata as ad
a = ad.AnnData(X=np.zeros((NC, 1), np.float32))
a.obsm['X_pca'] = P[:, :40]
sc.pp.neighbors(a, n_neighbors=15, use_rep='X_pca')
log('neighbors done')
sc.tl.umap(a, min_dist=0.3, random_state=0)
U = a.obsm['X_umap']
np.save('umap_coords.npy', U)
log('UMAP done')

# ---- 元数据对齐 ----
hdr = open('header.txt').read().strip().split(',')
ann = {r['nucleus_barcode']: r for r in csv.DictReader(open('GSE180759_annotation.txt'), delimiter='\t')}
obs = pd.DataFrame({
    'barcode': hdr,
    'cell_type': [ann[b]['cell_type'] for b in hdr],
    'pathology': [ann[b]['pathology'] for b in hdr],
    'donor': [ann[b]['NBB_case'] for b in hdr],
    'UMAP1': U[:, 0], 'UMAP2': U[:, 1],
    'IFIH1_lognorm': ifih1,
})
obs.to_csv('umap_obs.csv', index=False)
log('SAVED')
print(obs.pathology.value_counts().to_string())
print(obs.cell_type.value_counts().to_string())
