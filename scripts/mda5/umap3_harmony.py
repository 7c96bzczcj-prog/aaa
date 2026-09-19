"""对 GSE180759 的 PCA 做 Harmony 批次整合后重算 UMAP。

未整合的版本里，文库 18 单独形成一个少突胶质孤岛（5748 个核中 5166 个来自该文库），
其 IFIH1 阳性率 4.37% 对其余少突胶质的 1.29%，测序深度也更高。
按分期铺开的特征图会把这个单文库批次显示成"斑块周围白质"的生物学信号。
因此以文库为批次变量做整合。
"""
import numpy as np, pandas as pd, time
import scanpy as sc, anndata as ad
import harmonypy

t0 = time.time()
def log(*a): print(*a, '|', round(time.time() - t0), 's', flush=True)

P = np.load('umap_pca.npy')
obs = pd.read_csv('umap_obs.csv')
d = pd.read_parquet('df_180759.parquet').reset_index(drop=True)
assert len(d) == len(obs)
meta = pd.DataFrame({'libr': d.libr.astype(str).values, 'donor': d.donor.astype(str).values})
log('PCA', P.shape, '文库数', meta.libr.nunique())

ho = harmonypy.run_harmony(P[:, :40], meta, ['libr'], max_iter_harmony=20)
H = np.asarray(ho.Z_corr)
if H.shape[0] != P.shape[0]:          # 这一版返回 PCs x cells
    H = H.T
H = np.ascontiguousarray(H, dtype=np.float32)
assert H.shape[0] == P.shape[0], H.shape
np.save('umap_harmony.npy', H)
log('harmony done', H.shape)

a = ad.AnnData(X=np.zeros((H.shape[0], 1), np.float32))
a.obsm['X_h'] = H
sc.pp.neighbors(a, n_neighbors=15, use_rep='X_h')
log('neighbors done')
sc.tl.umap(a, min_dist=0.3, random_state=0)
U = a.obsm['X_umap']
log('umap done')

obs['UMAP1'] = U[:, 0]; obs['UMAP2'] = U[:, 1]
obs['libr'] = meta.libr.values
obs.to_csv('umap_obs_harmony.csv', index=False)
log('SAVED')

# 整合是否奏效：原来那团里文库 18 的占比
oli = obs[obs.cell_type == 'oligodendrocytes']
for lb, g in oli.groupby('libr'):
    if len(g) < 500: continue
    print(f'  文库 {lb:>3}  n={len(g):>5}  UMAP1 中位 {g.UMAP1.median():6.2f}  UMAP2 中位 {g.UMAP2.median():6.2f}')
