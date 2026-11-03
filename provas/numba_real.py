"""Igual ao numba_frio.py, mas com o PCA de verdade do Kang 2018 (24.562 células x 30 componentes) lido do resultado da análise.
Uso: NUMBA_CACHE_DIR=$(mktemp -d) python provas/numba_real.py dados/kang2018_resultado.h5ad   (de novo, com a mesma pasta, para ver o efeito do cache em disco)"""
import sys
import time

t0 = time.perf_counter()
import anndata as ad
import numpy as np
import scanpy as sc

t_imp = time.perf_counter() - t0
X = ad.read_h5ad(sys.argv[1], backed="r").obsm["X_pca"][:]
a = ad.AnnData(np.zeros((X.shape[0], 2), dtype="float32"))
a.obsm["X_pca"] = np.ascontiguousarray(X)
t = time.perf_counter(); sc.pp.neighbors(a, use_rep="X_pca", random_state=0); t_nb = time.perf_counter() - t
t = time.perf_counter(); sc.tl.umap(a, random_state=0); t_um = time.perf_counter() - t
print(f"python {sys.version.split()[0]} | {X.shape[0]} células, X_pca {X.dtype} | importar {t_imp:.1f} s | vizinhos {t_nb:.1f} s | umap {t_um:.1f} s")
