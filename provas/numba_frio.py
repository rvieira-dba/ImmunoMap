"""Custo de compilação do numba na primeira busca de vizinhos e no primeiro UMAP, em um processo novo.
Uso: NUMBA_CACHE_DIR=$(mktemp -d) python provas/numba_frio.py   (a pasta vazia é a "frio"; rode de novo, com a mesma pasta, para ver o efeito do cache em disco)"""
import sys
import time

t0 = time.perf_counter()
import anndata as ad
import numpy as np
import scanpy as sc

t_imp = time.perf_counter() - t0
rng = np.random.default_rng(0)
a = ad.AnnData(np.zeros((6000, 5), dtype="float32"))
a.obsm["X_pca"] = rng.normal(size=(6000, 30)).astype("float32")
t = time.perf_counter(); sc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca", random_state=0); t_nb = time.perf_counter() - t
t = time.perf_counter(); sc.tl.umap(a, random_state=0); t_um = time.perf_counter() - t
print(f"python {sys.version.split()[0]} | importar {t_imp:.1f} s | vizinhos {t_nb:.1f} s | umap {t_um:.1f} s")
