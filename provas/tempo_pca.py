"""Quanto tempo leva cada solver do PCA no Kang 2018, e quanto o resultado difere do `arpack` (o exato de referência).
Roda o mesmo caminho da análise (QC, normalização, genes variáveis, escala) e só então compara os três solvers.
Uso (na raiz do repositório, com o ImmunoMap instalado):  python provas/tempo_pca.py dados/kang2018.h5ad"""
import sys, time, warnings; warnings.filterwarnings("ignore")
import anndata as ad, scanpy as sc, numpy as np
from immunomap.config import Config
from immunomap.analysis import steps
from immunomap import pipeline
a = ad.read_h5ad(sys.argv[1]); cfg = Config(); pipeline.validate_input(a, cfg); a, _ = steps.quality_control(a, cfg)
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="sample", flavor="seurat")
h = a[:, a.var["highly_variable"]].copy(); sc.pp.scale(h, max_value=10)
print("default solver:", sc.settings.__dict__.get("_default_svd_solver", "?"), "| X:", type(h.X).__name__, h.X.dtype)
R = {}
for solver in ["arpack", "covariance_eigh", "randomized"]:
    g = h.copy(); t = time.perf_counter(); sc.tl.pca(g, n_comps=30, random_state=0, svd_solver=solver); R[solver] = (time.perf_counter() - t, g.obsm["X_pca"], g.uns["pca"]["variance_ratio"])
    print(f"{solver:16} {R[solver][0]:5.1f} s")
ref = R["arpack"][1]
for s in ["covariance_eigh", "randomized"]:
    X = R[s][1]; sign = np.sign((X * ref).sum(0)); d = np.abs(X * sign - ref).max(axis=0)
    print(f"{s:16} máx. dif. (após alinhar sinais) por componente: até {d.max():.2e} | sinais trocados: {(sign < 0).sum()} | var. explicada dif. {np.abs(R[s][2]-R['arpack'][2]).max():.1e}")
