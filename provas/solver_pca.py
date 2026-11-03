"""PBMC3k: quantas células mudam de rótulo se o PCA trocar o solver exato (covariance_eigh) por outro?
Uso: python -I provas/solver_pca.py . dados/pbmc3k.h5ad"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(sys.argv[1]) / "src"))
import pandas as pd
import scanpy as sc

from immunomap.analysis import steps
from immunomap.config import Config

raw = sc.read_h5ad(sys.argv[2])
raw.var_names_make_unique()
for key in ("sample", "condition", "donor"):
    raw.obs[key] = pd.Categorical(["pbmc3k"] * raw.n_obs)
orig = sc.tl.pca


def rotulos(solver, seed, n_pcs):
    def pca(*a, **k):
        k["svd_solver"] = solver
        return orig(*a, **k)

    sc.tl.pca = pca
    cfg = Config(random_state=seed, integration="none", n_pcs=n_pcs)
    a, _ = steps.quality_control(raw.copy(), cfg)
    a = steps.normalise_and_reduce(a, cfg)
    steps.embed_and_cluster(a, "X_pca", cfg, compute_umap=False)
    key = steps.annotation_clusters(a, cfg)
    steps.annotate(a, cfg, groupby=key)
    return a.obs["cell_type"].astype(str), a.obs["leiden"].nunique(), a.obs[key].nunique()


for n_pcs in (30, 10):
    for seed in range(3):
        base, g0, f0 = rotulos("covariance_eigh", seed, n_pcs)
        linha = f"{n_pcs:2} PCs · semente {seed} · covariance_eigh: {g0} grupos no mapa, {f0} finos"
        for solver in ("arpack", "randomized"):
            lab, g, f = rotulos(solver, seed, n_pcs)
            dif = int((lab != base).sum())
            linha += f" | {solver}: {dif} células com outro rótulo ({100 * dif / len(base):.2f}%), {g}/{f} grupos"
        print(linha, flush=True)
