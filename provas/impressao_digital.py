"""Impressão digital do pipeline no PBMC3k: roda os mesmos passos do ImmunoMap em 3 sementes e 2 números de componentes
e imprime um hash dos resultados. Mesma entrada + mesmo hash em outra pilha de bibliotecas = resultado idêntico.
Uso: python -I provas/impressao_digital.py . dados/pbmc3k.h5ad [--cod]   (--cod também imprime os rótulos da semente 0, em base64)"""
import base64
import hashlib
import importlib.metadata as m
import json
import sys
import warnings
import zlib
from pathlib import Path

warnings.filterwarnings("ignore")
repo, data = Path(sys.argv[1]), sys.argv[2]
sys.path.insert(0, str(repo / "src"))

import numpy as np
import pandas as pd
import scanpy as sc

from immunomap.analysis import steps
from immunomap.config import Config

sha = lambda xs: hashlib.sha256("\x1f".join(map(str, xs)).encode()).hexdigest()[:12]  # noqa: E731
raw = sc.read_h5ad(data)
raw.var_names_make_unique()
for key in ("sample", "condition", "donor"):
    raw.obs[key] = pd.Categorical(["pbmc3k"] * raw.n_obs)
X = raw.X
print("pilha:", {p: m.version(p) for p in ["scanpy", "anndata", "pandas", "numpy", "scipy", "scikit-learn", "igraph", "leidenalg", "numba"]},
      "| python", sys.version.split()[0])
print("entrada:", json.dumps({"n_obs": raw.n_obs, "n_vars": raw.n_vars, "soma_X": float(X.sum()), "obs": sha(raw.obs_names), "var": sha(raw.var_names)}))
tipos = {}
for n_pcs in (30, 10):
    for seed in range(3):
        cfg = Config(random_state=seed, integration="none", n_pcs=n_pcs)
        a, _ = steps.quality_control(raw.copy(), cfg)
        a = steps.normalise_and_reduce(a, cfg)
        steps.embed_and_cluster(a, "X_pca", cfg, compute_umap=False)
        key = steps.annotation_clusters(a, cfg)
        steps.annotate(a, cfg, groupby=key)
        lab = a.obs["cell_type"].astype(str)
        pca = np.abs(a.obsm["X_pca"])
        out = {"n_pcs": n_pcs, "seed": seed, "celulas": a.n_obs,
               "hvg": sha(a.var_names[a.var["highly_variable"]]),
               "pca_abs_5casas": hashlib.sha256(np.round(pca, 5).tobytes()).hexdigest()[:12], "pca_soma": round(float(pca.sum()), 6),
               "leiden_mapa": sha(a.obs["leiden"]), "leiden_fino": sha(a.obs[key]), "grupos": [int(a.obs["leiden"].nunique()), int(a.obs[key].nunique())],
               "rotulos": sha(lab), "contagens": lab.value_counts().sort_index().to_dict()}
        print("saida:", json.dumps(out, ensure_ascii=False), flush=True)
        if "--cod" in sys.argv and n_pcs == 30 and seed == 0:
            cats = sorted(lab.unique()); codes = np.array([cats.index(x) for x in lab], dtype=np.uint8)
            print("COD:", json.dumps({"cats": cats, "crc": zlib.crc32(codes.tobytes()), "b64": base64.b64encode(zlib.compress(codes.tobytes(), 9)).decode()}))
