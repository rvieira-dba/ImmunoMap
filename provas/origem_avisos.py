"""Roda o DESeq2 do teste pequeno e imprime de ONDE sai cada aviso (arquivo:linha), com as versões.
Uso: python -I provas/origem_avisos.py ."""
import importlib.metadata as m
import sys
import warnings
from pathlib import Path

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "src"))
sys.path.insert(0, str(repo / "tests"))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from immunomap.analysis import steps
from immunomap.config import Config

# mesma montagem de tests/conftest.py::toy_counts
rng = np.random.default_rng(0)
genes = [f"G{i}" for i in range(196)] + ["CD3E", "MS4A1", "ISG15", "MT-CO1"]
rows, obs = [], []
for donor in ["d1", "d2", "d3"]:
    for cond in ["ctrl", "stim"]:
        for kind in ["T", "B"]:
            lam = np.full(len(genes), 1.0)
            block = range(0, 15) if kind == "T" else range(15, 30)
            lam[list(block)] = 15
            lam[genes.index("CD3E" if kind == "T" else "MS4A1")] = 20
            if cond == "stim":
                lam[genes.index("ISG15")] = 25
            rows.append(rng.poisson(lam, size=(60, len(genes))))
            obs += [{"sample": f"{donor}_{cond}", "condition": cond, "donor": donor, "kind": kind}] * 60
X = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
a = ad.AnnData(X=X, obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(X.shape[0])]), var=pd.DataFrame(index=genes))
a.obs["cell_type"] = np.where(a.obs["kind"] == "T", "T cells", "B cells")

cfg = Config(data_dir="d", results_dir="r", figures_dir="f", min_genes_per_cell=10, max_genes_per_cell=1000,
             n_top_genes=100, n_pcs=10, min_cells_pseudobulk=5, n_jobs=1,
             marker_genes={"T cells": ["CD3E"], "B cells": ["MS4A1"]})
a.layers["counts"] = a.X.copy()
counts, meta = steps.pseudobulk(a, cfg)

print("python", sys.version.split()[0], "| pandas", m.version("pandas"), "| pydeseq2", m.version("pydeseq2"),
      "| anndata", m.version("anndata"), "| numpy", m.version("numpy"))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    try:
        de = steps.differential_expression(counts, meta, cfg)
        print("DE ok:", len(de), "linhas")
    except Exception as exc:  # noqa: BLE001
        print("DE FALHOU:", type(exc).__name__, str(exc)[:200])

# O worker devolve só o texto; para achar a origem, repete o que ele faz, sem engolir o local.
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ct = meta["cell_type"].unique()[0]
mm = meta[meta["cell_type"] == ct]
c = counts.loc[mm.index]; c = c.loc[:, c.sum(axis=0) >= 10]
md = mm[["condition", "donor"]].copy()
md["condition"] = pd.Categorical(md["condition"], categories=["ctrl", "stim"])
seen = {}
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    try:
        dds = DeseqDataSet(counts=c, metadata=md, design="~donor + condition", quiet=True, n_cpus=1)
        dds.deseq2()
        DeseqStats(dds, contrast=["condition", "stim", "ctrl"], quiet=True, n_cpus=1).summary()
        print("DESeq2 direto ok")
    except Exception as exc:  # noqa: BLE001
        import traceback
        print("DESeq2 direto FALHOU:", type(exc).__name__, str(exc)[:300])
        tb = traceback.extract_tb(exc.__traceback__)[-1]
        print("   em", tb.filename.split("site-packages/")[-1], "linha", tb.lineno, "→", tb.line)
for w in caught:
    key = (w.category.__name__, w.filename.split("site-packages/")[-1], w.lineno, str(w.message)[:70])
    seen[key] = seen.get(key, 0) + 1
for (cat, fn, ln, msg), n in seen.items():
    print(f"   {n}× {cat} em {fn}:{ln} — {msg}")
