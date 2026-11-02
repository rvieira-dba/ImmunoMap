"""
Concordância da anotação do ImmunoMap com rótulos de especialistas no PBMC3k,
variando a semente e a resolução usada para nomear os tipos.

    python benchmarks/anotacao_pbmc3k.py                 # baixa os dados na primeira vez
    python benchmarks/anotacao_pbmc3k.py --seeds 10 --resolutions 0.6 1.5
    python benchmarks/anotacao_pbmc3k.py --n-pcs 10          # como no Cell-Atlas

Referência: ``scanpy.datasets.pbmc3k_processed()``, coluna ``louvain`` (anotada à mão no
tutorial do Scanpy). CD4 T e CD8 T viram "T cells" e Megakaryocytes vira "Platelets",
para casar com os nomes do ImmunoMap. Os passos são os mesmos do pipeline, sem o UMAP.
"""
from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

import pandas as pd
import scanpy as sc

from immunomap.analysis import steps
from immunomap.config import Config
from immunomap.datasets import load_dataset

COARSE = {"CD4 T cells": "T cells", "CD8 T cells": "T cells", "Megakaryocytes": "Platelets"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=5, help="sementes 0..N-1 (padrão 5)")
    ap.add_argument("--resolutions", type=float, nargs="+", default=[0.6, 1.0, 1.5, 2.0])
    ap.add_argument("--n-pcs", type=int, default=Config().n_pcs, help="componentes principais (padrão do ImmunoMap)")
    ap.add_argument("--data-dir", type=Path, default=Path("dados"))
    ap.add_argument("--raw", type=Path, help="h5ad bruto do PBMC3k (padrão: baixa)")
    ap.add_argument("--expert", type=Path, help="pbmc3k_processed.h5ad (padrão: baixa)")
    ap.add_argument("--out", type=Path, default=Path("resultados/benchmark_anotacao_pbmc3k.csv"))
    args = ap.parse_args(argv)

    raw = sc.read_h5ad(args.raw) if args.raw else load_dataset("pbmc3k", args.data_dir)
    raw.var_names_make_unique()
    for key in ("sample", "condition", "donor"):
        if key not in raw.obs:
            raw.obs[key] = pd.Categorical(["pbmc3k"] * raw.n_obs)
    expert = (sc.read_h5ad(args.expert) if args.expert else sc.datasets.pbmc3k_processed()).obs["louvain"]
    expert = expert.astype(str).map(lambda x: COARSE.get(x, x))

    rows = []
    for seed in range(args.seeds):
        cfg = Config(random_state=seed, integration="none", n_pcs=args.n_pcs)
        a, _ = steps.quality_control(raw.copy(), cfg)
        a = steps.normalise_and_reduce(a, cfg)
        steps.embed_and_cluster(a, "X_pca", cfg, compute_umap=False)
        common = a.obs_names.intersection(expert.index)
        for res in args.resolutions:
            key = steps.annotation_clusters(a, dataclasses.replace(cfg, annotation_resolution=res))
            steps.annotate(a, cfg, groupby=key)
            pred = a.obs.loc[common, "cell_type"].astype(str)
            row = {"semente": seed, "resolucao": res, "grupos": a.obs[key].nunique(),
                   "concordancia_%": round(100 * float((pred == expert[common]).mean()), 1),
                   "unassigned_%": round(100 * float((a.obs["cell_type"].astype(str) == "Unassigned").mean()), 1)}
            rows.append(row)
            print("  ".join(f"{k}={v}" for k, v in row.items()), flush=True)

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    summary = df.groupby("resolucao").agg(grupos=("grupos", "mean"), media_pct=("concordancia_%", "mean"),
                                          pior_pct=("concordancia_%", "min"), melhor_pct=("concordancia_%", "max"))
    print("\n" + summary.round(1).to_string())
    print(f"\nTabela completa: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
