"""End-to-end run used by both the CLI and the app."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Optional

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

from .analysis import plots, steps
from .config import Config
from .io.write import write_h5ad

Progress = Optional[Callable[[str], None]]


@dataclass
class Result:
    adata: ad.AnnData
    qc: pd.DataFrame
    scores: pd.DataFrame
    composition: pd.DataFrame
    metrics: Dict[str, Dict[str, float]]
    de: pd.DataFrame
    pb_meta: pd.DataFrame
    rep: str
    figures: Dict[str, Path] = field(default_factory=dict)


def normalise_design(obs: pd.DataFrame, cfg: Config) -> Optional[str]:
    """
    Make ``obs`` describe one sample per (sample, condition), in place.

    * a missing donor column becomes the sample column;
    * a column that changes from cell to cell (``nCount_RNA``, a barcode…) is refused as sample or condition;
    * a sample seen in more than one condition — typical when the "sample" column is the patient,
      as in the raw Kang 2018 file — is split into ``<sample>_<condition>``, and the original value
      becomes the donor so the design stays paired.

    Returns a note for the user when samples were split, otherwise ``None``.
    """
    s, c, d = cfg.sample_key, cfg.condition_key, cfg.donor_key
    for key in (s, c):
        if key not in obs:
            raise ValueError(f"Coluna obrigatória ausente em obs: '{key}'.")
    n = len(obs)
    for key, what, limit in ((s, "amostra", max(50, n // 5)), (c, "condição", 20)):
        k = obs[key].nunique()
        if k > limit:
            raise ValueError(
                f"A coluna usada como {what} tem {k:,} valores diferentes para {n:,} células: parece uma medida "
                f"por célula, não o nome da {what}. Escolha a coluna que identifica a {what} "
                f"(ex.: paciente, ctrl/stim).".replace(",", "."))
    if d not in obs:
        obs[d] = obs[s]
    per_sample = obs.groupby(obs[s].astype(str), observed=True)[c].nunique()
    mixed = per_sample[per_sample > 1]
    if mixed.empty:
        return None
    donor_was_sample = (obs[d].astype(str) == obs[s].astype(str)).all()
    original = obs[s].astype(str)
    before = original.nunique()
    obs[s] = pd.Categorical(original + "_" + obs[c].astype(str))
    if donor_was_sample:
        obs[d] = pd.Categorical(original)
    return (f"{len(mixed)} de {before} amostras tinham mais de uma condição; separei cada uma por condição "
            f"({before} → {obs[s].nunique()} amostras)" + (" e usei o valor original como doador." if donor_was_sample else "."))


def validate_input(adata: ad.AnnData, cfg: Config, progress: Progress = None) -> None:
    note = normalise_design(adata.obs, cfg)
    if note and progress:
        progress(note)


def run(adata: ad.AnnData, cfg: Config, progress: Progress = None, save: bool = True,
        compute_de: bool = True, name: str = "analise") -> Result:
    validate_input(adata, cfg, progress)
    adata, qc = steps.quality_control(adata, cfg, progress)
    adata = steps.normalise_and_reduce(adata, cfg, progress)

    # Before integration: same graph on raw PCA, kept for the comparison.
    multi = adata.obs[cfg.sample_key].nunique() > 1
    metrics: Dict[str, Dict[str, float]] = {}
    if multi and cfg.integration != "none":
        steps.embed_and_cluster(adata, "X_pca", cfg, progress, umap_key="X_umap_sem_integracao")

    rep = steps.integrate(adata, cfg, progress)
    steps.embed_and_cluster(adata, rep, cfg, progress)
    fine = steps.annotation_clusters(adata, cfg, progress)
    scores = steps.annotate(adata, cfg, groupby=fine, progress=progress)
    if multi:  # same final cell-type labels for both embeddings: a fair before/after
        metrics["sem integração"] = steps.integration_metrics(adata, "X_pca", cfg.sample_key, seed=cfg.random_state)
        if rep != "X_pca":
            label = {"X_harmony": "Harmony", "X_scvi": "scVI"}[rep]
            metrics[label] = steps.integration_metrics(adata, rep, cfg.sample_key, seed=cfg.random_state)

    comp = steps.composition(adata, cfg)
    counts, pb_meta = steps.pseudobulk(adata, cfg)
    de = steps.differential_expression(counts, pb_meta, cfg, progress) if compute_de else pd.DataFrame()

    res = Result(adata, qc, scores, comp, metrics, de, pb_meta, rep)
    if save:
        _save_outputs(res, cfg, name, progress)
    return res


def _save_outputs(res: Result, cfg: Config, name: str, progress: Progress) -> None:
    if progress:
        progress("Gravando tabelas, figuras e o .h5ad final")
    cfg.ensure_dirs()
    out, fig = cfg.results_dir / name, cfg.figures_dir / name
    out.mkdir(parents=True, exist_ok=True)
    res.qc.to_csv(out / "qc_por_amostra.csv", index=False)
    res.scores.to_csv(out / "escores_marcadores.csv")
    res.composition.to_csv(out / "composicao_por_amostra.csv")
    res.pb_meta.to_csv(out / "pseudobulk_amostras.csv")
    if not res.de.empty:
        res.de.to_csv(out / "expressao_diferencial.csv", index=False)
    (out / "metricas_integracao.json").write_text(json.dumps(res.metrics, indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "parametros.json").write_text(json.dumps(cfg.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")

    a = res.adata
    colors = [cfg.sample_key, cfg.condition_key, "cell_type"]
    if "X_umap_sem_integracao" in a.obsm:
        a.obsm["X_umapsem"] = a.obsm["X_umap_sem_integracao"]
        res.figures["umap_sem_integracao"] = fig / "umap_sem_integracao.png"
        plt.close(plots.umap_panels(a, colors[:2], basis="umapsem", title="Antes da integração", path=res.figures["umap_sem_integracao"]))
        del a.obsm["X_umapsem"]
    res.figures["umap"] = fig / "umap_integrado.png"
    plt.close(plots.umap_panels(a, colors, title=f"Embedding: {res.rep}", path=res.figures["umap"]))
    res.figures["composicao"] = fig / "composicao_por_amostra.png"
    plt.close(plots.composition_bars(res.composition, cfg.condition_key, path=res.figures["composicao"]))
    if not res.de.empty:
        for ct in res.de["cell_type"].unique():
            p = fig / f"volcano_{ct.replace(' ', '_').replace('+', 'pos')}.png"
            plt.close(plots.volcano(res.de, ct, path=p))
            res.figures[f"volcano:{ct}"] = p
    write_h5ad(a, cfg.data_dir / f"{name}_resultado.h5ad")
