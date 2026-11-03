"""
The multi-sample workflow, one function per step.

Raw counts are kept in ``layers["counts"]`` from the start: pseudobulk needs them,
and normalisation happens on a copy in ``X``.
"""
from __future__ import annotations

import logging
import warnings
from typing import Callable, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ..config import Config

Progress = Optional[Callable[[str], None]]
UNASSIGNED = "Unassigned"


def _say(progress: Progress, msg: str) -> None:
    if progress:
        progress(msg)


# ------------------------------------------------------------------ QC
def quality_control(adata: ad.AnnData, cfg: Config, progress: Progress = None) -> Tuple[ad.AnnData, pd.DataFrame]:
    """Per-sample QC with fixed thresholds; returns the filtered data and a per-sample report."""
    _say(progress, "Controle de qualidade por amostra")
    adata.var["mt"] = np.asarray(adata.var_names.str.upper().str.startswith("MT-"), dtype=bool)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    has_mt = bool(adata.var["mt"].any())

    keep = (adata.obs["n_genes_by_counts"] >= cfg.min_genes_per_cell) & \
           (adata.obs["n_genes_by_counts"] <= cfg.max_genes_per_cell)
    if has_mt:
        keep &= adata.obs["pct_counts_mt"] <= cfg.max_mito_percent

    key = cfg.sample_key
    before = adata.obs.groupby(key, observed=True).size().rename("celulas_antes")
    after = adata.obs[keep].groupby(key, observed=True).size().rename("celulas_depois")
    report = pd.concat([before, after], axis=1).fillna(0).astype(int)
    report["mantidas_%"] = (100 * report["celulas_depois"] / report["celulas_antes"]).round(1)
    report["mediana_genes"] = adata.obs[keep].groupby(key, observed=True)["n_genes_by_counts"].median().round(0)
    if has_mt:
        report["mediana_mito_%"] = adata.obs[keep].groupby(key, observed=True)["pct_counts_mt"].median().round(2)

    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=cfg.min_cells_per_gene)
    report.attrs["has_mt"] = has_mt
    return adata, report.reset_index()


# ------------------------------------------------------------------ normalise / reduce
def normalise_and_reduce(adata: ad.AnnData, cfg: Config, progress: Progress = None) -> ad.AnnData:
    _say(progress, "Normalização, genes variáveis e PCA")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata  # full log-normalised transcriptome for markers and plots
    n_batches = adata.obs[cfg.sample_key].nunique()
    batch_key = cfg.sample_key if n_batches > 1 else None
    sc.pp.highly_variable_genes(adata, n_top_genes=min(cfg.n_top_genes, adata.n_vars - 1),
                                batch_key=batch_key, flavor="seurat")
    hvg = adata[:, adata.var["highly_variable"]].copy()
    with warnings.catch_warnings():
        # scale() zero-centres, so a sparse matrix becomes dense here by design; the warning is noise.
        warnings.filterwarnings("ignore", message="zero-centering a sparse array/matrix densifies it")
        sc.pp.scale(hvg, max_value=10)
    n_pcs = min(cfg.n_pcs, hvg.n_vars - 1, hvg.n_obs - 1)
    sc.tl.pca(hvg, n_comps=n_pcs, random_state=cfg.random_state)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]
    return adata


# ------------------------------------------------------------------ integration
def integrate(adata: ad.AnnData, cfg: Config, progress: Progress = None) -> str:
    """Return the name of the embedding to use downstream."""
    n_batches = adata.obs[cfg.sample_key].nunique()
    if cfg.integration == "none" or n_batches < 2:
        return "X_pca"
    if cfg.integration == "harmony":
        _say(progress, "Integração com Harmony")
        import harmonypy

        log = logging.getLogger("harmonypy")  # 0.1.x logs at INFO even with verbose=False
        level = log.level
        log.setLevel(logging.WARNING)
        try:
            ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [cfg.sample_key],
                                       random_state=cfg.random_state, verbose=False)
        finally:
            log.setLevel(level)
        Z = np.asarray(ho.Z_corr)
        adata.obsm["X_harmony"] = Z.T if Z.shape[0] != adata.n_obs else Z
        return "X_harmony"
    if cfg.integration == "scvi":
        _say(progress, "Integração com scVI (pode levar alguns minutos)")
        try:
            import scvi
        except ImportError as exc:
            raise RuntimeError("scVI não instalado: pip install -e \".[scvi]\"") from exc
        hv = adata[:, adata.var["highly_variable"]].copy()
        scvi.settings.seed = cfg.random_state
        scvi.model.SCVI.setup_anndata(hv, layer="counts", batch_key=cfg.sample_key)
        model = scvi.model.SCVI(hv)
        model.train(max_epochs=min(400, max(50, 20000 * 400 // max(hv.n_obs, 1) // 100)))
        adata.obsm["X_scvi"] = model.get_latent_representation()
        return "X_scvi"
    raise ValueError(cfg.integration)


def embed_and_cluster(adata: ad.AnnData, rep: str, cfg: Config, progress: Progress = None,
                      umap_key: str = "X_umap", compute_umap: bool = True) -> None:
    _say(progress, f"Vizinhança, UMAP e Leiden ({rep})")
    sc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, use_rep=rep, random_state=cfg.random_state)
    if compute_umap:
        sc.tl.umap(adata, random_state=cfg.random_state)
        if umap_key != "X_umap":
            adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, resolution=cfg.resolution, key_added="leiden", flavor="igraph",
                 n_iterations=2, random_state=cfg.random_state)


def annotation_clusters(adata: ad.AnnData, cfg: Config, progress: Progress = None,
                        key: str = "leiden_anotacao") -> str:
    """Finer Leiden partition on the current neighbour graph, used only to name cell types."""
    _say(progress, f"Grupos finos para a anotação (resolução {cfg.annotation_resolution:g})")
    sc.tl.leiden(adata, resolution=cfg.annotation_resolution, key_added=key, flavor="igraph",
                 n_iterations=2, random_state=cfg.random_state)
    return key


# ------------------------------------------------------------------ annotation
def _sorted_groups(groups) -> List[str]:
    groups = [str(g) for g in groups]
    try:
        return sorted(groups, key=int)
    except ValueError:
        return sorted(groups)


def annotate(adata: ad.AnnData, cfg: Config, groupby: str = "leiden", progress: Progress = None) -> pd.DataFrame:
    """Marker-score annotation (same rule as Cell-Atlas). Returns the per-cluster score table."""
    _say(progress, "Anotação por marcadores clássicos")
    source = adata.raw if adata.raw is not None else adata
    clusters = adata.obs[groupby].astype(str)
    order = _sorted_groups(clusters.unique())
    scores = {}
    for cell_type, genes in cfg.marker_genes.items():
        genes = [g for g in genes if g in source.var_names]
        if not genes:
            continue
        m = source[:, genes].X
        m = m.toarray() if sparse.issparse(m) else np.asarray(m)
        means = pd.DataFrame(m, columns=genes, index=adata.obs_names).groupby(clusters.values).mean().reindex(order)
        span = (means.max() - means.min()).replace(0, np.nan)
        scores[cell_type] = ((means - means.min()) / span).fillna(0.0).mean(axis=1)
    table = pd.DataFrame(scores, index=order)
    labels = {}
    for cl, row in table.iterrows():
        best = row.idxmax() if len(row) else UNASSIGNED
        labels[cl] = best if len(row) and row[best] >= cfg.min_annotation_score else UNASSIGNED
    adata.obs["cell_type"] = pd.Categorical(clusters.map(labels))
    table.index.name = groupby
    return table


# ------------------------------------------------------------------ metrics
def integration_metrics(adata: ad.AnnData, rep: str, batch_key: str, label_key: str = "cell_type",
                        k: int = 30, max_cells: int = 5000, seed: int = 0) -> Dict[str, float]:
    """
    Two numbers in [0, 1]:

    * ``mistura_lotes`` — mean normalised entropy of sample labels among each cell's k nearest
      neighbours (1 = samples perfectly mixed, 0 = each neighbourhood from one sample);
    * ``separacao_tipos`` — silhouette of cell types rescaled to [0, 1] (higher = types stay apart).
    """
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import NearestNeighbors

    X = np.asarray(adata.obsm[rep])
    batches = adata.obs[batch_key].astype(str).values
    labels = adata.obs[label_key].astype(str).values if label_key in adata.obs else None
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=min(max_cells, adata.n_obs), replace=False)
    nb = len(set(batches))
    out = {"mistura_lotes": float("nan"), "separacao_tipos": float("nan")}
    if nb > 1:
        nn = NearestNeighbors(n_neighbors=min(k + 1, adata.n_obs)).fit(X)
        _, ind = nn.kneighbors(X[idx])
        ent = []
        for row in ind[:, 1:]:
            _, counts = np.unique(batches[row], return_counts=True)
            p = counts / counts.sum()
            ent.append(float(-(p * np.log(p)).sum() / np.log(nb)))
        out["mistura_lotes"] = float(np.mean(ent))
    if labels is not None and len(set(labels[idx])) > 1:
        out["separacao_tipos"] = float((silhouette_score(X[idx], labels[idx]) + 1) / 2)
    return out


def composition(adata: ad.AnnData, cfg: Config) -> pd.DataFrame:
    """Share of each cell type per sample (%), with condition for plotting."""
    t = pd.crosstab(adata.obs[cfg.sample_key], adata.obs["cell_type"], normalize="index") * 100
    cond = adata.obs.groupby(cfg.sample_key, observed=True)[cfg.condition_key].first()
    t.insert(0, cfg.condition_key, cond.reindex(t.index).astype(str))
    return t.round(2)


# ------------------------------------------------------------------ pseudobulk
def pseudobulk(adata: ad.AnnData, cfg: Config, layer: str = "counts") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Sum raw counts per (sample, cell type). Returns ``(counts, metadata)`` where rows are
    ``sample|cell_type`` and metadata has sample, condition, donor, cell_type and n_cells.
    Groups with fewer than ``min_cells_pseudobulk`` cells are dropped.
    """
    X = adata.layers[layer] if layer in adata.layers else adata.X
    X = sparse.csr_matrix(X)
    keys = adata.obs[cfg.sample_key].astype(str) + "|" + adata.obs["cell_type"].astype(str)
    codes, uniques = pd.factorize(keys)
    indicator = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                                  shape=(len(uniques), len(codes)))
    summed = indicator @ X
    counts = pd.DataFrame(np.asarray(summed.todense()).astype(np.int64), index=uniques, columns=adata.var_names)
    meta = (adata.obs.assign(_k=keys.values)
            .groupby("_k", observed=True)
            .agg(sample=(cfg.sample_key, "first"), condition=(cfg.condition_key, "first"),
                 donor=(cfg.donor_key, "first") if cfg.donor_key in adata.obs else (cfg.sample_key, "first"),
                 cell_type=("cell_type", "first"), n_cells=("cell_type", "size"))
            .reindex(uniques))
    for c in ("sample", "condition", "donor", "cell_type"):
        meta[c] = meta[c].astype(str)
    keep = meta["n_cells"] >= cfg.min_cells_pseudobulk
    return counts[keep.values], meta[keep]


def _report_deseq_warnings(caught, cell_type: str, progress: Progress) -> None:
    """PyDESeq2 warnings as one short note per message (or re-raised as warnings without a progress sink)."""
    seen = set()
    for w in caught:
        msg = str(w.message)
        if msg in seen:
            continue
        seen.add(msg)
        if progress is None:
            warnings.warn(msg, w.category, stacklevel=2)
        elif "did not converge" in msg:
            progress(f"   nota ({cell_type}): a curva de dispersão não convergiu; o DESeq2 usou a média (comum em grupos pequenos)")
        elif "residual degrees of freedom is less than 3" in msg:
            progress(f"   nota ({cell_type}): poucas réplicas para o desenho; a dispersão fica menos precisa")
        else:
            progress(f"   nota do DESeq2 ({cell_type}): {msg}")


def differential_expression(counts: pd.DataFrame, meta: pd.DataFrame, cfg: Config,
                            progress: Progress = None) -> pd.DataFrame:
    """
    DESeq2 (PyDESeq2) per cell type on pseudobulk counts. Paired design
    ``~donor + condition`` when every donor appears in both conditions.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    conditions = sorted(meta["condition"].unique())
    if len(conditions) < 2:
        return pd.DataFrame()
    ref = cfg.reference_condition or conditions[0]
    others = [c for c in conditions if c != ref]
    results = []
    for cell_type in sorted(meta["cell_type"].unique()):
        if cell_type == UNASSIGNED:
            continue
        m = meta[meta["cell_type"] == cell_type]
        per_cond = m.groupby("condition").size()
        if (per_cond.reindex(conditions).fillna(0) < 2).any():
            continue  # need ≥2 replicates per condition
        c = counts.loc[m.index]
        c = c.loc[:, c.sum(axis=0) >= 10]
        paired = m.groupby("donor")["condition"].nunique().eq(len(conditions)).all() and m["donor"].nunique() > 1
        design = "~donor + condition" if paired else "~condition"
        _say(progress, f"Expressão diferencial: {cell_type} ({'pareado' if paired else 'não pareado'})")
        md = m[["condition", "donor"]].copy()
        md["condition"] = pd.Categorical(md["condition"], categories=[ref] + others)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dds = DeseqDataSet(counts=c, metadata=md, design=design, quiet=True, n_cpus=1)
            dds.deseq2()
            stats = {}
            for other in others:
                st = DeseqStats(dds, contrast=["condition", other, ref], quiet=True, n_cpus=1)
                st.summary()
                stats[other] = st
        _report_deseq_warnings(caught, cell_type, progress)
        for other, st in stats.items():
            r = st.results_df.copy()
            r.insert(0, "gene", r.index)
            r.insert(0, "contraste", f"{other} vs {ref}")
            r.insert(0, "cell_type", cell_type)
            r["desenho"] = "pareado" if paired else "não pareado"
            results.append(r.reset_index(drop=True))
    if not results:
        return pd.DataFrame()
    out = pd.concat(results, ignore_index=True)
    return out.sort_values(["cell_type", "padj"], na_position="last")
