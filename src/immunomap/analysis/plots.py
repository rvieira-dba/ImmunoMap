"""Figures saved as PNG (also returned as matplotlib figures for the app)."""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402


def _save(fig, path: Optional[Path]):
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
    return fig


def umap_panels(adata: ad.AnnData, colors: List[str], basis: str = "umap", title: str = "",
                path: Optional[Path] = None):
    n = len(colors)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.6), squeeze=False)
    for ax, color in zip(axes[0], colors):
        sc.pl.embedding(adata, basis=basis, color=color, ax=ax, show=False, frameon=False,
                        legend_fontsize=7, size=max(120000 / max(adata.n_obs, 1), 2))
    if title:
        fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return _save(fig, path)


def composition_bars(comp: pd.DataFrame, condition_key: str, path: Optional[Path] = None):
    data = comp.drop(columns=[condition_key]).sort_index()
    order = comp.sort_values(condition_key).index
    data = data.loc[order]
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(data) + 3), 4.2))
    bottom = np.zeros(len(data))
    cmap = plt.get_cmap("tab10")
    for i, col in enumerate(data.columns):
        ax.bar(range(len(data)), data[col].values, bottom=bottom, label=col, color=cmap(i % 10), width=0.85)
        bottom += data[col].values
    ax.set_xticks(range(len(data)), [f"{s}\n({comp.loc[s, condition_key]})" for s in data.index],
                  rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("% das células da amostra")
    ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.set_title("Composição de tipos celulares por amostra")
    fig.tight_layout()
    return _save(fig, path)


def volcano(de: pd.DataFrame, cell_type: str, path: Optional[Path] = None, top: int = 10):
    d = de[de["cell_type"] == cell_type].dropna(subset=["padj", "log2FoldChange"])
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    y = -np.log10(d["padj"].clip(lower=1e-300))
    sig = (d["padj"] < 0.05) & (d["log2FoldChange"].abs() > 1)
    ax.scatter(d["log2FoldChange"], y, s=6, c=np.where(sig, "#c2477a", "#9aa4ae"), linewidths=0)
    for _, r in d[sig].nsmallest(top, "padj").iterrows():
        ax.annotate(r["gene"], (r["log2FoldChange"], -np.log10(max(r["padj"], 1e-300))), fontsize=7)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#666")
    ax.axvline(-1, ls=":", lw=0.8, c="#666"); ax.axvline(1, ls=":", lw=0.8, c="#666")
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10(p ajustado)")
    ax.set_title(f"{cell_type} — {int(sig.sum())} genes com |log2FC| > 1 e padj < 0,05", fontsize=9)
    fig.tight_layout()
    return _save(fig, path)
