"""
Readers for raw count matrices.

Supported inputs (``fmt="auto"`` detects by extension/contents):

* ``h5ad``  — AnnData file; counts taken from ``layers["counts"]``, else ``.raw``, else ``X``
* ``10x-h5`` — Cell Ranger ``filtered_feature_bc_matrix.h5``
* ``10x-mtx`` — folder (or ``.zip``/``.tar.gz``) with ``matrix.mtx(.gz)``, ``barcodes.tsv(.gz)`` and
  ``features.tsv(.gz)``/``genes.tsv``
* ``csv``/``tsv`` — dense table; by default **genes in rows, cells in columns** (``transpose=False``),
  the usual layout of GEO supplementary tables.

Every reader returns an AnnData with raw integer-like counts in ``X`` and unique gene names.
"""
from __future__ import annotations

import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path
from typing import Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

FORMATS = ("auto", "h5ad", "10x-h5", "10x-mtx", "csv", "tsv")
PathLike = Union[str, Path]


def detect_format(path: PathLike) -> str:
    p = Path(path)
    name = p.name.lower()
    if p.is_dir():
        return "10x-mtx"
    if name.endswith(".h5ad"):
        return "h5ad"
    if name.endswith(".h5"):
        return "10x-h5"
    if name.endswith((".zip", ".tar.gz", ".tgz", ".tar")):
        return "10x-mtx"
    if name.endswith((".csv", ".csv.gz")):
        return "csv"
    if name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "tsv"
    raise ValueError(f"Formato não reconhecido: {p.name}. Use .h5ad, .h5 (10x), pasta/arquivo .mtx ou .csv/.tsv.")


def _counts_from_h5ad(adata: ad.AnnData) -> ad.AnnData:
    if "counts" in adata.layers:
        out = ad.AnnData(X=adata.layers["counts"], obs=adata.obs.copy(), var=adata.var.copy())
    elif adata.raw is not None and _looks_like_counts(adata.raw.X):
        out = ad.AnnData(X=adata.raw.X, obs=adata.obs.copy(), var=adata.raw.var.copy())
    else:
        out = ad.AnnData(X=adata.X, obs=adata.obs.copy(), var=adata.var.copy())
    return out


def _looks_like_counts(X) -> bool:
    """True when a sample of the values is non-negative and integer-valued."""
    data = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return True
    sample = data[: min(data.size, 200_000)]
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def _find_mtx_dir(root: Path) -> Path:
    for candidate in [root, *sorted(p for p in root.rglob("*") if p.is_dir())]:
        names = {f.name for f in candidate.iterdir()} if candidate.is_dir() else set()
        if any(n.startswith("matrix.mtx") for n in names):
            return candidate
    raise ValueError("Nenhum matrix.mtx(.gz) encontrado na pasta/arquivo enviado.")


def _read_mtx(path: Path) -> ad.AnnData:
    import scanpy as sc

    if path.is_dir():
        return sc.read_10x_mtx(_find_mtx_dir(path), var_names="gene_symbols", cache=False)
    tmp = Path(tempfile.mkdtemp(prefix="immunomap_"))
    try:
        if path.name.lower().endswith(".zip"):
            with zipfile.ZipFile(path) as zf:
                zf.extractall(tmp)
        else:
            with tarfile.open(path) as tf:
                tf.extractall(tmp, filter="data") if hasattr(tarfile, "data_filter") else tf.extractall(tmp)
        return sc.read_10x_mtx(_find_mtx_dir(tmp), var_names="gene_symbols", cache=False)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def read_counts(path: PathLike, fmt: str = "auto", transpose: bool = False,
                sample_id: Optional[str] = None) -> ad.AnnData:
    """Read one sample (or a pre-combined file) and return raw counts."""
    import scanpy as sc

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")
    fmt = detect_format(path) if fmt == "auto" else fmt
    if fmt not in FORMATS:
        raise ValueError(f"fmt deve ser um de {FORMATS}")

    if fmt == "h5ad":
        adata = _counts_from_h5ad(sc.read_h5ad(path))
    elif fmt == "10x-h5":
        adata = sc.read_10x_h5(path)
    elif fmt == "10x-mtx":
        adata = _read_mtx(path)
    else:
        sep = "," if fmt == "csv" else "\t"
        table = pd.read_csv(path, sep=sep, index_col=0)
        if not transpose:  # genes in rows → cells in rows
            table = table.T
        adata = ad.AnnData(X=sparse.csr_matrix(table.values.astype(np.float32)),
                           obs=pd.DataFrame(index=table.index.astype(str)),
                           var=pd.DataFrame(index=table.columns.astype(str)))

    adata.var_names_make_unique()
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if not _looks_like_counts(adata.X):
        raise ValueError(
            f"{path.name}: os valores não parecem contagens brutas (há negativos ou decimais). "
            "Envie a matriz de contagens, não dados já normalizados."
        )
    if sample_id:
        adata.obs_names = [f"{sample_id}_{b}" for b in adata.obs_names]
    return adata
