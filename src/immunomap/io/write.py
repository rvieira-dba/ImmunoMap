"""Atomic .h5ad writing (a crash never leaves a half-written file)."""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import pandas as pd


def _plain_strings(frame: pd.DataFrame) -> None:
    for c in frame.columns:
        if isinstance(frame[c].dtype, pd.StringDtype):
            frame[c] = frame[c].astype(object)
    if isinstance(frame.index.dtype, pd.StringDtype):
        frame.index = frame.index.astype(object)


def write_h5ad(adata: ad.AnnData, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _plain_strings(adata.obs)
    _plain_strings(adata.var)
    tmp = path.with_name(path.name + ".tmp")
    try:
        adata.write_h5ad(tmp)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return path
