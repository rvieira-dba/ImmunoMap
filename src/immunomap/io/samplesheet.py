"""
Sample sheets: which file is which sample, of which condition and donor.

CSV columns: ``sample_id, condition, path`` (required) and ``donor, format, transpose`` (optional).
Relative paths are resolved against the sheet's folder.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

import anndata as ad
import pandas as pd


@dataclass
class Sample:
    sample_id: str
    condition: str
    path: Path
    donor: str = ""
    fmt: str = "auto"
    transpose: bool = False


REQUIRED = ("sample_id", "condition", "path")


def load_samplesheet(path) -> List[Sample]:
    path = Path(path)
    table = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in REQUIRED if c not in table.columns]
    if missing:
        raise ValueError(f"Planilha sem as colunas obrigatórias: {', '.join(missing)}")
    if table["sample_id"].duplicated().any():
        raise ValueError("sample_id repetido na planilha: cada amostra precisa de um nome único.")
    samples = []
    for row in table.to_dict("records"):
        p = Path(row["path"])
        if not p.is_absolute():
            p = (path.parent / p).resolve()
        samples.append(Sample(
            sample_id=row["sample_id"].strip(),
            condition=row["condition"].strip(),
            path=p,
            donor=row.get("donor", "").strip(),
            fmt=(row.get("format", "") or "auto").strip(),
            transpose=str(row.get("transpose", "")).strip().lower() in {"1", "true", "sim", "yes"},
        ))
    return samples


def combine_samples(samples: Iterable[Sample], sample_key="sample", condition_key="condition",
                    donor_key="donor", progress: Optional[callable] = None) -> ad.AnnData:
    """Read each sample, annotate obs and concatenate on shared genes (inner join)."""
    from .readers import read_counts

    parts = []
    samples = list(samples)
    for i, s in enumerate(samples, 1):
        if progress:
            progress(f"Lendo {s.sample_id} ({i}/{len(samples)})")
        a = read_counts(s.path, fmt=s.fmt, transpose=s.transpose, sample_id=s.sample_id)
        a.obs[sample_key] = s.sample_id
        a.obs[condition_key] = s.condition
        a.obs[donor_key] = s.donor or s.sample_id
        parts.append(a)
    if not parts:
        raise ValueError("Nenhuma amostra informada.")
    combined = ad.concat(parts, join="inner", merge="same") if len(parts) > 1 else parts[0]
    if combined.n_vars == 0:
        raise ValueError("As amostras não têm genes em comum: confira se usam a mesma nomenclatura.")
    for key in (sample_key, condition_key, donor_key):
        combined.obs[key] = combined.obs[key].astype("category")
    return combined
