"""
Known public datasets, downloaded on demand into ``dados/`` and standardised to
``obs[sample, condition, donor]`` with raw counts in ``X``.

Downloads need internet access to the hosts listed in each entry. Sources are
the ones documented by each project; if a host changes, the error message says
which URL failed so it can be updated here.
"""
from __future__ import annotations

import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional

import anndata as ad
import pandas as pd


@dataclass(frozen=True)
class Dataset:
    key: str
    title: str
    summary: str
    question: str
    samples: str
    citation: str
    doi: str
    source: str
    loader: Callable[[Path], ad.AnnData]


def _fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_name(dest.name + ".part")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "immunomap/0.1 (+https://github.com/rubenslyra/ImmunoMap)"})
        with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
            while chunk := r.read(1 << 20):
                fh.write(chunk)
        tmp.replace(dest)
    except Exception as exc:  # noqa: BLE001
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Falha ao baixar {url}: {exc}") from exc
    return dest


def _standardise(a: ad.AnnData, sample, condition, donor) -> ad.AnnData:
    a.obs["sample"] = pd.Categorical(sample)
    a.obs["condition"] = pd.Categorical(condition)
    a.obs["donor"] = pd.Categorical(donor)
    a.var_names_make_unique()
    return a


# ---------------------------------------------------------------- loaders
def _load_pbmc3k(data_dir: Path) -> ad.AnnData:
    import scanpy as sc

    a = sc.datasets.pbmc3k()
    n = a.n_obs
    return _standardise(a, ["pbmc3k"] * n, ["saudavel"] * n, ["doador1"] * n)


def _pick(obs: pd.DataFrame, *names: str) -> pd.Series:
    for n in names:
        if n in obs.columns:
            return obs[n].astype(str)
    raise KeyError(f"Nenhuma das colunas {names} encontrada em obs: {list(obs.columns)}")


KANG_URL = "https://exampledata.scverse.org/pertpy/kang_2018.h5ad"  # o mesmo arquivo do pertpy.data.kang_2018()


def _load_kang2018(data_dir: Path) -> ad.AnnData:
    path = _fetch(KANG_URL, Path(data_dir) / "brutos" / "kang_2018.h5ad")
    a = ad.read_h5ad(path)
    from .io.readers import _counts_from_h5ad

    obs = a.obs.copy()
    a = _counts_from_h5ad(a)
    condition = _pick(obs, "label", "stim", "condition")
    donor = _pick(obs, "replicate", "ind", "patient", "donor")
    sample = donor + "_" + condition
    a.obs = obs
    if "cell_type" in obs.columns:
        a.obs["author_cell_type"] = obs["cell_type"].astype(str)
    return _standardise(a, sample.values, condition.values, donor.values)


_10X = {
    "pbmc8k_v2": "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz",
    "pbmc10k_v3": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5",
}


def _load_pbmc_v2_v3(data_dir: Path) -> ad.AnnData:
    import scanpy as sc
    from .io.readers import read_counts

    raw = data_dir / "brutos"
    v2 = _fetch(_10X["pbmc8k_v2"], raw / "pbmc8k_filtered_gene_bc_matrices.tar.gz")
    v3 = _fetch(_10X["pbmc10k_v3"], raw / "pbmc_10k_v3_filtered_feature_bc_matrix.h5")
    a2 = read_counts(v2, sample_id="v2")
    a3 = read_counts(v3, sample_id="v3")
    a2.obs["sample"], a3.obs["sample"] = "chromium_v2", "chromium_v3"
    combined = ad.concat([a2, a3], join="inner")
    n2, n3 = a2.n_obs, a3.n_obs
    return _standardise(combined, ["chromium_v2"] * n2 + ["chromium_v3"] * n3,
                        ["saudavel"] * (n2 + n3), ["doador_v2"] * n2 + ["doador_v3"] * n3)


REGISTRY: Dict[str, Dataset] = {
    "pbmc3k": Dataset(
        key="pbmc3k",
        title="PBMC3k (10x Genomics, 2016)",
        summary="~2.700 PBMC de um doador saudável. Uma amostra só: o caso de referência do Cell-Atlas.",
        question="O pipeline reencontra os tipos imunes conhecidos?",
        samples="1 amostra",
        citation="Zheng, G. X. Y. et al. Nat. Commun. 8, 14049 (2017).",
        doi="10.1038/ncomms14049",
        source="scanpy.datasets.pbmc3k()",
        loader=_load_pbmc3k,
    ),
    "kang2018": Dataset(
        key="kang2018",
        title="Kang et al. 2018 — lúpus, interferon-β",
        summary="PBMC de 8 pacientes com lúpus, cada um em controle e estimulado com IFN-β: 16 amostras pareadas.",
        question="Como cada tipo celular responde ao interferon? (integração + pseudobulk pareado por paciente)",
        samples="16 amostras (8 pacientes × 2 condições)",
        citation="Kang, H. M. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat. Biotechnol. 36, 89 (2018).",
        doi="10.1038/nbt.4042",
        source=KANG_URL + " (cópia do scverse/pertpy) · GEO GSE96583",
        loader=_load_kang2018,
    ),
    "pbmc_v2_v3": Dataset(
        key="pbmc_v2_v3",
        title="PBMC em duas químicas 10x (v2 × v3)",
        summary="O mesmo tipo de amostra sequenciado com Chromium v2 (8k células) e v3 (10k células).",
        question="Quanto da diferença entre amostras é técnica (efeito de lote)?",
        samples="2 amostras (uma por tecnologia)",
        citation="10x Genomics, conjuntos públicos pbmc8k (v2) e pbmc_10k_v3.",
        doi="",
        source=" · ".join(_10X.values()),
        loader=_load_pbmc_v2_v3,
    ),
}


def load_dataset(key: str, data_dir: Path, progress: Optional[Callable[[str], None]] = None) -> ad.AnnData:
    """Return the dataset, reusing ``dados/<key>.h5ad`` when it already exists."""
    from .io.write import write_h5ad

    if key not in REGISTRY:
        raise KeyError(f"Conjunto desconhecido: {key}. Opções: {', '.join(REGISTRY)}")
    cached = Path(data_dir) / f"{key}.h5ad"
    if cached.exists():
        if progress:
            progress(f"Usando cópia local {cached}")
        return ad.read_h5ad(cached)
    if progress:
        progress(f"Baixando {REGISTRY[key].title}…")
    a = REGISTRY[key].loader(Path(data_dir))
    write_h5ad(a, cached)
    return a
