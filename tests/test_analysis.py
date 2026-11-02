import numpy as np
import pandas as pd

from immunomap.analysis import steps
from immunomap.config import Config
from immunomap.pipeline import run


def _cfg(tmp_path, **kw):
    base = dict(data_dir=tmp_path / "d", results_dir=tmp_path / "r", figures_dir=tmp_path / "f",
                min_genes_per_cell=10, max_genes_per_cell=1000, n_top_genes=100, n_pcs=10, min_cells_pseudobulk=5,
                marker_genes={"T cells": ["CD3E"], "B cells": ["MS4A1"]})
    base.update(kw)
    return Config(**base)


def test_pseudobulk_sums_raw_counts(toy_counts, tmp_path):
    a = toy_counts.copy()
    a.layers["counts"] = a.X.copy()
    a.obs["cell_type"] = a.obs["kind"].map({"T": "T cells", "B": "B cells"})
    counts, meta = steps.pseudobulk(a, _cfg(tmp_path))
    assert len(counts) == 12  # 6 samples x 2 types
    key = "d1_ctrl|T cells"
    mask = ((a.obs["sample"] == "d1_ctrl") & (a.obs["cell_type"] == "T cells")).values
    assert counts.loc[key].sum() == int(a.X[mask].sum())
    assert meta.loc[key, "n_cells"] == 60


def test_integration_metrics_detect_mixing(toy_counts):
    a = toy_counts.copy()
    rng = np.random.default_rng(1)
    a.obs["cell_type"] = a.obs["kind"]
    mixed = rng.normal(size=(a.n_obs, 5)) + (a.obs["kind"] == "T").values[:, None] * 5
    split = mixed + pd.factorize(a.obs["sample"])[0][:, None] * 10
    a.obsm["mixed"], a.obsm["split"] = mixed, split
    good = steps.integration_metrics(a, "mixed", "sample")
    bad = steps.integration_metrics(a, "split", "sample")
    assert good["mistura_lotes"] > 0.8 > bad["mistura_lotes"]
    assert good["separacao_tipos"] > 0.6


def test_end_to_end_recovers_stimulated_gene(toy_counts, tmp_path):
    res = run(toy_counts.copy(), _cfg(tmp_path, resolution=0.3), save=True, name="toy")
    assert set(res.adata.obs["cell_type"]) <= {"T cells", "B cells", "Unassigned"}
    assert {"sem integração", "Harmony"} <= set(res.metrics)
    assert not res.de.empty
    top = res.de[res.de["padj"] < 0.05]
    assert "ISG15" in set(top["gene"])
    assert (top[top["gene"] == "ISG15"]["log2FoldChange"] > 0).all()
    assert set(res.de["desenho"]) == {"pareado"}
    assert (tmp_path / "r" / "toy" / "expressao_diferencial.csv").exists()

