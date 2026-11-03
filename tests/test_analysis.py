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



def test_cell_types_named_on_finer_partition(toy_counts, tmp_path):
    res = run(toy_counts.copy(), _cfg(tmp_path, resolution=0.3, annotation_resolution=2.0, integration="none"),
              save=False, compute_de=False)
    obs = res.adata.obs
    assert "leiden_anotacao" in obs
    assert obs["leiden_anotacao"].nunique() >= obs["leiden"].nunique()
    # one name per fine group, and the score table is indexed by the fine groups
    assert (obs.groupby("leiden_anotacao", observed=True)["cell_type"].nunique() == 1).all()
    assert res.scores.index.name == "leiden_anotacao"


def test_annotation_resolution_must_be_positive(tmp_path):
    import pytest

    with pytest.raises(ValueError):
        _cfg(tmp_path, annotation_resolution=0)


def test_sample_column_with_two_conditions_is_split(toy_counts, tmp_path):
    from immunomap.pipeline import normalise_design

    a = toy_counts.copy()
    a.obs["sample"] = a.obs["donor"].astype(str)  # "sample" = patient, as in the raw Kang file
    del a.obs["donor"]
    note = normalise_design(a.obs, _cfg(tmp_path))
    assert note and "separei" in note
    assert a.obs["sample"].nunique() == 6
    assert (a.obs.groupby("sample", observed=True)["condition"].nunique() == 1).all()
    assert a.obs["donor"].nunique() == 3  # the patient becomes the donor: paired design kept


def test_per_cell_column_refused_as_sample(toy_counts, tmp_path):
    import pytest
    from immunomap.pipeline import normalise_design

    a = toy_counts.copy()
    a.obs["sample"] = [f"c{i}" for i in range(a.n_obs)]  # e.g. a barcode or nCount_RNA
    with pytest.raises(ValueError, match="medida"):
        normalise_design(a.obs, _cfg(tmp_path))


def test_parallel_de_matches_one_at_a_time(toy_counts, tmp_path):
    a = toy_counts.copy()
    cfg1 = _cfg(tmp_path, resolution=0.3, integration="none", n_jobs=1)
    res = run(a, cfg1, save=False)
    counts, meta = steps.pseudobulk(res.adata, cfg1)
    serial = steps.differential_expression(counts, meta, cfg1)
    parallel = steps.differential_expression(counts, meta, _cfg(tmp_path, n_jobs=2))
    cols = ["cell_type", "gene", "log2FoldChange", "padj", "desenho"]
    pd.testing.assert_frame_equal(serial[cols].reset_index(drop=True), parallel[cols].reset_index(drop=True))


def test_n_jobs_must_not_be_negative(tmp_path):
    import pytest

    with pytest.raises(ValueError):
        _cfg(tmp_path, n_jobs=-1)
