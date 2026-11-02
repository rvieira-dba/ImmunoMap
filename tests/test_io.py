import gzip

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import io as spio, sparse

from immunomap.io import Sample, combine_samples, detect_format, load_samplesheet, read_counts


def test_detect_format(tmp_path):
    assert detect_format("a.h5ad") == "h5ad"
    assert detect_format("x/filtered_feature_bc_matrix.h5") == "10x-h5"
    assert detect_format("m.zip") == "10x-mtx"
    assert detect_format("t.csv") == "csv"
    assert detect_format(tmp_path) == "10x-mtx"
    with pytest.raises(ValueError):
        detect_format("foto.png")


def test_csv_genes_in_rows(tmp_path):
    p = tmp_path / "s.csv"
    pd.DataFrame({"cell1": [1, 0, 3], "cell2": [0, 2, 1]}, index=["A", "B", "C"]).to_csv(p)
    a = read_counts(p, sample_id="s1")
    assert a.shape == (2, 3)
    assert list(a.obs_names) == ["s1_cell1", "s1_cell2"]


def test_rejects_normalised_values(tmp_path):
    p = tmp_path / "n.csv"
    pd.DataFrame({"c1": [0.5, 1.2], "c2": [2.3, 0.1]}, index=["A", "B"]).to_csv(p)
    with pytest.raises(ValueError, match="contagens"):
        read_counts(p)


def test_mtx_folder(tmp_path):
    d = tmp_path / "mtx"
    d.mkdir()
    m = sparse.coo_matrix(np.array([[1, 0], [0, 4], [2, 2]]))  # genes x cells
    with gzip.open(d / "matrix.mtx.gz", "wb") as fh:
        spio.mmwrite(fh, m)
    with gzip.open(d / "barcodes.tsv.gz", "wt") as fh:
        fh.write("AAA-1\nCCC-1\n")
    with gzip.open(d / "features.tsv.gz", "wt") as fh:
        fh.write("E1\tCD3E\tGene Expression\nE2\tMS4A1\tGene Expression\nE3\tLYZ\tGene Expression\n")
    a = read_counts(d)
    assert a.shape == (2, 3)
    assert "CD3E" in a.var_names


def test_h5ad_prefers_counts_layer(tmp_path):
    a = ad.AnnData(X=sparse.csr_matrix(np.full((2, 2), 0.5, dtype=np.float32)))
    a.layers["counts"] = sparse.csr_matrix(np.array([[1, 2], [3, 4]], dtype=np.float32))
    p = tmp_path / "x.h5ad"
    a.write_h5ad(p)
    assert read_counts(p).X.sum() == 10


def test_samplesheet_and_combine(tmp_path):
    for name in ("a", "b"):
        pd.DataFrame({"c1": [1, 2, 0], "c2": [0, 1, 5]}, index=["G1", "G2", "G3"]).to_csv(tmp_path / f"{name}.csv")
    sheet = tmp_path / "amostras.csv"
    sheet.write_text("sample_id,condition,path,donor\nA,ctrl,a.csv,d1\nB,stim,b.csv,d1\n", encoding="utf-8")
    samples = load_samplesheet(sheet)
    assert [s.sample_id for s in samples] == ["A", "B"]
    combined = combine_samples(samples)
    assert combined.n_obs == 4
    assert set(combined.obs["condition"]) == {"ctrl", "stim"}
    assert combined.obs_names.is_unique


def test_samplesheet_requires_columns(tmp_path):
    sheet = tmp_path / "ruim.csv"
    sheet.write_text("nome,arquivo\nA,a.csv\n", encoding="utf-8")
    with pytest.raises(ValueError, match="obrigatórias"):
        load_samplesheet(sheet)
