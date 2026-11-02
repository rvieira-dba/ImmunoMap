import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture
def toy_counts():
    """3 donors x 2 conditions, 2 cell populations with distinct marker genes."""
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(196)] + ["CD3E", "MS4A1", "ISG15", "MT-CO1"]
    rows, obs = [], []
    for donor in ["d1", "d2", "d3"]:
        for cond in ["ctrl", "stim"]:
            for kind in ["T", "B"]:
                lam = np.full(len(genes), 1.0)
                block = range(0, 15) if kind == "T" else range(15, 30)  # a population-specific programme
                lam[list(block)] = 15
                lam[genes.index("CD3E" if kind == "T" else "MS4A1")] = 20
                if cond == "stim":
                    lam[genes.index("ISG15")] = 25
                m = rng.poisson(lam, size=(60, len(genes)))
                rows.append(m)
                obs += [{"sample": f"{donor}_{cond}", "condition": cond, "donor": donor, "kind": kind}] * 60
    X = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
    o = pd.DataFrame(obs, index=[f"c{i}" for i in range(X.shape[0])])
    return ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=genes))
