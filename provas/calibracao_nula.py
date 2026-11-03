"""Calibração sob hipótese nula: só amostras de CONTROLE do Kang 2018 (não há efeito de condição nenhum).
Sorteia 6 divisões das 8 amostras em 4 × 4 e testa 'A × B' em monócitos CD14+ de dois jeitos:
  (1) célula a célula (Wilcoxon, cada célula é uma réplica)  (2) pseudobulk + DESeq2 (cada amostra é uma réplica, o desenho do ImmunoMap)
Como não há efeito verdadeiro, todo gene 'significativo' é falso positivo.  Uso (depois de rodar o Kang 2018): python -I provas/calibracao_nula.py dados/kang2018_resultado.h5ad"""
import itertools
import sys
import time
import warnings

warnings.filterwarnings("ignore")
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed

CT, N_DIV = "CD14+ Monocytes", 6
a = ad.read_h5ad(sys.argv[1])
c = a[(a.obs["condition"] == "ctrl") & (a.obs["cell_type"] == CT)].copy()
del a
amostras = sorted(c.obs["sample"].unique())
print(f"{CT}: {c.n_obs} células de controle em {len(amostras)} amostras ({c.obs['sample'].value_counts().min()} a {c.obs['sample'].value_counts().max()} células por amostra)", flush=True)
divisoes = [s for s in itertools.combinations(amostras, 4) if amostras[0] in s]
escolhidas = [divisoes[i] for i in np.random.default_rng(1).choice(len(divisoes), size=N_DIV, replace=False)]
cont = c.layers["counts"]
pb = pd.DataFrame({s: np.asarray(cont[(c.obs["sample"] == s).to_numpy()].sum(axis=0)).ravel() for s in amostras}, index=c.var_names).T.astype(int)
pb = pb.loc[:, pb.sum(axis=0) >= 10]


def uma(i, A):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    warnings.filterwarnings("ignore")
    grp = c.obs["sample"].map(lambda s: "A" if s in A else "B").astype("category")
    cc = c.copy()
    cc.obs["grp"] = grp
    t = time.perf_counter()
    sc.tl.rank_genes_groups(cc, "grp", groups=["A"], reference="B", method="wilcoxon", use_raw=False)
    w = sc.get.rank_genes_groups_df(cc, group="A")
    t_w = time.perf_counter() - t
    md = pd.DataFrame({"grp": pd.Categorical(["A" if s in A else "B" for s in pb.index], categories=["B", "A"])}, index=pb.index)
    t = time.perf_counter()
    dds = DeseqDataSet(counts=pb, metadata=md, design="~grp", quiet=True, n_cpus=1)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["grp", "A", "B"], quiet=True, n_cpus=1)
    st.summary()
    d = st.results_df
    t_d = time.perf_counter() - t
    f = lambda p, l: (int((p < 0.05).sum()), int(((p < 0.05) & (l.abs() > 1)).sum()))  # noqa: E731
    return i, A, f(w["pvals_adj"], w["logfoldchanges"]), f(d["padj"], d["log2FoldChange"]), t_w, t_d


res = Parallel(n_jobs=2)(delayed(uma)(i, A) for i, A in enumerate(escolhidas))
print(f"{'divisão':8} {'amostras do grupo A':34} {'célula a célula: padj<0,05 | e |log2FC|>1':>46} {'pseudobulk DESeq2: padj<0,05 | e |log2FC|>1':>46}")
for i, A, w, d, tw, td in res:
    print(f"{i + 1:<8} {', '.join(x.replace('patient_', '').replace('_ctrl', '') for x in A):34} {w[0]:>22} | {w[1]:<21} {d[0]:>22} | {d[1]:<21}")
mw = np.mean([r[2][1] for r in res]); md_ = np.mean([r[3][1] for r in res])
print(f"média de falsos positivos com |log2FC|>1: célula a célula {mw:.0f} · pseudobulk {md_:.1f}  (genes testados: {pb.shape[1]})")
