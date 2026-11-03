"""Quais são os genes 'significativos' na divisão 4 do teste nulo (pseudobulk DESeq2)? Mesmo sorteio de calibracao_nula.py (semente 1).
Uso: python -I provas/calibracao_nula_divisao4.py dados/kang2018_resultado.h5ad"""
import itertools, sys, warnings
warnings.filterwarnings("ignore")
import anndata as ad, numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

CT = "CD14+ Monocytes"
a = ad.read_h5ad(sys.argv[1])
c = a[(a.obs["condition"] == "ctrl") & (a.obs["cell_type"] == CT)].copy()
amostras = sorted(c.obs["sample"].unique())
divisoes = [s for s in itertools.combinations(amostras, 4) if amostras[0] in s]
escolhidas = [divisoes[i] for i in np.random.default_rng(1).choice(len(divisoes), size=6, replace=False)]
A = escolhidas[3]
cont = c.layers["counts"]
pb = pd.DataFrame({s: np.asarray(cont[(c.obs["sample"] == s).to_numpy()].sum(axis=0)).ravel() for s in amostras}, index=c.var_names).T.astype(int)
pb = pb.loc[:, pb.sum(axis=0) >= 10]
md = pd.DataFrame({"grp": pd.Categorical(["A" if s in A else "B" for s in pb.index], categories=["B", "A"])}, index=pb.index)
dds = DeseqDataSet(counts=pb, metadata=md, design="~grp", quiet=True, n_cpus=1); dds.deseq2()
st = DeseqStats(dds, contrast=["grp", "A", "B"], quiet=True, n_cpus=1); st.summary()
d = st.results_df
sig = d[(d.padj < 0.05)].sort_values("padj")
print("grupo A:", A)
print("grupo B:", [s for s in amostras if s not in A])
print(f"{len(sig)} genes com padj<0,05; {int((sig.log2FoldChange.abs() > 1).sum())} com |log2FC|>1")
sexo = {"XIST", "TSIX", "RPS4Y1", "DDX3Y", "UTY", "ZFY", "KDM5D", "EIF1AY", "TXLNGY", "USP9Y", "NLGN4Y", "LINC00278", "TTTY15", "PRKY", "TMSB4Y"}
print("ligados ao sexo (cromossomos X/Y) entre eles:", sorted(set(sig.index) & sexo))
print(sig[["baseMean", "log2FoldChange", "padj"]].round(4).to_string())
