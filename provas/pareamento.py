"""Monócitos CD14+ do Kang 2018: DESeq2 pareado (`~donor + condition`) contra não pareado (`~condition`), com os mesmos 16 pseudobulks.
Mostra quantos genes mudam de classificação só por trocar o desenho (padj < 0,05 e |log2FC| > 1).
Uso (depois de rodar o Kang 2018): python -I provas/pareamento.py . dados/kang2018_resultado.h5ad"""
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(sys.argv[1]) / "src"))

import anndata as ad
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from immunomap.analysis import steps
from immunomap.config import Config

CT = "CD14+ Monocytes"
a = ad.read_h5ad(sys.argv[2])
counts, meta = steps.pseudobulk(a, Config())
m = meta[meta["cell_type"] == CT]
c = counts.loc[m.index]
c = c.loc[:, c.sum(axis=0) >= 10]
md = m[["condition", "donor"]].copy()
md["condition"] = pd.Categorical(md["condition"], categories=["ctrl", "stim"])
print(f"{CT}: {len(m)} pseudobulks ({m['donor'].nunique()} doadores), {c.shape[1]} genes com pelo menos 10 contagens", flush=True)


def rodar(design):
    t = time.perf_counter()
    dds = DeseqDataSet(counts=c, metadata=md, design=design, quiet=True, n_cpus=1)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "stim", "ctrl"], quiet=True, n_cpus=1)
    st.summary()
    return st.results_df, time.perf_counter() - t


res = {}
for nome, design in (("pareado     ~donor + condition", "~donor + condition"), ("não pareado ~condition", "~condition")):
    r, dt = rodar(design)
    res[nome] = r
    sig = (r.padj < 0.05) & (r.log2FoldChange.abs() > 1)
    print(f"{nome:30} {dt:5.1f} s · genes testados {r.padj.notna().sum()} · padj<0,05 e |log2FC|>1: {int(sig.sum())} · "
          f"mediana do erro padrão do log2FC {r.lfcSE.median():.3f} · ISG15 log2FC {r.loc['ISG15', 'log2FoldChange']:.2f} (padj {r.loc['ISG15', 'padj']:.1e})")
p, n = list(res.values())
sp, sn = (p.padj < 0.05) & (p.log2FoldChange.abs() > 1), (n.padj < 0.05) & (n.log2FoldChange.abs() > 1)
print(f"só no pareado: {int((sp & ~sn).sum())} · só no não pareado: {int((sn & ~sp).sum())} · nos dois: {int((sp & sn).sum())}")
