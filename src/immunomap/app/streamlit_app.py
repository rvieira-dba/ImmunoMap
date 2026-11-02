"""
ImmunoMap — interface no navegador (mobile-first: uma coluna, passos em sequência).

    immunomap app          # neste computador
    immunomap app --lan    # também pelo tablet/celular na mesma rede: http://<ip-do-pc>:8501
"""
from __future__ import annotations

import io
import sys
import tempfile
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

SRC = Path(__file__).resolve().parents[2]
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from immunomap import __version__  # noqa: E402
from immunomap.analysis import plots  # noqa: E402
from immunomap.config import INTEGRATIONS, Config  # noqa: E402
from immunomap.datasets import REGISTRY, load_dataset  # noqa: E402
from immunomap.io import Sample, combine_samples, read_counts  # noqa: E402
from immunomap.pipeline import run  # noqa: E402

st.set_page_config(page_title="ImmunoMap", layout="centered", initial_sidebar_state="collapsed")
st.markdown(
    """<style>
    .block-container {padding-top: 1.2rem; padding-bottom: 3rem; max-width: 760px;}
    div.stButton > button, div.stDownloadButton > button {min-height: 48px; width: 100%;}
    [data-testid="stMetricValue"] {font-size: 1.6rem;}
    </style>""",
    unsafe_allow_html=True,
)

INTEGRATION_LABELS = {"harmony": "Harmony (recomendado)", "none": "Sem integração", "scvi": "scVI (precisa de scvi-tools)"}
UPLOAD_DIR = Path(tempfile.gettempdir()) / "immunomap_uploads"

st.title("ImmunoMap")
st.caption(f"v{__version__} · scRNA-seq com várias amostras: QC → integração → tipos celulares → expressão diferencial")

# ------------------------------------------------------------------ 1. dados
st.header("1. Dados")
source = st.radio("De onde vêm os dados?", ["Conjunto conhecido", "Enviar arquivos (uma amostra por arquivo)", "Arquivo já combinado (.h5ad)"],
                  label_visibility="visible")

data_dir = Path(st.session_state.get("data_dir", "dados"))


@st.cache_resource(show_spinner=False)
def _known(key: str, data_dir: str):
    return load_dataset(key, Path(data_dir))


def _save_upload(f) -> Path:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    dest = UPLOAD_DIR / f.name
    dest.write_bytes(f.getbuffer())
    return dest


adata, run_name = None, None
if source == "Conjunto conhecido":
    key = st.selectbox("Conjunto", list(REGISTRY), format_func=lambda k: REGISTRY[k].title)
    d = REGISTRY[key]
    with st.container(border=True):
        st.markdown(f"**{d.samples}.** {d.summary}")
        st.markdown(f"**Pergunta:** {d.question}")
        st.caption(f"{d.citation}" + (f" · doi:{d.doi}" if d.doi else "") + f" · Fonte: {d.source}")
    if st.button("Carregar conjunto", type="primary", use_container_width=True):
        with st.status("Carregando (o primeiro uso baixa os dados)…", expanded=True) as s:
            try:
                st.session_state["input"] = _known(key, str(data_dir))
                st.session_state["input_name"] = key
                s.update(label="Conjunto carregado", state="complete")
            except Exception as exc:  # noqa: BLE001
                s.update(label="Falha ao carregar", state="error")
                st.error(str(exc))

elif source.startswith("Enviar arquivos"):
    st.write("Formatos: `.h5ad`, `.h5` da 10x, pasta `.mtx` compactada em `.zip`/`.tar.gz`, ou tabela `.csv`/`.tsv` (genes nas linhas).")
    files = st.file_uploader("Arquivos das amostras", accept_multiple_files=True,
                             type=["h5ad", "h5", "zip", "gz", "tgz", "tar", "csv", "tsv", "txt"])
    if files:
        st.write("Diga qual amostra é cada arquivo. Use o **mesmo doador** para amostras pareadas (ex.: antes e depois do tratamento).")
        default = pd.DataFrame({
            "arquivo": [f.name for f in files],
            "amostra": [Path(f.name).name.split(".")[0] for f in files],
            "condicao": ["controle"] * len(files),
            "doador": [""] * len(files),
        })
        table = st.data_editor(default, hide_index=True, use_container_width=True, disabled=["arquivo"],
                               key=f"sheet_{len(files)}")
        if st.button("Ler arquivos", type="primary", use_container_width=True):
            if table["amostra"].duplicated().any():
                st.error("Cada amostra precisa de um nome único.")
            else:
                with st.status("Lendo arquivos…", expanded=True) as s:
                    try:
                        by_name = {f.name: f for f in files}
                        samples = [Sample(sample_id=r.amostra.strip(), condition=r.condicao.strip(),
                                          donor=(r.doador or "").strip(), path=_save_upload(by_name[r.arquivo]))
                                   for r in table.itertuples()]
                        st.session_state["input"] = combine_samples(samples, progress=s.write)
                        st.session_state["input_name"] = "envio"
                        s.update(label="Arquivos lidos", state="complete")
                    except Exception as exc:  # noqa: BLE001
                        s.update(label="Falha ao ler", state="error")
                        st.error(str(exc))

else:
    f = st.file_uploader("Arquivo .h5ad com todas as amostras", type=["h5ad"])
    if f is not None:
        import scanpy as sc

        path = _save_upload(f)
        meta = sc.read_h5ad(path, backed="r").obs
        cols = list(meta.columns)
        if not cols:
            st.error("O arquivo não tem colunas em obs para indicar amostra e condição.")
        else:
            def _guess(options, names, fallback=0):
                low = [o.lower() for o in options]
                for n in names:
                    if n in low:
                        return low.index(n)
                return fallback

            sample_col = st.selectbox("Coluna da amostra", cols, index=_guess(cols, ["sample", "amostra", "sample_id", "orig.ident"]))
            condition_col = st.selectbox("Coluna da condição", cols,
                                         index=_guess(cols, ["condition", "condicao", "label", "stim", "group", "grupo"], min(1, len(cols) - 1)))
            donor_opts = ["(nenhuma)"] + cols
            donor_col = st.selectbox("Coluna do doador (opcional)", donor_opts,
                                     index=_guess(donor_opts, ["donor", "doador", "replicate", "patient", "paciente", "ind"]))
            if st.button("Usar este arquivo", type="primary", use_container_width=True):
                try:
                    a = read_counts(path)
                    a.obs = meta.reindex(a.obs_names).copy()
                    a.obs["sample"] = a.obs[sample_col].astype(str).astype("category")
                    a.obs["condition"] = a.obs[condition_col].astype(str).astype("category")
                    a.obs["donor"] = (a.obs[donor_col] if donor_col != "(nenhuma)" else a.obs["sample"]).astype(str).astype("category")
                    st.session_state["input"], st.session_state["input_name"] = a, Path(f.name).stem
                except Exception as exc:  # noqa: BLE001
                    st.error(str(exc))

adata = st.session_state.get("input")
if adata is not None:
    obs = adata.obs
    c1, c2, c3 = st.columns(3)
    c1.metric("Células", f"{adata.n_obs:,}".replace(",", "."))
    c2.metric("Amostras", obs["sample"].nunique())
    c3.metric("Condições", obs["condition"].nunique())
    with st.expander("Amostras carregadas"):
        st.dataframe(obs.groupby("sample", observed=True).agg(condicao=("condition", "first"), doador=("donor", "first"),
                                                              celulas=("condition", "size")), use_container_width=True)

# ------------------------------------------------------------------ 2. parâmetros
st.header("2. Parâmetros")
with st.expander("Ajustar (os padrões servem para PBMC)", expanded=False):
    integration = st.selectbox("Integração entre amostras", INTEGRATIONS, index=INTEGRATIONS.index("harmony"),
                               format_func=INTEGRATION_LABELS.get,
                               help="Com uma amostra só, a integração é pulada automaticamente.")
    resolution = st.slider("Resolução do Leiden (zoom dos grupos do mapa)", 0.1, 2.0, 0.6, 0.1)
    annotation_resolution = st.slider(
        "Resolução da anotação", 0.5, 3.0, 1.5, 0.1,
        help="Grupos mais finos, usados só para dar nome aos tipos celulares. Evita que dois tipos "
             "parecidos (ex.: T CD8 e NK) caiam no mesmo grupo e recebam um nome só.")
    n_pcs = st.slider("Componentes principais", 10, 50, 30, 5)
    min_genes = st.number_input("Mínimo de genes por célula", 50, 2000, 200, 50)
    max_genes = st.number_input("Máximo de genes por célula", 1000, 15000, 6000, 500)
    max_mito = st.number_input("Máximo de % mitocondrial", 1.0, 100.0, 15.0, 1.0)
    seed = st.number_input("Semente", 0, 10_000, 0, 1)
    do_de = st.toggle("Expressão diferencial entre condições (pseudobulk)", value=True)
    st.caption("A expressão diferencial precisa de pelo menos 2 amostras por condição em cada tipo celular.")

# ------------------------------------------------------------------ 3. rodar
st.header("3. Rodar")
if st.button("Rodar análise", type="primary", disabled=adata is None, use_container_width=True):
    cfg = Config(integration=integration, resolution=float(resolution),
                 annotation_resolution=float(annotation_resolution), n_pcs=int(n_pcs),
                 min_genes_per_cell=int(min_genes), max_genes_per_cell=int(max_genes),
                 max_mito_percent=float(max_mito), random_state=int(seed))
    with st.status("Analisando…", expanded=True) as s:
        try:
            res = run(adata.copy(), cfg, progress=s.write, compute_de=do_de,
                      name=st.session_state.get("input_name", "analise"))
            st.session_state["result"], st.session_state["cfg"] = res, cfg
            s.update(label="Análise concluída", state="complete", expanded=False)
        except Exception as exc:  # noqa: BLE001
            s.update(label="A análise falhou", state="error")
            st.error(str(exc))
if adata is None:
    st.caption("Carregue os dados no passo 1 para liberar este botão.")

# ------------------------------------------------------------------ 4. resultados
res = st.session_state.get("result")
if res is not None:
    cfg = st.session_state["cfg"]
    a = res.adata
    st.header("4. Resultados")
    tabs = st.tabs(["Resumo", "Mapa", "Composição", "Diferencial", "QC", "Baixar"])

    with tabs[0]:
        c1, c2 = st.columns(2)
        c1.metric("Células após QC", f"{a.n_obs:,}".replace(",", "."))
        c2.metric("Tipos celulares", a.obs["cell_type"].nunique())
        if res.metrics:
            st.subheader("Integração")
            st.dataframe(pd.DataFrame(res.metrics).T.rename(columns={"mistura_lotes": "mistura de amostras",
                                                                     "separacao_tipos": "separação de tipos"}).round(2),
                         use_container_width=True)
            st.caption("Ambos vão de 0 a 1. Boa integração aumenta a mistura de amostras sem derrubar a separação de tipos.")
        if not res.de.empty:
            sig = res.de[(res.de["padj"] < 0.05) & (res.de["log2FoldChange"].abs() > 1)]
            st.subheader("Genes diferenciais por tipo")
            st.dataframe(sig.groupby("cell_type").size().rename("genes (padj < 0,05; |log2FC| > 1)"), use_container_width=True)

    with tabs[1]:
        color = st.selectbox("Colorir por", ["cell_type", "sample", "condition", "donor", "leiden", "leiden_anotacao"],
                             format_func={"cell_type": "tipo celular", "sample": "amostra", "condition": "condição",
                                          "donor": "doador", "leiden": "grupo Leiden",
                                          "leiden_anotacao": "grupo fino (anotação)"}.get)
        before = "X_umap_sem_integracao" in a.obsm and st.toggle("Mostrar antes da integração")
        if before:
            a.obsm["X_umapsem"] = a.obsm["X_umap_sem_integracao"]
        fig = plots.umap_panels(a, [color], basis="umapsem" if before else "umap")
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
        st.caption("O mapa serve para ver vizinhanças. Distâncias entre ilhas e tamanhos de ilhas não são medidas.")

    with tabs[2]:
        fig = plots.composition_bars(res.composition, cfg.condition_key)
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
        st.dataframe(res.composition, use_container_width=True)

    with tabs[3]:
        if res.de.empty:
            st.info("Sem expressão diferencial: é preciso ter 2 condições e pelo menos 2 amostras por condição em cada tipo.")
        else:
            ct = st.selectbox("Tipo celular", sorted(res.de["cell_type"].unique()))
            fig = plots.volcano(res.de, ct)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)
            d = res.de[res.de["cell_type"] == ct]
            st.caption(f"{d['contraste'].iloc[0]} · desenho {d['desenho'].iloc[0]} · DESeq2 em pseudobulk (contagens somadas por amostra)")
            st.dataframe(d[["gene", "log2FoldChange", "padj", "baseMean"]].head(200).round(4),
                         use_container_width=True, hide_index=True)

    with tabs[4]:
        st.dataframe(res.qc, use_container_width=True, hide_index=True)

    with tabs[5]:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in (cfg.results_dir / st.session_state.get("input_name", "analise")).glob("*"):
                zf.write(p, f"resultados/{p.name}")
            for p in res.figures.values():
                if Path(p).exists():
                    zf.write(p, f"figuras/{Path(p).name}")
        st.download_button("Baixar tabelas e figuras (.zip)", buf.getvalue(), file_name="immunomap_resultados.zip",
                           mime="application/zip", use_container_width=True)
        st.caption(f"O objeto completo foi salvo em {cfg.data_dir}/{st.session_state.get('input_name', 'analise')}_resultado.h5ad")
