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
from immunomap.pipeline import normalise_design, run  # noqa: E402

ASSETS = Path(__file__).resolve().parent / "assets"
st.set_page_config(page_title="ImmunoMap", page_icon=str(ASSETS / "logo.png"), layout="centered",
                   initial_sidebar_state="collapsed")
st.logo(str(ASSETS / "logo.png"), size="large")


@st.cache_resource(show_spinner=False)
def _warm_up():
    """Compile numba's neighbour search once per server, in the background (about 30 s on 2 cores).

    The first neighbour search of each process pays this compilation; done while the user picks the
    data, the analysis itself starts warm. Needs a large input: a small one takes another code path.
    """
    import threading
    import warnings

    # numba enters and leaves warnings.catch_warnings() while it compiles, and catch_warnings() is not
    # thread-safe: it can undo the filter that steps.normalise_and_reduce sets for this expected warning
    # (6 leaks in 249 tries). A process-wide filter set before the thread starts survives (0 in 237).
    warnings.filterwarnings("ignore", message="zero-centering a sparse array/matrix densifies it")

    def work():
        try:
            import anndata as ad
            import numpy as np
            import scanpy as sc

            w = ad.AnnData(np.zeros((20000, 1), dtype="float32"))
            w.obsm["X_pca"] = np.random.default_rng(0).normal(size=(20000, 30)).astype("float32")
            sc.pp.neighbors(w, n_neighbors=15, use_rep="X_pca", random_state=0)
        except Exception:  # noqa: BLE001 — warm-up is only an optimisation
            pass

    thread = threading.Thread(target=work, name="immunomap-aquecimento", daemon=True)
    thread.start()
    return thread


_warm_up()
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
    if st.button("Carregar conjunto", type="primary", width="stretch"):
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
        table = st.data_editor(default, hide_index=True, width="stretch", disabled=["arquivo"],
                               key=f"sheet_{len(files)}")
        if st.button("Ler arquivos", type="primary", width="stretch"):
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

            sample_col = st.selectbox("Coluna da amostra", cols,
                                      index=_guess(cols, ["sample", "amostra", "sample_id", "orig.ident", "replicate",
                                                          "patient", "paciente", "donor", "doador", "ind"]),
                                      help="A que identifica cada amostra. Se for o paciente, e cada paciente aparecer "
                                           "nas duas condições, o ImmunoMap separa por condição sozinho.")
            condition_col = st.selectbox("Coluna da condição", cols,
                                         index=_guess(cols, ["condition", "condicao", "label", "stim", "group", "grupo"], min(1, len(cols) - 1)))
            donor_opts = ["(nenhuma)"] + cols
            donor_col = st.selectbox("Coluna do doador (opcional)", donor_opts,
                                     index=_guess(donor_opts, ["donor", "doador", "replicate", "patient", "paciente", "ind"]))
            if st.button("Usar este arquivo", type="primary", width="stretch"):
                try:
                    a = read_counts(path)
                    a.obs = meta.reindex(a.obs_names).copy()
                    a.obs["sample"] = a.obs[sample_col].astype(str).astype("category")
                    a.obs["condition"] = a.obs[condition_col].astype(str).astype("category")
                    a.obs["donor"] = (a.obs[donor_col] if donor_col != "(nenhuma)" else a.obs["sample"]).astype(str).astype("category")
                    note = normalise_design(a.obs, Config())
                    if note:
                        st.info(note)
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
                                                              celulas=("condition", "size")), width="stretch")

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
run_area = st.container()  # always present, so the elements below keep their place between reruns
if run_area.button("Rodar análise", type="primary", disabled=adata is None, width="stretch"):
    cfg = Config(integration=integration, resolution=float(resolution),
                 annotation_resolution=float(annotation_resolution), n_pcs=int(n_pcs),
                 min_genes_per_cell=int(min_genes), max_genes_per_cell=int(max_genes),
                 max_mito_percent=float(max_mito), random_state=int(seed))
    with run_area.status("Analisando…", expanded=True) as s:
        try:
            res = run(adata.copy(), cfg, progress=s.write, compute_de=do_de,
                      name=st.session_state.get("input_name", "analise"))
            st.session_state["result"], st.session_state["cfg"] = res, cfg
            st.session_state["result_id"] = st.session_state.get("result_id", 0) + 1
            st.session_state["cache"] = {}
            s.update(label="Análise concluída", state="complete", expanded=False)
        except Exception as exc:  # noqa: BLE001
            s.update(label="A análise falhou", state="error")
            run_area.error(str(exc))
if adata is None:
    run_area.caption("Carregue os dados no passo 1 para liberar este botão.")

# ------------------------------------------------------------------ 4. resultados
res = st.session_state.get("result")


def _cached(key, make):
    """Build once per result (figures as PNG, the zip): clicks then only swap images."""
    cache = st.session_state.setdefault("cache", {})
    if key not in cache:
        cache[key] = make()
    return cache[key]


def _png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


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
                         width="stretch")
            st.caption("Ambos vão de 0 a 1. Boa integração aumenta a mistura de amostras sem derrubar a separação de tipos.")
        if not res.de.empty:
            sig = res.de[(res.de["padj"] < 0.05) & (res.de["log2FoldChange"].abs() > 1)]
            st.subheader("Genes diferenciais por tipo")
            st.dataframe(sig.groupby("cell_type").size().rename("genes (padj < 0,05; |log2FC| > 1)"), width="stretch")

    @st.fragment
    def _tab_mapa():
        color = st.selectbox("Colorir por", ["cell_type", "sample", "condition", "donor", "leiden", "leiden_anotacao"],
                             format_func={"cell_type": "tipo celular", "sample": "amostra", "condition": "condição",
                                          "donor": "doador", "leiden": "grupo Leiden",
                                          "leiden_anotacao": "grupo fino (anotação)"}.get)
        before = "X_umap_sem_integracao" in a.obsm and st.toggle("Mostrar antes da integração")

        def _umap():
            if before:
                a.obsm["X_umapsem"] = a.obsm["X_umap_sem_integracao"]
            return _png(plots.umap_panels(a, [color], basis="umapsem" if before else "umap"))

        st.image(_cached(("umap", color, before), _umap), width="stretch")
        st.caption("O mapa serve para ver vizinhanças. Distâncias entre ilhas e tamanhos de ilhas não são medidas.")

    with tabs[1]:
        _tab_mapa()

    with tabs[2]:
        st.image(_cached("composicao", lambda: _png(plots.composition_bars(res.composition, cfg.condition_key))),
                 width="stretch")
        st.dataframe(res.composition, width="stretch")

    @st.fragment
    def _tab_diferencial():
        if res.de.empty:
            st.info("Sem expressão diferencial: é preciso ter 2 condições e pelo menos 2 amostras por condição em cada tipo.")
            return
        ct = st.selectbox("Tipo celular", sorted(res.de["cell_type"].unique()))
        st.image(_cached(("volcano", ct), lambda: _png(plots.volcano(res.de, ct))), width="stretch")
        d = res.de[res.de["cell_type"] == ct]
        st.caption(f"{d['contraste'].iloc[0]} · desenho {d['desenho'].iloc[0]} · DESeq2 em pseudobulk (contagens somadas por amostra)")
        st.dataframe(d[["gene", "log2FoldChange", "padj", "baseMean"]].head(200).round(4),
                     width="stretch", hide_index=True)

    with tabs[3]:
        _tab_diferencial()

    with tabs[4]:
        st.dataframe(res.qc, width="stretch", hide_index=True)

    with tabs[5]:
        def _zip():
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for p in (cfg.results_dir / st.session_state.get("input_name", "analise")).glob("*"):
                    zf.write(p, f"resultados/{p.name}")
                for p in res.figures.values():
                    if Path(p).exists():
                        zf.write(p, f"figuras/{Path(p).name}")
            return buf.getvalue()

        st.download_button("Baixar tabelas e figuras (.zip)", _cached("zip", _zip), file_name="immunomap_resultados.zip",
                           mime="application/zip", width="stretch")
        st.caption(f"O objeto completo foi salvo em {cfg.data_dir}/{st.session_state.get('input_name', 'analise')}_resultado.h5ad")
