# Versões: o que o ImmunoMap exige, por quê, e quem ainda usa as anteriores

**Regra do projeto: usar as versões mais novas que foram testadas, e tratar depreciação como erro.**
Um aviso de depreciação (`DeprecationWarning`, `PendingDeprecationWarning`, `FutureWarning`) derruba a suíte de testes (`filterwarnings` no `pyproject.toml`). Durante uma análise ele aparece como uma nota, "aviso de depreciação (…)", e a análise continua. Piso de versão só existe com motivo (uma função que usamos ou um erro conhecido); teto só existe com prova (`harmonypy`).

Cada afirmação traz a classe da prova:

| Marca | Significa |
| --- | --- |
| **M** | medido: o comando foi executado e a saída guardada (script em [`provas/`](../provas/LEIA-ME.md)) |
| **T** | coberto por teste automatizado (`pytest`) |
| **D** | documento ou metadado oficial (PyPI, PEP, SPEC, código-fonte da própria biblioteca), com o endereço |
| **A** | argumento ou hipótese: **não foi medido** |

## 1. O que o `pyproject.toml` exige

| Pacote | Exigência | Por que existe | Prova |
| --- | --- | --- | --- |
| Python | ≥ 3.12 | O 3.10 chegou ao fim em 2026-10-01. O 3.11 ainda recebe correções de segurança até outubro de 2027, mas o pip entrega nele o scanpy 1.11.5 e o anndata 0.12.19 (as últimas versões que aceitam 3.11). Com o matplotlib 3.11, a suíte **falha** no 3.11 (o scanpy 1.11.5 chama `Colormap.set_bad`, depreciado), e o anndata 0.12.19 exige `pandas < 3`, o que prende o 3.11 no pandas 2. A SPEC 0 já daria o 3.11 por encerrado em outubro de 2025 | **M** 3.11 com matplotlib 3.11.2: 27 testes passam e 1 falha com `PendingDeprecationWarning: The set_bad function will be deprecated`; com matplotlib 3.10.9: 28 passam (isola a causa). **D** [PEP 619](https://peps.python.org/pep-0619/) (3.10 em fim de vida em 2026-10-01); [PEP 664](https://peps.python.org/pep-0664/) (3.11.0 em 2022-10-24, suporte até cerca de outubro de 2027); `Requires-Dist: pandas<3` no METADATA do anndata 0.12.19 |
| `pydeseq2` | ≥ 0.5.4 | A 0.5.2 e a 0.5.3 **quebram no pandas 3.0.5** (`TypeError` em `pandas/core/internals/blocks.py`) e, no pandas 2.3.3, a 0.5.2 emite `FutureWarning: Setting an item of incompatible dtype`. Antes da 0.5.0 não existe o argumento `design=` | **M** `provas/origem_avisos.py`: no pandas 3.0.5, a 0.5.4 roda sem aviso de depreciação e a 0.5.2 e a 0.5.3 dão `TypeError` (`origem_avisos_container.txt`); no pandas 2.3.3, a 0.5.2 emite 3 `FutureWarning` e a 0.5.4 nenhum (`origem_avisos_pandas233.txt`). **T** `test_deseq_deprecation_is_reported_as_deprecation` |
| `scanpy` | ≥ 1.12 | A 1.11.5 (a última que roda no 3.11) chama `cmap.set_bad(...)` em `scanpy/plotting/_tools/scatterplots.py:184`, e o matplotlib 3.11 marcou essa função como depreciada. A 1.12 está na faixa que o `covariance_eigh` exige (≥ 1.11): PCA exato e cerca de 3 vezes mais rápido (11,8 → 3,6 s no Kang 2018) | **M** o teste do 3.11 acima; 0 células com rótulo diferente do `arpack` (`provas/solver_pca.py`); tempos em `provas/tempo_pca.py` |
| `scikit-learn` | ≥ 1.5 | Com versões anteriores o scanpy troca o `covariance_eigh` por `arpack` e avisa | **D** `scanpy/preprocessing/_pca/__init__.py` (1.12.5): `covariance_eigh` só entra na lista de solvers aceitos `if pkg_version("scikit-learn") >= Version("1.5")` |
| `harmonypy` | ≥ 0.1, **< 0.2** | A 0.1.x só usa NumPy, pandas, scikit-learn e SciPy. A 0.2.0 exige **PyTorch** (`import harmonypy` falha sem ele; o pacote do `torch` 2.15.0 tem 127 MB no Windows e 539 MB no Linux, fora as dependências dele). A 2.x é compilada e **não publica pacote pronto (wheel) para Windows**: o pip tentaria compilar C++ | **M** `ModuleNotFoundError: No module named 'torch'` ao importar a 0.2.0. **M** `pip download harmonypy --only-binary=:all: --platform win_amd64 --python-version 3.14` escolhe a 0.2.0, porque a 2.1.0 (2026-10-06) só tem wheels para Linux e macOS. **D** `Requires-Dist: torch` no METADATA da 0.2.0; tamanhos dos arquivos do `torch` no PyPI (2026-11-03) |
| `joblib` | ≥ 1.4 | `Parallel(return_as="generator_unordered")`, que entrega cada tipo celular assim que termina | **D** no código da 1.3.2 só há uma frase "futuras versões vão aceitar `generator_unordered`"; a 1.4.0 valida e implementa o valor. **T** `test_parallel_de_matches_one_at_a_time` |
| `streamlit` | ≥ 1.50 | `st.fragment`, `st.logo(size=)` e `width="stretch"`. O argumento `use_container_width` será removido (aviso do Streamlit: "after 2025-12-31") | **D** wheels do Streamlit inspecionadas na 0.3.1. **M** 27 avisos por sessão antes da correção. **T** `tests/test_app.py` |
| `pandas` | ≥ 2.3, **< 3.1** | Piso: o scanpy 1.12 e o anndata 0.13 pedem ≥ 2.3. Teto (provisório): com o pandas 3.1.0 (2026-10-22), o anndata 0.13.4 emite `Pandas4Warning` ao gravar um `.h5ad` (lê a opção `future.infer_string`) e ao criar uma visão de um `AnnData` (usa a opção `mode.chained_assignment`). O ImmunoMap não usa nenhuma das duas opções, e não há o que corrigir do nosso lado | **M** Python 3.12 com pandas 3.1.0: 13 dos 28 testes falham pela trava; com pandas 3.0.6: 28 passam. **D** `anndata/_io/specs/methods.py:1269` e `anndata/_core/anndata.py:1101` (0.13.4). Para tirar o teto, rode `pip install -e ".[dev]"`, depois `pip install --no-deps pandas==3.1.0` e `python -m pytest -q`; quando os 28 passarem, o teto sai |
| `anndata`, `matplotlib` | ≥ 0.11, ≥ 3.9 | Mesmo piso do `pydeseq2` 0.5.4 (o instalador já exigiria) | **D** METADATA do `pydeseq2` 0.5.4 |
| `leidenalg`, `igraph` | ≥ 0.10 | Herdados das versões anteriores do projeto | **A** não testados com versões menores |

O que sai de cada análise: `versoes.json` (Python, sistema e versão de cada biblioteca) ao lado de `parametros.json`. A mesma entrada com as mesmas versões dá o mesmo resultado; se as versões mudam, o arquivo mostra quais.

## 2. Quem fica nas versões anteriores

### 2.1 O que cada Python consegue instalar (**D**, PyPI, consulta de 2026-11-03)

Última versão de cada pacote que o `Requires-Python` do PyPI deixa instalar em cada Python (e a data do lançamento). Refazer: `python -I provas/versoes_por_python.py`.

| Pacote | Python 3.10 | Python 3.11 | Python 3.12, 3.13 e 3.14 |
| --- | --- | --- | --- |
| scanpy | 1.11.5 (2025-10-21) | 1.11.5 (2025-10-21) | 1.12.5 (2026-10-29) |
| anndata | 0.11.4 (2025-03-26) | 0.12.19 (2026-07-02) | 0.13.4 (2026-09-18) |
| pandas | 2.3.3 (2025-09-29) | 3.1.0 (2026-10-22) | 3.1.0 (2026-10-22) |
| numpy | 2.2.6 (2025-05-17) | 2.4.6 (2026-05-18) | 2.5.4 (2026-10-10) |
| scipy | 1.15.3 (2025-05-08) | 1.17.1 (2026-02-23) | 1.18.1 (2026-08-21) |
| scikit-learn | 1.7.2 (2025-09-09) | 1.9.1 (2026-09-10) | 1.9.1 (2026-09-10) |
| pydeseq2 | **0.5.2** (2025-06-26) | 0.5.4 (2026-01-23) | 0.5.4 (2026-01-23) |
| streamlit | 1.65.0 (2026-10-02) | 1.66.1 (2026-10-21) | 1.66.1 (2026-10-21) |
| matplotlib | 3.10.9 (2026-04-24) | 3.11.2 (2026-09-11) | 3.11.2 (2026-09-11) |

Quem está no Python 3.10 não escolhe ficar no `pydeseq2` 0.5.2: é a última versão que o instalador entrega, e é justamente a que acusa depreciação no pandas 2.3 e quebra no pandas 3.

No Python 3.11 acontece algo parecido com o pandas: existe pandas 3.1.0 para ele, mas o último anndata que roda no 3.11 (0.12.19) declara `pandas<3`, então uma instalação limpa resolve para **pandas 2.3.3** (**M**, `freeze_py311.txt`; **D** METADATA do anndata 0.12.19).

### 2.2 Quanto do uso vem de cada Python (**M**, pypistats.org, 27 de outubro a 2 de novembro de 2026)

| Python | scanpy | pydeseq2 |
| --- | --- | --- |
| 3.10 | 8,0% | 4,5% |
| 3.11 | 18,3% | 17,6% |
| 3.12 | 35,4% | 37,4% |
| 3.13 | 35,1% | 39,1% |
| 3.14 | 3,1% | 1,3% |

Contagem de **downloads, não de pessoas**: as linhas da API (`category`, `date`, `downloads`) não trazem rótulo de integração contínua nem de espelho, então não separam teste de análise. Os valores foram lidos da API do pypistats (7 dias, 228.888 downloads do scanpy e 153.250 do pydeseq2) e estão em `pypistats.txt` (o JSON bruto não foi guardado); as categorias raras (3.9, 3.15, nulo) ficaram de fora e a soma pode variar cerca de 0,1 ponto. Em números: 26,3% dos downloads do scanpy vêm de Python 3.10 ou 3.11, que **não instalam o scanpy 1.12**.

### 2.3 Por que ficam (e o que a prova mostra de cada motivo)

1. **Não têm escolha** (**D**, tabela 2.1). O instalador respeita o `Requires-Python`; quem está num Python antigo recebe as bibliotecas antigas. Este é o único motivo que a prova mostra por inteiro.
2. **Outro pacote limita a versão** (**M**, `provas/vizinhos_ecossistema.py`, metadado do PyPI de 2026-11-03). De 30 pacotes do ecossistema de célula única, só 2 têm, na última versão, uma faixa que exclui o que o ImmunoMap instala: o `omicverse` 2.3.2 exige `pandas < 3.0.0` (e o ImmunoMap instala o pandas 3.0.5), e o `popv` 0.6.1 exige `scikit-learn < 1.8.0` e `harmonypy < 0.1.0`. Quem instala um desses fica em pandas 2 (ou numa `harmonypy` antiga) mesmo que queira atualizar. O mesmo vale dentro da pilha: o anndata 0.12.19 exige `pandas < 3`, e o ImmunoMap agora faz o mesmo com o pandas 3.1 (seção 1). O teto prova o **efeito**; o motivo que os autores tiveram não foi lido (**A**).
3. **A maioria avançou, e deixou o Python antigo para trás** (**M**, mesma consulta). 14 dos 30 pacotes (contando o próprio scanpy e o anndata) exigem Python ≥ 3.12 na última versão, entre eles squidpy, scvi-tools, cellrank, scirpy, muon, liana e o `pertpy`, que exige ≥ 3.13. 13 ainda declaram aceitar o Python 3.10 ou anterior; 10 deles têm a última versão de 2025 ou antes (a mais antiga é o scrublet, de 2020) e 3 foram publicados em 2026 (scvelo, palantir, omicverse). Esse é o "legado" que se consegue medir: a data da última versão. O `Requires-Python` declarado não prova que o pacote funcione nas bibliotecas novas (**A**).
4. **Evitar mudança que quebra** (**D** + **M**). O pandas 3.0 trocou comportamentos padrão (tipo texto `str` no lugar de `object`; Copy-on-Write) e as notas da versão pedem passar antes pela 2.3 "and ensure your code is working without warnings" ([notas do pandas 3.0](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)). O risco é real aqui: o `pydeseq2` 0.5.2 e o 0.5.3 falham no pandas 3.0.5. Quem prefere não correr esse risco fixa `pandas<3`, como o `omicverse` fez (item 2); a motivação de cada pessoa não foi medida (**A**).
5. **Reproduzir um resultado já publicado** (**M**). Em cinco pilhas (`impressao_digital_py3*.txt`), os rótulos finais do PBMC3k saem idênticos nas 6 combinações. Quatro pilhas, de Python 3.11 a 3.14, são idênticas em tudo, mesmo com scanpy 1.11.5 e 1.12.5 e pandas 2.3.3 e 3.0. A mais antiga (Python 3.10, numpy 2.2.6, scikit-learn 1.7.2, medida em outra máquina) muda a partição do mapa em 1 das 6 combinações (9 grupos × 8; 30 componentes, semente 0). A causa, biblioteca ou máquina, não foi isolada. O que a prova mostra: mudar de pilha pode mudar um número do mapa; quem precisa reproduzir fixa as versões (**A** quanto à motivação).
6. **Bibliotecas que testam o próprio piso** (**D**). A [SPEC 0](https://scientific-python.org/specs/spec-0000/) recomenda dar suporte a cada versão do Python por 3 anos e a cada versão das bibliotecas-base (NumPy, SciPy, pandas…) por 2, "para limitar o custo de manutenção". Os projetos que a seguem testam os pisos, e por isso as versões antigas continuam sendo exercitadas. Que o scanpy e o `pydeseq2` façam isso em seus testes (**A**): o GitHub não pôde ser consultado daqui.
7. **Servidores, contêineres e sistemas de suporte longo** (**A**). Plausível, não medido; o PyPI não identifica quem baixa e não há como nomear usuários.

### 2.4 O que não foi medido

- Downloads **por versão do pacote** (por exemplo, quantos usam o `pydeseq2` ≤ 0.5.3). O pypistats não separa por versão; a consulta abaixo, escrita para o conjunto público do PyPI no BigQuery, **não foi executada** (pede um projeto no Google Cloud):

  ```sql
  SELECT file.version AS versao, COUNT(*) AS downloads
  FROM `bigquery-public-data.pypi.file_downloads`
  WHERE project = 'pydeseq2' AND timestamp >= TIMESTAMP_SUB(CURRENT_TIMESTAMP(), INTERVAL 7 DAY)
  GROUP BY versao ORDER BY downloads DESC LIMIT 15
  ```

- Quem são as pessoas ou instituições, e em que situação (teste ou análise) usam as versões antigas. O que se mediu é o lado dos pacotes (seção 2.3, itens 2 e 3), não o das pessoas.
- Windows com Python 3.14 em execução (só você tem essa máquina); o que foi conferido é a existência de wheels para todas as dependências.

## 3. As pilhas testadas (Linux, 2 núcleos; instalações limpas, com o instalador escolhendo as versões)

| | Python 3.14.6 | Python 3.13.18 | Python 3.12.3 | Python 3.11.18 (fora do piso) |
| --- | --- | --- | --- | --- |
| scanpy / anndata | 1.12.5 / 0.13.4 | 1.12.5 / 0.13.4 | 1.12.5 / 0.13.4 | 1.11.5 / 0.12.19 |
| pandas / numpy / scipy | 3.0.6 / 2.5.4 / 1.18.1 | 3.0.5 / 2.5.3 / 1.18.1 | 3.0.6 / 2.5.4 / 1.18.1 | 2.3.3 / 2.4.6 / 1.17.1 |
| scikit-learn / pydeseq2 | 1.9.1 / 0.5.4 | 1.9.1 / 0.5.4 | 1.9.1 / 0.5.4 | 1.9.1 / 0.5.4 |
| matplotlib / streamlit | 3.11.2 / 1.66.1 | 3.11.2 / 1.66.1 | 3.11.2 / 1.66.1 | 3.11.2 ou 3.10.9 / 1.66.1 |
| Testes (28), depreciação = erro | **28 passam** | **28 passam** | **28 passam** | matplotlib 3.11.2: 27 passam, 1 falha · matplotlib 3.10.9: 28 passam |

Mais duas execuções que explicam os tetos: Python 3.12 com **pandas 3.1.0** (13 dos 28 testes falham por `Pandas4Warning` do anndata 0.13.4) e a pilha antiga do Python 3.10 (scanpy 1.11.5, anndata 0.11.4, pandas 2.3.3, numpy 2.2.6, scikit-learn 1.7.2, `pydeseq2` 0.5.2, em outra máquina), usada só para comparar rótulos (`impressao_digital_py310_vm.txt`) e fora do piso. Dependências para Windows com Python 3.14: `uv pip compile pyproject.toml --python-version 3.14 --python-platform x86_64-pc-windows-msvc --only-binary :all:` resolve 87 pacotes, todos com wheel (**M**, `resolucao_windows_py314.txt`). O que **não** foi executado: Windows de fato, e o Python 3.12 e o 3.14 só nesta máquina Linux.

## 4. Como refazer

```bash
python -m pytest -q                                          # 28 testes; depreciação = erro
python -I provas/versoes_por_python.py                       # tabela 2.1 (precisa de rede)
python -I provas/vizinhos_ecossistema.py                     # quem limita as versões, seção 2.3 (precisa de rede)
python -I provas/rodas_harmonypy.py                          # harmonypy: wheels e tamanho do torch (precisa de rede)
python -I provas/impressao_digital.py . dados/pbmc3k.h5ad    # compara duas pilhas: rode nas duas e compare os hashes
python -I provas/origem_avisos.py .                          # de onde sai cada aviso do DESeq2, com as versões
python -I provas/pareamento.py . dados/kang2018_resultado.h5ad    # pareado × não pareado (monócitos CD14+)
NUMBA_CACHE_DIR=$(mktemp -d) python provas/numba_real.py dados/kang2018_resultado.h5ad   # duas vezes, na mesma pasta: a frio × com cache
uv pip compile pyproject.toml --python-version 3.14 --python-platform x86_64-pc-windows-msvc --only-binary :all:
```

As demais provas (PCA, calibração sob hipótese nula, tempo e memória) estão em [`provas/LEIA-ME.md`](../provas/LEIA-ME.md), com os resultados de referência. As saídas originais de 2026-11-03 estão em `provas/logs/2026-11-03/`.
