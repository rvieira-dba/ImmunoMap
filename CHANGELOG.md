# Changelog

Formato [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

## [0.3.1] - 2026-11-03

### Alterado
- **Análise mais rápida, com os mesmos resultados.** No Kang 2018 com 2 núcleos, a linha de comando caiu de cerca de 4,5 para 3,5 minutos, e o app, para cerca de 3. O Streamlit não era a causa: o tempo estava em três etapas.
  - **Expressão diferencial:** o DESeq2 roda um tipo celular por processo, em paralelo (147 → cerca de 95 s com 2 núcleos; com mais núcleos, o ganho cresce). Novo parâmetro `n_jobs` (`--jobs` na CLI): 0 usa todos os núcleos e 1 roda um tipo de cada vez. Se o paralelo falhar, a análise continua um tipo de cada vez.
  - **Compilação do numba:** a primeira busca de vizinhos de cada processo leva cerca de 30 s compilando código. O app faz isso em segundo plano assim que abre, enquanto você escolhe os dados.
  - **PCA:** passou a usar o método `covariance_eigh`. É exato como o `arpack` (diferença de cerca de 1e-11) e cerca de 4 vezes mais rápido (12,6 → 3,5 s).
- O mapa "antes da integração" não roda mais o Leiden, cujo resultado era descartado.
- App: as figuras e o `.zip` são gerados uma vez por análise e reaproveitados. Antes, cada clique redesenhava todas as abas; agora só a parte que mudou é atualizada (cerca de 0,4 s para trocar a cor do mapa).

### Corrigido
- App: trocar a cor do mapa ou o tipo celular do volcano fazia a página voltar para a aba "Resumo". Agora a aba continua aberta.
- App: o terminal mostrava "Please replace `use_container_width` with `width`" a cada clique (27 vezes numa sessão, no Streamlit 1.66). Os botões, imagens e tabelas passaram a usar `width="stretch"`.
- App: o aviso do scanpy sobre centralizar matriz esparsa ainda podia aparecer quando a análise começava durante a compilação em segundo plano. O `warnings.catch_warnings()` do Python não é seguro entre threads, e o numba entra e sai dele enquanto compila. O filtro agora vale para o processo todo, definido antes de a thread começar.
- Versões mínimas corrigidas: Streamlit ≥ 1.50 (o app usa `st.fragment`, `st.logo(size=…)` e `width="stretch"`) e scikit-learn ≥ 1.5 (com versões anteriores, o scanpy troca o `covariance_eigh` pelo `arpack` e avisa).

### Verificado
- Kang 2018 de ponta a ponta: as 57.908 linhas da expressão diferencial são idênticas às da 0.3.0 (log2FC e padj), os rótulos são os mesmos em 100% das células e as métricas de integração não mudam.
- PBMC3k: com `covariance_eigh` e com `arpack`, os rótulos saem idênticos, célula a célula, nas 40 combinações do benchmark (30 e 10 componentes × 5 sementes × 4 resoluções). A tabela de concordância do README não muda.
- 21 testes, entre eles um que compara o DESeq2 em paralelo com o de um tipo de cada vez.

## [0.3.0] - 2026-11-02

### Corrigido
- **Dados brutos com o paciente na coluna de amostra** (caso do arquivo original do Kang 2018, que só tem `replicate` e `label`): uma amostra com duas condições agora é separada em `<amostra>_<condição>`, e o paciente vira o doador, mantendo o desenho pareado. Antes, o pseudobulk somava controle e estimulado no mesmo grupo e a expressão diferencial saía errada sem aviso.
- Uma coluna que muda de célula para célula (como `nCount_RNA`) é recusada como amostra ou condição, com mensagem explicando o motivo.
- App: sugere `replicate`/`patient` como coluna da amostra quando não existe uma coluna `sample`.
- Saída do terminal sem ruído: some o aviso do scanpy sobre centralizar matriz esparsa (esperado nessa etapa) e os logs com data e hora do harmonypy 0.1.x. Avisos do PyDESeq2 viram uma nota curta em português por tipo celular.

### Adicionado
- Logo (pino de mapa com três grupos de células) no app (cabeçalho e ícone da aba), no README e num banner 1280×640 para a prévia do repositório no GitHub (`docs/assets/immunomap-banner.png`).
- Volcano destaca sempre os genes de interferon ISG15, IFIT1 e MX1.
- Figuras de referência do Kang 2018 em `docs/assets/`.

### Verificado
- Kang 2018 de ponta a ponta pelos dois caminhos, conjunto conhecido e arquivo bruto com as colunas escolhidas: resultados idênticos (24.562 de 24.673 células, 16 amostras, cerca de 4 a 5 minutos com 2 núcleos). O Harmony leva a mistura de amostras de 0,61 a 0,79 sem perder a separação de tipos (0,60 → 0,63). Os genes de interferon sobem em todos os tipos, mais nos monócitos CD14+ (2.185 genes; ISG15 com log2FC 7,8). Concordância com os rótulos dos autores: 91,1%.

## [0.2.0] - 2026-11-02

### Alterado
- Os tipos celulares são nomeados num agrupamento Leiden mais fino (`annotation_resolution`, padrão 1,5), separado da resolução do mapa (`resolution`, 0,6). Um grupo grosso podia juntar dois tipos parecidos, como T CD8 e NK, e cada grupo recebe um nome só. No PBMC3k (30 componentes principais, sementes 0 a 4), a concordância com os rótulos de especialistas passou de 82,6–96,7% para 96,9–97,8%.
- Novo parâmetro no app ("Resolução da anotação"), na CLI (`--annotation-resolution`) e no TOML. O mapa pode ser colorido pelo grupo fino.

### Adicionado
- `benchmarks/anotacao_pbmc3k.py`: reproduz a comparação de sementes × resolução (e `--n-pcs` para comparar com o Cell-Atlas).

### Corrigido
- App: integração padrão volta a ser Harmony, como na CLI e no README (estava "Sem integração").
- Licença declarada no formato SPDX (`license = "MIT"`), como pede o setuptools ≥ 77. Some o aviso `project.license as a TOML table is deprecated` na instalação.

## [0.1.1] - 2026-11-02

### Corrigido
- Instalação no Windows: `harmonypy` fixado em 0.1.x (NumPy puro). A série 2.x não publica wheel para Windows e tentava compilar C++ (Armadillo).
- Kang 2018 baixado direto do arquivo h5ad do scverse, sem o `pertpy` (dependência pesada que não resolve no Windows com Python 3.14). O extra `[kang]` deixa de existir.

## [0.1.0] - 2026-11-02

### Added

- Projeto derivado do Cell-Atlas v0.2.0, voltado a **várias amostras**.
- Entrada de dados: `.h5ad`, `.h5` da 10x, pasta/arquivo `.mtx` (inclusive `.zip`/`.tar.gz`) e `.csv`/`.tsv`; planilha de amostras (`sample_id, condition, path, donor, format`); validação de contagens brutas.
- Conjuntos conhecidos sob demanda: PBMC3k, Kang 2018 (lúpus, IFN-β) e PBMC em duas químicas 10x (v2 × v3).
- QC por amostra, genes variáveis por lote, integração (Harmony; scVI opcional), Leiden, anotação por marcadores.
- Métricas de integração antes/depois (mistura de amostras e separação de tipos).
- Pseudobulk por amostra × tipo celular e expressão diferencial com PyDESeq2 (desenho pareado `~donor + condition` quando possível).
- CLI `immunomap` (`datasets`, `download`, `run`, `app`) e app Streamlit mobile-first, acessível pelo tablet com `immunomap app --lan`.
