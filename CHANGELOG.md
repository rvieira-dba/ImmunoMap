# Changelog

Formato [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

## [0.3.2] - 2026-11-03

### Alterado
- **Versões mínimas com motivo e prova.** Python ≥ 3.12 (o 3.10 saiu de suporte em 2026-10-01, PEP 619; no 3.11 o scanpy para na 1.11.5, que chama `Colormap.set_bad`, depreciado no matplotlib 3.11, e o anndata para na 0.12.19, que exige `pandas < 3`), `scanpy` ≥ 1.12, `pydeseq2` ≥ 0.5.4, `anndata` ≥ 0.11, `pandas` ≥ 2.3 e < 3.1 e `matplotlib` ≥ 3.9. O teto do pandas é provisório: com o pandas 3.1.0, o anndata 0.13.4 emite `Pandas4Warning` e 13 dos 28 testes falham; o teto sai quando o anndata corrigir. O piso antigo do `pydeseq2` (≥ 0.4) estava errado: a 0.4 não tem o argumento `design=`, a 0.5.2 e a 0.5.3 quebram no pandas 3 (`TypeError`) e, no pandas 2.3, a 0.5.2 emite `FutureWarning`. O teto `harmonypy < 0.2` continua, agora com a prova (a 0.2.0 exige PyTorch; a 2.x não tem wheel para Windows). [`docs/VERSOES.md`](docs/VERSOES.md) explica cada piso, mostra a última versão que cada Python consegue instalar e quem ainda usa as versões anteriores, com a classe da prova de cada afirmação (medido, teste, documento ou argumento).
- **Depreciação não passa.** O `pyproject.toml` transforma `DeprecationWarning`, `PendingDeprecationWarning` e `FutureWarning` em erro nos testes. Numa análise, um aviso de depreciação do DESeq2 aparece como a nota "aviso de depreciação (…)" e a análise continua. Antes, o processo que roda o DESeq2 devolvia só o texto do aviso, e o aviso voltava como `UserWarning`: a depreciação ficava escondida.

### Adicionado
- `versoes.json` em cada análise (ao lado de `parametros.json`, e dentro do `.zip` do app): Python, sistema e versão de cada biblioteca. A mesma entrada com as mesmas versões dá o mesmo resultado.
- `provas/`: 13 scripts que refazem as medições citadas na documentação (PCA exato e seu tempo, calibração sob hipótese nula, pareado × não pareado, impressão digital entre pilhas, origem dos avisos, versões por Python, vizinhos do ecossistema, pacotes do `harmonypy`, compilação do numba, tempo e memória) e, em `provas/logs/2026-11-03/`, as saídas originais (caminhos temporários trocados por `<repositorio>`, `<venv>`, `<sobreposicao>` e `<trabalho>`).
- Testes: 28 (eram 21). Novos: a trava de depreciação, a categoria do aviso, o `versoes.json` e o app iniciando sem erro nem aviso de depreciação (`streamlit.testing`).
- `requirements.txt`: o ambiente exato do teste no Windows 11 com Python 3.14.7 (`pip install -r requirements.txt` instala tudo, inclusive o próprio ImmunoMap). As faixas aceitas continuam no `pyproject.toml`. O `wrapt` fica na 2.5.2 estável: o `formulaic` 1.2.2 pede `wrapt>=1.17.0rc1`, e esse `rc1` deixa o pip escolher uma versão de desenvolvimento (2.6.0.dev1).

### Removido
- `.agents/` e `.claude/` saíram do repositório e foram para o `.gitignore`: eram atalhos para dentro da `.venv` local e quebravam em qualquer outro computador.

### Corrigido
- Documentação de desempenho, refeita com medições de ponta a ponta pela linha de comando (Kang 2018, 2 núcleos, arquivo local): 191,5 s com 2 processos e 265,6 s com `--jobs 1` (a rodada final da 0.3.2, com a mesma configuração, levou 184,3 s). No mesmo roteiro de medição etapa por etapa, a 0.3.0 levou 268 s e a 0.3.1, 206 s; cerca de 12 s dessa diferença são importação e leitura a frio na medição da 0.3.0. O "cerca de 3,5 minutos" da 0.3.1 era uma estimativa. O PCA com `covariance_eigh` medido de novo: 11,8 → 3,6 s (cerca de 3 vezes; a 0.3.1 anunciou "cerca de 4" numa medição sem registro).
- Memória: o processo principal chega a 2,3 GB no Kang 2018 e cada processo extra, a 359 MB (o README dizia 300 MB por processo).
- CHANGELOG: a seção da 0.2.1 voltou; na 0.3.0 ela tinha sido renomeada e misturada à da 0.3.0.
- Tags do repositório: a `v0.2.0` daqui apontava para o commit da 0.1.0, com a mensagem do Cell-Atlas. Cada versão do ImmunoMap, da 0.1.0 à 0.3.2, passa a ter a própria tag.

### Verificado
- PCA: no PBMC3k, `arpack` e `covariance_eigh` dão os mesmos rótulos nas 6 combinações testadas; o solver `randomized` muda de 46 a 154 rótulos de 2.698 (1,7% a 5,7%).
- Hipótese nula (8 controles do Kang 2018 divididos 4 × 4 ao acaso, 6 vezes): o teste célula a célula chama de 78 a 183 genes "significativos" em todas as divisões; o pseudobulk com DESeq2, 0 genes em 5 das 6.
- Em cinco pilhas de bibliotecas, os rótulos finais do PBMC3k saem idênticos nas 6 combinações testadas (30 e 10 componentes × 3 sementes). Quatro pilhas, de Python 3.11 a 3.14 (scanpy 1.11.5 e 1.12.5, pandas 2.3.3 e 3.0, numpy 2.4 e 2.5, scikit-learn 1.9.1), são idênticas em todos os campos. A quinta, em Python 3.10 (numpy 2.2.6, scipy 1.15.3, scikit-learn 1.7.2, medida em outra máquina), difere em ruído de ponto flutuante no PCA e muda a partição do mapa em 1 das 6 combinações (9 grupos × 8); a causa, biblioteca ou máquina, não foi isolada.
- Kang 2018 de ponta a ponta pela linha de comando com a 0.3.2 (instalações limpas, 2 núcleos): as 57.908 linhas da expressão diferencial são idênticas às da 0.3.0 (diferença 0,0 em log2FC e padj), os rótulos são os mesmos em 100% das células e as métricas de integração não mudam, no Python 3.13 (184,3 s) e no 3.14.6 (253,2 s na primeira execução, com o cache do numba vazio; 190,7 s na segunda).
- Os 28 testes passam, com depreciação como erro, em Python 3.12, 3.13 e 3.14. No 3.11 com matplotlib 3.11.2, 1 falha (`set_bad` no scanpy 1.11.5); com matplotlib 3.10.9, os 28 passam. Com o pandas 3.1.0 (Python 3.12), 13 falham (`Pandas4Warning` no anndata 0.13.4). Dependências para Windows com Python 3.14: 87 pacotes, todos com wheel.
- Compilação do numba: a primeira busca de vizinhos do Kang 2018 leva cerca de 106 s com a pasta de cache vazia, cerca de 30 s com o cache em disco e menos de 1 s numa segunda chamada no mesmo processo (Python 3.13 e 3.14.6). A primeira análise depois de instalar leva, por isso, cerca de 77 s a mais, e a "lentidão do 3.14" da primeira medição era o cache vazio, não o Python.
- Pareado × não pareado (monócitos CD14+ do Kang 2018, os mesmos 16 pseudobulks): 2.185 genes no desenho pareado (`~donor + condition`) e 2.277 no não pareado; 35 só no pareado, 127 só no não pareado e 2.150 nos dois (`provas/pareamento.py`).
- Windows 11, Python 3.14.7, 4 núcleos, pandas 3.0.6: o Kang 2018 roda de ponta a ponta pelo app. Controle de qualidade (24.673 → 24.562 células) e mistura de amostras (0,61 → 0,79) iguais aos do Linux; separação de tipos e expressão diferencial próximas, não idênticas (0,61 → 0,64 contra 0,60 → 0,63; 57.862 contra 57.908 linhas; 2.201 contra 2.185 genes nos monócitos CD14+). A causa não foi isolada: compare resultados sempre na mesma máquina.

### Problemas conhecidos
- Com o pandas 3.1, fora da faixa suportada, o `Pandas4Warning` do anndata aparece como "nota do DESeq2" e com a dica errada de atualizar o pydeseq2. A correção entra na próxima versão.

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

### Adicionado
- Logo (pino de mapa com três grupos de células) no app (cabeçalho e ícone da aba), no README e num banner 1280×640 para a prévia do repositório no GitHub (`docs/assets/immunomap-banner.png`).
- Volcano destaca sempre os genes de interferon ISG15, IFIT1 e MX1.
- Figuras de referência do Kang 2018 em `docs/assets/`.

### Verificado
- Kang 2018 de ponta a ponta pelos dois caminhos, conjunto conhecido e arquivo bruto com as colunas escolhidas: resultados idênticos (24.562 de 24.673 células, 16 amostras, cerca de 4 a 5 minutos com 2 núcleos). O Harmony leva a mistura de amostras de 0,61 a 0,79 sem perder a separação de tipos (0,60 → 0,63). Os genes de interferon sobem em todos os tipos, mais nos monócitos CD14+ (2.185 genes; ISG15 com log2FC 7,8). Concordância com os rótulos dos autores: 91,1%.

## [0.2.1] - 2026-11-02

### Corrigido
- Saída do terminal sem ruído: some o aviso do scanpy sobre centralizar matriz esparsa (esperado nessa etapa) e os logs com data e hora do harmonypy 0.1.x.
- Avisos do PyDESeq2 viram uma nota curta em português por tipo celular (por exemplo, quando a curva de dispersão não converge e o DESeq2 usa a média).

### Verificado
- Kang 2018 de ponta a ponta (24.673 células, 16 amostras, scanpy 1.12.5): cerca de 5 minutos com 2 núcleos. Harmony leva a mistura de amostras de 0,61 a 0,79 sem perder a separação de tipos (0,60 → 0,63); os genes de interferon (ISG15, IFIT1, MX1) sobem em todos os tipos, mais nos monócitos CD14+ (2.185 genes diferenciais). Concordância com os rótulos dos autores: 91,1%.

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

[0.3.2]: https://github.com/rubenslyra/ImmunoMap/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/rubenslyra/ImmunoMap/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/rubenslyra/ImmunoMap/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/rubenslyra/ImmunoMap/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/rubenslyra/ImmunoMap/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/rubenslyra/ImmunoMap/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/rubenslyra/ImmunoMap/releases/tag/v0.1.0
