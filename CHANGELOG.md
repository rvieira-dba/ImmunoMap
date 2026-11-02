# Changelog

Formato [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

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
