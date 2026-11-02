# Changelog

Formato [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

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
