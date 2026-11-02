# ImmunoMap

**scRNA-seq com várias amostras: do arquivo bruto à resposta de cada tipo celular.**

O ImmunoMap é a continuação do [Cell-Atlas](https://github.com/rubenslyra/Cell-Atlas), que analisa uma amostra só. Aqui a pergunta muda: *o que muda nas células entre condições* (controle × estimulado, saudável × doente), separando o que é biologia do que é diferença técnica entre amostras.

```
arquivos (.h5ad · .h5 · .mtx · .csv) ──► QC por amostra ──► normalização + genes variáveis por lote ──► PCA
   └─ ou conjunto conhecido                                                                              │
                                                     Harmony (integração) ◄──────────────────────────────┘
                                                            │
                                   Leiden + UMAP ──► anotação por marcadores ──► pseudobulk ──► DESeq2
```

---

## Instalação

```powershell
cd D:\source\Academicos-Projetos\ImmunoMap
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
pytest                      # confere a instalação
# opcional:
pip install -e ".[scvi]"    # integração com scVI (instala PyTorch, ~2 GB)
```

Requer Python ≥ 3.10. As dependências foram conferidas para Windows com Python 3.14: todas têm pacote pronto (wheel), então nada precisa de compilador.

## Uso

### Pelo navegador (também no tablet)

```powershell
immunomap app          # abre em http://localhost:8501
immunomap app --lan    # aceita o tablet na mesma rede: http://<IP-do-PC>:8501
```

Para descobrir o IP do PC, rode `ipconfig` e use o "Endereço IPv4". Na primeira vez, o Windows pode pedir para liberar o Python no firewall: permita só em **redes privadas**.

O app tem quatro passos, numa coluna só:

1. **Dados.** Três formas de entrada:
   - **conjunto conhecido** (baixado e guardado em `dados/` no primeiro uso);
   - **arquivos enviados** (um por amostra; você informa amostra, condição e doador numa tabela);
   - **`.h5ad` já combinado** (você escolhe as colunas de amostra, condição e doador).
2. **Parâmetros:** integração, resolução, QC e semente.
3. **Rodar.**
4. **Resultados:**
   - Resumo e métricas de integração;
   - Mapa (antes e depois da integração);
   - Composição por amostra;
   - Expressão diferencial (volcano e tabela);
   - QC;
   - Baixar tudo em `.zip`.

### Pela linha de comando

```powershell
immunomap datasets                          # lista os conjuntos conhecidos
immunomap download kang2018                 # baixa e padroniza
immunomap run --dataset kang2018            # análise completa
immunomap run --samples amostras.csv        # seus arquivos, descritos numa planilha
immunomap run --input meus.h5ad --sample-col paciente --condition-col grupo --donor-col paciente
immunomap run --dataset pbmc_v2_v3 --integration none   # ver o efeito de lote sem correção
immunomap run --dataset pbmc3k --annotation-resolution 0.6 --seed 1   # ver o erro T CD8 → NK
python benchmarks/anotacao_pbmc3k.py                   # sementes × resolução da anotação
```

Formato da planilha ([`amostras.example.csv`](amostras.example.csv)):

| coluna | obrigatória | exemplo |
| --- | --- | --- |
| `sample_id` | sim | `P1_ctrl` |
| `condition` | sim | `ctrl` / `stim` |
| `path` | sim | arquivo `.h5ad`/`.h5`/`.csv` ou pasta `.mtx` (relativo à planilha) |
| `donor` | não | `P1`. Mesmo doador em duas condições ⇒ análise **pareada** |
| `format` | não | `auto`, `h5ad`, `10x-h5`, `10x-mtx`, `csv`, `tsv` |

O ImmunoMap exige **contagens brutas** (números inteiros). Dados já normalizados são recusados com uma mensagem explicando o motivo.

## Conjuntos conhecidos

| chave | conjunto | amostras | pergunta | fonte |
| --- | --- | --- | --- | --- |
| `pbmc3k` | PBMC3k, 10x Genomics (2016) | 1 | O pipeline reencontra os tipos imunes? | `scanpy.datasets.pbmc3k()` |
| `kang2018` | Kang et al. 2018, lúpus + IFN-β | 16 (8 pacientes × 2) | Como cada tipo responde ao interferon? | cópia h5ad do scverse (a mesma do `pertpy`) · GEO GSE96583 |
| `pbmc_v2_v3` | PBMC em Chromium v2 e v3 | 2 | Quanto é efeito técnico (lote)? | 10x Genomics (pbmc8k, pbmc_10k_v3) |

Os downloads precisam de internet. Se algum endereço mudar, a mensagem de erro mostra qual URL falhou; ela fica em `src/immunomap/datasets.py`.

## Saídas

Para cada análise `<nome>`:

| Arquivo | Conteúdo |
| --- | --- |
| `resultados/<nome>/qc_por_amostra.csv` | células antes/depois do QC, mediana de genes e % mitocondrial |
| `resultados/<nome>/metricas_integracao.json` | mistura de amostras e separação de tipos, antes e depois |
| `resultados/<nome>/composicao_por_amostra.csv` | % de cada tipo celular por amostra |
| `resultados/<nome>/pseudobulk_amostras.csv` | grupos amostra × tipo usados na estatística |
| `resultados/<nome>/expressao_diferencial.csv` | DESeq2 por tipo celular: log2FC, p, padj, desenho |
| `figuras/<nome>/*.png` | UMAP antes/depois, composição, volcano por tipo |
| `dados/<nome>_resultado.h5ad` | objeto completo para continuar no Python |

## Decisões de método

- **QC por amostra.** Cada amostra tem sua própria distribuição. O filtro de % mitocondrial só é aplicado se o arquivo tiver genes `MT-`.
- **Contagens brutas preservadas** em `layers["counts"]`: o pseudobulk soma contagens, não valores normalizados.
- **Genes variáveis por lote** (`batch_key=sample`): evita escolher genes que só diferenciam amostras.
- **Anotação em grupos finos.** O nome de cada tipo vem de um Leiden mais fino (`annotation_resolution` = 1,5) do que o do mapa (`resolution` = 0,6). Cada grupo recebe um nome só; se um grupo grosso junta T CD8 e NK, todos viram NK. Concordância com os rótulos de especialistas do PBMC3k, sementes 0 a 4 ([`benchmarks/anotacao_pbmc3k.py`](benchmarks/anotacao_pbmc3k.py)):

  | resolução da anotação | grupos | concordância média | pior semente |
  | --- | --- | --- | --- |
  | 0,6 | 7–9 | 87,9% | 82,6% |
  | 1,0 | 10 | 97,0% | 96,5% |
  | **1,5 (padrão)** | 12–13 | 97,4% | 96,9% |
  | 2,0 | 18–20 | 97,8% | 97,6% |

  A causa é a combinação de 30 componentes principais (padrão do ImmunoMap, pensado para conjuntos maiores) com a resolução 0,6: saem poucos grupos, e T CD8 e NK caem juntos. Com 10 componentes, como no Cell-Atlas, a resolução 0,6 fica estável (96,3–96,8%, `--n-pcs 10`). Anotando em 1,5, o resultado fica estável nos dois casos (96,7–97,8%).

  Lição: um número de acurácia só vale junto com os parâmetros e as sementes testadas. (Medido no Linux com scanpy 1.11; outras versões do igraph/scanpy mudam alguns décimos.)
- **Harmony** como padrão: rápido, roda na CPU e corrige o PCA sem mexer nas contagens. As métricas antes/depois mostram se ajudou.
- **Pseudobulk + DESeq2**: o teste compara **amostras**, não células. Assim, mil células de um paciente contam como um paciente, e não como mil réplicas. O desenho fica pareado (`~donor + condition`) quando cada doador aparece nas duas condições.

## Desenvolvimento

```powershell
pytest        # 16 testes: leitores, planilha, pseudobulk, métricas, pipeline ponta a ponta, CLI
```

## Referências

- Kang, H. M. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nat. Biotechnol.* 36, 89 (2018). doi:10.1038/nbt.4042
- Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* 16, 1289 (2019). doi:10.1038/s41592-019-0619-0
- Luecken, M. D. et al. Benchmarking atlas-level data integration in single-cell genomics. *Nat. Methods* 19, 41 (2022). doi:10.1038/s41592-021-01336-8
- Squair, J. W. et al. Confronting false discoveries in single-cell differential expression. *Nat. Commun.* 12, 5692 (2021). doi:10.1038/s41467-021-25960-2
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15, 550 (2014). doi:10.1186/s13059-014-0550-8
- Muzellec, B. et al. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39, btad547 (2023). doi:10.1093/bioinformatics/btad547
- Heumos, L. et al. Best practices for single-cell analysis across modalities. *Nat. Rev. Genet.* 24, 550 (2023). doi:10.1038/s41576-023-00586-w

Material de estudo comentado, com versões em português: [`docs/MATERIAL.md`](docs/MATERIAL.md).

---

MIT · Rubens Lyra
