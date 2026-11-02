# Material de estudo do ImmunoMap

Trilha para entender **cada passo** da análise com várias amostras. Cada tema tem: a ideia em uma frase, onde aparece no código, a leitura original (inglês) e um caminho em português.

> Ordem sugerida: 1 → 2 → 3 → 4 → 5. Os temas 1 e 2 vêm do Cell-Atlas; aqui só o que muda.

---

## 1. Por que várias amostras mudam tudo

**Ideia.** Com uma amostra você descreve *quais* células existem. Com várias você pergunta *o que muda* entre grupos. Mas cada amostra também traz ruído técnico próprio (dia, operador, lote de reagente, chip).

- No código: `io/samplesheet.py` (cada célula recebe `sample`, `condition`, `donor`).
- Leitura: Heumos et al. 2023, *Best practices for single-cell analysis across modalities*, seções "Data integration" e "Differential gene expression". doi:10.1038/s41576-023-00586-w
- Livro aberto (o mesmo grupo): [sc-best-practices.org](https://www.sc-best-practices.org), capítulos *Integration* e *Differential gene expression analysis*.
- Em português: busque "efeito de lote RNA-seq" no [Google Acadêmico](https://scholar.google.com.br/scholar?hl=pt-BR&q=efeito+de+lote+RNA-seq) e "single cell RNA-seq aula" no [YouTube](https://www.youtube.com/results?search_query=single+cell+RNA-seq+aula+portugu%C3%AAs).

**Analogia de infra.** Amostra = servidor; condição = versão do deploy. Se cada versão rodou num datacenter diferente, você não sabe se a latência mudou pelo código ou pelo datacenter. O desenho do experimento (cada doador nas duas condições) é o que separa as duas coisas.

## 2. Efeito de lote e integração (Harmony)

**Ideia.** Integração alinha as amostras no espaço do PCA para que *o mesmo tipo celular* fique junto, sem apagar diferenças biológicas reais.

- No código: `analysis/steps.py → integrate()` (Harmony gera `X_harmony`); `integration_metrics()` mede antes/depois.
- Métricas usadas:
  - **mistura de amostras** (entropia entre vizinhos): perto de 1 = amostras bem misturadas;
  - **separação de tipos** (silhueta reescalada 0–1): alta = tipos continuam distintos.
  - O objetivo é subir a primeira **sem** derrubar a segunda.
- Leitura original: Korsunsky et al. 2019, Harmony. doi:10.1038/s41592-019-0619-0
- Comparação de métodos: Luecken et al. 2022, *Benchmarking atlas-level data integration*. doi:10.1038/s41592-021-01336-8
- Para ver na prática: `immunomap run --dataset pbmc_v2_v3 --integration none` e depois com `harmony`, comparando os dois UMAPs.

**Cuidado.** Integrar demais é um erro real: se as condições forem confundidas com o lote (todo controle num dia, todo estimulado em outro), a integração pode "corrigir" a própria biologia.

## 3. Composição: proporção de tipos por amostra

**Ideia.** Uma doença pode mudar *quantas* células de cada tipo existem, não só o que elas expressam.

- No código: `composition()` → `composicao_por_amostra.csv` e o gráfico de barras.
- Proporções somam 100%: se um tipo sobe, outros "descem" sem mudar de verdade. Métodos próprios tratam isso, por exemplo scCODA (Büttner et al. 2021, doi:10.1038/s41467-021-27150-6), que fica como passo futuro.

## 4. Pseudobulk e expressão diferencial (DESeq2)

**Ideia.** Para testar diferença entre condições, some as contagens de cada **tipo celular × amostra** e compare amostras, não células.

- No código: `pseudobulk()` soma `layers["counts"]`; `differential_expression()` roda PyDESeq2 por tipo celular com desenho `~donor + condition` (pareado) ou `~condition`.
- Por que não testar célula a célula: células do mesmo paciente não são independentes. Tratar mil células como mil réplicas gera falsos positivos em massa. Squair et al. 2021, doi:10.1038/s41467-021-25960-2
- DESeq2: Love, Huber & Anders 2014, doi:10.1186/s13059-014-0550-8. Versão Python: Muzellec et al. 2023, doi:10.1093/bioinformatics/btad547
- Vocabulário mínimo:
  - **log2FC**: log2 da razão entre condições; +1 = dobrou, −1 = caiu pela metade;
  - **p**: chance de ver uma diferença assim se nada mudasse;
  - **padj**: p corrigido para milhares de genes testados (Benjamini–Hochberg); use padj < 0,05;
  - **pareado**: cada doador é comparado com ele mesmo, o que remove a variação entre pessoas.
- Em português: o [vignette do DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) é em inglês, mas há aulas em português buscando "DESeq2 tutorial português" no [YouTube](https://www.youtube.com/results?search_query=DESeq2+tutorial+portugu%C3%AAs).

## 5. Os conjuntos de dados

| conjunto | o que é | o que ensina | origem |
| --- | --- | --- | --- |
| **PBMC3k** | ~2.700 células de sangue de um doador saudável (10x, 2016) | ponto de partida: o pipeline reencontra os tipos? | 10x Genomics; usado no tutorial do Scanpy |
| **Kang 2018** | PBMCs de 8 pacientes com lúpus, cada um em controle e estimulado com **interferon-β** por 6 h | resposta por tipo celular; desenho pareado | Kang et al. 2018, doi:10.1038/nbt.4042 · GEO GSE96583 |
| **PBMC v2 × v3** | PBMCs de doadores saudáveis em duas versões da química Chromium | efeito puramente técnico (lote) | 10x Genomics, datasets públicos pbmc8k e pbmc_10k_v3 |

**O que esperar no Kang 2018.** O interferon-β liga os *genes estimulados por interferon* (ISGs) em quase todos os tipos: **ISG15, IFI6, IFIT1, IFIT3, MX1, OAS1**. Monócitos costumam ter a resposta mais forte. Se o volcano de monócitos CD14 não mostrar esses genes no alto à direita, algo está errado no pipeline. É o "teste de fumaça" biológico do projeto.

**Por que o lúpus é relevante.** No lúpus eritematoso sistêmico, muitos pacientes têm uma "assinatura de interferon" elevada no sangue. Esse é um dos motivos do uso clínico de anifrolumabe (bloqueio do receptor de interferon tipo I). Referência de revisão: Crow, M. K. *Type I interferon in the pathogenesis of lupus*. J. Immunol. 192, 5459 (2014).

## 6. Exercícios guiados

1. Rode `pbmc_v2_v3` sem integração. As amostras se separam no UMAP? Anote as duas métricas.
2. Rode de novo com Harmony. Qual métrica mudou, e quanto?
3. Rode `kang2018`. Na aba Diferencial, procure ISG15 nos monócitos. Qual o log2FC?
4. Mude a resolução de 0,6 para 1,2. Os tipos se dividem? A DE fica mais ou menos estável?
5. Explique em 3 frases, para um recrutador, por que o teste é por amostra e não por célula.

---

Referências completas no [README](../README.md#referências).
