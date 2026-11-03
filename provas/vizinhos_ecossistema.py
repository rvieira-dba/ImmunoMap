"""Quem fica preso a versões antigas por causa de OUTRO pacote? Para cada pacote do ecossistema de célula única, lê o
Requires-Dist da última versão publicada (metadado oficial do PyPI: https://pypi.org/pypi/<pacote>/json) e confere se a
faixa que ele aceita para scanpy, anndata, pandas, numpy, scikit-learn, pydeseq2 e harmonypy admite as versões que o
ImmunoMap 0.3.2 instala hoje (Python 3.13 limpo). Também mostra o Requires-Python.
Não diz quantas pessoas usam cada pacote nem por quê: mostra o MECANISMO (uma faixa que exclui a versão nova).
Uso: python -I provas/vizinhos_ecossistema.py   (precisa de rede e do pacote packaging)"""
import datetime as dt
import json
import urllib.error
import urllib.request

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version

# Versões que a instalação limpa do ImmunoMap 0.3.2 escolheu (Python 3.13, versoes.json da rodada final)
NOSSA_PILHA = {"scanpy": "1.12.5", "anndata": "0.13.4", "pandas": "3.0.5", "numpy": "2.5.3", "scikit-learn": "1.9.1",
               "pydeseq2": "0.5.4", "harmonypy": "0.1.0"}
ALVOS = list(NOSSA_PILHA)
PACOTES = ["squidpy", "scvi-tools", "cellrank", "scvelo", "decoupler", "pertpy", "scirpy", "muon", "mudata", "spatialdata",
           "celltypist", "scrublet", "bbknn", "scanorama", "liana", "omicverse", "scib", "scgen", "scarches", "rapids-singlecell",
           "infercnvpy", "palantir", "moscot", "ehrapy", "cell2location", "popv", "pydeseq2", "scanpy", "anndata", "sc-toolbox"]
ENV = {"extra": "", "python_version": "3.13", "python_full_version": "3.13.1", "sys_platform": "linux", "platform_system": "Linux",
       "os_name": "posix", "platform_machine": "x86_64", "implementation_name": "cpython", "platform_python_implementation": "CPython"}


def faixa(req_dist, alvo):
    out = []
    for r in req_dist or []:
        try:
            q = Requirement(r)
        except Exception:
            continue
        if q.name.lower().replace("_", "-") != alvo:
            continue
        if q.marker is not None:
            if "extra" in str(q.marker):
                continue  # dependência opcional: não conta
            try:
                if not q.marker.evaluate(ENV):
                    continue
            except Exception:
                pass
        out.append(str(q.specifier) or "(qualquer)")
    return out


hoje = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
print(f"consulta em {hoje} (UTC) · fonte: pypi.org/pypi/<pacote>/json (última versão de cada pacote, dependências sem 'extra')")
print("pilha do ImmunoMap 0.3.2 usada na conferência: " + ", ".join(f"{k} {v}" for k, v in NOSSA_PILHA.items()) + "\n")
print(f"{'pacote':18} {'versão':9} {'data':11} {'Requires-Python':16} faixas que o pacote aceita (✗ = exclui a versão do ImmunoMap)")
bloqueios = []
sem_pacote = []
for nome in PACOTES:
    try:
        j = json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{nome}/json", timeout=60))
    except urllib.error.HTTPError as e:
        sem_pacote.append(f"{nome} (HTTP {e.code})")
        continue
    info = j["info"]
    data = min((f["upload_time_iso_8601"] for f in j["urls"]), default="")[:10]
    partes = []
    for alvo in ALVOS:
        if alvo == nome:
            continue
        for sp in faixa(info.get("requires_dist"), alvo):
            ok = True if sp == "(qualquer)" else SpecifierSet(sp).contains(Version(NOSSA_PILHA[alvo]), prereleases=True)
            partes.append(f"{alvo} {sp} {'✓' if ok else '✗'}")
            if not ok:
                bloqueios.append((nome, info["version"], alvo, sp))
    print(f"{nome:18} {info['version']:9} {data:11} {(info.get('requires_python') or '—'):16} " + ("; ".join(partes) if partes else "(não declara faixa para esses pacotes)"))

print()
if sem_pacote:
    print("não encontrados no PyPI: " + ", ".join(sem_pacote))
print(f"pacotes cuja última versão EXCLUI algo da pilha do ImmunoMap: {len({b[0] for b in bloqueios})} de {len(PACOTES) - len(sem_pacote)} consultados")
for b in bloqueios:
    print(f"  {b[0]} {b[1]}: exige {b[2]} {b[3]}, e o ImmunoMap instala {b[2]} {NOSSA_PILHA[b[2]]}")
