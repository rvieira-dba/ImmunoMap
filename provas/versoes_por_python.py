"""Para cada pacote: a última versão que ainda instala em cada Python (pelo Requires-Python que o PyPI guarda) e a data.
Fonte: https://pypi.org/pypi/<pacote>/json  (metadado oficial). Uso: python -I provas/versoes_por_python.py   (precisa de rede e do pacote packaging)"""
import datetime as dt
import json
import urllib.request
from packaging.version import Version, InvalidVersion
from packaging.specifiers import SpecifierSet

PACOTES = ["scanpy", "anndata", "pandas", "numpy", "scipy", "scikit-learn", "pydeseq2", "harmonypy", "streamlit", "matplotlib"]
PYTHONS = ["3.10", "3.11", "3.12", "3.13", "3.14"]
hoje = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
print(f"consulta em {hoje} (UTC) · fonte: pypi.org/pypi/<pacote>/json\n")
print(f"{'pacote':13} {'mais nova':10} {'(data)':11} " + " ".join(f"{'py' + p:>20}" for p in PYTHONS))
for nome in PACOTES:
    j = json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{nome}/json", timeout=60))
    rel = []
    for v, files in j["releases"].items():
        try:
            ver = Version(v)
        except InvalidVersion:
            continue
        if ver.is_prerelease or ver.is_devrelease or not files or all(f.get("yanked") for f in files):
            continue
        rp = next((f.get("requires_python") for f in files if f.get("requires_python")), None)
        data = min(f["upload_time_iso_8601"] for f in files)[:10]
        rel.append((ver, data, rp))
    rel.sort(reverse=True)
    linha = f"{nome:13} {str(rel[0][0]):10} {rel[0][1]:11} "
    for p in PYTHONS:
        ok = next(((v, d) for v, d, rp in rel if (not rp) or SpecifierSet(rp).contains(p + ".0")), None)
        linha += f"{(str(ok[0]) + ' ' + ok[1]) if ok else '—':>20} "
    print(linha)
