"""Por que `harmonypy < 0.2`? Mostra, para cada versão do harmonypy, se existe pacote pronto (wheel) para Windows,
se ela exige PyTorch, e quanto pesa o `torch` no Windows e no Linux. Fonte: https://pypi.org/pypi/<pacote>/<versão>/json
Uso: python -I provas/rodas_harmonypy.py   (precisa de rede)"""
import datetime as dt
import json
import urllib.request


def pypi(caminho):
    return json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{caminho}/json", timeout=60))


print(f"consulta em {dt.datetime.now(dt.timezone.utc):%Y-%m-%d} (UTC) · fonte: pypi.org\n")
print(f"{'versão':8} {'data':11} {'wheels p/ Windows':18} {'tipos de arquivo':28} exige torch?")
for v in ("0.1.0", "0.2.0", "2.0.0", "2.0.2", "2.1.0"):
    j = pypi(f"harmonypy/{v}")
    arqs = j["urls"]
    win = [a["filename"] for a in arqs if "win" in a["filename"] and a["filename"].endswith(".whl")]
    puros = [a["filename"] for a in arqs if a["filename"].endswith("none-any.whl")]
    tipos = ", ".join(sorted({("wheel pura" if a["filename"] in puros else "wheel compilada" if a["filename"].endswith(".whl") else "código-fonte") for a in arqs}))
    torch = any(r.lower().startswith("torch") and "extra ==" not in r for r in (j["info"].get("requires_dist") or []))
    data = min(a["upload_time_iso_8601"] for a in arqs)[:10]
    print(f"{v:8} {data:11} {('sim' if (win or puros) else 'NÃO'):18} {tipos:28} {'SIM' if torch else 'não'}")

j = pypi("torch")
versao = j["info"]["version"]
print(f"\ntorch {versao}: tamanho do pacote por plataforma (Python 3.13)")
for a in j["urls"]:
    f = a["filename"]
    if "cp313-cp313-" in f and ("win_amd64" in f or "manylinux" in f and "x86_64" in f):
        print(f"  {f}  {a['size'] / 1e6:,.0f} MB")
