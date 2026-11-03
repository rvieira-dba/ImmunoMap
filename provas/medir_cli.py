"""Roda `immunomap run --dataset kang2018` de verdade (pelo main() da CLI), com hora em cada linha e memória medida.
Uso (na raiz do repositório; precisa de psutil):  python provas/medir_cli.py [--jobs 1 ...]"""
import io
import sys
import threading
import time

T0 = time.perf_counter()
import psutil

peak = {"main": 0, "filho_max": 0, "soma": 0, "filhos": 0}


def amostrar():
    p = psutil.Process()
    while True:
        try:
            main = p.memory_info().rss
            kids = [c.memory_info().rss for c in p.children(recursive=True)]
        except Exception:  # noqa: BLE001
            main, kids = 0, []
        peak["main"] = max(peak["main"], main)
        if kids:
            peak["filho_max"] = max(peak["filho_max"], max(kids))
            peak["filhos"] = max(peak["filhos"], len(kids))
        peak["soma"] = max(peak["soma"], main + sum(kids))
        time.sleep(0.25)


threading.Thread(target=amostrar, daemon=True).start()


class Carimbo(io.TextIOBase):
    def __init__(self, out):
        self.out, self.buf = out, ""

    def write(self, s):
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.out.write(f"{time.perf_counter() - T0:7.1f}  {line}\n")
            self.out.flush()
        return len(s)

    def flush(self):
        self.out.flush()


sys.stdout = Carimbo(sys.stdout)
from immunomap.cli.main import main  # noqa: E402

rc = main(["run", "--dataset", "kang2018", *sys.argv[1:]])
total = time.perf_counter() - T0
mb = lambda b: f"{b / 1e6:.0f} MB"  # noqa: E731
print(f"TOTAL {total:.1f} s · rc {rc} · memória: processo principal {mb(peak['main'])}, "
      f"maior processo filho {mb(peak['filho_max'])} ({peak['filhos']} filhos), pico somado {mb(peak['soma'])}")
