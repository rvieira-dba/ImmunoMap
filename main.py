#!/usr/bin/env python
"""Atalho na raiz: ``python main.py <comando>`` (sem instalar o pacote)."""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if SRC.is_dir() and str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from immunomap.cli.main import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
