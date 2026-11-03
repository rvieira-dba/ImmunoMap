"""Versions behind a result. Same input + same versions = same result; this file says which versions those were."""
from __future__ import annotations

import importlib.metadata as md
import os
import platform
from typing import Dict, Optional

from . import __version__

# The libraries whose version can change a number in the output (or make the run fail).
PACKAGES = ("scanpy", "anndata", "pandas", "numpy", "scipy", "scikit-learn", "igraph", "leidenalg", "numba",
            "pydeseq2", "harmonypy", "joblib", "matplotlib", "streamlit")


def _version(package: str) -> Optional[str]:
    try:
        return md.version(package)
    except md.PackageNotFoundError:
        return None


def environment_report() -> Dict[str, object]:
    """Python, system and library versions as a JSON-ready dict (written to ``versoes.json`` with each run)."""
    return {
        "immunomap": __version__,
        "python": platform.python_version(),
        "sistema": platform.system() + " " + platform.release(),
        "nucleos": os.cpu_count(),
        "bibliotecas": {p: _version(p) for p in PACKAGES},
    }
