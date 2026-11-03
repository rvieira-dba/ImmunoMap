"""The app starts clean. Streamlit prints its deprecations to the terminal, not to the page, so the test reads stderr."""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
APP = Path(__file__).resolve().parents[1] / "src" / "immunomap" / "app" / "streamlit_app.py"


def test_app_source_has_no_deprecated_streamlit_argument():
    # Streamlit announced the removal of use_container_width after 2025-12-31; width="stretch" replaces it.
    assert "use_container_width" not in APP.read_text(encoding="utf-8")


def test_app_starts_without_errors_or_deprecation_notices(monkeypatch, capfd):
    import scanpy as sc
    from streamlit.testing.v1 import AppTest

    monkeypatch.setattr(sc.pp, "neighbors", lambda *a, **k: None)  # keeps the numba warm-up thread out of the test
    at = AppTest.from_file(str(APP), default_timeout=60).run()
    assert not at.exception
    err = capfd.readouterr().err.lower()
    assert "will be removed" not in err and "deprecat" not in err
