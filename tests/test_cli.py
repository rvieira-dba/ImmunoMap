import pytest

from immunomap import __version__
from immunomap.cli.main import main
from immunomap.config import Config
from immunomap.datasets import REGISTRY


def test_version(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert __version__ in capsys.readouterr().out


def test_datasets_listed(capsys):
    assert main(["datasets"]) == 0
    out = capsys.readouterr().out
    for key in REGISTRY:
        assert key in out


def test_jobs_flag_reaches_config():
    from immunomap.cli.main import _config, build_parser

    args = build_parser().parse_args(["run", "--dataset", "kang2018", "--jobs", "1"])
    assert _config(args).n_jobs == 1
    assert _config(build_parser().parse_args(["run", "--dataset", "kang2018"])).n_jobs == 0


def test_unknown_config_key(tmp_path):
    p = tmp_path / "c.json"
    p.write_text('{"resolucao": 1}', encoding="utf-8")
    with pytest.raises(ValueError):
        Config.from_file(p)
