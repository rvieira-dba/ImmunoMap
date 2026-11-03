"""
ImmunoMap CLI.

    immunomap datasets                     # lista os conjuntos conhecidos
    immunomap download kang2018            # baixa e padroniza em dados/
    immunomap run --dataset kang2018       # análise completa
    immunomap run --samples amostras.csv   # seus arquivos, descritos numa planilha
    immunomap run --input meus.h5ad --sample-col paciente --condition-col grupo
    immunomap app                          # interface no navegador (Streamlit)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

from .. import __version__
from ..config import INTEGRATIONS, Config


def _config(args) -> Config:
    cfg = Config.from_file(args.config) if getattr(args, "config", None) else Config()
    over = {k: v for k, v in {
        "integration": getattr(args, "integration", None),
        "resolution": getattr(args, "resolution", None),
        "annotation_resolution": getattr(args, "annotation_resolution", None),
        "random_state": getattr(args, "seed", None),
        "n_jobs": getattr(args, "jobs", None),
        "data_dir": getattr(args, "data_dir", None),
    }.items() if v is not None}
    return Config.from_mapping({**cfg.__dict__, **over}) if over else cfg


def cmd_datasets(args) -> int:
    from ..datasets import REGISTRY

    for d in REGISTRY.values():
        print(f"{d.key:12} {d.title}\n{'':12} {d.samples} — {d.question}\n{'':12} {d.citation}\n")
    return 0


def cmd_download(args) -> int:
    from ..datasets import load_dataset

    cfg = _config(args)
    cfg.ensure_dirs()
    a = load_dataset(args.name, cfg.data_dir, progress=print)
    print(f"✅ {args.name}: {a.n_obs} células × {a.n_vars} genes, {a.obs['sample'].nunique()} amostras")
    return 0


def cmd_run(args) -> int:
    from ..pipeline import run

    cfg = _config(args)
    cfg.ensure_dirs()
    if args.dataset:
        from ..datasets import load_dataset
        adata, name = load_dataset(args.dataset, cfg.data_dir, progress=print), args.dataset
    elif args.samples:
        from ..io import combine_samples, load_samplesheet
        adata, name = combine_samples(load_samplesheet(args.samples), progress=print), Path(args.samples).stem
    elif args.input:
        from ..io import read_counts
        import scanpy as sc
        src = sc.read_h5ad(args.input) if str(args.input).endswith(".h5ad") else None
        adata = read_counts(args.input)
        if src is not None:
            adata.obs = src.obs.reindex(adata.obs_names)
        for col, key in ((args.sample_col, "sample"), (args.condition_col, "condition"), (args.donor_col, "donor")):
            if col:
                adata.obs[key] = adata.obs[col].astype(str).astype("category")
        name = Path(args.input).stem
    else:
        print("Informe --dataset, --samples ou --input.", file=sys.stderr)
        return 2
    res = run(adata, cfg, progress=print, compute_de=not args.no_de, name=name)
    print("\n📋 QC por amostra:\n" + res.qc.to_string(index=False))
    if res.metrics:
        print("\n📋 Integração (0–1, maior é melhor):")
        for k, v in res.metrics.items():
            print(f"  {k:15} mistura de amostras {v['mistura_lotes']:.2f} · separação de tipos {v['separacao_tipos']:.2f}")
    if not res.de.empty:
        sig = res.de[(res.de["padj"] < 0.05) & (res.de["log2FoldChange"].abs() > 1)]
        print(f"\n📋 Genes diferenciais (padj < 0,05, |log2FC| > 1) por tipo:\n{sig.groupby('cell_type').size().to_string()}")
    print(f"\n✅ Resultados em {cfg.results_dir / name} e {cfg.figures_dir / name}")
    return 0


def cmd_app(args) -> int:
    app = Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(args.port),
           "--server.maxUploadSize", "4096", "--browser.gatherUsageStats", "false",
           "--theme.primaryColor", "#0f6e6a", "--client.toolbarMode", "minimal"]
    if args.lan:
        cmd += ["--server.address", "0.0.0.0"]
    return subprocess.call(cmd)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="immunomap", description="scRNA-seq com várias amostras: integração, anotação e pseudobulk.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("datasets", help="lista os conjuntos públicos conhecidos").set_defaults(func=cmd_datasets)

    d = sub.add_parser("download", help="baixa um conjunto conhecido para dados/")
    d.add_argument("name"); d.add_argument("--data-dir", type=Path); d.add_argument("-c", "--config", type=Path)
    d.set_defaults(func=cmd_download)

    r = sub.add_parser("run", help="roda a análise completa")
    src = r.add_mutually_exclusive_group()
    src.add_argument("--dataset", help="nome de um conjunto conhecido (ver `immunomap datasets`)")
    src.add_argument("--samples", type=Path, help="planilha CSV: sample_id, condition, path[, donor, format]")
    src.add_argument("--input", type=Path, help="um arquivo já combinado (.h5ad, .h5, .csv)")
    r.add_argument("--sample-col"); r.add_argument("--condition-col"); r.add_argument("--donor-col")
    r.add_argument("--integration", choices=INTEGRATIONS)
    r.add_argument("--resolution", type=float, help="grupos do mapa (padrão 0.6)")
    r.add_argument("--annotation-resolution", type=float,
                   help="grupos finos usados só para nomear os tipos (padrão 1.5)")
    r.add_argument("--seed", type=int)
    r.add_argument("--jobs", type=int,
                   help="tipos celulares do DESeq2 em paralelo (padrão 0 = todos os núcleos; 1 = um de cada vez)")
    r.add_argument("--no-de", action="store_true", help="pula a expressão diferencial")
    r.add_argument("--data-dir", type=Path); r.add_argument("-c", "--config", type=Path)
    r.set_defaults(func=cmd_run)

    a = sub.add_parser("app", help="abre a interface no navegador")
    a.add_argument("--port", type=int, default=8501)
    a.add_argument("--lan", action="store_true", help="aceita conexões da rede local (abrir no tablet)")
    a.set_defaults(func=cmd_app)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"❌ {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
