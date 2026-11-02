"""Central configuration (defaults < TOML/JSON file < CLI flags or app inputs)."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, List, Mapping

INTEGRATIONS = ("none", "harmony", "scvi")


@dataclass
class Config:
    data_dir: Path = field(default_factory=lambda: Path("dados"))
    results_dir: Path = field(default_factory=lambda: Path("resultados"))
    figures_dir: Path = field(default_factory=lambda: Path("figuras"))
    random_state: int = 0

    # Column names in adata.obs
    sample_key: str = "sample"
    condition_key: str = "condition"
    donor_key: str = "donor"

    # Quality control (applied per sample)
    min_genes_per_cell: int = 200
    max_genes_per_cell: int = 6000
    max_mito_percent: float = 15.0
    min_cells_per_gene: int = 3

    # Normalization / dimensionality
    target_sum: float = 1e4
    n_top_genes: int = 2000
    n_pcs: int = 30

    # Integration and clustering
    integration: str = "harmony"
    n_neighbors: int = 15
    resolution: float = 0.6

    # Annotation
    min_annotation_score: float = 0.5
    marker_genes: Dict[str, List[str]] = field(default_factory=lambda: {
        "T cells": ["CD3E", "CD3D", "IL7R"],
        "B cells": ["MS4A1", "CD79A"],
        "CD14+ Monocytes": ["CD14", "LYZ"],
        "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
        "NK cells": ["GNLY", "NKG7"],
        "Dendritic cells": ["FCER1A", "CST3"],
        "Platelets": ["PPBP", "PF4"],
    })

    # Pseudobulk differential expression
    min_cells_pseudobulk: int = 10
    reference_condition: str = ""  # empty = first condition alphabetically

    def __post_init__(self) -> None:
        self.data_dir, self.results_dir, self.figures_dir = map(Path, (self.data_dir, self.results_dir, self.figures_dir))
        if self.integration not in INTEGRATIONS:
            raise ValueError(f"integration must be one of {INTEGRATIONS}")
        if self.min_genes_per_cell >= self.max_genes_per_cell:
            raise ValueError("min_genes_per_cell must be lower than max_genes_per_cell")
        if self.n_pcs < 2 or self.n_neighbors < 2 or self.resolution <= 0:
            raise ValueError("n_pcs and n_neighbors must be >= 2 and resolution > 0")

    def ensure_dirs(self) -> None:
        for d in (self.data_dir, self.results_dir, self.figures_dir):
            d.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in d.items()}

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "Config":
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(values) - known)
        if unknown:
            raise ValueError(f"Unknown configuration keys: {', '.join(unknown)}")
        return cls(**dict(values))

    @classmethod
    def from_file(cls, path) -> "Config":
        path = Path(path)
        if path.suffix.lower() == ".json":
            values = json.loads(path.read_text(encoding="utf-8"))
        elif path.suffix.lower() == ".toml":
            try:
                import tomllib
            except ModuleNotFoundError:  # pragma: no cover
                import tomli as tomllib  # type: ignore
            with path.open("rb") as fh:
                values = tomllib.load(fh)
        else:
            raise ValueError("Config must be .toml or .json")
        values = values.get("immunomap", values)
        return cls.from_mapping(values)
