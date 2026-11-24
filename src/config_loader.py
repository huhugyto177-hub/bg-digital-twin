"""Load baseline and benchmark YAML configs into typed models."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import yaml

from .schemas import (
    BaselinePreset,
    CircuitParams,
    MoleculePerturbation,
    SimulationConfig,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"


def load_yaml(path: Path | str) -> dict[str, Any]:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}


def load_baseline_params(preset: BaselinePreset | str) -> CircuitParams:
    preset = BaselinePreset(preset) if isinstance(preset, str) else preset
    if preset == BaselinePreset.PD_DEPLETED:
        path = CONFIGS / "disease_pd_depleted.yaml"
    else:
        path = CONFIGS / "baseline_healthy.yaml"
    data = load_yaml(path)
    return CircuitParams(**{k: v for k, v in data.items() if k in CircuitParams.model_fields})


def parse_perturbation(data: Optional[dict[str, Any]]) -> Optional[MoleculePerturbation]:
    if not data:
        return None
    return MoleculePerturbation(**data)


def parse_perturbation_yaml(text: str | dict[str, Any]) -> MoleculePerturbation:
    """Parse a MoleculePerturbation from YAML text or a benchmark-style wrapper.

    Accepts a YAML string, an already-loaded ``dict``, a bare perturbation
    mapping, or a document with a top-level ``perturbation:`` key (as in
    ``configs/benchmarks/*.yaml``).
    """
    if isinstance(text, dict):
        raw = text
    else:
        raw = yaml.safe_load(text) or {}
    if not isinstance(raw, dict):
        raise ValueError("YAML must be a mapping (MoleculePerturbation or benchmark wrapper).")
    data = raw.get("perturbation", raw)
    if not isinstance(data, dict) or "name" not in data or "pk" not in data:
        raise ValueError(
            "Expected MoleculePerturbation fields (name, pk, modulations) or a "
            "benchmark YAML with a top-level 'perturbation' key."
        )
    return MoleculePerturbation(**data)


def load_benchmark(path: Path | str) -> tuple[SimulationConfig, Optional[MoleculePerturbation], CircuitParams, dict]:
    raw = load_yaml(path)
    preset = BaselinePreset(raw.get("baseline_preset", "healthy"))
    sim = raw.get("simulation", {})
    config = SimulationConfig(
        baseline_preset=preset,
        duration_h=sim.get("duration_h", 8.0),
        dt_h=sim.get("dt_h", 0.05),
        noise_std=sim.get("noise_std", 0.0),
        seed=sim.get("seed", 42),
        n_sensitivity_runs=sim.get("n_sensitivity_runs", 30),
        analysis_window_start_h=sim.get("analysis_window_start_h"),
        analysis_window_end_h=sim.get("analysis_window_end_h"),
    )
    perturbation = parse_perturbation(raw.get("perturbation"))
    params = load_baseline_params(preset)
    return config, perturbation, params, raw
