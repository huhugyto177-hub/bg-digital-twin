"""Run examples/custom_candidate.yaml through the Python API.

Run from the repository root:

    python examples/custom_candidate.py

The YAML uses the illustrative study-form defaults (D1 agonist, max_effect 0.8,
onset 0 h, Tmax 1.5 h, half-life 2 h). It does not contain an EC50 or a Cmax.
The printed KPIs are dimensionless proxies.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.circuit_model import run_simulation
from src.config_loader import load_baseline_params, parse_perturbation_yaml
from src.presentation import fmt_sig
from src.readouts import analysis_window_for
from src.schemas import BaselinePreset, SimulationConfig
from src.uncertainty import run_uncertainty

YAML_PATH = Path(__file__).resolve().parent / "custom_candidate.yaml"
KPI_KEYS = (
    "beta_mean_window",
    "motor_mean_window",
    "beta_power_change_pct",
    "gpi_snr_output_rate",
    "direct_indirect_balance_index",
    "dyskinesia_risk_proxy",
    "stn_burstiness_proxy",
    "therapeutic_index_margin",
    "overdrive_peak",
)


def main() -> None:
    pert = parse_perturbation_yaml(YAML_PATH.read_text(encoding="utf-8"))
    cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=8.0,
        dt_h=0.05,
        noise_std=0.0,
        seed=42,
    )
    params = load_baseline_params(cfg.baseline_preset)
    result = run_simulation(cfg, pert, params)
    unc = run_uncertainty(cfg, pert, params, n_runs=8)
    start, end = analysis_window_for(cfg)

    print(f"candidate: {result.perturbation_name}")
    print(f"baseline: {result.baseline_preset}  seed: {result.seed}")
    print(f"analysis window: {start:g}-{end:g} model hours")
    print(f"notes: {pert.notes.strip()}")
    print("\nPoint estimate, then P10-P90 across 8 parameter draws:")
    for key in KPI_KEYS:
        print(
            f"  {key}: {fmt_sig(result.kpis.get(key))}  "
            f"[{fmt_sig(unc.p10.get(key))}, {fmt_sig(unc.p90.get(key))}]"
        )
    print("\nAssumptions recorded on the result:")
    for line in result.assumptions:
        print(f"  - {line}")


if __name__ == "__main__":
    main()
