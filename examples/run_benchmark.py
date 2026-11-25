"""Reproduce three shipped benchmarks and print dimensionless proxy KPIs.

Run from the repository root:

    python examples/run_benchmark.py

The printed values are this model's proxies for the scenario files. They are
not clinical measurements and they are not a claim of predictive accuracy.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.circuit_model import run_simulation
from src.config_loader import load_benchmark
from src.presentation import fmt_sig
from src.readouts import DEFAULT_ANALYSIS_WINDOW

BENCH = ROOT / "configs" / "benchmarks"
FILES = (
    "01_healthy_control.yaml",
    "02_untreated_pd.yaml",
    "06_ldopa_normalized_reference.yaml",
)
KPI_KEYS = (
    "beta_mean_window",
    "motor_mean_window",
    "beta_power_change_pct",
    "gpi_snr_output_rate",
    "direct_indirect_balance_index",
    "dyskinesia_risk_proxy",
    "stn_burstiness_proxy",
    "overdrive_peak",
    "on_time_h",
    "off_time_h",
)


def main() -> None:
    window = DEFAULT_ANALYSIS_WINDOW
    print(
        "Dimensionless proxies. Analysis window "
        f"{window[0]:g}-{window[1]:g} model hours. Seed is the file seed."
    )
    rows = []
    for name in FILES:
        config, pert, params, raw = load_benchmark(BENCH / name)
        result = run_simulation(config, pert, params)
        rows.append((raw.get("name", name), result))
        print(f"\n{raw.get('name', name)}  baseline={config.baseline_preset.value}")
        for key in KPI_KEYS:
            print(f"  {key}: {fmt_sig(result.kpis.get(key))}")

    healthy = rows[0][1].kpis
    depleted = rows[1][1].kpis
    print(
        "\nIn these files, the depleted baseline has a higher window-mean beta proxy "
        f"({fmt_sig(depleted['beta_mean_window'])} vs {fmt_sig(healthy['beta_mean_window'])}) "
        "and a lower window-mean motor-readiness proxy "
        f"({fmt_sig(depleted['motor_mean_window'])} vs {fmt_sig(healthy['motor_mean_window'])}) "
        "than the healthy baseline. That is a directional property of this model, "
        "not a measure of predictive accuracy."
    )


if __name__ == "__main__":
    main()
