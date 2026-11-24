"""Automated Client Decision Brief generator (Markdown + high-DPI PNG charts)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .circuit_model import run_simulation
from .config_loader import load_benchmark
from .evidence_registry import (
    cou_and_limitations_markdown,
    load_evidence_registry,
    references_markdown,
)
from .schemas import MoleculePerturbation, SimulationConfig
from .uncertainty import run_uncertainty

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"


def _plot_trajectories(
    result,
    out_path: Path,
    title: str,
) -> None:
    t = np.asarray(result.time_h)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), dpi=160)
    fig.suptitle(title, fontsize=12)

    ax = axes[0, 0]
    for key in ("d1", "d2", "stn", "gpi", "thal"):
        ax.plot(t, result.states[key], label=key)
    ax.set_title("Circuit states")
    ax.set_xlabel("Time (h)")
    ax.legend(fontsize=7, ncol=2)

    ax = axes[0, 1]
    ax.plot(t, result.readouts["beta_proxy"], color="C3", label="beta proxy")
    ax.plot(t, result.readouts["motor_readiness"], color="C2", label="motor readiness")
    ax.set_title("Primary readouts")
    ax.set_xlabel("Time (h)")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(t, result.readouts["direct_pathway"], label="direct")
    ax.plot(t, result.readouts["indirect_hyperdirect"], label="indirect/hyperdirect")
    ax.plot(t, result.readouts["pathway_balance"], label="balance", linestyle="--")
    ax.set_title("Pathway decomposition")
    ax.set_xlabel("Time (h)")
    ax.legend(fontsize=7)

    ax = axes[1, 1]
    ax.plot(t, result.readouts["overdrive_index"], color="C1")
    ax.axhline(1.15, color="k", linestyle=":", linewidth=1, label="overdrive threshold")
    ax.set_title("Over-drive / dyskinesia risk proxy")
    ax.set_xlabel("Time (h)")
    ax.legend(fontsize=7)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_uncertainty(unc, out_path: Path, title: str) -> None:
    keys = ["beta_mean_window", "motor_mean_window", "overdrive_peak"]
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    x = np.arange(len(keys))
    means = [unc.mean[k] for k in keys]
    p10 = [unc.p10[k] for k in keys]
    p90 = [unc.p90[k] for k in keys]
    yerr = np.vstack([np.array(means) - np.array(p10), np.array(p90) - np.array(means)])
    ax.bar(x, means, yerr=yerr, capsize=4, color="steelblue", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(keys, rotation=15, ha="right")
    ax.set_title(title)
    ax.set_ylabel("KPI (mean ± P10–P90)")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def generate_decision_brief(
    config: SimulationConfig,
    perturbation: Optional[MoleculePerturbation],
    params,
    out_dir: Path,
    scenario_name: str = "scenario",
    n_sensitivity_runs: Optional[int] = None,
) -> Path:
    """Write Markdown brief + PNGs + config JSON to out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    result = run_simulation(config, perturbation, params)

    traj_png = out_dir / f"{scenario_name}_trajectories.png"
    unc_png = out_dir / f"{scenario_name}_uncertainty.png"
    _plot_trajectories(result, traj_png, f"{scenario_name}: trajectories & readouts")

    n_runs = n_sensitivity_runs or min(config.n_sensitivity_runs, 30)
    unc = run_uncertainty(config, perturbation, params, n_runs=n_runs)
    _plot_uncertainty(unc, unc_png, f"{scenario_name}: uncertainty (n={n_runs})")

    # Save exact configuration JSON
    cfg_json = {
        "scenario_name": scenario_name,
        "simulation": config.model_dump(mode="json"),
        "perturbation": perturbation.model_dump(mode="json") if perturbation else None,
        "circuit_params": params.model_dump(mode="json"),
        "kpis": result.kpis,
        "uncertainty": {
            "mean": unc.mean,
            "p10": unc.p10,
            "p50": unc.p50,
            "p90": unc.p90,
            "sensitivity_ranking": unc.sensitivity_ranking,
        },
    }
    (out_dir / f"{scenario_name}_config.json").write_text(
        json.dumps(cfg_json, indent=2), encoding="utf-8"
    )

    registry = load_evidence_registry()
    md_lines = [
        f"# Client Decision Brief — {scenario_name}",
        "",
        cou_and_limitations_markdown(registry),
        "",
        "## Scenario",
        f"- **Name:** {scenario_name}",
        f"- **Baseline preset:** {config.baseline_preset.value}",
        f"- **Perturbation:** {perturbation.name if perturbation else 'None'}",
        f"- **Seed:** {config.seed}",
        f"- **Duration:** {config.duration_h} h (dt={config.dt_h})",
        "",
        "## Key Performance Indicators",
        "",
        "| KPI | Point estimate | P10 | P50 | P90 |",
        "|---|---:|---:|---:|---:|",
    ]
    for k in sorted(result.kpis.keys()):
        md_lines.append(
            f"| {k} | {result.kpis[k]:.4f} | {unc.p10.get(k, float('nan')):.4f} | "
            f"{unc.p50.get(k, float('nan')):.4f} | {unc.p90.get(k, float('nan')):.4f} |"
        )

    md_lines.extend(
        [
            "",
            "## Sensitivity ranking (OAT, motor_mean_window)",
            "",
        ]
    )
    for name, score in unc.sensitivity_ranking:
        md_lines.append(f"- **{name}**: {score:.3f}")

    md_lines.extend(["", "## Assumptions", ""])
    for a in result.assumptions:
        md_lines.append(f"- {a}")
    if perturbation is not None:
        md_lines.append(
            "- Custom/novel molecule inputs (affinity, BBB, toxicity, clinical efficacy) "
            "are user-supplied or unknown; the twin only maps stated target modulation + "
            "PK exposure → circuit readouts."
        )

    md_lines.extend(["", references_markdown(registry)])

    md_lines.extend(
        [
            "## Figures",
            f"- ![{scenario_name} trajectories]({traj_png.name})",
            f"- ![{scenario_name} uncertainty]({unc_png.name})",
            "",
            f"Configuration JSON: `{scenario_name}_config.json`",
            "",
        ]
    )

    md_path = out_dir / f"{scenario_name}_decision_brief.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    return md_path


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Generate BG Digital Twin decision brief.")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to benchmark/scenario YAML (relative to bg_digital_twin/ or absolute).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output directory (default: outputs/<scenario>).",
    )
    args = parser.parse_args(argv)

    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = ROOT / cfg_path
    if not cfg_path.exists():
        raise SystemExit(f"Config not found: {cfg_path}")

    config, perturbation, params, raw = load_benchmark(cfg_path)
    scenario = raw.get("name", cfg_path.stem)
    out_dir = Path(args.out) if args.out else OUTPUTS / scenario
    md = generate_decision_brief(config, perturbation, params, out_dir, scenario_name=scenario)
    print(f"Wrote decision brief: {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
