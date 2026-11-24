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
    ax.set_ylabel("(dimensionless proxy)")
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
    ax.set_ylabel("(dimensionless proxy)")
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


def preclinical_brief_markdown(study, registry: Optional[dict] = None) -> str:
    """One-page decision brief as Markdown. Neutral attribution only."""
    from .evidence_registry import evidence_caption
    from .presentation import fmt_sig

    reg = registry or load_evidence_registry()
    lines = [
        "# Preclinical Decision Brief",
        "",
        "Basal Ganglia Digital Twin v1",
        "",
        "## Study question and Context of Use",
        "",
        study.question,
        "",
        cou_and_limitations_markdown(reg),
        "",
        f"Scaling: {study.scaling_note}",
        "",
    ]
    if getattr(study, "illustration", False):
        lines.append(
            "This run used a mechanistic illustration preset, not a calibrated regimen. "
            "potency not specified; normalized exposure shape"
        )
        lines.append("")
    if getattr(study, "potency_note", ""):
        lines.append(study.potency_note)
        lines.append("")
    lines.extend(
        [
            "## Dose sweep and pathway comparison",
            "",
            "| Dose | Level | Beta change % | Motor window | Overdrive peak | Therapeutic index | In band |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in study.doses:
        k = row.kpis
        lines.append(
            "| {label} | {level} | {beta} | {motor} | {od} | {ti} | {band} |".format(
                label=row.label,
                level=fmt_sig(row.level),
                beta=fmt_sig(k.get("beta_power_change_pct")),
                motor=fmt_sig(k.get("motor_mean_window")),
                od=fmt_sig(k.get("overdrive_peak")),
                ti=fmt_sig(k.get("therapeutic_index_margin")),
                band="yes" if row.in_band else "no",
            )
        )
    lines.append("")
    if study.inverted_u:
        lines.append(
            "Inverted-U: a higher exposure has a lower therapeutic-index margin than an interior dose."
        )
        lines.append("")
    lines.extend(["## Uncertainty and sensitivity", ""])
    ref = study.doses[len(study.doses) // 2] if study.doses else None
    if ref is not None and ref.uncertainty is not None:
        unc = ref.uncertainty
        lines.append(f"P10–P90 at dose {ref.label} (dimensionless proxy unless noted):")
        lines.append("")
        for key in ("beta_mean_window", "motor_mean_window", "overdrive_peak", "therapeutic_index_margin"):
            lines.append(
                f"- {key}: mean {fmt_sig(unc.mean.get(key))} "
                f"(P10 {fmt_sig(unc.p10.get(key))}, P90 {fmt_sig(unc.p90.get(key))})"
            )
        lines.append("")
        lines.append("Sensitivity drivers (share of |Δ motor window|):")
        lines.append("")
        for name, score in (unc.sensitivity_ranking or [])[:5]:
            lines.append(f"- {name}: {fmt_sig(score)}")
    else:
        lines.append("Uncertainty was not run for this export.")
    lines.extend(
        [
            "",
            "## Suggested next wet-lab experiment",
            "",
            study.suggestion,
            "",
            "## Limitations, non-claims, and citations",
            "",
            "All readouts are dimensionless model proxies. This brief does not state mg/kg, "
            "stereotaxic coordinates, or a clinical dose.",
            "",
            evidence_caption(
                reg,
                parameter_keys=["beta_proxy", "motor_readiness", "off_time", "d2_indirect", "d1_direct"],
                evidence_ids=list(getattr(study, "evidence_ids", []) or []),
                limit=8,
            ),
            "",
        ]
    )
    return "\n".join(lines)


def write_preclinical_brief_pdf(study, path: Path, registry: Optional[dict] = None) -> Path:
    """One A4 page via matplotlib. No extra PDF dependency."""
    from matplotlib.backends.backend_pdf import PdfPages

    from .presentation import apply_mpl_style, fmt_sig

    apply_mpl_style()
    reg = registry or load_evidence_registry()
    fig = plt.figure(figsize=(8.27, 11.69), dpi=120)
    fig.patch.set_facecolor("white")
    fig.text(0.06, 0.965, "Basal Ganglia Digital Twin v1", fontsize=12, fontweight="bold", color="#1B4F72")
    fig.text(0.06, 0.945, "Preclinical Decision Brief", fontsize=10, color="#1C2833")
    question = (study.question or "")[:220]
    fig.text(0.06, 0.925, "Study question: " + question, fontsize=7.5)
    fig.text(
        0.06,
        0.905,
        "Context of Use: mechanistic circuit-level translation only. Not affinity, BBB, ADMET, or clinical efficacy.",
        fontsize=7,
    )

    ax_table = fig.add_axes([0.06, 0.62, 0.88, 0.26])
    ax_table.axis("off")
    ax_table.set_title("Dose sweep (dimensionless proxies)", loc="left", fontsize=9)
    cell = [["Dose", "Level", "Beta %", "Motor", "Overdrive", "TI", "Band"]]
    for row in study.doses:
        k = row.kpis
        cell.append(
            [
                str(row.label),
                fmt_sig(row.level),
                fmt_sig(k.get("beta_power_change_pct")),
                fmt_sig(k.get("motor_mean_window")),
                fmt_sig(k.get("overdrive_peak")),
                fmt_sig(k.get("therapeutic_index_margin")),
                "yes" if row.in_band else "no",
            ]
        )
    table = ax_table.table(cellText=cell, loc="upper left", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.25)

    ax = fig.add_axes([0.08, 0.36, 0.84, 0.22])
    for row in study.doses:
        ax.plot(row.time_h, row.readouts["beta_proxy"], label=f"{row.label} beta")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Beta (dimensionless proxy)")
    ax.set_title("Beta proxy by dose", fontsize=9)
    ax.legend(fontsize=6, ncol=min(5, max(1, len(study.doses))), frameon=False)

    ref = study.doses[len(study.doses) // 2] if study.doses else None
    sens = "Uncertainty not run."
    if ref is not None and ref.uncertainty is not None and ref.uncertainty.sensitivity_ranking:
        top = ", ".join(f"{n} {fmt_sig(s)}" for n, s in ref.uncertainty.sensitivity_ranking[:3])
        sens = f"P10–P90 stored for {ref.label}. Top sensitivity drivers: {top}."
    flag = "Inverted-U flagged." if study.inverted_u else "No inverted-U on therapeutic index."
    fig.text(0.06, 0.33, "Uncertainty: " + sens, fontsize=7)
    fig.text(0.06, 0.312, flag, fontsize=7)
    suggestion = (study.suggestion or "")[:700]
    fig.text(0.06, 0.292, "Next experiment (suggestion, not a protocol): " + suggestion[:420], fontsize=6.5)
    claims = " | ".join(reg.get("non_claims", [])[:5])
    fig.text(0.06, 0.08, "Non-claims: " + claims[:420], fontsize=6)
    fig.text(0.06, 0.035, "Basal Ganglia Digital Twin v1 — dimensionless proxies only.", fontsize=7, color="#1B4F72")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)
    return path


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
