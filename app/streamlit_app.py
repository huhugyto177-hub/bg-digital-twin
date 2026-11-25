"""Basal Ganglia Digital Twin — Streamlit walkthrough.

Primary workflow is a preclinical study design. Benchmark and custom-candidate
modes remain available. Attribution in the UI is the product name only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import __version__
from src.catalog import (
    BASELINES,
    ILLUSTRATIVE_DOSE_LABELS,
    ILLUSTRATIVE_DOSE_LEVELS,
    ILLUSTRATIVE_DURATION_H,
    ILLUSTRATIVE_HALF_LIFE_H,
    ILLUSTRATIVE_ONSET_H,
    ILLUSTRATIVE_PK_LABEL,
    ILLUSTRATIVE_REFERENCE,
    ILLUSTRATIVE_TMAX_H,
    MECHANISM_ILLUSTRATIONS,
    SESSION_DEFAULTS,
    baseline_by_preset,
    custom_form_perturbation,
    illustration_by_id,
    study_candidate_name,
    study_question,
)
from src.circuit_model import simulate
from src.config_loader import load_baseline_params, load_benchmark, parse_perturbation_yaml
from src.evidence_registry import (
    cou_and_limitations_markdown,
    evidence_caption,
    evidence_table_rows,
    load_evidence_registry,
)
from src.guide import getting_started_markdown
from src.presentation import (
    PRIMARY_KPIS,
    apply_mpl_style,
    fmt_sig,
    frame_to_csv,
    sensitivity_figure,
    trajectories_frame,
    window_label,
)
from src.readouts import READOUT_DOCS, analysis_window_for, compute_readouts, summarize_kpis
from src.report_generator import preclinical_brief_markdown, write_preclinical_brief_pdf
from src.schemas import (
    BaselinePreset,
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    SimulationConfig,
    TargetModulation,
    TargetNode,
)
from src.study_design import (
    ILLUSTRATION_BANNER,
    POTENCY_NORMALIZED,
    PRESET_SPECS,
    default_scenario_duration_h,
    format_field_problems,
    parse_study_design_inputs,
    preset_reference_perturbation,
    run_study_from_base,
)
from src.uncertainty import run_uncertainty

BENCH_DIR = ROOT / "configs" / "benchmarks"
COU_BANNER = (
    "**Context of Use:** Mechanistic circuit-level translation only. "
    "This tool does **not** predict de novo binding affinity, BBB permeability, "
    "peripheral ADMET/toxicity, or clinical efficacy without experimental validation. "
    "It does **not** invent EC50, Cmax, or half-life."
)
CUSTOM_DISCLAIMER = (
    "**Custom / newly discovered molecules:** Affinity, BBB penetration, and toxicity "
    "are *user inputs or unknown*. The twin only translates your stated target "
    "modulation (node + mechanism + max_effect / EC50 if you supply one) and an exposure "
    "curve into circuit readouts. It does **not** invent binding from SMILES or structures."
)
APP_CSS = """
<style>
.stApp { background: #F6F4F0; }
div.block-container { padding-top: 1.1rem; max-width: 1180px; }
section[data-testid="stSidebar"] {
  background: #E7F0F0;
  border-right: 1px solid #D5E2E4;
}
section[data-testid="stSidebar"] > div { padding-top: 1.2rem; }
.bgdt-hero {
  background: linear-gradient(115deg, #14233A 0%, #1B3A4B 58%, #0F6E6B 100%);
  border-radius: 16px;
  padding: 1.35rem 1.6rem 1.2rem;
  margin: 0 0 0.9rem 0;
}
.bgdt-kicker {
  color: #9FDED6;
  letter-spacing: 0.16em;
  text-transform: uppercase;
  font-size: 0.72rem;
  margin: 0;
}
.bgdt-hero h1 {
  color: #F7F5F2;
  font-size: 1.85rem;
  line-height: 1.15;
  margin: 0.2rem 0 0.4rem 0;
  letter-spacing: -0.02em;
}
.bgdt-hero p { color: #E6F2F1; font-size: 1.02rem; margin: 0; }
div[data-testid="stMetric"] {
  background: #FFFFFF;
  border: 1px solid #D9E3E6;
  border-left: 4px solid #0F6E6B;
  border-radius: 12px;
  padding: 0.55rem 0.8rem 0.45rem;
  box-shadow: 0 1px 2px rgba(20, 35, 58, 0.05);
}
div[data-testid="stNumberInput"] input:disabled {
  color: #14233A !important;
  -webkit-text-fill-color: #14233A !important;
  opacity: 1;
}
div[data-testid="stTabs"] button[data-baseweb="tab"] {
  border-radius: 8px 8px 0 0;
  padding-top: 0.45rem;
  padding-bottom: 0.45rem;
}
div[data-testid="stTabs"] button[aria-selected="true"] {
  border-bottom: 2px solid #0F6E6B;
}
.bgdt-footer {
  color: #5C6B73;
  font-size: 0.86rem;
  padding-top: 0.2rem;
}
</style>
"""
HERO_HTML = """
<div class="bgdt-hero">
  <p class="bgdt-kicker">Preclinical CBGT circuit simulator</p>
  <h1>Basal Ganglia Digital Twin</h1>
  <p>Map a stated mechanism and an exposure curve to cortico–basal ganglia–thalamic readouts.</p>
</div>
"""
PK_HELP = (
    ILLUSTRATIVE_PK_LABEL
    + " This value is a shared scenario timebase, not a measured half-life, Tmax, or onset for a named drug."
)
CARD_KEYS = [
    ("motor_mean_window", "Motor readiness", "dimensionless proxy"),
    ("beta_power_change_pct", "Beta power change", "percent vs untreated reference"),
    ("gpi_snr_output_rate", "GPi/SNr output rate", "dimensionless rate proxy"),
    ("direct_indirect_balance_index", "Direct/indirect balance", "index, −1 to +1"),
    ("dyskinesia_risk_proxy", "Dyskinesia-risk proxy", "fraction of the analysis window"),
    ("stn_burstiness_proxy", "STN burstiness proxy", "dimensionless proxy"),
    ("off_time_h", "OFF-time proxy", "model hours in the analysis window"),
    ("therapeutic_index_margin", "Therapeutic index margin", "dimensionless proxy"),
]


def _list_benchmarks() -> dict[str, Path]:
    return {p.stem: p for p in sorted(BENCH_DIR.glob("*.yaml"))}


def _bench_label(stem: str) -> str:
    number, _, rest = stem.partition("_")
    if not rest:
        return stem
    return f"{number} · {rest.replace('_', ' ')}"


def _baseline_options() -> list[str]:
    return [item.preset for item in BASELINES]


def _baseline_title(preset: str) -> str:
    return baseline_by_preset(preset).title


def _finite_values(rows) -> list[float]:
    values: list[float] = []
    for row in rows:
        for key in ("y", "p10", "p90"):
            series = row.get(key)
            if series is None:
                continue
            for item in series:
                value = float(item)
                if value == value and value not in (float("inf"), float("-inf")):
                    values.append(value)
    return values


def _calm_ylim(ax, values: list[float]) -> None:
    """Keep a flat proxy looking flat, and print tick labels without a 1e-6 offset."""
    if not values:
        return
    lo, hi = min(values), max(values)
    pad = max((hi - lo) * 0.18, 0.02)
    ax.set_ylim(lo - pad, hi + pad)
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")


def _plot_series(study_or_rows, ylabel: str, title: str):
    apply_mpl_style()
    fig, ax = plt.subplots(figsize=(7.6, 3.5))
    rows = list(study_or_rows)
    for row in rows:
        ax.plot(row["t"], row["y"], label=row["label"])
        if row.get("p10") is not None:
            ax.fill_between(row["t"], row["p10"], row["p90"], alpha=0.16)
    _calm_ylim(ax, _finite_values(rows))
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    return fig


def _readout_panel(t, readouts, title: str):
    apply_mpl_style()
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.5), sharex=True, constrained_layout=True)
    specs = [
        (axes[0, 0], "beta_proxy", "Beta proxy", "dimensionless"),
        (axes[0, 1], "motor_readiness", "Motor readiness", "dimensionless"),
        (axes[1, 0], "gpi_output", "GPi/SNr output", "dimensionless rate"),
        (axes[1, 1], "direct_indirect_balance_index", "Direct/indirect balance", "index (−1 to +1)"),
    ]
    for ax, key, title_text, ylabel in specs:
        ax.plot(t, readouts[key], color="#0F6E6B")
        _calm_ylim(ax, _finite_values([{"y": readouts[key]}]))
        ax.set_title(title_text, fontsize=10)
        ax.set_ylabel(ylabel)
    axes[1, 0].set_xlabel("Time (h)")
    axes[1, 1].set_xlabel("Time (h)")
    fig.suptitle(title, fontsize=12)
    return fig


def _dose_response_figure(study):
    apply_mpl_style()
    levels = [d.level for d in study.doses]
    band = [d.level for d in study.doses if d.in_band]
    panels = [
        ("beta_power_change_pct", "Beta change", "percent vs untreated", "#0F6E6B"),
        ("motor_mean_window", "Motor readiness", "dimensionless", "#1B3A4B"),
        ("overdrive_peak", "Overdrive peak", "dimensionless", "#C47B4A"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.1), sharex=True)
    for ax, (key, title, unit, color) in zip(axes, panels):
        ax.plot(levels, [d.kpis[key] for d in study.doses], marker="o", color=color)
        if band:
            ax.axvspan(min(band), max(band), color="#0F6E6B", alpha=0.12)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(unit)
        ax.set_xlabel("Dose level (user scale)")
    fig.suptitle("Exposure sweep", fontsize=12)
    fig.tight_layout()
    return fig


def _metric_cards(kpis: dict, config, reference_label: str | None = None):
    window = window_label(config)
    for start in (0, 4):
        cols = st.columns(4)
        for col, (key, label, unit) in zip(cols, CARD_KEYS[start : start + 4]):
            col.metric(label, fmt_sig(kpis.get(key)))
            suffix = f" Reference dose {reference_label}." if reference_label else ""
            col.caption(f"{unit}. Window {window}.{suffix}")


def _kpi_table(study) -> pd.DataFrame:
    rows = []
    for dose in study.doses:
        unc = dose.uncertainty
        for key, unit in PRIMARY_KPIS:
            if key not in dose.kpis:
                continue
            rows.append(
                {
                    "dose": dose.label,
                    "metric": key,
                    "point": fmt_sig(dose.kpis[key]),
                    "mean": fmt_sig(unc.mean.get(key)) if unc else "n/a",
                    "p10": fmt_sig(unc.p10.get(key)) if unc else "n/a",
                    "p50": fmt_sig(unc.p50.get(key)) if unc else "n/a",
                    "p90": fmt_sig(unc.p90.get(key)) if unc else "n/a",
                    "unit": unit,
                    "window": window_label(study.config),
                    "definition": READOUT_DOCS.get(key, ""),
                }
            )
    return pd.DataFrame(rows)


def _config_from_session(session: dict, baseline: str, duration_h: float, n_sensitivity: int) -> SimulationConfig:
    return SimulationConfig(
        baseline_preset=BaselinePreset(baseline),
        duration_h=duration_h,
        dt_h=0.05,
        noise_std=session["noise_std"],
        seed=session["seed"],
        n_sensitivity_runs=n_sensitivity,
        analysis_window_start_h=session["analysis_window_start_h"],
        analysis_window_end_h=session["analysis_window_end_h"],
    )


def _window_ok(session: dict) -> bool:
    return session["analysis_window_end_h"] > session["analysis_window_start_h"]


def _show_study(study, evidence_ids: list[str] | None = None):
    if study.illustration:
        st.info(ILLUSTRATION_BANNER)
    if study.potency_note and study.potency_note not in (study.scaling_note or ""):
        st.caption(study.potency_note)
    st.write(study.scaling_note)
    if study.inverted_u:
        st.warning(
            "Inverted-U: therapeutic index is lower at a higher exposure than at an interior dose."
        )
    ref = next((d for d in study.doses if abs(d.level - 1.0) < 1e-9), study.doses[len(study.doses) // 2])
    _metric_cards(ref.kpis, study.config, ref.label)
    st.dataframe(_kpi_table(study), width="stretch", hide_index=True)
    st.pyplot(_dose_response_figure(study), clear_figure=True)
    st.pyplot(
        _readout_panel(ref.time_h, ref.readouts, f"Reference dose {ref.label}"),
        clear_figure=True,
    )
    series = []
    for dose in study.doses:
        item = {"label": dose.label, "t": dose.time_h, "y": dose.readouts["beta_proxy"]}
        if dose.uncertainty is not None and dose.uncertainty.trajectory_p10:
            item["p10"] = dose.uncertainty.trajectory_p10["beta_proxy"]
            item["p90"] = dose.uncertainty.trajectory_p90["beta_proxy"]
        series.append(item)
    st.pyplot(
        _plot_series(series, "Beta proxy (dimensionless)", "Beta proxy with P10–P90 where sampled"),
        clear_figure=True,
    )
    st.markdown(
        evidence_caption(
            load_evidence_registry(),
            parameter_keys=["beta_proxy", "motor_readiness", "d1_direct", "d2_indirect", "off_time"],
            evidence_ids=evidence_ids or study.evidence_ids,
        )
    )
    if ref.uncertainty is not None:
        st.pyplot(
            sensitivity_figure(ref.uncertainty.sensitivity_ranking, f"Sensitivity at {ref.label}"),
            clear_figure=True,
        )
    traj = trajectories_frame(
        [
            {"name": d.label, "t": d.time_h, "states": d.states, "readouts": d.readouts}
            for d in study.doses
        ]
    )
    kpi_csv = _kpi_table(study).to_csv(index=False)
    cfg = {
        "question": study.question,
        "candidate": study.candidate_name,
        "scaling_note": study.scaling_note,
        "illustration": study.illustration,
        "potency_note": study.potency_note,
        "simulation": study.config.model_dump(mode="json") if study.config else None,
        "doses": [
            {"label": d.label, "level": d.level, "perturbation": d.perturbation.model_dump(mode="json")}
            for d in study.doses
        ],
    }
    c1, c2, c3 = st.columns(3)
    c1.download_button("Trajectories CSV", frame_to_csv(traj), file_name="trajectories.csv", mime="text/csv")
    c2.download_button("KPI table CSV", kpi_csv, file_name="kpis.csv", mime="text/csv")
    c3.download_button(
        "Config JSON",
        json.dumps(cfg, indent=2),
        file_name="study_config.json",
        mime="application/json",
    )
    md = preclinical_brief_markdown(study)
    st.download_button("Decision brief (Markdown)", md, file_name="preclinical_decision_brief.md", mime="text/markdown")
    pdf_path = ROOT / "outputs" / "preclinical_decision_brief.pdf"
    write_preclinical_brief_pdf(study, pdf_path)
    st.download_button(
        "Decision brief (PDF)",
        pdf_path.read_bytes(),
        file_name="preclinical_decision_brief.pdf",
        mime="application/pdf",
    )
    with st.expander("Suggested next experiment and limitations"):
        st.write(study.suggestion)
        st.markdown(cou_and_limitations_markdown())


def _run_preset(preset_id: str, session: dict):
    spec, config, pert, params, raw = preset_reference_perturbation(preset_id)
    config = config.model_copy(
        update={
            "seed": session["seed"],
            "noise_std": session["noise_std"],
        }
    )
    study = run_study_from_base(
        config,
        pert,
        params,
        labels=["0.25×", "0.5×", "1×", "2×", "4×"],
        levels=[0.25, 0.5, 1.0, 2.0, 4.0],
        reference_level=1.0,
        use_hill=False,
        question=f"Illustration sweep for {spec['title']}: how does relative exposure change circuit proxies?",
        n_uncertainty=6,
        illustration=True,
        evidence_ids=list(raw.get("evidence_ids") or []),
        potency_note=str(raw.get("potency_note") or POTENCY_NORMALIZED),
    )
    st.session_state["study"] = study
    st.session_state["study_evidence"] = list(raw.get("evidence_ids") or [])


def _study_design_tab(session: dict):
    st.subheader("Preclinical study design")
    st.write(
        "The fields below are prefilled so Run works immediately. "
        + ILLUSTRATIVE_PK_LABEL
        + " Assay PK stays off until you turn it on and type values you already have."
    )
    st.markdown("**Mechanistic illustration presets, not calibrated.** One click runs a relative-exposure sweep.")
    cols = st.columns(len(PRESET_SPECS))
    for col, spec in zip(cols, PRESET_SPECS):
        col.caption(spec["blurb"])
        if col.button(spec["title"], key=f"preset_{spec['id']}", help="Runs this illustration with its file. Not measured PK."):
            with st.spinner("Running illustration sweep…"):
                _run_preset(spec["id"], session)
    st.divider()
    st.markdown("**Or use the prefilled candidate.** " + ILLUSTRATIVE_PK_LABEL)
    choice = st.selectbox(
        "Mechanism illustration",
        [item.id for item in MECHANISM_ILLUSTRATIONS],
        format_func=lambda item_id: illustration_by_id(item_id).title,
        help="Fills the target, the circuit mechanism, and an illustrative max_effect. Not a measured potency.",
    )
    spec = illustration_by_id(choice)
    st.caption(spec.blurb)
    name = st.text_input(
        "Candidate name",
        value=study_candidate_name(),
        help="A label only. The default is not a development code for a real molecule.",
    )
    question = st.text_input(
        "Study question",
        value=study_question(),
        help="Shown on the decision brief. The default asks about beta versus overdrive.",
    )
    nodes = [n.value for n in TargetNode]
    mechanisms = [m.value for m in Mechanism]
    c1, c2, c3 = st.columns(3)
    target = c1.selectbox(
        "Target node",
        nodes,
        index=nodes.index(spec.target_node),
        key=f"sd_target_{choice}",
        help="Circuit node that receives the gain. See Getting started for the mapping.",
    )
    mechanism = c2.selectbox(
        "Mechanism",
        mechanisms,
        index=mechanisms.index(spec.mechanism),
        key=f"sd_mech_{choice}",
        help="Model field that sets the sign of the gain. The catalog blurb states the pharmacological class.",
    )
    max_effect = c3.text_input(
        "max_effect at the reference dose",
        value=spec.max_effect,
        key=f"sd_me_{choice}",
        help="Unitless gain at the reference dose. Illustrative default, not a measured potency or EC50.",
    )
    ia = st.text_input(
        "Intrinsic activity (required for partial agonist)",
        value="0.5",
        help="Fraction of max_effect, greater than 0 and at most 1. Illustrative scenario fraction, not a measured potency. Used when the mechanism is partial agonist.",
    )
    st.caption(ILLUSTRATIVE_PK_LABEL + " Onset, Tmax, half-life, and duration are a scenario timebase.")
    p1, p2, p3 = st.columns(3)
    onset = p1.text_input("Onset (h)", value=ILLUSTRATIVE_ONSET_H, help=PK_HELP + " Unit: hours. Default 0.")
    tmax = p2.text_input("Tmax (h)", value=ILLUSTRATIVE_TMAX_H, help=PK_HELP + " Unit: hours after onset. Default 1.5.")
    half_life = p3.text_input(
        "Half-life (h)",
        value=ILLUSTRATIVE_HALF_LIFE_H,
        help=PK_HELP + " Unit: hours. Default 2.",
    )
    use_assay = st.checkbox(
        "Supply assay PK (EC50, reference Cmax, ka)",
        value=False,
        help="Off by default. The normalized path does not use EC50, Cmax, or ka. Turn this on only to type values you already measured. Nothing is invented.",
    )
    if use_assay:
        e1, e2, e3 = st.columns(3)
        ec50 = e1.text_input(
            "EC50 (nM)",
            value="",
            help="Nanomolar EC50 from an assay you supply. Required when assay PK is on. Left blank on purpose.",
        )
        cmax_ref = e2.text_input(
            "Reference Cmax (nM)",
            value="",
            help="Nanomolar Cmax you supply for the reference dose. Required with an EC50. Not invented.",
        )
        ka = e3.text_input(
            "ka (1/h, optional)",
            value="",
            help="Absorption rate in 1/h. Optional. Blank derives ka from Tmax. Not inferred from structure.",
        )
        reference = ILLUSTRATIVE_REFERENCE
    else:
        ec50 = ""
        cmax_ref = ""
        ka = ""
        reference = st.text_input(
            "Reference dose level",
            value=ILLUSTRATIVE_REFERENCE,
            help="User-scale level that maps onto max_effect. Default 1. Not a milligram dose.",
        )
        st.caption("Normalized path. EC50, reference Cmax, and ka stay off.")
    st.caption("Dose levels on the user scale. The first three are filled. Leave unused rows blank.")
    dose_labels, dose_levels = [], []
    for i in range(5):
        a, b = st.columns(2)
        lab = a.text_input(
            f"Dose {i + 1} label",
            value=ILLUSTRATIVE_DOSE_LABELS[i],
            key=f"dose_lab_{i}",
            help="Short label. Rows 4 and 5 are optional and start blank.",
        )
        lev = b.text_input(
            f"Dose {i + 1} level",
            value=ILLUSTRATIVE_DOSE_LEVELS[i],
            key=f"dose_lvl_{i}",
            help="Positive user-scale level. Defaults 0.5, 1, and 2. Not a milligram dose. Illustrative default, not measured PK.",
        )
        dose_labels.append(lab)
        dose_levels.append(lev)
    options = _baseline_options()
    baseline = st.selectbox(
        "Baseline",
        options,
        index=options.index(SESSION_DEFAULTS["baseline"]),
        format_func=_baseline_title,
        key="sd_base",
        help="Circuit state before the exposure is applied. Moderate PD depletion is the prefilled reference.",
    )
    st.caption(baseline_by_preset(baseline).summary)
    duration = st.text_input(
        "Duration (h)",
        value=ILLUSTRATIVE_DURATION_H,
        help="Scenario length in hours. Default 8. Not a clinical treatment interval. Illustrative default, not measured PK.",
    )
    default_duration = default_scenario_duration_h()
    if default_duration is None:
        st.caption("Enter Duration (h). No default scenario duration is set.")
    else:
        st.caption(
            f"Leave Duration (h) blank to use the default scenario duration of {default_duration:g} h. "
            + ILLUSTRATIVE_PK_LABEL
        )
    if st.button("Run entered study", type="primary", help="Runs the values currently in the form, including the illustrative defaults."):
        if not _window_ok(session):
            st.warning("Analysis window end (h) must be greater than the start. The study was not run.")
        else:
            parsed = parse_study_design_inputs(
                max_effect=max_effect,
                onset=onset,
                tmax=tmax,
                half_life=half_life,
                ka=ka,
                ec50=ec50,
                cmax_ref=cmax_ref,
                reference=reference,
                intrinsic_activity=ia,
                mechanism=mechanism,
                dose_labels=dose_labels,
                dose_levels=dose_levels,
                duration=duration,
            )
            if not parsed.ok:
                st.warning(format_field_problems(parsed.problems))
            elif not name.strip():
                st.warning("Enter a candidate name. The study was not run.")
            elif not (3 <= len(parsed.labels) <= 5):
                st.warning("Enter 3 to 5 doses (label and level). The study was not run.")
            else:
                if parsed.used_default_duration and parsed.duration_h is not None:
                    st.info(
                        f"Duration (h) was blank, so this run uses the default scenario duration of {parsed.duration_h:g} h."
                    )
                try:
                    pk_kwargs = {
                        "onset_h": parsed.onset_h,
                        "tmax_h": parsed.tmax_h,
                        "half_life_h": parsed.half_life_h,
                        "n_doses": 1,
                    }
                    if parsed.ka_per_h is not None:
                        pk_kwargs["ka_per_h"] = parsed.ka_per_h
                    mod = {
                        "target_node": target,
                        "mechanism": mechanism,
                        "max_effect": parsed.max_effect,
                        "hill_coeff": 1.0,
                    }
                    if parsed.use_hill:
                        mod["ec50_nM"] = parsed.ec50_nM
                        pk_kwargs["cmax_nM"] = parsed.cmax_nM
                    if parsed.intrinsic_activity is not None and mechanism == Mechanism.PARTIAL_AGONIST.value:
                        mod["intrinsic_activity"] = parsed.intrinsic_activity
                    pert = MoleculePerturbation(
                        name=name.strip(),
                        notes=POTENCY_NORMALIZED if not parsed.use_hill else "EC50 and Cmax are user-supplied inputs, not predictions.",
                        pk=PKProfile(**pk_kwargs),
                        modulations=[TargetModulation(**mod)],
                    )
                    cfg = _config_from_session(session, baseline, parsed.duration_h, 6)
                    params = load_baseline_params(cfg.baseline_preset)
                    study = run_study_from_base(
                        cfg,
                        pert,
                        params,
                        labels=parsed.labels,
                        levels=parsed.levels,
                        reference_level=parsed.reference_level,
                        use_hill=parsed.use_hill,
                        question=question.strip() or f"How does exposure change circuit proxies for {name.strip()}?",
                        n_uncertainty=6,
                        illustration=False,
                        evidence_ids=list(spec.evidence_ids),
                        potency_note="" if parsed.use_hill else POTENCY_NORMALIZED,
                    )
                    st.session_state["study"] = study
                    st.session_state["study_evidence"] = list(spec.evidence_ids)
                except Exception:
                    st.warning("The study was not run. Check the entered values and try again.")
    study = st.session_state.get("study")
    if study is not None:
        _show_study(study, st.session_state.get("study_evidence"))


def _run_benchmark(path: Path, use_file_duration: bool, duration_override: float):
    config, pert, params, raw = load_benchmark(path)
    if not use_file_duration:
        config = config.model_copy(update={"duration_h": duration_override})
    t, states, _ = simulate(config, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(
        t,
        readouts,
        therapeutic_window=analysis_window_for(config),
        pk=pert.pk if pert else None,
    )
    unc = run_uncertainty(config.model_copy(update={"n_sensitivity_runs": 8}), pert, params, n_runs=8)
    return config, pert, params, raw, t, states, readouts, kpis, unc


def _benchmark_tab():
    st.subheader("Mechanistic illustration benchmarks")
    st.caption(
        "Files 01–17 are unchanged. Each one is a mechanistic illustration. "
        "Pre-existing files 01–05 keep their original numbers."
    )
    benches = _list_benchmarks()
    choice = st.selectbox(
        "Benchmark",
        list(benches),
        format_func=_bench_label,
        help="Default is the first file, healthy control. IDs 01–17 stay fixed.",
    )
    use_file = st.checkbox(
        "Use each benchmark’s configured duration",
        value=True,
        help="On by default. The file’s own duration and analysis window are used.",
    )
    duration = st.number_input(
        "Duration override (h)",
        min_value=1.0,
        value=8.0,
        disabled=use_file,
        help="Used only when the checkbox above is off. Default 8 h. Not a clinical interval.",
    )
    if st.button(
        "Run benchmark",
        key="run_bench",
        type="primary",
        help="Runs the selected file with its stored perturbation.",
    ):
        with st.spinner("Running benchmark…"):
            st.session_state["bench"] = _run_benchmark(benches[choice], use_file, duration)
            st.session_state["bench_name"] = choice
    packed = st.session_state.get("bench")
    if not packed:
        return
    config, pert, params, raw, t, states, readouts, kpis, unc = packed
    if raw.get("illustration") == "mechanistic" or raw.get("potency_note"):
        st.info(ILLUSTRATION_BANNER)
    if "03_ldopa" in st.session_state.get("bench_name", ""):
        st.caption(
            "File 03 keeps the pre-existing illustration numbers, including EC50 and Cmax. "
            "Those values were not re-verified as published PK."
        )
    if raw.get("pk_timing_notes"):
        st.caption(str(raw.get("pk_timing_notes")))
    _metric_cards(kpis, config)
    rows = []
    for key, unit in PRIMARY_KPIS:
        if key not in kpis:
            continue
        rows.append(
            {
                "metric": key,
                "point": fmt_sig(kpis[key]),
                "mean": fmt_sig(unc.mean.get(key)),
                "p10": fmt_sig(unc.p10.get(key)),
                "p50": fmt_sig(unc.p50.get(key)),
                "p90": fmt_sig(unc.p90.get(key)),
                "unit": unit,
                "window": window_label(config),
                "definition": READOUT_DOCS.get(key, ""),
            }
        )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    chart_title = str(raw.get("name", "benchmark")).replace("_", " ")
    st.pyplot(_readout_panel(t, readouts, chart_title), clear_figure=True)
    item = {"label": chart_title, "t": t, "y": readouts["beta_proxy"]}
    if unc.trajectory_p10:
        item["p10"] = unc.trajectory_p10["beta_proxy"]
        item["p90"] = unc.trajectory_p90["beta_proxy"]
    st.pyplot(_plot_series([item], "Beta proxy (dimensionless)", "Beta proxy"), clear_figure=True)
    st.markdown(
        evidence_caption(
            load_evidence_registry(),
            parameter_keys=["beta_proxy", "motor_readiness", "dopamine_tone"],
            evidence_ids=list(raw.get("evidence_ids") or []),
        )
    )
    st.pyplot(sensitivity_figure(unc.sensitivity_ranking, "Sensitivity ranking"), clear_figure=True)
    traj = trajectories_frame([{"name": raw.get("name", "benchmark"), "t": t, "states": states, "readouts": readouts}])
    c1, c2, c3 = st.columns(3)
    c1.download_button("Trajectories CSV", frame_to_csv(traj), file_name="benchmark_trajectories.csv", mime="text/csv", key="b_traj")
    c2.download_button("KPI table CSV", pd.DataFrame(rows).to_csv(index=False), file_name="benchmark_kpis.csv", mime="text/csv", key="b_kpi")
    c3.download_button(
        "Config JSON",
        json.dumps(
            {
                "simulation": config.model_dump(mode="json"),
                "perturbation": pert.model_dump(mode="json") if pert else None,
                "kpis": kpis,
            },
            indent=2,
        ),
        file_name="benchmark_config.json",
        mime="application/json",
        key="b_json",
    )


def _custom_tab(session: dict):
    st.subheader("Custom candidate")
    st.markdown(CUSTOM_DISCLAIMER)
    st.caption(ILLUSTRATIVE_PK_LABEL + " The form is prefilled. Assay PK is off.")
    mode = st.radio(
        "Input",
        ["Form", "YAML"],
        horizontal=True,
        help="Form uses the illustrative defaults. YAML is the same scenario in text.",
    )
    options = _baseline_options()
    baseline = st.selectbox(
        "Baseline",
        options,
        index=options.index(SESSION_DEFAULTS["baseline"]),
        format_func=_baseline_title,
        key="cu_base",
        help="Circuit state for this single candidate. Default is moderate PD depletion.",
    )
    st.caption(baseline_by_preset(baseline).summary)
    if mode == "YAML":
        from src.catalog import ILLUSTRATIVE_CUSTOM_YAML

        text = st.text_area(
            "MoleculePerturbation YAML",
            value=ILLUSTRATIVE_CUSTOM_YAML,
            height=280,
            help="Prefilled illustrative scenario. No EC50 or Cmax. Edit only with values you supply.",
        )
        payload = text
        form_spec = None
    else:
        choice = st.selectbox(
            "Mechanism illustration",
            [item.id for item in MECHANISM_ILLUSTRATIONS],
            format_func=lambda item_id: illustration_by_id(item_id).title,
            key="cu_ill",
            help="Fills target, mechanism, and illustrative max_effect.",
        )
        form_spec = illustration_by_id(choice)
        st.caption(form_spec.blurb)
        name = st.text_input(
            "Name",
            value="Illustrative custom candidate",
            help="Label only. Illustrative default, not a real compound name.",
        )
        nodes = [n.value for n in TargetNode]
        mechanisms = [m.value for m in Mechanism]
        target = st.selectbox(
            "Target",
            nodes,
            index=nodes.index(form_spec.target_node),
            key=f"cu_target_{choice}",
            help="Circuit node for the gain.",
        )
        mechanism = st.selectbox(
            "Mechanism",
            mechanisms,
            index=mechanisms.index(form_spec.mechanism),
            key=f"cu_mech_{choice}",
            help="Circuit-sign field. Read the blurb for the pharmacological class.",
        )
        max_effect = st.text_input(
            "max_effect",
            value=form_spec.max_effect,
            key=f"cu_me_{choice}",
            help="Unitless gain. Illustrative default, not a measured potency.",
        )
        onset = st.text_input("Onset h", value=ILLUSTRATIVE_ONSET_H, key="cu_on", help=PK_HELP + " Unit: hours.")
        tmax = st.text_input("Tmax h", value=ILLUSTRATIVE_TMAX_H, key="cu_tm", help=PK_HELP + " Unit: hours.")
        half_life = st.text_input(
            "Half-life h",
            value=ILLUSTRATIVE_HALF_LIFE_H,
            key="cu_hl",
            help=PK_HELP + " Unit: hours. Default 2.",
        )
        ia = st.text_input(
            "Intrinsic activity",
            value="0.5",
            key="cu_ia",
            help="Scenario fraction in (0, 1]. Required for a partial agonist. Not a measured potency.",
        )
        duration = st.text_input(
            "Duration (h)",
            value=ILLUSTRATIVE_DURATION_H,
            key="cu_dur",
            help="Scenario length in hours. Default 8. Illustrative default, not measured PK.",
        )
        use_assay = st.checkbox(
            "Supply assay PK (EC50, reference Cmax, ka)",
            value=False,
            key="cu_assay",
            help="Off by default. Turn on only to enter assay values you have. EC50 and Cmax are not invented.",
        )
        ec50 = cmax = ka = ""
        if use_assay:
            ec50 = st.text_input("EC50 nM", value="", help="Nanomolar EC50 you supply. Required when this box is on.")
            cmax = st.text_input("Cmax nM", value="", help="Nanomolar Cmax you supply. Required when this box is on.")
            ka = st.text_input("ka (1/h, optional)", value="", help="Optional absorption rate in 1/h.")
        payload = {
            "name": name,
            "target_node": target,
            "mechanism": mechanism,
            "max_effect": max_effect,
            "onset": onset,
            "tmax": tmax,
            "half_life": half_life,
            "intrinsic_activity": ia,
            "ec50": ec50,
            "cmax": cmax,
            "ka": ka,
            "use_assay": use_assay,
            "duration": duration,
        }
    if st.button(
        "Run custom candidate",
        type="primary",
        help="Runs the prefilled candidate, or your edits, on the normalized path unless assay PK is on.",
    ):
        if not _window_ok(session):
            st.warning("Analysis window end (h) must be greater than the start. The candidate was not run.")
        else:
            try:
                if isinstance(payload, dict):
                    pert, problems = custom_form_perturbation(
                        name=payload["name"],
                        target_node=payload["target_node"],
                        mechanism=payload["mechanism"],
                        max_effect=payload["max_effect"],
                        onset=payload["onset"],
                        tmax=payload["tmax"],
                        half_life=payload["half_life"],
                        intrinsic_activity=payload["intrinsic_activity"],
                        ec50=payload["ec50"],
                        cmax=payload["cmax"],
                        ka=payload["ka"],
                        use_assay=payload["use_assay"],
                    )
                    if problems or pert is None:
                        st.warning("Check the custom candidate form: " + "; ".join(problems) + ". The candidate was not run.")
                        pert = None
                    duration_text = payload["duration"]
                else:
                    pert = parse_perturbation_yaml(payload)
                    duration_text = ILLUSTRATIVE_DURATION_H
                if pert is not None:
                    from src.study_design import resolve_scenario_duration

                    resolved = resolve_scenario_duration(duration_text)
                    if resolved.problem is not None:
                        st.warning(format_field_problems([resolved.problem]))
                    else:
                        cfg = _config_from_session(session, baseline, resolved.hours, 8)
                        params = load_baseline_params(cfg.baseline_preset)
                        t, states, _ = simulate(cfg, pert, params)
                        readouts = compute_readouts(t, states, params)
                        kpis = summarize_kpis(
                            t, readouts, therapeutic_window=analysis_window_for(cfg), pk=pert.pk
                        )
                        unc = run_uncertainty(cfg, pert, params, n_runs=8)
                        evidence = list(form_spec.evidence_ids) if form_spec is not None else []
                        st.session_state["custom"] = (cfg, pert, t, states, readouts, kpis, unc, evidence)
            except Exception as exc:  # noqa: BLE001
                st.error(str(exc))
    packed = st.session_state.get("custom")
    if not packed:
        return
    cfg, pert, t, states, readouts, kpis, unc, evidence = packed
    if pert.pk.cmax_nM is None and not any(m.ec50_nM for m in pert.modulations):
        st.caption(POTENCY_NORMALIZED + " " + ILLUSTRATIVE_PK_LABEL)
    _metric_cards(kpis, cfg)
    st.pyplot(_readout_panel(t, readouts, pert.name), clear_figure=True)
    item = {
        "label": pert.name,
        "t": t,
        "y": readouts["beta_proxy"],
        "p10": unc.trajectory_p10["beta_proxy"],
        "p90": unc.trajectory_p90["beta_proxy"],
    }
    st.pyplot(_plot_series([item], "Beta proxy (dimensionless)", pert.name), clear_figure=True)
    st.markdown(
        evidence_caption(
            load_evidence_registry(),
            parameter_keys=["beta_proxy", "d1_direct", "dopamine_tone"],
            evidence_ids=evidence,
        )
    )
    traj = trajectories_frame([{"name": pert.name, "t": t, "states": states, "readouts": readouts}])
    st.download_button("Trajectories CSV", frame_to_csv(traj), file_name="custom_trajectories.csv", mime="text/csv", key="c_traj")
    st.download_button(
        "Config JSON",
        json.dumps(
            {"simulation": cfg.model_dump(mode="json"), "perturbation": pert.model_dump(mode="json"), "kpis": kpis},
            indent=2,
        ),
        file_name="custom_config.json",
        mime="application/json",
        key="c_json",
    )


def _getting_started_tab():
    st.markdown(getting_started_markdown())


def _limitations_tab():
    st.markdown(cou_and_limitations_markdown())
    st.subheader("Readout definitions")
    for key, text in READOUT_DOCS.items():
        st.markdown(f"**{key}.** {text}")
    st.subheader("Evidence registry")
    st.dataframe(pd.DataFrame(evidence_table_rows()), width="stretch", hide_index=True)


def _sidebar() -> dict:
    with st.sidebar:
        st.markdown("**Basal Ganglia Digital Twin**")
        st.caption(f"Version {__version__}")
        st.markdown(
            "Mechanistic circuit translation. Fields in the main forms start filled. "
            + ILLUSTRATIVE_PK_LABEL
        )
        seed = st.number_input(
            "Random seed",
            min_value=0,
            value=SESSION_DEFAULTS["seed"],
            step=1,
            help="Integer seed. Default 42. The same seed reproduces the same trajectory.",
        )
        noise = st.number_input(
            "Observation noise (dimensionless SD)",
            min_value=0.0,
            value=SESSION_DEFAULTS["noise_std"],
            step=0.01,
            help="Standard deviation added to the observed states. Default 0. Dimensionless.",
        )
        start = st.number_input(
            "Analysis window start (h)",
            min_value=0.0,
            value=SESSION_DEFAULTS["analysis_window_start_h"],
            step=0.1,
            help="Start of the window used for means, ON/OFF time, and the dyskinesia-risk fraction. Unit: model hours. Default 0.8.",
        )
        end = st.number_input(
            "Analysis window end (h)",
            min_value=0.1,
            value=SESSION_DEFAULTS["analysis_window_end_h"],
            step=0.1,
            help="End of that window. Unit: model hours. Default 3.5. Must be greater than the start.",
        )
        st.caption(
            "Seed, noise, and the analysis window apply to Study design and Custom candidate. "
            "Benchmarks keep the window stored in each file."
        )
        st.markdown("---")
        st.caption(
            "Non-claims: no affinity from structure, no BBB prediction, no ADMET or toxicity "
            "prediction, no clinical-efficacy claim, and no invented EC50, Cmax, or half-life. "
            "The full text is on Getting started and the last tab."
        )
    return {
        "seed": int(seed),
        "noise_std": noise,
        "analysis_window_start_h": start,
        "analysis_window_end_h": end,
    }


def main():
    st.set_page_config(page_title="Basal Ganglia Digital Twin", layout="wide")
    apply_mpl_style()
    st.markdown(APP_CSS, unsafe_allow_html=True)
    session = _sidebar()
    st.markdown(HERO_HTML, unsafe_allow_html=True)
    st.markdown(COU_BANNER)
    getting_tab, study_tab, bench_tab, custom_tab, limit_tab = st.tabs(
        [
            "Getting started",
            "Study design",
            "Benchmarks",
            "Custom candidate",
            "Context of Use & Limitations",
        ]
    )
    with getting_tab:
        _getting_started_tab()
    with study_tab:
        _study_design_tab(session)
    with bench_tab:
        _benchmark_tab()
    with custom_tab:
        _custom_tab(session)
    with limit_tab:
        _limitations_tab()
    st.divider()
    st.markdown(
        f'<div class="bgdt-footer">Basal Ganglia Digital Twin · version {__version__} · '
        "Dimensionless circuit proxies. Context of Use and the non-claims are on "
        "Getting started and the last tab.</div>",
        unsafe_allow_html=True,
    )


main()
