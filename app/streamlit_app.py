"""Basal Ganglia Digital Twin v1 — Streamlit walkthrough.

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
from src.circuit_model import simulate
from src.config_loader import load_baseline_params, load_benchmark, parse_perturbation_yaml
from src.evidence_registry import (
    cou_and_limitations_markdown,
    evidence_caption,
    evidence_table_rows,
    load_evidence_registry,
)
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
    preset_reference_perturbation,
    run_study_from_base,
)
from src.uncertainty import run_uncertainty

BENCH_DIR = ROOT / "configs" / "benchmarks"
COU_BANNER = (
    "**Context of Use:** Mechanistic circuit-level translation only. "
    "This tool does **not** predict de novo binding affinity, BBB permeability, "
    "peripheral ADMET/toxicity, or clinical efficacy without experimental validation."
)
CUSTOM_DISCLAIMER = (
    "**Custom / newly discovered molecules:** Affinity, BBB penetration, and toxicity "
    "are *user inputs or unknown*. The twin only translates your stated target "
    "modulation (node + mechanism + max_effect / EC50 if known) and PK exposure "
    "(onset, tmax, half-life, cmax, dosing) into circuit readouts. "
    "It does **not** invent binding from SMILES or structures."
)
EXAMPLE_PERTURBATION_YAML = """\
# Example placeholders only — not measured PK or potency.
name: My novel candidate
description: User-defined mechanistic scenario (assay/PK assumptions required).
notes: Affinity/BBB/toxicity are user inputs or unknown — not predicted by the twin. potency not specified; normalized exposure shape
pk:
  onset_h: 0.0
  tmax_h: 1.0
  half_life_h: 2.0
  n_doses: 1
modulations:
  - target_node: d1_direct
    mechanism: agonist
    max_effect: 0.7
    hill_coeff: 1.0
"""

CARD_KEYS = [
    ("motor_mean_window", "Motor readiness", "dimensionless proxy"),
    ("beta_power_change_pct", "Beta power change", "percent vs untreated reference"),
    ("off_time_h", "OFF-time proxy", "model hours in the analysis window"),
    ("therapeutic_index_margin", "Therapeutic index margin", "dimensionless proxy"),
]


def _list_benchmarks() -> dict[str, Path]:
    return {p.stem: p for p in sorted(BENCH_DIR.glob("*.yaml"))}


def _plot_series(study_or_rows, ylabel: str, key: str, title: str):
    apply_mpl_style()
    fig, ax = plt.subplots(figsize=(7.4, 3.4))
    rows = study_or_rows
    for row in rows:
        ax.plot(row["t"], row["y"], label=row["label"])
        if row.get("p10") is not None:
            ax.fill_between(row["t"], row["p10"], row["p90"], alpha=0.18)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    return fig


def _dose_response_figure(study):
    apply_mpl_style()
    levels = [d.level for d in study.doses]
    fig, ax = plt.subplots(figsize=(7.4, 3.4))
    ax.plot(levels, [d.kpis["beta_power_change_pct"] for d in study.doses], marker="o", label="Beta change %")
    ax.plot(levels, [d.kpis["motor_mean_window"] for d in study.doses], marker="o", label="Motor readiness")
    ax.plot(levels, [d.kpis["overdrive_peak"] for d in study.doses], marker="o", label="Overdrive peak")
    band = [d.level for d in study.doses if d.in_band]
    if band:
        ax.axvspan(min(band), max(band), color="#1F7A6B", alpha=0.12, label="Therapeutic-window band")
    ax.set_xlabel("Dose level (user scale)")
    ax.set_ylabel("(dimensionless proxy)")
    ax.set_title("Exposure sweep")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    return fig


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
                }
            )
    return pd.DataFrame(rows)


def _show_study(study, evidence_ids: list[str] | None = None):
    if study.illustration:
        st.info(ILLUSTRATION_BANNER)
    if study.potency_note:
        st.caption(study.potency_note)
    st.write(study.scaling_note)
    if study.inverted_u:
        st.warning(
            "Inverted-U: therapeutic index is lower at a higher exposure than at an interior dose."
        )
    ref = next((d for d in study.doses if abs(d.level - 1.0) < 1e-9), study.doses[len(study.doses) // 2])
    cols = st.columns(len(CARD_KEYS))
    window = window_label(study.config)
    for col, (key, label, unit) in zip(cols, CARD_KEYS):
        col.metric(label, fmt_sig(ref.kpis.get(key)))
        col.caption(f"{unit}. Window {window}. Reference dose {ref.label}.")
    st.dataframe(_kpi_table(study), width="stretch", hide_index=True)
    st.pyplot(
        _dose_response_figure(study),
        clear_figure=True,
    )
    series = []
    for dose in study.doses:
        item = {
            "label": dose.label,
            "t": dose.time_h,
            "y": dose.readouts["beta_proxy"],
        }
        if dose.uncertainty is not None and dose.uncertainty.trajectory_p10:
            item["p10"] = dose.uncertainty.trajectory_p10["beta_proxy"]
            item["p90"] = dose.uncertainty.trajectory_p90["beta_proxy"]
        series.append(item)
    st.pyplot(
        _plot_series(series, "Beta (dimensionless proxy)", "beta_proxy", "Beta proxy with P10–P90 where sampled"),
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
            {
                "name": d.label,
                "t": d.time_h,
                "states": d.states,
                "readouts": d.readouts,
            }
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


def _run_preset(preset_id: str):
    spec, config, pert, params, raw = preset_reference_perturbation(preset_id)
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


def _float_field(label: str, key: str, placeholder: str) -> str:
    return st.text_input(label, value="", key=key, placeholder=placeholder)


def _study_design_tab():
    st.subheader("Preclinical study design")
    st.write(
        "Enter 3–5 doses and the PK you supply, or start from a mechanistic illustration preset. "
        "Empty fields are intentional: this page does not fill in a real PK or potency for a novel molecule."
    )
    st.markdown("**Mechanistic illustration presets, not calibrated.** One click runs a relative-exposure sweep.")
    cols = st.columns(len(PRESET_SPECS))
    for col, spec in zip(cols, PRESET_SPECS):
        col.caption(spec["blurb"])
        if col.button(spec["title"], key=f"preset_{spec['id']}"):
            with st.spinner("Running illustration sweep…"):
                _run_preset(spec["id"])
    st.divider()
    st.markdown("**Or enter your own candidate.** Placeholders are examples, not measured values.")
    name = st.text_input("Candidate name", value="", placeholder="e.g. Candidate A — label only")
    question = st.text_input(
        "Study question",
        value="",
        placeholder="e.g. Which relative exposure lowers beta without overdrive?",
    )
    c1, c2, c3 = st.columns(3)
    target = c1.selectbox("Target node", [n.value for n in TargetNode], key="sd_target")
    mechanism = c2.selectbox("Mechanism", [m.value for m in Mechanism], key="sd_mech")
    max_effect = c3.text_input("max_effect at the reference dose", value="", placeholder="e.g. 0.8 — example only")
    ia = st.text_input(
        "Intrinsic activity (required for partial agonist)",
        value="",
        placeholder="e.g. 0.5 — scenario fraction, not a measured potency",
    )
    p1, p2, p3, p4 = st.columns(4)
    onset = p1.text_input("Onset (h)", value="", placeholder="e.g. 0 — example only")
    tmax = p2.text_input("Tmax (h)", value="", placeholder="e.g. 1.5 — example only, not a measured value")
    half_life = p3.text_input("Half-life (h)", value="", placeholder="e.g. 2 — example only, not a measured value")
    ka = p4.text_input("ka (1/h, optional)", value="", placeholder="leave blank to derive from Tmax")
    e1, e2, e3 = st.columns(3)
    ec50 = e1.text_input("EC50 (nM, optional)", value="", placeholder="leave blank if potency is unknown")
    cmax_ref = e2.text_input("Reference Cmax (nM)", value="", placeholder="only if you supplied an EC50")
    reference = e3.text_input("Reference dose level", value="", placeholder="e.g. 1 — the level that maps to max_effect")
    st.caption("Dose levels. Fill 3 to 5. Leave unused rows blank.")
    labels, levels = [], []
    for i in range(5):
        a, b = st.columns(2)
        lab = a.text_input(f"Dose {i+1} label", value="", key=f"dose_lab_{i}", placeholder="e.g. low")
        lev = b.text_input(f"Dose {i+1} level", value="", key=f"dose_lvl_{i}", placeholder="e.g. 0.5 — example only")
        if lab.strip() and lev.strip():
            labels.append(lab.strip())
            levels.append(lev.strip())
    baseline = st.selectbox("Baseline", [p.value for p in BaselinePreset], index=1, key="sd_base")
    duration = st.text_input("Duration (h)", value="", placeholder="e.g. 8 — scenario length, not a clinical interval")
    if st.button("Run entered study", type="primary"):
        try:
            if not name.strip():
                raise ValueError("Enter a candidate name.")
            if not (3 <= len(labels) <= 5):
                raise ValueError("Enter 3 to 5 doses (label and level).")
            me = float(max_effect)
            pk_kwargs = {
                "onset_h": float(onset),
                "tmax_h": float(tmax),
                "half_life_h": float(half_life),
                "n_doses": 1,
            }
            if ka.strip():
                pk_kwargs["ka_per_h"] = float(ka)
            mod = {
                "target_node": target,
                "mechanism": mechanism,
                "max_effect": me,
                "hill_coeff": 1.0,
            }
            use_hill = bool(ec50.strip())
            if use_hill:
                mod["ec50_nM"] = float(ec50)
                if not cmax_ref.strip():
                    raise ValueError("Hill path needs a reference Cmax that you supply.")
                pk_kwargs["cmax_nM"] = float(cmax_ref)
                ref_level = float(cmax_ref)
            else:
                if not reference.strip():
                    raise ValueError("Enter the reference dose level (the level that maps to max_effect).")
                ref_level = float(reference)
            if mechanism == Mechanism.PARTIAL_AGONIST.value:
                if not ia.strip():
                    raise ValueError("Partial agonist needs an intrinsic-activity fraction you supply.")
                mod["intrinsic_activity"] = float(ia)
            elif ia.strip():
                mod["intrinsic_activity"] = float(ia)
            pert = MoleculePerturbation(
                name=name.strip(),
                notes=POTENCY_NORMALIZED if not use_hill else "EC50 and Cmax are user-supplied inputs, not predictions.",
                pk=PKProfile(**pk_kwargs),
                modulations=[TargetModulation(**mod)],
            )
            cfg = SimulationConfig(
                baseline_preset=BaselinePreset(baseline),
                duration_h=float(duration) if duration.strip() else 8.0,
                dt_h=0.05,
                noise_std=0.0,
                seed=42,
                n_sensitivity_runs=6,
            )
            params = load_baseline_params(cfg.baseline_preset)
            study = run_study_from_base(
                cfg,
                pert,
                params,
                labels=labels,
                levels=[float(x) for x in levels],
                reference_level=ref_level,
                use_hill=use_hill,
                question=question.strip() or f"How does exposure change circuit proxies for {name.strip()}?",
                n_uncertainty=6,
                illustration=False,
                potency_note="" if use_hill else POTENCY_NORMALIZED,
            )
            st.session_state["study"] = study
            st.session_state["study_evidence"] = []
        except Exception as exc:  # noqa: BLE001 — show the validation message
            st.error(str(exc))
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
    st.caption("Labeled as mechanistic illustrations. Pre-existing files 01–05 keep their original numbers.")
    benches = _list_benchmarks()
    choice = st.selectbox("Benchmark", list(benches))
    use_file = st.checkbox("Use each benchmark’s configured duration", value=True)
    duration = st.number_input("Duration override (h)", min_value=1.0, value=8.0, disabled=use_file)
    if st.button("Run benchmark", key="run_bench"):
        with st.spinner("Running benchmark…"):
            st.session_state["bench"] = _run_benchmark(benches[choice], use_file, float(duration))
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
    cols = st.columns(len(CARD_KEYS))
    for col, (key, label, unit) in zip(cols, CARD_KEYS):
        col.metric(label, fmt_sig(kpis.get(key)))
        col.caption(f"{unit}. Window {window_label(config)}.")
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
    item = {"label": raw.get("name", "benchmark"), "t": t, "y": readouts["beta_proxy"]}
    if unc.trajectory_p10:
        item["p10"] = unc.trajectory_p10["beta_proxy"]
        item["p90"] = unc.trajectory_p90["beta_proxy"]
    st.pyplot(_plot_series([item], "Beta (dimensionless proxy)", "beta", "Beta proxy"), clear_figure=True)
    motor = {"label": "motor", "t": t, "y": readouts["motor_readiness"]}
    if unc.trajectory_p10:
        motor["p10"] = unc.trajectory_p10["motor_readiness"]
        motor["p90"] = unc.trajectory_p90["motor_readiness"]
    st.pyplot(
        _plot_series([motor], "Motor readiness (dimensionless proxy)", "motor", "Motor readiness"),
        clear_figure=True,
    )
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


def _custom_tab():
    st.subheader("Custom candidate")
    st.markdown(CUSTOM_DISCLAIMER)
    mode = st.radio("Input", ["YAML", "Form"], horizontal=True)
    baseline = st.selectbox("Baseline", [p.value for p in BaselinePreset], index=1, key="cu_base")
    if mode == "YAML":
        text = st.text_area("MoleculePerturbation YAML", value=EXAMPLE_PERTURBATION_YAML, height=280)
        payload = text
    else:
        st.caption("Fields start empty. Placeholders are examples, not measured values.")
        name = st.text_input("Name", value="", placeholder="Candidate label")
        target = st.selectbox("Target", [n.value for n in TargetNode], key="cu_target")
        mechanism = st.selectbox("Mechanism", [m.value for m in Mechanism], key="cu_mech")
        max_effect = st.text_input("max_effect", value="", placeholder="e.g. 0.7 — example only")
        ec50 = st.text_input("EC50 nM (optional)", value="", placeholder="leave blank if unknown")
        onset = st.text_input("Onset h", value="", key="cu_on", placeholder="e.g. 0 — example only")
        tmax = st.text_input("Tmax h", value="", key="cu_tm", placeholder="e.g. 1 — example only")
        half_life = st.text_input("Half-life h", value="", key="cu_hl", placeholder="e.g. 2 — example only")
        cmax = st.text_input("Cmax nM (optional)", value="", placeholder="only with an EC50 you supply")
        ia = st.text_input("Intrinsic activity", value="", key="cu_ia", placeholder="required for partial agonist")
        payload = {
            "name": name,
            "notes": POTENCY_NORMALIZED,
            "pk": {
                "onset_h": onset,
                "tmax_h": tmax,
                "half_life_h": half_life,
                "n_doses": 1,
            },
            "modulations": [
                {
                    "target_node": target,
                    "mechanism": mechanism,
                    "max_effect": max_effect,
                    "hill_coeff": 1.0,
                }
            ],
        }
        if ec50.strip():
            payload["modulations"][0]["ec50_nM"] = ec50
        if cmax.strip():
            payload["pk"]["cmax_nM"] = cmax
        if ia.strip():
            payload["modulations"][0]["intrinsic_activity"] = ia
    if st.button("Run custom candidate"):
        try:
            if isinstance(payload, dict):
                pk = payload["pk"]
                for key in ("onset_h", "tmax_h", "half_life_h"):
                    if str(pk[key]).strip() == "":
                        raise ValueError("Enter onset, Tmax, and half-life. They are not filled in for you.")
                    pk[key] = float(pk[key])
                if "cmax_nM" in pk:
                    pk["cmax_nM"] = float(pk["cmax_nM"])
                mod = payload["modulations"][0]
                if str(mod["max_effect"]).strip() == "":
                    raise ValueError("Enter max_effect.")
                mod["max_effect"] = float(mod["max_effect"])
                if "ec50_nM" in mod:
                    mod["ec50_nM"] = float(mod["ec50_nM"])
                if "intrinsic_activity" in mod:
                    mod["intrinsic_activity"] = float(mod["intrinsic_activity"])
                if not str(payload["name"]).strip():
                    raise ValueError("Enter a candidate name.")
                pert = MoleculePerturbation(**payload)
            else:
                pert = parse_perturbation_yaml(payload)
            cfg = SimulationConfig(
                baseline_preset=BaselinePreset(baseline),
                duration_h=8.0,
                dt_h=0.05,
                noise_std=0.0,
                seed=42,
                n_sensitivity_runs=8,
            )
            params = load_baseline_params(cfg.baseline_preset)
            t, states, _ = simulate(cfg, pert, params)
            readouts = compute_readouts(t, states, params)
            kpis = summarize_kpis(t, readouts, therapeutic_window=analysis_window_for(cfg), pk=pert.pk)
            unc = run_uncertainty(cfg, pert, params, n_runs=8)
            st.session_state["custom"] = (cfg, pert, t, states, readouts, kpis, unc)
        except Exception as exc:  # noqa: BLE001
            st.error(str(exc))
    packed = st.session_state.get("custom")
    if not packed:
        return
    cfg, pert, t, states, readouts, kpis, unc = packed
    if pert.pk.cmax_nM is None and not any(m.ec50_nM for m in pert.modulations):
        st.caption(POTENCY_NORMALIZED)
    cols = st.columns(len(CARD_KEYS))
    for col, (key, label, unit) in zip(cols, CARD_KEYS):
        col.metric(label, fmt_sig(kpis.get(key)))
        col.caption(unit)
    item = {"label": pert.name, "t": t, "y": readouts["beta_proxy"], "p10": unc.trajectory_p10["beta_proxy"], "p90": unc.trajectory_p90["beta_proxy"]}
    st.pyplot(_plot_series([item], "Beta (dimensionless proxy)", "beta", pert.name), clear_figure=True)
    st.markdown(evidence_caption(load_evidence_registry(), parameter_keys=["beta_proxy", "d1_direct", "dopamine_tone"]))
    traj = trajectories_frame([{"name": pert.name, "t": t, "states": states, "readouts": readouts}])
    st.download_button("Trajectories CSV", frame_to_csv(traj), file_name="custom_trajectories.csv", mime="text/csv", key="c_traj")
    st.download_button(
        "Config JSON",
        json.dumps({"simulation": cfg.model_dump(mode="json"), "perturbation": pert.model_dump(mode="json"), "kpis": kpis}, indent=2),
        file_name="custom_config.json",
        mime="application/json",
        key="c_json",
    )


def _limitations_tab():
    st.markdown(cou_and_limitations_markdown())
    st.subheader("Readout definitions")
    for key, text in READOUT_DOCS.items():
        st.markdown(f"**{key}.** {text}")
    st.subheader("Evidence registry")
    st.dataframe(pd.DataFrame(evidence_table_rows()), width="stretch", hide_index=True)


def main():
    st.set_page_config(page_title="Basal Ganglia Digital Twin v1", layout="wide")
    apply_mpl_style()
    st.title("Basal Ganglia Digital Twin v1")
    st.markdown(COU_BANNER)
    study_tab, bench_tab, custom_tab, limit_tab = st.tabs(
        ["Study design", "Benchmarks", "Custom candidate", "Context of Use & Limitations"]
    )
    with study_tab:
        _study_design_tab()
    with bench_tab:
        _benchmark_tab()
    with custom_tab:
        _custom_tab()
    with limit_tab:
        _limitations_tab()
    st.divider()
    st.caption(
        f"Basal Ganglia Digital Twin v1 · version {__version__} · "
        "Context of Use & Limitations is the last tab. Dimensionless proxies only."
    )


main()
