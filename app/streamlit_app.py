"""Basal Ganglia Digital Twin — interactive Streamlit walkthrough UI."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.circuit_model import simulate
from src.config_loader import load_baseline_params, load_benchmark, parse_perturbation_yaml
from src.evidence_registry import (
    cou_and_limitations_markdown,
    evidence_table_rows,
    load_evidence_registry,
    references_markdown,
)
from src.readouts import compute_readouts, summarize_kpis
from src.report_generator import generate_decision_brief
from src.schemas import (
    BaselinePreset,
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    SimulationConfig,
    TargetModulation,
    TargetNode,
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
name: My novel candidate
description: User-defined mechanistic scenario (assay/PK assumptions required).
notes: Affinity/BBB/toxicity are user inputs or unknown — not predicted by the twin.
pk:
  onset_h: 0.2
  tmax_h: 1.0
  half_life_h: 2.5
  cmax_nM: 150.0
  n_doses: 1
modulations:
  - target_node: d1_direct
    mechanism: agonist
    max_effect: 0.7
    ec50_nM: 50.0
    hill_coeff: 1.0
"""


def _list_benchmarks() -> dict[str, Path]:
    return {p.stem: p for p in sorted(BENCH_DIR.glob("*.yaml"))}


def _run_scenario(config: SimulationConfig, pert, params):
    t, states, gains = simulate(config, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(t, readouts)
    return t, states, readouts, kpis, gains


def _candidate_from_sidebar(prefix: str) -> MoleculePerturbation | None:
    enabled = st.sidebar.checkbox(f"Enable candidate {prefix}", value=(prefix == "A"), key=f"{prefix}_en")
    if not enabled:
        return None
    name = st.sidebar.text_input(f"{prefix} name", value=f"Candidate {prefix}", key=f"{prefix}_name")
    target = st.sidebar.selectbox(
        f"{prefix} target_node",
        options=[n.value for n in TargetNode],
        index=0 if prefix == "A" else 2,
        key=f"{prefix}_tgt",
    )
    mech = st.sidebar.selectbox(
        f"{prefix} mechanism",
        options=[m.value for m in Mechanism],
        index=0,
        key=f"{prefix}_mech",
    )
    max_effect = st.sidebar.slider(
        f"{prefix} max_effect", 0.0, 2.0, 0.8 if prefix == "A" else 0.6, 0.05, key=f"{prefix}_me"
    )
    use_ec50 = st.sidebar.checkbox(f"{prefix} supply EC50 (assay)", value=True, key=f"{prefix}_use_ec50")
    ec50 = None
    if use_ec50:
        ec50 = float(
            st.sidebar.number_input(f"{prefix} ec50_nM", 1.0, 2000.0, 80.0, 5.0, key=f"{prefix}_ec50")
        )
    st.sidebar.caption(f"{prefix} PKProfile (empirical exposure — not ADMET prediction)")
    onset = st.sidebar.number_input(f"{prefix} onset_h", 0.0, 4.0, 0.2, 0.05, key=f"{prefix}_onset")
    tmax = st.sidebar.number_input(f"{prefix} tmax_h", 0.1, 6.0, 1.0, 0.1, key=f"{prefix}_tmax")
    half = st.sidebar.number_input(f"{prefix} half_life_h", 0.1, 12.0, 2.0, 0.1, key=f"{prefix}_hl")
    cmax = st.sidebar.number_input(f"{prefix} cmax_nM", 1.0, 5000.0, 200.0, 10.0, key=f"{prefix}_cmax")
    return MoleculePerturbation(
        name=name,
        description=f"User-defined candidate {prefix}",
        notes=(
            "UI-specified empirical PK/potency — not structure-derived. "
            "Affinity/BBB/toxicity are user inputs or unknown."
        ),
        pk=PKProfile(
            onset_h=float(onset),
            tmax_h=float(tmax),
            half_life_h=float(half),
            cmax_nM=float(cmax),
            n_doses=1,
        ),
        modulations=[
            TargetModulation(
                target_node=TargetNode(target),
                mechanism=Mechanism(mech),
                max_effect=float(max_effect),
                ec50_nM=ec50,
                hill_coeff=1.0,
            )
        ],
    )


def _load_yaml_candidate(label: str) -> MoleculePerturbation | None:
    st.sidebar.markdown(f"**{label}**")
    uploaded = st.sidebar.file_uploader(
        f"{label} YAML file", type=["yaml", "yml"], key=f"{label}_upload"
    )
    paste = st.sidebar.text_area(
        f"{label} YAML paste",
        value="" if label == "YAML B" else EXAMPLE_PERTURBATION_YAML,
        height=180 if label == "YAML A" else 100,
        key=f"{label}_paste",
    )
    text = None
    if uploaded is not None:
        text = uploaded.read().decode("utf-8")
    elif paste.strip():
        text = paste
    if not text:
        return None
    try:
        return parse_perturbation_yaml(text)
    except Exception as exc:  # noqa: BLE001 — show validation to user
        st.sidebar.error(f"{label} parse error: {exc}")
        return None


def _reference_runs(
    preset: str, seed: int, n_sens: int, duration: float
) -> list[tuple[str, SimulationConfig, object, object]]:
    """Healthy / untreated PD / L-DOPA reference arms for custom comparisons."""
    refs: list[tuple[str, SimulationConfig, object, object]] = []
    healthy = SimulationConfig(
        baseline_preset=BaselinePreset.HEALTHY,
        duration_h=float(duration),
        dt_h=0.05,
        noise_std=0.0,
        seed=int(seed),
        n_sensitivity_runs=int(n_sens),
    )
    refs.append(("Healthy (ref)", healthy, None, load_baseline_params(BaselinePreset.HEALTHY)))

    pd_cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=float(duration),
        dt_h=0.05,
        noise_std=0.0,
        seed=int(seed),
        n_sensitivity_runs=int(n_sens),
    )
    refs.append(("Untreated PD (ref)", pd_cfg, None, load_baseline_params(BaselinePreset.PD_DEPLETED)))

    ld_path = BENCH_DIR / "03_ldopa_standard.yaml"
    if ld_path.exists():
        cfg, pert, params, raw = load_benchmark(ld_path)
        cfg = cfg.model_copy(
            update={
                "seed": int(seed),
                "n_sensitivity_runs": int(n_sens),
                "duration_h": float(duration),
            }
        )
        refs.append((raw.get("name", "L-DOPA (ref)"), cfg, pert, params))
    return refs


def main() -> None:
    st.set_page_config(page_title="BG Digital Twin v1", layout="wide")
    st.title("Basal Ganglia Digital Twin (v1)")
    st.caption(
        "Preclinical therapeutic & mechanism simulator — CBGT twin grounded in "
        "classical + contemporary basal ganglia / PD / dopaminergic literature"
    )
    st.warning(COU_BANNER)

    benches = _list_benchmarks()

    with st.sidebar:
        st.header("Simulation controls")
        mode = st.radio(
            "Input mode",
            [
                "Benchmark presets",
                "Custom 1–2 candidates (form)",
                "Custom YAML upload / paste",
            ],
            index=0,
        )
        preset = st.selectbox(
            "Baseline preset (for custom runs)",
            options=[BaselinePreset.HEALTHY.value, BaselinePreset.PD_DEPLETED.value],
            index=1,
        )
        n_sens = st.slider("Uncertainty runs", 5, 60, 20, 5)
        seed = st.number_input("Random seed", 0, 10_000, 42, 1)
        duration = st.slider("Duration (h)", 2.0, 24.0, 8.0, 0.5)

        candidates: list[tuple[str, object, object, object]] = []
        if mode == "Benchmark presets":
            pick = st.multiselect(
                "Select 1–2 benchmarks to compare",
                options=list(benches.keys()),
                default=["03_ldopa_standard", "04_selective_d1_agonist"][:2],
                max_selections=2,
            )
            for key in pick:
                cfg, pert, params, raw = load_benchmark(benches[key])
                cfg = cfg.model_copy(
                    update={
                        "seed": int(seed),
                        "n_sensitivity_runs": int(n_sens),
                        "duration_h": float(duration),
                    }
                )
                candidates.append((raw.get("name", key), cfg, pert, params))
        else:
            st.info(CUSTOM_DISCLAIMER)
            include_refs = st.checkbox(
                "Also run Healthy / Untreated PD / L-DOPA references",
                value=True,
                help="Compare your custom candidate(s) against standard reference arms.",
            )
            params = load_baseline_params(preset)
            cfg = SimulationConfig(
                baseline_preset=BaselinePreset(preset),
                duration_h=float(duration),
                dt_h=0.05,
                noise_std=0.0,
                seed=int(seed),
                n_sensitivity_runs=int(n_sens),
            )
            if mode == "Custom 1–2 candidates (form)":
                st.subheader("Candidate A — mechanistic + PK assumptions")
                pert_a = _candidate_from_sidebar("A")
                st.subheader("Candidate B (optional)")
                pert_b = _candidate_from_sidebar("B")
                if pert_a:
                    candidates.append((pert_a.name, cfg, pert_a, params))
                if pert_b:
                    candidates.append((pert_b.name, cfg, pert_b, params))
            else:
                st.subheader("YAML MoleculePerturbation")
                st.caption(
                    "Paste or upload a bare perturbation YAML, or a full benchmark "
                    "file with a top-level `perturbation:` key."
                )
                pert_a = _load_yaml_candidate("YAML A")
                pert_b = _load_yaml_candidate("YAML B")
                if pert_a:
                    candidates.append((pert_a.name, cfg, pert_a, params))
                if pert_b:
                    candidates.append((pert_b.name, cfg, pert_b, params))
            if include_refs:
                for name, rcfg, rpert, rparams in _reference_runs(preset, seed, n_sens, duration):
                    candidates.append((name, rcfg, rpert, rparams))

        run = st.button("Run simulation", type="primary")

    if not run:
        st.info(
            "Configure the sidebar and click **Run simulation**. "
            "For novel molecules, supply target_node / mechanism / max_effect "
            "(+ EC50 if known) and a PKProfile — do not expect SMILES→affinity."
        )
        st.markdown(cou_and_limitations_markdown())
        with st.expander("References / Literature Basis (always shown)", expanded=True):
            st.markdown(references_markdown())
            st.dataframe(evidence_table_rows(load_evidence_registry()), width="stretch")
        with st.expander("Example MoleculePerturbation YAML"):
            st.code(EXAMPLE_PERTURBATION_YAML, language="yaml")
        return

    if len(candidates) < 1:
        st.error("Select at least one candidate, YAML perturbation, or enable reference arms.")
        return

    results = []
    for name, cfg, pert, params in candidates:
        t, states, readouts, kpis, gains = _run_scenario(cfg, pert, params)
        unc = run_uncertainty(cfg, pert, params, n_runs=int(n_sens))
        results.append(
            {
                "name": name,
                "config": cfg,
                "pert": pert,
                "params": params,
                "t": t,
                "states": states,
                "readouts": readouts,
                "kpis": kpis,
                "unc": unc,
            }
        )

    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "Trajectory & Readout Comparison",
            "Pathway Decomposition",
            "Sensitivity & Uncertainty",
            "Assumptions / Citations / Export",
        ]
    )

    with tab1:
        st.subheader("Primary readouts")
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
        for res in results:
            axes[0].plot(res["t"], res["readouts"]["beta_proxy"], label=res["name"])
            axes[1].plot(res["t"], res["readouts"]["motor_readiness"], label=res["name"])
        axes[0].set_title("LFP beta proxy")
        axes[1].set_title("Motor readiness")
        for ax in axes:
            ax.set_xlabel("Time (h)")
            ax.legend(fontsize=7)
        st.pyplot(fig)
        plt.close(fig)

        cols = st.columns(min(len(results), 4))
        for col, res in zip(cols, results):
            with col:
                st.metric(f"{res['name']} β window", f"{res['kpis']['beta_mean_window']:.3f}")
                st.metric(f"{res['name']} motor window", f"{res['kpis']['motor_mean_window']:.3f}")
                st.metric(f"{res['name']} overdrive flag", int(res["kpis"]["overdrive_flag"]))

    with tab2:
        st.subheader("Direct vs indirect/hyperdirect")
        fig, ax = plt.subplots(figsize=(8, 3.5))
        for res in results:
            ax.plot(res["t"], res["readouts"]["pathway_balance"], label=f"{res['name']} balance")
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Pathway balance")
        ax.legend(fontsize=7)
        st.pyplot(fig)
        plt.close(fig)

        fig2, axes = plt.subplots(1, 2, figsize=(10, 3.5))
        for res in results:
            axes[0].plot(res["t"], res["readouts"]["direct_pathway"], label=res["name"])
            axes[1].plot(res["t"], res["readouts"]["indirect_hyperdirect"], label=res["name"])
        axes[0].set_title("Direct pathway drive")
        axes[1].set_title("Indirect/hyperdirect drive")
        for ax in axes:
            ax.legend(fontsize=7)
            ax.set_xlabel("Time (h)")
        st.pyplot(fig2)
        plt.close(fig2)

    with tab3:
        st.subheader("P10–P90 uncertainty")
        for res in results:
            unc = res["unc"]
            st.markdown(f"**{res['name']}**")
            st.write(
                {
                    "motor_mean_window": {
                        "mean": unc.mean["motor_mean_window"],
                        "P10": unc.p10["motor_mean_window"],
                        "P90": unc.p90["motor_mean_window"],
                    },
                    "beta_mean_window": {
                        "mean": unc.mean["beta_mean_window"],
                        "P10": unc.p10["beta_mean_window"],
                        "P90": unc.p90["beta_mean_window"],
                    },
                }
            )
            st.write("Sensitivity ranking:", unc.sensitivity_ranking)

    with tab4:
        st.warning(CUSTOM_DISCLAIMER)
        st.markdown(cou_and_limitations_markdown())
        st.subheader("Assumptions & Literature Basis")
        st.markdown(references_markdown())
        st.dataframe(evidence_table_rows(load_evidence_registry()), width="stretch")

        st.subheader("Export decision brief")
        for res in results:
            safe = res["name"].replace(" ", "_").replace("(", "").replace(")", "")
            if st.button(f"Export brief: {res['name']}", key=f"export_{safe}"):
                out = ROOT / "outputs" / safe
                md = generate_decision_brief(
                    res["config"],
                    res["pert"],
                    res["params"],
                    out,
                    scenario_name=safe,
                    n_sensitivity_runs=int(n_sens),
                )
                st.success(f"Wrote {md}")
                cfg_path = out / f"{safe}_config.json"
                if cfg_path.exists():
                    st.download_button(
                        f"Download config JSON ({res['name']})",
                        data=cfg_path.read_text(encoding="utf-8"),
                        file_name=cfg_path.name,
                        mime="application/json",
                        key=f"dl_{safe}",
                    )


if __name__ == "__main__":
    main()
