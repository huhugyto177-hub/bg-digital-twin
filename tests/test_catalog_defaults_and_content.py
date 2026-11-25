"""Defaults, new baselines, new mechanisms, new readouts, and non-claims."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from src.catalog import (
    BASELINES,
    ILLUSTRATIVE_CUSTOM_YAML,
    ILLUSTRATIVE_PK_LABEL,
    MECHANISM_ILLUSTRATIONS,
    NEW_READOUT_KEYS,
    SESSION_DEFAULTS,
    custom_form_perturbation,
    default_illustration,
    run_illustration_on_baseline,
    run_illustrative_custom,
    run_illustrative_study,
    study_form_kwargs,
)
from src.circuit_model import simulate
from src.config_loader import load_baseline_params, load_benchmark, parse_perturbation_yaml
from src.evidence_registry import load_evidence_registry
from src.guide import getting_started_markdown
from src.pk_model import compute_perturbation_gains
from src.presentation import fmt_sig
from src.readouts import analysis_window_for, compute_readouts, summarize_kpis
from src.schemas import (
    BaselinePreset,
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    SimulationConfig,
    TargetModulation,
    TargetNode,
)
from src.study_design import parse_study_design_inputs
from tests.test_study_and_readouts import EXACT_NON_CLAIMS

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "configs" / "benchmarks"


def _untreated(preset: str):
    params = load_baseline_params(preset)
    cfg = SimulationConfig(
        baseline_preset=BaselinePreset(preset),
        duration_h=6.0,
        dt_h=0.05,
        noise_std=0.0,
        seed=42,
    )
    t, states, _ = simulate(cfg, None, params)
    kpis = summarize_kpis(t, compute_readouts(t, states, params), therapeutic_window=analysis_window_for(cfg))
    return params, states, kpis


def _finite_new(kpis: dict) -> None:
    for key in NEW_READOUT_KEYS:
        assert key in kpis
        assert math.isfinite(kpis[key]), key


def test_study_and_custom_forms_run_on_illustrative_defaults():
    parsed = parse_study_design_inputs(**study_form_kwargs())
    assert parsed.ok
    assert not parsed.use_hill
    assert parsed.ec50_nM is None
    assert parsed.cmax_nM is None
    assert parsed.ka_per_h is None
    assert parsed.labels == ["low", "mid", "high"]
    assert parsed.levels == pytest.approx([0.5, 1.0, 2.0])
    assert parsed.onset_h == pytest.approx(0.0)
    assert parsed.tmax_h == pytest.approx(1.5)
    assert parsed.half_life_h == pytest.approx(2.0)
    assert parsed.duration_h == pytest.approx(8.0)
    assert parsed.reference_level == pytest.approx(1.0)

    study = run_illustrative_study(n_uncertainty=0)
    assert len(study.doses) == 3
    assert POTENCY_NOTE_PRESENT(study.scaling_note)
    for dose in study.doses:
        _finite_new(dose.kpis)
        assert all(np.isfinite(series).all() for series in dose.states.values())

    spec = default_illustration()
    pert, problems = custom_form_perturbation(
        name="Illustrative custom candidate",
        target_node=spec.target_node,
        mechanism=spec.mechanism,
        max_effect=spec.max_effect,
        onset="0",
        tmax="1.5",
        half_life="2",
        intrinsic_activity="0.5",
        ec50="",
        cmax="",
        ka="",
        use_assay=False,
    )
    assert problems == []
    assert pert is not None
    assert pert.pk.cmax_nM is None
    assert pert.pk.ka_per_h is None
    assert all(mod.ec50_nM is None for mod in pert.modulations)

    _cfg, yaml_pert, _t, states, _readouts, kpis, _gains = run_illustrative_custom()
    assert yaml_pert.pk.cmax_nM is None
    assert all(mod.ec50_nM is None for mod in yaml_pert.modulations)
    _finite_new(kpis)
    assert all(np.isfinite(series).all() for series in states.values())
    assert ILLUSTRATIVE_PK_LABEL.split(",")[0] in ILLUSTRATIVE_CUSTOM_YAML or "not measured PK" in ILLUSTRATIVE_CUSTOM_YAML

    blocked, blocked_problems = custom_form_perturbation(
        name="Needs assay values",
        target_node="d1_direct",
        mechanism="agonist",
        max_effect="0.8",
        onset="0",
        tmax="1.5",
        half_life="2",
        intrinsic_activity="0.5",
        ec50="",
        cmax="",
        ka="",
        use_assay=True,
    )
    assert blocked is None
    assert blocked_problems
    assert parse_perturbation_yaml(ILLUSTRATIVE_CUSTOM_YAML).pk.half_life_h == pytest.approx(2.0)


def POTENCY_NOTE_PRESENT(note: str) -> bool:
    return "normalized" in note


def test_every_mechanism_illustration_parses_and_runs():
    for spec in MECHANISM_ILLUSTRATIONS:
        parsed = parse_study_design_inputs(**study_form_kwargs(spec.id))
        assert parsed.ok, spec.id
        assert not parsed.use_hill
        out = run_illustration_on_baseline(spec, "pd_depleted")
        pert = out["perturbation"]
        assert pert.pk.cmax_nM is None
        assert all(mod.ec50_nM is None for mod in pert.modulations)
        _finite_new(out["kpis"])
        assert all(np.isfinite(series).all() for series in out["states"].values())
        gain = out["exposure_gains"][spec.target_node]
        if spec.mechanism in {Mechanism.NAM.value, Mechanism.ANTAGONIST.value, Mechanism.CHANNEL_BLOCKER.value}:
            assert float(np.min(gain)) < 0.0
        else:
            assert float(np.max(gain)) > 0.0


def test_mechanism_directions_match_the_mapping():
    d1 = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "d1_agonist"), "pd_depleted")
    assert d1["kpis"]["direct_mean_window"] > d1["base"]["direct_mean_window"]

    mglu = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "mglur5_nam"), "pd_depleted")
    assert mglu["kpis"]["beta_mean_window"] < mglu["base"]["beta_mean_window"]

    htr = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "htr1a_agonist"), "pd_depleted")
    assert htr["kpis"]["direct_mean_window"] < htr["base"]["direct_mean_window"]

    m4 = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "m4_pam"), "pd_depleted")
    assert m4["kpis"]["direct_mean_window"] < m4["base"]["direct_mean_window"]

    gaba = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "gaba_a_gpi"), "pd_depleted")
    assert gaba["kpis"]["gpi_snr_output_rate"] < gaba["base"]["gpi_snr_output_rate"]

    a2a = run_illustration_on_baseline(next(s for s in MECHANISM_ILLUSTRATIONS if s.id == "a2a_antagonist"), "pd_depleted")
    assert a2a["kpis"]["beta_mean_window"] < a2a["base"]["beta_mean_window"]


def test_new_baselines_are_finite_and_directionally_ordered():
    rows = {}
    for info in BASELINES:
        params, states, kpis = _untreated(info.preset)
        rows[info.preset] = kpis
        _finite_new(kpis)
        for series in states.values():
            assert np.all(np.isfinite(series))
            assert float(np.min(series)) >= params.state_floor - 1e-8
            assert float(np.max(series)) <= params.state_ceil + 1e-8
        assert -1.05 <= kpis["direct_indirect_balance_index"] <= 1.05

    healthy, mild, depleted, advanced = (
        rows["healthy"],
        rows["pd_mild"],
        rows["pd_depleted"],
        rows["pd_advanced"],
    )
    assert healthy["beta_mean_window"] < mild["beta_mean_window"] < depleted["beta_mean_window"]
    assert depleted["beta_mean_window"] < advanced["beta_mean_window"]
    assert healthy["motor_mean_window"] > mild["motor_mean_window"] > depleted["motor_mean_window"]
    assert advanced["gpi_snr_output_rate"] > depleted["gpi_snr_output_rate"] > 0.0

    hd = rows["hd_indirect_loss"]
    assert hd["motor_mean_window"] > healthy["motor_mean_window"]
    assert hd["indirect_mean_window"] < healthy["indirect_mean_window"]
    assert hd["beta_mean_window"] < healthy["beta_mean_window"]

    dys = rows["dystonia_surround"]
    assert dys["motor_mean_window"] > healthy["motor_mean_window"]
    assert dys["indirect_mean_window"] < healthy["indirect_mean_window"]
    assert dys["beta_mean_window"] < depleted["beta_mean_window"]
    assert dys["overdrive_peak"] < 1.15

    pulse = MoleculePerturbation(
        name="illustrative dopamine-tone pulse",
        notes=ILLUSTRATIVE_PK_LABEL,
        pk=PKProfile(onset_h=0.0, tmax_h=1.5, half_life_h=2.0, n_doses=1),
        modulations=[
            TargetModulation(
                target_node=TargetNode.DOPAMINE_TONE,
                mechanism=Mechanism.AGONIST,
                max_effect=0.8,
            )
        ],
    )
    cfg = SimulationConfig(duration_h=8.0, dt_h=0.05, noise_std=0.0, seed=42)

    def overdrive(preset: str) -> float:
        params = load_baseline_params(preset)
        t, states, _ = simulate(cfg, pulse, params)
        kpis = summarize_kpis(t, compute_readouts(t, states, params), pk=pulse.pk)
        _finite_new(kpis)
        return kpis["overdrive_peak"]

    assert overdrive("lid_prone") > overdrive("pd_depleted")
    assert overdrive("lid_prone") >= 1.15
    gains = compute_perturbation_gains(
        np.linspace(0.0, 8.0, 40),
        pulse,
        endogenous_tone=0.24,
    )
    assert float(np.max(gains["dopamine_tone"])) > 0.0
    assert "essential_tremor" not in {item.value for item in BaselinePreset}


def test_benchmark_ids_01_through_17_are_unchanged_and_the_default_file_runs():
    names = sorted(path.name for path in BENCH.glob("*.yaml"))
    numbered = [name for name in names if name[:2].isdigit()]
    assert numbered[0].startswith("01_")
    assert any(name.startswith("17_") for name in numbered)
    assert len(numbered) == 17
    assert [int(name[:2]) for name in numbered] == list(range(1, 18))
    config, pert, params, _raw = load_benchmark(BENCH / numbered[0])
    t, states, _ = simulate(config, pert, params)
    kpis = summarize_kpis(t, compute_readouts(t, states, params), therapeutic_window=analysis_window_for(config))
    _finite_new(kpis)


def test_nonclaims_and_no_owner_name_on_new_surfaces():
    guide = getting_started_markdown()
    for claim in EXACT_NON_CLAIMS:
        assert claim in guide
    app = (ROOT / "app" / "streamlit_app.py").read_text(encoding="utf-8")
    assert "Getting started" in app
    assert ILLUSTRATIVE_PK_LABEL in app
    assert "Supply assay PK" in app
    registry = load_evidence_registry()
    assert registry["non_claims"] == EXACT_NON_CLAIMS
    blobs = [
        guide,
        app,
        (ROOT / "README.md").read_text(encoding="utf-8"),
        (ROOT / "src" / "catalog.py").read_text(encoding="utf-8"),
        (ROOT / "src" / "guide.py").read_text(encoding="utf-8"),
        (ROOT / "configs" / "evidence_registry.yaml").read_text(encoding="utf-8"),
    ]
    for info in BASELINES:
        blobs.append((ROOT / "configs" / _yaml_for(info.preset)).read_text(encoding="utf-8"))
    for blob in blobs:
        assert "Hao Li" not in blob
        assert "Li & Jin" not in blob
    for blob in (guide, app):
        assert "university" not in blob.lower()
        assert "@" not in blob
    assert SESSION_DEFAULTS["seed"] == 42
    assert SESSION_DEFAULTS["noise_std"] == 0.0
    assert SESSION_DEFAULTS["analysis_window_start_h"] == pytest.approx(0.8)
    assert SESSION_DEFAULTS["analysis_window_end_h"] == pytest.approx(3.5)
    assert "value=SESSION_DEFAULTS" in app


def test_display_collapses_numerical_dust():
    assert fmt_sig(5.32e-66) == "0"
    assert fmt_sig(-4.41e-08) == "0"
    assert fmt_sig(0.747) == "0.747"
    assert fmt_sig(float("nan")) == "n/a"


def _yaml_for(preset: str) -> str:
    return {
        "healthy": "baseline_healthy.yaml",
        "pd_mild": "baseline_pd_mild.yaml",
        "pd_depleted": "disease_pd_depleted.yaml",
        "pd_advanced": "baseline_pd_advanced.yaml",
        "lid_prone": "baseline_lid_prone.yaml",
        "hd_indirect_loss": "baseline_hd_indirect_loss.yaml",
        "dystonia_surround": "baseline_dystonia_surround.yaml",
    }[preset]
