"""Study design, new benchmarks, steady-state init, and neutral attribution."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.circuit_model import STATE_NAMES, simulate, untreated_steady_state
from src.config_loader import load_benchmark
from src.dose_response import flag_inverted_u, run_dose_response
from src.evidence_registry import (
    DEFAULT_EVIDENCE,
    evidence_caption,
    load_evidence_registry,
    references_markdown,
)
from src.pk_model import compute_perturbation_gains, exposure_curve, time_vector
from src.readouts import analysis_window_for, compute_readouts, summarize_kpis
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
    POTENCY_NORMALIZED,
    PRESET_SPECS,
    expand_dose_levels,
    preset_reference_perturbation,
    run_study_from_base,
)

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "configs" / "benchmarks"
EXACT_NON_CLAIMS = [
    "Does not predict de novo chemical binding affinity from SMILES/structures.",
    "Does not predict blood-brain barrier permeability.",
    "Does not predict peripheral ADMET/toxicity.",
    "Does not claim clinical efficacy without experimental validation.",
    "Does not invent EC50/Cmax/half-life for newly discovered molecules; user must supply mechanistic and PK assumptions.",
]


def _kpis(path: Path):
    config, pert, params, raw = load_benchmark(path)
    t, states, _ = simulate(config, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(
        t,
        readouts,
        therapeutic_window=analysis_window_for(config),
        pk=pert.pk if pert else None,
    )
    t0, s0, _ = simulate(config, None, params)
    base = summarize_kpis(
        t0,
        compute_readouts(t0, s0, params),
        therapeutic_window=analysis_window_for(config),
    )
    return config, pert, params, raw, t, states, readouts, kpis, base


def test_untreated_start_matches_steady_state_and_is_flat():
    config, _, params, _ = load_benchmark(BENCH / "02_untreated_pd.yaml")
    y0 = untreated_steady_state(params)
    t, states, _ = simulate(config, None, params)
    for i, name in enumerate(STATE_NAMES):
        assert states[name][0] == pytest.approx(float(y0[i]), abs=1e-6)
        assert float(np.max(np.abs(states[name] - states[name][0]))) < 1e-3
    assert t[0] == 0.0


def test_on_and_off_partition_the_analysis_window():
    config, pert, params, raw, t, states, readouts, kpis, _ = _kpis(
        BENCH / "06_ldopa_normalized_reference.yaml"
    )
    start, end = analysis_window_for(config)
    mask = (t >= start) & (t <= end)
    dt = float(np.median(np.diff(t)))
    expected = float(np.sum(mask) * dt)
    assert kpis["on_time_h"] + kpis["off_time_h"] == pytest.approx(expected, rel=1e-6, abs=1e-6)


def test_new_benchmarks_have_no_asserted_ec50_or_cmax():
    for path in sorted(BENCH.glob("*.yaml")):
        if path.name[:2] < "06":
            continue
        _, pert, _, raw, _, _, _, _, _ = _kpis(path)
        assert pert is not None
        assert pert.pk.cmax_nM is None
        assert all(m.ec50_nM is None for m in pert.modulations)
        assert POTENCY_NORMALIZED in (raw.get("potency_note") or "")


def test_negative_control_d2_antagonist_worsens_beta():
    _, _, _, _, _, _, _, kpis, base = _kpis(BENCH / "08_d2_antagonist_haloperidol_like.yaml")
    assert kpis["beta_mean_window"] > base["beta_mean_window"]
    assert kpis["motor_mean_window"] <= base["motor_mean_window"] + 1e-9


def test_d2_agonist_and_a2a_and_gpr6_lower_beta():
    for name in (
        "07_d2_agonist_pramipexole_like.yaml",
        "11_a2a_istradefylline_like.yaml",
        "17_gpr6_inverse_agonist.yaml",
        "12_amantadine_nmda_stn.yaml",
        "13_stn_channel_block_proxy.yaml",
    ):
        _, _, _, _, _, _, _, kpis, base = _kpis(BENCH / name)
        assert kpis["beta_mean_window"] < base["beta_mean_window"]


def test_comt_envelope_lengthens_horizon_on_time_on_the_same_record():
    cfg, long_pert, params, _ = load_benchmark(BENCH / "09_comt_opicapone_like_adjunct.yaml")
    _, short_pert, _, _ = load_benchmark(BENCH / "06_ldopa_normalized_reference.yaml")

    def horizon(pert):
        t, states, _ = simulate(cfg, pert, params)
        k = summarize_kpis(t, compute_readouts(t, states, params), therapeutic_window=analysis_window_for(cfg), pk=pert.pk)
        return k["on_time_horizon_h"]

    assert horizon(long_pert) > horizon(short_pert)


def test_repeat_dosing_exposure_exceeds_one_pulse():
    _, pert, _, _, t, _, _, _, _ = _kpis(BENCH / "14_ldopa_repeat_dosing.yaml")
    trap = getattr(np, "trapezoid", None) or np.trapz
    multi = exposure_curve(
        t, pert.pk.onset_h, pert.pk.tmax_h, pert.pk.half_life_h,
        n_doses=pert.pk.n_doses, dosing_interval_h=pert.pk.dosing_interval_h,
    )
    single = exposure_curve(
        t, pert.pk.onset_h, pert.pk.tmax_h, pert.pk.half_life_h, n_doses=1,
    )
    assert float(trap(multi, t)) > float(trap(single, t))


def test_partial_agonist_and_pam_are_not_the_same_trajectory():
    _, partial, params, _, t, _, r_partial, _, _ = _kpis(BENCH / "15_d1_partial_agonist_tavapadon_like.yaml")
    cfg, pam, _, _ = load_benchmark(BENCH / "16_d1_pam_glovadalen_like.yaml")
    _, states, _ = simulate(cfg, pam, params)
    r_pam = compute_readouts(t, states, params) if len(t) == len(states["d1"]) else compute_readouts(*simulate(cfg, pam, params)[:2], params)
    # Recompute on each file's own time base.
    t_p, s_p, _ = simulate(load_benchmark(BENCH / "15_d1_partial_agonist_tavapadon_like.yaml")[0], partial, params)
    t_m, s_m, _ = simulate(cfg, pam, params)
    assert not np.allclose(s_p["d1"], s_m["d1"])
    assert partial.modulations[0].intrinsic_activity == pytest.approx(0.5)
    gains = compute_perturbation_gains(
        time_vector(4.0, 0.1), pam, endogenous_tone=params.dopamine_tone
    )
    agonist = pam.model_copy(deep=True)
    agonist.modulations[0].mechanism = Mechanism.AGONIST
    gains_ag = compute_perturbation_gains(
        time_vector(4.0, 0.1), agonist, endogenous_tone=params.dopamine_tone
    )
    assert float(np.max(gains["d1_direct"])) < float(np.max(gains_ag["d1_direct"]))


def test_ka_changes_the_exposure_shape():
    t = time_vector(6.0, 0.05)
    slow = exposure_curve(t, 0.0, 1.0, 2.0, ka_per_h=0.4)
    fast = exposure_curve(t, 0.0, 1.0, 2.0, ka_per_h=3.0)
    assert not np.allclose(slow, fast)


def test_dose_response_normalized_path_and_inverted_u_flag():
    cfg, pert, params, _ = load_benchmark(BENCH / "06_ldopa_normalized_reference.yaml")
    result = run_dose_response(cfg, pert, params, [0.25, 0.5, 1.0, 2.0])
    assert result.scaling_mode == "normalized_max_effect"
    assert any(POTENCY_NORMALIZED in note for note in result.notes)
    assert result.metrics["beta_suppression_pct"].shape == (4,)
    assert flag_inverted_u(np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.4, 0.12]))
    assert not flag_inverted_u(np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.2, 0.3]))


def test_study_design_accepts_three_to_five_doses_and_writes_a_brief(tmp_path):
    spec, config, pert, params, raw = preset_reference_perturbation("ldopa")
    with pytest.raises(ValueError):
        expand_dose_levels(pert, ["a", "b"], [1, 2], 1.0, False)
    study = run_study_from_base(
        config,
        pert,
        params,
        labels=["low", "mid", "high"],
        levels=[0.5, 1.0, 2.0],
        reference_level=1.0,
        use_hill=False,
        question="Which relative exposure is inside the illustrated band?",
        n_uncertainty=4,
        illustration=True,
        evidence_ids=list(raw.get("evidence_ids") or []),
        potency_note=POTENCY_NORMALIZED,
    )
    assert len(study.doses) == 3
    assert POTENCY_NORMALIZED in study.scaling_note
    md = preclinical_brief_markdown(study)
    assert "Context of Use & Limitations" in md
    assert "Basal Ganglia Digital Twin v1" in md
    for claim in EXACT_NON_CLAIMS:
        assert claim in md
    assert "Suggestion, not a protocol" in md
    assert "Hao Li" not in md and "Li & Jin" not in md and "Li H" not in md
    pdf = write_preclinical_brief_pdf(study, tmp_path / "brief.pdf")
    blob = pdf.read_bytes()
    assert blob.startswith(b"%PDF")
    assert b"Hao Li" not in blob


def test_registry_size_identifiers_and_exact_non_claims():
    reg = load_evidence_registry()
    citations = reg["citations"]
    assert len(citations) >= 30
    assert reg["non_claims"] == EXACT_NON_CLAIMS
    assert DEFAULT_EVIDENCE["non_claims"] == EXACT_NON_CLAIMS
    for item in citations:
        assert item.get("doi") or item.get("pmid")
        blob = " ".join(str(item.get(k, "")) for k in ("citation", "authors", "id", "assumption"))
        assert "Hao Li" not in blob
        assert "Li & Jin" not in blob
        assert "Li H" not in blob
    refs = references_markdown(reg)
    assert "elife_action_selection_2023" in refs
    assert "Hao Li" not in refs
    caption = evidence_caption(reg, evidence_ids=["elife_action_selection_2023"])
    assert "10.7554/eLife.87644" in caption
    assert "Hao Li" not in caption
    assert "Li H" not in caption


def test_rendered_surfaces_omit_owner_identity():
    paths = [
        ROOT / "README.md",
        ROOT / "app" / "streamlit_app.py",
        ROOT / "src" / "report_generator.py",
        ROOT / "src" / "evidence_registry.py",
        ROOT / "configs" / "evidence_registry.yaml",
    ]
    for path in paths:
        text = path.read_text(encoding="utf-8")
        assert "Hao Li" not in text, path
        assert "Li & Jin" not in text, path
        assert "li_jin_2023" not in text, path
        assert "Li H," not in text and "Li H." not in text, path


def test_presets_are_labeled_illustrations_without_ec50():
    assert {s["id"] for s in PRESET_SPECS} >= {"ldopa", "d1_pam", "d1_partial", "d2_modulator", "gpr6"}
    for spec in PRESET_SPECS:
        _, _, pert, _, raw = preset_reference_perturbation(spec["id"])
        assert POTENCY_NORMALIZED in (raw.get("potency_note") or "")
        assert pert.pk.cmax_nM is None
        assert all(m.ec50_nM is None for m in pert.modulations)


def test_user_hill_path_keeps_ec50_fixed():
    base = MoleculePerturbation(
        name="user hill",
        pk=PKProfile(onset_h=0.0, tmax_h=1.0, half_life_h=2.0, cmax_nM=100.0),
        modulations=[
            TargetModulation(
                target_node=TargetNode.D1_DIRECT,
                mechanism=Mechanism.AGONIST,
                max_effect=0.8,
                ec50_nM=40.0,
            )
        ],
    )
    rows = expand_dose_levels(base, ["low", "mid", "high"], [50, 100, 200], 100, True)
    assert rows[0][2].pk.cmax_nM == pytest.approx(50)
    assert rows[2][2].pk.cmax_nM == pytest.approx(200)
    assert all(r[2].modulations[0].ec50_nM == pytest.approx(40) for r in rows)
    cfg = SimulationConfig(baseline_preset=BaselinePreset.PD_DEPLETED, duration_h=4.0, dt_h=0.1, n_sensitivity_runs=2)
    from src.config_loader import load_baseline_params

    study = run_study_from_base(
        cfg,
        base,
        load_baseline_params(BaselinePreset.PD_DEPLETED),
        labels=["low", "mid", "high"],
        levels=[50, 100, 200],
        reference_level=100,
        use_hill=True,
        n_uncertainty=0,
    )
    assert "Hill path" in study.scaling_note
    assert "Hao Li" not in study.suggestion
