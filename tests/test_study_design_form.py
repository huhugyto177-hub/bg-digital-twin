"""Numeric parsing for the Study design form, including a blank duration."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.config_loader import load_baseline_params
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
    default_scenario_duration_h,
    format_field_problems,
    parse_numeric_field,
    parse_study_design_inputs,
    parse_user_levels,
    resolve_scenario_duration,
    run_study_from_base,
)

ROOT = Path(__file__).resolve().parents[1]
RAW_FLOAT_ERROR = "could not convert string to float"


def _entered_study(**overrides):
    payload = dict(
        max_effect="0.8",
        onset="0",
        tmax="1.5",
        half_life="2",
        ka="",
        ec50="",
        cmax_ref="",
        reference="1",
        intrinsic_activity="",
        mechanism="agonist",
        dose_labels=["low", "mid", "high", "", ""],
        dose_levels=["0.5", "1", "2", "", ""],
        duration="",
    )
    payload.update(overrides)
    return parse_study_design_inputs(**payload)


def test_parse_numeric_field_blank_whitespace_and_invalid():
    blank = parse_numeric_field("")
    assert blank.blank and blank.value is None and not blank.invalid
    spaces = parse_numeric_field("  \t")
    assert spaces.blank and spaces.value is None
    assert parse_numeric_field(None).blank
    assert parse_numeric_field("   ").blank

    number = parse_numeric_field("  1.5 ")
    assert number.value == pytest.approx(1.5)
    assert not number.blank and not number.invalid
    assert parse_numeric_field("1e2").value == pytest.approx(100.0)

    for raw in ("abc", "8h", "1,5", "nan", "inf", "+inf", "-inf"):
        parsed = parse_numeric_field(raw)
        assert parsed.invalid, raw
        assert parsed.value is None and not parsed.blank

    with pytest.raises(ValueError) as exc:
        parse_user_levels(["a", "b", "c"], ["0.5", "", "2"], "", "1")
    assert RAW_FLOAT_ERROR not in str(exc.value)


def test_blank_duration_uses_existing_default_and_invalid_does_not():
    expected = SimulationConfig().duration_h
    assert expected == default_scenario_duration_h()
    assert expected == 8.0
    loader = (ROOT / "src" / "config_loader.py").read_text(encoding="utf-8")
    assert f'sim.get("duration_h", {expected})' in loader

    blank = resolve_scenario_duration("")
    assert blank.used_default and blank.problem is None
    assert blank.hours == expected
    whitespace = resolve_scenario_duration("   ")
    assert whitespace.used_default and whitespace.hours == expected

    explicit = resolve_scenario_duration(" 6.5 ")
    assert explicit.hours == pytest.approx(6.5)
    assert not explicit.used_default and explicit.problem is None

    invalid = resolve_scenario_duration("abc")
    assert invalid.hours is None and not invalid.used_default
    assert invalid.problem is not None
    assert invalid.problem.field == "Duration (h)" and invalid.problem.reason == "invalid"
    zero = resolve_scenario_duration("0")
    assert zero.hours is None and zero.problem is not None and zero.problem.reason == "positive"
    message = format_field_problems([invalid.problem, zero.problem])
    assert "Duration (h)" in message
    assert RAW_FLOAT_ERROR not in message

    entered = _entered_study(duration="")
    assert entered.ok and entered.used_default_duration
    assert entered.duration_h == expected
    assert entered.labels == ["low", "mid", "high"]
    assert entered.levels == pytest.approx([0.5, 1.0, 2.0])
    assert not entered.use_hill and entered.reference_level == pytest.approx(1.0)

    whitespace_form = _entered_study(duration=" \n ")
    assert whitespace_form.ok and whitespace_form.used_default_duration
    assert whitespace_form.duration_h == expected

    kept = _entered_study(duration="12")
    assert kept.ok and not kept.used_default_duration
    assert kept.duration_h == pytest.approx(12.0)

    bad_duration = _entered_study(duration="nope")
    assert not bad_duration.ok
    assert bad_duration.duration_h is None
    assert not bad_duration.used_default_duration
    assert [item.field for item in bad_duration.problems] == ["Duration (h)"]
    assert RAW_FLOAT_ERROR not in format_field_problems(bad_duration.problems)

    cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=entered.duration_h,
        dt_h=0.05,
        noise_std=0.0,
        seed=42,
        n_sensitivity_runs=2,
    )
    assert cfg.duration_h == expected
    assert cfg.baseline_preset == BaselinePreset.PD_DEPLETED
    pert = MoleculePerturbation(
        name="Candidate A",
        notes=POTENCY_NORMALIZED,
        pk=PKProfile(
            onset_h=entered.onset_h,
            tmax_h=entered.tmax_h,
            half_life_h=entered.half_life_h,
            n_doses=1,
        ),
        modulations=[
            TargetModulation(
                target_node=TargetNode.D1_DIRECT,
                mechanism=Mechanism.AGONIST,
                max_effect=entered.max_effect,
                hill_coeff=1.0,
            )
        ],
    )
    study = run_study_from_base(
        cfg,
        pert,
        load_baseline_params(BaselinePreset.PD_DEPLETED),
        labels=entered.labels,
        levels=entered.levels,
        reference_level=entered.reference_level,
        use_hill=False,
        n_uncertainty=0,
    )
    assert study.config is not None
    assert study.config.duration_h == expected
    assert len(study.doses) == 3


def test_blank_duration_without_a_default_names_the_field(monkeypatch):
    monkeypatch.setattr("src.study_design.default_scenario_duration_h", lambda: None)
    resolved = resolve_scenario_duration("")
    assert resolved.hours is None and not resolved.used_default
    assert resolved.problem is not None
    assert resolved.problem.field == "Duration (h)" and resolved.problem.reason == "blank"
    message = format_field_problems([resolved.problem])
    assert "Duration (h)" in message and "blank" in message
    assert RAW_FLOAT_ERROR not in message


def test_invalid_or_blank_numbers_name_fields_and_do_not_run():
    missing_level = _entered_study(
        dose_labels=["low", "mid", "high", "extra", ""],
        dose_levels=["0.5", "", "2", "abc", ""],
    )
    assert not missing_level.ok
    named = [item.field for item in missing_level.problems]
    assert named == ["Dose 2 level", "Dose 4 level"]
    assert "Duration (h)" not in named
    warning = format_field_problems(missing_level.problems)
    assert "Dose 2 level is blank" in warning
    assert "Dose 4 level is not a number" in warning
    assert RAW_FLOAT_ERROR not in warning

    missing_required = _entered_study(max_effect="  ", onset="soon", half_life="-1", ka="fast")
    required_names = [item.field for item in missing_required.problems]
    assert "max_effect at the reference dose" in required_names
    assert "Onset (h)" in required_names
    assert "Half-life (h)" in required_names
    assert "ka (1/h, optional)" in required_names
    assert RAW_FLOAT_ERROR not in format_field_problems(missing_required.problems)

    partial = _entered_study(mechanism="partial_agonist", intrinsic_activity="")
    assert [item.field for item in partial.problems] == [
        "Intrinsic activity (required for partial agonist)"
    ]

    hill_missing_cmax = _entered_study(ec50="40", cmax_ref="", reference="")
    assert not hill_missing_cmax.ok
    assert [item.field for item in hill_missing_cmax.problems] == ["Reference Cmax (nM)"]

    unused_optional = _entered_study(ka="  ", ec50="  ", cmax_ref="  ")
    assert unused_optional.ok
    assert unused_optional.ka_per_h is None
    assert not unused_optional.use_hill


def test_study_design_tab_parses_numbers_without_raw_float():
    source = (ROOT / "app" / "streamlit_app.py").read_text(encoding="utf-8")
    start = source.index("def _study_design_tab")
    end = source.index("\ndef _run_benchmark")
    body = source[start:end]
    assert "float(" not in body
    assert "st.error" not in body
    assert "parse_study_design_inputs" in body
    assert "default scenario duration" in body
    assert RAW_FLOAT_ERROR not in body
    assert "Hao Li" not in body
