"""Evidence registry coverage + custom MoleculePerturbation smoke tests."""

from __future__ import annotations

from pathlib import Path

from src.circuit_model import run_simulation
from src.config_loader import load_baseline_params, parse_perturbation_yaml
from src.evidence_registry import (
    cou_and_limitations_markdown,
    evidence_table_rows,
    load_evidence_registry,
    references_markdown,
)
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

ROOT = Path(__file__).resolve().parents[1]


def test_evidence_registry_has_multiple_broad_citations():
    reg = load_evidence_registry()
    rows = evidence_table_rows(reg)
    assert len(rows) >= 10, "Expected a broad multi-citation literature set"
    ids = {r["id"] for r in rows}
    # Classical CBGT / PD + beta + Triple-Control among many
    assert "albin_young_penney_1989" in ids or "delong_1990" in ids
    assert "brown_2003_beta" in ids or "hammond_2007_beta_osc" in ids
    assert "li_jin_2023" in ids
    assert "gerfen_surmeier_2011" in ids
    assert reg.get("calibration_note")
    assert any("affinity" in nc.lower() or "SMILES" in nc for nc in reg.get("non_claims", []))


def test_references_and_cou_markdown_include_literature():
    md = references_markdown()
    assert "## References / Literature Basis" in md
    assert "Citation" in md
    assert "li_jin_2023" in md
    cou = cou_and_limitations_markdown()
    assert "Context of Use" in cou
    assert "Calibration" in cou or "directional" in cou.lower()


def test_custom_molecule_perturbation_runs():
    pert = MoleculePerturbation(
        name="Custom test agonist",
        description="Programmatic novel-molecule-style assumptions",
        notes="Affinity/BBB unknown — user-supplied max_effect + PK only.",
        pk=PKProfile(onset_h=0.2, tmax_h=1.0, half_life_h=2.0, cmax_nM=100.0, n_doses=1),
        modulations=[
            TargetModulation(
                target_node=TargetNode.D1_DIRECT,
                mechanism=Mechanism.AGONIST,
                max_effect=0.6,
                ec50_nM=40.0,
                hill_coeff=1.0,
            )
        ],
    )
    cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=4.0,
        dt_h=0.1,
        noise_std=0.0,
        seed=7,
        n_sensitivity_runs=5,
    )
    params = load_baseline_params(BaselinePreset.PD_DEPLETED)
    result = run_simulation(cfg, pert, params)
    assert result.perturbation_name == "Custom test agonist"
    assert len(result.time_h) > 10
    assert "beta_mean_window" in result.kpis
    assert "motor_mean_window" in result.kpis


def test_parse_perturbation_yaml_and_run():
    text = """
name: YAML novel candidate
notes: User-supplied mechanistic + PK assumptions.
pk:
  onset_h: 0.1
  tmax_h: 0.8
  half_life_h: 1.5
  cmax_nM: 80.0
  n_doses: 1
modulations:
  - target_node: dopamine_tone
    mechanism: agonist
    max_effect: 0.5
    ec50_nM: 100.0
"""
    pert = parse_perturbation_yaml(text)
    assert pert.name == "YAML novel candidate"
    cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=3.0,
        dt_h=0.1,
        seed=3,
        n_sensitivity_runs=3,
    )
    params = load_baseline_params(BaselinePreset.PD_DEPLETED)
    result = run_simulation(cfg, pert, params)
    assert result.kpis["motor_mean_window"] >= 0.0


def test_decision_brief_includes_references_section(tmp_path):
    pert = MoleculePerturbation(
        name="Brief custom mol",
        pk=PKProfile(onset_h=0.2, tmax_h=1.0, half_life_h=2.0, cmax_nM=120.0),
        modulations=[
            TargetModulation(
                target_node=TargetNode.D2_INDIRECT,
                mechanism=Mechanism.AGONIST,
                max_effect=0.4,
                ec50_nM=60.0,
            )
        ],
    )
    cfg = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=3.0,
        dt_h=0.1,
        seed=11,
        n_sensitivity_runs=4,
    )
    params = load_baseline_params(BaselinePreset.PD_DEPLETED)
    md_path = generate_decision_brief(
        cfg, pert, params, tmp_path, scenario_name="custom_brief_test", n_sensitivity_runs=4
    )
    text = md_path.read_text(encoding="utf-8")
    assert "## References / Literature Basis" in text
    assert "li_jin_2023" in text
    assert "Context of Use" in text
    assert "user-supplied or unknown" in text.lower() or "Custom/novel" in text
