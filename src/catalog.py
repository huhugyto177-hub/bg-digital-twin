"""Baseline and mechanism catalogs, plus the illustrative form defaults.

Defaults are a shared scenario timebase. They are not measured PK, EC50, or
Cmax for any named molecule. Assay fields stay blank unless the caller turns
the assay path on and supplies values.
"""

from __future__ import annotations

from dataclasses import dataclass

from .circuit_model import simulate
from .config_loader import load_baseline_params
from .pk_model import compute_perturbation_gains
from .readouts import analysis_window_for, compute_readouts, summarize_kpis
from .schemas import (
    BaselinePreset,
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    SimulationConfig,
    TargetModulation,
    TargetNode,
)
from .study_design import parse_numeric_field, parse_study_design_inputs, run_study_from_base

ILLUSTRATIVE_PK_LABEL = "Illustrative default, not measured PK."
POTENCY_NOT_ASSERTED = "potency not specified; normalized exposure shape"

# Shared scenario timebase already used as placeholders in the study form.
ILLUSTRATIVE_ONSET_H = "0"
ILLUSTRATIVE_TMAX_H = "1.5"
ILLUSTRATIVE_HALF_LIFE_H = "2"
ILLUSTRATIVE_DURATION_H = "8"
ILLUSTRATIVE_REFERENCE = "1"
ILLUSTRATIVE_DOSE_LABELS = ["low", "mid", "high", "", ""]
ILLUSTRATIVE_DOSE_LEVELS = ["0.5", "1", "2", "", ""]

NEW_READOUT_KEYS = (
    "dyskinesia_risk_proxy",
    "gpi_snr_output_rate",
    "stn_burstiness_proxy",
    "direct_indirect_balance_index",
)

SESSION_DEFAULTS = {
    "seed": 42,
    "noise_std": 0.0,
    "analysis_window_start_h": 0.8,
    "analysis_window_end_h": 3.5,
    "baseline": BaselinePreset.PD_DEPLETED.value,
}


@dataclass(frozen=True)
class BaselineInfo:
    preset: str
    title: str
    summary: str
    evidence_ids: tuple[str, ...]


BASELINES: tuple[BaselineInfo, ...] = (
    BaselineInfo(
        preset="healthy",
        title="Healthy",
        summary=(
            "Reference circuit with intact dopamine tone and both striatal pathways. "
            "Directional prior from the classical direct/indirect skeleton."
        ),
        evidence_ids=("albin_1989", "delong_wichmann_2007"),
    ),
    BaselineInfo(
        preset="pd_mild",
        title="Mild PD",
        summary=(
            "Partial dopamine loss. Beta sits between healthy and the moderate depleted "
            "reference, and some motor readiness remains. Scenario levels, not measured "
            "concentrations. Directional prior: dopamine depletion tracks severity."
        ),
        evidence_ids=("bernheimer_1973", "kish_1988", "fearnley_1991", "delong_1990"),
    ),
    BaselineInfo(
        preset="pd_depleted",
        title="PD depleted (moderate reference)",
        summary=(
            "The moderate dopamine-depleted state used by benchmarks 01–17. "
            "Kept stable so those files do not change."
        ),
        evidence_ids=("delong_1990", "delong_wichmann_2007"),
    ),
    BaselineInfo(
        preset="pd_advanced",
        title="Advanced PD",
        summary=(
            "More severe dopamine loss than the moderate reference: higher GPi/SNr "
            "output and a higher beta proxy. Scenario levels ordered by the same "
            "severity prior, not a fitted staging scale."
        ),
        evidence_ids=("fearnley_1991", "kish_1988", "bernheimer_1973", "delong_1990"),
    ),
    BaselineInfo(
        preset="lid_prone",
        title="LID-prone",
        summary=(
            "Depleted dopamine tone with a sensitized direct pathway, so the same "
            "dopamine-tone pulse reaches the overdrive threshold more readily than "
            "on the moderate PD reference. Not a clinical dyskinesia rating."
        ),
        evidence_ids=("cenci_2007_lid", "calabresi_2010", "picconi_2003", "jenner_2008"),
    ),
    BaselineInfo(
        preset="hd_indirect_loss",
        title="Huntington-like indirect-pathway loss",
        summary=(
            "Early chorea-stage sketch: indirect-pathway (GPe-projecting) drive is "
            "reduced and dopamine tone is not the lesion. Motor readiness rises and "
            "the indirect term falls relative to healthy. Not a cell-count model."
        ),
        evidence_ids=("reiner_1988", "albin_1990_hd", "albin_1989"),
    ),
    BaselineInfo(
        preset="dystonia_surround",
        title="Dystonia-like lost surround",
        summary=(
            "Weaker indirect-pathway surround and a modestly easier thalamic release, "
            "without a parkinsonian beta elevation. A rate model cannot represent the "
            "pattern and somatosensory changes that a rate-only dystonia account leaves out. "
            "Essential tremor is not offered: it is not a cortico–basal rate imbalance."
        ),
        evidence_ids=("mink_2003_surround", "vitek_2002_dystonia"),
    ),
)


@dataclass(frozen=True)
class MechanismIllustration:
    id: str
    title: str
    target_node: str
    mechanism: str
    max_effect: str
    blurb: str
    evidence_ids: tuple[str, ...]
    intrinsic_activity: str = ""


MECHANISM_ILLUSTRATIONS: tuple[MechanismIllustration, ...] = (
    MechanismIllustration(
        id="d1_agonist",
        title="D1 direct-pathway agonist",
        target_node="d1_direct",
        mechanism="agonist",
        max_effect="0.8",
        blurb=(
            "Positive gain on the direct-pathway node. The 0.8 magnitude is an "
            "illustrative default, not a measured potency."
        ),
        evidence_ids=("gerfen_surmeier_2011", "surmeier_2007", "kravitz_2010"),
    ),
    MechanismIllustration(
        id="mglur5_nam",
        title="mGluR5 NAM",
        target_node="stn_hyperdirect",
        mechanism="nam",
        max_effect="0.6",
        blurb=(
            "mGluR5 negative allosteric modulation mapped as reduced glutamatergic "
            "drive on the STN/hyperdirect node. Illustrative gain, not a binding prediction."
        ),
        evidence_ids=("conn_2005_mglu",),
    ),
    MechanismIllustration(
        id="htr1a_agonist",
        title="5-HT1A agonist",
        target_node="dopamine_tone",
        mechanism="nam",
        max_effect="0.35",
        blurb=(
            "5-HT1A agonism is mapped as a negative dopamine-tone gain, because published "
            "work links it to less dopamine release from serotonin terminals. The NAM "
            "field is the circuit sign, not a claim that the ligand is a dopamine NAM."
        ),
        evidence_ids=("carta_2007_5ht",),
    ),
    MechanismIllustration(
        id="m4_pam",
        title="Muscarinic M4 PAM",
        target_node="d1_direct",
        mechanism="nam",
        max_effect="0.45",
        blurb=(
            "M4 positive allosteric modulation on direct-pathway neurons is mapped as "
            "a negative D1-node gain. The PAM pharmacology is described here; the model "
            "field carries the circuit sign."
        ),
        evidence_ids=("shen_2015_m4",),
    ),
    MechanismIllustration(
        id="gaba_a_gpi",
        title="GABA-A facilitation in GPi/SNr",
        target_node="gpi_snr",
        mechanism="nam",
        max_effect="0.5",
        blurb=(
            "A GABA-A agonist that inhibits GPi/SNr is a negative gain on that node, "
            "the same direction as focal muscimol inactivation. Illustrative gain, "
            "not a concentration and not a surgical plan."
        ),
        evidence_ids=("baron_2002_gaba_gpi",),
    ),
    MechanismIllustration(
        id="a2a_antagonist",
        title="Adenosine A2A antagonist",
        target_node="d2_indirect",
        mechanism="agonist",
        max_effect="0.7",
        blurb=(
            "A2A antagonism is illustrated on the indirect-pathway node with the same "
            "gain sign as benchmark 11. Potency is not specified. No published half-life "
            "is copied into the timebase."
        ),
        evidence_ids=("schwarzschild_2006", "mizuno_2013"),
    ),
)


def baseline_by_preset(preset: str) -> BaselineInfo:
    return next(item for item in BASELINES if item.preset == preset)


def illustration_by_id(illustration_id: str) -> MechanismIllustration:
    return next(item for item in MECHANISM_ILLUSTRATIONS if item.id == illustration_id)


def default_illustration() -> MechanismIllustration:
    return MECHANISM_ILLUSTRATIONS[0]


def study_form_kwargs(illustration_id: str | None = None, **overrides) -> dict:
    """String payload for ``parse_study_design_inputs``. Assay fields stay blank."""
    spec = illustration_by_id(illustration_id) if illustration_id else default_illustration()
    payload = dict(
        max_effect=spec.max_effect,
        onset=ILLUSTRATIVE_ONSET_H,
        tmax=ILLUSTRATIVE_TMAX_H,
        half_life=ILLUSTRATIVE_HALF_LIFE_H,
        ka="",
        ec50="",
        cmax_ref="",
        reference=ILLUSTRATIVE_REFERENCE,
        intrinsic_activity=spec.intrinsic_activity or "0.5",
        mechanism=spec.mechanism,
        dose_labels=list(ILLUSTRATIVE_DOSE_LABELS),
        dose_levels=list(ILLUSTRATIVE_DOSE_LEVELS),
        duration=ILLUSTRATIVE_DURATION_H,
    )
    payload.update(overrides)
    return payload


def study_candidate_name() -> str:
    return "Illustrative candidate"


def study_question() -> str:
    return "Which relative exposure lowers the beta proxy without raising the overdrive proxy?"


ILLUSTRATIVE_CUSTOM_YAML = """\
# Illustrative default, not measured PK. No EC50 or Cmax is asserted.
name: Illustrative custom candidate
description: User-defined mechanistic scenario. Assay and PK values are optional and off.
notes: Affinity, BBB, and toxicity are user inputs or unknown. potency not specified; normalized exposure shape
pk:
  onset_h: 0.0
  tmax_h: 1.5
  half_life_h: 2.0
  n_doses: 1
modulations:
  - target_node: d1_direct
    mechanism: agonist
    max_effect: 0.8
    hill_coeff: 1.0
"""


def simulation_config(
    baseline: str,
    duration_h: float,
    *,
    seed: int = SESSION_DEFAULTS["seed"],
    noise_std: float = SESSION_DEFAULTS["noise_std"],
    analysis_window_start_h: float = SESSION_DEFAULTS["analysis_window_start_h"],
    analysis_window_end_h: float = SESSION_DEFAULTS["analysis_window_end_h"],
    n_sensitivity_runs: int = 6,
) -> SimulationConfig:
    return SimulationConfig(
        baseline_preset=BaselinePreset(baseline),
        duration_h=float(duration_h),
        dt_h=0.05,
        noise_std=float(noise_std),
        seed=int(seed),
        n_sensitivity_runs=int(n_sensitivity_runs),
        analysis_window_start_h=float(analysis_window_start_h),
        analysis_window_end_h=float(analysis_window_end_h),
    )


def perturbation_for_illustration(spec: MechanismIllustration) -> MoleculePerturbation:
    """Reference-dose perturbation on the normalized path. No EC50 or Cmax."""
    modulation = {
        "target_node": TargetNode(spec.target_node),
        "mechanism": Mechanism(spec.mechanism),
        "max_effect": float(spec.max_effect),
        "hill_coeff": 1.0,
    }
    if spec.intrinsic_activity:
        modulation["intrinsic_activity"] = float(spec.intrinsic_activity)
    return MoleculePerturbation(
        name=spec.title,
        description=spec.blurb,
        notes=f"{ILLUSTRATIVE_PK_LABEL} {POTENCY_NOT_ASSERTED}",
        pk=PKProfile(
            onset_h=float(ILLUSTRATIVE_ONSET_H),
            tmax_h=float(ILLUSTRATIVE_TMAX_H),
            half_life_h=float(ILLUSTRATIVE_HALF_LIFE_H),
            n_doses=1,
        ),
        modulations=[TargetModulation(**modulation)],
    )


def run_illustrative_study(n_uncertainty: int = 0, illustration_id: str | None = None):
    """Run the Study design form exactly as it is prefilled."""
    spec = illustration_by_id(illustration_id) if illustration_id else default_illustration()
    parsed = parse_study_design_inputs(**study_form_kwargs(spec.id))
    if not parsed.ok:
        raise RuntimeError("Illustrative study defaults did not parse.")
    pert = perturbation_for_illustration(spec)
    cfg = simulation_config(SESSION_DEFAULTS["baseline"], parsed.duration_h, n_sensitivity_runs=max(n_uncertainty, 1))
    params = load_baseline_params(cfg.baseline_preset)
    return run_study_from_base(
        cfg,
        pert,
        params,
        labels=parsed.labels,
        levels=parsed.levels,
        reference_level=parsed.reference_level,
        use_hill=parsed.use_hill,
        question=study_question(),
        n_uncertainty=n_uncertainty,
        illustration=True,
        evidence_ids=list(spec.evidence_ids),
        potency_note=POTENCY_NOT_ASSERTED,
    )


def _read_number(
    raw: str,
    field: str,
    problems: list[str],
    *,
    required: bool,
    positive: bool,
) -> float | None:
    parsed = parse_numeric_field(raw)
    if parsed.blank:
        if required:
            problems.append(f"{field} is blank.")
        return None
    if parsed.invalid or parsed.value is None:
        problems.append(f"{field} is not a number.")
        return None
    if positive and parsed.value <= 0:
        problems.append(f"{field} must be a positive number.")
        return None
    if not positive and parsed.value < 0:
        problems.append(f"{field} must be zero or positive.")
        return None
    return parsed.value


def custom_form_perturbation(
    *,
    name: str,
    target_node: str,
    mechanism: str,
    max_effect: str,
    onset: str,
    tmax: str,
    half_life: str,
    intrinsic_activity: str,
    ec50: str,
    cmax: str,
    ka: str,
    use_assay: bool,
) -> tuple[MoleculePerturbation | None, list[str]]:
    """Build one custom candidate. Assay PK is ignored unless ``use_assay`` is true."""
    problems: list[str] = []
    if not str(name).strip():
        problems.append("Enter a candidate name.")
    effect = _read_number(max_effect, "max_effect", problems, required=True, positive=True)
    onset_h = _read_number(onset, "Onset (h)", problems, required=True, positive=False)
    tmax_h = _read_number(tmax, "Tmax (h)", problems, required=True, positive=True)
    half_life_h = _read_number(half_life, "Half-life (h)", problems, required=True, positive=True)
    ia = None
    needs_ia = mechanism == Mechanism.PARTIAL_AGONIST.value
    if needs_ia or str(intrinsic_activity).strip():
        ia = _read_number(
            intrinsic_activity,
            "Intrinsic activity",
            problems,
            required=needs_ia,
            positive=True,
        )
        if ia is not None and not (0.0 < ia <= 1.0):
            problems.append("Intrinsic activity must be greater than 0 and at most 1.")
            ia = None
    ec50_nM = cmax_nM = ka_per_h = None
    if use_assay:
        ec50_nM = _read_number(ec50, "EC50 (nM)", problems, required=True, positive=True)
        cmax_nM = _read_number(cmax, "Reference Cmax (nM)", problems, required=True, positive=True)
        ka_per_h = _read_number(ka, "ka (1/h)", problems, required=False, positive=True)
        if ec50_nM is None and not any("EC50" in item for item in problems):
            problems.append("EC50 (nM) is required when assay PK is on. It is not filled in for you.")
        if cmax_nM is None and not any("Cmax" in item for item in problems):
            problems.append("Reference Cmax (nM) is required when assay PK is on. It is not filled in for you.")
    if problems or effect is None or onset_h is None or tmax_h is None or half_life_h is None:
        return None, problems
    pk_kwargs = {
        "onset_h": onset_h,
        "tmax_h": tmax_h,
        "half_life_h": half_life_h,
        "n_doses": 1,
    }
    if ka_per_h is not None:
        pk_kwargs["ka_per_h"] = ka_per_h
    if cmax_nM is not None:
        pk_kwargs["cmax_nM"] = cmax_nM
    mod = {
        "target_node": TargetNode(target_node),
        "mechanism": Mechanism(mechanism),
        "max_effect": effect,
        "hill_coeff": 1.0,
    }
    if ec50_nM is not None:
        mod["ec50_nM"] = ec50_nM
    if ia is not None and mechanism == Mechanism.PARTIAL_AGONIST.value:
        mod["intrinsic_activity"] = ia
    notes = (
        "EC50 and Cmax are user-supplied inputs, not predictions."
        if use_assay
        else f"{ILLUSTRATIVE_PK_LABEL} {POTENCY_NOT_ASSERTED}"
    )
    return (
        MoleculePerturbation(
            name=name.strip(),
            notes=notes,
            pk=PKProfile(**pk_kwargs),
            modulations=[TargetModulation(**mod)],
        ),
        [],
    )


def run_illustrative_custom(n_uncertainty: int = 0):
    """Run the Custom candidate form exactly as it is prefilled (normalized path)."""
    from .config_loader import parse_perturbation_yaml

    pert = parse_perturbation_yaml(ILLUSTRATIVE_CUSTOM_YAML)
    cfg = simulation_config(SESSION_DEFAULTS["baseline"], 8.0, n_sensitivity_runs=max(n_uncertainty, 1))
    params = load_baseline_params(cfg.baseline_preset)
    t, states, gains = simulate(cfg, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(t, readouts, therapeutic_window=analysis_window_for(cfg), pk=pert.pk)
    return cfg, pert, t, states, readouts, kpis, gains


def run_illustration_on_baseline(spec: MechanismIllustration, baseline: str):
    """One normalized-path run of a catalog mechanism on one baseline."""
    pert = perturbation_for_illustration(spec)
    cfg = simulation_config(baseline, 8.0, n_sensitivity_runs=1)
    params = load_baseline_params(baseline)
    t, states, gains = simulate(cfg, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(t, readouts, therapeutic_window=analysis_window_for(cfg), pk=pert.pk)
    t0, s0, _ = simulate(cfg, None, params)
    base = summarize_kpis(t0, compute_readouts(t0, s0, params), therapeutic_window=analysis_window_for(cfg))
    return {
        "perturbation": pert,
        "kpis": kpis,
        "base": base,
        "gains": gains,
        "states": states,
        "readouts": readouts,
        "exposure_gains": compute_perturbation_gains(t, pert, endogenous_tone=float(params.dopamine_tone)),
    }
