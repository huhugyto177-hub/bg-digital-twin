"""Preclinical study-design sweeps.

Users supply PK and dose levels. Illustration presets are labeled as such and
do not present unpublished numbers as measured PK or potency.

Normalized path (no EC50): dose levels scale ``max_effect``. The note is
``potency not specified; normalized exposure shape``.
Hill path (user-supplied EC50 and Cmax): dose levels scale Cmax and EC50 is fixed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .circuit_model import simulate
from .dose_response import flag_inverted_u
from .readouts import (
    MOTOR_ON_DELTA,
    OVERDRIVE_THRESHOLD,
    analysis_window_for,
    compute_readouts,
    summarize_kpis,
)
from .schemas import (
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    SimulationConfig,
    TargetModulation,
    TargetNode,
)
from .uncertainty import run_uncertainty

POTENCY_NORMALIZED = "potency not specified; normalized exposure shape"
ILLUSTRATION_BANNER = (
    "Mechanistic illustration preset, not a calibrated regimen. "
    + POTENCY_NORMALIZED
)


@dataclass
class DoseRow:
    label: str
    level: float
    perturbation: MoleculePerturbation
    kpis: dict[str, float]
    in_band: bool
    time_h: np.ndarray
    readouts: dict[str, np.ndarray]
    states: dict[str, np.ndarray]
    uncertainty: object | None = None


@dataclass
class StudyRun:
    question: str
    candidate_name: str
    doses: list[DoseRow]
    inverted_u: bool
    suggestion: str
    scaling_note: str
    illustration: bool
    evidence_ids: list[str] = field(default_factory=list)
    potency_note: str = ""
    config: Optional[SimulationConfig] = None


def _scale_perturbation(base: MoleculePerturbation, factor: float, cmax: float | None) -> MoleculePerturbation:
    data = base.model_dump()
    if cmax is not None and any(m.get("ec50_nM") for m in data["modulations"]):
        data["pk"]["cmax_nM"] = float(cmax)
    else:
        data["pk"]["cmax_nM"] = None
        for mod in data["modulations"]:
            mod["max_effect"] = float(mod["max_effect"]) * float(factor)
            mod["ec50_nM"] = None
    return MoleculePerturbation(**data)


def expand_dose_levels(
    base: MoleculePerturbation,
    labels: list[str],
    levels: list[float],
    reference_level: float,
    use_hill: bool,
) -> list[tuple[str, float, MoleculePerturbation]]:
    """Build one perturbation per dose. ``levels`` are relative doses or Cmax values."""
    if len(labels) != len(levels):
        raise ValueError("Each dose needs a label and a level.")
    if not 3 <= len(labels) <= 5:
        raise ValueError("A study design uses 3 to 5 doses.")
    if reference_level <= 0 or any(lv <= 0 for lv in levels):
        raise ValueError("Dose levels must be positive.")
    out = []
    for label, level in zip(labels, levels):
        factor = float(level) / float(reference_level)
        cmax = float(level) if use_hill else None
        pert = _scale_perturbation(base, factor if not use_hill else 1.0, cmax)
        pert = pert.model_copy(update={"name": f"{base.name} {label}"})
        out.append((label, float(level), pert))
    return out


def in_therapeutic_band(kpis: dict[str, float], baseline_motor: float) -> bool:
    benefit = float(kpis["motor_mean_window"]) - float(baseline_motor)
    return (
        benefit >= MOTOR_ON_DELTA
        and float(kpis["overdrive_peak"]) < OVERDRIVE_THRESHOLD
        and float(kpis["beta_power_change_pct"]) < 0.0
    )


def suggest_next_experiment(rows: list[dict]) -> str:
    """Rule-based wet-lab suggestion. Not a protocol and not an efficacy claim."""
    ordered = sorted(rows, key=lambda r: r["level"])
    band = [r["label"] for r in ordered if r["in_band"]]
    levels = np.array([r["level"] for r in ordered], dtype=float)
    ti = np.array([r["ti"] for r in ordered], dtype=float)
    inverted = flag_inverted_u(levels, ti)
    sentences = [
        "Suggestion, not a protocol and not a claim of efficacy.",
    ]
    if band:
        sentences.append(
            "The illustrated therapeutic-index band (motor benefit without crossing "
            f"the overdrive threshold, and lower beta than the untreated reference) "
            f"covers: {', '.join(band)}."
        )
        sentences.append(
            "A reasonable next wet-lab comparison is those doses in a unilateral "
            "6-OHDA lesion preparation, reading LFP beta and an LID rating. "
            "Registry entries ungerstedt_1970 and cenci_lundblad_2007 describe that "
            "class of experiment; this note does not reproduce their methods."
        )
    else:
        sentences.append(
            "No dose in this sweep sat inside the illustrated therapeutic-index band. "
            "Confirm target engagement before a lesion LFP/LID study, and do not "
            "treat the highest entered dose as preferable."
        )
    if inverted and len(ordered):
        peak_i = int(np.argmax(ti))
        sentences.append(
            f"Therapeutic index is lower at a higher exposure than at {ordered[peak_i]['label']} "
            "(inverted-U). Include higher doses only as a comparison, not as the preferred range."
        )
    sentences.append(
        "Do not copy these dimensionless gains into mg/kg. Exposure and potency "
        "remain user-supplied or unknown."
    )
    return " ".join(sentences)


def run_study_from_base(
    config: SimulationConfig,
    base: MoleculePerturbation,
    params,
    labels: list[str],
    levels: list[float],
    reference_level: float,
    use_hill: bool,
    question: str = "",
    n_uncertainty: int = 8,
    illustration: bool = False,
    evidence_ids: list[str] | None = None,
    potency_note: str = "",
) -> StudyRun:
    window = analysis_window_for(config)
    cfg = config.model_copy(update={"noise_std": 0.0})
    t0, s0, _ = simulate(cfg, None, params)
    base_k = summarize_kpis(t0, compute_readouts(t0, s0, params), therapeutic_window=window)
    doses_spec = expand_dose_levels(base, labels, levels, reference_level, use_hill)
    rows: list[DoseRow] = []
    for label, level, pert in doses_spec:
        t, states, _ = simulate(cfg, pert, params)
        readouts = compute_readouts(t, states, params)
        kpis = summarize_kpis(t, readouts, therapeutic_window=window, pk=pert.pk)
        unc = None
        if n_uncertainty and n_uncertainty > 0:
            unc = run_uncertainty(
                cfg.model_copy(update={"n_sensitivity_runs": int(n_uncertainty)}),
                pert,
                params,
                n_runs=int(n_uncertainty),
            )
        rows.append(
            DoseRow(
                label=label,
                level=level,
                perturbation=pert,
                kpis=kpis,
                in_band=in_therapeutic_band(kpis, base_k["motor_mean_window"]),
                time_h=t,
                readouts=readouts,
                states=states,
                uncertainty=unc,
            )
        )
    ordered = sorted(rows, key=lambda r: r.level)
    inverted = flag_inverted_u(
        np.array([r.level for r in ordered]),
        np.array([r.kpis["therapeutic_index_margin"] for r in ordered]),
    )
    suggestion_rows = [
        {
            "label": r.label,
            "level": r.level,
            "ti": r.kpis["therapeutic_index_margin"],
            "in_band": r.in_band,
        }
        for r in ordered
    ]
    if use_hill:
        scaling = (
            "Hill path: user-supplied EC50 is fixed and each dose sets Cmax. "
            "Those values are inputs, not predictions."
        )
    else:
        scaling = (
            "Normalized path: dose levels scale max_effect and no EC50 or Cmax is asserted. "
            + POTENCY_NORMALIZED
        )
    return StudyRun(
        question=question or f"How does exposure level change circuit proxies for {base.name}?",
        candidate_name=base.name,
        doses=ordered,
        inverted_u=inverted,
        suggestion=suggest_next_experiment(suggestion_rows),
        scaling_note=scaling,
        illustration=illustration,
        evidence_ids=list(evidence_ids or []),
        potency_note=potency_note,
        config=cfg,
    )


def parse_user_levels(
    labels: list[str],
    raw_levels: list[str],
    ec50_text: str,
    reference_text: str,
) -> tuple[bool, float, list[float], Optional[float]]:
    """Parse blank-or-number form fields. Raises ValueError with a user-facing message."""
    labels = [lb.strip() for lb in labels]
    if any(not lb for lb in labels):
        raise ValueError("Each dose needs a label.")
    try:
        levels = [float(x) for x in raw_levels]
    except ValueError as exc:
        raise ValueError("Dose levels must be numbers you supply. Fields start empty on purpose.") from exc
    ec50 = None
    if ec50_text.strip():
        ec50 = float(ec50_text)
        if ec50 <= 0:
            raise ValueError("EC50 must be positive if you supply one.")
    if not reference_text.strip():
        raise ValueError("Enter the reference dose (the level that maps to max_effect, or the reference Cmax).")
    reference = float(reference_text)
    use_hill = ec50 is not None
    return use_hill, reference, levels, ec50


def apply_user_pk(base_kwargs: dict, ec50: float | None) -> MoleculePerturbation:
    """Build the reference perturbation from user-supplied fields only."""
    return MoleculePerturbation(**base_kwargs)


PRESET_SPECS: list[dict] = [
    {
        "id": "ldopa",
        "title": "L-DOPA reference",
        "file": "06_ldopa_normalized_reference.yaml",
        "blurb": "Combined dopamine-tone, D1, and D2 illustration. Not a calibrated regimen.",
    },
    {
        "id": "d1_pam",
        "title": "D1 PAM",
        "file": "16_d1_pam_glovadalen_like.yaml",
        "blurb": "D1 positive allosteric modulator. Gain scales with baseline dopamine tone.",
    },
    {
        "id": "d1_partial",
        "title": "D1 partial agonist",
        "file": "15_d1_partial_agonist_tavapadon_like.yaml",
        "blurb": "Partial agonist. Intrinsic activity is a scenario fraction, not a published value.",
    },
    {
        "id": "d2_modulator",
        "title": "D2/D3 modulator",
        "file": "07_d2_agonist_pramipexole_like.yaml",
        "blurb": "D2-node agonist. Half-life 12.9 h is from Wright 1997; Tmax is the illustration timebase.",
    },
    {
        "id": "gpr6",
        "title": "GPR6 inverse agonist",
        "file": "17_gpr6_inverse_agonist.yaml",
        "blurb": "Reduced D2-MSN/indirect drive. Not calibrated to a published potency.",
    },
]


def preset_reference_perturbation(preset_id: str):
    """Load the illustration perturbation that corresponds to relative dose 1."""
    from .config_loader import load_benchmark
    from pathlib import Path

    spec = next(s for s in PRESET_SPECS if s["id"] == preset_id)
    path = Path(__file__).resolve().parents[1] / "configs" / "benchmarks" / spec["file"]
    config, pert, params, raw = load_benchmark(path)
    if pert is None:
        raise ValueError(f"Preset {preset_id} has no perturbation.")
    return spec, config, pert, params, raw
