"""Exposure–response sweep for one mechanistic candidate.

The sweep multiplies exposure and records dimensionless proxy readouts.
It does not fit a clinical dose, invent an EC50, or claim efficacy.

Scaling
-------
* If every used modulation has an EC50 and the PK profile has ``cmax_nM``,
  each multiple scales ``cmax_nM`` (a Cmax/EC50 multiple). EC50 is unchanged.
* Otherwise the candidate is on the normalized path (potency not specified).
  Each multiple scales ``max_effect``. Intrinsic activity is not scaled.

Therapeutic-window band
    A multiple is inside the band when, relative to a drug-free run of the
    same baseline:

    * window-mean motor readiness exceeds the drug-free window mean by
      ``MOTOR_ON_DELTA``, and
    * overdrive peak is below ``OVERDRIVE_THRESHOLD``, and
    * beta suppression versus the run's own untreated baseline is positive.

Inverted-U
    Flagged for a metric when the maximum is strictly interior and exceeds
    both endpoints by ``INVERTED_U_MARGIN`` (dimensionless proxy).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .circuit_model import simulate
from .readouts import (
    MOTOR_ON_DELTA,
    OVERDRIVE_THRESHOLD,
    analysis_window_for,
    compute_readouts,
    summarize_kpis,
)
from .schemas import CircuitParams, MoleculePerturbation, SimulationConfig

INVERTED_U_MARGIN = 0.02


@dataclass
class DoseResponseResult:
    multiples: np.ndarray
    metrics: dict[str, np.ndarray]
    in_therapeutic_window: np.ndarray
    inverted_u: dict[str, bool]
    baseline_motor: float
    baseline_beta: float
    scaling_mode: str
    notes: list[str] = field(default_factory=list)


def scaling_mode(perturbation: MoleculePerturbation) -> str:
    pk = perturbation.pk
    hill = pk.cmax_nM is not None and any(m.ec50_nM is not None for m in perturbation.modulations)
    if hill:
        return "cmax_over_ec50"
    return "normalized_max_effect"


def scale_perturbation(perturbation: MoleculePerturbation, multiple: float) -> MoleculePerturbation:
    """Return a copy whose exposure is multiplied. Intrinsic activity is kept."""
    if multiple <= 0:
        raise ValueError("exposure multiple must be positive")
    data = perturbation.model_dump()
    mode = scaling_mode(perturbation)
    if mode == "cmax_over_ec50":
        data["pk"]["cmax_nM"] = float(data["pk"]["cmax_nM"]) * float(multiple)
    else:
        for mod in data["modulations"]:
            mod["max_effect"] = float(mod["max_effect"]) * float(multiple)
    data["notes"] = (
        (data.get("notes") or "")
        + f" Exposure multiple {multiple:g} ({mode})."
    ).strip()
    return MoleculePerturbation(**data)


def flag_inverted_u(
    multiples: np.ndarray,
    values: np.ndarray,
    margin: float = INVERTED_U_MARGIN,
) -> bool:
    """True when the series peaks inside the sweep, above both endpoints."""
    y = np.asarray(values, dtype=float)
    x = np.asarray(multiples, dtype=float)
    if y.size < 3 or not np.all(np.isfinite(y)):
        return False
    order = np.argsort(x)
    y = y[order]
    i_peak = int(np.argmax(y))
    if i_peak == 0 or i_peak == len(y) - 1:
        return False
    peak = float(y[i_peak])
    return (peak - float(y[0]) > margin) and (peak - float(y[-1]) > margin)


def run_dose_response(
    config: SimulationConfig,
    perturbation: MoleculePerturbation,
    params: CircuitParams,
    multiples: np.ndarray | list[float],
) -> DoseResponseResult:
    """Sweep exposure multiples for one candidate. Nominal (noise-free) only."""
    grid = np.asarray(list(multiples), dtype=float)
    if grid.size < 2:
        raise ValueError("dose-response sweep needs at least two multiples")
    if np.any(grid <= 0):
        raise ValueError("exposure multiples must be positive")

    window = analysis_window_for(config)
    mode = scaling_mode(perturbation)
    cfg = config.model_copy(update={"noise_std": 0.0})

    t_u, states_u, _ = simulate(cfg, None, params)
    base_k = summarize_kpis(
        t_u, compute_readouts(t_u, states_u, params), therapeutic_window=window
    )

    keys = (
        "beta_suppression_pct",
        "motor_mean_window",
        "overdrive_peak",
        "overdrive_mean_window",
        "therapeutic_window_score",
        "off_time_h",
        "on_time_h",
        "on_time_horizon_h",
    )
    store = {k: [] for k in keys}
    for multiple in grid:
        pert = scale_perturbation(perturbation, float(multiple))
        t, states, _ = simulate(cfg, pert, params)
        kpis = summarize_kpis(
            t, compute_readouts(t, states, params), therapeutic_window=window
        )
        for key in keys:
            store[key].append(float(kpis[key]))
    metrics = {k: np.asarray(v, dtype=float) for k, v in store.items()}

    in_band = (
        (metrics["motor_mean_window"] >= base_k["motor_mean_window"] + MOTOR_ON_DELTA)
        & (metrics["overdrive_peak"] < OVERDRIVE_THRESHOLD)
        & (metrics["beta_suppression_pct"] > 0.0)
    )
    inverted = {
        "motor_mean_window": flag_inverted_u(grid, metrics["motor_mean_window"]),
        "beta_suppression_pct": flag_inverted_u(grid, metrics["beta_suppression_pct"]),
        "overdrive_peak": flag_inverted_u(grid, metrics["overdrive_peak"]),
    }
    notes = [
        "Mechanistic exposure sweep. Not a clinical dose-finding result.",
        (
            "Scaling mode: Cmax/EC50 multiple (EC50 held fixed)."
            if mode == "cmax_over_ec50"
            else "Scaling mode: normalized max_effect multiple. potency not specified; normalized exposure shape."
        ),
        (
            "Therapeutic-window band: motor window mean at least "
            f"{MOTOR_ON_DELTA} above the drug-free window mean, overdrive peak "
            f"< {OVERDRIVE_THRESHOLD}, and positive beta suppression versus the untreated baseline."
        ),
    ]
    return DoseResponseResult(
        multiples=grid,
        metrics=metrics,
        in_therapeutic_window=np.asarray(in_band, dtype=bool),
        inverted_u=inverted,
        baseline_motor=float(base_k["motor_mean_window"]),
        baseline_beta=float(base_k["beta_mean_window"]),
        scaling_mode=mode,
        notes=notes,
    )
