"""Virtual preclinical readouts from CBGT state trajectories.

All activity readouts are dimensionless proxies of a mean-field rate model.
Time readouts use the simulation clock that the PK envelope is defined on.
They are not clinical hours, clinical ON/OFF time, or measured LFP power.

Definitions
-----------
Analysis window
    Default ``(0.8, 3.5)`` model hours, overridable on ``SimulationConfig``.
    Window means, ON/OFF time, and overdrive time-in-window are computed
    only on samples inside this window. If the window misses the time base,
    the full record is used.

Untreated baseline
    Value of a readout at t = 0. Simulations are initialized at the drug-free
    steady state of the chosen baseline, so t = 0 is that untreated state.

Motor threshold
    ``motor_readiness[0] + MOTOR_ON_DELTA`` with ``MOTOR_ON_DELTA = 0.05``
    (dimensionless proxy).

ON-time proxy (``on_time_h``)
    Model hours inside the analysis window during which motor readiness is
    at or above the motor threshold.

OFF-time proxy (``off_time_h``)
    Model hours inside the analysis window during which motor readiness is
    below the motor threshold. With the discrete sample count, ON-time and
    OFF-time partition the window.

Horizon ON-time (``on_time_horizon_h``)
    Same threshold, counted over the full simulated record. Used to compare
    long exposure envelopes that outlast the default window.

Time-to-onset (``time_to_motor_onset_h``)
    First sample at which motor readiness reaches the motor threshold.
    ``NaN`` and ``motor_onset_reached = 0`` if it never does. Searched on
    the full record, not only the analysis window.

Beta suppression (``beta_suppression_pct``)
    ``100 * (beta[0] - mean_window(beta)) / beta[0]``. Percent of the
    untreated baseline beta proxy. Positive means the window mean is lower
    than the untreated baseline. Not a clinical percent change in LFP power.

Overdrive time-in-window (``overdrive_time_in_window``)
    Model hours inside the window during which ``overdrive_index`` exceeds
    ``OVERDRIVE_THRESHOLD`` (1.15, dimensionless proxy).

Pathway balance
    ``direct - indirect_hyperdirect``, where direct is ``w_d1_gpi * d1`` and
    indirect/hyperdirect is ``w_stn_gpi * stn - w_gpe_gpi * gpe``.

GPi output (``gpi_output``)
    The GPi/SNr state itself (dimensionless proxy of inhibitory output).

STN–GPe coupling (``stn_gpe_coupling``)
    ``stn * (1 + 0.5 * max(0.4 - gpe, 0))``. The same term that enters the
    beta proxy. Not a measured coherence.

Therapeutic-window score
    Window mean of ``max(motor - motor[0], 0) * (overdrive < threshold)``.
    Motor benefit above the untreated baseline on samples that are not in
    overdrive. Dimensionless proxy. Zero when benefit and overdrive coincide
    or when there is no benefit.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .schemas import CircuitParams

# Documented constants. Do not retune to force a benchmark result.
MOTOR_ON_DELTA = 0.05
OVERDRIVE_THRESHOLD = 1.15
DEFAULT_ANALYSIS_WINDOW = (0.8, 3.5)

READOUT_DOCS: dict[str, str] = {
    "pathway_balance": (
        "Direct minus indirect/hyperdirect drive onto GPi (dimensionless proxy). "
        "Higher means more thalamic release."
    ),
    "beta_proxy": (
        "Weighted STN–GPe coupling, GPi output, and low pathway balance "
        "(dimensionless proxy). Not a measured 13–30 Hz spectrum."
    ),
    "motor_readiness": (
        "Thalamic facilitation plus pathway balance, reduced by beta "
        "(dimensionless proxy). Not a clinical motor score."
    ),
    "overdrive_index": (
        "Excess direct-pathway and thalamic drive with low GPi "
        "(dimensionless proxy). Not a dyskinesia rating scale."
    ),
    "direct_pathway": "w_d1_gpi * D1-MSN state (dimensionless proxy).",
    "indirect_hyperdirect": (
        "w_stn_gpi * STN − w_gpe_gpi * GPe (dimensionless proxy)."
    ),
    "gpi_output": "GPi/SNr state (dimensionless proxy of inhibitory output).",
    "stn_gpe_coupling": (
        "stn * (1 + 0.5 * max(0.4 − gpe, 0)) (dimensionless proxy). "
        "Same coupling term used inside the beta proxy."
    ),
    "therapeutic_window_trace": (
        "max(motor − motor(t=0), 0) on samples below the overdrive threshold; "
        "otherwise 0 (dimensionless proxy)."
    ),
    "beta_mean_window": "Mean beta proxy inside the analysis window (dimensionless proxy).",
    "motor_mean_window": "Mean motor readiness inside the analysis window (dimensionless proxy).",
    "on_time_h": (
        "Model hours inside the analysis window at or above the motor threshold "
        "(untreated motor + 0.05). ON-time proxy, not clinical ON time."
    ),
    "off_time_h": (
        "Model hours inside the analysis window below the motor threshold. "
        "OFF-time proxy, not clinical OFF time."
    ),
    "on_time_horizon_h": (
        "Model hours at or above the motor threshold over the full simulation."
    ),
    "time_to_motor_onset_h": (
        "First time motor readiness reaches the motor threshold. "
        "n/a if it never does."
    ),
    "beta_suppression_pct": (
        "Percent reduction of window-mean beta versus the untreated (t=0) beta proxy."
    ),
    "beta_power_change_pct": (
        "Baseline-normalized percent change in the beta proxy versus the untreated "
        "reference (t=0 of this baseline, which is the untreated disease steady state "
        "when the preset is pd_depleted). Negative means lower beta than that reference. "
        "Dimensionless model proxy, not measured beta power."
    ),
    "therapeutic_index_margin": (
        "Motor-readiness benefit above the untreated reference minus any overdrive "
        "peak above the 1.15 threshold. Dimensionless proxy. Falls when a higher "
        "exposure buys less net benefit."
    ),
    "motor_at_tmax": "Motor readiness nearest the user-supplied Tmax (dimensionless proxy).",
    "motor_at_trough": "Motor readiness at the end of the dosing interval, or at the end of the record (dimensionless proxy).",
    "beta_at_tmax": "Beta proxy nearest Tmax (dimensionless proxy).",
    "beta_at_trough": "Beta proxy at the trough sample (dimensionless proxy).",
    "overdrive_at_tmax": "Overdrive index nearest Tmax (dimensionless proxy).",
    "overdrive_at_trough": "Overdrive index at the trough sample (dimensionless proxy).",
    "overdrive_time_in_window": (
        "Model hours inside the window with overdrive_index ≥ 1.15."
    ),
    "therapeutic_window_score": (
        "Window mean of motor benefit on samples that are not in overdrive "
        "(dimensionless proxy)."
    ),
}


def analysis_window_for(config) -> tuple[float, float]:
    """Return the configured analysis window, or the module default."""
    start = getattr(config, "analysis_window_start_h", None)
    end = getattr(config, "analysis_window_end_h", None)
    if start is None or end is None:
        return DEFAULT_ANALYSIS_WINDOW
    return (float(start), float(end))


def compute_readouts(
    t: np.ndarray,
    states: dict[str, np.ndarray],
    params: Optional[CircuitParams] = None,
) -> dict[str, np.ndarray]:
    """Derive pathway, beta, motor, overdrive, and node traces."""
    params = params or CircuitParams()
    d1 = states["d1"]
    d2 = states["d2"]
    gpe = states["gpe"]
    stn = states["stn"]
    gpi = states["gpi"]
    thal = states["thal"]
    ctx = states["ctx"]

    direct = params.w_d1_gpi * d1
    indirect = params.w_stn_gpi * stn - params.w_gpe_gpi * gpe
    pathway_balance = direct - indirect

    stn_gpe_coupling = stn * (1.0 + 0.5 * np.maximum(0.4 - gpe, 0.0))
    beta_proxy = 0.55 * stn_gpe_coupling + 0.35 * gpi + 0.15 * (1.0 - np.clip(pathway_balance, -1, 1))
    beta_proxy = np.clip(beta_proxy, 0.0, 3.0)

    motor_readiness = np.clip(
        0.7 * thal + 0.4 * pathway_balance + 0.15 * ctx - 0.25 * beta_proxy,
        0.0,
        3.0,
    )
    overdrive_index = np.clip(
        0.6 * d1 + 0.35 * thal - 0.4 * gpi - 0.2 * beta_proxy + 0.25 * np.maximum(pathway_balance, 0),
        0.0,
        3.0,
    )
    benefit = np.maximum(motor_readiness - float(motor_readiness[0]), 0.0)
    therapeutic_window_trace = benefit * (overdrive_index < OVERDRIVE_THRESHOLD)

    return {
        "pathway_balance": pathway_balance,
        "beta_proxy": beta_proxy,
        "motor_readiness": motor_readiness,
        "overdrive_index": overdrive_index,
        "direct_pathway": direct,
        "indirect_hyperdirect": indirect,
        "gpi_output": np.asarray(gpi, dtype=float),
        "stn_gpe_coupling": np.asarray(stn_gpe_coupling, dtype=float),
        "therapeutic_window_trace": np.asarray(therapeutic_window_trace, dtype=float),
    }


def _trap(y: np.ndarray, t: np.ndarray) -> float:
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(y, t))


def _dt(t: np.ndarray) -> float:
    if len(t) < 2:
        return 0.0
    return float(t[1] - t[0])


def _nearest(t: np.ndarray, y: np.ndarray, t_query: float) -> float:
    return float(y[int(np.argmin(np.abs(t - t_query)))])


def exposure_sample_times(pk, duration_h: float) -> tuple[float, float]:
    """Clock times for the Tmax sample and the end-of-interval trough sample.

    Tmax is ``onset + tmax`` on the simulation clock (the Bateman peak is near
    tmax after onset). The trough is the end of the dosing interval when one
    is set, otherwise the end of the record. Both are clipped to the record.
    """
    peak_t = min(float(pk.onset_h) + float(pk.tmax_h), float(duration_h))
    if pk.dosing_interval_h:
        trough_t = min(float(pk.onset_h) + float(pk.dosing_interval_h), float(duration_h))
    else:
        trough_t = float(duration_h)
    return peak_t, max(trough_t, 0.0)


def summarize_kpis(
    t: np.ndarray,
    readouts: dict[str, np.ndarray],
    therapeutic_window: tuple[float, float] | None = None,
    overdrive_threshold: float = OVERDRIVE_THRESHOLD,
    motor_on_delta: float = MOTOR_ON_DELTA,
    pk=None,
    reference_beta: float | None = None,
    reference_motor: float | None = None,
) -> dict[str, float]:
    """AUC / peak / window KPIs, including ON/OFF and therapeutic-window score."""
    if therapeutic_window is None:
        therapeutic_window = DEFAULT_ANALYSIS_WINDOW
    t0, t1 = therapeutic_window
    mask = (t >= t0) & (t <= t1)
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)

    def mean_window(y: np.ndarray) -> float:
        return float(np.mean(y[mask]))

    beta = readouts["beta_proxy"]
    motor = readouts["motor_readiness"]
    od = readouts["overdrive_index"]
    dt = _dt(t)

    time_overdrive = float(np.sum(od[mask] > overdrive_threshold) * dt)
    overdrive_flag = 1.0 if float(np.max(od[mask])) >= overdrive_threshold else 0.0

    motor_threshold = float(motor[0]) + float(motor_on_delta)
    above = motor >= motor_threshold
    on_time = float(np.sum(above[mask]) * dt)
    off_time = float(np.sum(~above[mask]) * dt)
    on_time_horizon = float(np.sum(above) * dt)
    if np.any(above):
        time_to_onset = float(t[int(np.argmax(above))])
        onset_reached = 1.0
    else:
        time_to_onset = float("nan")
        onset_reached = 0.0

    beta0 = float(beta[0])
    beta_mean = mean_window(beta)
    if abs(beta0) > 1e-9:
        beta_suppression = 100.0 * (beta0 - beta_mean) / beta0
    else:
        beta_suppression = 0.0

    tw = readouts.get("therapeutic_window_trace")
    if tw is None:
        benefit_trace = np.maximum(motor - float(motor[0]), 0.0)
        tw = benefit_trace * (od < overdrive_threshold)

    beta_ref = float(beta[0] if reference_beta is None else reference_beta)
    motor_ref = float(motor[0] if reference_motor is None else reference_motor)
    if abs(beta_ref) > 1e-9:
        beta_power_change = 100.0 * (beta_mean - beta_ref) / beta_ref
    else:
        beta_power_change = 0.0
    motor_benefit = mean_window(motor) - motor_ref
    overdrive_excess = max(0.0, float(np.max(od)) - float(overdrive_threshold))
    # Benefit above the untreated reference, reduced by any overdrive past threshold.
    therapeutic_index_margin = float(motor_benefit - overdrive_excess)

    if pk is not None and len(t):
        peak_t, trough_t = exposure_sample_times(pk, float(t[-1]))
        sampled = {
            "beta_at_tmax": _nearest(t, beta, peak_t),
            "beta_at_trough": _nearest(t, beta, trough_t),
            "motor_at_tmax": _nearest(t, motor, peak_t),
            "motor_at_trough": _nearest(t, motor, trough_t),
            "overdrive_at_tmax": _nearest(t, od, peak_t),
            "overdrive_at_trough": _nearest(t, od, trough_t),
            "sample_tmax_h": float(peak_t),
            "sample_trough_h": float(trough_t),
        }
    else:
        sampled = {
            "beta_at_tmax": float("nan"),
            "beta_at_trough": float("nan"),
            "motor_at_tmax": float("nan"),
            "motor_at_trough": float("nan"),
            "overdrive_at_tmax": float("nan"),
            "overdrive_at_trough": float("nan"),
            "sample_tmax_h": float("nan"),
            "sample_trough_h": float("nan"),
        }

    return {
        "beta_auc": _trap(beta, t),
        "beta_peak": float(np.max(beta)),
        "beta_mean_window": beta_mean,
        "motor_auc": _trap(motor, t),
        "motor_peak": float(np.max(motor)),
        "motor_mean_window": mean_window(motor),
        "overdrive_peak": float(np.max(od)),
        "overdrive_mean_window": mean_window(od),
        "overdrive_time_in_window": time_overdrive,
        "overdrive_flag": overdrive_flag,
        "pathway_balance_mean_window": mean_window(readouts["pathway_balance"]),
        "direct_mean_window": mean_window(readouts["direct_pathway"]),
        "indirect_mean_window": mean_window(readouts["indirect_hyperdirect"]),
        "gpi_mean_window": mean_window(readouts["gpi_output"]) if "gpi_output" in readouts else mean_window(np.zeros_like(motor)),
        "stn_gpe_mean_window": (
            mean_window(readouts["stn_gpe_coupling"])
            if "stn_gpe_coupling" in readouts
            else float("nan")
        ),
        "time_to_motor_onset_h": time_to_onset,
        "motor_onset_reached": onset_reached,
        "on_time_h": on_time,
        "off_time_h": off_time,
        "on_time_horizon_h": on_time_horizon,
        "beta_suppression_pct": float(beta_suppression),
        "beta_power_change_pct": float(beta_power_change),
        "therapeutic_window_score": mean_window(tw),
        "therapeutic_index_margin": therapeutic_index_margin,
        **sampled,
    }
