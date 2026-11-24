"""Cortico–Basal–Ganglia–Thalamic dynamical engine (Triple-Control inspired).

Maps dopaminergic tone and receptor/pathway modulations to circuit state
trajectories. Equations are separated from I/O; integration via solve_ivp.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from .pk_model import compute_perturbation_gains, time_vector
from .schemas import (
    BaselinePreset,
    CircuitParams,
    MoleculePerturbation,
    SimResult,
    SimulationConfig,
)

STATE_NAMES = ("ctx", "d1", "d2", "gpe", "stn", "gpi", "thal")


def _relu(x: float) -> float:
    return x if x > 0.0 else 0.0


def _clip(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else hi if x > hi else x


def d2_inverted_u(d2_activity: float, center: float, width: float) -> float:
    """Triple-Control style inverted-U: surround efficacy peaks at moderate D2."""
    w = max(width, 1e-6)
    return math_exp(-((d2_activity - center) ** 2) / (2.0 * w * w))


def math_exp(x: float) -> float:
    return float(np.exp(np.clip(x, -50.0, 50.0)))


def circuit_rhs(
    t: float,
    y: np.ndarray,
    params: CircuitParams,
    gain_fns: dict[str, callable],
) -> list[float]:
    """Mean-field rate equations for CBGT loop."""
    ctx, d1, d2, gpe, stn, gpi, thal = [float(v) for v in y]

    g_da = float(gain_fns.get("dopamine_tone", lambda _t: 0.0)(t))
    g_d1 = float(gain_fns.get("d1_direct", lambda _t: 0.0)(t))
    g_d2 = float(gain_fns.get("d2_indirect", lambda _t: 0.0)(t))
    g_stn = float(gain_fns.get("stn_hyperdirect", lambda _t: 0.0)(t))
    g_gpe = float(gain_fns.get("gpe", lambda _t: 0.0)(t))
    g_gpi = float(gain_fns.get("gpi_snr", lambda _t: 0.0)(t))

    da = _clip(params.dopamine_tone + g_da, 0.0, 2.0)
    # Direct pathway facilitated by DA / D1 agonists
    d1_gain = 0.4 + 0.6 * da + g_d1
    # g_d2 adds to D2-MSN drive (agonist > 0, antagonist < 0). Dopamine tone
    # itself still lowers this gain. GPe inhibition uses d2 * inverted-U
    # surround, so the circuit effect of a D2-node gain is not monotone.
    d2_gain = 0.7 + 0.3 * (1.0 - 0.5 * da) + g_d2
    surround = d2_inverted_u(d2, params.d2_inverted_u_center, params.d2_inverted_u_width)

    # Cortical drive + thalamic feedback
    i_ctx = params.drive_ctx + params.w_thal_ctx * thal
    # D1-MSN (direct)
    i_d1 = params.drive_d1 + d1_gain * params.w_ctx_d1 * ctx
    # D2-MSN (indirect / surround)
    i_d2 = params.drive_d2 + d2_gain * params.w_ctx_d2 * ctx
    # GPe: inhibited by D2, excited by STN, plus modulation
    i_gpe = (
        params.drive_gpe
        - params.w_d2_gpe * d2 * surround
        + params.w_stn_gpe * stn
        + g_gpe
    )
    # STN: hyperdirect from cortex, inhibited by GPe
    i_stn = (
        params.drive_stn
        + (params.w_ctx_stn + g_stn) * ctx
        - params.w_gpe_stn * gpe
    )
    # GPi/SNr: inhibited by D1 (direct), excited by STN, inhibited by GPe
    i_gpi = (
        params.drive_gpi
        - params.w_d1_gpi * d1
        + params.w_stn_gpi * stn
        - params.w_gpe_gpi * gpe
        + g_gpi
    )
    # Thalamus: inhibited by GPi
    i_thal = params.drive_thal - params.w_gpi_thal * gpi + 0.3 * ctx

    dctx = (-ctx + _relu(i_ctx)) / params.tau_ctx
    dd1 = (-d1 + _relu(i_d1)) / params.tau_d1
    dd2 = (-d2 + _relu(i_d2)) / params.tau_d2
    dgpe = (-gpe + _relu(i_gpe)) / params.tau_gpe
    dstn = (-stn + _relu(i_stn)) / params.tau_stn
    dgpi = (-gpi + _relu(i_gpi)) / params.tau_gpi
    dthal = (-thal + _relu(i_thal)) / params.tau_thal
    return [dctx, dd1, dd2, dgpe, dstn, dgpi, dthal]


def default_params_for_preset(preset: BaselinePreset) -> CircuitParams:
    if preset == BaselinePreset.PD_DEPLETED:
        return CircuitParams(
            dopamine_tone=0.22,
            drive_d1=0.12,
            drive_d2=0.32,
            drive_stn=0.42,
            drive_gpi=0.55,
            w_ctx_stn=1.05,
            w_stn_gpi=1.45,
        )
    return CircuitParams()


def load_circuit_params(path: Path | str) -> CircuitParams:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    return CircuitParams(**{k: v for k, v in data.items() if k in CircuitParams.model_fields})


def _make_gain_interpolators(
    t_grid: np.ndarray,
    gains: dict[str, np.ndarray],
) -> dict[str, callable]:
    fns: dict[str, callable] = {}
    for key, series in gains.items():
        if key.startswith("_"):
            continue
        # Closure over arrays
        tt = t_grid
        yy = series

        def _fn(t: float, _tt=tt, _yy=yy) -> float:
            return float(np.interp(t, _tt, _yy))

        fns[key] = _fn
    return fns


def _zero_gain_fns() -> dict[str, callable]:
    keys = (
        "dopamine_tone",
        "d1_direct",
        "d2_indirect",
        "stn_hyperdirect",
        "gpe",
        "gpi_snr",
    )
    return {k: (lambda _t: 0.0) for k in keys}


_SS_CACHE: dict[tuple, np.ndarray] = {}


def _params_cache_key(params: CircuitParams) -> tuple:
    return tuple(sorted((k, float(v)) for k, v in params.model_dump().items()))


def untreated_steady_state(params: CircuitParams, warmup_h: float = 25.0) -> np.ndarray:
    """Drug-free fixed point of the chosen baseline.

    A short drug-free integration supplies the guess; a root polish keeps
    t = 0 on the untreated equilibrium so the opening transient is not an
    artifact of the heuristic initial guess.
    """
    key = _params_cache_key(params)
    cached = _SS_CACHE.get(key)
    if cached is not None:
        return cached.copy()

    zero = _zero_gain_fns()

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        return circuit_rhs(_t, y, params, zero)

    guess = initial_state(params)
    sol = solve_ivp(
        rhs,
        (0.0, float(warmup_h)),
        guess,
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"Drug-free warm-up failed: {sol.message}")
    y = np.clip(sol.y[:, -1], params.state_floor, params.state_ceil)

    from scipy.optimize import root

    polished = root(lambda yy: np.asarray(rhs(0.0, yy), dtype=float), y, tol=1e-11)
    if polished.success and np.all(np.isfinite(polished.x)):
        y2 = np.clip(np.asarray(polished.x, dtype=float), params.state_floor, params.state_ceil)
        r1 = float(np.linalg.norm(rhs(0.0, y)))
        r2 = float(np.linalg.norm(rhs(0.0, y2)))
        if r2 <= r1:
            y = y2
    _SS_CACHE[key] = np.asarray(y, dtype=float)
    return _SS_CACHE[key].copy()


def initial_state(params: CircuitParams) -> np.ndarray:
    """Heuristic resting guess used only to start the drug-free fixed-point solve."""
    da = params.dopamine_tone
    return np.array(
        [
            0.55,
            0.35 + 0.25 * da,
            0.45 + 0.1 * (1.0 - da),
            0.5,
            0.4 + 0.25 * (1.0 - da),
            0.55 + 0.3 * (1.0 - da),
            0.4 + 0.2 * da,
        ],
        dtype=float,
    )


def simulate(
    config: SimulationConfig,
    perturbation: Optional[MoleculePerturbation] = None,
    params: Optional[CircuitParams] = None,
    return_arrays: bool = True,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Integrate CBGT ODEs. Returns (time, states, gains). Deterministic given seed."""
    rng = np.random.default_rng(config.seed)
    params = params or default_params_for_preset(config.baseline_preset)
    t = time_vector(config.duration_h, config.dt_h)
    gains = compute_perturbation_gains(
        t, perturbation, endogenous_tone=float(params.dopamine_tone)
    )
    gain_fns = _make_gain_interpolators(t, gains)

    y0 = untreated_steady_state(params)
    # Small deterministic jitter from seed (reproducible)
    if config.noise_std > 0:
        y0 = y0 + rng.normal(0.0, config.noise_std, size=y0.shape)

    def rhs(tt: float, yy: np.ndarray) -> list[float]:
        dydt = circuit_rhs(tt, yy, params, gain_fns)
        if config.noise_std > 0:
            # Process noise scaled small; use hash of time bin for determinism via rng stream
            pass
        return dydt

    sol = solve_ivp(
        rhs,
        (float(t[0]), float(t[-1])),
        y0,
        t_eval=t,
        method="RK45",
        rtol=1e-6,
        atol=1e-8,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    # Bound states post-hoc for physiological ceiling
    y = np.clip(sol.y, params.state_floor, params.state_ceil)
    # Optional observation noise (deterministic)
    if config.noise_std > 0:
        y = y + rng.normal(0.0, config.noise_std * 0.1, size=y.shape)
        y = np.clip(y, params.state_floor, params.state_ceil)

    states = {name: y[i] for i, name in enumerate(STATE_NAMES)}
    if return_arrays:
        return t, states, gains
    return t, states, gains


def run_simulation(
    config: SimulationConfig,
    perturbation: Optional[MoleculePerturbation] = None,
    params: Optional[CircuitParams] = None,
) -> SimResult:
    """High-level runner returning a SimResult (readouts filled by caller or here)."""
    from .readouts import analysis_window_for, compute_readouts, summarize_kpis

    t, states, gains = simulate(config, perturbation, params)
    readouts = compute_readouts(t, states, params or default_params_for_preset(config.baseline_preset))
    kpis = summarize_kpis(
        t,
        readouts,
        therapeutic_window=analysis_window_for(config),
        pk=perturbation.pk if perturbation is not None else None,
    )
    assumptions = [
        "Context of Use: mechanistic circuit-level translation only.",
        "Does NOT predict de novo binding affinity, BBB permeability, peripheral ADMET, or clinical efficacy.",
        "PK curves are user-specified empirical envelopes, not predicted from chemistry.",
        "Mean-field circuit abstraction; not a full spiking network.",
        f"Deterministic given seed={config.seed}.",
    ]
    return SimResult(
        config=config,
        baseline_preset=config.baseline_preset.value,
        perturbation_name=perturbation.name if perturbation else None,
        seed=config.seed,
        time_h=t.tolist(),
        states={k: v.tolist() for k, v in states.items()},
        readouts={k: v.tolist() for k, v in readouts.items()},
        kpis=kpis,
        assumptions=assumptions,
        metadata={"gains_peak": {k: float(np.max(np.abs(v))) for k, v in gains.items() if not k.startswith("_")}},
    )
