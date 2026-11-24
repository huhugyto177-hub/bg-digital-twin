"""Virtual preclinical readouts from CBGT state trajectories."""

from __future__ import annotations

from typing import Optional

import numpy as np

from .schemas import CircuitParams


def compute_readouts(
    t: np.ndarray,
    states: dict[str, np.ndarray],
    params: Optional[CircuitParams] = None,
) -> dict[str, np.ndarray]:
    """Derive pathway balance, beta proxy, motor readiness, over-drive index."""
    params = params or CircuitParams()
    d1 = states["d1"]
    d2 = states["d2"]
    gpe = states["gpe"]
    stn = states["stn"]
    gpi = states["gpi"]
    thal = states["thal"]
    ctx = states["ctx"]

    # Pathway balance: direct vs indirect+hyperdirect drive onto GPi
    direct = params.w_d1_gpi * d1
    indirect = params.w_stn_gpi * stn - params.w_gpe_gpi * gpe
    pathway_balance = direct - indirect  # higher → more release of thalamus

    # LFP beta proxy (13–30 Hz abstraction): elevated when STN–GPe loop + high GPi
    # Pathological beta associated with DA depletion / STN–GPe resonance proxy
    stn_gpe_couple = stn * (1.0 + 0.5 * np.maximum(0.4 - gpe, 0.0))
    beta_proxy = 0.55 * stn_gpe_couple + 0.35 * gpi + 0.15 * (1.0 - np.clip(pathway_balance, -1, 1))
    beta_proxy = np.clip(beta_proxy, 0.0, 3.0)

    # Motor readiness / action-selection proxy: thalamic facilitation + pathway balance
    motor_readiness = np.clip(0.7 * thal + 0.4 * pathway_balance + 0.15 * ctx - 0.25 * beta_proxy, 0.0, 3.0)

    # Over-drive / dyskinesia risk: excessive direct + low GPi + high thal relative to beta
    overdrive_index = np.clip(
        0.6 * d1 + 0.35 * thal - 0.4 * gpi - 0.2 * beta_proxy + 0.25 * np.maximum(pathway_balance, 0),
        0.0,
        3.0,
    )

    return {
        "pathway_balance": pathway_balance,
        "beta_proxy": beta_proxy,
        "motor_readiness": motor_readiness,
        "overdrive_index": overdrive_index,
        "direct_pathway": direct,
        "indirect_hyperdirect": indirect,
    }


def summarize_kpis(
    t: np.ndarray,
    readouts: dict[str, np.ndarray],
    therapeutic_window: tuple[float, float] | None = None,
    overdrive_threshold: float = 1.15,
) -> dict[str, float]:
    """AUC / peak / time-in-window KPIs."""
    if therapeutic_window is None:
        # Default: around expected PK peak region for 8h sims
        therapeutic_window = (0.8, 3.5)
    t0, t1 = therapeutic_window
    mask = (t >= t0) & (t <= t1)
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)

    def auc(y: np.ndarray) -> float:
        trap = getattr(np, "trapezoid", None) or np.trapz
        return float(trap(y, t))

    def peak(y: np.ndarray) -> float:
        return float(np.max(y))

    def mean_window(y: np.ndarray) -> float:
        return float(np.mean(y[mask]))

    beta = readouts["beta_proxy"]
    motor = readouts["motor_readiness"]
    od = readouts["overdrive_index"]

    time_overdrive = float(np.sum(od[mask] > overdrive_threshold) * (t[1] - t[0] if len(t) > 1 else 0.0))
    overdrive_flag = 1.0 if float(np.max(od[mask])) >= overdrive_threshold else 0.0

    return {
        "beta_auc": auc(beta),
        "beta_peak": peak(beta),
        "beta_mean_window": mean_window(beta),
        "motor_auc": auc(motor),
        "motor_peak": peak(motor),
        "motor_mean_window": mean_window(motor),
        "overdrive_peak": peak(od),
        "overdrive_mean_window": mean_window(od),
        "overdrive_time_in_window": time_overdrive,
        "overdrive_flag": overdrive_flag,
        "pathway_balance_mean_window": mean_window(readouts["pathway_balance"]),
        "direct_mean_window": mean_window(readouts["direct_pathway"]),
        "indirect_mean_window": mean_window(readouts["indirect_hyperdirect"]),
    }
