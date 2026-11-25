"""Monte Carlo / LHS uncertainty and one-at-a-time sensitivity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .circuit_model import default_params_for_preset, simulate
from .readouts import analysis_window_for, compute_readouts, summarize_kpis
from .schemas import CircuitParams, MoleculePerturbation, SimulationConfig


@dataclass
class UncertaintyResult:
    kpi_names: list[str]
    mean: dict[str, float]
    p10: dict[str, float]
    p50: dict[str, float]
    p90: dict[str, float]
    samples: dict[str, np.ndarray]
    sensitivity_ranking: list[tuple[str, float]]
    time_h: Optional[np.ndarray] = None
    trajectory_p10: Optional[dict[str, np.ndarray]] = None
    trajectory_p50: Optional[dict[str, np.ndarray]] = None
    trajectory_p90: Optional[dict[str, np.ndarray]] = None


def _latin_hypercube(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    """Simple LHS in [0, 1]^d."""
    cut = np.linspace(0, 1, n + 1)
    u = rng.uniform(size=(n, d))
    pts = np.zeros((n, d))
    for j in range(d):
        pts[:, j] = cut[:-1] + u[:, j] * (cut[1] - cut[0])
        rng.shuffle(pts[:, j])
    return pts


def _perturb_params(
    base: CircuitParams,
    scale_da: float,
    scale_potency: float,
) -> CircuitParams:
    data = base.model_dump()
    data["dopamine_tone"] = float(np.clip(base.dopamine_tone * scale_da, 0.05, 2.0))
    # Potency/exposure proxy: scale key synaptic weights lightly with exposure uncertainty
    for key in ("w_d1_gpi", "w_stn_gpi", "w_ctx_d1", "w_ctx_d2"):
        data[key] = float(data[key] * scale_potency)
    return CircuitParams(**data)


def _perturb_molecule(
    mol: Optional[MoleculePerturbation],
    scale_exposure: float,
) -> Optional[MoleculePerturbation]:
    if mol is None:
        return None
    data = mol.model_dump()
    mods = []
    for m in data["modulations"]:
        m = dict(m)
        m["max_effect"] = float(m["max_effect"] * scale_exposure)
        if m.get("ec50_nM") is not None:
            m["ec50_nM"] = float(m["ec50_nM"] / max(scale_exposure, 1e-6))
        mods.append(m)
    data["modulations"] = mods
    pk = dict(data["pk"])
    if pk.get("cmax_nM") is not None:
        pk["cmax_nM"] = float(pk["cmax_nM"] * scale_exposure)
    data["pk"] = pk
    return MoleculePerturbation(**data)


def run_uncertainty(
    config: SimulationConfig,
    perturbation: Optional[MoleculePerturbation] = None,
    params: Optional[CircuitParams] = None,
    n_runs: Optional[int] = None,
    relative_spread: float = 0.15,
    method: str = "lhs",
) -> UncertaintyResult:
    """Gaussian MC or LHS over ±relative_spread on DA tone and potency/exposure."""
    n = n_runs or config.n_sensitivity_runs
    rng = np.random.default_rng(config.seed + 17)
    base = params or default_params_for_preset(config.baseline_preset)

    if method == "lhs":
        unit = _latin_hypercube(n, 2, rng)
        # Map [0,1] → [1-s, 1+s]
        scales = 1.0 + (unit - 0.5) * 2.0 * relative_spread
    else:
        scales = rng.normal(1.0, relative_spread / 1.96, size=(n, 2))
        scales = np.clip(scales, 1.0 - relative_spread, 1.0 + relative_spread)

    window = analysis_window_for(config)
    band_keys = (
        "beta_proxy",
        "motor_readiness",
        "overdrive_index",
        "pathway_balance",
        "gpi_output",
        "stn_gpe_coupling",
        "stn_burstiness_proxy",
        "direct_indirect_balance_index",
        "therapeutic_window_trace",
    )
    kpi_rows: list[dict[str, float]] = []
    band_stack: dict[str, list[np.ndarray]] = {k: [] for k in band_keys}
    t_ref: Optional[np.ndarray] = None
    for i in range(n):
        scale_da, scale_pot = float(scales[i, 0]), float(scales[i, 1])
        p = _perturb_params(base, scale_da, scale_pot)
        mol = _perturb_molecule(perturbation, scale_pot)
        cfg = config.model_copy(update={"seed": config.seed + i, "noise_std": 0.0})
        t, states, _ = simulate(cfg, mol, p)
        readouts = compute_readouts(t, states, p)
        kpi_rows.append(
            summarize_kpis(
                t,
                readouts,
                therapeutic_window=window,
                pk=mol.pk if mol is not None else None,
            )
        )
        if t_ref is None:
            t_ref = t
        for key in band_keys:
            band_stack[key].append(np.asarray(readouts[key], dtype=float))

    names = sorted(kpi_rows[0].keys())
    samples = {k: np.array([row[k] for row in kpi_rows], dtype=float) for k in names}

    def pct(arr: np.ndarray, q: float) -> float:
        if arr.size == 0 or np.all(np.isnan(arr)):
            return float("nan")
        return float(np.nanpercentile(arr, q))

    def _mean(arr: np.ndarray) -> float:
        if arr.size == 0 or np.all(np.isnan(arr)):
            return float("nan")
        return float(np.nanmean(arr))

    mean = {k: _mean(v) for k, v in samples.items()}
    p10 = {k: pct(v, 10) for k, v in samples.items()}
    p50 = {k: pct(v, 50) for k, v in samples.items()}
    p90 = {k: pct(v, 90) for k, v in samples.items()}

    def band(stack: list[np.ndarray], q: float) -> np.ndarray:
        arr = np.vstack(stack)
        return np.nanpercentile(arr, q, axis=0)

    trajectory_p10 = {k: band(v, 10) for k, v in band_stack.items()}
    trajectory_p50 = {k: band(v, 50) for k, v in band_stack.items()}
    trajectory_p90 = {k: band(v, 90) for k, v in band_stack.items()}

    ranking = one_at_a_time_sensitivity(config, perturbation, base, relative_spread)
    return UncertaintyResult(
        kpi_names=names,
        mean=mean,
        p10=p10,
        p50=p50,
        p90=p90,
        samples=samples,
        sensitivity_ranking=ranking,
        time_h=t_ref,
        trajectory_p10=trajectory_p10,
        trajectory_p50=trajectory_p50,
        trajectory_p90=trajectory_p90,
    )


def one_at_a_time_sensitivity(
    config: SimulationConfig,
    perturbation: Optional[MoleculePerturbation],
    base: CircuitParams,
    relative_spread: float = 0.15,
) -> list[tuple[str, float]]:
    """OAT sensitivity on motor_mean_window absolute change."""
    factors = {
        "dopamine_tone": "dopamine_tone",
        "w_d1_gpi": "w_d1_gpi",
        "w_stn_gpi": "w_stn_gpi",
        "w_ctx_stn": "w_ctx_stn",
        "exposure_potency": None,
    }

    window = analysis_window_for(config)

    def metric(params: CircuitParams, mol: Optional[MoleculePerturbation], seed: int) -> float:
        cfg = config.model_copy(update={"seed": seed, "noise_std": 0.0})
        t, states, _ = simulate(cfg, mol, params)
        r = compute_readouts(t, states, params)
        return summarize_kpis(
            t,
            r,
            therapeutic_window=window,
            pk=mol.pk if mol is not None else None,
        )["motor_mean_window"]

    baseline = metric(base, perturbation, config.seed)
    scores: list[tuple[str, float]] = []
    for name, field in factors.items():
        if field is None:
            mol_hi = _perturb_molecule(perturbation, 1.0 + relative_spread)
            mol_lo = _perturb_molecule(perturbation, 1.0 - relative_spread)
            hi = metric(base, mol_hi, config.seed)
            lo = metric(base, mol_lo, config.seed)
        else:
            data_hi = base.model_dump()
            data_lo = base.model_dump()
            data_hi[field] = float(data_hi[field] * (1.0 + relative_spread))
            data_lo[field] = float(data_lo[field] * (1.0 - relative_spread))
            hi = metric(CircuitParams(**data_hi), perturbation, config.seed)
            lo = metric(CircuitParams(**data_lo), perturbation, config.seed)
        scores.append((name, abs(hi - lo)))
    scores.sort(key=lambda x: x[1], reverse=True)
    # Normalize
    total = sum(s for _, s in scores) or 1.0
    return [(n, s / total) for n, s in scores]
