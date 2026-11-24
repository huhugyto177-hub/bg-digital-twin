"""Pharmacokinetic exposure and target-occupancy curves.

Context of Use: maps user-specified PK parameters and potencies to
time-varying modulation envelopes. Does NOT predict BBB permeability,
de novo affinity from structure, or clinical ADMET.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from .schemas import Mechanism, MoleculePerturbation, TargetModulation, TargetNode


def _single_dose_concentration(
    t: np.ndarray,
    onset_h: float,
    tmax_h: float,
    half_life_h: float,
    cmax: float = 1.0,
    ka_per_h: Optional[float] = None,
) -> np.ndarray:
    """Bateman-like single-dose curve.

    If ``ka_per_h`` is omitted, the absorption rate is chosen so the peak lies
    near ``tmax_h`` after onset. A user-supplied ka is used as given and is not
    refit to tmax.
    """
    t_rel = np.maximum(t - onset_h, 0.0)
    if ka_per_h is not None:
        ka = float(ka_per_h)
    else:
        ka = math.log(2.0) / max(tmax_h, 1e-6)
    ke = math.log(2.0) / max(half_life_h, 1e-6)
    if abs(ka - ke) < 1e-9:
        ke = ke * 1.01
    # Normalize so max ≈ cmax
    raw = (ka / (ka - ke)) * (np.exp(-ke * t_rel) - np.exp(-ka * t_rel))
    raw = np.where(t_rel > 0, raw, 0.0)
    peak = float(np.max(raw)) if raw.size else 1.0
    if peak <= 0:
        return np.zeros_like(t)
    return cmax * (raw / peak)


def exposure_curve(
    t: np.ndarray,
    onset_h: float,
    tmax_h: float,
    half_life_h: float,
    cmax_nM: Optional[float] = None,
    dosing_interval_h: Optional[float] = None,
    n_doses: int = 1,
    ka_per_h: Optional[float] = None,
) -> np.ndarray:
    """Single- or multi-dose exposure. If cmax_nM is None, returns unitless 0–1 curve."""
    cmax = float(cmax_nM) if cmax_nM is not None else 1.0
    conc = np.zeros_like(t, dtype=float)
    interval = dosing_interval_h if dosing_interval_h is not None else max(tmax_h + half_life_h, 1.0)
    for i in range(max(n_doses, 1)):
        dose_onset = onset_h + i * interval
        conc = conc + _single_dose_concentration(
            t, dose_onset, tmax_h, half_life_h, cmax=cmax, ka_per_h=ka_per_h
        )
    return conc


def hill_occupancy(
    concentration: np.ndarray,
    ec50_nM: float,
    hill_coeff: float = 1.0,
    max_effect: float = 1.0,
) -> np.ndarray:
    """Hill equation occupancy / effect."""
    c = np.maximum(concentration, 0.0)
    n = max(hill_coeff, 1e-9)
    ec50 = max(ec50_nM, 1e-12)
    frac = (c**n) / (ec50**n + c**n)
    return max_effect * frac


def mechanism_sign(mechanism: Mechanism) -> float:
    if mechanism in (Mechanism.AGONIST, Mechanism.PARTIAL_AGONIST, Mechanism.PAM):
        return 1.0
    if mechanism in (Mechanism.ANTAGONIST, Mechanism.NAM, Mechanism.CHANNEL_BLOCKER):
        return -1.0
    return 1.0


def modulation_effect_curve(
    t: np.ndarray,
    modulation: TargetModulation,
    pk_conc: np.ndarray,
    cmax_nM: Optional[float],
    endogenous_tone: float = 1.0,
) -> np.ndarray:
    """Time-varying signed effect for one target modulation.

    Partial agonist: the curve is the full-agonist envelope multiplied by
    ``intrinsic_activity`` (user-supplied fraction of ``max_effect``).
    PAM: the positive envelope is multiplied by baseline ``endogenous_tone``
    clipped to [0, 1], a dimensionless proxy for dependence on endogenous
    dopamine. Agonist envelopes are not scaled by tone.
    """
    sign = mechanism_sign(modulation.mechanism)
    if modulation.ec50_nM is not None and cmax_nM is not None:
        # Use absolute concentration with Hill occupancy
        occ = hill_occupancy(
            pk_conc,
            ec50_nM=modulation.ec50_nM,
            hill_coeff=modulation.hill_coeff,
            max_effect=abs(modulation.max_effect),
        )
        effect = sign * occ
    else:
        # Normalized envelope: scale unitless exposure by max_effect.
        # Potency is not asserted on this path.
        if cmax_nM is not None and np.max(pk_conc) > 0:
            norm = pk_conc / np.max(pk_conc)
        else:
            peak = np.max(pk_conc) if pk_conc.size else 1.0
            norm = pk_conc / peak if peak > 0 else pk_conc
        effect = sign * abs(modulation.max_effect) * norm
    if modulation.mechanism == Mechanism.PARTIAL_AGONIST:
        ia = float(modulation.intrinsic_activity if modulation.intrinsic_activity is not None else 1.0)
        effect = effect * ia
    elif modulation.mechanism == Mechanism.PAM:
        tone = float(np.clip(endogenous_tone, 0.0, 1.0))
        effect = effect * tone
    return effect


def compute_perturbation_gains(
    t: np.ndarray,
    perturbation: Optional[MoleculePerturbation],
    endogenous_tone: float = 1.0,
) -> dict[str, np.ndarray]:
    """Return time series of additive gains keyed by TargetNode value.

    Absent perturbation → all-zero gains.
    """
    keys = [n.value for n in TargetNode]
    gains = {k: np.zeros_like(t, dtype=float) for k in keys}
    if perturbation is None:
        return gains

    pk = perturbation.pk
    conc = exposure_curve(
        t,
        onset_h=pk.onset_h,
        tmax_h=pk.tmax_h,
        half_life_h=pk.half_life_h,
        cmax_nM=pk.cmax_nM,
        dosing_interval_h=pk.dosing_interval_h,
        n_doses=pk.n_doses,
        ka_per_h=pk.ka_per_h,
    )
    for mod in perturbation.modulations:
        eff = modulation_effect_curve(
            t, mod, conc, pk.cmax_nM, endogenous_tone=endogenous_tone
        )
        gains[mod.target_node.value] = gains[mod.target_node.value] + eff
    gains["_exposure"] = conc
    return gains


def time_vector(duration_h: float, dt_h: float) -> np.ndarray:
    n = int(math.floor(duration_h / dt_h)) + 1
    return np.linspace(0.0, duration_h, n)
