"""Unit tests for PK curve bounds & timing."""

from __future__ import annotations

import numpy as np

from src.pk_model import exposure_curve, hill_occupancy, compute_perturbation_gains, time_vector
from src.schemas import (
    Mechanism,
    MoleculePerturbation,
    PKProfile,
    TargetModulation,
    TargetNode,
)


def test_exposure_nonnegative_and_peaks_near_tmax():
    t = time_vector(8.0, 0.05)
    conc = exposure_curve(t, onset_h=0.5, tmax_h=1.0, half_life_h=2.0, cmax_nM=100.0)
    assert np.all(conc >= -1e-12)
    # Peak after onset, near onset+tmax
    t_peak = t[int(np.argmax(conc))]
    assert 0.5 <= t_peak <= 2.5
    assert abs(float(np.max(conc)) - 100.0) < 5.0


def test_multi_dose_higher_auc_than_single():
    t = time_vector(12.0, 0.05)
    single = exposure_curve(t, 0.0, 1.0, 2.0, cmax_nM=50.0, n_doses=1)
    multi = exposure_curve(
        t, 0.0, 1.0, 2.0, cmax_nM=50.0, dosing_interval_h=4.0, n_doses=3
    )
    trap = getattr(np, "trapezoid", None) or np.trapz
    assert float(trap(multi, t)) > float(trap(single, t))


def test_hill_occupancy_bounds():
    c = np.array([0.0, 10.0, 100.0, 1e6])
    occ = hill_occupancy(c, ec50_nM=100.0, hill_coeff=1.0, max_effect=0.8)
    assert occ[0] == 0.0
    assert 0.0 < occ[1] < 0.4
    assert abs(occ[2] - 0.4) < 0.05  # half max at EC50
    assert occ[-1] <= 0.8 + 1e-9


def test_perturbation_gains_keys_and_sign():
    t = time_vector(6.0, 0.1)
    mol = MoleculePerturbation(
        name="test",
        description="",
        pk=PKProfile(onset_h=0.0, tmax_h=1.0, half_life_h=2.0, cmax_nM=200.0),
        modulations=[
            TargetModulation(
                target_node=TargetNode.D1_DIRECT,
                mechanism=Mechanism.AGONIST,
                max_effect=0.5,
                ec50_nM=100.0,
            ),
            TargetModulation(
                target_node=TargetNode.D2_INDIRECT,
                mechanism=Mechanism.ANTAGONIST,
                max_effect=0.4,
                ec50_nM=100.0,
            ),
        ],
    )
    gains = compute_perturbation_gains(t, mol)
    assert "d1_direct" in gains and "d2_indirect" in gains
    assert float(np.max(gains["d1_direct"])) > 0
    assert float(np.min(gains["d2_indirect"])) < 0
