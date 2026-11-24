"""Stability, determinism, and physiological bounds tests."""

from __future__ import annotations

import numpy as np

from src.circuit_model import default_params_for_preset, simulate
from src.config_loader import load_baseline_params
from src.schemas import BaselinePreset, SimulationConfig


def test_states_remain_bounded():
    config = SimulationConfig(
        baseline_preset=BaselinePreset.PD_DEPLETED,
        duration_h=8.0,
        dt_h=0.05,
        noise_std=0.0,
        seed=7,
    )
    params = load_baseline_params(BaselinePreset.PD_DEPLETED)
    t, states, _ = simulate(config, None, params)
    for name, series in states.items():
        assert np.all(np.isfinite(series)), name
        assert float(np.min(series)) >= params.state_floor - 1e-9
        assert float(np.max(series)) <= params.state_ceil + 1e-9


def test_determinism_identical_seeds():
    config = SimulationConfig(
        baseline_preset=BaselinePreset.HEALTHY,
        duration_h=5.0,
        dt_h=0.05,
        noise_std=0.02,
        seed=123,
    )
    params = default_params_for_preset(BaselinePreset.HEALTHY)
    t1, s1, _ = simulate(config, None, params)
    t2, s2, _ = simulate(config, None, params)
    np.testing.assert_array_equal(t1, t2)
    for k in s1:
        np.testing.assert_array_equal(s1[k], s2[k])


def test_different_seeds_can_differ_with_noise():
    base = dict(
        baseline_preset=BaselinePreset.HEALTHY,
        duration_h=5.0,
        dt_h=0.05,
        noise_std=0.05,
    )
    params = default_params_for_preset(BaselinePreset.HEALTHY)
    _, s1, _ = simulate(SimulationConfig(seed=1, **base), None, params)
    _, s2, _ = simulate(SimulationConfig(seed=2, **base), None, params)
    # With observation noise, trajectories should not be identical
    assert not np.array_equal(s1["thal"], s2["thal"])


def test_healthy_more_motor_ready_than_pd():
    from src.readouts import compute_readouts, summarize_kpis

    cfg = SimulationConfig(duration_h=6.0, dt_h=0.05, seed=42, noise_std=0.0)
    th, sh, _ = simulate(
        cfg.model_copy(update={"baseline_preset": BaselinePreset.HEALTHY}),
        None,
        load_baseline_params(BaselinePreset.HEALTHY),
    )
    tp, sp, _ = simulate(
        cfg.model_copy(update={"baseline_preset": BaselinePreset.PD_DEPLETED}),
        None,
        load_baseline_params(BaselinePreset.PD_DEPLETED),
    )
    kh = summarize_kpis(th, compute_readouts(th, sh, load_baseline_params(BaselinePreset.HEALTHY)))
    kp = summarize_kpis(tp, compute_readouts(tp, sp, load_baseline_params(BaselinePreset.PD_DEPLETED)))
    assert kh["motor_mean_window"] > kp["motor_mean_window"]
    assert kp["beta_mean_window"] > kh["beta_mean_window"]
