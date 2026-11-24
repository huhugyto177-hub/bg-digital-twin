"""Verification that all 5 benchmark scenarios reproduce expected directional outcomes."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from src.circuit_model import simulate
from src.config_loader import load_benchmark
from src.readouts import compute_readouts, summarize_kpis

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "configs" / "benchmarks"


def _run(path: Path):
    config, pert, params, raw = load_benchmark(path)
    t, states, gains = simulate(config, pert, params)
    readouts = compute_readouts(t, states, params)
    kpis = summarize_kpis(t, readouts)
    return config, pert, params, raw, t, states, readouts, kpis, gains


def test_01_and_02_untreated_pd_higher_beta_lower_motor_vs_healthy():
    _, _, _, _, _, _, _, k_h, _ = _run(BENCH / "01_healthy_control.yaml")
    _, _, _, _, _, _, _, k_pd, _ = _run(BENCH / "02_untreated_pd.yaml")
    assert k_pd["beta_mean_window"] > k_h["beta_mean_window"]
    assert k_pd["motor_mean_window"] < k_h["motor_mean_window"]


def test_03_ldopa_restores_motor_and_suppresses_beta():
    _, _, _, _, _, _, _, k_pd, _ = _run(BENCH / "02_untreated_pd.yaml")
    _, _, _, _, _, _, _, k_ld, _ = _run(BENCH / "03_ldopa_standard.yaml")
    # Therapeutic window at peak: motor up, beta down vs untreated PD
    assert k_ld["motor_mean_window"] > k_pd["motor_mean_window"]
    assert k_ld["beta_mean_window"] < k_pd["beta_mean_window"]


def test_04_selective_d1_distinct_pathway_vs_ldopa():
    _, _, _, _, _, _, r_ld, k_ld, _ = _run(BENCH / "03_ldopa_standard.yaml")
    _, _, _, _, _, _, r_d1, k_d1, _ = _run(BENCH / "04_selective_d1_agonist.yaml")
    # Distinct pathway dynamics: D1-selective should emphasize direct relative to L-DOPA's broader DA
    # Compare ratio of direct to absolute indirect engagement in therapeutic window
    t = np.asarray(_run(BENCH / "03_ldopa_standard.yaml")[4])
    # Use KPI pathway components
    assert abs(k_d1["direct_mean_window"] - k_ld["direct_mean_window"]) > 1e-3 or abs(
        k_d1["pathway_balance_mean_window"] - k_ld["pathway_balance_mean_window"]
    ) > 1e-3
    # Trajectories should not be identical
    assert not np.allclose(r_d1["direct_pathway"], r_ld["direct_pathway"])
    assert not np.allclose(r_d1["pathway_balance"], r_ld["pathway_balance"])


def test_05_overdrive_triggers_flag():
    _, _, _, _, _, _, _, k_od, _ = _run(BENCH / "05_overdrive_dyskinesia_case.yaml")
    assert k_od["overdrive_flag"] == 1.0
    assert k_od["overdrive_peak"] >= 1.15


def test_identical_seeds_bit_for_bit():
    a = _run(BENCH / "04_selective_d1_agonist.yaml")
    b = _run(BENCH / "04_selective_d1_agonist.yaml")
    np.testing.assert_array_equal(a[4], b[4])  # time
    for k in a[5]:
        np.testing.assert_array_equal(a[5][k], b[5][k])
    for k in a[6]:
        np.testing.assert_array_equal(a[6][k], b[6][k])
