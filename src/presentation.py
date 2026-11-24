"""Tables, CSV exports, and the shared matplotlib style.

Displayed numbers use three significant figures. CSV files keep full precision.
"""

from __future__ import annotations

import io
import math
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

from .readouts import DEFAULT_ANALYSIS_WINDOW, READOUT_DOCS, analysis_window_for

# Restrained scientific palette. Shared by the app and the decision-brief figures.
_COLORS = ["#1B4F72", "#B85C38", "#1F7A6B", "#5C4A72", "#8A6A3B", "#3D5A80"]

PRIMARY_KPIS: list[tuple[str, str]] = [
    ("motor_mean_window", "dimensionless proxy"),
    ("beta_mean_window", "dimensionless proxy"),
    ("beta_suppression_pct", "percent of untreated baseline beta proxy"),
    ("on_time_h", "model hours inside the analysis window"),
    ("off_time_h", "model hours inside the analysis window"),
    ("on_time_horizon_h", "model hours over the full simulation"),
    ("time_to_motor_onset_h", "model hours"),
    ("overdrive_time_in_window", "model hours inside the analysis window"),
    ("overdrive_peak", "dimensionless proxy"),
    ("therapeutic_window_score", "dimensionless proxy"),
    ("therapeutic_index_margin", "dimensionless proxy"),
    ("beta_power_change_pct", "percent change vs untreated reference"),
    ("motor_at_tmax", "dimensionless proxy at Tmax"),
    ("motor_at_trough", "dimensionless proxy at trough"),
    ("beta_at_tmax", "dimensionless proxy at Tmax"),
    ("beta_at_trough", "dimensionless proxy at trough"),
    ("pathway_balance_mean_window", "dimensionless proxy"),
    ("gpi_mean_window", "dimensionless proxy"),
    ("stn_gpe_mean_window", "dimensionless proxy"),
]


def apply_mpl_style() -> None:
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.28,
            "grid.linewidth": 0.6,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#2C3E50",
            "axes.labelcolor": "#1C2833",
            "xtick.color": "#1C2833",
            "ytick.color": "#1C2833",
            "axes.prop_cycle": cycler(color=_COLORS),
            "lines.linewidth": 1.8,
        }
    )


def fmt_sig(value: float | None, sig: int = 3) -> str:
    """Format a number to ``sig`` significant figures. Non-finite values are n/a."""
    if value is None:
        return "n/a"
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "n/a"
    if math.isnan(x) or math.isinf(x):
        return "n/a"
    return f"{x:.{sig}g}"


def window_label(config=None) -> str:
    if config is None:
        start, end = DEFAULT_ANALYSIS_WINDOW
    else:
        start, end = analysis_window_for(config)
    return f"{fmt_sig(start)}–{fmt_sig(end)} h"


def kpi_comparison_frame(results: Iterable[dict]) -> pd.DataFrame:
    """Long table of mean / P10 / P50 / P90 for the primary KPIs."""
    rows = []
    for res in results:
        unc = res["unc"]
        window = window_label(res["config"])
        for key, unit in PRIMARY_KPIS:
            if key not in unc.mean:
                continue
            rows.append(
                {
                    "scenario": res["name"],
                    "metric": key,
                    "definition": READOUT_DOCS.get(key, ""),
                    "mean": unc.mean[key],
                    "p10": unc.p10[key],
                    "p50": unc.p50[key],
                    "p90": unc.p90[key],
                    "unit": unit,
                    "window": window,
                }
            )
    return pd.DataFrame(rows)


def format_kpi_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for col in ("mean", "p10", "p50", "p90"):
        if col in out.columns:
            out[col] = out[col].map(lambda v: fmt_sig(v))
    return out


def trajectories_frame(results: Iterable[dict]) -> pd.DataFrame:
    """Long CSV-ready frame of time, states, and readouts."""
    frames = []
    for res in results:
        t = np.asarray(res["t"], dtype=float)
        data = {"scenario": res["name"], "time_h": t}
        for name, series in res["states"].items():
            data[f"state_{name}"] = np.asarray(series, dtype=float)
        for name, series in res["readouts"].items():
            data[name] = np.asarray(series, dtype=float)
        frames.append(pd.DataFrame(data))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def frame_to_csv(frame: pd.DataFrame) -> str:
    buf = io.StringIO()
    frame.to_csv(buf, index=False)
    return buf.getvalue()


def sensitivity_figure(ranking: list[tuple[str, float]], title: str):
    """Horizontal bar chart of normalized OAT scores."""
    apply_mpl_style()
    labels = [name.replace("_", " ") for name, _ in ranking][::-1]
    scores = [score for _, score in ranking][::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.barh(labels, scores, color="#1B4F72")
    ax.set_xlabel("Share of |Δ motor window| (dimensionless proxy)")
    ax.set_title(title)
    fig.tight_layout()
    return fig
