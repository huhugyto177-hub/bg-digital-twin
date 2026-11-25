"""Pydantic v2 data models for the Basal Ganglia Digital Twin."""

from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


class TargetNode(str, Enum):
    DOPAMINE_TONE = "dopamine_tone"
    D1_DIRECT = "d1_direct"
    D2_INDIRECT = "d2_indirect"
    STN_HYPERDIRECT = "stn_hyperdirect"
    GPE = "gpe"
    GPI_SNR = "gpi_snr"


class Mechanism(str, Enum):
    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    PARTIAL_AGONIST = "partial_agonist"
    PAM = "pam"
    NAM = "nam"
    CHANNEL_BLOCKER = "channel_blocker"


class TargetModulation(BaseModel):
    """Receptor/channel/pathway modulation specification."""

    target_node: TargetNode
    mechanism: Mechanism
    max_effect: float = Field(..., description="Peak signed effect magnitude (unitless gain).")
    ec50_nM: Optional[float] = Field(None, gt=0)
    hill_coeff: float = Field(1.0, gt=0)
    intrinsic_activity: Optional[float] = Field(
        None,
        gt=0,
        le=1.0,
        description=(
            "Partial-agonist ceiling as a fraction of max_effect. "
            "User-supplied mechanistic assumption, not a measured potency."
        ),
    )

    model_config = {"extra": "forbid"}

    @model_validator(mode="after")
    def _partial_agonist_needs_intrinsic_activity(self):
        if self.mechanism == Mechanism.PARTIAL_AGONIST and self.intrinsic_activity is None:
            raise ValueError(
                "partial_agonist requires intrinsic_activity in (0, 1] "
                "(user-supplied fraction of max_effect; not a measured potency)."
            )
        return self


class PKProfile(BaseModel):
    """Pharmacokinetic exposure profile (empirical, not ADMET-predictive)."""

    onset_h: float = Field(0.0, ge=0)
    tmax_h: float = Field(..., gt=0)
    half_life_h: float = Field(..., gt=0)
    cmax_nM: Optional[float] = Field(None, gt=0)
    dosing_interval_h: Optional[float] = Field(None, gt=0)
    n_doses: int = Field(1, ge=1)
    ka_per_h: Optional[float] = Field(
        None,
        gt=0,
        description="Optional user-supplied absorption rate (1/h). Not inferred from structure.",
    )

    model_config = {"extra": "forbid"}


class MoleculePerturbation(BaseModel):
    """Named therapeutic/experimental perturbation candidate."""

    name: str
    description: str = ""
    modulations: list[TargetModulation] = Field(default_factory=list)
    pk: PKProfile
    notes: str = ""

    model_config = {"extra": "forbid"}


class BaselinePreset(str, Enum):
    HEALTHY = "healthy"
    PD_MILD = "pd_mild"
    PD_DEPLETED = "pd_depleted"
    PD_ADVANCED = "pd_advanced"
    LID_PRONE = "lid_prone"
    HD_INDIRECT_LOSS = "hd_indirect_loss"
    DYSTONIA_SURROUND = "dystonia_surround"


class SimulationConfig(BaseModel):
    """Run configuration for a deterministic simulation."""

    baseline_preset: BaselinePreset = BaselinePreset.HEALTHY
    duration_h: float = Field(8.0, gt=0)
    dt_h: float = Field(0.05, gt=0)
    noise_std: float = Field(0.0, ge=0)
    seed: int = 42
    n_sensitivity_runs: int = Field(30, ge=1)
    # Optional analysis window (model hours). None keeps the readout default.
    analysis_window_start_h: Optional[float] = Field(None, ge=0)
    analysis_window_end_h: Optional[float] = Field(None, gt=0)

    model_config = {"extra": "forbid"}

    @model_validator(mode="after")
    def _window_ordered(self):
        start = self.analysis_window_start_h
        end = self.analysis_window_end_h
        if (start is None) ^ (end is None):
            raise ValueError("Set both analysis_window_start_h and analysis_window_end_h, or neither.")
        if start is not None and end is not None and end <= start:
            raise ValueError("analysis_window_end_h must be greater than analysis_window_start_h.")
        return self


class CircuitParams(BaseModel):
    """CBGT mean-field parameters (Triple-Control inspired)."""

    # Baseline dopaminergic tone (0–1 scale; ~1 healthy, ~0.2 PD depleted)
    dopamine_tone: float = Field(1.0, ge=0, le=2.0)

    # Pathway gains
    w_ctx_d1: float = 1.0
    w_ctx_d2: float = 1.0
    w_ctx_stn: float = 0.8  # hyperdirect
    w_d1_gpi: float = 1.2  # direct: D1 inhibits GPi
    w_d2_gpe: float = 1.0
    w_gpe_stn: float = 1.0
    w_gpe_gpi: float = 0.6
    w_stn_gpe: float = 1.0
    w_stn_gpi: float = 1.2
    w_gpi_thal: float = 1.4
    w_thal_ctx: float = 0.5

    # Time constants (hours of simulation time; abstract rate units)
    tau_ctx: float = 0.15
    tau_d1: float = 0.2
    tau_d2: float = 0.2
    tau_gpe: float = 0.18
    tau_stn: float = 0.12
    tau_gpi: float = 0.15
    tau_thal: float = 0.2

    # Basal drives
    drive_ctx: float = 0.6
    drive_d1: float = 0.2
    drive_d2: float = 0.25
    drive_gpe: float = 0.35
    drive_stn: float = 0.3
    drive_gpi: float = 0.4
    drive_thal: float = 0.25

    # Triple-Control D2 inverted-U curvature
    d2_inverted_u_center: float = 0.55
    d2_inverted_u_width: float = 0.35

    # Noise / bounds
    state_floor: float = 0.0
    state_ceil: float = 2.5

    model_config = {"extra": "forbid"}


class SimResult(BaseModel):
    """Serialized simulation output summary (arrays stored separately as dicts)."""

    config: SimulationConfig
    baseline_preset: str
    perturbation_name: Optional[str] = None
    seed: int
    time_h: list[float]
    states: dict[str, list[float]]
    readouts: dict[str, list[float]]
    kpis: dict[str, float]
    assumptions: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("time_h")
    @classmethod
    def _nonempty(cls, v: list[float]) -> list[float]:
        if not v:
            raise ValueError("time_h must be non-empty")
        return v
