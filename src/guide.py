"""Getting-started copy for the Streamlit app and for tests.

The five non-claims are included verbatim from the evidence registry so the
first screen states the same exclusions as the decision brief.
"""

from __future__ import annotations

from .catalog import (
    BASELINES,
    ILLUSTRATIVE_PK_LABEL,
    MECHANISM_ILLUSTRATIONS,
    NEW_READOUT_KEYS,
)
from .evidence_registry import cou_and_limitations_markdown
from .readouts import READOUT_DOCS


def getting_started_markdown() -> str:
    baseline_lines = "\n".join(
        f"- **{item.title}** (`{item.preset}`). {item.summary}" for item in BASELINES
    )
    mechanism_lines = "\n".join(
        f"- **{item.title}.** Node `{item.target_node}`, circuit field `{item.mechanism}`, "
        f"illustrative max_effect {item.max_effect}. {item.blurb}"
        for item in MECHANISM_ILLUSTRATIONS
    )
    readout_lines = "\n".join(
        f"- **{key}.** {READOUT_DOCS[key]}" for key in NEW_READOUT_KEYS
    )
    core_lines = "\n".join(
        f"- **{key}.** {READOUT_DOCS[key]}"
        for key in (
            "beta_proxy",
            "beta_power_change_pct",
            "motor_readiness",
            "on_time_h",
            "off_time_h",
            "overdrive_index",
            "therapeutic_index_margin",
        )
    )
    return f"""## What this tool does

Basal Ganglia Digital Twin maps a **stated** circuit modulation and a **stated** exposure curve onto a cortico–basal ganglia–thalamic rate model. It returns dimensionless proxies: beta, motor readiness, GPi/SNr output, pathway balance, and an overdrive-based dyskinesia-risk fraction.

It does not design molecules. Affinity, blood–brain barrier penetration, ADMET, and clinical efficacy stay outside the model.

## Three-step quick start

1. Stay on **Study design**. The form is already filled with an illustrative D1 direct-pathway agonist (max_effect 0.8), onset 0 h, Tmax 1.5 h, half-life 2 h, doses low/mid/high at 0.5/1/2, reference dose 1, duration 8 h, and the moderate PD baseline.
2. Leave **Supply assay PK** off. That keeps EC50, reference Cmax, and ka empty. The run uses the normalized exposure shape and does not invent a potency.
3. Click **Run entered study**. Read the cards and the KPI table for the reference dose (mid). Then change one thing, such as the baseline or the mechanism illustration, and run again.

## What each input means

| Input | Unit | How to read the default |
|---|---|---|
| Target node | circuit node | Where the gain is applied. |
| Mechanism | model field | Sets the sign of the gain. A catalog blurb states when the pharmacological class and the circuit sign differ. |
| max_effect | unitless gain at the reference dose | Illustrative magnitude. Not an EC50. |
| Intrinsic activity | fraction of max_effect, (0, 1] | Used for a partial agonist. 0.5 is a scenario fraction, not a published value. |
| Onset | hours on the model clock | {ILLUSTRATIVE_PK_LABEL} Default 0. |
| Tmax | hours after onset | {ILLUSTRATIVE_PK_LABEL} Default 1.5. |
| Half-life | hours | {ILLUSTRATIVE_PK_LABEL} Default 2. Not a published half-life for any named drug. |
| Reference dose | user scale | The level that maps to max_effect. Default 1. |
| Dose label and level | user scale | Three doses are filled (0.5, 1, 2). Two extra rows are optional and blank. |
| Duration | hours | Scenario length. Default 8. Not a clinical treatment interval. |
| Baseline | circuit state | Moderate PD depletion is the prefilled comparison. |
| Random seed | integer | Default 42. Same seed, same trajectory. |
| Observation noise | dimensionless SD | Default 0. |
| Analysis window | model hours | Default 0.8–3.5. Window means and ON/OFF time use this interval. |
| EC50 | nM | Off until you enable assay PK and type a value you have. |
| Reference Cmax | nM | Off until you enable assay PK. Required only with an EC50 you supply. |
| ka | 1/h | Off until you enable assay PK. Blank derives the absorption rate from Tmax. |

## How to read the main proxies

{core_lines}

### Added readouts

{readout_lines}

Time is the simulation clock shared with the exposure curve. It is not clinical clock time. A negative **beta power change** means the window-mean beta proxy is lower than the untreated reference of the selected baseline.

## Worked example

Leave every Study design field at its illustrative default and click **Run entered study**. The candidate is an illustrative D1 agonist on the moderate PD baseline. The mid dose is the reference (level 1) and uses max_effect 0.8. Low and high scale that gain to 0.5× and 2× because assay PK is off.

On the result, compare across low, mid, and high:

- **Beta power change** (percent versus the untreated reference). Negative means lower beta than that reference.
- **Motor readiness** (dimensionless) inside 0.8–3.5 h.
- **Dyskinesia-risk proxy** (fraction of that window at or above the overdrive threshold of 1.15).
- **GPi/SNr output rate** and the **direct/indirect balance index**.

The therapeutic-index margin falls when overdrive past 1.15 eats the motor benefit. That pattern is a model flag, not a dose recommendation. Do not copy the unitless gains into mg/kg.

## Baselines

{baseline_lines}

Benchmarks 01–17 stay on `healthy` or `pd_depleted`. New baselines do not renumber those files.

## Mechanism illustrations

Each row below is a circuit mapping with a citation in the evidence registry. max_effect values are illustrative. None of them carries an EC50, Cmax, or a published half-life.

{mechanism_lines}

## Context of Use and non-claims

{cou_and_limitations_markdown()}
"""
