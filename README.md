# Basal Ganglia Digital Twin (v1)

Preclinical therapeutic & mechanism simulator for cortico–basal–ganglia–thalamic (CBGT) circuit translation. Grounded in classical and contemporary basal ganglia / PD / dopaminergic / LFP-beta literature (Albin–DeLong–Wichmann, Gerfen/Surmeier, STN–GPe beta work, L-DOPA PK/PD, QSP framing). Triple-Control (Li & Jin) is **one** supporting source among many — not the sole foundation.

## Context of Use

Mechanistic circuit-level translation: map **defined** receptor/channel/pathway modulations and **user-specified** PK exposure curves to predicted basal ganglia circuit activity, proxy LFP beta-band shifts, and motor-readiness / over-drive readouts.

The twin is calibrated to *directional* consensus from the evidence registry. It does **not** invent binding affinity, BBB permeability, toxicity, or clinical efficacy for novel chemical structures.

## Non-claims (strict)

This system does **NOT**:

- predict de novo chemical binding affinity from SMILES/structures
- predict blood–brain barrier permeability
- predict peripheral ADMET/toxicity
- claim clinical efficacy without experimental validation
- invent EC50 / Cmax / half-life for newly discovered molecules (user must supply mechanistic + PK assumptions)

## Setup

```bash
cd bg_digital_twin
python -m pip install -r requirements.txt
```

## Quickstart

### Run tests

```bash
python -m pytest -v
```

### Generate a client decision brief

Run from `bg_digital_twin/` (so the `src` package is importable):

```bash
python -m src.report_generator --config configs/benchmarks/04_selective_d1_agonist.yaml
```

Outputs Markdown + high-DPI PNGs + config JSON under `outputs/<scenario>/`. Every brief includes **Context of Use**, assumptions, and a full **References / Literature Basis** table from the evidence registry.

### Launch Streamlit UI

```bash
python -m streamlit run app/streamlit_app.py
```

## Deploy on Streamlit Community Cloud

This repository root **is** the app root (`bg_digital_twin`). Config and data paths are resolved from the package via `Path(__file__)`, so they work on Streamlit Cloud without absolute local paths.

1. Push this repo to GitHub (public is simplest for free Community Cloud).
2. Open [share.streamlit.io](https://share.streamlit.io) → **New app**.
3. Connect the GitHub repo and set:

| Setting | Value |
|---------|--------|
| Branch | `main` |
| Main file path | `app/streamlit_app.py` |
| Root directory | *(leave blank / repo root)* |

4. Deploy. No secrets are required for the default app.

Local entry (same as Cloud): from this directory, `streamlit run app/streamlit_app.py`.

## Literature tracking (evidence registry)

| File | Role |
|------|------|
| `configs/evidence_registry.yaml` | Authoritative citation list + CoU / non-claims / calibration note |
| `src/evidence_registry.py` | Loader, defaults, `references_markdown()`, table rows for UI/briefs |

Decision Briefs (`src/report_generator.py`) and the Streamlit **Assumptions / Citations** tab always render the full References table from this registry.

### How to add a citation when calibrating

1. Edit `configs/evidence_registry.yaml` (preferred) — add a block under `citations:`:

```yaml
  - id: my_new_citation_yyyy
    parameter_keys: [dopamine_tone, beta_proxy]   # twin params this supports
    citation: "Author et al. Title. Journal. Year;..."
    assumption: "One sentence: what directional prior this paper justifies."
```

2. Optionally mirror the same entry in `DEFAULT_EVIDENCE` inside `src/evidence_registry.py` so fresh installs without a YAML file still carry the citation.
3. Re-run a brief or open Streamlit — the new row appears automatically. No code change is required for rendering.

## Custom / newly discovered molecules

**Important:** There is no structure→affinity / SMILES path. For a self-defined or novel molecule you supply *mechanistic assumptions* and an empirical `PKProfile`. The twin only maps those stated inputs → circuit readouts.

### Required user inputs

| Field | Meaning |
|-------|---------|
| `target_node` | Circuit target (`dopamine_tone`, `d1_direct`, `d2_indirect`, `stn_hyperdirect`, `gpe`, `gpi_snr`) |
| `mechanism` | `agonist` / `antagonist` / `pam` / `nam` / `channel_blocker` |
| `max_effect` | Peak signed gain magnitude (unitless) |
| `ec50_nM` | Optional; from assay if known |
| `PKProfile` | `onset_h`, `tmax_h`, `half_life_h`, `cmax_nM`, optional multi-dose |

Affinity, BBB, and toxicity remain **user inputs or unknown**.

### Streamlit steps

1. Launch the app (`python -m streamlit run app/streamlit_app.py`).
2. Choose **Custom 1–2 candidates (form)** *or* **Custom YAML upload / paste**.
3. Set baseline (typically `pd_depleted`) and optionally keep **Also run Healthy / Untreated PD / L-DOPA references** checked.
4. Form mode: enable Candidate A/B; set target, mechanism, max_effect, optional EC50, and PK fields.
5. YAML mode: paste or upload a `MoleculePerturbation` document (example below).
6. Click **Run simulation** → compare trajectories / pathway / uncertainty tabs.
7. Export a Decision Brief from **Assumptions / Citations / Export** (includes full literature table).

### YAML example

```yaml
name: My novel candidate
description: User-defined mechanistic scenario (assay/PK assumptions required).
notes: Affinity/BBB/toxicity are user inputs or unknown — not predicted by the twin.
pk:
  onset_h: 0.2
  tmax_h: 1.0
  half_life_h: 2.5
  cmax_nM: 150.0
  n_doses: 1
modulations:
  - target_node: d1_direct
    mechanism: agonist
    max_effect: 0.7
    ec50_nM: 50.0
    hill_coeff: 1.0
```

A full benchmark file with a top-level `perturbation:` key is also accepted. Programmatically:

```python
from src.config_loader import parse_perturbation_yaml, load_baseline_params
from src.circuit_model import run_simulation
from src.schemas import SimulationConfig, BaselinePreset

pert = parse_perturbation_yaml(open("my_mol.yaml").read())
cfg = SimulationConfig(baseline_preset=BaselinePreset.PD_DEPLETED, seed=42)
params = load_baseline_params(BaselinePreset.PD_DEPLETED)
result = run_simulation(cfg, pert, params)
```

Or CLI brief from a custom benchmark YAML under `configs/benchmarks/`.

## Architecture

| Module | Role |
|--------|------|
| `src/schemas.py` | Pydantic v2 models (`MoleculePerturbation`, `PKProfile`, …) |
| `src/pk_model.py` | Exposure & Hill occupancy envelopes |
| `src/circuit_model.py` | CBGT ODE engine (`solve_ivp`) |
| `src/readouts.py` | Beta proxy, motor readiness, over-drive KPIs |
| `src/uncertainty.py` | LHS/MC P10–P90 + OAT sensitivity |
| `src/evidence_registry.py` | Citation / assumption map + References markdown |
| `src/report_generator.py` | Decision brief CLI (always includes literature table) |
| `src/config_loader.py` | Baseline/benchmark loaders + `parse_perturbation_yaml` |
| `app/streamlit_app.py` | Interactive comparison UI (benchmarks + custom form/YAML) |

## Benchmarks

Five reference scenarios under `configs/benchmarks/` (healthy, untreated PD, L-DOPA, selective D1 agonist, overdrive/dyskinesia case). Molecules in benchmarks are plain `MoleculePerturbation` YAML blocks — the same schema used for custom candidates.

## Reproducibility

Every run is deterministic given its random seed; decision briefs save exact configuration JSON and P10–P90 intervals.
