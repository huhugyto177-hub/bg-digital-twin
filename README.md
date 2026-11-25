# Basal Ganglia Digital Twin (v1)

Preclinical therapeutic & mechanism simulator for cortico–basal–ganglia–thalamic (CBGT) circuit translation. Grounded in classical and contemporary basal ganglia / PD / dopaminergic / LFP-beta literature (Albin–DeLong–Wichmann, Gerfen/Surmeier, STN–GPe beta work, L-DOPA PK/PD, QSP framing). A 2023 eLife paper on direct and indirect pathway interactions (doi:10.7554/eLife.87644) is **one** supporting source among many — not the sole foundation.

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
| `mechanism` | `agonist` / `partial_agonist` / `antagonist` / `pam` / `nam` / `channel_blocker` |
| `max_effect` | Peak signed gain magnitude (unitless) |
| `ec50_nM` | Optional; from assay if known |
| `PKProfile` | `onset_h`, `tmax_h`, `half_life_h`, `cmax_nM`, optional multi-dose |

Affinity, BBB, and toxicity remain **user inputs or unknown**.

### Streamlit steps

1. Launch the app (`python -m streamlit run app/streamlit_app.py`).
2. **Getting started** is the first tab: what the tool does, a 3-step quick start, input units, how to read each readout, a worked example, and the Context of Use / non-claims.
3. **Study design** is prefilled. Onset 0 h, Tmax 1.5 h, half-life 2 h, duration 8 h, doses 0.5 / 1 / 2, and reference dose 1 are an illustrative timebase, not measured PK. EC50, reference Cmax, and ka stay off until you enable assay PK and type values you have.
4. Or click a mechanistic illustration preset (L-DOPA reference, D1 PAM, D1 partial agonist, D2/D3 modulator, GPR6 inverse agonist). Each one is labeled not calibrated.
5. **Benchmarks** runs a file under `configs/benchmarks/`. Leave “Use each benchmark’s configured duration” checked unless you mean to override it. IDs 01–17 are unchanged.
6. **Custom candidate** accepts a prefilled form or YAML on the same illustrative timebase.
7. Download trajectories CSV, the KPI table, config JSON, and a one-page decision brief (Markdown and PDF).
8. **Context of Use & Limitations** lists the five non-claims and the evidence registry. The same non-claims are on Getting started, and a short Context of Use line stays above the tabs.

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
| `src/readouts.py` | Dimensionless proxies: beta, motor, ON/OFF, therapeutic index |
| `src/uncertainty.py` | LHS/MC P10–P90 + OAT sensitivity |
| `src/dose_response.py` | Exposure sweep and inverted-U flag |
| `src/study_design.py` | 3–5 dose preclinical sweeps and illustration presets |
| `src/catalog.py` | Baseline and mechanism catalogs; illustrative form defaults |
| `src/guide.py` | Getting started text, including the non-claims |
| `src/presentation.py` | Tables, CSV, and the shared matplotlib style |
| `src/evidence_registry.py` | Citation / assumption map + References markdown |
| `src/report_generator.py` | Decision brief Markdown and one-page PDF |
| `src/config_loader.py` | Baseline/benchmark loaders + `parse_perturbation_yaml` |
| `app/streamlit_app.py` | Study design (primary), benchmarks, and custom candidates |

## Readouts

All activity readouts are **dimensionless proxies**. Time readouts use the simulation clock. They are not clinical hours or measured LFP power. The default analysis window is 0.8–3.5 model hours.

| Readout | Definition |
|---------|------------|
| `beta_proxy` | STN–GPe coupling proxy, high when the depleted circuit is imbalanced |
| `motor_readiness` | Thalamic release proxy |
| `beta_power_change_pct` | Percent change of window-mean beta versus the untreated reference (negative means a reduction) |
| `beta_suppression_pct` | Percent drop of window-mean beta versus beta at t = 0 |
| `on_time_h` / `off_time_h` | Model hours inside the analysis window above / below the motor threshold (`motor[0] + 0.05`). They partition the window |
| `on_time_horizon_h` | Same threshold counted over the full record |
| `time_to_motor_onset_h` | First time the motor threshold is crossed; missing if it never is |
| `overdrive_index` | Excess-drive proxy; the flag uses a peak of 1.15 inside the window |
| `therapeutic_index_margin` | Motor benefit above the untreated reference, reduced by overdrive past 1.15 |
| `pathway_balance`, `gpi_output`, `stn_gpe_coupling` | Direct-versus-indirect balance, GPi rate, and the STN–GPe term |
| `dyskinesia_risk_proxy` | Fraction of analysis-window samples with overdrive index ≥ 1.15. Not a clinical rating |
| `gpi_snr_output_rate` | Window mean of the GPi/SNr state (dimensionless rate proxy) |
| `stn_burstiness_proxy` | Window mean of STN × max(0.4 − GPe, 0). Not a spike-burst index |
| `direct_indirect_balance_index` | Window mean of (direct − indirect) / (|direct| + |indirect|), from −1 to +1 |

Simulations start at the drug-free steady state of the selected baseline, so t = 0 is not a startup transient.

## Benchmarks

Files under `configs/benchmarks/` are mechanistic illustrations. `01`–`05` keep their original numbers. `06`–`17` do not assert an EC50 or Cmax; the note on each file is “potency not specified; normalized exposure shape”. Timing is the shared illustration timebase unless a field cites a paper (Wright 1997 half-life 12.9 h for the pramipexole-like file; Ferreira 2022 “at least doubled” applied to that timebase for the COMT file). The STN file is a channel-block proxy, not a DBS module.

| File | Illustration | Source class |
|------|----------------|--------------|
| `06_ldopa_normalized_reference` | Combined dopamine-tone, D1, and D2 | Contin 2010; Nutt 2008; Nutt & Holford 1996 |
| `07_d2_agonist_pramipexole_like` | D2/D3 agonist | Wright 1997 (half-life); Kaye 2000; Rascol 2000 |
| `08_d2_antagonist_haloperidol_like` | Negative control, beta rises | Kapur 2000; Farde 1992 |
| `09_comt_opicapone_like_adjunct` | Longer L-DOPA-like envelope | Ferreira 2022; Männistö 1999 |
| `10_maob_safinamide_like_adjunct` | Scenario-prolonged envelope, not a published half-life | Youdim 2006; Finberg 2016; Loprete 2016 |
| `11_a2a_istradefylline_like` | Indirect-pathway illustration | Schwarzschild 2006; Mizuno 2013 |
| `12_amantadine_nmda_stn` | STN channel block | Verhagen Metman 1998; Aoki 1988 |
| `13_stn_channel_block_proxy` | Stimulation-on envelope, not a drug half-life | Benazzouz 2000; McIntyre 2010 |
| `14_ldopa_repeat_dosing` | Three pulses; interval is a scenario choice | Contin 2010; Nutt 2008 |
| `15_d1_partial_agonist_tavapadon_like` | Partial agonist; intrinsic activity is a scenario fraction | Bezard 2024; Young 2020 |
| `16_d1_pam_glovadalen_like` | PAM scaled by baseline dopamine tone | Rose 2022 (method framing) |
| `17_gpr6_inverse_agonist` | Reduced indirect-pathway drive | Sun 2021; Barekatain 2024 |

The Streamlit **Study design** tab can run L-DOPA reference, D1 PAM, D1 partial agonist, D2/D3 modulator, and GPR6 inverse-agonist presets in one click. Each preset is labeled “mechanistic illustration preset, not calibrated”.

## Baselines

`healthy` and `pd_depleted` are unchanged. `pd_depleted` remains the moderate reference behind benchmarks 01–17. Added circuit states are directional sketches, not fitted concentrations:

| Preset | What it sketches | Registry ids |
|--------|------------------|--------------|
| `pd_mild` | Partial dopamine loss between healthy and the moderate reference | `bernheimer_1973`, `kish_1988`, `fearnley_1991` |
| `pd_advanced` | More severe dopamine loss than `pd_depleted` (higher GPi output and beta) | `fearnley_1991`, `kish_1988`, `bernheimer_1973` |
| `lid_prone` | Depleted tone with a sensitized direct pathway | `cenci_2007_lid`, `calabresi_2010`, `picconi_2003`, `jenner_2008` |
| `hd_indirect_loss` | Early Huntington-like loss of indirect-pathway drive | `reiner_1988`, `albin_1990_hd` |
| `dystonia_surround` | Weaker surround inhibition, not a parkinsonian beta state | `mink_2003_surround`, `vitek_2002_dystonia` |

Essential tremor is not a baseline. The rate skeleton here does not represent a cerebello–thalamic oscillator, and the dystonia file does not claim that a rate change is sufficient (Vitek 2002).

## Mechanism illustrations

These are extra mappings onto the existing nodes. Magnitudes are illustrative. None of them asserts EC50, Cmax, or a published half-life. The shared timebase remains onset 0 h, Tmax 1.5 h, half-life 2 h.

| Illustration | Node and circuit field | Why that sign | Registry id |
|--------------|------------------------|---------------|-------------|
| D1 direct-pathway agonist | `d1_direct` agonist | Positive direct-pathway gain | `gerfen_surmeier_2011` |
| mGluR5 NAM | `stn_hyperdirect` NAM | Less glutamatergic drive at STN | `conn_2005_mglu` |
| 5-HT1A agonist | `dopamine_tone` NAM | Less dopamine release from serotonin terminals. The NAM field is the circuit sign, not the receptor class | `carta_2007_5ht` |
| Muscarinic M4 PAM | `d1_direct` NAM | M4 PAM reduces direct-pathway potentiation. The model field carries that negative gain | `shen_2015_m4` |
| GABA-A facilitation in GPi/SNr | `gpi_snr` NAM | Focal GABA-A agonism lowers GPi output, so the node gain is negative | `baron_2002_gaba_gpi` |
| Adenosine A2A antagonist | `d2_indirect` agonist | Same gain sign as benchmark 11 | `schwarzschild_2006` |

## Reproducibility

Every run is deterministic given its random seed; decision briefs save exact configuration JSON and P10–P90 intervals.
