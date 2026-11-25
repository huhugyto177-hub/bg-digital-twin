# Changelog

This file summarizes the releases visible in git history. It does not add claims beyond those commits. Dates are the commit dates.

## 1.2.0 — 2026-11-24

Commit `142e9af` (pull request #3). `src.__version__` is `1.2.0`. The Streamlit footer prints that string.

- Getting started is the first tab. It states what the tool does, a three-step quick start, the input units, how to read the proxies, and the non-claims.
- Study design and the custom candidate open with illustrative defaults: a D1 direct-pathway agonist at max_effect 0.8, onset 0 h, Tmax 1.5 h, half-life 2 h, doses 0.5 / 1 / 2, duration 8 h, and the moderate PD baseline. Those values are a scenario timebase. EC50, Cmax, and ka stay off until assay PK is enabled.
- Baselines added beside the existing `healthy` and `pd_depleted` pair: `pd_mild`, `pd_advanced`, `lid_prone`, `hd_indirect_loss`, and `dystonia_surround`. They are directional sketches. Benchmarks `01`–`17` are unchanged, and `pd_depleted` remains their moderate reference.
- Mechanism illustrations added on the existing nodes: mGluR5 NAM on `stn_hyperdirect`; 5-HT1A as a negative `dopamine_tone` gain; muscarinic M4 PAM as a negative `d1_direct` gain; GABA-A facilitation as a negative `gpi_snr` gain; and an A2A antagonist illustration with the same sign as benchmark `11`. Magnitudes are illustrative. None of them asserts an EC50, a Cmax, or a published half-life.
- Readouts added: `dyskinesia_risk_proxy`, `gpi_snr_output_rate`, `stn_burstiness_proxy`, and `direct_indirect_balance_index`.
- Pull request #2 (`0fb8f71`) is in this release's history. A blank Duration uses the 8 h default, and invalid free-text numbers warn by field name and do not run. That commit left the package string at `1.1.0`. Pull request #3 set it to `1.2.0`.

## 1.1.0 — 2026-11-24

Commit `9159d52` (pull request #1). Study design, illustration benchmarks, and the citation registry. `src.__version__` moved from `1.0.0` to `1.1.0`.

- Partial agonists, an optional absorption rate, and an analysis window are user inputs.
- Simulations start at the drug-free steady state of the selected baseline.
- ON-time, OFF-time, beta change, and therapeutic-index values are defined as dimensionless model quantities.
- The evidence registry lists citations that carry a DOI or PMID. The 2023 eLife circuit paper is cited without a personal name.
- Benchmarks `06`–`17` were added. New scenarios stay on a normalized exposure path unless a published timing value is named in the file.
- Study design compares 3–5 doses and can export a Markdown and PDF decision brief.
- The Streamlit app leads with study design. Presets, benchmarks, and custom candidates show proxy units, P10–P90 bands, evidence captions, and CSV downloads.

## 1.0.0 — 2026-11-23

Commit `f2b8b4d`. Initial simulator, laid out for a Streamlit run from this repository root.

- Mean-field cortico–basal ganglia–thalamic model, with healthy and dopamine-depleted baselines.
- Benchmarks `01`–`05`.
- Empirical PK envelope, circuit integration, readouts, and uncertainty intervals.
- Decision-brief generator and the first Streamlit app.
- Evidence registry and the pytest suite.
