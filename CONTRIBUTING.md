# Contributing

Thanks for improving Basal Ganglia Digital Twin. This is research software. Changes to the simulator need the same care as a methods section: a reader should be able to see what was assumed, what was taken from a paper, and what the model does not claim.

## Development setup

Use Python 3.11 or 3.12. From a clone of the repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the app from the repository root:

```bash
streamlit run app/streamlit_app.py
```

The `src` package is imported from the repository root. `pytest.ini` puts that root on the Python path. You do not need an editable install.

## Tests

```bash
python -m pytest
```

GitHub Actions runs that command on Python 3.11 and 3.12 for every push and pull request (`.github/workflows/tests.yml`). Please run the suite locally before opening a pull request.

Tests check directional behavior of the shipped scenarios, the evidence registry, and that rendered project text does not carry a personal name in the citation that is stored without one. Do not weaken those checks to make a new scenario pass. Do not retune circuit constants to force a benchmark result.

## Style

Match the surrounding code.

- Python, with `from __future__ import annotations` in modules that use modern annotations.
- Four-space indentation.
- Pydantic models in `src/schemas.py` for inputs. Prefer `extra = "forbid"` so a misspelled field fails loudly.
- Keep scientific constants named and commented where they already are (`MOTOR_ON_DELTA`, `OVERDRIVE_THRESHOLD`, the analysis window).
- Readouts stay dimensionless proxies. Docstrings should say what a number is, and that it is not a clinical score or a measured spectrum.
- Do not add a dependency for a small helper that the standard library or the existing requirements already cover.

There is no separate formatter configuration in the repository. A pull request that only reformats files is harder to review than one that changes behavior on purpose.

## Pull requests

1. Branch from `main`.
2. Keep the change focused. A mechanism, a readout, and a UI redesign are three reviews.
3. Run `python -m pytest`.
4. Fill in `.github/PULL_REQUEST_TEMPLATE.md`, including the scientific-claims checklist.
5. Do not commit generated files under `outputs/`, virtual environments, or `.env` files. `outputs/` is gitignored except for `outputs/.gitkeep`.
6. Do not put personal names, institutional affiliations, employer names, email addresses, or phone numbers in documentation, UI strings, issue templates, or citations that this repository renders. The GitHub username already present in the repository URL is fine. Literature citations keep the bibliographic form already used in `configs/evidence_registry.yaml`.

## Scientific-claims rule

Every new mechanism or baseline needs a verifiable citation. A verifiable citation has a DOI or a PMID that a reviewer can open. Put it in `configs/evidence_registry.yaml` (and, if you want a fallback when that file is absent, in the default list in `src/evidence_registry.py`). State, in one sentence, the directional prior the paper supports.

Do not invent numbers. In particular, do not invent:

- EC50, IC50, Ki, or any other potency
- Cmax, half-life, Tmax, onset, bioavailability, or any other PK value
- a percent effect, an accuracy, a correlation, or a sample size

If a paper reports a timing value and the scenario uses that value, the benchmark YAML must say so in `pk_timing_notes`, the way `07_d2_agonist_pramipexole_like.yaml` names the half-life source and says that Tmax is the shared illustration timebase. If potency is not taken from a paper, omit EC50 and Cmax and keep the note `potency not specified; normalized exposure shape`.

A partial agonist needs a user-supplied `intrinsic_activity` in `(0, 1]`. Label a scenario fraction as a scenario fraction.

This project does not add:

- SMILES-to-affinity or structure-to-activity prediction
- blood–brain barrier prediction
- ADMET or toxicity prediction
- clinical-efficacy claims
- medical advice

A pull request that needs any of those does not belong in this repository. Open a discussion only if you are proposing a different tool, and say that plainly.

## Where to start

- Bug reports, feature requests, and new mechanisms each have a form under `.github/ISSUE_TEMPLATE/`.
- Behavior questions about an existing readout are often answered by `src/readouts.py` and the Context of Use tab in the app.
- Security reports go through GitHub private vulnerability reporting. See `SECURITY.md`.
