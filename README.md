# KwanTube

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19744600.svg)](https://doi.org/10.5281/zenodo.19744600)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)

**KwanTube** (v3.5.1) is a reproducible computational framework for the study of quantum dynamics and decoherence in microtubule-inspired open quantum systems. It provides the full numerical stack required to reproduce all claims presented in the associated manuscript (PRX Life / bioRxiv).

The Python package maintains the legacy import namespace `qmc_mt` for backward compatibility with established analysis pipelines.

---

## Public Data Repository (L4 Empirical Layer)

The complete companion public-data repository (~1.52 GB — RCSB, PubChem, OpenAlex,
CrossRef, Europe PMC) is too large for GitHub and is freely available at:

> 🔗 **[Google Drive — KwanTube\_public\_data\_repo](https://drive.google.com/drive/folders/1DmXBlcdwP7gY-k56RKdNGRvGSRwd5ZzU?usp=sharing)**

No login required. See [`PUBLIC_DATA.md`](PUBLIC_DATA.md) for the full contents,
provenance, licences, and usage instructions.

The pipeline re-fetches all data automatically from the original open-access APIs:

```bash
python src/scripts/public_data/fetch_public_data.py
```

---

## Reproducibility Status

| Artifact | Path |
|:---------|:-----|
| Canonical entry point | `src/scripts/validation/reproduce_paper_results.py` |
| Validation ledger | `outputs_data/raw_json/structural/validation_report.json` |
| Integrity report | `outputs_data/raw_json/structural/outputs_validation_report.json` |
| Living SI | `LIVING_SI.md` |
| Pipeline map | `PIPELINE_MAP.md` |
| Reviewer quickstart | `QUICKSTART_REPRODUCE.md` |

---

## Quick Start

For reviewers and readers: the manuscript results and Living Supplementary can be reproduced in approximately 5-15 minutes using [`QUICKSTART_REPRODUCE.md`](QUICKSTART_REPRODUCE.md).

### One-shot reproduction (Phases 6–9)

The launchers below create or activate a local `.venv`, install dependencies, and run the canonical workflow:

| Platform | Launcher |
|:---------|:---------|
| Windows | `reproduce_results.bat` |
| Linux / macOS | `reproduce_results.sh` |
| macOS Finder | `reproduce_results.command` |
| Linux desktop | `reproduce_results_linux.desktop` |

On Linux/macOS, restore executable permissions after cloning:

```bash
chmod +x reproduce_results.sh reproduce_results.command reproduce_results_linux.desktop
```

---

## Universal Access — Zero-Friction Reproduction

KwanTube is designed so that **any user — regardless of operating system or technical background** — can reproduce all manuscript results.

### Option A — Docker (recommended for reviewers and non-technical users)

No Python installation required. Docker must be installed ([get Docker](https://docs.docker.com/get-docker/)).

```bash
git clone https://github.com/facundofirmenich/KwanTube.git
cd KwanTube
docker compose up --build
```

All figures and validation reports appear in your local `outputs_data/` folder when the container finishes (~5–15 min on a modern laptop).

### Option B — Makefile (Linux / macOS / WSL)

```bash
make reproduce     # Full pipeline (Phases 6–9)
make validate      # Validation only
make figures       # Figures only
make test          # Full test suite
make help          # List all targets
```

### Option C — Platform launchers (double-click)

| Platform | File |
|:---------|:-----|
| Windows | `reproduce_results.bat` |
| Linux / macOS | `reproduce_results.sh` |
| macOS Finder | `reproduce_results.command` |
| Linux desktop | `reproduce_results_linux.desktop` |

### Option D — Manual (Python users)

```bash
pip install -e .
python src/scripts/validation/reproduce_paper_results.py --mode paper
```

---

### Manual canonical entry point

```bash
python src/scripts/validation/reproduce_paper_results.py --mode paper
```

This regenerates `validation_report.json`, `LIVING_SI.md`, and all manuscript-level figures.

---

## Pipeline Overview

The full pipeline runs in **10 sequential phases**. See [`PIPELINE_MAP.md`](PIPELINE_MAP.md) for the complete reference.

| Phase | Scripts | Key outputs |
|:------|:--------|:------------|
| 0 | `pip install -e .` | Environment |
| 1 | `fetch_public_data.py` | Raw structural / spectroscopic data |
| 2 | `curate_compact.py`, `build_registry.py` | `raw_csv/compact/` |
| 3 | `pdb_tubulin_analysis.py`, `build_hamiltonian.py` | `H_1JFF.npz`, `H_6DPU.npz` |
| 4 | `heom_production_driver.py` ⚠️ | `production/window_*.npz` (96 × 30 ps) |
| 5 | `fit_heom_kww_relaxation.py` [NUEVO] | KWW fits (beta~0.44 purity) |
| 6 | `bayesian_heom_hierarchy_v2.py`, `sensitivity.py`, `sbc_report.py` | `raw_json/metrics/` |
| 7 | `compute_subradiant_decay_spectrum.py` [NUEVO] | Optical subradiance modes |
| 8 | `heom_pade_convergence.py`, `export_claim_traceability.py` | `raw_txt+md/reports/` |
| 9 | `audit_lineage.py`, **`reproduce_paper_results.py`** | `validation_report.json`, `LIVING_SI.md` |
| 10 | `seal_outputs.py`, `validate_outputs.py` | SHA-256 sidecars, Integrity check |

> Phase 4 (HEOM production) is computationally expensive. Pre-computed artifacts are included in the repository; re-execution is only required after Hamiltonian changes.

---

## Output Directory Architecture (v3.5.1)

```
outputs_data/
├── figures_final/
├── production/
├── raw_csv/
│   ├── bath/                   # Empirical bath parameters
│   ├── compact/                # Curated compact tables
│   ├── flags/                  # Pipeline completion sentinels
│   ├── heom_+bayesian_analysis/# Bayesian HEOM hierarchy outputs
│   └── theory/                 # Fisher information barrier
├── raw_json/
│   ├── audit/                  # Lineage audit
│   ├── metrics/                # Physical diagnostics (Redfield, Sobol, SBC, ...)
│   ├── progress/               # Execution heartbeats
│   └── structural/             # Validation reports
├── raw_npz/                    # Hamiltonians and HEOM trajectories
├── raw_pkl/                    # Padé checkpoints
├── raw_txt+md/
│   ├── logs/                   # Execution logs
│   └── reports/                # Convergence reports, traceability matrix
└── verification/
```

---

## Bayesian HEOM Hierarchy (v2)

KwanTube includes a hierarchical contraction layer for HEOM convergence meta-validation.

- **Engine**: `src/scripts/analysis/bayesian_heom_hierarchy_v2.py`
- **Input**: `outputs_data/raw_csv/heom_+bayesian_analysis/heom_bayes_input_current.csv`
- **Outputs**: `outputs_data/raw_csv/heom_+bayesian_analysis/` (CSV), `outputs_data/raw_txt+md/reports/diagnostics_v2.txt`

---

## Manuscript Figures

Figures are regenerated automatically by the launchers. To generate independently:

```bash
python src/scripts/figures/generate_paper_figures.py
python src/scripts/analysis/extract_heom_production_figure.py
```

Output artifacts are saved to `outputs_data/figures_final/`.

---

## Integrity Verification

To verify or regenerate all SHA-256 sidecars:

```bash
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
# Expected: validated=96 bad=0
```

---

## Testing

```bash
pytest tests/ -v
python src/qmc_mt/test_ns_consistency.py
```

CI is defined in `.github/workflows/ci.yml` and covers installation, unit tests, NS consistency, Bayesian smoke tests, and figure generation integrity.

---

## Repository Architecture

See [`PIPELINE_MAP.md`](PIPELINE_MAP.md) for the full annotated repository tree and sequential execution guide.

- `src/qmc_mt/` — Core physical and statistical implementations
- `src/scripts/` — Executable pipeline scripts (public_data/, analysis/, figures/, heom/, validation/)
- `config/` — Physical and numerical parameters (`physics_params.yaml`, `numerical_params.yaml`)
- `outputs_data/` — Validated numerical artifacts (hierarchical v3.5.1 structure)
- `tests/` — Automated test suite

---

## License

This project is licensed under the **GNU GPLv3** (`LICENSE`).

## Community and Support

- Contribution guide: `CONTRIBUTING.md`
- Code of conduct: `CODE_OF_CONDUCT.md`
- Security policy: `SECURITY.md`
- Support channels: `SUPPORT.md`
- Open maintenance plan: `OPEN_DEVELOPMENT_ROADMAP.md`
