# KwanTube

**KwanTube** is a reproducible computational framework and research repository for the study of quantum dynamics and decoherence in microtubule-inspired open quantum systems. It provides the full numerical stack required to reproduce the claims presented in the associated manuscript (PRX Life / bioRxiv).

The Python package maintains the legacy import namespace `qmc_mt` for backward compatibility with established analysis pipelines.

## Reproducibility Status

- **Canonical Entry Point**: `reproduce_paper_results.py`
- **Release Runbook**: `RELEASE_RUNBOOK.md`
- **Primary Artifacts**:
  - `validation_report.json`: Machine-auditable numerical ledger (SHA-256 stamped).
  - `LIVING_SI.md`: Automated, human-readable Supplementary Information.

## Quick Start

Execute the following commands from the `git_repo/` directory.

### Fast Validation (CI Smoke Test)
```bash
python src/reproduce_paper_results.py --fast
```

### Full Repository Reproduction (Canonical)
```bash
python src/reproduce_paper_results.py
```

### Extended Detection-Power / ROC Sensitivity Sweep
```bash
python src/reproduce_paper_results.py --full-roc
```
This legacy flag generates the fixed-threshold detection-power surface used in the manuscript.

## Manuscript Figures

Generate all canonical manuscript figures using the dedicated visualization engine:
```bash
python src/generate_paper_figures.py
```
Output artifacts are saved to `figures_final/`.

Additional diagnostic figures used in the manuscript/SI may be produced by their corresponding result scripts and stored in `figures_final/`.

## Continuous Integration and Quality Assurance

The repository includes an automated CI workflow defined in `.github/workflows/ci.yml`. The validation suite ensures:

1. Successful package installation in an editable environment.
2. Unit test coverage (`tests/`).
3. Numerical consistency of the Nested Sampling engine (`src/test_ns_consistency.py`).
4. Stability of the Bayesian HEOM hierarchy (`src/bayesian_heom_hierarchy_v2.py`).
5. End-to-end reproducibility smoke tests.
6. Figure generation integrity.

## Bayesian HEOM Hierarchy (v2)

KwanTube includes a specialized hierarchical contraction layer for HEOM convergence validation in the small-$N$ regime.

- **Engine**: `src/bayesian_heom_hierarchy_v2.py`
- **Input Dataset**: `src/heom_bayes_input_current.csv`
- **Output Directory**: `heom_bayes_out_v2/`

This module applies log-scale contraction modeling and hierarchical shrinkage to summarize existing HEOM evidence. It is designed for meta-validation of convergence ledgers and does not replace the requirement for full-system HEOM production runs.

## Repository Architecture

- `src/qmc_mt/`: Core physical and statistical implementations.
- `src/reproduce_paper_results.py`: End-to-end validation pipeline.
- `src/generate_paper_figures.py`: Manuscript figure generation engine.
- `outputs_data/`: Validated numerical artifacts (.json, .npz, .pkl).
- `heom_acceptance_criteria.md`: Pre-registered validation thresholds for HEOM integration.

## License

This project is licensed under the **GNU GPLv3** (`LICENSE`).
