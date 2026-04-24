---
title: 'qmc_mt: A calibrated framework for tracing quantum decoherence in microtubule proteins'
tags:
  - Python
  - biophysics
  - quantum biology
  - microtubules
  - decoherence
  - open quantum systems
authors:
  - name: Facundo Firmenich
    orcid: 0009-0002-6578-3811
    affiliation: "1, 2"
  - name: Pau Firmenich
    affiliation: "1"
  - name: León Firmenich
    affiliation: "1"
affiliations:
  - name: "CEDESUR Research Group, Barcelona, Spain"
    index: 1
  - name: "Independent Research Collaboration in Quantum Biophysics, Barcelona, Spain"
    index: 2
date: 26 March 2026
bibliography: paper.bib
---

# Summary

`qmc_mt` is an open-source Python package for reproducible modelling and validation of decoherence-focused calculations in microtubule-inspired open quantum systems. The software provides a practical computational stack for rate calculations, parameter inversion, sensitivity analysis, model comparison, and reproducible report generation. It is designed to support transparent, testable workflows rather than narrative-only claims.

In the current release (`v3.5.0`), the repository also includes a small-$N$ Bayesian convergence layer for HEOM ledgers (`src/bayesian_heom_hierarchy_v2.py`). This layer summarizes already-computed convergence evidence through hierarchical contraction modelling and is intended for methods/supplementary validation.

# Statement of Need

Computational studies in this domain often suffer from one or more of the following issues: weak reproducibility, fragmented scripts, limited automated validation, and inconsistent reporting across numerical experiments. `qmc_mt` addresses this software gap by providing a coherent, test-backed pipeline that makes it straightforward to:

1. compute calibrated decoherence-related quantities from shared parameter sets,
2. run inversion and sensitivity modules under controlled seeds and settings,
3. compare model alternatives with explicit quantitative criteria,
4. regenerate machine-readable and human-readable outputs in a single command path.

# Key Features

## 🔬 Core Physics Engine
The software implements a TLS-oriented open-system workflow with modular components. It includes:
- **Vibrational Dephasing ($\Gamma_{\text{vib}}$)**: Integrates Ohmic, Drude-Lorentz, and super-Ohmic spectral densities.
- **Ionic Regulation**: Debye-Hückel fluctuations corrected by microtubule wall shielding.
- **QED Cavity Shielding**: Collective Rabi splitting mediated by conformational transition dipoles within the MT lumen.

## 📊 Analytical Infrastructure
- **Global Sensitivity**: Sobol-style decomposition for key parameters.
- **Multi-temp Inversion**: Parameter recovery from temperature-dependent synthetic/benchmark data.
- **Automated Validation**: Generation of `validation_report.json` and `LIVING_SI.md` through canonical reproduction commands.
- **HEOM Bayes v2 (small-$N$)**: stable contraction-ratio inference from existing HEOM convergence ledgers.

# Comparison to State of the Art

General-purpose quantum toolkits are essential but usually do not provide domain-specific validation pipelines for this problem class. `qmc_mt` contributes a domain-focused layer: reproducibility scripts, regression tests, release runbook, and convergence-validation utilities tailored to this repository’s numerical workflow.

# Reproducibility and Quality Control

The project is distributed as a script-first reproducible stack. Current release quality controls include:

- canonical entrypoint: `python reproduce_paper_results.py --fast`,
- unit tests in `tests/`,
- dedicated smoke test for Bayesian HEOM v2,
- CI workflow executing tests and smoke pipelines on push/PR,
- release runbook (`RELEASE_RUNBOOK.md`) for deterministic pre-release checks.

# Scope and Guardrails

`qmc_mt` is a software and reproducibility contribution. Its Bayesian HEOM v2 component is explicitly a **validation summary layer** over existing convergence outputs. It does **not** generate new HEOM dynamics and does **not** replace new full-system HEOM production runs.

# Acknowledgements

We thank collaborators and colleagues for technical discussions on open quantum systems, HEOM benchmarking, and reproducibility workflows.

# References
