---
title: 'KwanTube: A calibrated framework for tracing quantum decoherence in microtubule proteins'
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
    affiliation: "1"
  - name: Pau Firmenich
    affiliation: "1"
  - name: León Firmenich
    affiliation: "1"
affiliations:
  - name: "Centro de Estudios del Sur (CEDESUR) & Universidad Nacional Arturo Jauretche (UNAJ), Argentina"
    index: 1
date: 26 March 2026
bibliography: paper.bib
---

# Summary

`KwanTube` is an open-source Python package for reproducible modelling and validation of decoherence-focused calculations in microtubule-inspired open quantum systems. The software provides a practical computational stack for rate calculations, parameter inversion, sensitivity analysis, model comparison, and reproducible report generation. It is designed to support transparent, testable workflows rather than narrative-only claims.

In the current release (`v3.5.0`), the repository also includes a small-$N$ Bayesian convergence layer for HEOM ledgers (`src/scripts/bayesian_heom_hierarchy_v2.py`). This layer summarizes already-computed convergence evidence through hierarchical contraction modelling and is intended for methods/supplementary validation.

# Statement of Need

Computational studies in this domain often suffer from one or more of the following issues: weak reproducibility, fragmented scripts, limited automated validation, and inconsistent reporting across numerical experiments. `qmc_mt` addresses this software gap by providing a coherent, test-backed pipeline that makes it straightforward to:

1. compute calibrated decoherence-related quantities from shared parameter sets,
2. run inversion and sensitivity modules under controlled seeds and settings,
3. compare model alternatives with explicit quantitative criteria,
4. regenerate machine-readable and human-readable outputs in a single command path.

# State of the field

The Python scientific ecosystem already provides robust general-purpose numerical foundations, including NumPy [@harris2020numpy], SciPy [@virtanen2020scipy], Matplotlib [@hunter2007matplotlib], and pandas [@mckinney2010pandas]. For open quantum dynamics, QuTiP is a widely used toolkit [@johansson2012qutip; @johansson2013qutip2].

`qmc_mt` is not positioned as a replacement for those libraries. Instead, it contributes a domain-specific reproducibility and validation layer for microtubule-inspired open-system studies, integrating calibration scripts, inversion/sensitivity workflows, release-grade ledgers, and HEOM-oriented contraction diagnostics in a single, test-backed pipeline [@qmcmt_zenodo_2026].

The build-vs-contribute decision was driven by scope and interface constraints. Existing general-purpose frameworks prioritize broad solver coverage and generic APIs, while this project requires tightly coupled, repository-native reproducibility components (structured validation ledgers, manuscript-aligned regeneration scripts, and pre-registered HEOM acceptance checks) that are specific to this research workflow and release process. For this reason, the software is implemented as an integration layer on top of standard scientific Python foundations rather than as a fork or extension of a single upstream toolkit.

# Software design

## Core physics engine
The software implements a TLS-oriented open-system workflow with modular components, aligned with standard open-quantum-system formalism [@breuer2002open] and HEOM-oriented validation practice [@tanimura2020heom]. It includes:
- **Vibrational Dephasing ($\Gamma_{\text{vib}}$)**: Integrates Ohmic, Drude-Lorentz, and super-Ohmic spectral densities.
- **Ionic Regulation**: Debye-Hückel fluctuations corrected by microtubule wall shielding.
- **QED Cavity Shielding**: Collective Rabi splitting mediated by conformational transition dipoles within the MT lumen.

## Analytical infrastructure
- **Global Sensitivity**: Sobol-style decomposition for key parameters.
- **Multi-temp Inversion**: Parameter recovery from temperature-dependent synthetic/benchmark data.
- **Automated Validation**: Generation of `validation_report.json` and `LIVING_SI.md` through canonical reproduction commands.
- **HEOM Bayes v2 (small-$N$)**: stable contraction-ratio inference from existing HEOM convergence ledgers.

# Research impact statement

The repository is organized around executable reproducibility rather than narrative-only reporting. The current release provides a deterministic reproduction path, CI-backed verification, and machine-auditable artifacts intended for transparent reuse, re-analysis, and extension across related computational studies [@qmcmt_zenodo_2026].

Near-term significance is supported by reviewer-ready operational evidence: a public tagged release with archival DOI [@qmcmt_zenodo_2026], automated CI checks for core pipeline paths, deterministic script entrypoints for regeneration of manuscript artifacts, and explicit guardrails that separate validation summaries from full HEOM production dynamics. This evidence is intended to make independent verification and method transfer straightforward for other computational groups.

# Reproducibility and quality control

The project is distributed as a script-first reproducible stack. Current release quality controls include:

- canonical entrypoint: `python src/scripts/reproduce_paper_results.py --fast`,
- unit tests in `tests/`,
- dedicated smoke test for Bayesian HEOM v2,
- CI workflow executing tests and smoke pipelines on push/PR, with `pytest` as the primary unit-test harness [@pytest],
- release runbook (`RELEASE_RUNBOOK.md`) for deterministic pre-release checks.

# AI usage disclosure

No generative AI system was used to generate scientific results, numerical ledgers, or benchmark outputs reported by this software paper. AI assistance was limited, when used, to editorial support; all technical claims, code paths, and references were verified by the authors before release.

# Scope and guardrails

`qmc_mt` is a software and reproducibility contribution. Its Bayesian HEOM v2 component is explicitly a **validation summary layer** over existing convergence outputs. It does **not** generate new HEOM dynamics and does **not** replace new full-system HEOM production runs.

# Acknowledgements

We thank collaborators and colleagues for technical discussions on open quantum systems, HEOM benchmarking, and reproducibility workflows.

# References
