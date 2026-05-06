# Repository Release Runbook

This runbook specifies the canonical sequence for validating the KwanTube repository prior to an official release or manuscript submission.

## 1. Environment Setup

Initialize the environment in editable mode to ensure all package-level dependencies are resolved.

```bash
pip install -e .
```

## 2. Reproducibility Smoke Test (Fast Mode)

Execute the fast-mode validation pipeline to verify the integrity of the numerical ledger.

```bash
python scripts/reproduce_paper_results.py --fast
```

**Expected Outcomes**:
- Successful regeneration of `outputs_data/validation_report.json`.
- Successful regeneration of `LIVING_SI.md`.
- All summary checks marked as `[OK]`.

## 3. Manuscript Figure Generation

Verify that the visualization engine produces all required manuscript figures.

```bash
python scripts/generate_paper_figures.py
```

**Expected Artifacts** (`git_repo/figures_final/`):
- `fig1_landscape.(pdf|png)`
- `fig2_signatures.(pdf|png)`
- `fig3_frohlich.(pdf|png)`
- `fig4_scaling.(pdf|png)`

## 4. Comprehensive Test Suite

Execute the full unit and consistency test batteries.

```bash
python -m pytest tests/ -v
python scripts/test_ns_consistency.py
```

## 5. Bayesian HEOM Hierarchy Validation (v2)

Validate the Bayesian contraction layer for existing HEOM convergence evidence.

```bash
python src/scripts/bayesian_heom_hierarchy_v2.py src/heom_bayes_input_current.csv --output-dir heom_bayes_out_v2_ci
```

**Core Deliverables** (`outputs_data/heom_bayes_out_v2_ci/`):
- `hierarchy_global_contraction.csv`
- `hierarchy_group_shrinkage.csv`
- `diagnostics_v2.txt`
- `posterior_plots_v2.png`

## 5b. Empirical Bath-Parameter Audit (L4 Public-Data Layer)

Run the institutional public-data pipeline to extract empirical tubulin-specific
bath coupling constraints from RCSB PDB and generate the Fisher information barrier
for Experiment 4 (cavity doublet detection).

```bash
python src/scripts/scripts_data_real/run_pipeline_vscode.py
```

**Expected outputs** (`outputs_data/analysis/`):
- `structures_compact.csv` — 507 RCSB tubulin structures with Wilson B-factors
- `bath_params_empirical.csv` — `eta_proxy=0.6026`, `H_s=0.53` (N=364 structures)
- `fisher_barrier.csv` — Cramér-Rao SNR requirements for Experiment 4 (36-point grid)
- `comparative_panels_compact.csv` — U_phys per mechanism (equilibrium, Fröhlich, QED cavity, subradiance)
- `claim_traceability_matrix_v2.md` — L4 traceability ledger for PRX Life reviewers

**Physical interpretation**: The Wilson B-factor median of 48.2 Å² (vs. FMO reference
12 Å²) yields a tubulin-specific bath coupling proxy `eta_proxy ≈ 0.60`, confirming
that the conservative range `eta ∈ [0.1, 1.0]` adopted in §2.1.1 of the manuscript
is empirically grounded, not arbitrary.

## 6. Editorial Consistency Checks

Manual verification of metadata and document integrity:
- [ ] Ensure `paper.md` scope aligns with v3.5.0 features.
- [ ] Confirm `paper.bib` coverage for all software citations.
- [ ] Verify version synchronization between `CITATION.cff` and `pyproject.toml`.
- [ ] Confirm `LIVING_SI.md` hash matches `heom_acceptance_criteria.md`.

## 7. Extended Reproducibility (Optional)

For final submission, execute the high-resolution sweeps:

```bash
python scripts/reproduce_paper_results.py
python scripts/reproduce_paper_results.py --full-roc
```
