# Repository Release Runbook

This runbook specifies the canonical sequence for validating the Quantubulin repository prior to an official release or manuscript submission.

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
python -m unittest discover -s tests -v
python -m unittest src/test_ns_consistency.py -v
```

## 5. Bayesian HEOM Hierarchy Validation (v2)

Validate the Bayesian contraction layer for existing HEOM convergence evidence.

```bash
python src/bayesian_heom_hierarchy_v2.py src/heom_bayes_input_current.csv --output-dir heom_bayes_out_v2
```

**Core Deliverables**:
- `heom_bayes_out_v2/hierarchy_global_contraction.csv`
- `heom_bayes_out_v2/hierarchy_group_shrinkage.csv`
- `heom_bayes_out_v2/diagnostics_v2.txt`

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
