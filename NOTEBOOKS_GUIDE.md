# Notebooks Guide (polish / reproducibility)

This repository is currently script-first for deterministic CI execution.

If notebooks are added, keep these rules:

1. Store notebooks under `notebooks/`.
2. Keep source-of-truth logic in `src/qmc_mt/` and scripts, not in notebook-only cells.
3. Ensure every notebook has:
   - fixed random seeds,
   - explicit dependency cell,
   - exported artifact path under `figures_final/` or `reports/`.
4. Add a smoke check in CI for notebook execution only after runtime is stable.

Suggested notebook set:
- `01_geometry_and_hamiltonian_audit.ipynb`
- `02_redfield_baseline_checks.ipynb`
- `03_heom_convergence_visual_audit.ipynb` (non-hierarchical Bayes scope)
- `04_sensitivity_and_model_selection.ipynb`
