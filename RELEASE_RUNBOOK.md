# Release Runbook (non-Bayesian hierarchy scope)

This runbook documents the canonical order to validate the repository before release.

## 1) Environment

```bash
pip install -e .
```

## 2) Fast reproducibility smoke

```bash
python reproduce_paper_results.py --fast
```

Expected:
- `validation_report.json` regenerated
- `LIVING_SI.md` regenerated
- all checks in summary marked `[OK]`

## 3) Figure generation

```bash
python generate_paper_figures.py
```

Expected outputs under `git_repo/figures_final/`:
- `fig1_landscape.(pdf|png)`
- `fig2_signatures.(pdf|png)`
- `fig3_frohlich.(pdf|png)`
- `fig4_scaling.(pdf|png)`

## 4) Test suite

```bash
python -m unittest discover -s tests -v
python -m unittest test_ns_consistency.py -v
```

## 4b) Bayesian HEOM hierarchy v2 (small-N contraction layer)

```bash
python src/bayesian_heom_hierarchy_v2.py src/heom_bayes_input_current.csv --output-dir heom_bayes_out_v2
```

Expected core outputs:
- `heom_bayes_out_v2/hierarchy_global_contraction.csv`
- `heom_bayes_out_v2/hierarchy_group_shrinkage.csv`
- `heom_bayes_out_v2/diagnostics_v2.txt`

Interpretation guardrail:
- This layer summarizes existing HEOM convergence evidence via random-effects contraction.
- It does **not** replace new full-system HEOM production runs.

## 5) Editorial consistency checks

Manual checks before release:
- `paper.md` coherent with current scope and version.
- `paper.bib` contains entries for all citations used in `paper/paper.tex`.
- `CITATION.cff` version aligned with `pyproject.toml`.

## 6) Optional heavier runs

```bash
python reproduce_paper_results.py
python reproduce_paper_results.py --full-roc
```
