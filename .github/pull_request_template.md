## Summary

Briefly explain what this PR changes and why.

## Type of change

- [ ] Bug fix
- [ ] Feature
- [ ] Refactor
- [ ] Documentation
- [ ] Test/CI

## Validation

List exact commands run and key outputs.

```bash
python -m unittest discover -s tests -v
python src/qmc_mt/test_ns_consistency.py
python src/scripts/reproduce_paper_results.py --fast
```

## Reproducibility impact

- [ ] No impact on reproducibility artifacts
- [ ] Updates reproducibility artifacts (describe below)

Artifacts changed:

## Checklist

- [ ] Tests added/updated where needed
- [ ] Documentation updated where needed
- [ ] No hardcoded secrets
- [ ] Version/paths/names remain consistent
