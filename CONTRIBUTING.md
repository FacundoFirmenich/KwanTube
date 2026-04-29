# Contributing to KwanTube

Thanks for contributing to `KwanTube`.

## Scope

- Bug fixes
- Reproducibility improvements
- Documentation improvements
- Tests and CI hardening
- Performance and numerical stability improvements

## Development workflow

1. Fork and create a feature branch from `main`.
2. Keep changes small and atomic.
3. Add or update tests for behavior changes.
4. Run local checks before opening a PR:

```bash
pip install -e .
python -m unittest discover -s tests -v
python src/scripts/test_ns_consistency.py
python src/scripts/reproduce_paper_results.py --fast
```

5. Open a PR with:
   - clear motivation,
   - summary of changes,
   - validation evidence (commands and outputs),
   - impact on reproducibility artifacts (if any).

## Coding expectations

- Follow existing style and naming conventions.
- Prefer explicit, deterministic behavior.
- Avoid hidden randomness; use fixed seeds where applicable.
- Keep manuscript claims traceable to scripts and artifacts.

## Reporting issues

Use GitHub Issues and include:

- environment details,
- minimal reproduction steps,
- expected vs observed behavior,
- relevant logs or stack traces.

## Security

Do not report security vulnerabilities in public issues. See `SECURITY.md`.
