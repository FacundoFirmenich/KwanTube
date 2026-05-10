# QUICKSTART_REPRODUCE.md

Reviewer-oriented reproduction guide for KwanTube v3.5.1.1.

## Purpose

The quickstart reproduction regenerates the validation ledger and Living Supplementary from archived artifacts. The full 30 ps HEOM production trajectory is provided as a frozen benchmark artifact with SHA-256 verification due to its multi-day computational cost on commodity hardware.

## Expected runtime

- Quickstart validation: approximately 5-15 minutes on a modern laptop.
- Full 30 ps HEOM production: multi-day commodity-hardware benchmark; use archived artifacts for routine review.

## One-command launchers

### Windows

```bat
reproduce_results.bat
```

### macOS/Linux

```bash
chmod +x reproduce_results.sh reproduce_results.command reproduce_results_linux.desktop
./reproduce_results.sh
```

### Makefile

```bash
make reproduce
```

### Canonical Python entry point

```bash
python src/scripts/validation/reproduce_paper_results.py --mode paper --full-roc
```

## Expected outputs

- `LIVING_SI.md`
- `outputs_data/raw_json/structural/validation_report.json`
- `outputs_data/raw_json/structural/outputs_validation_report.json`
- final figures in `outputs_data/figures_final/`
- raw CSV/JSON/TXT diagnostics under `outputs_data/`

## Validation scope

The quickstart validates:

- equilibrium decoherence baselines
- Redfield/HEOM archived-artifact consistency
- HEOM structured-relaxation diagnostics
- Fröhlich dimensional gating audit
- Bayesian evidence and SBC calibration
- lattice radiative-family checks
- Living SI regeneration and ledger integrity

## Not regenerated from scratch by default

- The full 30 ps HEOM production trajectory is not recomputed in quickstart mode because it requires multi-day wall time on commodity hardware.
- Large public-data mirrors are not bundled in GitHub; see `PUBLIC_DATA.md` for the public-data repository.

## Release identifiers

- Version: 3.5.1.1
- Zenodo DOI: `10.5281/zenodo.19744600`
- License: GPL-3.0
- Canonical validation command: `python src/scripts/validation/reproduce_paper_results.py --mode paper --full-roc`
