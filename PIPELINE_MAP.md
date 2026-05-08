# PIPELINE_MAP.md - KwanTube v3.5.1

> **Last updated**: 2026-05-08 | **Architecture**: v3.5.1 (pre-freeze)

---

## 1. Repository Tree

```
KwanTube/
├── config/
│   ├── numerical_params.yaml
│   └── physics_params.yaml
├── outputs_data/
│   ├── figures_final/              # Publication-ready PNG/PDF figures
│   ├── production/                 # HEOM windows (000-029 production, 030-095 staged v2)
│   ├── raw_csv/
│   │   ├── bath/                   # bath_params_empirical.csv, bath_params_proxy.csv
│   │   ├── compact/                # *_compact.csv, data_registry.csv
│   │   ├── flags/                  # _done_*.flag
│   │   ├── heom_+bayesian_analysis/# hierarchy_*.csv, heom_bayes_input_current.csv
│   │   ├── theory/                 # fisher_barrier.csv
│   │   └── subradiant_modes*.csv    # radiative spectra N130/N260/N520
│   ├── raw_json/
│   │   ├── audit/                  # lineage_audit.json
│   │   ├── metrics/                # KWW, HEOM/Redfield, mean-force, subradiance, SBC, Sobol
│   │   ├── nonequilibrium/          # frohlich_universal_gating_audit.json
│   │   ├── progress/               # _progress_*.json
│   │   └── structural/             # validation_report.json, outputs_validation_report.json
│   ├── interactive/                 # epistemic_graph.html
│   ├── raw_npz/                    # H_*.npz, heom_*.npz, master_results.npz
│   ├── raw_pkl/                    # heom_*.pkl, pade_ckpt_*.pkl
│   ├── raw_txt+md/
│   │   ├── logs/                   # execution_memory.log.txt
│   │   └── reports/                # diagnostics_v2.txt, claim_traceability_matrix_v2.md
│   └── verification/
├── src/
│   ├── kt_utils/                   # Cross-cutting utilities v3.5.1
│   │   ├── __init__.py
│   │   ├── paths.py                # Canonical path resolution
│   │   └── logging.py              # Centralized logger with RUN_AUDIT
│   ├── qmc_mt/                     # Main package (physics + auditing)
│   │   ├── lattice.py               # B-lattice family N130/N260/N520
│   │   ├── run_audit.py            # Centralized execution log
│   │   └── validate_integrity.py   # SHA-256 cryptographic engine
│   └── scripts/
│       ├── analysis/               # Inference and post-processing engines
│       │   ├── assemble_master_results.py
│       │   ├── bayesian_heom_hierarchy_v2.py
│       │   ├── build_hamiltonian.py
│       │   ├── build_epistemic_graph.py
│       │   ├── analyze_heom_structured_relaxation.py
│       │   ├── compute_detectability_metrics.py
│       │   ├── compute_subradiant_decay_spectrum.py [family-aware]
│       │   ├── export_claim_traceability.py
│       │   ├── fit_heom_kww_relaxation.py
│       │   ├── frohlich_universal_gating_audit.py
│       │   └── heom_pade_convergence.py
│       ├── figures/                # Paper plot generation
│       │   └── generate_paper_figures.py
│       ├── heom/                   # HEOM production drivers (expensive)
│       │   └── heom_production_driver.py
│       ├── public_data/            # External data ingestion and curation
│       │   ├── fetch_public_data.py
│       │   ├── curate_compact.py
│       │   ├── build_registry.py
│       │   └── run_pipeline_vscode.py
│       └── validation/             # Technical auditing and reproduction
│           ├── audit_lineage.py
│           ├── reproduce_paper_results.py
│           ├── seal_outputs.py
│           └── validate_outputs.py
├── tests/                          # CI validation suite
├── CITATION.cff, README.md, LIVING_SI.md
└── pyproject.toml                  # Package configuration
```

---

## 2. Sequential Execution Pipeline

### PHASE 1 — Ingestion (Public Data)
```bash
python src/scripts/public_data/fetch_public_data.py
```
> **Note**: Raw data (~1.5 GB) is downloaded via API and saved locally.
> See [`PUBLIC_DATA.md`](PUBLIC_DATA.md) for the Drive mirror link.

### PHASE 2 — Curation and Registry
```bash
python src/scripts/public_data/curate_compact.py
python src/scripts/public_data/build_registry.py
```

### PHASE 3 — Structure and Hamiltonian
```bash
python src/qmc_mt/pdb_tubulin_analysis.py
python src/scripts/analysis/build_hamiltonian.py
```

### PHASE 4 — HEOM Production (Only if necessary)
```bash
python src/scripts/heom/heom_production_driver.py
```

### PHASE 5 — Relaxation Dynamics (KWW)
```bash
python src/scripts/analysis/fit_heom_kww_relaxation.py
python src/scripts/analysis/analyze_heom_structured_relaxation.py
```
*   **Results**: $\beta \approx 0.44$ (quantum purity) and subunitary clustering $\beta=0.370$--$0.462$ across six population/purity/entropy observables, supporting structured non-Markovian distributed relaxation over a finite window without asserting thermodynamic glass transition.

### PHASE 6 — Bayesian Inference and Sensitivity
```bash
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py
python src/qmc_mt/sensitivity.py
python src/qmc_mt/sbc_report.py
```
*   **Canonical Sobol**: `src/qmc_mt/sensitivity.py` generates `outputs_data/raw_json/metrics/sensitivity_sobol_final.json` with Saltelli base `N=50000`, bootstrap `n=200` and `CI=0.95` intervals. The reproduction ledger validates this precision via `sobol_canonical_precision`.

### PHASE 7 — Optical Coherence and Subradiance
```bash
python src/qmc_mt/lattice.py
python src/scripts/analysis/compute_subradiant_decay_spectrum.py --n-layers 10
python src/scripts/analysis/compute_subradiant_decay_spectrum.py
python src/scripts/analysis/compute_subradiant_decay_spectrum.py --n-layers 40
python src/scripts/analysis/frohlich_universal_gating_audit.py
```
*   **Results**: 13-protofilament family `N130/N260/N520`, with convergent excitonic gaps, increasing `lowest_mode_ipr` and monotonic increase of free-space subradiant fraction (`70.8% -> 77.3% -> 84.8%`).
*   **Artifact contract**: `subradiant_decay_spectrum_N130.json`, `subradiant_decay_spectrum.json` (`N260`), `subradiant_decay_spectrum_N520.json` and their associated `subradiant_modes*.csv`.
*   **Fröhlich audit**: separates `L_omega=v_g/(2f_F)` from `L_gamma=v_g/(2 gamma_Hz)`; for MT at 0.1 THz yields `L_omega≈10 nm` and `L_gamma≈10 um` only if `gamma_Hz≈1e8`.

### PHASE 8 — Traceability and Convergence
```bash
python src/scripts/analysis/export_claim_traceability.py
python src/scripts/analysis/heom_pade_convergence.py
python src/scripts/analysis/build_epistemic_graph.py
```
*   **Interactive graph**: generates `outputs_data/raw_json/structural/epistemic_graph.json` and `outputs_data/interactive/epistemic_graph.html`, synchronizing claims, constraints, scope boundaries, TUR, Fröhlich linewidth-conditional and structured non-Markovian diagnostics.

### PHASE 9 — Forensic Audit and Reproduction
```bash
python src/scripts/validation/audit_lineage.py
python src/scripts/validation/reproduce_paper_results.py --mode paper
```
*   Generates `LIVING_SI.md` and `outputs_data/raw_json/structural/validation_report.json` synchronized with available data.
*   `--mode paper` is maintained as a compatible alias for manuscript reproduction mode; `--fast` executes the same contract in smoke test mode.
*   Total macro-checks and SHA-256 of the report are release-specific and must be read from `validation_report.json`, not manually copied to the manuscript before freeze.
*   The ledger includes a metacheck (`validation_ledger_self_consistent`) that validates schema, unique names, boolean statuses, non-empty details and minimum domain coverage of the validation set itself.

### PHASE 10 — Final Sealing (Cryptographic Integrity)
```bash
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
```
*   **Release target**: `bad=0` in `validate_outputs.py` after sealing. The number `validated=N` is the count of sealed binary artifacts (`*.npz`, `*.pkl`) and may change when new artifacts covered by that contract are added.
*   **Do not confuse**: macro-checks from `reproduce_paper_results.py` validate claims, tables and JSON/CSV/figure artifacts required by the manuscript; `validate_outputs.py` validates cryptographic integrity of sealed binary outputs.

---

## 3. Convenience Orchestrators

| Launcher | Use |
|:---------|:----|
| `src/scripts/public_data/run_pipeline_vscode.py` | Isolated execution of Phases 1-3 |
| `reproduce_results.sh / .bat` | Execution of analysis block (5-9) |

---

## 4. Test Suite (GitHub Actions)

The repository maintains a CI (`ci.yml`) that validates each push:
1.  **Unit Tests**: `pytest tests/`
2.  **SBC Check**: Statistical calibration of the Bayesian engine.
3.  **Smoke Tests**: Fast execution of `reproduce_paper_results.py --fast`.
4.  **Integrity**: SHA-256 signature verification on critical artifacts.
