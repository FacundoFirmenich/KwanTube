# PIPELINE_MAP.md — KwanTube v3.5.1

> **Última actualización**: 2026-05-05 | **Arquitectura**: v3.5.1

---

## 1. Árbol del Repositorio

```
KwanTube/
├── config/
│   ├── numerical_params.yaml
│   └── physics_params.yaml
├── outputs_data/
│   ├── figures_final/
│   ├── production/
│   ├── raw_csv/
│   │   ├── bath/                       # bath_params_empirical.csv, bath_params_proxy.csv
│   │   ├── compact/                    # *_compact.csv, data_registry.csv
│   │   ├── flags/                      # _done_*.flag
│   │   ├── heom_+bayesian_analysis/    # hierarchy_*.csv, extrapolated_jumps.csv,
│   │   │                               # heom_bayes_input_current.csv, level_reference_checks.csv
│   │   └── theory/                     # fisher_barrier.csv
│   ├── raw_json/
│   │   ├── audit/                      # lineage_audit.json
│   │   ├── metrics/                    # redfield_summary.json, pdb_tubulin_analysis.json,
│   │   │                               # heom_convergence_summary.json, sensitivity_sobol_final.json,
│   │   │                               # prior_sensitivity.json, sbc_report.json,
│   │   │                               # meanforce_diagnosis.json, redfield_tau_coh.json
│   │   ├── progress/                   # _progress_*.json, _done_fetch_public_data.flag
│   │   └── structural/                 # outputs_validation_report.json, validation_report.json
│   ├── raw_npz/                        # H_*.npz, heom_*.npz, window_*.npz
│   ├── raw_pkl/                        # heom_*.pkl, pade_ckpt_*.pkl, pade_refined_*.pkl
│   ├── raw_txt+md/
│   │   ├── logs/                       # execution_memory.log.txt
│   │   └── reports/                    # heom_convergence_report.txt, diagnostics_v2.txt,
│   │                                   # claim_traceability_matrix_v2.md
│   └── verification/
├── src/
│   ├── qmc_mt/                         # Paquete principal (importable)
│   │   ├── [EJECUTABLES]
│   │   │   ├── diagnose_ss_and_meanforce.py
│   │   │   ├── pdb_tubulin_analysis.py
│   │   │   ├── rank_pairs.py
│   │   │   ├── recompute_tau_coh.py
│   │   │   ├── redfield_tubulin.py
│   │   │   ├── sbc_report.py
│   │   │   ├── sensitivity.py
│   │   │   └── sensitivity_priors.py
│   │   ├── [INFRAESTRUCTURA]
│   │   │   ├── run_audit.py            # Decorador de auditoría y log
│   │   │   └── validate_integrity.py   # Motor SHA-256
│   │   └── [LIBRERÍA INTERNA]
│   │       ├── bootstrap.py, control_variates.py, core.py
│   │       ├── inversion.py, lattice.py, meta.py, model_selection.py
│   │       ├── nested_sampling.py, noneq.py, open_system.py
│   │       ├── primary_data.py, qmc.py, richardson_fragment.py
│   │       ├── robust_interp.py, roc.py, rqmc.py, sbc.py
│   │       └── test_ns_consistency.py  # [TEST]
│   └── scripts/
│       ├── analysis/
│       │   ├── assemble_master_results.py
│       │   ├── bayesian_heom_hierarchy_v2.py
│       │   └── extract_heom_production_figure.py
│       ├── data/
│       │   ├── _audit_data.py          # [INTERNO]
│       │   ├── build_hamiltonian.py
│       │   ├── build_registry.py
│       │   ├── compute_detectability_metrics.py
│       │   ├── curate_compact.py
│       │   ├── export_claim_traceability.py
│       │   ├── fetch_public_data.py
│       │   ├── run_comparative_panels.py
│       │   └── run_pipeline_vscode.py  # [ORQUESTADOR VSCode]
│       ├── figures/
│       │   ├── extract_vector_figure.py
│       │   └── generate_paper_figures.py
│       ├── heom/
│       │   ├── heom_acceptance_criteria.md
│       │   └── heom_production_driver.py
│       └── validation/
│           ├── audit_lineage.py
│           ├── heom_pade_convergence.py
│           ├── reproduce_paper_results.py  # PUNTO CANÓNICO
│           ├── seal_outputs.py
│           └── validate_outputs.py
├── tests/
│   ├── test_bayesian_heom_hierarchy_v2_smoke.py
│   ├── test_interpolation_integrity.py
│   ├── test_inversion.py
│   ├── test_parameter_traceability.py
│   └── test_physics.py
├── CITATION.cff, CONTRIBUTING.md, LICENSE
├── LIVING_SI.md                        # SI auto-generado
├── PIPELINE_MAP.md                     # Este documento
├── README.md
├── RELEASE_RUNBOOK.md
├── paper.bib / paper.md
├── pyproject.toml                      # v3.5.1
├── reproduce_results.bat / .sh / .command / _linux.desktop
└── requirements.txt
```

### Archivos legacy — pendientes de eliminación

| Archivo | Motivo |
|:--------|:-------|
| `push_biorxiv_final.ps1` | Script de push ad hoc; reemplazado por workflow Git |
| `do_push_kwantube.ps1` | Ídem; redundante con `git push` estándar |

---

## 2. Pipeline de Ejecución Secuencial

### FASE 0 — Instalación

```bash
pip install -e .
# o: pip install -r requirements.txt
```

---

### FASE 1 — Adquisición de Datos Públicos

```bash
python src/scripts/data/fetch_public_data.py
```

| Salida | Destino |
|--------|---------|
| Progress heartbeat | `raw_json/progress/_progress_fetch_public_data.json` |
| Done flag | `raw_json/progress/_done_fetch_public_data.flag` |

> **Nota sobre el repositorio de datos públicos**: Los datos descargados (~1.52 GB)
> no pueden subirse a GitHub por restricciones de tamaño. Están disponibles
> libremente en:
> **https://drive.google.com/drive/folders/1DmXBlcdwP7gY-k56RKdNGRvGSRwd5ZzU?usp=sharing**
>
> El script re-descarga todo desde las APIs originales (RCSB, PubChem, OpenAlex,
> CrossRef, Europe PMC) de forma completamente reproducible.
> Ver [`PUBLIC_DATA.md`](PUBLIC_DATA.md) para detalles.

---

### FASE 2 — Curación de Tablas Compactas

```bash
python src/scripts/data/curate_compact.py
python src/scripts/data/build_registry.py
```

| Salida | Destino |
|--------|---------|
| `structures_compact.csv`, `studies_compact.csv`, `spectral_compact.csv` | `raw_csv/compact/` |
| `data_registry.csv` | `raw_csv/compact/` |

---

### FASE 3 — Análisis PDB y Hamiltoniano

```bash
python src/qmc_mt/pdb_tubulin_analysis.py
python src/scripts/data/build_hamiltonian.py
```

| Salida | Destino |
|--------|---------|
| `pdb_tubulin_analysis.json` | `raw_json/metrics/` |
| `H_1JFF.npz`, `H_6DPU.npz` | `raw_npz/` |

---

### FASE 4 — Producción HEOM ⚠️ costosa

```bash
python src/scripts/heom/heom_production_driver.py
```

| Salida | Destino |
|--------|---------|
| `window_000.npz` … `window_095.npz` | `production/` |

> Los artefactos de producción ya están presentes. No re-ejecutar salvo cambio de Hamiltoniano.

---

### FASE 5 — Diagnósticos Físicos

```bash
python src/qmc_mt/redfield_tubulin.py
python src/qmc_mt/recompute_tau_coh.py
python src/qmc_mt/diagnose_ss_and_meanforce.py
python src/qmc_mt/rank_pairs.py
```

| Salida | Destino |
|--------|---------|
| `redfield_summary.json` | `raw_json/metrics/` |
| `redfield_tau_coh.json` | `raw_json/metrics/` |
| `meanforce_diagnosis.json` | `raw_json/metrics/` |

---

### FASE 6 — Bayesiano y Sensibilidad

```bash
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py
python src/qmc_mt/sensitivity.py
python src/qmc_mt/sensitivity_priors.py
python src/qmc_mt/sbc_report.py
```

| Salida | Destino |
|--------|---------|
| `hierarchy_global_contraction.csv`, `hierarchy_group_shrinkage.csv`, `extrapolated_jumps.csv`, `level_reference_checks.csv` | `raw_csv/heom_+bayesian_analysis/` |
| `diagnostics_v2.txt` | `raw_txt+md/reports/` |
| `sensitivity_sobol_final.json`, `prior_sensitivity.json`, `sbc_report.json` | `raw_json/metrics/` |

---

### FASE 7 — Detectabilidad y Trazabilidad

```bash
python src/scripts/data/compute_detectability_metrics.py
python src/scripts/data/run_comparative_panels.py
python src/scripts/data/export_claim_traceability.py
```

| Salida | Destino |
|--------|---------|
| `bath_params_empirical.csv` | `raw_csv/bath/` |
| `metrics_compact.csv`, `metrics_compact_summary.csv`, `comparative_panels_compact.csv` | `raw_csv/compact/` |
| `fisher_barrier.csv` | `raw_csv/theory/` |
| `_done_compute_detectability_metrics.flag` | `raw_csv/flags/` |
| `claim_traceability_matrix_v2.md` | `raw_txt+md/reports/` |

---

### FASE 8 — Validación Técnica (punto de entrada canónico)

```bash
python src/scripts/validation/heom_pade_convergence.py
python src/scripts/validation/audit_lineage.py
python src/scripts/validation/reproduce_paper_results.py --mode paper
```

| Salida | Destino |
|--------|---------|
| `heom_convergence_summary.json` | `raw_json/metrics/` |
| `heom_convergence_report.txt` | `raw_txt+md/reports/` |
| `lineage_audit.json` | `raw_json/audit/` |
| `validation_report.json` | `raw_json/structural/` |
| `LIVING_SI.md` | raíz del repo |

---

### FASE 9 — Ensamblado y Figuras

```bash
python src/scripts/analysis/assemble_master_results.py
python src/scripts/analysis/extract_heom_production_figure.py
python src/scripts/figures/generate_paper_figures.py
python src/scripts/figures/extract_vector_figure.py
```

| Salida | Destino |
|--------|---------|
| `master_results.npz` | `production/` |
| Todas las figuras del manuscrito | `figures_final/` |

---

### FASE 10 — Sellado de Integridad

```bash
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
# Resultado esperado: validated=96 bad=0
```

| Salida | Destino |
|--------|---------|
| `*.sha256` (sidecars) | junto a cada `.npz` / `.pkl` |
| `outputs_validation_report.json` | `raw_json/structural/` |

---

## 3. Orquestadores de Conveniencia

| Launcher | Plataforma | Fases cubiertas |
|:---------|:-----------|:----------------|
| `reproduce_results.sh` | Linux / macOS | 6 → 9 |
| `reproduce_results.bat` | Windows | 6 → 9 |
| `reproduce_results.command` | macOS Finder | 6 → 9 |
| `reproduce_results_linux.desktop` | Linux DE | 6 → 9 |
| `src/scripts/data/run_pipeline_vscode.py` | VSCode (cualquier OS) | 1 → 7 |

---

## 4. Suite de Tests

```bash
pytest tests/ -v
python src/qmc_mt/test_ns_consistency.py
```

| Test | Cobertura |
|:-----|:----------|
| `test_bayesian_heom_hierarchy_v2_smoke.py` | Smoke test motor Bayesiano |
| `test_interpolation_integrity.py` | Interpolación robusta |
| `test_inversion.py` | Inversión de fase cuántica |
| `test_parameter_traceability.py` | Trazabilidad de parámetros físicos |
| `test_physics.py` | Constantes y fórmulas físicas |
| `test_ns_consistency.py` | Motor Nested Sampling |

---

*Para regenerar `LIVING_SI.md` y `validation_report.json`:*
```bash
python src/scripts/validation/reproduce_paper_results.py --mode paper
```
