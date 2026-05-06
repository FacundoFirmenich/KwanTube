# PIPELINE_MAP.md — KwanTube v3.5.1

> **Última actualización**: 2026-05-06 | **Arquitectura**: v3.5.1 (Final bioRxiv)

---

## 1. Árbol del Repositorio

```
KwanTube/
├── config/
│   ├── numerical_params.yaml
│   └── physics_params.yaml
├── outputs_data/
│   ├── figures_final/              # Figuras PNG/PDF publication-ready
│   ├── production/                 # Ventanas HEOM (000-095)
│   ├── raw_csv/
│   │   ├── bath/                   # bath_params_empirical.csv, bath_params_proxy.csv
│   │   ├── compact/                # *_compact.csv, data_registry.csv
│   │   ├── flags/                  # _done_*.flag
│   │   ├── heom_+bayesian_analysis/# hierarchy_*.csv, heom_bayes_input_current.csv
│   │   └── theory/                 # fisher_barrier.csv
│   ├── raw_json/
│   │   ├── audit/                  # lineage_audit.json
│   │   ├── metrics/                # KWW_fit.json, subradiant_spectrum.json, SBC, Sobol
│   │   ├── progress/               # _progress_*.json
│   │   └── structural/             # outputs_validation_report.json
│   ├── raw_npz/                    # H_*.npz, heom_*.npz, master_results.npz
│   ├── raw_pkl/                    # heom_*.pkl, pade_ckpt_*.pkl
│   ├── raw_txt+md/
│   │   ├── logs/                   # execution_memory.log.txt
│   │   └── reports/                # diagnostics_v2.txt, claim_traceability_matrix_v2.md
│   └── verification/
├── src/
│   ├── qmc_mt/                     # Paquete principal (Lógica física y auditoría)
│   │   ├── run_audit.py            # Log de ejecución centralizado
│   │   └── validate_integrity.py   # Motor criptográfico SHA-256
│   └── scripts/
│       ├── analysis/               # Motores de inferencia y post-procesamiento
│       │   ├── assemble_master_results.py
│       │   ├── bayesian_heom_hierarchy_v2.py
│       │   ├── build_hamiltonian.py
│       │   ├── compute_detectability_metrics.py
│       │   ├── compute_subradiant_decay_spectrum.py [NUEVO]
│       │   ├── export_claim_traceability.py
│       │   ├── fit_heom_kww_relaxation.py [NUEVO]
│       │   └── heom_pade_convergence.py
│       ├── figures/                # Generación de plots para el paper
│       │   └── generate_paper_figures.py
│       ├── heom/                   # Drivers de producción HEOM (costosos)
│       │   └── heom_production_driver.py
│       ├── public_data/            # Ingestión y curación de datos externos
│       │   ├── fetch_public_data.py
│       │   ├── curate_compact.py
│       │   ├── build_registry.py
│       │   └── run_pipeline_vscode.py
│       └── validation/             # Auditoría técnica y reproducción
│           ├── audit_lineage.py
│           ├── reproduce_paper_results.py
│           ├── seal_outputs.py
│           └── validate_outputs.py
├── tests/                          # Suite de validación CI
├── CITATION.cff, README.md, LIVING_SI.md
└── paper.bib / paper.md            # Manuscrito fuente
```

---

## 2. Pipeline de Ejecución Secuencial

### FASE 1 — Ingestión (Datos Públicos)
```bash
python src/scripts/public_data/fetch_public_data.py
```
> **Nota**: Los datos brutos (~1.5 GB) se descargan vía API y se guardan localmente.
> Ver [`PUBLIC_DATA.md`](PUBLIC_DATA.md) para el link al mirror en Drive.

### FASE 2 — Curación y Registro
```bash
python src/scripts/public_data/curate_compact.py
python src/scripts/public_data/build_registry.py
```

### FASE 3 — Estructura y Hamiltoniano
```bash
python src/qmc_mt/pdb_tubulin_analysis.py
python src/scripts/analysis/build_hamiltonian.py
```

### FASE 4 — Producción HEOM (Solo si es necesario)
```bash
python src/scripts/heom/heom_production_driver.py
```

### FASE 5 — Dinámica de Relajación (KWW)
```bash
python src/scripts/analysis/fit_heom_kww_relaxation.py
```
*   **Resultados**: $\beta \approx 0.44$ (pureza cuántica), indicando desviación Markoviana.

### FASE 6 — Inferencia Bayesiana y Sensibilidad
```bash
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py
python src/qmc_mt/sensitivity.py
python src/qmc_mt/sbc_report.py
```

### FASE 7 — Coherencia Óptica y Subradiancia
```bash
python src/scripts/analysis/compute_subradiant_decay_spectrum.py
```
*   **Resultados**: Identificación de modos protegidos ($>77\%$ subradiantes).

### FASE 8 — Trazabilidad y Convergencia
```bash
python src/scripts/analysis/export_claim_traceability.py
python src/scripts/analysis/heom_pade_convergence.py
```

### FASE 9 — Auditoría Forense y Reproducción
```bash
python src/scripts/validation/audit_lineage.py
python src/scripts/validation/reproduce_paper_results.py --mode paper
```
*   Genera el `LIVING_SI.md` definitivo sincronizado con los datos.

### FASE 10 — Sellado Final (Integridad Criptográfica)
```bash
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
```
*   **Meta**: `validated=96 bad=0`. Este es el estado de "Release Ready".

---

## 3. Orquestadores de Conveniencia

| Launcher | Uso |
|:---------|:----|
| `src/scripts/public_data/run_pipeline_vscode.py` | Ejecución aislada de Fases 1-3 |
| `reproduce_results.sh / .bat` | Ejecución del bloque de análisis (5-9) |

---

## 4. Suite de Tests (GitHub Actions)

El repositorio mantiene un CI (`ci.yml`) que valida cada push:
1.  **Unit Tests**: `pytest tests/`
2.  **SBC Check**: Calibración estadística del motor Bayesiano.
3.  **Smoke Tests**: Ejecución rápida de `reproduce_paper_results.py --fast`.
4.  **Integrity**: Verificación de firmas SHA-256 en artefactos críticos.
