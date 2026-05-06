# EXECUTION_ORDER.md — KwanTube v3.5.1
# Orden Canónico de Ejecución de Scripts

> **Versión**: 3.5.1 | **Actualizado**: 2026-05-05  
> Documento de referencia oficial. Refleja la reorganización de carpetas efectuada en v3.5.1+.

---

## Prefacio

Los scripts están organizados en **grupos funcionales** que respetan su dependencia de datos:
los outputs de cada fase son inputs de la siguiente.  
Los tres últimos scripts tienen **prioridad inversamente proporcional** a su orden:
`reproduce_paper_results` → `seal_outputs` → `validate_outputs`.

**Todas las rutas son relativas a `KwanTube/`.**

---

## GRUPO 0 — Instalación (una sola vez)

```bash
pip install -e .
# o: pip install -r requirements.txt
```

---

## GRUPO 1 — Datos Públicos (fetch + curación)

> Ejecutar **en este orden exacto**.  
> `fetch_public_data.py` solo si la data no está ya en `data_downloaded_public_repos/`.

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 1.1 | `fetch_public_data.py` | `src/scripts/public_data/` | ❌ NO (data ya presente en `data_downloaded_public_repos/`) |
| 1.2 | `_audit_data.py` | `src/scripts/public_data/` | ✅ SÍ |
| 1.3 | `curate_compact.py` | `src/scripts/public_data/` | ✅ SÍ |
| 1.4 | `build_registry.py` | `src/scripts/public_data/` | ✅ SÍ |
| 1.5 | `run_comparative_panels.py` | `src/scripts/public_data/` | ✅ SÍ |

---

## GRUPO 2 — Física y Hamiltoniano

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 2.1 | `pdb_tubulin_analysis.py` | `src/qmc_mt/` | ✅ SÍ |
| 2.2 | `build_hamiltonian.py` | `src/scripts/analysis/` | ✅ SÍ |

---

## GRUPO 3 — Producción HEOM (⚠️ muy costoso: ~58-59 h)

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 3.1 | `heom_production_driver.py` | `src/scripts/heom/` | ❌ NO (artefactos ya presentes) |

---

## GRUPO 4 — Diagnósticos Físicos Rápidos

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 4.1 | `redfield_tubulin.py` | `src/qmc_mt/` | ✅ SÍ |
| 4.2 | `recompute_tau_coh.py` | `src/qmc_mt/` | ✅ SÍ |
| 4.3 | `diagnose_ss_and_meanforce.py` | `src/qmc_mt/` | ✅ SÍ |
| 4.4 | `rank_pairs.py` | `src/qmc_mt/` | ✅ SÍ |

---

## GRUPO 5 — Ensamblado de Resultados HEOM

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 5.1 | `assemble_master_results.py` | `src/scripts/analysis/` | ✅ SÍ |

---

## GRUPO 6 — Bayesiano, Sensibilidad y Calibración

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 6.1 | `bayesian_heom_hierarchy_v2.py` | `src/scripts/analysis/` | ✅ SÍ |
| 6.2 | `heom_pade_convergence.py` | `src/scripts/analysis/` | ✅ SÍ |
| 6.3 | `sensitivity.py` | `src/qmc_mt/` | ✅ SÍ |
| 6.4 | `sensitivity_priors.py` | `src/qmc_mt/` | ✅ SÍ |
| 6.5 | `sbc_report.py` | `src/qmc_mt/` | ✅ SÍ |

---

## GRUPO 7 — Detectabilidad, Métricas y Trazabilidad

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 7.1 | `compute_detectability_metrics.py` | `src/scripts/analysis/` | ✅ SÍ |
| 7.2 | `export_claim_traceability.py` | `src/scripts/analysis/` | ✅ SÍ |

---

## GRUPO 8 — Figuras

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 8.1 | `generate_paper_figures.py` | `src/scripts/figures/` | ✅ SÍ |
| 8.2 | `extract_vector_figure.py` | `src/scripts/figures/` | ❌ NO (dependía de script HEOM benchmark deprecado) |

---

## GRUPO 9 — Validación de Linaje y Resultados

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 9.1 | `audit_lineage.py` | `src/scripts/validation/` | ✅ SÍ |

---

## GRUPO 10 — Tres últimos: reproducibilidad, sellado e integridad

> **Ejecutar SIEMPRE en este orden exacto. Sin excepciones.**

| Orden | Script | Carpeta | Ejecutar |
|:-----:|:-------|:--------|:--------:|
| 10.1 — **ANTEPENÚLTIMO** | `reproduce_paper_results.py` | `src/scripts/validation/` | ✅ SÍ |
| 10.2 — **PENÚLTIMO** | `seal_outputs.py` | `src/scripts/validation/` | ✅ SÍ |
| 10.3 — **ÚLTIMO** | `validate_outputs.py` | `src/scripts/validation/` | ✅ SÍ |

---

## Scripts NO ejecutar (en ninguna circunstancia rutinaria)

| Script | Carpeta | Razón |
|:-------|:--------|:------|
| `fetch_public_data.py` | `public_data/` | Data ya en `data_downloaded_public_repos/` (1.52 GB) |
| `heom_production_driver.py` | `heom/` | ~59 h de cómputo; artefactos ya presentes |
| `extract_vector_figure.py` | `figures/` | Script HEOM benchmark deprecado |
| `run_pipeline_vscode.py` | `public_data/` | Orquestador de conveniencia; no reemplaza el orden arriba |

---

## Comandos de ejecución completa (flujo mínimo sin HEOM ni fetch)

```bash
# GRUPO 1 — Curación (sin fetch)
python src/scripts/public_data/_audit_data.py
python src/scripts/public_data/curate_compact.py
python src/scripts/public_data/build_registry.py
python src/scripts/public_data/run_comparative_panels.py

# GRUPO 2 — Hamiltoniano
python src/qmc_mt/pdb_tubulin_analysis.py
python src/scripts/analysis/build_hamiltonian.py

# GRUPO 4 — Diagnósticos físicos
python src/qmc_mt/redfield_tubulin.py
python src/qmc_mt/recompute_tau_coh.py
python src/qmc_mt/diagnose_ss_and_meanforce.py
python src/qmc_mt/rank_pairs.py

# GRUPO 5 — Ensamblado
python src/scripts/analysis/assemble_master_results.py

# GRUPO 6 — Bayesiano y sensibilidad
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py
python src/scripts/analysis/heom_pade_convergence.py
python src/qmc_mt/sensitivity.py
python src/qmc_mt/sensitivity_priors.py
python src/qmc_mt/sbc_report.py

# GRUPO 7 — Detectabilidad y trazabilidad
python src/scripts/analysis/compute_detectability_metrics.py
python src/scripts/analysis/export_claim_traceability.py

# GRUPO 8 — Figuras
python src/scripts/figures/generate_paper_figures.py

# GRUPO 9 — Linaje
python src/scripts/validation/audit_lineage.py

# GRUPO 10 — Reproducibilidad + sellado + validación (SIEMPRE ESTE ORDEN)
python src/scripts/validation/reproduce_paper_results.py --mode paper
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
```

---

## Estructura de carpetas `src/scripts/` (post-reorganización v3.5.1)

```
src/scripts/
├── analysis/                          ← análisis físico y bayesiano
│   ├── assemble_master_results.py
│   ├── bayesian_heom_hierarchy_v2.py
│   ├── build_hamiltonian.py           ← movido desde data/
│   ├── compute_detectability_metrics.py ← movido desde data/
│   ├── export_claim_traceability.py   ← movido desde data/
│   ├── extract_heom_production_figure.py
│   └── heom_pade_convergence.py       ← movido desde validation/
├── figures/
│   ├── extract_vector_figure.py       ← DEPRECADO
│   └── generate_paper_figures.py
├── heom/
│   ├── heom_acceptance_criteria.md
│   └── heom_production_driver.py
├── public_data/                       ← renombrado desde data/
│   ├── _audit_data.py
│   ├── build_registry.py
│   ├── curate_compact.py
│   ├── fetch_public_data.py
│   ├── run_comparative_panels.py
│   └── run_pipeline_vscode.py
└── validation/
    ├── audit_lineage.py
    ├── reproduce_paper_results.py
    ├── seal_outputs.py
    └── validate_outputs.py
```
