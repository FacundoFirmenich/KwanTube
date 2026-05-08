# EXECUTION_ORDER.md — KwanTube v3.5.1
# Canonical Script Execution Order

> **Version**: 3.5.1 | **Updated**: 2026-05-08  
> Official reference document. Reflects folder reorganization introduced in v3.5.1+.  
> **Current ledger**: 22/22 checks | SHA-256: `228cfd088821c5ba489b6dfc78fc5af448cf5b26855b35e1008c464d23b09273`

---

## Preface

Scripts are organized into **functional groups** respecting data dependencies:
outputs from each phase become inputs for the next.  
The last three scripts have **inverse priority** to their order:
`reproduce_paper_results` → `seal_outputs` → `validate_outputs`.

**All paths are relative to `KwanTube/`.**

---

## GROUP 0 — Installation (one-time)

```bash
pip install -e .
# or: pip install -r requirements.txt
```

---

## GROUP 1 — Public Data (fetch + curation)

> Execute **in this exact order**.  
> `fetch_public_data.py` only if data is not already in `data_downloaded_public_repos/`.

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 1.1 | `fetch_public_data.py` | `src/scripts/public_data/` | ❌ NO (data already present in `data_downloaded_public_repos/`) |
| 1.2 | `_audit_data.py` | `src/scripts/public_data/` | ✅ YES |
| 1.3 | `curate_compact.py` | `src/scripts/public_data/` | ✅ YES |
| 1.4 | `build_registry.py` | `src/scripts/public_data/` | ✅ YES |
| 1.5 | `run_comparative_panels.py` | `src/scripts/public_data/` | ✅ YES |

---

## GROUP 2 — Physics and Hamiltonian

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 2.1 | `pdb_tubulin_analysis.py` | `src/qmc_mt/` | ✅ YES |
| 2.2 | `build_hamiltonian.py` | `src/scripts/analysis/` | ✅ YES |

---

## GROUP 3 — HEOM Production (⚠️ expensive: ~58-59 h)

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 3.1 | `heom_production_driver.py` | `src/scripts/heom/` | ❌ NO (artifacts already present) |

---

## GROUP 4 — Fast Physics Diagnostics

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 4.1 | `redfield_tubulin.py` | `src/qmc_mt/` | ✅ YES |
| 4.2 | `recompute_tau_coh.py` | `src/qmc_mt/` | ✅ YES |
| 4.3 | `diagnose_ss_and_meanforce.py` | `src/qmc_mt/` | ✅ YES |
| 4.4 | `rank_pairs.py` | `src/qmc_mt/` | ✅ YES |

---

## GROUP 5 — HEOM Assembly and Relaxation

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 5.1 | `assemble_master_results.py` | `src/scripts/analysis/` | ✅ YES |
| 5.2 | `fit_heom_kww_relaxation.py` | `src/scripts/analysis/` | ✅ YES |
| 5.3 | `analyze_heom_structured_relaxation.py` | `src/scripts/analysis/` | ✅ YES |

---

## GROUP 6 — Bayesian, Sensitivity and Calibration

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 6.1 | `bayesian_heom_hierarchy_v2.py` | `src/scripts/analysis/` | ✅ YES |
| 6.2 | `heom_pade_convergence.py` | `src/scripts/analysis/` | ✅ YES |
| 6.3 | `sensitivity.py` | `src/qmc_mt/` | ✅ YES |
| 6.4 | `sensitivity_priors.py` | `src/qmc_mt/` | ✅ YES |
| 6.5 | `sbc_report.py` | `src/qmc_mt/` | ✅ YES |

---

## GROUP 7 — Detectability, Metrics, Traceability and Propagation

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 7.1 | `compute_detectability_metrics.py` | `src/scripts/analysis/` | ✅ YES |
| 7.2 | `export_claim_traceability.py` | `src/scripts/analysis/` | ✅ YES |
| 7.3 | `frohlich_universal_gating_audit.py` | `src/scripts/analysis/` | ✅ YES |
| 7.4 | `build_epistemic_graph.py` | `src/scripts/analysis/` | ✅ YES |
| 7.5 | `propagations.py` | `src/scripts/analysis/` | ✅ YES |

---

## GROUP 8 — Figures

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 8.1 | `generate_paper_figures.py` | `src/scripts/figures/` | ✅ YES |
| 8.2 | `extract_vector_figure.py` | `src/scripts/figures/` | ❌ NO (deprecated HEOM benchmark script) |

---

## GROUP 9 — Lineage and Results Validation

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 9.1 | `audit_lineage.py` | `src/scripts/validation/` | ✅ YES |

---

## GROUP 10 — Last three: reproducibility, sealing and integrity

> **Execute ALWAYS in this exact order. No exceptions.**

| Order | Script | Folder | Execute |
|:-----:|:-------|:-------|:-------:|
| 10.1 — **THIRD-TO-LAST** | `reproduce_paper_results.py` | `src/scripts/validation/` | ✅ YES |
| 10.2 — **SECOND-TO-LAST** | `seal_outputs.py` | `src/scripts/validation/` | ✅ YES |
| 10.3 — **LAST** | `validate_outputs.py` | `src/scripts/validation/` | ✅ YES |

---

## Scripts NOT to execute (under any routine circumstances)

| Script | Folder | Reason |
|:-------|:-------|:-------|
| `fetch_public_data.py` | `public_data/` | Data already in `data_downloaded_public_repos/` (1.52 GB) |
| `heom_production_driver.py` | `heom/` | ~59 h computation; artifacts already present |
| `extract_vector_figure.py` | `figures/` | Deprecated HEOM benchmark script |
| `run_pipeline_vscode.py` | `public_data/` | Convenience orchestrator; does not replace the order above |

---

## Complete execution commands (minimal flow without HEOM or fetch)

```bash
# GROUP 1 — Curation (without fetch)
python src/scripts/public_data/_audit_data.py
python src/scripts/public_data/curate_compact.py
python src/scripts/public_data/build_registry.py
python src/scripts/public_data/run_comparative_panels.py

# GROUP 2 — Hamiltonian
python src/qmc_mt/pdb_tubulin_analysis.py
python src/scripts/analysis/build_hamiltonian.py

# GROUP 4 — Physics diagnostics
python src/qmc_mt/redfield_tubulin.py
python src/qmc_mt/recompute_tau_coh.py
python src/qmc_mt/diagnose_ss_and_meanforce.py
python src/qmc_mt/rank_pairs.py

# GROUP 5 — Assembly and Relaxation
python src/scripts/analysis/assemble_master_results.py
python src/scripts/analysis/fit_heom_kww_relaxation.py
python src/scripts/analysis/analyze_heom_structured_relaxation.py

# GROUP 6 — Bayesian and sensitivity
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py
python src/scripts/analysis/heom_pade_convergence.py
python src/qmc_mt/sensitivity.py
python src/qmc_mt/sensitivity_priors.py
python src/qmc_mt/sbc_report.py

# GROUP 7 — Detectability, traceability and propagation
python src/scripts/analysis/compute_detectability_metrics.py
python src/scripts/analysis/export_claim_traceability.py
python src/scripts/analysis/frohlich_universal_gating_audit.py
python src/scripts/analysis/build_epistemic_graph.py
python src/scripts/analysis/propagations.py

# GROUP 8 — Figures
python src/scripts/figures/generate_paper_figures.py

# GROUP 9 — Lineage
python src/scripts/validation/audit_lineage.py

# GROUP 10 — Reproducibility + sealing + validation (ALWAYS THIS ORDER)
python src/scripts/validation/reproduce_paper_results.py --mode paper
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py
```

---

## Folder structure `src/` (post-reorganization v3.5.1)

```
src/
├── kt_utils/                        ← v3.5.1: cross-cutting utilities
│   ├── __init__.py
│   ├── paths.py                     # Canonical path resolution
│   └── logging.py                   # Centralized logger with RUN_AUDIT
├── qmc_mt/                          ← Main package (physics + auditing)
│   ├── core.py
│   ├── lattice.py                   # B-lattice family N130/N260/N520
│   ├── sensitivity.py
│   ├── sbc_report.py
│   ├── run_audit.py
│   └── validate_integrity.py        # SHA-256 cryptographic engine
└── scripts/
    ├── analysis/                    ← Physical and Bayesian analysis
    │   ├── analyze_heom_structured_relaxation.py  ← v3.5.1: structured non-Markovian diagnostics
    │   ├── assemble_master_results.py
    │   ├── bayesian_heom_hierarchy_v2.py
    │   ├── build_epistemic_graph.py           ← v3.5.1: interactive epistemic graph
    │   ├── build_hamiltonian.py
    │   ├── compute_detectability_metrics.py
    │   ├── export_claim_traceability.py
    │   ├── extract_heom_production_figure.py
    │   ├── fit_heom_kww_relaxation.py        ← v3.5.1: KWW fit for HEOM
    │   ├── frohlich_universal_gating_audit.py ← v3.5.1: Fröhlich dimensional audit
    │   ├── heom_pade_convergence.py
    │   └── propagations.py                    ← v3.5.1: MC error propagation
    ├── figures/
    │   ├── extract_vector_figure.py       ← DEPRECATED
    │   └── generate_paper_figures.py
    ├── heom/
    │   ├── heom_acceptance_criteria.md
    │   └── heom_production_driver.py
    ├── public_data/                       ← renamed from data/
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
