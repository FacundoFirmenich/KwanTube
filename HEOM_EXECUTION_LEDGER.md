# HEOM EXECUTION LEDGER (Complete Audit Trail)
**Session Date:** 2026-04-22  
**Audit Window:** > 12:00 AM (Local Time)

## 1. Pipeline Sanitation & Path Normalization (Boilerplate Injection)
| Target Script | Action | Status | Rationale |
|:---|:---|:---|:---|
| `check_hierarchy_convergence.py` | Route Sync | FIXED | Portability across project subdirs |
| `compare_heom_vs_redfield_1jff.py`| Fix / Refactor | FIXED | HEOM-Redfield benchmark integrity |
| `heom_pade_convergence.py` | Fix / Refactor | FIXED | Hierarchical sweep validation |
| `heom_production_1jff.py` | Fix / Checkpoint | FIXED | 51h production script stability |
| `reproduce_paper.py` | Boilerplate | FIXED | Main reproduction entry point |

## 2. Core Scientific Validations (Tier-1 Results)
| Script | Outcome | Metric | Verdict |
|:---|:---|:---|:---|
| `check_redfield_vals.py` | Audited | P_init(500fs)=0.62 | Discrepancy detected (~26%) |
| `sensitivity.py` | SUCCESS | S1(eta) = 1.0216 | **One-parameter-dominance** |
| `diagnose_ss_and_meanforce.py` | SUCCESS | KL = 0.327 nats | **Non-thermal biological state** |
| `observable_level_conv_refined.py` | SUCCESS | dFrob < 4e-3 | **NC=7 / Nk=1 Standard** |
| `jff_calibration_check.py` | SUCCESS | Ratio = 0.549 | Scalability to 8 sites confirmed |
| `sensitivity_priors.py` | SUCCESS | BF > 10 | Independent of Prior Choice |

## 3. Production & Documentation Artefacts
| Script | Outcome | Result | Description |
|:---|:---|:---|:---|
| `generate_paper_figures.py` | SUCCESS | Figs 1, 2, 3, 4 | 600 DPI High-Res Export |
| `generate_living_si.py` | SUCCESS | LIVING_SI.md | Consolidated Multi-Audit Report |
| `jff_short_validation.py` | **RUNNING** | Pending | 200fs Stress-Test (NC=8 validation) |

## 4. Technical Audit Notes
*   **Boilerplate Absolute:** Se eliminaron todas las rutas relativas `../` a favor de `PROJECT_ROOT` dinámico.
*   **Checkpoints:** Sincronización de volcados `.pkl` entre los scripts de convergencia y el auditor refinado.
*   **Liouvillian Construction:** Optimización en la construcción de ADOs para sistemas de 8 sitios (QuTiP 5.x).

---
*Certified by antigravity-agentic-coding*
