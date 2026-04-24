"""
sensitivity_priors.py -- Robustness analysis for Bayesian evidence.
Varies prior_sd across wide ranges to verify that BF10 remains
decisive/strong regardless of the weakly-informative prior choice.
"""
import numpy as np
import json
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qmc_mt.meta import _study_bf, StudyRecord
from qmc_mt.primary_data import BABCOCK_2024, KALRA_2024

def scan_study(study: StudyRecord, priors: list[float]) -> list[float]:
    """Helper for the orchestrator: scan BF10 over a list of prior_sd."""
    bfs = []
    for psd in priors:
        res = _study_bf(study, prior_sd=psd, n_live=400, seed=42)
        bfs.append(res.BF10)
    return bfs

def run_sensitivity():
    print("Iniciando Análisis de Sensibilidad de Priors...")
    
    # --- Babcock Sensitivity (Log-ratio scale) ---
    babcock_priors = np.logspace(-1, 1, 15)
    babcock_bfs = scan_study(BABCOCK_2024, babcock_priors.tolist())
        
    # --- Kalra Sensitivity (Seconds scale) ---
    kalra_priors = np.linspace(10, 600, 15)
    kalra_bfs = scan_study(KALRA_2024, kalra_priors.tolist())

    # --- Plotting ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot Babcock
    ax = axes[0]
    ax.plot(babcock_priors, babcock_bfs, 'o-', color='navy', lw=2)
    ax.axhline(100, color='red', ls='--', alpha=0.5, label='Decisivo (100)')
    ax.axhline(10, color='orange', ls='--', alpha=0.5, label='Fuerte (10)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Prior SD (log-ratio units)")
    ax.set_ylabel("BF10")
    ax.set_title("Sensibilidad: Babcock 2024")
    ax.legend()
    ax.grid(True, which="both", ls="-", alpha=0.2)
    
    # Plot Kalra
    ax = axes[1]
    ax.plot(kalra_priors, kalra_bfs, 's-', color='darkgreen', lw=2)
    ax.axhline(100, color='red', ls='--', alpha=0.5)
    ax.axhline(10, color='orange', ls='--', alpha=0.5)
    ax.set_yscale('log')
    ax.set_xlabel("Prior SD (seconds)")
    ax.set_ylabel("BF10")
    ax.set_title("Sensibilidad: Kalra 2024")
    ax.grid(True, which="both", ls="-", alpha=0.2)

    fig.tight_layout()
    out_dir = Path("figures_final")
    out_dir.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"prior_sensitivity.{ext}", dpi=600)
    plt.close()

    # --- Results JSON ---
    summary = {
        "babcock": {
            "priors": babcock_priors.tolist(),
            "bfs": babcock_bfs,
            "min_bf": min(babcock_bfs),
            "max_bf": max(babcock_bfs),
            "remains_decisive": all(b > 100 for b, p in zip(babcock_bfs, babcock_priors) if p >= 0.5)
        },
        "kalra": {
            "priors": kalra_priors.tolist(),
            "bfs": kalra_bfs,
            "min_bf": min(kalra_bfs),
            "max_bf": max(kalra_bfs),
            "remains_decisive_kalra": all(b > 10 for b, p in zip(kalra_bfs, kalra_priors) if 50 <= p <= 220)
        }
    }
    
    with open("prior_sensitivity.json", "w") as f:
        json.dump(summary, f, indent=2)
        
    return summary

if __name__ == "__main__":
    report = run_sensitivity()
    print("\n=== RESUMEN DE SENSIBILIDAD DEL PRIOR ===")
    print(f"Babcock 2024: BF10 min={report['babcock']['min_bf']:.1f}, max={report['babcock']['max_bf']:.1f}")
    print(f"   -> ¿Es robusto (>100)?: {'SÍ' if report['babcock']['remains_decisive'] else 'NO'}")
    
    print(f"Kalra 2024:   BF10 min={report['kalra']['min_bf']:.1f}, max={report['kalra']['max_bf']:.1f}")
    print(f"   -> ¿Es robusto (>10)?:  {'SÍ' if report['kalra']['remains_decisive_kalra'] else 'NO'}")
    
    print("\nResultados guardados en: prior_sensitivity.json y figures_final/prior_sensitivity.png")
