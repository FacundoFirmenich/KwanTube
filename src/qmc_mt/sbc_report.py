"""
sbc_report.py -- Simulation-Based Calibration for the Nested Sampling engine.
Validates that the evidence engine (nested_sampling.py) used for the
Bayesian meta-analysis is statistically well-calibrated.
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

from qmc_mt.nested_sampling import nested_sample
from qmc_mt.sbc import simulation_based_calibration
from qmc_mt.meta import _norm_ppf  # Sincronización exacta de la transformación de prior

# ---------- constantes del modelo de calibración ----------
SIGMA_DATA = 1.0
PRIOR_SD   = 2.0

def prior_sampler(rng: np.random.Generator) -> float:
    return float(rng.normal(0, PRIOR_SD))

def data_sampler(theta: float, rng: np.random.Generator) -> np.ndarray:
    return rng.normal(theta, SIGMA_DATA, size=1)

def posterior_sampler_ns(y: np.ndarray, L: int, rng: np.random.Generator) -> np.ndarray:
    """
    Uses Nested Sampling to draw L samples from the posterior p(theta | y).
    Sincronizado con meta.py (Acklam PPF + Likelihood constante).
    """
    val_y = float(y[0])
    s = SIGMA_DATA
    
    def loglike(theta: np.ndarray) -> float:
        t = float(theta[0])
        # Incluye la constante de normalización para match exacto con meta.py
        return -0.5 * ((val_y - t) / s) ** 2 - 0.5 * np.log(2 * np.pi * s**2)
        
    def prior_transform(u: np.ndarray) -> np.ndarray:
        return np.array([PRIOR_SD * _norm_ppf(float(u[0]))])

    # 2. Ejecutar Nested Sampling
    res = nested_sample(loglike, prior_transform, ndim=1, n_live=100, seed=rng.integers(0, 1e6))
    
    # 3. Re-muestreo por importancia para obtener L muestras de la posterior
    idx = rng.choice(len(res.samples), size=L, p=res.weights)
    return res.samples[idx].ravel()

def run_sbc_report():
    print("Iniciando Simulation-Based Calibration (SBC) para Nested Sampling...")
    print(f"Configuración: n_sim=400, L=99, prior_sd={PRIOR_SD}, sigma_data={SIGMA_DATA}")
    
    # Ejecutar SBC (esto tomará ~1-2 minutos dado que corre NS 400 veces)
    res = simulation_based_calibration(
        prior_sampler=prior_sampler,
        data_sampler=data_sampler,
        posterior_sampler=posterior_sampler_ns,
        n_sim=400,
        L=99,
        n_bins=20,
        seed=42
    )
    
    # Generar Histograma de Rangos (SI Figure)
    plt.figure(figsize=(8, 5))
    plt.bar(res.bin_edges[:-1], res.bin_counts, width=np.diff(res.bin_edges), 
            align='edge', color='skyblue', edgecolor='navy', alpha=0.7)
    plt.axhline(res.n_sim / len(res.bin_counts), color='red', linestyle='--', label='Uniformidad Ideal')
    plt.xlabel("Rango (Rank)")
    plt.ylabel("Frecuencia")
    plt.title(f"SBC Rank Histogram: Nested Sampling Engine\n(p-value = {res.p_value:.4f})")
    plt.legend()
    
    out_dir = Path("figures_final")
    out_dir.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        plt.savefig(out_dir / f"sbc_calibration_ns.{ext}", dpi=600)
    plt.close()
    
    report = {
        "engine": "nested_sampling.py",
        "n_simulations": res.n_sim,
        "n_post_samples_per_sim": res.L,
        "chi2_stat": res.chi2,
        "p_value": res.p_value,
        "is_calibrated_at_0.05": res.p_value > 0.05,
        "interpretation": "Muestreo estadísticamente consistente" if res.p_value > 0.05 else "Muestreo SESGADO (Mala calibración)"
    }
    
    with open("sbc_report.json", "w") as f:
        json.dump(report, f, indent=2)
        
    return report

if __name__ == "__main__":
    report = run_sbc_report()
    print("\n=== REPORTE DE CALIBRACIÓN SBC ===")
    print(json.dumps(report, indent=2))
    print("\nFigura guardada en: figures_final/sbc_calibration_ns.png")
