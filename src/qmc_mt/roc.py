"""Neyman-Pearson Detection Power Surface (Fixed alpha = 0.05).

This tool computes the statistical power (P_D = 1 - beta) of detecting a
spectral doublet (Rabi splitting) against a null hypothesis of a single
Lorentzian peak, across a grid of signal strengths and noise levels.

Ref Results/SI:
The detection grid should be interpreted as a Neyman--Pearson power surface 
rather than as a full ROC curve. Rows correspond to the signal split (dl), 
columns to log10(SNR), and each cell reports the detection probability 
P_D = 1 - beta at fixed false-alarm rate alpha. Thus, the matrix quantifies 
operational detectability across acquisition conditions rather than 
threshold-swept classifier performance.
"""
import numpy as np
import sys
from pathlib import Path
from scipy.stats import norm

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

def neyman_pearson_power_surface(dl_grid, snr_exp_grid, n_mc: int = 100,
                                 seed: int = 42, alpha: float = 0.05) -> dict:
    """
    Calcula la superficie de potencia de deteccion de Neyman-Pearson.

    DEFINICIONES EXPLICITAS PARA EL MANUSCRITO:
    - FILAS (Axis 0): Intensidad/separacion de senal (dl).
    - COLUMNAS (Axis 1): Nivel de ruido (log10 SNR).
    - METRICA: Probabilidad de Deteccion (P_D) = Potencia Estadistica (1 - beta).
    - CRITERIO: Detector tipo Neyman-Pearson con umbral critico z_alpha (alpha=0.05).

    Ref SI:
    The detection grid should be interpreted as a Neyman--Pearson power surface 
    rather than as a full ROC curve. Rows correspond to the signal split (dl), 
    columns to log10(SNR), and each cell reports the detection probability 
    P_D = 1 - beta at fixed false-alarm rate alpha. Thus, the matrix quantifies 
    operational detectability across acquisition conditions rather than 
    threshold-swept classifier performance.
    """
    rng = np.random.default_rng(seed)
    dl_grid      = np.asarray(dl_grid, float)
    snr_exp_grid = np.asarray(snr_exp_grid, float)
    n_dl, n_snr  = len(dl_grid), len(snr_exp_grid)

    z_alpha  = norm.isf(alpha)
    P_D_mc   = np.zeros((n_dl, n_snr))
    P_D_ana  = np.zeros((n_dl, n_snr))

    for i, dl in enumerate(dl_grid):
        for j, snr_exp in enumerate(snr_exp_grid):
            snr   = 10.0 ** snr_exp
            mu1   = dl * np.sqrt(snr / 1.0e3)
            # Monte Carlo bajo H1 (Deteccion de senal real)
            x = rng.normal(mu1, 1.0, n_mc)
            P_D_mc[i, j]  = float(np.mean(x > z_alpha))
            P_D_ana[i, j] = float(norm.sf(z_alpha - mu1))

    return {
        "dl_grid":        dl_grid.tolist(),
        "snr_exp_grid":   snr_exp_grid.tolist(),
        "P_D_grid":       P_D_mc.tolist(),
        "P_D_analytic":   P_D_ana.tolist(),
        "alpha":          alpha,
        "n_mc":           n_mc,
    }


def roc_surface(dl_grid, snr_exp_grid, n_mc: int = 100,
                seed: int = 42, alpha: float = 0.05) -> dict:
    """
    Backward-compatible alias used by reproduction scripts.
    """
    return neyman_pearson_power_surface(
        dl_grid=dl_grid,
        snr_exp_grid=snr_exp_grid,
        n_mc=n_mc,
        seed=seed,
        alpha=alpha,
    )


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    dl_vals = [0.5, 1.0, 1.5, 2.0]
    snr_exp_vals = [2.0, 2.5, 3.0, 3.5]
    
    res = neyman_pearson_power_surface(dl_vals, snr_exp_vals, n_mc=200)
    
    print("\n=== NEYMAN-PEARSON DETECTION POWER SURFACE (alpha=0.05) ===")
    print("Metrica: Probabilidad de Deteccion (1 - beta)")
    print("\nEje X: log10(SNR) ->", "  ".join([f"{v:6.1f}" for v in snr_exp_vals]))
    print("-" * 55)
    for i, dl in enumerate(dl_vals):
        row = "  ".join([f"{p:6.3f}" for p in res["P_D_grid"][i]])
        print(f"dl={dl:3.1f} | {row}")
    print("-" * 55)
    print("Nota: dl es la intensidad de senal; log10(SNR) es el nivel de ruido.")
